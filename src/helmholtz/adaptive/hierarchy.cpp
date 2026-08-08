#include "helmholtz/adaptive/hierarchy.h"
#include "lod/quasi_interp.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace lod2d::helmholtz::adaptive {
namespace {

Eigen::SparseMatrix<double> identity_sparse(int size) {
    Eigen::SparseMatrix<double> result(size, size);
    result.setIdentity();
    return result;
}

Eigen::SparseMatrix<double> cg_to_dg(const TriMesh &mesh) {
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(3 * mesh.elems.size());
    for (int element = 0; element < static_cast<int>(mesh.elems.size()); ++element) {
        for (int local = 0; local < 3; ++local) {
            triplets.emplace_back(
                3 * element + local, mesh.elems[element][local], 1.0);
        }
    }
    Eigen::SparseMatrix<double> result(
        3 * static_cast<int>(mesh.elems.size()),
        static_cast<int>(mesh.nodes.size()));
    result.setFromTriplets(triplets.begin(), triplets.end());
    return result;
}

int refinement_increment(double parent_area, double child_area) {
    if (!(parent_area > 0.0) || !(child_area > 0.0))
        throw std::runtime_error("NVB hierarchy encountered a nonpositive element area");
    const double raw = std::log2(parent_area / child_area);
    const int increment = static_cast<int>(std::llround(raw));
    if (increment < 0 || std::abs(raw - increment) > 1e-9)
        throw std::runtime_error("NVB child area is not a dyadic fraction of its parent");
    return increment;
}

using CoordinateKey = std::pair<long long, long long>;

CoordinateKey coordinate_key(const Point2 &point) {
    return {std::llround(point.x() * 1e13), std::llround(point.y() * 1e13)};
}

void align_completion_to_reference(
    NestedFineMesh &completion,
    const TriMesh &reference) {
    const TriMesh &candidate = completion.refinement.mesh;
    if (candidate.nodes.size() != reference.nodes.size()
        || candidate.elems.size() != reference.elems.size()) {
        throw std::runtime_error(
            "adaptive coarse refinement changed the fixed fine mesh size");
    }

    std::map<CoordinateKey, int> candidate_nodes;
    for (int node = 0; node < static_cast<int>(candidate.nodes.size()); ++node) {
        if (!candidate_nodes.emplace(coordinate_key(candidate.nodes[node]), node).second)
            throw std::runtime_error("fixed fine mesh contains duplicate coordinates");
    }
    std::vector<int> reference_to_candidate_node(reference.nodes.size(), -1);
    for (int node = 0; node < static_cast<int>(reference.nodes.size()); ++node) {
        const auto found = candidate_nodes.find(coordinate_key(reference.nodes[node]));
        if (found == candidate_nodes.end()
            || (candidate.nodes[found->second] - reference.nodes[node]).norm() > 1e-12) {
            throw std::runtime_error(
                "adaptive coarse refinement changed the fixed fine mesh geometry");
        }
        reference_to_candidate_node[node] = found->second;
    }

    std::map<std::array<int, 3>, int> candidate_elements;
    for (int element = 0; element < static_cast<int>(candidate.elems.size()); ++element) {
        std::array<int, 3> key = candidate.elems[element];
        std::sort(key.begin(), key.end());
        if (!candidate_elements.emplace(key, element).second)
            throw std::runtime_error("fixed fine mesh contains duplicate elements");
    }
    std::vector<int> reference_to_candidate_element(reference.elems.size(), -1);
    for (int element = 0; element < static_cast<int>(reference.elems.size()); ++element) {
        std::array<int, 3> key{};
        for (int local = 0; local < 3; ++local)
            key[local] = reference_to_candidate_node[reference.elems[element][local]];
        std::sort(key.begin(), key.end());
        const auto found = candidate_elements.find(key);
        if (found == candidate_elements.end())
            throw std::runtime_error(
                "adaptive coarse refinement changed the fixed fine triangulation");
        reference_to_candidate_element[element] = found->second;
    }

    std::vector<Eigen::Triplet<double>> node_triplets;
    node_triplets.reserve(reference.nodes.size());
    for (int node = 0; node < static_cast<int>(reference.nodes.size()); ++node)
        node_triplets.emplace_back(node, reference_to_candidate_node[node], 1.0);
    Eigen::SparseMatrix<double> node_permutation(
        reference.nodes.size(), candidate.nodes.size());
    node_permutation.setFromTriplets(node_triplets.begin(), node_triplets.end());

    std::vector<Eigen::Triplet<double>> element_triplets;
    std::vector<Eigen::Triplet<double>> dg_triplets;
    element_triplets.reserve(reference.elems.size());
    dg_triplets.reserve(3 * reference.elems.size());
    std::vector<int> aligned_levels(reference.elems.size());
    for (int element = 0; element < static_cast<int>(reference.elems.size()); ++element) {
        const int candidate_element = reference_to_candidate_element[element];
        element_triplets.emplace_back(element, candidate_element, 1.0);
        aligned_levels[element] = completion.element_levels[candidate_element];
        for (int local = 0; local < 3; ++local) {
            const int candidate_node =
                reference_to_candidate_node[reference.elems[element][local]];
            const Triangle &triangle = candidate.elems[candidate_element];
            const auto position = std::find(
                triangle.begin(), triangle.end(), candidate_node);
            if (position == triangle.end())
                throw std::runtime_error("failed to align a fine DG basis function");
            const int candidate_local =
                static_cast<int>(std::distance(triangle.begin(), position));
            dg_triplets.emplace_back(
                3 * element + local,
                3 * candidate_element + candidate_local,
                1.0);
        }
    }

    Eigen::SparseMatrix<double> element_permutation(
        reference.elems.size(), candidate.elems.size());
    element_permutation.setFromTriplets(
        element_triplets.begin(), element_triplets.end());
    Eigen::SparseMatrix<double> dg_permutation(
        3 * reference.elems.size(), 3 * candidate.elems.size());
    dg_permutation.setFromTriplets(dg_triplets.begin(), dg_triplets.end());

    completion.refinement.P_node =
        node_permutation * completion.refinement.P_node;
    completion.refinement.P_elem =
        element_permutation * completion.refinement.P_elem;
    completion.refinement.P_dg =
        dg_permutation * completion.refinement.P_dg;
    completion.refinement.mesh = reference;
    completion.element_levels = std::move(aligned_levels);
}

} // namespace

std::vector<int> fine_element_parents(
    const Eigen::SparseMatrix<double> &prolongation,
    int fine_element_count,
    int coarse_element_count) {
    if (prolongation.rows() != fine_element_count
        || prolongation.cols() != coarse_element_count) {
        throw std::invalid_argument("element prolongation dimensions do not match the meshes");
    }
    std::vector<int> parents(fine_element_count, -1);
    for (int coarse = 0; coarse < prolongation.outerSize(); ++coarse) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(prolongation, coarse); it; ++it) {
            if (std::abs(it.value()) <= 1e-14) continue;
            if (std::abs(it.value() - 1.0) > 1e-12 || parents[it.row()] >= 0)
                throw std::runtime_error("element prolongation is not a unique parent map");
            parents[it.row()] = coarse;
        }
    }
    if (std::find(parents.begin(), parents.end(), -1) != parents.end())
        throw std::runtime_error("element prolongation leaves a fine element without a parent");
    return parents;
}

std::vector<int> refinement_child_levels(
    const TriMesh &parent_mesh,
    const std::vector<int> &parent_levels,
    const RefineOutput &refinement) {
    if (parent_levels.size() != parent_mesh.elems.size())
        throw std::invalid_argument("parent level count must match parent elements");
    const std::vector<int> parents = fine_element_parents(
        refinement.P_elem,
        static_cast<int>(refinement.mesh.elems.size()),
        static_cast<int>(parent_mesh.elems.size()));
    const std::vector<double> parent_areas = compute_area(parent_mesh);
    const std::vector<double> child_areas = compute_area(refinement.mesh);
    std::vector<int> levels(refinement.mesh.elems.size());
    for (int child = 0; child < static_cast<int>(levels.size()); ++child) {
        const int parent = parents[child];
        levels[child] = parent_levels[parent]
                      + refinement_increment(parent_areas[parent], child_areas[child]);
    }
    return levels;
}

AdaptiveMeshHierarchy::AdaptiveMeshHierarchy(
    const TriMesh &initial_mesh,
    int initial_coarse_level,
    int fine_level)
    : initial_mesh_(initial_mesh), fine_level_(fine_level) {
    if (initial_mesh.nodes.empty() || initial_mesh.elems.empty())
        throw std::invalid_argument("adaptive hierarchy initial mesh must not be empty");
    if (initial_coarse_level < 0 || fine_level <= initial_coarse_level)
        throw std::invalid_argument("adaptive levels must satisfy 0 <= H_initial < h");

    RefineOutput initial_coarse = refine_mesh_nvb(initial_mesh_, initial_coarse_level);
    coarse_mesh_ = std::move(initial_coarse.mesh);
    coarse_levels_.assign(coarse_mesh_.elems.size(), initial_coarse_level);
    coarse_element_ids_.resize(coarse_mesh_.elems.size());
    coarse_parent_ids_.assign(coarse_mesh_.elems.size(), 0);
    for (std::uint64_t &id : coarse_element_ids_) id = next_element_id_++;

    refresh_fine_embedding();
    cert_audit_mesh_ = fine_completion_.refinement.mesh;
    cert_audit_element_levels_ = fine_completion_.element_levels;
    const int fine_nodes = static_cast<int>(cert_audit_mesh_.nodes.size());
    const int fine_elements = static_cast<int>(cert_audit_mesh_.elems.size());
    fine_to_cert_audit_.mesh = cert_audit_mesh_;
    fine_to_cert_audit_.P_node = identity_sparse(fine_nodes);
    fine_to_cert_audit_.P_elem = identity_sparse(fine_elements);
    fine_to_cert_audit_.P_dg = identity_sparse(3 * fine_elements);
    refresh_coarse_to_cert_audit();
}

NestedFineMesh complete_to_fine_level(
    const TriMesh &coarse_mesh,
    const std::vector<int> &coarse_levels,
    int fine_level) {
    if (coarse_levels.size() != coarse_mesh.elems.size())
        throw std::invalid_argument("coarse level count must match coarse elements");
    if (coarse_levels.empty() || fine_level < 0)
        throw std::invalid_argument("fine completion requires nonempty levels and a valid target");

    NestedFineMesh result;
    result.refinement.mesh = coarse_mesh;
    result.refinement.P_node = identity_sparse(static_cast<int>(coarse_mesh.nodes.size()));
    result.refinement.P_elem = identity_sparse(static_cast<int>(coarse_mesh.elems.size()));
    result.refinement.P_dg = identity_sparse(3 * static_cast<int>(coarse_mesh.elems.size()));
    result.element_levels = coarse_levels;

    const int max_iterations = 4 * (fine_level + 1);
    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        std::vector<int> marked;
        for (int element = 0; element < static_cast<int>(result.element_levels.size()); ++element) {
            if (result.element_levels[element] < fine_level) marked.push_back(element);
        }
        if (marked.empty()) break;

        const TriMesh parent_mesh = result.refinement.mesh;
        const std::vector<int> parent_levels = result.element_levels;
        RefineOutput step = bisect_newest_vertex(parent_mesh, marked);
        result.element_levels = refinement_child_levels(parent_mesh, parent_levels, step);
        if (*std::max_element(result.element_levels.begin(), result.element_levels.end()) > fine_level)
            throw std::runtime_error("NVB closure refined beyond the fixed fine level");
        result.refinement.P_node = step.P_node * result.refinement.P_node;
        result.refinement.P_elem = step.P_elem * result.refinement.P_elem;
        result.refinement.P_dg = step.P_dg * result.refinement.P_dg;
        result.refinement.mesh = std::move(step.mesh);
    }

    if (result.element_levels.empty()
        || *std::min_element(result.element_levels.begin(), result.element_levels.end()) != fine_level
        || *std::max_element(result.element_levels.begin(), result.element_levels.end()) != fine_level) {
        throw std::runtime_error("failed to complete adaptive coarse mesh to the fixed fine level");
    }
    return result;
}

NestedFineMesh AdaptiveMeshHierarchy::build_nested_fine_mesh() const {
    return fine_completion_;
}

void AdaptiveMeshHierarchy::refresh_fine_embedding() {
    NestedFineMesh updated = complete_to_fine_level(
        coarse_mesh_, coarse_levels_, fine_level_);
    if (!fine_completion_.refinement.mesh.nodes.empty())
        align_completion_to_reference(updated, fine_completion_.refinement.mesh);
    fine_completion_ = std::move(updated);
    if (fine_to_cert_audit_.P_node.rows() != 0) refresh_coarse_to_cert_audit();
}

void AdaptiveMeshHierarchy::refresh_coarse_to_cert_audit() {
    coarse_to_cert_audit_ =
        fine_to_cert_audit_.P_node * fine_completion_.refinement.P_node;
    coarse_elements_to_cert_audit_ =
        fine_to_cert_audit_.P_elem * fine_completion_.refinement.P_elem;
    coarse_dg_to_cert_audit_ =
        fine_to_cert_audit_.P_dg * fine_completion_.refinement.P_dg;
    coarse_to_cert_audit_.prune(1e-15);
    coarse_to_cert_audit_.makeCompressed();
    coarse_elements_to_cert_audit_.prune(1e-15);
    coarse_elements_to_cert_audit_.makeCompressed();
    coarse_dg_to_cert_audit_.prune(1e-15);
    coarse_dg_to_cert_audit_.makeCompressed();

    fine_quasi_interpolation_ = build_quasi_interp(
        coarse_mesh_, fine_completion_.refinement.mesh,
        fine_completion_.refinement.P_dg,
        cg_to_dg(fine_completion_.refinement.mesh),
        static_cast<int>(fine_completion_.refinement.mesh.nodes.size()),
        static_cast<int>(coarse_mesh_.nodes.size()));
    cert_audit_quasi_interpolation_ = build_quasi_interp(
        coarse_mesh_, cert_audit_mesh_, coarse_dg_to_cert_audit_,
        cg_to_dg(cert_audit_mesh_),
        static_cast<int>(cert_audit_mesh_.nodes.size()),
        static_cast<int>(coarse_mesh_.nodes.size()));
}

void AdaptiveMeshHierarchy::refine(const std::vector<int> &marked_elements) {
    if (marked_elements.empty()) return;
    for (int element : marked_elements) {
        if (element < 0 || element >= static_cast<int>(coarse_mesh_.elems.size()))
            throw std::out_of_range("adaptive marked element index is out of range");
        if (coarse_levels_[element] >= fine_level_)
            throw std::invalid_argument("cannot refine a coarse element at the fixed fine level");
    }

    const TriMesh parent_mesh = coarse_mesh_;
    const std::vector<int> parent_levels = coarse_levels_;
    const std::vector<std::uint64_t> parent_ids = coarse_element_ids_;
    const std::vector<std::uint64_t> parent_parent_ids = coarse_parent_ids_;
    RefineOutput refinement = bisect_newest_vertex(parent_mesh, marked_elements);
    std::vector<int> new_levels = refinement_child_levels(parent_mesh, parent_levels, refinement);
    if (*std::max_element(new_levels.begin(), new_levels.end()) > fine_level_)
        throw std::runtime_error("adaptive NVB closure would exceed the fixed fine level");
    const std::vector<int> parents = fine_element_parents(
        refinement.P_elem,
        static_cast<int>(refinement.mesh.elems.size()),
        static_cast<int>(parent_mesh.elems.size()));

    std::vector<std::uint64_t> new_ids(new_levels.size());
    std::vector<std::uint64_t> new_parent_ids(new_levels.size());
    for (int child = 0; child < static_cast<int>(new_levels.size()); ++child) {
        const int parent = parents[child];
        if (new_levels[child] == parent_levels[parent]) {
            new_ids[child] = parent_ids[parent];
            new_parent_ids[child] = parent_parent_ids[parent];
        } else {
            new_ids[child] = next_element_id_++;
            new_parent_ids[child] = parent_ids[parent];
        }
    }

    coarse_mesh_ = std::move(refinement.mesh);
    coarse_levels_ = std::move(new_levels);
    coarse_element_ids_ = std::move(new_ids);
    coarse_parent_ids_ = std::move(new_parent_ids);
    refresh_fine_embedding();
}

void AdaptiveMeshHierarchy::refine_cert_audit_from_fine_elements(
    const std::vector<int> &marked_fine_elements) {
    if (marked_fine_elements.empty()) return;
    const int fine_element_count =
        static_cast<int>(fine_completion_.refinement.mesh.elems.size());
    std::vector<char> selected(fine_element_count, false);
    for (int element : marked_fine_elements) {
        if (element < 0 || element >= fine_element_count)
            throw std::out_of_range("cert-audit fine element index is out of range");
        selected[element] = true;
    }

    const std::vector<int> fine_parents = fine_element_parents(
        fine_to_cert_audit_.P_elem,
        static_cast<int>(cert_audit_mesh_.elems.size()),
        fine_element_count);
    std::vector<int> marked_audit_elements;
    marked_audit_elements.reserve(cert_audit_mesh_.elems.size());
    for (int element = 0; element < static_cast<int>(fine_parents.size()); ++element) {
        if (selected[fine_parents[element]]) marked_audit_elements.push_back(element);
    }
    if (marked_audit_elements.empty()) return;

    const TriMesh parent_mesh = cert_audit_mesh_;
    const std::vector<int> parent_levels = cert_audit_element_levels_;
    RefineOutput step = bisect_newest_vertex(parent_mesh, marked_audit_elements);
    cert_audit_element_levels_ = refinement_child_levels(
        parent_mesh, parent_levels, step);
    fine_to_cert_audit_.P_node =
        step.P_node * fine_to_cert_audit_.P_node;
    fine_to_cert_audit_.P_elem =
        step.P_elem * fine_to_cert_audit_.P_elem;
    fine_to_cert_audit_.P_dg =
        step.P_dg * fine_to_cert_audit_.P_dg;
    cert_audit_mesh_ = std::move(step.mesh);
    fine_to_cert_audit_.mesh = cert_audit_mesh_;
    refresh_coarse_to_cert_audit();
}

} // namespace lod2d::helmholtz::adaptive

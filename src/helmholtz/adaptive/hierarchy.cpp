#include "helmholtz/adaptive/hierarchy.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace lod2d::helmholtz::adaptive {
namespace {

Eigen::SparseMatrix<double> identity_sparse(int size) {
    Eigen::SparseMatrix<double> result(size, size);
    result.setIdentity();
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
    return complete_to_fine_level(coarse_mesh_, coarse_levels_, fine_level_);
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
}

} // namespace lod2d::helmholtz::adaptive

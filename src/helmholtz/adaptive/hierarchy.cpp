#include "helmholtz/adaptive/hierarchy.h"

#include "helmholtz/boundary.h"
#include "lod/quasi_interp.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

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
        for (int local = 0; local < 3; ++local)
            triplets.emplace_back(3 * element + local, mesh.elems[element][local], 1.0);
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

double cross(const Point2 &first, const Point2 &second) {
    return first.x() * second.y() - first.y() * second.x();
}

std::array<double, 3> barycentric_coordinates(
    const Point2 &point,
    const Point2 &first,
    const Point2 &second,
    const Point2 &third) {
    const Point2 edge_one = second - first;
    const Point2 edge_two = third - first;
    const double denominator = cross(edge_one, edge_two);
    if (std::abs(denominator) <= 1e-15)
        throw std::runtime_error("nested mesh contains a degenerate parent triangle");
    const Point2 relative = point - first;
    const double second_weight = cross(relative, edge_two) / denominator;
    const double third_weight = cross(edge_one, relative) / denominator;
    return {1.0 - second_weight - third_weight, second_weight, third_weight};
}

bool inside_triangle(const std::array<double, 3> &weights, double tolerance) {
    return weights[0] >= -tolerance
        && weights[1] >= -tolerance
        && weights[2] >= -tolerance
        && weights[0] <= 1.0 + tolerance
        && weights[1] <= 1.0 + tolerance
        && weights[2] <= 1.0 + tolerance;
}

std::array<double, 3> clean_weights(std::array<double, 3> weights) {
    for (double &weight : weights) {
        if (std::abs(weight) <= 1e-13) weight = 0.0;
        if (std::abs(weight - 1.0) <= 1e-13) weight = 1.0;
    }
    const double sum = weights[0] + weights[1] + weights[2];
    if (std::abs(sum) <= 1e-15)
        throw std::runtime_error("nested mesh interpolation produced zero barycentric weight");
    for (double &weight : weights) weight /= sum;
    return weights;
}

double total_area(const TriMesh &mesh) {
    double result = 0.0;
    for (double area : compute_area(mesh)) result += std::abs(area);
    return result;
}

double triangle_diameter(const TriMesh &mesh, const Triangle &triangle) {
    return std::max({
        (mesh.nodes[triangle[0]] - mesh.nodes[triangle[1]]).norm(),
        (mesh.nodes[triangle[1]] - mesh.nodes[triangle[2]]).norm(),
        (mesh.nodes[triangle[2]] - mesh.nodes[triangle[0]]).norm()});
}

void validate_indices(
    const std::vector<int> &indices,
    int upper_bound,
    const char *description) {
    std::set<int> unique;
    for (int index : indices) {
        if (index < 0 || index >= upper_bound)
            throw std::out_of_range(std::string(description) + " index is out of range");
        if (!unique.insert(index).second)
            throw std::invalid_argument(std::string(description) + " contains a duplicate index");
    }
}

} // namespace

RefineOutput build_nested_mesh_embedding(
    const TriMesh &parent_mesh,
    const TriMesh &child_mesh) {
    if (parent_mesh.nodes.empty() || parent_mesh.elems.empty()
        || child_mesh.nodes.empty() || child_mesh.elems.empty()) {
        throw std::invalid_argument("nested mesh embedding requires two nonempty meshes");
    }
    validate_boundary_tags(parent_mesh);
    validate_boundary_tags(child_mesh);

    const double parent_area = total_area(parent_mesh);
    const double child_area = total_area(child_mesh);
    const double area_scale = std::max({1.0, parent_area, child_area});
    if (std::abs(parent_area - child_area) > 1e-10 * area_scale)
        throw std::invalid_argument("nested meshes do not cover the same domain");
    if (!parent_mesh.boundary_edges.empty()) {
        if (child_mesh.boundary_edges.empty())
            throw std::invalid_argument("nested child mesh lost explicit boundary tags");
        for (BoundaryTag tag : {BoundaryTag::Dirichlet, BoundaryTag::Robin}) {
            const double parent_measure = boundary_measure(parent_mesh, tag);
            const double child_measure = boundary_measure(child_mesh, tag);
            const double scale = std::max({1.0, parent_measure, child_measure});
            if (std::abs(parent_measure - child_measure) > 1e-10 * scale)
                throw std::invalid_argument("nested meshes have inconsistent boundary tags");
        }
    }

    const int parent_nodes = static_cast<int>(parent_mesh.nodes.size());
    const int parent_elements = static_cast<int>(parent_mesh.elems.size());
    const int child_nodes = static_cast<int>(child_mesh.nodes.size());
    const int child_elements = static_cast<int>(child_mesh.elems.size());
    constexpr double tolerance = 2e-11;

    std::vector<int> element_parents(child_elements, -1);
    std::vector<std::array<std::array<double, 3>, 3>> element_weights(child_elements);
    for (int child = 0; child < child_elements; ++child) {
        const Triangle &child_triangle = child_mesh.elems[child];
        const Point2 centroid = (
            child_mesh.nodes[child_triangle[0]]
            + child_mesh.nodes[child_triangle[1]]
            + child_mesh.nodes[child_triangle[2]]) / 3.0;
        for (int parent = 0; parent < parent_elements; ++parent) {
            const Triangle &parent_triangle = parent_mesh.elems[parent];
            const auto centroid_weights = barycentric_coordinates(
                centroid,
                parent_mesh.nodes[parent_triangle[0]],
                parent_mesh.nodes[parent_triangle[1]],
                parent_mesh.nodes[parent_triangle[2]]);
            if (!inside_triangle(centroid_weights, tolerance)) continue;

            bool all_inside = true;
            std::array<std::array<double, 3>, 3> weights{};
            for (int local = 0; local < 3; ++local) {
                weights[local] = barycentric_coordinates(
                    child_mesh.nodes[child_triangle[local]],
                    parent_mesh.nodes[parent_triangle[0]],
                    parent_mesh.nodes[parent_triangle[1]],
                    parent_mesh.nodes[parent_triangle[2]]);
                if (!inside_triangle(weights[local], tolerance)) {
                    all_inside = false;
                    break;
                }
                weights[local] = clean_weights(weights[local]);
            }
            if (!all_inside) continue;
            if (element_parents[child] >= 0)
                throw std::invalid_argument("nested child triangle has more than one parent");
            element_parents[child] = parent;
            element_weights[child] = weights;
        }
        if (element_parents[child] < 0)
            throw std::invalid_argument("child mesh is not a conforming refinement of parent mesh");
    }

    std::vector<Eigen::Triplet<double>> element_triplets;
    std::vector<Eigen::Triplet<double>> dg_triplets;
    element_triplets.reserve(child_elements);
    dg_triplets.reserve(9 * child_elements);
    std::vector<std::map<int, double>> nodal_rows(child_nodes);
    for (int child = 0; child < child_elements; ++child) {
        const int parent = element_parents[child];
        const Triangle &parent_triangle = parent_mesh.elems[parent];
        const Triangle &child_triangle = child_mesh.elems[child];
        element_triplets.emplace_back(child, parent, 1.0);
        for (int child_local = 0; child_local < 3; ++child_local) {
            std::map<int, double> row;
            for (int parent_local = 0; parent_local < 3; ++parent_local) {
                const double weight = element_weights[child][child_local][parent_local];
                if (std::abs(weight) <= 1e-14) continue;
                row[parent_triangle[parent_local]] += weight;
                dg_triplets.emplace_back(
                    3 * child + child_local,
                    3 * parent + parent_local,
                    weight);
            }
            const int child_node = child_triangle[child_local];
            if (nodal_rows[child_node].empty()) {
                nodal_rows[child_node] = row;
            } else {
                std::set<int> columns;
                for (const auto &[column, value] : nodal_rows[child_node]) {
                    (void)value;
                    columns.insert(column);
                }
                for (const auto &[column, value] : row) {
                    (void)value;
                    columns.insert(column);
                }
                for (int column : columns) {
                    const double old_value = nodal_rows[child_node].contains(column)
                        ? nodal_rows[child_node].at(column) : 0.0;
                    const double new_value = row.contains(column) ? row.at(column) : 0.0;
                    if (std::abs(old_value - new_value) > 2e-10)
                        throw std::invalid_argument(
                            "nested nodal interpolation is inconsistent across a parent edge");
                }
            }
        }
    }

    std::vector<Eigen::Triplet<double>> node_triplets;
    node_triplets.reserve(3 * child_nodes);
    for (int child_node = 0; child_node < child_nodes; ++child_node) {
        if (nodal_rows[child_node].empty())
            throw std::invalid_argument("nested child mesh contains an unused node");
        Point2 reconstructed = Point2::Zero();
        for (const auto &[parent_node, weight] : nodal_rows[child_node]) {
            node_triplets.emplace_back(child_node, parent_node, weight);
            reconstructed += weight * parent_mesh.nodes[parent_node];
        }
        if ((reconstructed - child_mesh.nodes[child_node]).norm() > 2e-10)
            throw std::runtime_error("nested nodal prolongation does not reproduce coordinates");
    }

    Eigen::SparseMatrix<double> node_prolongation(child_nodes, parent_nodes);
    node_prolongation.setFromTriplets(node_triplets.begin(), node_triplets.end());
    Eigen::SparseMatrix<double> element_prolongation(child_elements, parent_elements);
    element_prolongation.setFromTriplets(element_triplets.begin(), element_triplets.end());
    Eigen::SparseMatrix<double> dg_prolongation(3 * child_elements, 3 * parent_elements);
    dg_prolongation.setFromTriplets(dg_triplets.begin(), dg_triplets.end());
    return {
        child_mesh,
        std::move(node_prolongation),
        std::move(element_prolongation),
        std::move(dg_prolongation)};
}

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

Eigen::SparseMatrix<double> restrict_constraint_columns(
    const Eigen::SparseMatrix<double> &constraints,
    const std::vector<int> &ambient_dofs) {
    validate_indices(
        ambient_dofs,
        static_cast<int>(constraints.cols()),
        "kernel constraint degree of freedom");
    std::vector<Eigen::Triplet<double>> triplets;
    for (int local = 0; local < static_cast<int>(ambient_dofs.size()); ++local) {
        const int ambient = ambient_dofs[local];
        for (Eigen::SparseMatrix<double>::InnerIterator it(constraints, ambient); it; ++it)
            triplets.emplace_back(it.row(), local, it.value());
    }
    Eigen::SparseMatrix<double> result(constraints.rows(), ambient_dofs.size());
    result.setFromTriplets(triplets.begin(), triplets.end());
    return result;
}

AdaptiveMeshHierarchy::AdaptiveMeshHierarchy(
    const TriMesh &initial_mesh,
    int initial_coarse_level,
    int fine_level)
    : initial_mesh_(initial_mesh), fine_level_(fine_level) {
    if (initial_mesh.nodes.empty() || initial_mesh.elems.empty())
        throw std::invalid_argument("adaptive hierarchy initial mesh must not be empty");
    if (initial_coarse_level < 0 || fine_level < initial_coarse_level)
        throw std::invalid_argument("adaptive levels must satisfy 0 <= H_initial <= h");

    RefineOutput initial_coarse = refine_mesh_nvb(initial_mesh_, initial_coarse_level);
    coarse_mesh_ = std::move(initial_coarse.mesh);
    coarse_levels_.assign(coarse_mesh_.elems.size(), initial_coarse_level);
    coarse_element_ids_.resize(coarse_mesh_.elems.size());
    coarse_parent_ids_.assign(coarse_mesh_.elems.size(), 0);
    for (std::uint64_t &id : coarse_element_ids_) id = next_element_id_++;

    fine_completion_ = complete_to_fine_level(coarse_mesh_, coarse_levels_, fine_level_);
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

NestedFineMesh AdaptiveMeshHierarchy::build_nested_fine_mesh() const {
    return fine_completion_;
}

void AdaptiveMeshHierarchy::refresh_fine_embedding() {
    if (fine_completion_.refinement.mesh.nodes.empty()) {
        fine_completion_ = complete_to_fine_level(coarse_mesh_, coarse_levels_, fine_level_);
    } else {
        RefineOutput embedding = build_nested_mesh_embedding(
            coarse_mesh_, fine_completion_.refinement.mesh);
        fine_completion_.refinement = std::move(embedding);
    }
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

    fine_parent_coarse_elements_ = fine_element_parents(
        fine_completion_.refinement.P_elem,
        static_cast<int>(fine_completion_.refinement.mesh.elems.size()),
        static_cast<int>(coarse_mesh_.elems.size()));
    cert_audit_parent_fine_elements_ = fine_element_parents(
        fine_to_cert_audit_.P_elem,
        static_cast<int>(cert_audit_mesh_.elems.size()),
        static_cast<int>(fine_completion_.refinement.mesh.elems.size()));

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
    validate_current_embeddings();
}

void AdaptiveMeshHierarchy::validate_current_embeddings() const {
    constexpr double composition_tolerance = 1e-10;
    constexpr double right_inverse_tolerance = 1e-9;
    constexpr double coordinate_tolerance = 2e-10;

    validate_boundary_tags(coarse_mesh_);
    validate_boundary_tags(fine_completion_.refinement.mesh);
    validate_boundary_tags(cert_audit_mesh_);

    const Eigen::SparseMatrix<double> composed_nodes =
        fine_to_cert_audit_.P_node * fine_completion_.refinement.P_node;
    const Eigen::SparseMatrix<double> composed_elements =
        fine_to_cert_audit_.P_elem * fine_completion_.refinement.P_elem;
    const Eigen::SparseMatrix<double> composed_dg =
        fine_to_cert_audit_.P_dg * fine_completion_.refinement.P_dg;
    if ((composed_nodes - coarse_to_cert_audit_).norm() > composition_tolerance
        || (composed_elements - coarse_elements_to_cert_audit_).norm()
               > composition_tolerance
        || (composed_dg - coarse_dg_to_cert_audit_).norm()
               > composition_tolerance) {
        throw std::runtime_error(
            "three-level prolongations fail the production composition check");
    }

    auto coordinates = [](const TriMesh &mesh) {
        Eigen::MatrixXd result(mesh.nodes.size(), 2);
        for (int node = 0; node < static_cast<int>(mesh.nodes.size()); ++node)
            result.row(node) = mesh.nodes[node].transpose();
        return result;
    };
    auto maximum_absolute_entry = [](const Eigen::MatrixXd &matrix) {
        return matrix.size() == 0 ? 0.0 : matrix.cwiseAbs().maxCoeff();
    };
    const Eigen::MatrixXd coarse_coordinates = coordinates(coarse_mesh_);
    const Eigen::MatrixXd fine_coordinates =
        coordinates(fine_completion_.refinement.mesh);
    const Eigen::MatrixXd audit_coordinates = coordinates(cert_audit_mesh_);
    if (maximum_absolute_entry(
            fine_completion_.refinement.P_node * coarse_coordinates
            - fine_coordinates) > coordinate_tolerance
        || maximum_absolute_entry(
            fine_to_cert_audit_.P_node * fine_coordinates
            - audit_coordinates) > coordinate_tolerance
        || maximum_absolute_entry(
            coarse_to_cert_audit_ * coarse_coordinates
            - audit_coordinates) > coordinate_tolerance) {
        throw std::runtime_error(
            "three-level prolongations fail the production coordinate check");
    }

    std::vector<char> is_dirichlet(coarse_mesh_.nodes.size(), false);
    for (int node : dirichlet_nodes(coarse_mesh_)) is_dirichlet[node] = true;
    std::vector<Eigen::Triplet<double>> expected_triplets;
    expected_triplets.reserve(coarse_mesh_.nodes.size());
    for (int node = 0; node < static_cast<int>(coarse_mesh_.nodes.size()); ++node) {
        if (!is_dirichlet[node]) expected_triplets.emplace_back(node, node, 1.0);
    }
    Eigen::SparseMatrix<double> expected_right_inverse(
        coarse_mesh_.nodes.size(), coarse_mesh_.nodes.size());
    expected_right_inverse.setFromTriplets(
        expected_triplets.begin(), expected_triplets.end());
    const Eigen::SparseMatrix<double> fine_right_inverse =
        fine_quasi_interpolation_ * fine_completion_.refinement.P_node;
    const Eigen::SparseMatrix<double> audit_right_inverse =
        cert_audit_quasi_interpolation_ * coarse_to_cert_audit_;
    if ((fine_right_inverse - expected_right_inverse).norm()
            > right_inverse_tolerance
        || (audit_right_inverse - expected_right_inverse).norm()
            > right_inverse_tolerance) {
        throw std::runtime_error(
            "three-level interpolation fails the production right-inverse check");
    }
}

void AdaptiveMeshHierarchy::refine(const std::vector<int> &marked_elements) {
    if (marked_elements.empty()) return;
    validate_indices(
        marked_elements,
        static_cast<int>(coarse_mesh_.elems.size()),
        "adaptive marked element");

    const TriMesh parent_mesh = coarse_mesh_;
    const std::vector<int> parent_levels = coarse_levels_;
    const std::vector<std::uint64_t> parent_ids = coarse_element_ids_;
    const std::vector<std::uint64_t> parent_parent_ids = coarse_parent_ids_;
    RefineOutput refinement = bisect_newest_vertex(parent_mesh, marked_elements);
    std::vector<int> new_levels = refinement_child_levels(parent_mesh, parent_levels, refinement);
    const std::vector<int> parents = fine_element_parents(
        refinement.P_elem,
        static_cast<int>(refinement.mesh.elems.size()),
        static_cast<int>(parent_mesh.elems.size()));

    // Keep at least one local NVB generation between a changed coarse parent
    // and all of its fine descendants.  Previously the old fine mesh was used
    // as a hard capacity ceiling: once H caught h, the embedding rebuild threw
    // and adaptive refinement stopped.  Extend the master fine mesh here;
    // refine_fine_elements also advances the cert-audit mesh as needed.
    std::vector<int> required_fine_level(parent_mesh.elems.size(), -1);
    for (int child = 0; child < static_cast<int>(new_levels.size()); ++child) {
        const int parent = parents[child];
        if (new_levels[child] > parent_levels[parent]) {
            required_fine_level[parent] = std::max(
                required_fine_level[parent], new_levels[child] + 1);
        }
    }
    for (int iteration = 0; iteration < 64; ++iteration) {
        std::vector<int> marked_fine;
        for (int fine = 0;
             fine < static_cast<int>(fine_parent_coarse_elements_.size());
             ++fine) {
            const int parent = fine_parent_coarse_elements_[fine];
            if (required_fine_level[parent] >= 0
                && fine_completion_.element_levels[fine]
                       < required_fine_level[parent]) {
                marked_fine.push_back(fine);
            }
        }
        if (marked_fine.empty()) break;
        refine_fine_elements(marked_fine);
        if (iteration == 63) {
            throw std::runtime_error(
                "master fine and cert-audit meshes failed to expand with H");
        }
    }

    RefineOutput coarse_to_existing_fine = build_nested_mesh_embedding(
        refinement.mesh, fine_completion_.refinement.mesh);

    std::uint64_t candidate_next_id = next_element_id_;
    std::vector<std::uint64_t> new_ids(new_levels.size());
    std::vector<std::uint64_t> new_parent_ids(new_levels.size());
    for (int child = 0; child < static_cast<int>(new_levels.size()); ++child) {
        const int parent = parents[child];
        if (new_levels[child] == parent_levels[parent]) {
            new_ids[child] = parent_ids[parent];
            new_parent_ids[child] = parent_parent_ids[parent];
        } else {
            new_ids[child] = candidate_next_id++;
            new_parent_ids[child] = parent_ids[parent];
        }
    }

    coarse_mesh_ = std::move(refinement.mesh);
    coarse_levels_ = std::move(new_levels);
    coarse_element_ids_ = std::move(new_ids);
    coarse_parent_ids_ = std::move(new_parent_ids);
    next_element_id_ = candidate_next_id;
    fine_completion_.refinement = std::move(coarse_to_existing_fine);
    ++coarse_mesh_version_;
    ++interpolation_version_;
    ++boundary_version_;
    ++corrector_space_version_;
    refresh_coarse_to_cert_audit();
}

void AdaptiveMeshHierarchy::refine_fine_elements(
    const std::vector<int> &marked_fine_elements) {
    if (marked_fine_elements.empty()) return;
    const int old_fine_element_count =
        static_cast<int>(fine_completion_.refinement.mesh.elems.size());
    validate_indices(marked_fine_elements, old_fine_element_count, "fine marked element");

    const TriMesh old_fine_mesh = fine_completion_.refinement.mesh;
    const std::vector<int> old_fine_levels = fine_completion_.element_levels;
    RefineOutput fine_step = bisect_newest_vertex(old_fine_mesh, marked_fine_elements);
    std::vector<int> new_fine_levels = refinement_child_levels(
        old_fine_mesh, old_fine_levels, fine_step);
    const std::vector<int> new_fine_parents = fine_element_parents(
        fine_step.P_elem,
        static_cast<int>(fine_step.mesh.elems.size()),
        old_fine_element_count);

    std::vector<int> target_level = old_fine_levels;
    std::vector<char> changed(old_fine_element_count, false);
    for (int child = 0; child < static_cast<int>(new_fine_levels.size()); ++child) {
        const int parent = new_fine_parents[child];
        target_level[parent] = std::max(target_level[parent], new_fine_levels[child]);
        changed[parent] = changed[parent] || new_fine_levels[child] > old_fine_levels[parent];
    }

    TriMesh updated_audit_mesh = cert_audit_mesh_;
    std::vector<int> updated_audit_levels = cert_audit_element_levels_;
    RefineOutput old_fine_to_updated_audit = fine_to_cert_audit_;
    bool audit_changed = false;
    for (int iteration = 0; iteration < 64; ++iteration) {
        const std::vector<int> audit_parents = fine_element_parents(
            old_fine_to_updated_audit.P_elem,
            static_cast<int>(updated_audit_mesh.elems.size()),
            old_fine_element_count);
        std::vector<int> marked_audit;
        for (int element = 0; element < static_cast<int>(audit_parents.size()); ++element) {
            if (updated_audit_levels[element] < target_level[audit_parents[element]])
                marked_audit.push_back(element);
        }
        if (marked_audit.empty()) break;
        const TriMesh audit_parent_mesh = updated_audit_mesh;
        const std::vector<int> audit_parent_levels = updated_audit_levels;
        RefineOutput audit_step = bisect_newest_vertex(audit_parent_mesh, marked_audit);
        updated_audit_levels = refinement_child_levels(
            audit_parent_mesh, audit_parent_levels, audit_step);
        old_fine_to_updated_audit.P_node =
            audit_step.P_node * old_fine_to_updated_audit.P_node;
        old_fine_to_updated_audit.P_elem =
            audit_step.P_elem * old_fine_to_updated_audit.P_elem;
        old_fine_to_updated_audit.P_dg =
            audit_step.P_dg * old_fine_to_updated_audit.P_dg;
        updated_audit_mesh = std::move(audit_step.mesh);
        old_fine_to_updated_audit.mesh = updated_audit_mesh;
        audit_changed = true;
        if (iteration == 63)
            throw std::runtime_error("cert-audit mesh failed to catch up with local fine refinement");
    }

    RefineOutput new_fine_to_audit;
    bool nested = false;
    for (int attempt = 0; attempt < 8 && !nested; ++attempt) {
        try {
            new_fine_to_audit = build_nested_mesh_embedding(fine_step.mesh, updated_audit_mesh);
            nested = true;
        } catch (const std::invalid_argument &) {
            const std::vector<int> audit_parents = fine_element_parents(
                old_fine_to_updated_audit.P_elem,
                static_cast<int>(updated_audit_mesh.elems.size()),
                old_fine_element_count);
            std::vector<int> marked_audit;
            for (int element = 0; element < static_cast<int>(audit_parents.size()); ++element) {
                if (changed[audit_parents[element]]) marked_audit.push_back(element);
            }
            if (marked_audit.empty() || attempt == 7) throw;
            const TriMesh audit_parent_mesh = updated_audit_mesh;
            const std::vector<int> audit_parent_levels = updated_audit_levels;
            RefineOutput audit_step = bisect_newest_vertex(audit_parent_mesh, marked_audit);
            updated_audit_levels = refinement_child_levels(
                audit_parent_mesh, audit_parent_levels, audit_step);
            old_fine_to_updated_audit.P_node =
                audit_step.P_node * old_fine_to_updated_audit.P_node;
            old_fine_to_updated_audit.P_elem =
                audit_step.P_elem * old_fine_to_updated_audit.P_elem;
            old_fine_to_updated_audit.P_dg =
                audit_step.P_dg * old_fine_to_updated_audit.P_dg;
            updated_audit_mesh = std::move(audit_step.mesh);
            old_fine_to_updated_audit.mesh = updated_audit_mesh;
            audit_changed = true;
        }
    }
    if (!nested)
        throw std::runtime_error("local fine refinement is not nested in cert-audit mesh");

    fine_completion_.refinement = build_nested_mesh_embedding(coarse_mesh_, fine_step.mesh);
    fine_completion_.element_levels = std::move(new_fine_levels);
    cert_audit_mesh_ = std::move(updated_audit_mesh);
    cert_audit_element_levels_ = std::move(updated_audit_levels);
    fine_to_cert_audit_ = std::move(new_fine_to_audit);
    ++fine_mesh_version_;
    if (audit_changed) ++cert_audit_mesh_version_;
    ++interpolation_version_;
    ++boundary_version_;
    ++corrector_space_version_;
    refresh_coarse_to_cert_audit();
}

std::vector<int> AdaptiveMeshHierarchy::fine_elements_in_coarse_patch(
    const std::vector<int> &coarse_patch_elements) const {
    if (coarse_patch_elements.empty()) return {};
    validate_indices(
        coarse_patch_elements,
        static_cast<int>(coarse_mesh_.elems.size()),
        "coarse corrector patch element");
    std::vector<char> selected(coarse_mesh_.elems.size(), false);
    for (int coarse : coarse_patch_elements) selected[coarse] = true;
    std::vector<int> fine_elements;
    for (int fine = 0; fine < static_cast<int>(fine_parent_coarse_elements_.size()); ++fine) {
        if (selected[fine_parent_coarse_elements_[fine]])
            fine_elements.push_back(fine);
    }
    if (fine_elements.empty())
        throw std::runtime_error("coarse corrector patch has no fine descendants");
    return fine_elements;
}

void AdaptiveMeshHierarchy::refine_fine_in_coarse_patch(
    const std::vector<int> &coarse_patch_elements) {
    refine_fine_elements(fine_elements_in_coarse_patch(coarse_patch_elements));
}

void AdaptiveMeshHierarchy::refine_cert_audit_from_fine_elements(
    const std::vector<int> &marked_fine_elements) {
    if (marked_fine_elements.empty()) return;
    const int fine_element_count =
        static_cast<int>(fine_completion_.refinement.mesh.elems.size());
    validate_indices(marked_fine_elements, fine_element_count, "cert-audit fine element");
    std::vector<char> selected(fine_element_count, false);
    for (int element : marked_fine_elements) selected[element] = true;

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
    cert_audit_element_levels_ = refinement_child_levels(parent_mesh, parent_levels, step);
    cert_audit_mesh_ = std::move(step.mesh);
    fine_to_cert_audit_ = build_nested_mesh_embedding(
        fine_completion_.refinement.mesh, cert_audit_mesh_);
    ++cert_audit_mesh_version_;
    ++interpolation_version_;
    ++boundary_version_;
    refresh_coarse_to_cert_audit();
}

Eigen::SparseMatrix<double> AdaptiveMeshHierarchy::fine_kernel_constraints(
    const std::vector<int> &fine_dofs) const {
    return restrict_constraint_columns(fine_quasi_interpolation_, fine_dofs);
}

Eigen::SparseMatrix<double> AdaptiveMeshHierarchy::cert_audit_kernel_constraints(
    const std::vector<int> &cert_audit_dofs) const {
    return restrict_constraint_columns(cert_audit_quasi_interpolation_, cert_audit_dofs);
}

ReferenceEpochHierarchy::ReferenceEpochHierarchy(
    const TriMesh &initial_mesh,
    int initial_coarse_level,
    int reference_level)
    : initial_mesh_(initial_mesh), reference_level_(reference_level) {
    if (initial_mesh.nodes.empty() || initial_mesh.elems.empty()) {
        throw std::invalid_argument(
            "reference-epoch hierarchy initial mesh must not be empty");
    }
    if (initial_coarse_level < 0 || reference_level < initial_coarse_level) {
        throw std::invalid_argument(
            "reference-epoch levels must satisfy 0 <= H_initial <= h_ref");
    }

    RefineOutput initial_coarse =
        refine_mesh_nvb(initial_mesh_, initial_coarse_level);
    coarse_mesh_ = std::move(initial_coarse.mesh);
    coarse_levels_.assign(coarse_mesh_.elems.size(), initial_coarse_level);
    reference_completion_ = complete_to_fine_level(
        coarse_mesh_, coarse_levels_, reference_level_);
    reference_to_ambient_.mesh = reference_completion_.refinement.mesh;
    ambient_element_levels_ = reference_completion_.element_levels;
    refresh_embeddings();
}

void ReferenceEpochHierarchy::refresh_embeddings() {
    const TriMesh reference = reference_completion_.refinement.mesh;
    const TriMesh ambient = reference_to_ambient_.mesh;
    reference_completion_.refinement = build_nested_mesh_embedding(
        coarse_mesh_, reference);
    reference_to_ambient_ = build_nested_mesh_embedding(reference, ambient);
    coarse_to_ambient_ = build_nested_mesh_embedding(coarse_mesh_, ambient);

    reference_parent_coarse_elements_ = fine_element_parents(
        reference_completion_.refinement.P_elem,
        static_cast<int>(reference.elems.size()),
        static_cast<int>(coarse_mesh_.elems.size()));
    ambient_parent_coarse_elements_ = fine_element_parents(
        coarse_to_ambient_.P_elem,
        static_cast<int>(ambient.elems.size()),
        static_cast<int>(coarse_mesh_.elems.size()));
    ambient_parent_reference_elements_ = fine_element_parents(
        reference_to_ambient_.P_elem,
        static_cast<int>(ambient.elems.size()),
        static_cast<int>(reference.elems.size()));

    reference_quasi_interpolation_ = build_quasi_interp(
        coarse_mesh_, reference,
        reference_completion_.refinement.P_dg,
        cg_to_dg(reference),
        static_cast<int>(reference.nodes.size()),
        static_cast<int>(coarse_mesh_.nodes.size()));
    ambient_quasi_interpolation_ = build_quasi_interp(
        coarse_mesh_, ambient,
        coarse_to_ambient_.P_dg,
        cg_to_dg(ambient),
        static_cast<int>(ambient.nodes.size()),
        static_cast<int>(coarse_mesh_.nodes.size()));
    validate_current_embeddings();
}

void ReferenceEpochHierarchy::validate_current_embeddings() const {
    constexpr double composition_tolerance = 1e-10;
    constexpr double coordinate_tolerance = 2e-10;
    constexpr double right_inverse_tolerance = 1e-9;

    validate_boundary_tags(coarse_mesh_);
    validate_boundary_tags(reference_mesh());
    validate_boundary_tags(ambient_mesh());
    if (reference_completion_.element_levels.size()
            != reference_mesh().elems.size()
        || ambient_element_levels_.size() != ambient_mesh().elems.size()) {
        throw std::runtime_error(
            "reference-epoch element levels do not match their meshes");
    }

    const Eigen::SparseMatrix<double> composed_nodes =
        reference_to_ambient_.P_node
        * reference_completion_.refinement.P_node;
    const Eigen::SparseMatrix<double> composed_elements =
        reference_to_ambient_.P_elem
        * reference_completion_.refinement.P_elem;
    const Eigen::SparseMatrix<double> composed_dg =
        reference_to_ambient_.P_dg
        * reference_completion_.refinement.P_dg;
    if ((composed_nodes - coarse_to_ambient_.P_node).norm()
            > composition_tolerance
        || (composed_elements - coarse_to_ambient_.P_elem).norm()
            > composition_tolerance
        || (composed_dg - coarse_to_ambient_.P_dg).norm()
            > composition_tolerance) {
        throw std::runtime_error(
            "reference-epoch prolongations fail the composition check");
    }

    const auto coordinates = [](const TriMesh &mesh) {
        Eigen::MatrixXd result(mesh.nodes.size(), 2);
        for (int node = 0; node < static_cast<int>(mesh.nodes.size()); ++node)
            result.row(node) = mesh.nodes[node].transpose();
        return result;
    };
    const auto maximum_absolute_entry = [](const Eigen::MatrixXd &matrix) {
        return matrix.size() == 0 ? 0.0 : matrix.cwiseAbs().maxCoeff();
    };
    const Eigen::MatrixXd coarse_coordinates = coordinates(coarse_mesh_);
    const Eigen::MatrixXd reference_coordinates = coordinates(reference_mesh());
    const Eigen::MatrixXd ambient_coordinates = coordinates(ambient_mesh());
    if (maximum_absolute_entry(
            coarse_to_reference() * coarse_coordinates
            - reference_coordinates) > coordinate_tolerance
        || maximum_absolute_entry(
            reference_to_ambient() * reference_coordinates
            - ambient_coordinates) > coordinate_tolerance
        || maximum_absolute_entry(
            coarse_to_ambient() * coarse_coordinates
            - ambient_coordinates) > coordinate_tolerance) {
        throw std::runtime_error(
            "reference-epoch prolongations fail the coordinate check");
    }

    std::vector<char> is_dirichlet(coarse_mesh_.nodes.size(), false);
    for (int node : dirichlet_nodes(coarse_mesh_)) is_dirichlet[node] = true;
    std::vector<Eigen::Triplet<double>> expected_triplets;
    expected_triplets.reserve(coarse_mesh_.nodes.size());
    for (int node = 0; node < static_cast<int>(coarse_mesh_.nodes.size()); ++node) {
        if (!is_dirichlet[node])
            expected_triplets.emplace_back(node, node, 1.0);
    }
    Eigen::SparseMatrix<double> expected_right_inverse(
        coarse_mesh_.nodes.size(), coarse_mesh_.nodes.size());
    expected_right_inverse.setFromTriplets(
        expected_triplets.begin(), expected_triplets.end());
    const Eigen::SparseMatrix<double> reference_right_inverse =
        reference_quasi_interpolation_ * coarse_to_reference();
    const Eigen::SparseMatrix<double> ambient_right_inverse =
        ambient_quasi_interpolation_ * coarse_to_ambient();
    if ((reference_right_inverse - expected_right_inverse).norm()
            > right_inverse_tolerance
        || (ambient_right_inverse - expected_right_inverse).norm()
            > right_inverse_tolerance) {
        throw std::runtime_error(
            "reference-epoch interpolation fails the right-inverse check");
    }
}

ReferenceEpochRefinementResult
ReferenceEpochHierarchy::refine_coarse_preserving_reference(
    const std::vector<int> &marked_elements) {
    ReferenceEpochRefinementResult result;
    result.previous_element_count = coarse_mesh_.elems.size();
    result.current_element_count = coarse_mesh_.elems.size();
    if (marked_elements.empty()) {
        result.detail = "no coarse elements were marked";
        return result;
    }
    validate_indices(
        marked_elements,
        static_cast<int>(coarse_mesh_.elems.size()),
        "reference-epoch marked coarse element");

    const TriMesh parent_mesh = coarse_mesh_;
    const std::vector<int> parent_levels = coarse_levels_;
    RefineOutput candidate = bisect_newest_vertex(parent_mesh, marked_elements);
    std::vector<int> candidate_levels = refinement_child_levels(
        parent_mesh, parent_levels, candidate);
    RefineOutput candidate_to_reference;
    try {
        candidate_to_reference = build_nested_mesh_embedding(
            candidate.mesh, reference_mesh());
    } catch (const std::invalid_argument &error) {
        result.status =
            ReferenceEpochRefinementStatus::ReferenceRefreshRequired;
        result.detail = std::string(
            "candidate coarse mesh is not contained in the fixed reference mesh: ")
            + error.what();
        return result;
    }

    coarse_mesh_ = std::move(candidate.mesh);
    coarse_levels_ = std::move(candidate_levels);
    reference_completion_.refinement = std::move(candidate_to_reference);
    ++coarse_mesh_version_;
    ++interpolation_version_;
    ++boundary_version_;
    ++corrector_space_version_;
    refresh_embeddings();

    result.status = ReferenceEpochRefinementStatus::Refined;
    result.current_element_count = coarse_mesh_.elems.size();
    result.detail = "coarse mesh refined while preserving the reference epoch";
    return result;
}

bool ReferenceEpochHierarchy::reference_embedding_holds() const {
    try {
        (void)build_nested_mesh_embedding(coarse_mesh_, reference_mesh());
        return true;
    } catch (const std::invalid_argument &) {
        return false;
    }
}

double ReferenceEpochHierarchy::ambient_ratio() const {
    if (ambient_parent_coarse_elements_.size()
        != ambient_mesh().elems.size()) {
        throw std::runtime_error(
            "ambient/coarse parent map is stale or incomplete");
    }
    std::vector<double> coarse_diameters(coarse_mesh_.elems.size());
    for (int element = 0;
         element < static_cast<int>(coarse_mesh_.elems.size()); ++element) {
        coarse_diameters[element] = triangle_diameter(
            coarse_mesh_, coarse_mesh_.elems[element]);
        if (!(coarse_diameters[element] > 0.0))
            throw std::runtime_error("coarse mesh contains a degenerate element");
    }
    double result = 0.0;
    for (int element = 0;
         element < static_cast<int>(ambient_mesh().elems.size()); ++element) {
        const int parent = ambient_parent_coarse_elements_[element];
        const double diameter = triangle_diameter(
            ambient_mesh(), ambient_mesh().elems[element]);
        result = std::max(result, diameter / coarse_diameters[parent]);
    }
    return result;
}

AmbientRatioEnforcementResult
ReferenceEpochHierarchy::enforce_ambient_ratio(double rho_star) {
    if (!std::isfinite(rho_star) || !(rho_star > 0.0) || !(rho_star < 1.0)) {
        throw std::invalid_argument(
            "ambient ratio target must be finite and lie in (0,1)");
    }
    AmbientRatioEnforcementResult result;
    const std::size_t initial_elements = ambient_mesh().elems.size();
    TriMesh candidate_ambient = ambient_mesh();
    std::vector<int> candidate_levels = ambient_element_levels_;
    constexpr int maximum_refinement_steps = 128;
    const double tolerance = 1e-12 * std::max(1.0, rho_star);
    for (int iteration = 0; iteration < maximum_refinement_steps; ++iteration) {
        const RefineOutput coarse_to_candidate = build_nested_mesh_embedding(
            coarse_mesh_, candidate_ambient);
        const std::vector<int> candidate_parents = fine_element_parents(
            coarse_to_candidate.P_elem,
            static_cast<int>(candidate_ambient.elems.size()),
            static_cast<int>(coarse_mesh_.elems.size()));
        std::vector<double> coarse_diameters(coarse_mesh_.elems.size());
        for (int element = 0;
             element < static_cast<int>(coarse_mesh_.elems.size()); ++element) {
            coarse_diameters[element] = triangle_diameter(
                coarse_mesh_, coarse_mesh_.elems[element]);
        }
        double maximum_ratio = 0.0;
        std::vector<int> marked_ambient;
        for (int element = 0;
             element < static_cast<int>(candidate_ambient.elems.size()); ++element) {
            const int parent = candidate_parents[element];
            const double ratio = triangle_diameter(
                candidate_ambient, candidate_ambient.elems[element])
                / coarse_diameters[parent];
            maximum_ratio = std::max(maximum_ratio, ratio);
            if (ratio > rho_star + tolerance)
                marked_ambient.push_back(element);
        }
        if (maximum_ratio <= rho_star + tolerance) {
            result.maximum_ratio = maximum_ratio;
            result.refined_elements =
                candidate_ambient.elems.size() - initial_elements;
            if (result.changed) {
                reference_to_ambient_.mesh = std::move(candidate_ambient);
                ambient_element_levels_ = std::move(candidate_levels);
                ++ambient_mesh_version_;
                ++interpolation_version_;
                ++boundary_version_;
                refresh_embeddings();
            }
            return result;
        }
        if (marked_ambient.empty()) {
            throw std::runtime_error(
                "ambient ratio exceeds the target but no violating element was found");
        }

        const TriMesh parent_mesh = candidate_ambient;
        const std::vector<int> parent_levels = candidate_levels;
        RefineOutput step = bisect_newest_vertex(parent_mesh, marked_ambient);
        candidate_levels = refinement_child_levels(
            parent_mesh, parent_levels, step);
        candidate_ambient = std::move(step.mesh);
        result.changed = true;
        ++result.refinement_steps;
    }
    throw std::runtime_error(
        "ambient shadow failed to satisfy the mesh-ratio target");
}

void ReferenceEpochHierarchy::refresh_reference_from_ambient() {
    reference_completion_.refinement.mesh = ambient_mesh();
    reference_completion_.element_levels = ambient_element_levels_;
    reference_to_ambient_.mesh = reference_completion_.refinement.mesh;
    if (!ambient_element_levels_.empty()) {
        reference_level_ = *std::min_element(
            ambient_element_levels_.begin(), ambient_element_levels_.end());
    }
    ++reference_epoch_;
    ++reference_mesh_version_;
    ++ambient_mesh_version_;
    ++interpolation_version_;
    ++boundary_version_;
    ++corrector_space_version_;
    refresh_embeddings();
}

} // namespace lod2d::helmholtz::adaptive

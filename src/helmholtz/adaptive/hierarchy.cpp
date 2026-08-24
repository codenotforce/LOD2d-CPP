#include "helmholtz/adaptive/hierarchy.h"

#include "helmholtz/boundary.h"
#include "lod/quasi_interp.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
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
        for (BoundaryTag tag : {
                 BoundaryTag::Dirichlet,
                 BoundaryTag::Neumann,
                 BoundaryTag::Robin}) {
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

RefineOutput update_nested_mesh_embedding_after_parent_refinement(
    const TriMesh &old_parent_mesh,
    const RefineOutput &parent_refinement,
    const TriMesh &child_mesh,
    const std::vector<int> &old_child_parent_elements) {
    const TriMesh &new_parent_mesh = parent_refinement.mesh;
    if (old_parent_mesh.elems.empty() || new_parent_mesh.elems.empty()
        || child_mesh.elems.empty()) {
        throw std::invalid_argument(
            "incremental nested embedding requires nonempty meshes");
    }
    if (old_child_parent_elements.size() != child_mesh.elems.size()) {
        throw std::invalid_argument(
            "incremental nested embedding parent map has the wrong size");
    }
    const std::vector<int> new_parent_old_parents = fine_element_parents(
        parent_refinement.P_elem,
        static_cast<int>(new_parent_mesh.elems.size()),
        static_cast<int>(old_parent_mesh.elems.size()));
    std::vector<std::vector<int>> descendants(old_parent_mesh.elems.size());
    for (int element = 0;
         element < static_cast<int>(new_parent_old_parents.size()); ++element) {
        descendants[new_parent_old_parents[element]].push_back(element);
    }

    constexpr double tolerance = 2e-11;
    const int parent_nodes = static_cast<int>(new_parent_mesh.nodes.size());
    const int parent_elements = static_cast<int>(new_parent_mesh.elems.size());
    const int child_nodes = static_cast<int>(child_mesh.nodes.size());
    const int child_elements = static_cast<int>(child_mesh.elems.size());
    std::vector<int> element_parents(child_elements, -1);
    std::vector<std::array<std::array<double, 3>, 3>> element_weights(
        child_elements);
    for (int child = 0; child < child_elements; ++child) {
        const int old_parent = old_child_parent_elements[child];
        if (old_parent < 0
            || old_parent >= static_cast<int>(descendants.size())) {
            throw std::out_of_range(
                "incremental nested embedding old parent is out of range");
        }
        const Triangle &child_triangle = child_mesh.elems[child];
        const Point2 centroid = (
            child_mesh.nodes[child_triangle[0]]
            + child_mesh.nodes[child_triangle[1]]
            + child_mesh.nodes[child_triangle[2]]) / 3.0;
        for (int parent : descendants[old_parent]) {
            const Triangle &parent_triangle = new_parent_mesh.elems[parent];
            if (!inside_triangle(
                    barycentric_coordinates(
                        centroid,
                        new_parent_mesh.nodes[parent_triangle[0]],
                        new_parent_mesh.nodes[parent_triangle[1]],
                        new_parent_mesh.nodes[parent_triangle[2]]),
                    tolerance)) {
                continue;
            }
            std::array<std::array<double, 3>, 3> weights{};
            bool all_inside = true;
            for (int local = 0; local < 3; ++local) {
                weights[local] = barycentric_coordinates(
                    child_mesh.nodes[child_triangle[local]],
                    new_parent_mesh.nodes[parent_triangle[0]],
                    new_parent_mesh.nodes[parent_triangle[1]],
                    new_parent_mesh.nodes[parent_triangle[2]]);
                if (!inside_triangle(weights[local], tolerance)) {
                    all_inside = false;
                    break;
                }
                weights[local] = clean_weights(weights[local]);
            }
            if (!all_inside) continue;
            if (element_parents[child] >= 0) {
                throw std::invalid_argument(
                    "incremental nested child triangle has more than one parent");
            }
            element_parents[child] = parent;
            element_weights[child] = weights;
        }
        if (element_parents[child] < 0) {
            throw std::invalid_argument(
                "refined parent mesh is not contained in the fixed child mesh");
        }
    }

    std::vector<Eigen::Triplet<double>> element_triplets;
    std::vector<Eigen::Triplet<double>> dg_triplets;
    element_triplets.reserve(child_elements);
    dg_triplets.reserve(9 * child_elements);
    std::vector<std::map<int, double>> nodal_rows(child_nodes);
    for (int child = 0; child < child_elements; ++child) {
        const int parent = element_parents[child];
        const Triangle &parent_triangle = new_parent_mesh.elems[parent];
        const Triangle &child_triangle = child_mesh.elems[child];
        element_triplets.emplace_back(child, parent, 1.0);
        for (int child_local = 0; child_local < 3; ++child_local) {
            std::map<int, double> row;
            for (int parent_local = 0; parent_local < 3; ++parent_local) {
                const double weight =
                    element_weights[child][child_local][parent_local];
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
                    const double new_value = row.contains(column)
                        ? row.at(column) : 0.0;
                    if (std::abs(old_value - new_value) > 2e-10) {
                        throw std::invalid_argument(
                            "incremental nodal interpolation is inconsistent across an edge");
                    }
                }
            }
        }
    }

    std::vector<Eigen::Triplet<double>> node_triplets;
    node_triplets.reserve(3 * child_nodes);
    for (int child_node = 0; child_node < child_nodes; ++child_node) {
        if (nodal_rows[child_node].empty()) {
            throw std::invalid_argument(
                "incremental nested child mesh contains an unused node");
        }
        Point2 reconstructed = Point2::Zero();
        for (const auto &[parent_node, weight] : nodal_rows[child_node]) {
            node_triplets.emplace_back(child_node, parent_node, weight);
            reconstructed += weight * new_parent_mesh.nodes[parent_node];
        }
        if ((reconstructed - child_mesh.nodes[child_node]).norm() > 2e-10) {
            throw std::runtime_error(
                "incremental nodal prolongation does not reproduce coordinates");
        }
    }

    Eigen::SparseMatrix<double> node_prolongation(child_nodes, parent_nodes);
    node_prolongation.setFromTriplets(node_triplets.begin(), node_triplets.end());
    Eigen::SparseMatrix<double> element_prolongation(
        child_elements, parent_elements);
    element_prolongation.setFromTriplets(
        element_triplets.begin(), element_triplets.end());
    Eigen::SparseMatrix<double> dg_prolongation(
        3 * child_elements, 3 * parent_elements);
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
    int reference_level,
    const std::uint64_t initial_reference_epoch)
    : initial_mesh_(initial_mesh), reference_level_(reference_level),
      reference_epoch_(initial_reference_epoch) {
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
    coarse_to_ambient_ = reference_completion_.refinement;
    const int reference_nodes = static_cast<int>(
        reference_completion_.refinement.mesh.nodes.size());
    const int reference_elements = static_cast<int>(
        reference_completion_.refinement.mesh.elems.size());
    reference_to_ambient_.mesh = reference_completion_.refinement.mesh;
    reference_to_ambient_.P_node = identity_sparse(reference_nodes);
    reference_to_ambient_.P_elem = identity_sparse(reference_elements);
    reference_to_ambient_.P_dg = identity_sparse(3 * reference_elements);
    ambient_element_levels_ = reference_completion_.element_levels;
    refresh_embedding_metadata();
}

void ReferenceEpochHierarchy::refresh_embeddings() {
    const TriMesh reference = reference_completion_.refinement.mesh;
    const TriMesh ambient = reference_to_ambient_.mesh;
    reference_completion_.refinement = build_nested_mesh_embedding(
        coarse_mesh_, reference);
    reference_to_ambient_ = build_nested_mesh_embedding(reference, ambient);
    coarse_to_ambient_ = build_nested_mesh_embedding(coarse_mesh_, ambient);

    refresh_embedding_metadata();
}

void ReferenceEpochHierarchy::refresh_embedding_metadata() {
    const TriMesh &reference = reference_completion_.refinement.mesh;
    const TriMesh &ambient = reference_to_ambient_.mesh;

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
    const std::vector<int> old_reference_parents =
        reference_parent_coarse_elements_;
    const std::vector<int> old_ambient_parents =
        ambient_parent_coarse_elements_;
    RefineOutput candidate = bisect_newest_vertex(parent_mesh, marked_elements);
    std::vector<int> candidate_levels = refinement_child_levels(
        parent_mesh, parent_levels, candidate);
    RefineOutput candidate_to_reference;
    RefineOutput candidate_to_ambient;
    try {
        candidate_to_reference =
            update_nested_mesh_embedding_after_parent_refinement(
                parent_mesh, candidate, reference_mesh(),
                old_reference_parents);
        candidate_to_ambient =
            update_nested_mesh_embedding_after_parent_refinement(
                parent_mesh, candidate, ambient_mesh(), old_ambient_parents);
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
    coarse_to_ambient_ = std::move(candidate_to_ambient);
    ++coarse_mesh_version_;
    ++interpolation_version_;
    ++boundary_version_;
    ++corrector_space_version_;
    refresh_embedding_metadata();

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
    std::vector<int> candidate_parents = ambient_parent_coarse_elements_;
    RefineOutput old_ambient_to_candidate;
    old_ambient_to_candidate.mesh = candidate_ambient;
    old_ambient_to_candidate.P_node = identity_sparse(
        static_cast<int>(candidate_ambient.nodes.size()));
    old_ambient_to_candidate.P_elem = identity_sparse(
        static_cast<int>(candidate_ambient.elems.size()));
    old_ambient_to_candidate.P_dg = identity_sparse(
        3 * static_cast<int>(candidate_ambient.elems.size()));
    constexpr int maximum_refinement_steps = 128;
    const double tolerance = 1e-12 * std::max(1.0, rho_star);
    for (int iteration = 0; iteration < maximum_refinement_steps; ++iteration) {
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
                reference_to_ambient_.P_node =
                    old_ambient_to_candidate.P_node
                    * reference_to_ambient_.P_node;
                reference_to_ambient_.P_elem =
                    old_ambient_to_candidate.P_elem
                    * reference_to_ambient_.P_elem;
                reference_to_ambient_.P_dg =
                    old_ambient_to_candidate.P_dg
                    * reference_to_ambient_.P_dg;
                coarse_to_ambient_.P_node =
                    old_ambient_to_candidate.P_node
                    * coarse_to_ambient_.P_node;
                coarse_to_ambient_.P_elem =
                    old_ambient_to_candidate.P_elem
                    * coarse_to_ambient_.P_elem;
                coarse_to_ambient_.P_dg =
                    old_ambient_to_candidate.P_dg
                    * coarse_to_ambient_.P_dg;
                reference_to_ambient_.mesh = candidate_ambient;
                coarse_to_ambient_.mesh = std::move(candidate_ambient);
                ambient_element_levels_ = std::move(candidate_levels);
                ++ambient_mesh_version_;
                ++interpolation_version_;
                ++boundary_version_;
                refresh_embedding_metadata();
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
        const std::vector<int> step_parents = fine_element_parents(
            step.P_elem,
            static_cast<int>(step.mesh.elems.size()),
            static_cast<int>(parent_mesh.elems.size()));
        std::vector<int> next_candidate_parents(step_parents.size());
        for (int child = 0;
             child < static_cast<int>(step_parents.size()); ++child) {
            next_candidate_parents[child] = candidate_parents[step_parents[child]];
        }
        candidate_levels = refinement_child_levels(
            parent_mesh, parent_levels, step);
        old_ambient_to_candidate.P_node =
            step.P_node * old_ambient_to_candidate.P_node;
        old_ambient_to_candidate.P_elem =
            step.P_elem * old_ambient_to_candidate.P_elem;
        old_ambient_to_candidate.P_dg =
            step.P_dg * old_ambient_to_candidate.P_dg;
        candidate_ambient = std::move(step.mesh);
        old_ambient_to_candidate.mesh = candidate_ambient;
        candidate_parents = std::move(next_candidate_parents);
        result.changed = true;
        ++result.refinement_steps;
    }
    throw std::runtime_error(
        "ambient shadow failed to satisfy the mesh-ratio target");
}

void ReferenceEpochHierarchy::refresh_reference_from_ambient() {
    reference_completion_.refinement = coarse_to_ambient_;
    reference_completion_.element_levels = ambient_element_levels_;
    const int ambient_nodes = static_cast<int>(ambient_mesh().nodes.size());
    const int ambient_elements = static_cast<int>(ambient_mesh().elems.size());
    reference_to_ambient_.mesh = ambient_mesh();
    reference_to_ambient_.P_node = identity_sparse(ambient_nodes);
    reference_to_ambient_.P_elem = identity_sparse(ambient_elements);
    reference_to_ambient_.P_dg = identity_sparse(3 * ambient_elements);
    reference_parent_coarse_elements_ = ambient_parent_coarse_elements_;
    ambient_parent_reference_elements_.resize(ambient_elements);
    std::iota(
        ambient_parent_reference_elements_.begin(),
        ambient_parent_reference_elements_.end(), 0);
    reference_quasi_interpolation_ = ambient_quasi_interpolation_;
    if (proposed_coarse_) {
        proposed_coarse_->to_reference = proposed_coarse_->to_candidate;
    }
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
    // The promoted reference is exactly the existing candidate.  All maps
    // above are therefore known algebraically; rebuilding them by geometric
    // search would duplicate the most expensive part of an epoch switch.
    validate_incremental_candidate_metadata();
}

void ReferenceEpochHierarchy::validate_incremental_candidate_metadata() const {
    validate_boundary_tags(candidate_mesh());
    const int coarse_nodes = static_cast<int>(coarse_mesh_.nodes.size());
    const int reference_nodes = static_cast<int>(reference_mesh().nodes.size());
    const int candidate_nodes = static_cast<int>(candidate_mesh().nodes.size());
    const int coarse_elements = static_cast<int>(coarse_mesh_.elems.size());
    const int reference_elements = static_cast<int>(reference_mesh().elems.size());
    const int candidate_elements = static_cast<int>(candidate_mesh().elems.size());
    if (reference_to_ambient_.P_node.rows() != candidate_nodes
        || reference_to_ambient_.P_node.cols() != reference_nodes
        || coarse_to_ambient_.P_node.rows() != candidate_nodes
        || coarse_to_ambient_.P_node.cols() != coarse_nodes
        || reference_to_ambient_.P_elem.rows() != candidate_elements
        || reference_to_ambient_.P_elem.cols() != reference_elements
        || coarse_to_ambient_.P_elem.rows() != candidate_elements
        || coarse_to_ambient_.P_elem.cols() != coarse_elements
        || reference_to_ambient_.P_dg.rows() != 3 * candidate_elements
        || reference_to_ambient_.P_dg.cols() != 3 * reference_elements
        || coarse_to_ambient_.P_dg.rows() != 3 * candidate_elements
        || coarse_to_ambient_.P_dg.cols() != 3 * coarse_elements
        || reference_completion_.refinement.P_node.rows() != reference_nodes
        || reference_completion_.refinement.P_node.cols() != coarse_nodes
        || reference_completion_.refinement.P_elem.rows() != reference_elements
        || reference_completion_.refinement.P_elem.cols() != coarse_elements
        || reference_completion_.refinement.P_dg.rows() != 3 * reference_elements
        || reference_completion_.refinement.P_dg.cols() != 3 * coarse_elements
        || reference_quasi_interpolation_.rows() != coarse_nodes
        || reference_quasi_interpolation_.cols() != reference_nodes
        || ambient_quasi_interpolation_.rows() != coarse_nodes
        || ambient_quasi_interpolation_.cols() != candidate_nodes
        || ambient_element_levels_.size() != candidate_mesh().elems.size()
        || ambient_parent_coarse_elements_.size()
            != candidate_mesh().elems.size()
        || ambient_parent_reference_elements_.size()
            != candidate_mesh().elems.size()) {
        throw std::runtime_error(
            "incremental candidate embedding metadata has inconsistent dimensions");
    }
    for (int parent : ambient_parent_coarse_elements_) {
        if (parent < 0 || parent >= coarse_elements)
            throw std::runtime_error(
                "incremental candidate/coarse parent is out of range");
    }
    for (int parent : ambient_parent_reference_elements_) {
        if (parent < 0 || parent >= reference_elements)
            throw std::runtime_error(
                "incremental candidate/reference parent is out of range");
    }
}

void ReferenceEpochHierarchy::replace_candidate(
    TriMesh mesh,
    std::vector<int> element_levels) {
    if (element_levels.size() != mesh.elems.size()) {
        throw std::invalid_argument(
            "candidate element levels do not match the candidate mesh");
    }
    reference_to_ambient_ = build_nested_mesh_embedding(
        reference_mesh(), mesh);
    coarse_to_ambient_ = build_nested_mesh_embedding(coarse_mesh_, mesh);
    reference_to_ambient_.mesh = mesh;
    coarse_to_ambient_.mesh = std::move(mesh);
    ambient_element_levels_ = std::move(element_levels);
    ++ambient_mesh_version_;
    ++interpolation_version_;
    ++boundary_version_;
    refresh_embedding_metadata();
}

void ReferenceEpochHierarchy::begin_reference_epoch() {
    proposed_coarse_.reset();
    bool same_mesh = reference_mesh().elems == candidate_mesh().elems
        && reference_mesh().nodes.size() == candidate_mesh().nodes.size();
    for (std::size_t node = 0;
         same_mesh && node < reference_mesh().nodes.size(); ++node) {
        same_mesh = reference_mesh().nodes[node].isApprox(
            candidate_mesh().nodes[node], 0.0);
    }
    if (same_mesh) return;
    const TriMesh reference = reference_mesh();
    replace_candidate(reference, reference_completion_.element_levels);
}

ReferenceEpochCoarseRefinementPreview
ReferenceEpochHierarchy::preview_coarse_refinement(
    const std::vector<int> &marked_elements,
    const ReferenceEpochRefinementGuard &resource_guard) const {
    if (proposed_coarse_) {
        throw std::logic_error(
            "a coarse refinement preview requires no pending proposal");
    }
    if (marked_elements.empty()) {
        throw std::invalid_argument(
            "a coarse refinement preview requires marked elements");
    }
    validate_indices(
        marked_elements,
        static_cast<int>(coarse_mesh_.elems.size()),
        "reference-epoch preview marked coarse element");

    ReferenceEpochCoarseRefinementPreview preview;
    preview.marked_elements = marked_elements;
    preview.coarse_mesh_version = coarse_mesh_version_;
    preview.reference_mesh_version = reference_mesh_version_;
    if (resource_guard) {
        resource_guard(
            ReferenceEpochRefinementGuardPoint::BeforeNvb, coarse_mesh_);
    }
    const auto nvb_begin = std::chrono::steady_clock::now();
    preview.refinement = bisect_newest_vertex(coarse_mesh_, marked_elements);
    if (resource_guard) {
        resource_guard(
            ReferenceEpochRefinementGuardPoint::AfterNvb,
            preview.refinement.mesh);
    }
    preview.element_levels = refinement_child_levels(
        coarse_mesh_, coarse_levels_, preview.refinement);
    preview.time_nvb_refine = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - nvb_begin).count();

    const auto embedding_begin = std::chrono::steady_clock::now();
    try {
        RefineOutput embedding =
            update_nested_mesh_embedding_after_parent_refinement(
                coarse_mesh_, preview.refinement, reference_mesh(),
                reference_parent_coarse_elements_);
        ReferenceEpochCoarseRefinementPreview::CachedEmbedding cached;
        cached.P_node = std::move(embedding.P_node);
        cached.P_elem = std::move(embedding.P_elem);
        cached.P_dg = std::move(embedding.P_dg);
        preview.to_reference = std::move(cached);
    } catch (const std::invalid_argument &) {
        preview.to_reference.reset();
    }
    preview.time_reference_embedding_update =
        std::chrono::duration<double>(
            std::chrono::steady_clock::now() - embedding_begin).count();
    return preview;
}

ReferenceEpochRefinementResult
ReferenceEpochHierarchy::propose_coarse_refinement(
    const std::vector<int> &marked_elements,
    const ReferenceEpochRefinementGuard &resource_guard) {
    ReferenceEpochRefinementResult result;
    result.previous_element_count = coarse_mesh_.elems.size();
    result.current_element_count = coarse_mesh_.elems.size();
    if (proposed_coarse_) {
        throw std::logic_error(
            "a coarse refinement proposal is already pending");
    }
    if (marked_elements.empty()) {
        result.detail = "no coarse elements were marked";
        return result;
    }
    validate_indices(
        marked_elements,
        static_cast<int>(coarse_mesh_.elems.size()),
        "reference-epoch marked coarse element");

    ProposedCoarseState proposed;
    if (resource_guard) {
        resource_guard(
            ReferenceEpochRefinementGuardPoint::BeforeNvb, coarse_mesh_);
    }
    const auto nvb_begin = std::chrono::steady_clock::now();
    proposed.refinement = bisect_newest_vertex(coarse_mesh_, marked_elements);
    if (resource_guard) {
        resource_guard(
            ReferenceEpochRefinementGuardPoint::AfterNvb,
            proposed.refinement.mesh);
    }
    proposed.levels = refinement_child_levels(
        coarse_mesh_, coarse_levels_, proposed.refinement);
    result.time_nvb_refine = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - nvb_begin).count();
    const auto cache = [](RefineOutput embedding) {
        ProposedCoarseState::CachedEmbedding result;
        result.P_node = std::move(embedding.P_node);
        result.P_elem = std::move(embedding.P_elem);
        result.P_dg = std::move(embedding.P_dg);
        return result;
    };
    const auto embedding_begin = std::chrono::steady_clock::now();
    try {
        proposed.to_reference = cache(
            update_nested_mesh_embedding_after_parent_refinement(
                coarse_mesh_, proposed.refinement, reference_mesh(),
                reference_parent_coarse_elements_));
    } catch (const std::invalid_argument &) {
        proposed.to_reference.reset();
    }
    try {
        proposed.to_candidate = cache(
            update_nested_mesh_embedding_after_parent_refinement(
                coarse_mesh_, proposed.refinement, candidate_mesh(),
                ambient_parent_coarse_elements_));
    } catch (const std::invalid_argument &) {
        proposed.to_candidate.reset();
    }
    result.time_proposed_embedding_update = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - embedding_begin).count();
    result.status = ReferenceEpochRefinementStatus::Refined;
    result.current_element_count = proposed.refinement.mesh.elems.size();
    result.detail = "coarse refinement proposed without mutating the committed mesh";
    proposed_coarse_ = std::move(proposed);
    return result;
}

ReferenceEpochRefinementResult
ReferenceEpochHierarchy::propose_coarse_refinement(
    ReferenceEpochCoarseRefinementPreview preview,
    const ReferenceEpochRefinementGuard &resource_guard) {
    ReferenceEpochRefinementResult result;
    result.previous_element_count = coarse_mesh_.elems.size();
    result.current_element_count = coarse_mesh_.elems.size();
    if (proposed_coarse_) {
        throw std::logic_error(
            "a coarse refinement proposal is already pending");
    }
    if (preview.marked_elements.empty()) {
        throw std::invalid_argument(
            "an adopted coarse refinement preview has no marked elements");
    }
    if (preview.coarse_mesh_version != coarse_mesh_version_
        || preview.reference_mesh_version != reference_mesh_version_) {
        throw std::logic_error(
            "coarse refinement preview is stale for the current hierarchy");
    }
    validate_indices(
        preview.marked_elements,
        static_cast<int>(coarse_mesh_.elems.size()),
        "adopted preview marked coarse element");
    const int coarse_nodes = static_cast<int>(coarse_mesh_.nodes.size());
    const int coarse_elements = static_cast<int>(coarse_mesh_.elems.size());
    const int proposed_nodes =
        static_cast<int>(preview.refinement.mesh.nodes.size());
    const int proposed_elements =
        static_cast<int>(preview.refinement.mesh.elems.size());
    if (preview.refinement.P_node.rows() != proposed_nodes
        || preview.refinement.P_node.cols() != coarse_nodes
        || preview.refinement.P_elem.rows() != proposed_elements
        || preview.refinement.P_elem.cols() != coarse_elements
        || preview.refinement.P_dg.rows() != 3 * proposed_elements
        || preview.refinement.P_dg.cols() != 3 * coarse_elements
        || preview.element_levels.size()
            != preview.refinement.mesh.elems.size()) {
        throw std::invalid_argument(
            "adopted coarse refinement preview has inconsistent dimensions");
    }
    if (preview.to_reference
        && (preview.to_reference->P_node.rows()
                != static_cast<int>(reference_mesh().nodes.size())
            || preview.to_reference->P_node.cols() != proposed_nodes
            || preview.to_reference->P_elem.rows()
                != static_cast<int>(reference_mesh().elems.size())
            || preview.to_reference->P_elem.cols() != proposed_elements
            || preview.to_reference->P_dg.rows()
                != 3 * static_cast<int>(reference_mesh().elems.size())
            || preview.to_reference->P_dg.cols() != 3 * proposed_elements)) {
        throw std::invalid_argument(
            "adopted coarse refinement preview has an invalid reference embedding");
    }
    if (resource_guard) {
        resource_guard(
            ReferenceEpochRefinementGuardPoint::BeforeNvb, coarse_mesh_);
        resource_guard(
            ReferenceEpochRefinementGuardPoint::AfterNvb,
            preview.refinement.mesh);
    }

    ProposedCoarseState proposed;
    proposed.refinement = std::move(preview.refinement);
    proposed.levels = std::move(preview.element_levels);
    if (preview.to_reference) {
        ProposedCoarseState::CachedEmbedding cached;
        cached.P_node = std::move(preview.to_reference->P_node);
        cached.P_elem = std::move(preview.to_reference->P_elem);
        cached.P_dg = std::move(preview.to_reference->P_dg);
        proposed.to_reference = std::move(cached);
    }
    const auto embedding_begin = std::chrono::steady_clock::now();
    try {
        RefineOutput embedding =
            update_nested_mesh_embedding_after_parent_refinement(
                coarse_mesh_, proposed.refinement, candidate_mesh(),
                ambient_parent_coarse_elements_);
        ProposedCoarseState::CachedEmbedding cached;
        cached.P_node = std::move(embedding.P_node);
        cached.P_elem = std::move(embedding.P_elem);
        cached.P_dg = std::move(embedding.P_dg);
        proposed.to_candidate = std::move(cached);
    } catch (const std::invalid_argument &) {
        proposed.to_candidate.reset();
    }
    result.time_proposed_embedding_update = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - embedding_begin).count();
    result.status = ReferenceEpochRefinementStatus::Refined;
    result.current_element_count = proposed.refinement.mesh.elems.size();
    result.detail =
        "cached coarse preview adopted without repeating NVB/reference embedding";
    proposed_coarse_ = std::move(proposed);
    return result;
}

ReferenceEpochRefinementResult
ReferenceEpochHierarchy::propose_identity_coarse() {
    ReferenceEpochRefinementResult result;
    result.previous_element_count = coarse_mesh_.elems.size();
    result.current_element_count = coarse_mesh_.elems.size();
    if (proposed_coarse_) {
        throw std::logic_error(
            "an identity coarse proposal cannot replace a pending proposal");
    }
    ProposedCoarseState proposed;
    proposed.refinement.mesh = coarse_mesh_;
    proposed.refinement.P_node = identity_sparse(
        static_cast<int>(coarse_mesh_.nodes.size()));
    proposed.refinement.P_elem = identity_sparse(
        static_cast<int>(coarse_mesh_.elems.size()));
    proposed.refinement.P_dg = identity_sparse(
        3 * static_cast<int>(coarse_mesh_.elems.size()));
    proposed.levels = coarse_levels_;
    ProposedCoarseState::CachedEmbedding to_reference;
    to_reference.P_node = reference_completion_.refinement.P_node;
    to_reference.P_elem = reference_completion_.refinement.P_elem;
    to_reference.P_dg = reference_completion_.refinement.P_dg;
    proposed.to_reference = std::move(to_reference);
    ProposedCoarseState::CachedEmbedding to_candidate;
    to_candidate.P_node = coarse_to_ambient_.P_node;
    to_candidate.P_elem = coarse_to_ambient_.P_elem;
    to_candidate.P_dg = coarse_to_ambient_.P_dg;
    proposed.to_candidate = std::move(to_candidate);
    proposed.identity = true;
    proposed_coarse_ = std::move(proposed);
    result.status = ReferenceEpochRefinementStatus::Refined;
    result.detail = "identity coarse proposal opened";
    return result;
}

ReferenceEpochRefinementResult
ReferenceEpochHierarchy::refine_proposed_coarse(
    const std::vector<int> &marked_proposed_elements,
    const ReferenceEpochRefinementGuard &resource_guard) {
    ReferenceEpochRefinementResult result;
    if (!proposed_coarse_) {
        throw std::logic_error(
            "refining a proposed coarse mesh requires a pending proposal");
    }
    const ProposedCoarseState &old = *proposed_coarse_;
    result.previous_element_count = old.refinement.mesh.elems.size();
    result.current_element_count = old.refinement.mesh.elems.size();
    if (marked_proposed_elements.empty()) {
        result.detail = "no proposed coarse elements were marked";
        return result;
    }
    validate_indices(
        marked_proposed_elements,
        static_cast<int>(old.refinement.mesh.elems.size()),
        "marked proposed coarse element");
    if (resource_guard) {
        resource_guard(
            ReferenceEpochRefinementGuardPoint::BeforeNvb,
            old.refinement.mesh);
    }

    // Cache the fine-to-old-proposal parent maps before refining the parent.
    // If an embedding was already known to be absent, keep it absent: a
    // transactional matching closure must not re-enter the global geometry
    // search merely to probe whether containment has changed.
    std::vector<int> old_reference_parents;
    std::vector<int> old_candidate_parents;
    const auto parent_begin = std::chrono::steady_clock::now();
    if (old.to_reference) {
        old_reference_parents = fine_element_parents(
            old.to_reference->P_elem,
            static_cast<int>(reference_mesh().elems.size()),
            static_cast<int>(old.refinement.mesh.elems.size()));
    }
    if (old.to_candidate) {
        old_candidate_parents = fine_element_parents(
            old.to_candidate->P_elem,
            static_cast<int>(candidate_mesh().elems.size()),
            static_cast<int>(old.refinement.mesh.elems.size()));
    }
    result.time_parent_map_update = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - parent_begin).count();

    const auto nvb_begin = std::chrono::steady_clock::now();
    RefineOutput step = bisect_newest_vertex(
        old.refinement.mesh, marked_proposed_elements);
    if (resource_guard) {
        resource_guard(
            ReferenceEpochRefinementGuardPoint::AfterNvb, step.mesh);
    }
    std::vector<int> next_levels = refinement_child_levels(
        old.refinement.mesh, old.levels, step);
    result.time_nvb_refine = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - nvb_begin).count();

    ProposedCoarseState next;
    next.identity = false;
    const auto composition_begin = std::chrono::steady_clock::now();
    next.refinement.P_node = step.P_node * old.refinement.P_node;
    next.refinement.P_elem = step.P_elem * old.refinement.P_elem;
    next.refinement.P_dg = step.P_dg * old.refinement.P_dg;
    next.refinement.P_node.makeCompressed();
    next.refinement.P_elem.makeCompressed();
    next.refinement.P_dg.makeCompressed();
    result.time_embedding_composition = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - composition_begin).count();

    const auto embedding_begin = std::chrono::steady_clock::now();
    const auto cache = [](RefineOutput embedding) {
        ProposedCoarseState::CachedEmbedding cached;
        cached.P_node = std::move(embedding.P_node);
        cached.P_elem = std::move(embedding.P_elem);
        cached.P_dg = std::move(embedding.P_dg);
        return cached;
    };
    if (old.to_reference) {
        try {
            next.to_reference = cache(
                update_nested_mesh_embedding_after_parent_refinement(
                    old.refinement.mesh, step, reference_mesh(),
                    old_reference_parents));
        } catch (const std::invalid_argument &) {
            next.to_reference.reset();
        }
    }
    if (old.to_candidate) {
        try {
            next.to_candidate = cache(
                update_nested_mesh_embedding_after_parent_refinement(
                    old.refinement.mesh, step, candidate_mesh(),
                    old_candidate_parents));
        } catch (const std::invalid_argument &) {
            next.to_candidate.reset();
        }
    }
    result.time_proposed_embedding_update = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - embedding_begin).count();

    next.refinement.mesh = std::move(step.mesh);
    next.levels = std::move(next_levels);
    result.status = ReferenceEpochRefinementStatus::Refined;
    result.current_element_count = next.refinement.mesh.elems.size();
    result.detail =
        "pending coarse proposal refined without mutating the committed mesh";
    proposed_coarse_ = std::move(next);
    return result;
}

ReferenceEpochRefinementResult ReferenceEpochHierarchy::enrich_candidate(
    const std::vector<int> &marked_candidate_elements,
    const ReferenceEpochRefinementGuard &resource_guard) {
    ReferenceEpochRefinementResult result;
    result.previous_element_count = candidate_mesh().elems.size();
    result.current_element_count = candidate_mesh().elems.size();
    if (marked_candidate_elements.empty()) {
        result.detail = "no candidate elements were marked";
        return result;
    }
    validate_indices(
        marked_candidate_elements,
        static_cast<int>(candidate_mesh().elems.size()),
        "marked candidate element");
    if (resource_guard) {
        resource_guard(
            ReferenceEpochRefinementGuardPoint::BeforeNvb,
            candidate_mesh());
    }
    const TriMesh parent = candidate_mesh();
    const std::vector<int> parent_levels = ambient_element_levels_;
    const auto nvb_begin = std::chrono::steady_clock::now();
    RefineOutput refined = bisect_newest_vertex(
        parent, marked_candidate_elements);
    if (resource_guard) {
        resource_guard(
            ReferenceEpochRefinementGuardPoint::AfterNvb, refined.mesh);
    }
    std::vector<int> levels = refinement_child_levels(
        parent, parent_levels, refined);
    result.time_nvb_refine = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - nvb_begin).count();

    const auto parent_begin = std::chrono::steady_clock::now();
    const std::vector<int> step_parents = fine_element_parents(
        refined.P_elem,
        static_cast<int>(refined.mesh.elems.size()),
        static_cast<int>(parent.elems.size()));
    std::vector<int> next_coarse_parents(step_parents.size());
    std::vector<int> next_reference_parents(step_parents.size());
    for (int child = 0; child < static_cast<int>(step_parents.size()); ++child) {
        const int old_parent = step_parents[child];
        next_coarse_parents[child] =
            ambient_parent_coarse_elements_[old_parent];
        next_reference_parents[child] =
            ambient_parent_reference_elements_[old_parent];
    }
    result.time_parent_map_update = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - parent_begin).count();

    const auto composition_begin = std::chrono::steady_clock::now();
    Eigen::SparseMatrix<double> next_reference_nodes =
        refined.P_node * reference_to_ambient_.P_node;
    Eigen::SparseMatrix<double> next_reference_elements =
        refined.P_elem * reference_to_ambient_.P_elem;
    Eigen::SparseMatrix<double> next_reference_dg =
        refined.P_dg * reference_to_ambient_.P_dg;
    Eigen::SparseMatrix<double> next_coarse_nodes =
        refined.P_node * coarse_to_ambient_.P_node;
    Eigen::SparseMatrix<double> next_coarse_elements =
        refined.P_elem * coarse_to_ambient_.P_elem;
    Eigen::SparseMatrix<double> next_coarse_dg =
        refined.P_dg * coarse_to_ambient_.P_dg;
    std::optional<ProposedCoarseState::CachedEmbedding>
        next_proposed_to_candidate;
    if (proposed_coarse_ && proposed_coarse_->to_candidate) {
        ProposedCoarseState::CachedEmbedding embedding;
        embedding.P_node = refined.P_node
            * proposed_coarse_->to_candidate->P_node;
        embedding.P_elem = refined.P_elem
            * proposed_coarse_->to_candidate->P_elem;
        embedding.P_dg = refined.P_dg
            * proposed_coarse_->to_candidate->P_dg;
        embedding.P_node.makeCompressed();
        embedding.P_elem.makeCompressed();
        embedding.P_dg.makeCompressed();
        next_proposed_to_candidate = std::move(embedding);
    }
    next_reference_nodes.makeCompressed();
    next_reference_elements.makeCompressed();
    next_reference_dg.makeCompressed();
    next_coarse_nodes.makeCompressed();
    next_coarse_elements.makeCompressed();
    next_coarse_dg.makeCompressed();
    result.time_embedding_composition = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - composition_begin).count();

    TriMesh next_candidate = std::move(refined.mesh);
    reference_to_ambient_.P_node = std::move(next_reference_nodes);
    reference_to_ambient_.P_elem = std::move(next_reference_elements);
    reference_to_ambient_.P_dg = std::move(next_reference_dg);
    coarse_to_ambient_.P_node = std::move(next_coarse_nodes);
    coarse_to_ambient_.P_elem = std::move(next_coarse_elements);
    coarse_to_ambient_.P_dg = std::move(next_coarse_dg);
    reference_to_ambient_.mesh = next_candidate;
    coarse_to_ambient_.mesh = std::move(next_candidate);
    ambient_element_levels_ = std::move(levels);
    ambient_parent_coarse_elements_ = std::move(next_coarse_parents);
    ambient_parent_reference_elements_ = std::move(next_reference_parents);
    if (proposed_coarse_ && next_proposed_to_candidate) {
        proposed_coarse_->to_candidate =
            std::move(next_proposed_to_candidate);
    }

    const auto interpolation_begin = std::chrono::steady_clock::now();
    ambient_quasi_interpolation_ = build_quasi_interp(
        coarse_mesh_, candidate_mesh(), coarse_to_ambient_.P_dg,
        cg_to_dg(candidate_mesh()),
        static_cast<int>(candidate_mesh().nodes.size()),
        static_cast<int>(coarse_mesh_.nodes.size()));
    if (proposed_coarse_ && !proposed_coarse_->to_candidate) {
        try {
            RefineOutput embedding = build_nested_mesh_embedding(
                proposed_coarse_->refinement.mesh, candidate_mesh());
            ProposedCoarseState::CachedEmbedding cached;
            cached.P_node = std::move(embedding.P_node);
            cached.P_elem = std::move(embedding.P_elem);
            cached.P_dg = std::move(embedding.P_dg);
            proposed_coarse_->to_candidate = std::move(cached);
        } catch (const std::invalid_argument &) {
            // Candidate closure may need another NVB step.  The negative
            // result is cached until that step changes the candidate.
        }
    }
    result.time_candidate_quasi_interpolation =
        std::chrono::duration<double>(
            std::chrono::steady_clock::now() - interpolation_begin).count();

    ++ambient_mesh_version_;
    ++interpolation_version_;
    ++boundary_version_;
    const auto validation_begin = std::chrono::steady_clock::now();
    validate_incremental_candidate_metadata();
    result.time_embedding_validation = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - validation_begin).count();
    result.status = ReferenceEpochRefinementStatus::Refined;
    result.current_element_count = candidate_mesh().elems.size();
    result.detail = "candidate mesh enriched while the reference remained fixed";
    return result;
}

bool ReferenceEpochHierarchy::candidate_contains_proposed_coarse() const {
    if (!proposed_coarse_) return true;
    return proposed_coarse_->to_candidate.has_value();
}

bool ReferenceEpochHierarchy::reference_contains_proposed_coarse() const {
    if (!proposed_coarse_) return true;
    return proposed_coarse_->to_reference.has_value();
}

int ReferenceEpochHierarchy::minimum_proposed_reference_level_gap() const {
    if (!proposed_coarse_ || !proposed_coarse_->to_reference)
        return std::numeric_limits<int>::max();
    const std::vector<int> parents = fine_element_parents(
        proposed_coarse_->to_reference->P_elem,
        static_cast<int>(reference_mesh().elems.size()),
        static_cast<int>(proposed_coarse_->refinement.mesh.elems.size()));
    int gap = std::numeric_limits<int>::max();
    for (int child = 0; child < static_cast<int>(parents.size()); ++child) {
        gap = std::min(
            gap,
            reference_completion_.element_levels[child]
                - proposed_coarse_->levels[parents[child]]);
    }
    return gap;
}

int ReferenceEpochHierarchy::minimum_proposed_reference_level_gap(
    const std::vector<char> &included_proposed_elements) const {
    if (!proposed_coarse_ || !proposed_coarse_->to_reference)
        return std::numeric_limits<int>::max();
    if (included_proposed_elements.size()
        != proposed_coarse_->refinement.mesh.elems.size()) {
        throw std::invalid_argument(
            "proposed reference gap mask does not match the proposed mesh");
    }
    const std::vector<int> parents = fine_element_parents(
        proposed_coarse_->to_reference->P_elem,
        static_cast<int>(reference_mesh().elems.size()),
        static_cast<int>(proposed_coarse_->refinement.mesh.elems.size()));
    int gap = std::numeric_limits<int>::max();
    for (int child = 0; child < static_cast<int>(parents.size()); ++child) {
        if (!included_proposed_elements[parents[child]]) continue;
        const int local_gap = reference_completion_.element_levels[child]
            - proposed_coarse_->levels[parents[child]];
        gap = std::min(gap, local_gap);
    }
    return gap;
}

int ReferenceEpochHierarchy::minimum_proposed_candidate_level_gap() const {
    std::vector<char> included;
    if (proposed_coarse_) {
        included.assign(proposed_coarse_->refinement.mesh.elems.size(), true);
    }
    return minimum_proposed_candidate_level_gap(included);
}

int ReferenceEpochHierarchy::minimum_proposed_candidate_level_gap(
    const std::vector<char> &included_proposed_elements) const {
    if (!proposed_coarse_ || !proposed_coarse_->to_candidate)
        return std::numeric_limits<int>::max();
    if (included_proposed_elements.size()
        != proposed_coarse_->refinement.mesh.elems.size()) {
        throw std::invalid_argument(
            "proposed candidate gap mask does not match the proposed mesh");
    }
    const std::vector<int> parents = fine_element_parents(
        proposed_coarse_->to_candidate->P_elem,
        static_cast<int>(candidate_mesh().elems.size()),
        static_cast<int>(proposed_coarse_->refinement.mesh.elems.size()));
    int gap = std::numeric_limits<int>::max();
    for (int child = 0; child < static_cast<int>(parents.size()); ++child) {
        if (!included_proposed_elements[parents[child]]) continue;
        gap = std::min(
            gap,
            ambient_element_levels_[child]
                - proposed_coarse_->levels[parents[child]]);
    }
    return gap;
}

int ReferenceEpochHierarchy::minimum_proposed_reference_level_gap_margin(
    const std::vector<int> &target_level_gaps) const {
    if (!proposed_coarse_ || !proposed_coarse_->to_reference)
        return std::numeric_limits<int>::max();
    if (target_level_gaps.size()
        != proposed_coarse_->refinement.mesh.elems.size()) {
        throw std::invalid_argument(
            "proposed reference gap targets do not match the proposed mesh");
    }
    if (std::any_of(
            target_level_gaps.begin(), target_level_gaps.end(),
            [](const int target) { return target < 0; })) {
        throw std::invalid_argument(
            "proposed reference gap targets must be nonnegative");
    }
    const std::vector<int> parents = fine_element_parents(
        proposed_coarse_->to_reference->P_elem,
        static_cast<int>(reference_mesh().elems.size()),
        static_cast<int>(proposed_coarse_->refinement.mesh.elems.size()));
    int margin = std::numeric_limits<int>::max();
    for (int child = 0; child < static_cast<int>(parents.size()); ++child) {
        const int parent = parents[child];
        margin = std::min(
            margin,
            reference_completion_.element_levels[child]
                - proposed_coarse_->levels[parent]
                - target_level_gaps[parent]);
    }
    return margin;
}

int ReferenceEpochHierarchy::minimum_proposed_candidate_level_gap_margin(
    const std::vector<int> &target_level_gaps) const {
    if (!proposed_coarse_ || !proposed_coarse_->to_candidate)
        return std::numeric_limits<int>::max();
    if (target_level_gaps.size()
        != proposed_coarse_->refinement.mesh.elems.size()) {
        throw std::invalid_argument(
            "proposed candidate gap targets do not match the proposed mesh");
    }
    if (std::any_of(
            target_level_gaps.begin(), target_level_gaps.end(),
            [](const int target) { return target < 0; })) {
        throw std::invalid_argument(
            "proposed candidate gap targets must be nonnegative");
    }
    const std::vector<int> parents = fine_element_parents(
        proposed_coarse_->to_candidate->P_elem,
        static_cast<int>(candidate_mesh().elems.size()),
        static_cast<int>(proposed_coarse_->refinement.mesh.elems.size()));
    int margin = std::numeric_limits<int>::max();
    for (int child = 0; child < static_cast<int>(parents.size()); ++child) {
        const int parent = parents[child];
        margin = std::min(
            margin,
            ambient_element_levels_[child]
                - proposed_coarse_->levels[parent]
                - target_level_gaps[parent]);
    }
    return margin;
}

const TriMesh &ReferenceEpochHierarchy::proposed_coarse_mesh() const {
    if (!proposed_coarse_) {
        throw std::logic_error("there is no pending coarse refinement proposal");
    }
    return proposed_coarse_->refinement.mesh;
}

const std::vector<int> &
ReferenceEpochHierarchy::proposed_coarse_levels() const {
    if (!proposed_coarse_) {
        throw std::logic_error("there is no pending coarse refinement proposal");
    }
    return proposed_coarse_->levels;
}

const Eigen::SparseMatrix<double> &
ReferenceEpochHierarchy::coarse_elements_to_proposed_coarse() const {
    if (!proposed_coarse_) {
        throw std::logic_error("there is no pending coarse refinement proposal");
    }
    return proposed_coarse_->refinement.P_elem;
}

const Eigen::SparseMatrix<double> &
ReferenceEpochHierarchy::coarse_nodes_to_proposed_coarse() const {
    if (!proposed_coarse_) {
        throw std::logic_error("there is no pending coarse refinement proposal");
    }
    return proposed_coarse_->refinement.P_node;
}

const Eigen::SparseMatrix<double> &
ReferenceEpochHierarchy::proposed_coarse_elements_to_reference() const {
    if (!proposed_coarse_) {
        throw std::logic_error("there is no pending coarse refinement proposal");
    }
    if (!proposed_coarse_->to_reference) {
        throw std::logic_error(
            "the pending coarse proposal is not contained in the reference");
    }
    return proposed_coarse_->to_reference->P_elem;
}

const Eigen::SparseMatrix<double> &
ReferenceEpochHierarchy::proposed_coarse_elements_to_candidate() const {
    if (!proposed_coarse_) {
        throw std::logic_error("there is no pending coarse refinement proposal");
    }
    if (!proposed_coarse_->to_candidate) {
        throw std::logic_error(
            "the pending coarse proposal is not contained in the candidate");
    }
    return proposed_coarse_->to_candidate->P_elem;
}

std::vector<int>
ReferenceEpochHierarchy::reference_parent_proposed_coarse_elements() const {
    const Eigen::SparseMatrix<double> &embedding =
        proposed_coarse_elements_to_reference();
    return fine_element_parents(
        embedding,
        static_cast<int>(reference_mesh().elems.size()),
        static_cast<int>(proposed_coarse_mesh().elems.size()));
}

std::vector<int>
ReferenceEpochHierarchy::candidate_parent_proposed_coarse_elements() const {
    const Eigen::SparseMatrix<double> &embedding =
        proposed_coarse_elements_to_candidate();
    return fine_element_parents(
        embedding,
        static_cast<int>(candidate_mesh().elems.size()),
        static_cast<int>(proposed_coarse_mesh().elems.size()));
}

ReferenceEpochRefinementResult
ReferenceEpochHierarchy::close_candidate_over_proposed_coarse(
    const ReferenceEpochRefinementGuard &resource_guard) {
    ReferenceEpochRefinementResult result;
    result.previous_element_count = candidate_mesh().elems.size();
    result.current_element_count = candidate_mesh().elems.size();
    if (!proposed_coarse_) {
        result.detail = "there is no pending coarse refinement proposal";
        return result;
    }
    if (candidate_contains_proposed_coarse()) {
        result.detail = "candidate already contains the proposed coarse mesh";
        return result;
    }

    std::vector<int> refined_coarse_parents;
    const Eigen::SparseMatrix<double> &element_map =
        proposed_coarse_->refinement.P_elem;
    for (int parent = 0; parent < element_map.outerSize(); ++parent) {
        int descendants = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(
                 element_map, parent); it; ++it) {
            if (it.value() != 0.0) ++descendants;
        }
        if (descendants > 1) refined_coarse_parents.push_back(parent);
    }
    if (refined_coarse_parents.empty()) {
        throw std::runtime_error(
            "a non-nested coarse proposal has no refined parent elements");
    }

    constexpr int maximum_closure_steps = 32;
    for (int step = 0; step < maximum_closure_steps; ++step) {
        std::vector<char> refine_parent(coarse_mesh_.elems.size(), false);
        for (int parent : refined_coarse_parents) refine_parent[parent] = true;
        std::vector<int> marked;
        for (int element = 0;
             element < static_cast<int>(candidate_mesh().elems.size());
             ++element) {
            const int parent = ambient_parent_coarse_elements_[element];
            if (refine_parent[parent]) marked.push_back(element);
        }
        if (marked.empty()) {
            throw std::runtime_error(
                "candidate closure could not locate the proposed coarse region");
        }
        const ReferenceEpochRefinementResult enriched =
            enrich_candidate(marked, resource_guard);
        result.time_nvb_refine += enriched.time_nvb_refine;
        result.time_embedding_composition +=
            enriched.time_embedding_composition;
        result.time_parent_map_update += enriched.time_parent_map_update;
        result.time_candidate_quasi_interpolation +=
            enriched.time_candidate_quasi_interpolation;
        result.time_embedding_validation +=
            enriched.time_embedding_validation;
        if (candidate_contains_proposed_coarse()) {
            result.status = ReferenceEpochRefinementStatus::Refined;
            result.current_element_count = candidate_mesh().elems.size();
            result.detail =
                "candidate closure now contains the proposed coarse mesh";
            return result;
        }
    }
    throw std::runtime_error(
        "candidate closure did not contain the proposed coarse mesh");
}

ReferenceEpochRefinementResult
ReferenceEpochHierarchy::deepen_candidate_over_proposed_coarse(
    const int minimum_level_gap,
    const ReferenceEpochRefinementGuard &resource_guard) {
    std::vector<char> included;
    if (proposed_coarse_) {
        included.assign(proposed_coarse_->refinement.mesh.elems.size(), true);
    }
    return deepen_candidate_over_proposed_coarse(
        minimum_level_gap, included, resource_guard);
}

ReferenceEpochRefinementResult
ReferenceEpochHierarchy::deepen_candidate_over_proposed_coarse(
    const int minimum_level_gap,
    const std::vector<char> &included_proposed_elements,
    const ReferenceEpochRefinementGuard &resource_guard) {
    ReferenceEpochRefinementResult result;
    result.previous_element_count = candidate_mesh().elems.size();
    result.current_element_count = candidate_mesh().elems.size();
    if (minimum_level_gap <= 0) {
        result.detail = "candidate depth target is disabled";
        return result;
    }
    if (!proposed_coarse_ || !proposed_coarse_->to_candidate) {
        throw std::logic_error(
            "candidate depth target requires a contained coarse proposal");
    }
    if (included_proposed_elements.size()
        != proposed_coarse_->refinement.mesh.elems.size()) {
        throw std::invalid_argument(
            "candidate depth mask does not match the proposed mesh");
    }
    std::vector<int> target_level_gaps(
        included_proposed_elements.size(), 0);
    for (int element = 0;
         element < static_cast<int>(included_proposed_elements.size());
         ++element) {
        if (included_proposed_elements[element])
            target_level_gaps[element] = minimum_level_gap;
    }
    return deepen_candidate_over_proposed_coarse(
        target_level_gaps, resource_guard);
}

ReferenceEpochCandidateDeepeningProbe
ReferenceEpochHierarchy::probe_candidate_deepening_over_proposed_coarse(
    const std::vector<int> &target_level_gaps,
    const std::vector<char> &protected_proposed_elements,
    const std::vector<char> &full_target_proposed_elements,
    const ReferenceEpochRefinementGuard &resource_guard) const {
    if (!proposed_coarse_ || !proposed_coarse_->to_candidate) {
        throw std::logic_error(
            "candidate depth probe requires a contained coarse proposal");
    }
    const std::size_t proposed_size =
        proposed_coarse_->refinement.mesh.elems.size();
    if (target_level_gaps.size() != proposed_size
        || protected_proposed_elements.size() != proposed_size
        || full_target_proposed_elements.size() != proposed_size) {
        throw std::invalid_argument(
            "candidate depth probe data do not match the proposed mesh");
    }
    if (std::any_of(
            target_level_gaps.begin(), target_level_gaps.end(),
            [](const int target) { return target < 0; })) {
        throw std::invalid_argument(
            "candidate depth probe targets must be nonnegative");
    }

    TriMesh mesh = candidate_mesh();
    std::vector<int> levels = ambient_element_levels_;
    std::vector<int> proposed_parents =
        candidate_parent_proposed_coarse_elements();
    ReferenceEpochCandidateDeepeningProbe result;
    constexpr int maximum_deepening_steps = 32;
    for (int step = 0; step <= maximum_deepening_steps; ++step) {
        std::vector<int> marked;
        marked.reserve(mesh.elems.size());
        for (int child = 0; child < static_cast<int>(mesh.elems.size());
             ++child) {
            const int parent = proposed_parents[child];
            if (levels[child] - proposed_coarse_->levels[parent]
                < target_level_gaps[parent]) {
                marked.push_back(child);
            }
        }
        if (marked.empty()) {
            result.target_satisfied = true;
            break;
        }
        if (step == maximum_deepening_steps) break;
        if (resource_guard) {
            resource_guard(
                ReferenceEpochRefinementGuardPoint::BeforeNvb, mesh);
        }
        RefineOutput refined = bisect_newest_vertex(mesh, marked);
        if (resource_guard) {
            resource_guard(
                ReferenceEpochRefinementGuardPoint::AfterNvb,
                refined.mesh);
        }
        std::vector<int> next_levels = refinement_child_levels(
            mesh, levels, refined);
        const std::vector<int> step_parents = fine_element_parents(
            refined.P_elem, static_cast<int>(refined.mesh.elems.size()),
            static_cast<int>(mesh.elems.size()));
        std::vector<int> next_proposed_parents(step_parents.size());
        for (int child = 0; child < static_cast<int>(step_parents.size());
             ++child) {
            next_proposed_parents[child] =
                proposed_parents[step_parents[child]];
        }
        mesh = std::move(refined.mesh);
        levels = std::move(next_levels);
        proposed_parents = std::move(next_proposed_parents);
        ++result.refinement_steps;
    }

    std::vector<std::size_t> child_counts(proposed_size, 0);
    result.minimum_gap_margin = std::numeric_limits<int>::max();
    result.minimum_full_target_gap = std::numeric_limits<int>::max();
    for (int child = 0; child < static_cast<int>(mesh.elems.size()); ++child) {
        const int parent = proposed_parents[child];
        ++child_counts[parent];
        const int gap = levels[child] - proposed_coarse_->levels[parent];
        result.minimum_gap_margin = std::min(
            result.minimum_gap_margin, gap - target_level_gaps[parent]);
        if (full_target_proposed_elements[parent]) {
            result.minimum_full_target_gap = std::min(
                result.minimum_full_target_gap, gap);
        }
    }
    for (int parent = 0; parent < static_cast<int>(proposed_size); ++parent) {
        if (protected_proposed_elements[parent]
            && child_counts[parent] != 1) {
            ++result.protected_parent_spill;
        }
    }
    result.final_element_count = mesh.elems.size();
    result.final_node_count = mesh.nodes.size();
    return result;
}

ReferenceEpochRefinementResult
ReferenceEpochHierarchy::deepen_candidate_over_proposed_coarse(
    const std::vector<int> &target_level_gaps,
    const ReferenceEpochRefinementGuard &resource_guard) {
    ReferenceEpochRefinementResult result;
    result.previous_element_count = candidate_mesh().elems.size();
    result.current_element_count = candidate_mesh().elems.size();
    if (!proposed_coarse_ || !proposed_coarse_->to_candidate) {
        throw std::logic_error(
            "candidate depth target requires a contained coarse proposal");
    }
    if (target_level_gaps.size()
        != proposed_coarse_->refinement.mesh.elems.size()) {
        throw std::invalid_argument(
            "candidate depth targets do not match the proposed mesh");
    }
    if (std::any_of(
            target_level_gaps.begin(), target_level_gaps.end(),
            [](const int target) { return target < 0; })) {
        throw std::invalid_argument(
            "candidate depth targets must be nonnegative");
    }
    if (std::none_of(
            target_level_gaps.begin(), target_level_gaps.end(),
            [](const int target) { return target > 0; })) {
        result.detail = "candidate depth target is disabled";
        return result;
    }
    constexpr int maximum_deepening_steps = 32;
    for (int step = 0; step < maximum_deepening_steps; ++step) {
        const std::vector<int> parents = fine_element_parents(
            proposed_coarse_->to_candidate->P_elem,
            static_cast<int>(candidate_mesh().elems.size()),
            static_cast<int>(proposed_coarse_->refinement.mesh.elems.size()));
        std::vector<int> marked;
        marked.reserve(parents.size());
        for (int child = 0; child < static_cast<int>(parents.size()); ++child) {
            const int target = target_level_gaps[parents[child]];
            if (target > 0
                && ambient_element_levels_[child]
                    - proposed_coarse_->levels[parents[child]]
                < target) {
                marked.push_back(child);
            }
        }
        if (marked.empty()) {
            result.current_element_count = candidate_mesh().elems.size();
            result.detail =
                "candidate satisfies the graded post-refresh level-gap targets";
            return result;
        }
        const ReferenceEpochRefinementResult enriched =
            enrich_candidate(marked, resource_guard);
        result.status = ReferenceEpochRefinementStatus::Refined;
        result.current_element_count = enriched.current_element_count;
        result.time_nvb_refine += enriched.time_nvb_refine;
        result.time_embedding_composition +=
            enriched.time_embedding_composition;
        result.time_parent_map_update += enriched.time_parent_map_update;
        result.time_candidate_quasi_interpolation +=
            enriched.time_candidate_quasi_interpolation;
        result.time_embedding_validation +=
            enriched.time_embedding_validation;
    }
    throw std::runtime_error(
        "candidate failed to reach the post-refresh level-gap target");
}

ReferenceEpochRefinementResult
ReferenceEpochHierarchy::commit_coarse_refinement() {
    ReferenceEpochRefinementResult result;
    result.previous_element_count = coarse_mesh_.elems.size();
    result.current_element_count = coarse_mesh_.elems.size();
    if (!proposed_coarse_) {
        result.detail = "there is no pending coarse refinement proposal";
        return result;
    }
    if (!proposed_coarse_->to_reference
        || !proposed_coarse_->to_candidate) {
        result.status =
            ReferenceEpochRefinementStatus::ReferenceRefreshRequired;
        result.detail =
            "cached proposed embedding is not contained in reference/candidate";
        return result;
    }

    if (proposed_coarse_->identity) {
        proposed_coarse_.reset();
        result.status = ReferenceEpochRefinementStatus::Refined;
        result.detail = "identity coarse proposal committed";
        return result;
    }

    ProposedCoarseState proposed = std::move(*proposed_coarse_);
    proposed_coarse_.reset();
    coarse_mesh_ = std::move(proposed.refinement.mesh);
    coarse_levels_ = std::move(proposed.levels);
    reference_completion_.refinement.P_node =
        std::move(proposed.to_reference->P_node);
    reference_completion_.refinement.P_elem =
        std::move(proposed.to_reference->P_elem);
    reference_completion_.refinement.P_dg =
        std::move(proposed.to_reference->P_dg);
    coarse_to_ambient_.P_node =
        std::move(proposed.to_candidate->P_node);
    coarse_to_ambient_.P_elem =
        std::move(proposed.to_candidate->P_elem);
    coarse_to_ambient_.P_dg =
        std::move(proposed.to_candidate->P_dg);
    reference_parent_coarse_elements_ = fine_element_parents(
        reference_completion_.refinement.P_elem,
        static_cast<int>(reference_mesh().elems.size()),
        static_cast<int>(coarse_mesh_.elems.size()));
    ambient_parent_coarse_elements_ = fine_element_parents(
        coarse_to_ambient_.P_elem,
        static_cast<int>(candidate_mesh().elems.size()),
        static_cast<int>(coarse_mesh_.elems.size()));
    reference_quasi_interpolation_ = build_quasi_interp(
        coarse_mesh_, reference_mesh(),
        reference_completion_.refinement.P_dg,
        cg_to_dg(reference_mesh()),
        static_cast<int>(reference_mesh().nodes.size()),
        static_cast<int>(coarse_mesh_.nodes.size()));
    ambient_quasi_interpolation_ = build_quasi_interp(
        coarse_mesh_, candidate_mesh(), coarse_to_ambient_.P_dg,
        cg_to_dg(candidate_mesh()),
        static_cast<int>(candidate_mesh().nodes.size()),
        static_cast<int>(coarse_mesh_.nodes.size()));
    ++coarse_mesh_version_;
    ++interpolation_version_;
    ++boundary_version_;
    ++corrector_space_version_;
    // Cached proposed embeddings are exact.  Production commit validates
    // their dimensions/parents without repeating global geometric searches
    // or the full sparse composition audit exercised by hierarchy tests.
    validate_incremental_candidate_metadata();
    result.status = ReferenceEpochRefinementStatus::Refined;
    result.current_element_count = coarse_mesh_.elems.size();
    result.detail = "proposed coarse refinement committed";
    return result;
}

void ReferenceEpochHierarchy::refresh_reference_from_candidate() {
    refresh_reference_from_ambient();
}

} // namespace lod2d::helmholtz::adaptive

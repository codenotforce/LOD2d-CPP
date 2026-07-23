#include "helmholtz/hp_interpolation.h"
#include "mesh/refine.h"

#include <cmath>
#include <stdexcept>
#include <unordered_map>

namespace lod2d::helmholtz {
namespace {

Eigen::Vector3d barycentric_coordinates(
    const TriMesh &mesh,
    const Triangle &triangle,
    const Point2 &point) {
    Eigen::Matrix2d jacobian;
    jacobian.col(0) =
        mesh.nodes[triangle[1]] - mesh.nodes[triangle[0]];
    jacobian.col(1) =
        mesh.nodes[triangle[2]] - mesh.nodes[triangle[0]];
    if (std::abs(jacobian.determinant()) <= 1e-15)
        throw std::invalid_argument(
            "barycentric evaluation encountered a degenerate triangle");
    const Eigen::Vector2d reference =
        jacobian.inverse() * (point - mesh.nodes[triangle[0]]);
    return {1.0 - reference.x() - reference.y(),
            reference.x(), reference.y()};
}

std::vector<int> fine_parents(
    const Eigen::SparseMatrix<double> &prolongation,
    int fine_elements,
    int coarse_elements) {
    if (prolongation.rows() != fine_elements
        || prolongation.cols() != coarse_elements)
        throw std::invalid_argument(
            "fine element prolongation has incompatible dimensions");
    std::vector<int> parents(fine_elements, -1);
    for (int coarse = 0; coarse < prolongation.outerSize(); ++coarse) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(
                 prolongation, coarse); it; ++it) {
            if (it.value() == 0.0) continue;
            if (parents[it.row()] >= 0)
                throw std::invalid_argument(
                    "fine element belongs to multiple coarse elements");
            parents[it.row()] = coarse;
        }
    }
    for (int parent : parents)
        if (parent < 0)
            throw std::invalid_argument(
                "fine element is missing its coarse parent");
    return parents;
}

} // namespace

HelmholtzHpInterpolation build_helmholtz_hp_interpolation(
    const TriMesh &coarse,
    const HpTriSpace &fine_space,
    const Eigen::SparseMatrix<double> &fine_element_prolongation) {
    const TriMesh &fine = fine_space.mesh();
    const int coarse_elements = static_cast<int>(coarse.elems.size());
    const int coarse_nodes = static_cast<int>(coarse.nodes.size());
    const int fine_elements = static_cast<int>(fine.elems.size());
    const int hp_dofs = fine_space.dof_count();
    const int local_size = fine_space.local_dof_count();
    const auto parents = fine_parents(
        fine_element_prolongation, fine_elements, coarse_elements);

    std::vector<std::unordered_map<int, double>> injection_rows(hp_dofs);
    std::vector<int> injection_visits(hp_dofs, 0);
    for (int element = 0; element < fine_elements; ++element) {
        const int parent = parents[element];
        const Triangle &coarse_triangle = coarse.elems[parent];
        const auto &dofs = fine_space.element_dofs()[element];
        for (int local = 0; local < local_size; ++local) {
            const int dof = dofs[local];
            const Eigen::Vector3d coarse_barycentric =
                barycentric_coordinates(
                    coarse, coarse_triangle,
                    fine_space.dof_points()[dof]);
            for (int i = 0; i < 3; ++i)
                injection_rows[dof][coarse_triangle[i]] +=
                    coarse_barycentric(i);
            ++injection_visits[dof];
        }
    }

    std::vector<Eigen::Triplet<double>> injection_triplets;
    for (int dof = 0; dof < hp_dofs; ++dof) {
        if (injection_visits[dof] == 0)
            throw std::runtime_error("hp DOF is not incident to any fine element");
        for (const auto &[coarse_node, accumulated] : injection_rows[dof]) {
            const double value = accumulated / injection_visits[dof];
            if (std::abs(value) > 1e-14)
                injection_triplets.emplace_back(dof, coarse_node, value);
        }
    }
    Eigen::SparseMatrix<double> injection(hp_dofs, coarse_nodes);
    injection.setFromTriplets(
        injection_triplets.begin(), injection_triplets.end());

    std::vector<Eigen::Triplet<double>> moment_triplets;
    const auto quadrature =
        triangle_gauss_quadrature(fine_space.degree() + 2);
    for (int element = 0; element < fine_elements; ++element) {
        const Triangle &fine_triangle = fine.elems[element];
        const int parent = parents[element];
        const Triangle &coarse_triangle = coarse.elems[parent];
        Eigen::Matrix2d jacobian;
        jacobian.col(0) =
            fine.nodes[fine_triangle[1]] - fine.nodes[fine_triangle[0]];
        jacobian.col(1) =
            fine.nodes[fine_triangle[2]] - fine.nodes[fine_triangle[0]];
        const double determinant = std::abs(jacobian.determinant());
        const auto &dofs = fine_space.element_dofs()[element];
        for (const auto &point : quadrature) {
            Point2 physical = Point2::Zero();
            for (int i = 0; i < 3; ++i)
                physical += point.barycentric(i)
                    * fine.nodes[fine_triangle[i]];
            const Eigen::Vector3d coarse_barycentric =
                barycentric_coordinates(coarse, coarse_triangle, physical);
            const auto hp_evaluation =
                fine_space.basis().evaluate(point.barycentric);
            const double weight = determinant * point.weight;
            for (int coarse_local = 0; coarse_local < 3; ++coarse_local) {
                for (int hp_local = 0; hp_local < local_size; ++hp_local) {
                    const double value = weight
                        * coarse_barycentric(coarse_local)
                        * hp_evaluation.values(hp_local);
                    if (std::abs(value) > 1e-15)
                        moment_triplets.emplace_back(
                            3 * parent + coarse_local,
                            dofs[hp_local], value);
                }
            }
        }
    }
    Eigen::SparseMatrix<double> moments(
        3 * coarse_elements, hp_dofs);
    moments.setFromTriplets(
        moment_triplets.begin(), moment_triplets.end());

    const auto coarse_areas = compute_area(coarse);
    static constexpr double inverse_mass[3][3] = {
        {9.0, -3.0, -3.0},
        {-3.0, 9.0, -3.0},
        {-3.0, -3.0, 9.0}
    };
    std::vector<Eigen::Triplet<double>> inverse_mass_triplets;
    for (int element = 0; element < coarse_elements; ++element) {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                inverse_mass_triplets.emplace_back(
                    3 * element + i, 3 * element + j,
                    inverse_mass[i][j] / coarse_areas[element]);
    }
    Eigen::SparseMatrix<double> inverse_dg_mass(
        3 * coarse_elements, 3 * coarse_elements);
    inverse_dg_mass.setFromTriplets(
        inverse_mass_triplets.begin(), inverse_mass_triplets.end());

    std::vector<int> coarse_incidence(coarse_nodes, 0);
    for (const Triangle &triangle : coarse.elems)
        for (int vertex : triangle) ++coarse_incidence[vertex];
    std::vector<Eigen::Triplet<double>> averaging_triplets;
    for (int element = 0; element < coarse_elements; ++element) {
        for (int local = 0; local < 3; ++local) {
            const int vertex = coarse.elems[element][local];
            averaging_triplets.emplace_back(
                vertex, 3 * element + local,
                1.0 / coarse_incidence[vertex]);
        }
    }
    Eigen::SparseMatrix<double> averaging(
        coarse_nodes, 3 * coarse_elements);
    averaging.setFromTriplets(
        averaging_triplets.begin(), averaging_triplets.end());

    Eigen::SparseMatrix<double> quasi_interpolation =
        averaging * inverse_dg_mass * moments;
    quasi_interpolation.makeCompressed();
    const Eigen::MatrixXd reproduction =
        Eigen::MatrixXd(quasi_interpolation * injection);
    const Eigen::MatrixXd identity =
        Eigen::MatrixXd::Identity(coarse_nodes, coarse_nodes);
    if ((reproduction - identity).norm() > 1e-9)
        throw std::runtime_error(
            "hp quasi-interpolation does not reproduce coarse P1");

    return {std::move(injection), std::move(quasi_interpolation)};
}

} // namespace lod2d::helmholtz

#include "helmholtz/hp_operators.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace lod2d::helmholtz {
namespace {

std::vector<double> resolve_coefficients(
    const std::vector<double> &values,
    int count,
    const char *name) {
    if (values.empty()) return std::vector<double>(count, 1.0);
    if (static_cast<int>(values.size()) != count)
        throw std::invalid_argument(std::string(name) + " count must match elements");
    return values;
}

Eigen::Vector3d edge_barycentric(int local_edge, double coordinate) {
    switch (local_edge) {
    case 0:
        return {1.0 - coordinate, coordinate, 0.0};
    case 1:
        return {0.0, 1.0 - coordinate, coordinate};
    case 2:
        return {coordinate, 0.0, 1.0 - coordinate};
    default:
        throw std::out_of_range("invalid local triangle edge");
    }
}

} // namespace

HelmholtzHpOperators assemble_helmholtz_hp_operators(
    const HpTriSpace &space,
    double wavenumber,
    const std::vector<double> &diffusion,
    const std::vector<double> &refractive_index,
    double boundary_beta) {
    if (wavenumber <= 0.0)
        throw std::invalid_argument("Helmholtz wavenumber must be positive");
    const TriMesh &mesh = space.mesh();
    const int element_count = static_cast<int>(mesh.elems.size());
    const int local_size = space.local_dof_count();
    const auto diffusion_values =
        resolve_coefficients(diffusion, element_count, "diffusion");
    const auto refractive_values =
        resolve_coefficients(refractive_index, element_count, "refractive index");
    const auto triangle_quadrature =
        triangle_gauss_quadrature(space.degree() + 2);
    const auto edge_quadrature = edge_gauss_quadrature(space.degree() + 2);

    std::vector<Eigen::Triplet<double>> stiffness_triplets;
    std::vector<Eigen::Triplet<double>> mass_triplets;
    std::vector<Eigen::Triplet<double>> robin_triplets;
    const std::size_t local_entries =
        static_cast<std::size_t>(element_count) * local_size * local_size;
    stiffness_triplets.reserve(local_entries);
    mass_triplets.reserve(local_entries);
    robin_triplets.reserve(local_entries / 4 + 1);

    HelmholtzHpOperators result;
    result.wavenumber = wavenumber;
    result.boundary_beta = boundary_beta;
    result.diffusion = diffusion_values;
    result.refractive_index = refractive_values;
    result.element_blocks.resize(element_count);

    static constexpr int local_edges[3][2] = {{0, 1}, {1, 2}, {2, 0}};
    for (int element = 0; element < element_count; ++element) {
        const Triangle &triangle = mesh.elems[element];
        Eigen::Matrix2d jacobian;
        jacobian.col(0) = mesh.nodes[triangle[1]] - mesh.nodes[triangle[0]];
        jacobian.col(1) = mesh.nodes[triangle[2]] - mesh.nodes[triangle[0]];
        const double determinant = jacobian.determinant();
        if (std::abs(determinant) <= 1e-15)
            throw std::invalid_argument(
                "hp assembly encountered a degenerate triangle");
        const double absolute_determinant = std::abs(determinant);
        const Eigen::Matrix2d inverse_transpose =
            jacobian.inverse().transpose();

        Eigen::MatrixXd local_stiffness =
            Eigen::MatrixXd::Zero(local_size, local_size);
        Eigen::MatrixXd local_mass =
            Eigen::MatrixXd::Zero(local_size, local_size);
        Eigen::MatrixXd local_robin =
            Eigen::MatrixXd::Zero(local_size, local_size);
        for (const auto &quadrature : triangle_quadrature) {
            const auto evaluation =
                space.basis().evaluate(quadrature.barycentric);
            const Eigen::MatrixXd physical_gradients =
                evaluation.reference_gradients
                * inverse_transpose.transpose();
            const double weight =
                absolute_determinant * quadrature.weight;
            local_stiffness.noalias() += diffusion_values[element] * weight
                * physical_gradients * physical_gradients.transpose();
            local_mass.noalias() += refractive_values[element] * weight
                * evaluation.values * evaluation.values.transpose();
        }

        for (int local_edge = 0; local_edge < 3; ++local_edge) {
            const int a = triangle[local_edges[local_edge][0]];
            const int b = triangle[local_edges[local_edge][1]];
            const int edge = space.edge_index(a, b);
            if (!space.boundary_edges()[edge]) continue;
            const double length = (mesh.nodes[a] - mesh.nodes[b]).norm();
            for (const auto &quadrature : edge_quadrature) {
                const auto evaluation = space.basis().evaluate(
                    edge_barycentric(local_edge, quadrature.coordinate));
                local_robin.noalias() +=
                    boundary_beta * length * quadrature.weight
                    * evaluation.values * evaluation.values.transpose();
            }
        }

        const auto &dofs = space.element_dofs()[element];
        for (int i = 0; i < local_size; ++i) {
            for (int j = 0; j < local_size; ++j) {
                if (local_stiffness(i, j) != 0.0)
                    stiffness_triplets.emplace_back(
                        dofs[i], dofs[j], local_stiffness(i, j));
                if (local_mass(i, j) != 0.0)
                    mass_triplets.emplace_back(
                        dofs[i], dofs[j], local_mass(i, j));
                if (local_robin(i, j) != 0.0)
                    robin_triplets.emplace_back(
                        dofs[i], dofs[j], local_robin(i, j));
            }
        }
        result.element_blocks[element] =
            local_stiffness.cast<Complex>()
            - (wavenumber * wavenumber) * local_mass.cast<Complex>()
            - Complex(0.0, wavenumber) * local_robin.cast<Complex>();
    }

    result.stiffness.resize(space.dof_count(), space.dof_count());
    result.stiffness.setFromTriplets(
        stiffness_triplets.begin(), stiffness_triplets.end());
    result.mass.resize(space.dof_count(), space.dof_count());
    result.mass.setFromTriplets(mass_triplets.begin(), mass_triplets.end());
    result.boundary_mass.resize(space.dof_count(), space.dof_count());
    result.boundary_mass.setFromTriplets(
        robin_triplets.begin(), robin_triplets.end());
    result.system = result.stiffness.cast<Complex>();
    result.system -=
        wavenumber * wavenumber * result.mass.cast<Complex>();
    result.system -= Complex(0.0, wavenumber)
        * result.boundary_mass.cast<Complex>();
    result.system.makeCompressed();
    return result;
}

} // namespace lod2d::helmholtz

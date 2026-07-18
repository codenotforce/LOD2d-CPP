#include "helmholtz/operators.h"

#include <Eigen/SparseLU>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>

namespace lod2d::helmholtz {
namespace {

struct QuadraturePoint {
    std::array<double, 3> barycentric;
    double weight;
};

constexpr std::array<QuadraturePoint, 7> kTriangleQuadrature{{
    {{{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}}, 0.225000000000000},
    {{{0.059715871789770, 0.470142064105115, 0.470142064105115}}, 0.132394152788506},
    {{{0.470142064105115, 0.059715871789770, 0.470142064105115}}, 0.132394152788506},
    {{{0.470142064105115, 0.470142064105115, 0.059715871789770}}, 0.132394152788506},
    {{{0.797426985353087, 0.101286507323456, 0.101286507323456}}, 0.125939180544827},
    {{{0.101286507323456, 0.797426985353087, 0.101286507323456}}, 0.125939180544827},
    {{{0.101286507323456, 0.101286507323456, 0.797426985353087}}, 0.125939180544827}
}};

std::uint64_t edge_key(int a, int b) {
    const auto lo = static_cast<std::uint32_t>(std::min(a, b));
    const auto hi = static_cast<std::uint32_t>(std::max(a, b));
    return (static_cast<std::uint64_t>(lo) << 32U) | hi;
}

double triangle_geometry(
    const TriMesh &mesh,
    const Triangle &tri,
    std::array<Eigen::Vector2d, 3> &gradients) {
    const Point2 &a = mesh.nodes[tri[0]];
    const Point2 &b = mesh.nodes[tri[1]];
    const Point2 &c = mesh.nodes[tri[2]];
    const double det = (b.x() - a.x()) * (c.y() - a.y())
                     - (b.y() - a.y()) * (c.x() - a.x());
    if (std::abs(det) <= 1e-15)
        throw std::invalid_argument("Helmholtz assembly encountered a degenerate triangle");

    gradients[0] = Eigen::Vector2d(b.y() - c.y(), c.x() - b.x()) / det;
    gradients[1] = Eigen::Vector2d(c.y() - a.y(), a.x() - c.x()) / det;
    gradients[2] = Eigen::Vector2d(a.y() - b.y(), b.x() - a.x()) / det;
    return 0.5 * std::abs(det);
}

std::vector<double> resolved_coefficients(
    const std::vector<double> &values,
    int element_count,
    const char *name) {
    if (values.empty()) return std::vector<double>(element_count, 1.0);
    if (static_cast<int>(values.size()) != element_count)
        throw std::invalid_argument(std::string(name) + " coefficient count must match element count");
    return values;
}

} // namespace

HelmholtzOperators assemble_helmholtz_operators(
    const TriMesh &mesh,
    double wavenumber,
    const std::vector<double> &diffusion,
    const std::vector<double> &refractive_index,
    double boundary_beta) {
    if (wavenumber <= 0.0)
        throw std::invalid_argument("Helmholtz wavenumber must be positive");
    if (mesh.nodes.empty() || mesh.elems.empty())
        throw std::invalid_argument("Helmholtz mesh must not be empty");

    const int node_count = static_cast<int>(mesh.nodes.size());
    const int element_count = static_cast<int>(mesh.elems.size());
    const auto diffusion_values = resolved_coefficients(diffusion, element_count, "diffusion");
    const auto refractive_values = resolved_coefficients(refractive_index, element_count, "refractive index");

    std::unordered_map<std::uint64_t, int> edge_counts;
    edge_counts.reserve(static_cast<size_t>(3 * element_count));
    for (const auto &tri : mesh.elems) {
        ++edge_counts[edge_key(tri[0], tri[1])];
        ++edge_counts[edge_key(tri[1], tri[2])];
        ++edge_counts[edge_key(tri[2], tri[0])];
    }

    std::vector<Eigen::Triplet<double>> stiffness_triplets;
    std::vector<Eigen::Triplet<double>> mass_triplets;
    std::vector<Eigen::Triplet<double>> boundary_triplets;
    stiffness_triplets.reserve(static_cast<size_t>(9 * element_count));
    mass_triplets.reserve(static_cast<size_t>(9 * element_count));
    boundary_triplets.reserve(static_cast<size_t>(4 * std::sqrt(element_count) + 16));

    HelmholtzOperators result;
    result.wavenumber = wavenumber;
    result.diffusion = diffusion_values;
    result.refractive_index = refractive_values;
    result.boundary_beta = boundary_beta;
    result.element_blocks.resize(element_count, Eigen::Matrix3cd::Zero());

    static constexpr double mass_pattern[3][3] = {
        {2.0, 1.0, 1.0},
        {1.0, 2.0, 1.0},
        {1.0, 1.0, 2.0}
    };
    static constexpr int local_edges[3][2] = {{0, 1}, {1, 2}, {2, 0}};

    for (int element = 0; element < element_count; ++element) {
        const Triangle &tri = mesh.elems[element];
        std::array<Eigen::Vector2d, 3> gradients;
        const double area = triangle_geometry(mesh, tri, gradients);
        Eigen::Matrix3d local_stiffness;
        Eigen::Matrix3d local_mass;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                local_stiffness(i, j) = diffusion_values[element] * area
                                      * gradients[i].dot(gradients[j]);
                local_mass(i, j) = refractive_values[element] * area
                                 * mass_pattern[i][j] / 12.0;
                stiffness_triplets.emplace_back(tri[i], tri[j], local_stiffness(i, j));
                mass_triplets.emplace_back(tri[i], tri[j], local_mass(i, j));
            }
        }

        Eigen::Matrix3d local_boundary = Eigen::Matrix3d::Zero();
        for (const auto &local_edge : local_edges) {
            const int i = local_edge[0];
            const int j = local_edge[1];
            if (edge_counts.at(edge_key(tri[i], tri[j])) != 1) continue;
            const double length = (mesh.nodes[tri[i]] - mesh.nodes[tri[j]]).norm();
            const double diagonal = boundary_beta * length / 3.0;
            const double off_diagonal = boundary_beta * length / 6.0;
            local_boundary(i, i) += diagonal;
            local_boundary(j, j) += diagonal;
            local_boundary(i, j) += off_diagonal;
            local_boundary(j, i) += off_diagonal;
            boundary_triplets.emplace_back(tri[i], tri[i], diagonal);
            boundary_triplets.emplace_back(tri[j], tri[j], diagonal);
            boundary_triplets.emplace_back(tri[i], tri[j], off_diagonal);
            boundary_triplets.emplace_back(tri[j], tri[i], off_diagonal);
        }

        result.element_blocks[element] = local_stiffness.cast<Complex>()
            - (wavenumber * wavenumber) * local_mass.cast<Complex>()
            - Complex(0.0, wavenumber) * local_boundary.cast<Complex>();
    }

    result.stiffness.resize(node_count, node_count);
    result.stiffness.setFromTriplets(stiffness_triplets.begin(), stiffness_triplets.end());
    result.mass.resize(node_count, node_count);
    result.mass.setFromTriplets(mass_triplets.begin(), mass_triplets.end());
    result.boundary_mass.resize(node_count, node_count);
    result.boundary_mass.setFromTriplets(boundary_triplets.begin(), boundary_triplets.end());

    result.system = result.stiffness.cast<Complex>();
    result.system -= (wavenumber * wavenumber) * result.mass.cast<Complex>();
    result.system -= Complex(0.0, wavenumber) * result.boundary_mass.cast<Complex>();
    result.system.makeCompressed();
    return result;
}

ComplexVector assemble_helmholtz_load(
    const TriMesh &mesh,
    const ComplexFunction &source) {
    if (!source) throw std::invalid_argument("Helmholtz source function is empty");
    ComplexVector load = ComplexVector::Zero(static_cast<int>(mesh.nodes.size()));

    for (const Triangle &tri : mesh.elems) {
        std::array<Eigen::Vector2d, 3> gradients;
        const double area = triangle_geometry(mesh, tri, gradients);
        for (const auto &quadrature : kTriangleQuadrature) {
            Point2 point = Point2::Zero();
            for (int i = 0; i < 3; ++i)
                point += quadrature.barycentric[i] * mesh.nodes[tri[i]];
            const Complex value = source(point);
            for (int i = 0; i < 3; ++i) {
                load(tri[i]) += area * quadrature.weight
                              * value * quadrature.barycentric[i];
            }
        }
    }
    return load;
}

ComplexVector solve_helmholtz_fem(
    const HelmholtzOperators &operators,
    const ComplexVector &load) {
    if (operators.system.rows() != load.size())
        throw std::invalid_argument("Helmholtz load size does not match the system matrix");
    Eigen::SparseLU<ComplexSparseMatrix> solver;
    solver.analyzePattern(operators.system);
    solver.factorize(operators.system);
    if (solver.info() != Eigen::Success)
        throw std::runtime_error("Helmholtz sparse LU factorization failed");
    ComplexVector solution = solver.solve(load);
    if (solver.info() != Eigen::Success || !solution.allFinite())
        throw std::runtime_error("Helmholtz sparse LU solve failed");
    return solution;
}

HelmholtzError compute_helmholtz_error(
    const TriMesh &mesh,
    const ComplexVector &solution,
    double wavenumber,
    const ComplexFunction &exact,
    const ComplexGradientFunction &exact_gradient) {
    if (solution.size() != static_cast<int>(mesh.nodes.size()))
        throw std::invalid_argument("Helmholtz solution size does not match mesh nodes");
    if (!exact || !exact_gradient)
        throw std::invalid_argument("Helmholtz exact solution callbacks must not be empty");

    double l2_squared = 0.0;
    double gradient_squared = 0.0;
    for (const Triangle &tri : mesh.elems) {
        std::array<Eigen::Vector2d, 3> gradients;
        const double area = triangle_geometry(mesh, tri, gradients);
        Eigen::Vector2cd discrete_gradient = Eigen::Vector2cd::Zero();
        for (int i = 0; i < 3; ++i)
            discrete_gradient += solution(tri[i]) * gradients[i].cast<Complex>();

        for (const auto &quadrature : kTriangleQuadrature) {
            Point2 point = Point2::Zero();
            Complex discrete_value = 0.0;
            for (int i = 0; i < 3; ++i) {
                point += quadrature.barycentric[i] * mesh.nodes[tri[i]];
                discrete_value += quadrature.barycentric[i] * solution(tri[i]);
            }
            const Complex value_error = exact(point) - discrete_value;
            const Eigen::Vector2cd gradient_error = exact_gradient(point) - discrete_gradient;
            l2_squared += area * quadrature.weight * std::norm(value_error);
            gradient_squared += area * quadrature.weight * gradient_error.squaredNorm();
        }
    }

    HelmholtzError error;
    error.l2 = std::sqrt(std::max(0.0, l2_squared));
    error.energy = std::sqrt(std::max(0.0, gradient_squared
        + wavenumber * wavenumber * l2_squared));
    return error;
}

double max_element_diameter(const TriMesh &mesh) {
    double diameter = 0.0;
    for (const Triangle &tri : mesh.elems) {
        diameter = std::max(diameter, (mesh.nodes[tri[0]] - mesh.nodes[tri[1]]).norm());
        diameter = std::max(diameter, (mesh.nodes[tri[1]] - mesh.nodes[tri[2]]).norm());
        diameter = std::max(diameter, (mesh.nodes[tri[2]] - mesh.nodes[tri[0]]).norm());
    }
    return diameter;
}

} // namespace lod2d::helmholtz

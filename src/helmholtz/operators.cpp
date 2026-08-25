#include "helmholtz/operators.h"
#include "helmholtz/boundary.h"

#include <Eigen/SparseLU>
#if defined(LOD2D_HAVE_UMFPACK)
#include <Eigen/UmfPackSupport>
#endif

#include <chrono>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <exception>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace lod2d::helmholtz {

namespace {
using FemSolveClock = std::chrono::steady_clock;
double fem_seconds_since(const FemSolveClock::time_point start) {
    return std::chrono::duration<double>(FemSolveClock::now() - start).count();
}
} // namespace

const char *helmholtz_fem_solver_kind_name(const HelmholtzFemSolverKind kind) {
    switch (kind) {
    case HelmholtzFemSolverKind::SparseLu: return "sparse_lu";
    case HelmholtzFemSolverKind::Umfpack: return "umfpack";
    }
    throw std::invalid_argument("unknown Helmholtz FEM solver kind");
}

bool helmholtz_fem_solver_available(const HelmholtzFemSolverKind kind) {
    if (kind == HelmholtzFemSolverKind::SparseLu) return true;
#if defined(LOD2D_HAVE_UMFPACK)
    if (kind == HelmholtzFemSolverKind::Umfpack) return true;
#endif
    return false;
}
namespace {

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
    validate_boundary_tags(mesh);

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
    result.dirichlet_nodes = dirichlet_nodes(mesh);
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
            if (boundary_tag(mesh, canonical_edge(tri[i], tri[j]))
                != BoundaryTag::Robin) continue;
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
    const ComplexFunction &source,
    const QuadraturePolicy &quadrature,
    const QuadratureContext &quadrature_context) {
    if (!source) throw std::invalid_argument("Helmholtz source function is empty");
    const int element_count = static_cast<int>(mesh.elems.size());
    std::vector<std::array<Complex, 3>> element_loads(
        element_count, std::array<Complex, 3>{});
    std::vector<std::exception_ptr> element_errors(element_count);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if(element_count >= 64)
#endif
    for (int element = 0; element < element_count; ++element) {
        try {
            for (const auto &point : triangle_quadrature_points(
                     mesh, element, quadrature, quadrature_context)) {
                const Complex value = source(point.point);
                for (int i = 0; i < 3; ++i) {
                    element_loads[element][i] +=
                        point.weight * value * point.barycentric[i];
                }
            }
        } catch (...) {
            element_errors[element] = std::current_exception();
        }
    }
    for (const std::exception_ptr &error : element_errors)
        if (error) std::rethrow_exception(error);

    // Scatter in element order. This preserves bitwise reproducibility across
    // OpenMP thread counts while parallelizing the expensive quadrature and
    // manufactured-source evaluation.
    ComplexVector load = ComplexVector::Zero(static_cast<int>(mesh.nodes.size()));
    for (int element = 0; element < element_count; ++element) {
        const Triangle &tri = mesh.elems[element];
        for (int i = 0; i < 3; ++i)
            load(tri[i]) += element_loads[element][i];
    }
    return load;
}

ComplexVector solve_helmholtz_fem(
    const HelmholtzOperators &operators,
    const ComplexVector &load,
    const HelmholtzFemSolverKind solver_kind,
    HelmholtzFemSolveTimings *timings) {
    const FemSolveClock::time_point total_start = FemSolveClock::now();
    if (timings) *timings = {};
    if (!helmholtz_fem_solver_available(solver_kind))
        throw std::runtime_error(
            std::string("requested Helmholtz FEM solver is unavailable: ")
            + helmholtz_fem_solver_kind_name(solver_kind));
    if (operators.system.rows() != load.size())
        throw std::invalid_argument("Helmholtz load size does not match the system matrix");
    std::vector<char> is_dirichlet(load.size(), false);
    for (int node : operators.dirichlet_nodes) {
        if (node < 0 || node >= load.size())
            throw std::invalid_argument("Dirichlet node index is out of range");
        is_dirichlet[node] = true;
    }
    std::vector<int> free_nodes;
    std::vector<int> global_to_free(load.size(), -1);
    for (int node = 0; node < load.size(); ++node) {
        if (is_dirichlet[node]) continue;
        global_to_free[node] = static_cast<int>(free_nodes.size());
        free_nodes.push_back(node);
    }
    if (free_nodes.empty())
        throw std::invalid_argument("Helmholtz FEM has no unconstrained degrees of freedom");

    const FemSolveClock::time_point reduction_start = FemSolveClock::now();
    std::vector<ComplexTriplet> triplets;
    for (int global_col : free_nodes) {
        const int local_col = global_to_free[global_col];
        for (ComplexSparseMatrix::InnerIterator it(operators.system, global_col); it; ++it) {
            const int local_row = global_to_free[it.row()];
            if (local_row >= 0) triplets.emplace_back(local_row, local_col, it.value());
        }
    }
    ComplexSparseMatrix reduced(free_nodes.size(), free_nodes.size());
    reduced.setFromTriplets(triplets.begin(), triplets.end());
    ComplexVector reduced_load(free_nodes.size());
    for (int local = 0; local < static_cast<int>(free_nodes.size()); ++local)
        reduced_load(local) = load(free_nodes[local]);
    reduced.makeCompressed();
    if (timings) timings->reduction_seconds = fem_seconds_since(reduction_start);

    ComplexVector reduced_solution;
    if (solver_kind == HelmholtzFemSolverKind::SparseLu) {
        Eigen::SparseLU<ComplexSparseMatrix> solver;
        const FemSolveClock::time_point analysis_start = FemSolveClock::now();
        solver.analyzePattern(reduced);
        if (timings) timings->analysis_seconds = fem_seconds_since(analysis_start);
        const FemSolveClock::time_point factorization_start = FemSolveClock::now();
        solver.factorize(reduced);
        if (timings)
            timings->factorization_seconds = fem_seconds_since(factorization_start);
        if (solver.info() != Eigen::Success)
            throw std::runtime_error("Helmholtz sparse LU factorization failed");
        const FemSolveClock::time_point solve_start = FemSolveClock::now();
        reduced_solution = solver.solve(reduced_load);
        if (timings) timings->solve_seconds = fem_seconds_since(solve_start);
        if (solver.info() != Eigen::Success || !reduced_solution.allFinite())
            throw std::runtime_error("Helmholtz sparse LU solve failed");
    } else if (solver_kind == HelmholtzFemSolverKind::Umfpack) {
#if defined(LOD2D_HAVE_UMFPACK)
        Eigen::UmfPackLU<ComplexSparseMatrix> solver;
        const FemSolveClock::time_point analysis_start = FemSolveClock::now();
        solver.analyzePattern(reduced);
        if (timings) timings->analysis_seconds = fem_seconds_since(analysis_start);
        const FemSolveClock::time_point factorization_start = FemSolveClock::now();
        solver.factorize(reduced);
        if (timings)
            timings->factorization_seconds = fem_seconds_since(factorization_start);
        if (solver.info() != Eigen::Success)
            throw std::runtime_error("Helmholtz UMFPACK factorization failed");
        const FemSolveClock::time_point solve_start = FemSolveClock::now();
        reduced_solution = solver.solve(reduced_load);
        if (timings) timings->solve_seconds = fem_seconds_since(solve_start);
        if (solver.info() != Eigen::Success || !reduced_solution.allFinite())
            throw std::runtime_error("Helmholtz UMFPACK solve failed");
#else
        throw std::runtime_error("Helmholtz UMFPACK support was not compiled");
#endif
    } else {
        throw std::invalid_argument("unknown Helmholtz FEM solver kind");
    }
    ComplexVector solution = ComplexVector::Zero(load.size());
    for (int local = 0; local < static_cast<int>(free_nodes.size()); ++local)
        solution(free_nodes[local]) = reduced_solution(local);
    if (timings) timings->total_seconds = fem_seconds_since(total_start);
    return solution;
}

HelmholtzError compute_helmholtz_error(
    const TriMesh &mesh,
    const ComplexVector &solution,
    double wavenumber,
    const ComplexFunction &exact,
    const ComplexGradientFunction &exact_gradient,
    const QuadraturePolicy &quadrature,
    const QuadratureContext &quadrature_context) {
    if (solution.size() != static_cast<int>(mesh.nodes.size()))
        throw std::invalid_argument("Helmholtz solution size does not match mesh nodes");
    if (!exact || !exact_gradient)
        throw std::invalid_argument("Helmholtz exact solution callbacks must not be empty");

    const int element_count = static_cast<int>(mesh.elems.size());
    std::vector<double> element_l2_squared(element_count, 0.0);
    std::vector<double> element_gradient_squared(element_count, 0.0);
    std::vector<std::exception_ptr> element_errors(element_count);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if(element_count >= 64)
#endif
    for (int element = 0; element < element_count; ++element) {
        try {
        const Triangle &tri = mesh.elems[element];
        std::array<Eigen::Vector2d, 3> gradients;
        (void)triangle_geometry(mesh, tri, gradients);
        Eigen::Vector2cd discrete_gradient = Eigen::Vector2cd::Zero();
        for (int i = 0; i < 3; ++i)
            discrete_gradient += solution(tri[i]) * gradients[i].cast<Complex>();

        for (const auto &point : triangle_quadrature_points(
                 mesh, element, quadrature, quadrature_context)) {
            Complex discrete_value = 0.0;
            for (int i = 0; i < 3; ++i) {
                discrete_value += point.barycentric[i] * solution(tri[i]);
            }
            const Complex value_error = exact(point.point) - discrete_value;
            const Eigen::Vector2cd gradient_error = exact_gradient(point.point) - discrete_gradient;
            element_l2_squared[element] +=
                point.weight * std::norm(value_error);
            element_gradient_squared[element] +=
                point.weight * gradient_error.squaredNorm();
        }
        } catch (...) {
            element_errors[element] = std::current_exception();
        }
    }

    // Exceptions may not escape an OpenMP work-sharing region.  Store one
    // exception per independent element and rethrow the first element-index
    // failure afterwards, preserving the serial failure order.
    for (const std::exception_ptr &error : element_errors)
        if (error) std::rethrow_exception(error);

    // Reduce in element order so results do not depend on the OpenMP thread
    // count or schedule.
    double l2_squared = 0.0;
    double gradient_squared = 0.0;
    for (int element = 0; element < element_count; ++element) {
        l2_squared += element_l2_squared[element];
        gradient_squared += element_gradient_squared[element];
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

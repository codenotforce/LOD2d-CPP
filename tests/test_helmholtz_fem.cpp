#include "helmholtz/model.h"
#include "helmholtz/manufactured.h"
#include "helmholtz/operators.h"
#include "mesh/refine.h"

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace lod2d;
using namespace lod2d::helmholtz;

namespace {

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

void verify_parallel_exact_error_exception_order() {
    // Use disconnected, well-shaped triangles so each callback can identify
    // its element from the x-coordinate.  Sixty-four elements exercise the
    // OpenMP path when it is available.
    TriMesh mesh;
    constexpr int element_count = 64;
    mesh.nodes.reserve(3 * element_count);
    mesh.elems.reserve(element_count);
    for (int element = 0; element < element_count; ++element) {
        const double x = 10.0 * element;
        const int node = static_cast<int>(mesh.nodes.size());
        mesh.nodes.emplace_back(x, 0.0);
        mesh.nodes.emplace_back(x + 1.0, 0.0);
        mesh.nodes.emplace_back(x, 1.0);
        mesh.elems.push_back({node, node + 1, node + 2});
    }
    const ComplexVector zero = ComplexVector::Zero(mesh.nodes.size());
    const ComplexFunction throwing_exact = [](const Point2 &point) -> Complex {
        const int element = static_cast<int>(std::floor(point.x() / 10.0));
        throw std::runtime_error("exact-element-" + std::to_string(element));
    };
    const ComplexGradientFunction unused_gradient = [](const Point2 &) {
        return Eigen::Vector2cd::Zero();
    };

    try {
        (void)compute_helmholtz_error(
            mesh, zero, 2.0, throwing_exact, unused_gradient);
    } catch (const std::runtime_error &error) {
        require(std::string(error.what()) == "exact-element-0",
                "parallel exact-error evaluation changed exception order");
        return;
    }
    throw std::runtime_error(
        "parallel exact-error evaluation swallowed callback exception");
}

void verify_parallel_load_exception_order() {
    TriMesh mesh;
    constexpr int element_count = 64;
    for (int element = 0; element < element_count; ++element) {
        const double x = 10.0 * element;
        const int node = static_cast<int>(mesh.nodes.size());
        mesh.nodes.emplace_back(x, 0.0);
        mesh.nodes.emplace_back(x + 1.0, 0.0);
        mesh.nodes.emplace_back(x, 1.0);
        mesh.elems.push_back({node, node + 1, node + 2});
    }
    const ComplexFunction throwing_source = [](const Point2 &point) -> Complex {
        const int element = static_cast<int>(std::floor(point.x() / 10.0));
        throw std::runtime_error("load-element-" + std::to_string(element));
    };
    try {
        (void)assemble_helmholtz_load(mesh, throwing_source);
    } catch (const std::runtime_error &error) {
        require(std::string(error.what()) == "load-element-0",
                "parallel load assembly changed exception order");
        return;
    }
    throw std::runtime_error("parallel load assembly swallowed callback exception");
}

} // namespace

int main() {
    try {
        verify_parallel_exact_error_exception_order();
        verify_parallel_load_exception_order();
        const double k = 2.0;
        const TriMesh initial = make_helmholtz_unit_square_mesh();
        const TriMesh boundary_mesh = refine_mesh_nvb(initial, 3).mesh;
        const HelmholtzOperators boundary_operators = assemble_helmholtz_operators(boundary_mesh, k);

        if (helmholtz_fem_solver_available(HelmholtzFemSolverKind::Umfpack)) {
            const ComplexVector probe_load = ComplexVector::Ones(
                boundary_mesh.nodes.size());
            HelmholtzFemSolveTimings sparse_timings;
            HelmholtzFemSolveTimings umfpack_timings;
            const ComplexVector sparse_solution = solve_helmholtz_fem(
                boundary_operators, probe_load,
                HelmholtzFemSolverKind::SparseLu, &sparse_timings);
            const ComplexVector umfpack_solution = solve_helmholtz_fem(
                boundary_operators, probe_load,
                HelmholtzFemSolverKind::Umfpack, &umfpack_timings);
            require(
                (sparse_solution - umfpack_solution).norm()
                    <= 1e-10 * std::max(1.0, sparse_solution.norm()),
                "UMFPACK and SparseLU Helmholtz solutions disagree");
            require(sparse_timings.total_seconds >= 0.0
                        && umfpack_timings.total_seconds >= 0.0,
                    "Helmholtz solver timings are invalid");
        }

        const Eigen::VectorXd ones = Eigen::VectorXd::Ones(boundary_mesh.nodes.size());
        const double boundary_measure = ones.dot(boundary_operators.boundary_mass * ones);
        require(std::abs(boundary_measure - 4.0) < 1e-12,
                "Robin boundary mass does not integrate the unit-square perimeter");

        const Eigen::MatrixXcd dense_system(boundary_operators.system);
        require((dense_system - dense_system.transpose()).norm() < 1e-12,
                "Helmholtz matrix must be complex symmetric for real coefficients");
        require((dense_system - dense_system.adjoint()).norm() > 1e-3,
                "Helmholtz matrix was accidentally assembled as Hermitian");
        const double manufactured_k = 4.0;
        const HelmholtzManufacturedSolution manufactured =
            make_polynomial_plane_wave_solution(manufactured_k);
        for (double t : {0.0, 0.37, 1.0}) {
            const std::array<std::pair<Point2, Eigen::Vector2d>, 4> boundary_points{{
                {Point2(0.0, t), Eigen::Vector2d(-1.0, 0.0)},
                {Point2(1.0, t), Eigen::Vector2d(1.0, 0.0)},
                {Point2(t, 0.0), Eigen::Vector2d(0.0, -1.0)},
                {Point2(t, 1.0), Eigen::Vector2d(0.0, 1.0)}}};
            for (const auto &[point, normal] : boundary_points) {
                const Complex value = manufactured.value(point);
                const Complex normal_derivative =
                    normal.cast<Complex>().dot(manufactured.gradient(point));
                require(std::abs(normal_derivative
                            - Complex(0.0, manufactured_k) * value) < 1e-12,
                        "polynomial plane wave violates homogeneous Robin data");
            }
        }
        std::vector<double> manufactured_errors;
        for (int level : {4, 6, 8}) {
            const TriMesh mesh = refine_mesh_nvb(initial, level).mesh;
            const HelmholtzOperators operators =
                assemble_helmholtz_operators(mesh, manufactured_k);
            const ComplexVector load =
                assemble_helmholtz_load(mesh, manufactured.source);
            const ComplexVector solution = solve_helmholtz_fem(operators, load);
            manufactured_errors.push_back(compute_helmholtz_error(
                mesh, solution, manufactured_k,
                manufactured.value, manufactured.gradient).energy);
        }
        require(manufactured_errors[2] < 0.4 * manufactured_errors[0],
                "polynomial plane-wave energy error did not converge");


        const Complex b = (k * k + Complex(0.0, 2.0 * k)) / Complex(2.0, -k);
        auto g = [=](double x) { return Complex(1.0, 0.0) - Complex(0.0, k) * x + b * x * x; };
        auto gp = [=](double x) { return -Complex(0.0, k) + 2.0 * b * x; };
        auto exact = [=](const Point2 &point) { return g(point.x()) * g(point.y()); };
        auto exact_gradient = [=](const Point2 &point) {
            Eigen::Vector2cd gradient;
            gradient << gp(point.x()) * g(point.y()), g(point.x()) * gp(point.y());
            return gradient;
        };
        auto source = [=](const Point2 &point) {
            return -2.0 * b * (g(point.x()) + g(point.y()))
                 - k * k * g(point.x()) * g(point.y());
        };

        const Complex left_robin = -gp(0.0) - Complex(0.0, k) * g(0.0);
        const Complex right_robin = gp(1.0) - Complex(0.0, k) * g(1.0);
        require(std::abs(left_robin) < 1e-12 && std::abs(right_robin) < 1e-12,
                "manufactured solution does not satisfy homogeneous Robin data");

        std::vector<double> diameters;
        std::vector<double> energy_errors;
        std::vector<double> l2_errors;
        for (int level : {3, 5, 7}) {
            const TriMesh mesh = refine_mesh_nvb(initial, level).mesh;
            const HelmholtzOperators operators = assemble_helmholtz_operators(mesh, k);
            const ComplexVector load = assemble_helmholtz_load(mesh, source);
            const ComplexVector solution = solve_helmholtz_fem(operators, load);
            const HelmholtzError error = compute_helmholtz_error(
                mesh, solution, k, exact, exact_gradient);
            diameters.push_back(max_element_diameter(mesh));
            energy_errors.push_back(error.energy);
            l2_errors.push_back(error.l2);
        }

        require(energy_errors[2] < 0.45 * energy_errors[0],
                "Helmholtz energy error did not decrease under global NVB");
        require(l2_errors[2] < 0.25 * l2_errors[0],
                "Helmholtz L2 error did not decrease under global NVB");
        const double energy_rate = std::log(energy_errors[0] / energy_errors[2])
                                 / std::log(diameters[0] / diameters[2]);
        const double l2_rate = std::log(l2_errors[0] / l2_errors[2])
                             / std::log(diameters[0] / diameters[2]);
        require(energy_rate > 0.65, "Helmholtz energy convergence rate is too low");
        require(l2_rate > 1.35, "Helmholtz L2 convergence rate is too low");

        std::cout << "Helmholtz FEM: energy rate=" << energy_rate
                  << " L2 rate=" << l2_rate << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_fem failed: " << error.what() << '\n';
        return 1;
    }
}

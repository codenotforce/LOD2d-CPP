#include "solver/right_gmres.h"

#include <Eigen/LU>
#include <algorithm>
#include <complex>
#include <iostream>
#include <stdexcept>

namespace {

using lod2d::solver::ComplexVector;
using lod2d::solver::RightGmresConfig;
using lod2d::solver::RightGmresResult;

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

Eigen::MatrixXcd nonnormal_matrix() {
    using Complex = std::complex<double>;
    Eigen::MatrixXcd matrix = Eigen::MatrixXcd::Zero(5, 5);
    matrix.diagonal() << Complex(4.0, 0.2), Complex(3.0, -0.3),
        Complex(2.5, 0.7), Complex(1.8, -0.4), Complex(1.2, 0.5);
    matrix(0, 1) = Complex(2.0, -0.5);
    matrix(1, 2) = Complex(-1.5, 0.8);
    matrix(2, 3) = Complex(1.0, 1.2);
    matrix(3, 4) = Complex(0.7, -0.9);
    matrix(0, 4) = Complex(-0.4, 0.3);
    return matrix;
}

} // namespace

int main() {
    try {
        const Eigen::MatrixXcd matrix = nonnormal_matrix();
        ComplexVector exact(5);
        exact << std::complex<double>(1.0, 0.2),
            std::complex<double>(-0.5, 0.8),
            std::complex<double>(0.3, -1.1),
            std::complex<double>(1.4, 0.0),
            std::complex<double>(-0.7, 0.4);
        const ComplexVector rhs = matrix * exact;
        const auto apply_matrix = [&](const ComplexVector &value) {
            return ComplexVector(matrix * value);
        };
        const auto identity = [](const ComplexVector &value) { return value; };

        RightGmresConfig config;
        config.restart = 4;
        config.max_iterations = 100;
        config.relative_tolerance = 1e-10;
        const RightGmresResult restarted =
            lod2d::solver::solve_right_preconditioned_gmres(
                5, apply_matrix, identity, rhs, config);
        require(restarted.converged, "restarted right GMRES did not converge");
        require(
            (restarted.solution - exact).norm() / exact.norm() < 1e-8,
            "restarted right GMRES returned the wrong solution");
        require(
            (rhs - matrix * restarted.solution).norm() / rhs.norm() < 1e-10,
            "restarted right GMRES true residual is too large");

        Eigen::FullPivLU<Eigen::MatrixXcd> exact_factorization(matrix);
        const auto exact_right_preconditioner = [&](const ComplexVector &value) {
            return ComplexVector(exact_factorization.solve(value));
        };
        config.restart = 5;
        config.max_iterations = 5;
        config.relative_tolerance = 1e-12;
        const RightGmresResult exact_preconditioned =
            lod2d::solver::solve_right_preconditioned_gmres(
                5, apply_matrix, exact_right_preconditioner, rhs, config);
        require(exact_preconditioned.converged, "exactly preconditioned GMRES failed");
        require(
            exact_preconditioned.iterations == 1,
            "exact right preconditioner should converge in one iteration");
        require(
            exact_preconditioned.relative_residual < 1e-12,
            "exactly preconditioned true residual is too large");

        const RightGmresResult zero =
            lod2d::solver::solve_right_preconditioned_gmres(
                5, apply_matrix, identity, ComplexVector::Zero(5), config);
        require(zero.converged && zero.iterations == 0 && zero.solution.norm() == 0.0,
            "zero right-hand side handling failed");

        int flexible_calls = 0;
        const auto variable_preconditioner =
            [&](const ComplexVector &value) {
                ++flexible_calls;
                ComplexVector result = value;
                for (int index = 0; index < result.size(); ++index) {
                    const double scale = flexible_calls % 2 == 0
                        ? 0.7 + 0.1 * index
                        : 1.3 - 0.1 * index;
                    result(index) *= scale;
                }
                return result;
            };
        config.restart = 5;
        config.max_iterations = 30;
        config.relative_tolerance = 1e-10;
        const RightGmresResult flexible =
            lod2d::solver::solve_right_preconditioned_fgmres(
                5, apply_matrix, variable_preconditioner, rhs, config);
        require(flexible.converged,
                "right-preconditioned FGMRES did not converge");
        require(flexible_calls == flexible.iterations,
                "FGMRES did not apply one variable preconditioner per iteration");
        require(
            (rhs - matrix * flexible.solution).norm() / rhs.norm() < 1e-10,
            "FGMRES true residual is too large");

        std::cout << "Right-preconditioned GMRES/FGMRES tests passed.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_right_gmres failed: " << error.what() << '\n';
        return 1;
    }
}

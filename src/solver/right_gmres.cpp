#include "solver/right_gmres.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace lod2d::solver {
namespace {

void validate_vector(const ComplexVector &vector, int dimension, const char *name) {
    if (vector.size() != dimension)
        throw std::runtime_error(std::string(name) + " returned the wrong vector size");
    if (!vector.allFinite())
        throw std::runtime_error(std::string(name) + " returned non-finite values");
}

} // namespace

RightGmresResult solve_right_preconditioned_gmres(
    int dimension,
    const ComplexLinearOperator &apply_operator,
    const ComplexLinearOperator &apply_preconditioner,
    const ComplexVector &right_hand_side,
    const RightGmresConfig &config) {
    if (dimension <= 0) throw std::invalid_argument("GMRES dimension must be positive");
    if (right_hand_side.size() != dimension)
        throw std::invalid_argument("GMRES right-hand side has the wrong size");
    if (!right_hand_side.allFinite())
        throw std::invalid_argument("GMRES right-hand side contains non-finite values");
    if (!apply_operator || !apply_preconditioner)
        throw std::invalid_argument("GMRES operators must not be empty");
    if (config.restart <= 0 || config.max_iterations <= 0)
        throw std::invalid_argument("GMRES iteration limits must be positive");
    if (!(config.relative_tolerance >= 0.0)
        || !(config.absolute_tolerance >= 0.0)
        || !std::isfinite(config.relative_tolerance)
        || !std::isfinite(config.absolute_tolerance)) {
        throw std::invalid_argument("GMRES tolerances must be finite and nonnegative");
    }

    RightGmresResult result;
    result.solution = ComplexVector::Zero(dimension);
    const double rhs_norm = right_hand_side.norm();
    const double denominator = std::max(rhs_norm, 1e-30);
    const double target = std::max(
        config.absolute_tolerance, config.relative_tolerance * rhs_norm);
    if (rhs_norm == 0.0) {
        result.converged = true;
        result.relative_residual = 0.0;
        result.message = "zero right-hand side";
        return result;
    }

    const int restart = std::min(config.restart, dimension);
    while (result.iterations < config.max_iterations) {
        ComplexVector residual = right_hand_side - apply_operator(result.solution);
        validate_vector(residual, dimension, "GMRES operator");
        const double beta = residual.norm();
        result.relative_residual = beta / denominator;
        if (beta <= target) {
            result.converged = true;
            result.message = "converged before restart";
            return result;
        }

        const int cycle_iterations = std::min(
            restart, config.max_iterations - result.iterations);
        Eigen::MatrixXcd basis = Eigen::MatrixXcd::Zero(
            dimension, cycle_iterations + 1);
        Eigen::MatrixXcd preconditioned = Eigen::MatrixXcd::Zero(
            dimension, cycle_iterations);
        Eigen::MatrixXcd hessenberg = Eigen::MatrixXcd::Zero(
            cycle_iterations + 1, cycle_iterations);
        basis.col(0) = residual / beta;

        ComplexVector best_solution = result.solution;
        double best_residual = beta;
        bool breakdown = false;
        for (int j = 0; j < cycle_iterations; ++j) {
            ComplexVector z = apply_preconditioner(basis.col(j));
            validate_vector(z, dimension, "GMRES preconditioner");
            preconditioned.col(j) = z;

            ComplexVector work = apply_operator(z);
            validate_vector(work, dimension, "GMRES operator");
            for (int i = 0; i <= j; ++i) {
                const std::complex<double> coefficient = basis.col(i).dot(work);
                hessenberg(i, j) += coefficient;
                work.noalias() -= coefficient * basis.col(i);
            }
            if (config.reorthogonalize) {
                for (int i = 0; i <= j; ++i) {
                    const std::complex<double> correction = basis.col(i).dot(work);
                    hessenberg(i, j) += correction;
                    work.noalias() -= correction * basis.col(i);
                }
            }

            const double next_norm = work.norm();
            hessenberg(j + 1, j) = next_norm;
            const double breakdown_scale = std::max(1.0, hessenberg.col(j).norm());
            breakdown = next_norm <= 64.0 * std::numeric_limits<double>::epsilon()
                * breakdown_scale;
            if (!breakdown) basis.col(j + 1) = work / next_norm;

            Eigen::VectorXcd least_squares_rhs =
                Eigen::VectorXcd::Zero(j + 2);
            least_squares_rhs(0) = beta;
            const Eigen::MatrixXcd small_hessenberg =
                hessenberg.block(0, 0, j + 2, j + 1);
            const Eigen::VectorXcd coefficients =
                small_hessenberg.colPivHouseholderQr().solve(least_squares_rhs);
            ComplexVector candidate = result.solution
                + preconditioned.leftCols(j + 1) * coefficients;
            validate_vector(candidate, dimension, "GMRES update");

            ComplexVector true_residual = right_hand_side - apply_operator(candidate);
            validate_vector(true_residual, dimension, "GMRES operator");
            const double true_norm = true_residual.norm();
            ++result.iterations;
            if (true_norm < best_residual) {
                best_residual = true_norm;
                best_solution = std::move(candidate);
            }
            result.relative_residual = true_norm / denominator;
            if (true_norm <= target) {
                result.solution = std::move(best_solution);
                result.converged = true;
                result.happy_breakdown = breakdown;
                result.message = breakdown ? "converged at happy breakdown" : "converged";
                return result;
            }
            if (breakdown) break;
        }

        result.solution = std::move(best_solution);
        ++result.restarts;
        if (breakdown) {
            result.happy_breakdown = true;
            result.message = "GMRES breakdown before reaching the true-residual tolerance";
            return result;
        }
    }

    const ComplexVector final_residual =
        right_hand_side - apply_operator(result.solution);
    validate_vector(final_residual, dimension, "GMRES operator");
    result.relative_residual = final_residual.norm() / denominator;
    result.message = "GMRES reached the maximum iteration count";
    return result;
}

RightGmresResult solve_right_preconditioned_fgmres(
    int dimension,
    const ComplexLinearOperator &apply_operator,
    const ComplexLinearOperator &apply_preconditioner,
    const ComplexVector &right_hand_side,
    const RightGmresConfig &config) {
    return solve_right_preconditioned_gmres(
        dimension,
        apply_operator,
        apply_preconditioner,
        right_hand_side,
        config);
}

} // namespace lod2d::solver

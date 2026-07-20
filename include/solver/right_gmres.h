#pragma once

#include <Eigen/Dense>
#include <functional>
#include <string>

namespace lod2d::solver {

using ComplexVector = Eigen::VectorXcd;
using ComplexLinearOperator =
    std::function<ComplexVector(const ComplexVector &)>;

struct RightGmresConfig {
    int restart = 30;
    int max_iterations = 200;
    double relative_tolerance = 1e-12;
    double absolute_tolerance = 0.0;
    bool reorthogonalize = true;
};

struct RightGmresResult {
    ComplexVector solution;
    int iterations = 0;
    int restarts = 0;
    double relative_residual = 0.0;
    bool converged = false;
    bool happy_breakdown = false;
    std::string message;
};

RightGmresResult solve_right_preconditioned_gmres(
    int dimension,
    const ComplexLinearOperator &apply_operator,
    const ComplexLinearOperator &apply_preconditioner,
    const ComplexVector &right_hand_side,
    const RightGmresConfig &config = {});

RightGmresResult solve_right_preconditioned_fgmres(
    int dimension,
    const ComplexLinearOperator &apply_operator,
    const ComplexLinearOperator &apply_preconditioner,
    const ComplexVector &right_hand_side,
    const RightGmresConfig &config = {});

} // namespace lod2d::solver

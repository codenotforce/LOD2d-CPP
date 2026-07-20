#include "helmholtz/patch_multigrid.h"

#include "helmholtz/shifted_laplacian.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace lod2d::helmholtz {

HelmholtzPatchVcycle::HelmholtzPatchVcycle(
    const HelmholtzPatchSystem &system,
    double epsilon,
    int pre_smooth,
    int post_smooth,
    int coarse_max_dofs,
    double jacobi_weight)
    : HelmholtzPatchVcycle(
          build_shifted_helmholtz_operator(system, epsilon),
          system.geometric_prolongations,
          pre_smooth,
          post_smooth,
          coarse_max_dofs,
          jacobi_weight) {}

HelmholtzPatchVcycle::HelmholtzPatchVcycle(
    const ComplexSparseMatrix &shifted_operator,
    const std::vector<Eigen::SparseMatrix<double>> &geometric_prolongations,
    int pre_smooth,
    int post_smooth,
    int coarse_max_dofs,
    double jacobi_weight)
    : pre_smooth_(pre_smooth),
      post_smooth_(post_smooth),
      jacobi_weight_(jacobi_weight) {
    if (shifted_operator.rows() <= 0
        || shifted_operator.rows() != shifted_operator.cols())
        throw std::invalid_argument(
            "V-cycle shifted operator must be square and nonempty");
    if (geometric_prolongations.empty())
        throw std::invalid_argument(
            "V-cycle requires at least one geometric prolongation");
    if (pre_smooth_ < 0 || post_smooth_ < 0)
        throw std::invalid_argument("V-cycle smoothing counts must be nonnegative");
    if (coarse_max_dofs <= 0)
        throw std::invalid_argument("V-cycle coarse maximum DOFs must be positive");
    if (!(jacobi_weight_ > 0.0) || !std::isfinite(jacobi_weight_))
        throw std::invalid_argument("V-cycle Jacobi weight must be finite and positive");

    std::vector<ComplexSparseMatrix> all_prolongations;
    all_prolongations.reserve(geometric_prolongations.size());
    int expected_rows = shifted_operator.rows();
    for (int level = static_cast<int>(geometric_prolongations.size()) - 1;
         level >= 0; --level) {
        const auto &prolongation = geometric_prolongations[level];
        if (prolongation.rows() != expected_rows
            || prolongation.cols() <= 0)
            throw std::invalid_argument(
                "V-cycle geometric prolongation dimensions are inconsistent");
        expected_rows = prolongation.cols();
    }
    for (const auto &prolongation : geometric_prolongations) {
        ComplexSparseMatrix complex_prolongation =
            prolongation.cast<Complex>();
        complex_prolongation.makeCompressed();
        all_prolongations.push_back(std::move(complex_prolongation));
    }

    std::vector<ComplexSparseMatrix> all_operators(
        all_prolongations.size() + 1);
    all_operators.back() = shifted_operator;
    all_operators.back().makeCompressed();
    for (int level = static_cast<int>(all_prolongations.size()) - 1;
         level >= 0; --level) {
        all_operators[level] =
            all_prolongations[level].adjoint()
            * all_operators[level + 1]
            * all_prolongations[level];
        all_operators[level].prune(Complex(0.0, 0.0), 1e-14);
        all_operators[level].makeCompressed();
    }

    int first_level = 0;
    for (int level = 0;
         level + 1 < static_cast<int>(all_operators.size()); ++level) {
        if (all_operators[level].rows() <= coarse_max_dofs)
            first_level = level;
    }
    operators_.assign(
        all_operators.begin() + first_level, all_operators.end());
    prolongations_.assign(
        all_prolongations.begin() + first_level,
        all_prolongations.end());

    inverse_diagonals_.reserve(operators_.size());
    for (const ComplexSparseMatrix &op : operators_) {
        ComplexVector inverse(op.rows());
        for (int row = 0; row < op.rows(); ++row) {
            const Complex diagonal = op.coeff(row, row);
            if (std::abs(diagonal) <= 1e-14)
                throw std::runtime_error(
                    "V-cycle shifted operator has a zero diagonal");
            inverse(row) = Complex(1.0, 0.0) / diagonal;
        }
        inverse_diagonals_.push_back(std::move(inverse));
    }

    coarse_solver_.analyzePattern(operators_.front());
    coarse_solver_.factorize(operators_.front());
    if (coarse_solver_.info() != Eigen::Success)
        throw std::runtime_error(
            "V-cycle coarse shifted Helmholtz factorization failed");
}

void HelmholtzPatchVcycle::smooth(
    int level,
    const ComplexVector &right_hand_side,
    int steps,
    ComplexVector &iterate) const {
    for (int step = 0; step < steps; ++step) {
        const ComplexVector residual =
            right_hand_side - operators_[level] * iterate;
        iterate.noalias() += jacobi_weight_
            * inverse_diagonals_[level].cwiseProduct(residual);
    }
}

ComplexVector HelmholtzPatchVcycle::cycle(
    int level,
    const ComplexVector &right_hand_side) const {
    if (level == 0) {
        ComplexVector solution = coarse_solver_.solve(right_hand_side);
        if (!solution.allFinite())
            throw std::runtime_error(
                "V-cycle coarse solve returned non-finite values");
        return solution;
    }

    ComplexVector iterate = ComplexVector::Zero(right_hand_side.size());
    smooth(level, right_hand_side, pre_smooth_, iterate);
    const ComplexVector residual =
        right_hand_side - operators_[level] * iterate;
    const ComplexVector coarse_rhs =
        prolongations_[level - 1].adjoint() * residual;
    iterate.noalias() +=
        prolongations_[level - 1] * cycle(level - 1, coarse_rhs);
    smooth(level, right_hand_side, post_smooth_, iterate);
    return iterate;
}

ComplexVector HelmholtzPatchVcycle::apply(
    const ComplexVector &right_hand_side) const {
    if (right_hand_side.size() != finest_dofs())
        throw std::invalid_argument(
            "V-cycle right-hand side has the wrong size");
    if (!right_hand_side.allFinite())
        throw std::invalid_argument(
            "V-cycle right-hand side contains non-finite values");
    return cycle(levels() - 1, right_hand_side);
}

double HelmholtzPatchVcycle::relative_residual(
    const ComplexVector &right_hand_side) const {
    const double scale = std::max(1e-30, right_hand_side.norm());
    const ComplexVector approximation = apply(right_hand_side);
    return (right_hand_side - operators_.back() * approximation).norm()
        / scale;
}

} // namespace lod2d::helmholtz

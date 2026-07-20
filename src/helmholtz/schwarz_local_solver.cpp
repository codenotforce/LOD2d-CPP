#include "helmholtz/schwarz_local_solver.h"

#include "helmholtz/patch_multigrid.h"

#include <Eigen/SparseLU>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace lod2d::helmholtz {
namespace {

const std::vector<Eigen::SparseMatrix<double>> &empty_prolongations() {
    static const std::vector<Eigen::SparseMatrix<double>> value;
    return value;
}

} // namespace

struct HelmholtzSchwarzLocalSolver::Impl {
    HelmholtzSchwarzLocalSolverConfig config;
    double epsilon = 0.0;
    int dimension = 0;
    ComplexSparseMatrix matrix;
    Eigen::SparseLU<ComplexSparseMatrix> factorization;
    std::unique_ptr<HelmholtzPatchVcycle> vcycle;
};

HelmholtzSchwarzLocalSolver::HelmholtzSchwarzLocalSolver(
    const ComplexSparseMatrix &matrix,
    const Eigen::SparseMatrix<double> &mass,
    double wavenumber,
    HelmholtzSchwarzLocalSolverConfig config)
    : HelmholtzSchwarzLocalSolver(
          matrix,
          mass,
          empty_prolongations(),
          wavenumber,
          config) {}

HelmholtzSchwarzLocalSolver::HelmholtzSchwarzLocalSolver(
    const ComplexSparseMatrix &matrix,
    const Eigen::SparseMatrix<double> &mass,
    const std::vector<Eigen::SparseMatrix<double>> &geometric_prolongations,
    double wavenumber,
    HelmholtzSchwarzLocalSolverConfig config)
    : impl_(std::make_unique<Impl>()) {
    if (matrix.rows() <= 0 || matrix.rows() != matrix.cols())
        throw std::invalid_argument(
            "Schwarz local Helmholtz matrix must be square and nonempty");
    if (config.kind == HelmholtzSchwarzLocalSolverKind::ShiftedGmres
        && (mass.rows() != matrix.rows() || mass.cols() != matrix.cols()))
        throw std::invalid_argument(
            "Schwarz local mass matrix has the wrong size");
    if (!(wavenumber > 0.0) || !std::isfinite(wavenumber))
        throw std::invalid_argument(
            "Schwarz local wavenumber must be finite and positive");
    if (!(config.shift_alpha >= 0.0) || !std::isfinite(config.shift_alpha))
        throw std::invalid_argument(
            "Schwarz local shift alpha must be finite and nonnegative");
    if (config.kind != HelmholtzSchwarzLocalSolverKind::SparseLu
        && config.kind != HelmholtzSchwarzLocalSolverKind::ShiftedGmres)
        throw std::invalid_argument("unknown Schwarz local solver kind");
    if (config.shifted_inverse != HelmholtzSchwarzShiftedInverseKind::SparseLu
        && config.shifted_inverse
            != HelmholtzSchwarzShiftedInverseKind::GeometricVcycle)
        throw std::invalid_argument("unknown Schwarz shifted inverse kind");
    if (config.kind == HelmholtzSchwarzLocalSolverKind::ShiftedGmres) {
        if (config.gmres.restart <= 0 || config.gmres.max_iterations <= 0)
            throw std::invalid_argument(
                "Schwarz local GMRES iteration limits must be positive");
        if (!(config.gmres.relative_tolerance >= 0.0)
            || !(config.gmres.absolute_tolerance >= 0.0)
            || !std::isfinite(config.gmres.relative_tolerance)
            || !std::isfinite(config.gmres.absolute_tolerance))
            throw std::invalid_argument(
                "Schwarz local GMRES tolerances must be finite and nonnegative");
        if (config.shifted_inverse
                == HelmholtzSchwarzShiftedInverseKind::GeometricVcycle
            && geometric_prolongations.empty())
            throw std::invalid_argument(
                "Schwarz geometric V-cycle requires patch prolongations");
    }

    impl_->config = config;
    impl_->dimension = matrix.rows();
    if (config.kind == HelmholtzSchwarzLocalSolverKind::SparseLu) {
        ComplexSparseMatrix factor_matrix = matrix;
        factor_matrix.makeCompressed();
        impl_->factorization.analyzePattern(factor_matrix);
        impl_->factorization.factorize(factor_matrix);
    } else {
        impl_->epsilon = config.shift_alpha * wavenumber * wavenumber;
        impl_->matrix = matrix;
        ComplexSparseMatrix shifted = matrix
            - Complex(0.0, impl_->epsilon) * mass.cast<Complex>();
        shifted.makeCompressed();
        if (config.shifted_inverse
            == HelmholtzSchwarzShiftedInverseKind::SparseLu) {
            impl_->factorization.analyzePattern(shifted);
            impl_->factorization.factorize(shifted);
        } else {
            impl_->vcycle = std::make_unique<HelmholtzPatchVcycle>(
                shifted,
                geometric_prolongations,
                config.vcycle_pre_smooth,
                config.vcycle_post_smooth,
                config.vcycle_coarse_max_dofs,
                config.vcycle_jacobi_weight);
        }
    }
    if (!impl_->vcycle && impl_->factorization.info() != Eigen::Success)
        throw std::runtime_error(
            "Schwarz local Helmholtz factorization failed");
}

HelmholtzSchwarzLocalSolver::~HelmholtzSchwarzLocalSolver() = default;
HelmholtzSchwarzLocalSolver::HelmholtzSchwarzLocalSolver(
    HelmholtzSchwarzLocalSolver &&) noexcept = default;
HelmholtzSchwarzLocalSolver &HelmholtzSchwarzLocalSolver::operator=(
    HelmholtzSchwarzLocalSolver &&) noexcept = default;

HelmholtzSchwarzLocalSolveResult HelmholtzSchwarzLocalSolver::solve(
    const ComplexVector &right_hand_side) const {
    if (right_hand_side.size() != impl_->dimension)
        throw std::invalid_argument(
            "Schwarz local right-hand side has the wrong size");
    if (!right_hand_side.allFinite())
        throw std::invalid_argument(
            "Schwarz local right-hand side is non-finite");

    HelmholtzSchwarzLocalSolveResult result;
    if (impl_->config.kind == HelmholtzSchwarzLocalSolverKind::SparseLu) {
        result.solution = impl_->factorization.solve(right_hand_side);
        if (impl_->factorization.info() != Eigen::Success
            || !result.solution.allFinite())
            throw std::runtime_error(
                "Schwarz local SparseLU solve failed");
        return result;
    }

    const auto apply_operator = [&](const solver::ComplexVector &vector) {
        return solver::ComplexVector(impl_->matrix * vector);
    };
    const auto apply_preconditioner =
        [&](const solver::ComplexVector &vector) {
            if (impl_->vcycle)
                return solver::ComplexVector(impl_->vcycle->apply(vector));
            solver::ComplexVector value =
                impl_->factorization.solve(vector);
            if (impl_->factorization.info() != Eigen::Success
                || !value.allFinite())
                throw std::runtime_error(
                    "Schwarz shifted preconditioner solve failed");
            return value;
        };
    const solver::RightGmresResult gmres =
        solver::solve_right_preconditioned_gmres(
            impl_->dimension,
            apply_operator,
            apply_preconditioner,
            right_hand_side,
            impl_->config.gmres);
    if (!gmres.converged)
        throw std::runtime_error(
            "Schwarz local shifted GMRES failed: " + gmres.message);
    result.solution = gmres.solution;
    result.iterations = gmres.iterations;
    result.restarts = gmres.restarts;
    result.relative_residual = gmres.relative_residual;
    return result;
}

ComplexMatrix HelmholtzSchwarzLocalSolver::solve_block(
    const ComplexMatrix &right_hand_sides) const {
    if (impl_->config.kind != HelmholtzSchwarzLocalSolverKind::SparseLu)
        throw std::logic_error(
            "Schwarz block solve is available only for SparseLU");
    if (right_hand_sides.rows() != impl_->dimension)
        throw std::invalid_argument(
            "Schwarz local block right-hand side has the wrong size");
    if (!right_hand_sides.allFinite())
        throw std::invalid_argument(
            "Schwarz local block right-hand side is non-finite");
    if (right_hand_sides.cols() == 0)
        return ComplexMatrix(impl_->dimension, 0);

    ComplexMatrix solutions = impl_->factorization.solve(right_hand_sides);
    if (impl_->factorization.info() != Eigen::Success
        || !solutions.allFinite())
        throw std::runtime_error("Schwarz local SparseLU block solve failed");
    return solutions;
}

const HelmholtzSchwarzLocalSolverConfig &
HelmholtzSchwarzLocalSolver::config() const {
    return impl_->config;
}

double HelmholtzSchwarzLocalSolver::epsilon() const {
    return impl_->epsilon;
}

int HelmholtzSchwarzLocalSolver::vcycle_levels() const {
    return impl_->vcycle ? impl_->vcycle->levels() : 0;
}

int HelmholtzSchwarzLocalSolver::vcycle_coarse_dofs() const {
    return impl_->vcycle ? impl_->vcycle->coarse_dofs() : 0;
}

int HelmholtzSchwarzLocalSolver::vcycle_finest_dofs() const {
    return impl_->vcycle ? impl_->vcycle->finest_dofs() : 0;
}

} // namespace lod2d::helmholtz

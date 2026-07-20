#pragma once

#include "helmholtz/types.h"
#include "solver/right_gmres.h"

#include <Eigen/Sparse>
#include <memory>
#include <vector>

namespace lod2d::helmholtz {

enum class HelmholtzSchwarzLocalSolverKind {
    SparseLu,
    ShiftedGmres
};

enum class HelmholtzSchwarzShiftedInverseKind {
    SparseLu,
    GeometricVcycle
};

struct HelmholtzSchwarzLocalSolverConfig {
    HelmholtzSchwarzLocalSolverKind kind =
        HelmholtzSchwarzLocalSolverKind::SparseLu;
    HelmholtzSchwarzShiftedInverseKind shifted_inverse =
        HelmholtzSchwarzShiftedInverseKind::SparseLu;
    double shift_alpha = 0.2;
    int vcycle_pre_smooth = 2;
    int vcycle_post_smooth = 2;
    int vcycle_coarse_max_dofs = 256;
    double vcycle_jacobi_weight = 0.6;
    solver::RightGmresConfig gmres;
};

struct HelmholtzSchwarzLocalSolveResult {
    ComplexVector solution;
    int iterations = 0;
    int restarts = 0;
    double relative_residual = 0.0;
};

class HelmholtzSchwarzLocalSolver {
public:
    HelmholtzSchwarzLocalSolver(
        const ComplexSparseMatrix &matrix,
        const Eigen::SparseMatrix<double> &mass,
        double wavenumber,
        HelmholtzSchwarzLocalSolverConfig config = {});
    HelmholtzSchwarzLocalSolver(
        const ComplexSparseMatrix &matrix,
        const Eigen::SparseMatrix<double> &mass,
        const std::vector<Eigen::SparseMatrix<double>> &geometric_prolongations,
        double wavenumber,
        HelmholtzSchwarzLocalSolverConfig config = {});
    ~HelmholtzSchwarzLocalSolver();

    HelmholtzSchwarzLocalSolver(HelmholtzSchwarzLocalSolver &&) noexcept;
    HelmholtzSchwarzLocalSolver &operator=(
        HelmholtzSchwarzLocalSolver &&) noexcept;
    HelmholtzSchwarzLocalSolver(
        const HelmholtzSchwarzLocalSolver &) = delete;
    HelmholtzSchwarzLocalSolver &operator=(
        const HelmholtzSchwarzLocalSolver &) = delete;

    HelmholtzSchwarzLocalSolveResult solve(
        const ComplexVector &right_hand_side) const;
    ComplexMatrix solve_block(const ComplexMatrix &right_hand_sides) const;

    const HelmholtzSchwarzLocalSolverConfig &config() const;
    double epsilon() const;
    int vcycle_levels() const;
    int vcycle_coarse_dofs() const;
    int vcycle_finest_dofs() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace lod2d::helmholtz

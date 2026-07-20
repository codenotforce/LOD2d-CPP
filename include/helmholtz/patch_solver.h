#pragma once

#include "helmholtz/patch_system.h"
#include "solver/right_gmres.h"

namespace lod2d::helmholtz {

enum class HelmholtzPatchSolverKind {
    DirectSaddle,
    DirectSchur,
    ShiftedGmres
};

enum class HelmholtzShiftRule {
    KappaSquared,
    PatchScaled,
    Absolute
};

enum class HelmholtzShiftedInverseKind {
    Identity,
    SparseLu,
    GeometricVcycle
};

struct HelmholtzShiftedLaplacianConfig {
    HelmholtzShiftRule rule = HelmholtzShiftRule::KappaSquared;
    double alpha = 0.2;
    double absolute_epsilon = 0.0;
    HelmholtzShiftedInverseKind inverse =
        HelmholtzShiftedInverseKind::SparseLu;
    int pre_smooth = 2;
    int post_smooth = 2;
    int coarse_max_dofs = 200;
    double jacobi_weight = 0.6;
};

struct HelmholtzPatchSolverConfig {
    HelmholtzPatchSolverKind kind =
        HelmholtzPatchSolverKind::DirectSaddle;
    solver::RightGmresConfig gmres;
    HelmholtzShiftedLaplacianConfig shifted;
    bool fallback_to_direct = false;
};

struct HelmholtzPatchSolveDiagnostics {
    double primal_residual = 0.0;
    double adjoint_residual = 0.0;
    double constraint_residual = 0.0;
    double schur_residual = 0.0;
    double schur_rcond = 1.0;
    double max_gmres_residual = 0.0;
    double vcycle_relative_residual = 0.0;
    int vcycle_levels = 0;
    int vcycle_coarse_dofs = 0;
    int vcycle_finest_dofs = 0;
    int gmres_right_hand_sides = 0;
    int gmres_total_iterations = 0;
    int gmres_max_iterations = 0;
    int gmres_restarts = 0;
    bool symbolic_reused = false;
    bool direct_fallback = false;
};

struct HelmholtzPatchSolveResult {
    ComplexMatrix corrector;
    ComplexMatrix multipliers;
    HelmholtzPatchSolveDiagnostics diagnostics;
};

double helmholtz_shift_epsilon(
    const HelmholtzPatchSystem &system,
    const HelmholtzShiftedLaplacianConfig &config);

HelmholtzPatchSolveResult solve_helmholtz_patch(
    const HelmholtzPatchSystem &system,
    const HelmholtzPatchSolverConfig &config = {});

} // namespace lod2d::helmholtz

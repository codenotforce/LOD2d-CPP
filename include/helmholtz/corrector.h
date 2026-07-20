#pragma once

#include "helmholtz/operators.h"
#include "helmholtz/patch_solver.h"
#include <Eigen/Sparse>
#include <vector>

namespace lod2d::helmholtz {

struct HelmholtzCorrectorEntry {
    int row = -1;
    int local_coarse_vertex = -1;
    Complex value = 0.0;
};

using HelmholtzElementCorrector = std::vector<HelmholtzCorrectorEntry>;

struct HelmholtzCorrectorDiagnostics {
    double max_primal_residual = 0.0;
    double max_adjoint_residual = 0.0;
    double max_constraint_residual = 0.0;
    int patch_count = 0;
    int patches_touching_physical_boundary = 0;
    int parallel_threads = 1;
    int symbolic_analyses = 0;
    int symbolic_reuses = 0;
    int direct_fallbacks = 0;
    int gmres_right_hand_sides = 0;
    int gmres_iterations = 0;
    double max_vcycle_relative_residual = 0.0;
    int max_vcycle_levels = 0;
    int max_vcycle_coarse_dofs = 0;
    int max_vcycle_finest_dofs = 0;
    int gmres_max_iterations = 0;
    int gmres_restarts = 0;
    double max_gmres_relative_residual = 0.0;
    double max_schur_residual = 0.0;
    double min_schur_reciprocal_condition = 1.0;
};

struct HelmholtzCorrectorResult {
    std::vector<HelmholtzElementCorrector> primal;
    std::vector<HelmholtzElementCorrector> adjoint;
    HelmholtzCorrectorDiagnostics diagnostics;
};

HelmholtzCorrectorResult build_helmholtz_correctors(
    const TriMesh &coarse,
    const TriMesh &fine,
    const Eigen::SparseMatrix<double> &fine_element_prolongation,
    const Eigen::SparseMatrix<double> &fine_dg_prolongation,
    const Eigen::SparseMatrix<double> &quasi_interpolation,
    const Eigen::SparseMatrix<double> &patches,
    const std::vector<TriMesh> &hierarchy_meshes,
    const std::vector<Eigen::SparseMatrix<double>> &node_level_prolongations,
    const std::vector<Eigen::SparseMatrix<double>> &element_level_prolongations,
    const HelmholtzOperators &operators,
    const HelmholtzPatchSolverConfig &solver_config = {});

ComplexSparseMatrix build_helmholtz_corrector_matrix(
    const TriMesh &coarse,
    int fine_node_count,
    const std::vector<HelmholtzElementCorrector> &correctors);

ComplexSparseMatrix build_helmholtz_corrected_basis(
    const Eigen::SparseMatrix<double> &coarse_to_fine,
    const TriMesh &coarse,
    int fine_node_count,
    const std::vector<HelmholtzElementCorrector> &correctors);

} // namespace lod2d::helmholtz

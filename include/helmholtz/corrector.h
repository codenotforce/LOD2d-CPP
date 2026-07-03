#pragma once

#include "helmholtz/operators.h"
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
    const HelmholtzElementBlocks &element_blocks);

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

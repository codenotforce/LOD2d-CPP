#pragma once

#include "helmholtz/operators.h"
#include "helmholtz/patch_solver.h"
#include <Eigen/Sparse>
#include <cstddef>
#include <limits>
#include <memory>
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
    int skipped_patch_count = 0;
    std::size_t skipped_patch_work_units = 0;
    int patches_touching_physical_boundary = 0;
    int parallel_threads = 1;
    int symbolic_analyses = 0;
    int symbolic_reuses = 0;
    int factorization_reuses = 0;
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
    int patch_cache_hits = 0;
    int patch_cache_misses = 0;
    // Misses caused solely by the configured cache-size guard. These are
    // reported separately from mathematically dirty patches so E1 can decide
    // whether a larger bounded cache is useful before spending memory on it.
    int patch_cache_oversized_misses = 0;
    int patch_cache_budget_rejections = 0;
    // Sums of per-patch worker time. They expose where parallel corrector
    // work is spent; correctors_ms in HelmholtzBuildTimings remains wall time.
    double patch_assembly_work_seconds = 0.0;
    double patch_solve_work_seconds = 0.0;
    double patch_pack_work_seconds = 0.0;
    int maximum_patch_dofs = 0;
    int maximum_patch_constraints = 0;
    int maximum_patch_rhs = 0;
};

struct HelmholtzCorrectorResult {
    std::vector<HelmholtzElementCorrector> primal;
    HelmholtzCorrectorDiagnostics diagnostics;
};

// Conservative in-memory cache for complete local patch solves. Hashes only
// locate candidates; a hit additionally requires exact comparison of the
// assembled patch system and every solver option.
class HelmholtzCorrectorPatchCache {
public:
    struct Statistics {
        std::size_t hits = 0;
        std::size_t misses = 0;
        std::size_t stores = 0;
        std::size_t evictions = 0;
        std::size_t oversized_misses = 0;
        std::size_t budget_rejections = 0;
        std::size_t entries = 0;
        std::size_t current_bytes = 0;
        std::size_t peak_bytes = 0;
    };

    explicit HelmholtzCorrectorPatchCache(
        std::size_t maximum_entries = 4096,
        std::size_t maximum_patch_dofs = 4096,
        std::size_t maximum_bytes = std::numeric_limits<std::size_t>::max());
    ~HelmholtzCorrectorPatchCache();
    HelmholtzCorrectorPatchCache(HelmholtzCorrectorPatchCache &&) noexcept;
    HelmholtzCorrectorPatchCache &operator=(HelmholtzCorrectorPatchCache &&) noexcept;
    HelmholtzCorrectorPatchCache(const HelmholtzCorrectorPatchCache &) = delete;
    HelmholtzCorrectorPatchCache &operator=(const HelmholtzCorrectorPatchCache &) = delete;

    void clear();
    Statistics statistics() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;

    friend HelmholtzCorrectorResult build_helmholtz_correctors(
        const TriMesh &,
        const TriMesh &,
        const Eigen::SparseMatrix<double> &,
        const Eigen::SparseMatrix<double> &,
        const Eigen::SparseMatrix<double> &,
        const Eigen::SparseMatrix<double> &,
        const std::vector<TriMesh> &,
        const std::vector<Eigen::SparseMatrix<double>> &,
        const std::vector<Eigen::SparseMatrix<double>> &,
        const HelmholtzOperators &,
        const HelmholtzPatchSolverConfig &,
        HelmholtzCorrectorPatchCache *,
        const std::vector<int> &);
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
    const HelmholtzPatchSolverConfig &solver_config = {},
    HelmholtzCorrectorPatchCache *cache = nullptr,
    const std::vector<int> &skipped_coarse_elements = {});

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

#include "helmholtz/corrector.h"

#include "helmholtz/patch_system.h"

#include <algorithm>
#include <atomic>
#include <exception>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lod2d::helmholtz {
namespace {

struct PatchResult {
    HelmholtzElementCorrector primal;
    HelmholtzPatchSolveDiagnostics diagnostics;
    bool touches_physical_boundary = false;
};

PatchResult solve_and_pack_patch(
    const HelmholtzPatchAssembler &assembler,
    int target,
    const HelmholtzPatchSolverConfig &solver_config) {
    const HelmholtzPatchSystem system = assembler.assemble(target);
    HelmholtzPatchSolveResult solved =
        solve_helmholtz_patch(system, solver_config);

    PatchResult result;
    result.diagnostics = solved.diagnostics;
    result.touches_physical_boundary = system.touches_physical_boundary;
    const int local_size = static_cast<int>(system.local_vertices.size());
    result.primal.reserve(static_cast<std::size_t>(local_size) * solved.corrector.cols());
    for (int row = 0; row < local_size; ++row) {
        for (int column = 0; column < solved.corrector.cols(); ++column) {
            const Complex value = solved.corrector(row, column);
            if (std::abs(value) <= 1e-14) continue;
            result.primal.push_back({system.local_vertices[row], column, value});
        }
    }
    return result;
}

void accumulate_diagnostics(
    const PatchResult &patch,
    HelmholtzCorrectorDiagnostics &diagnostics) {
    const auto &local = patch.diagnostics;
    diagnostics.max_primal_residual =
        std::max(diagnostics.max_primal_residual, local.primal_residual);
    diagnostics.max_adjoint_residual =
        std::max(diagnostics.max_adjoint_residual, local.adjoint_residual);
    diagnostics.max_constraint_residual =
        std::max(diagnostics.max_constraint_residual, local.constraint_residual);
    diagnostics.max_gmres_relative_residual =
        std::max(diagnostics.max_gmres_relative_residual, local.max_gmres_residual);
    diagnostics.max_schur_residual =
        std::max(diagnostics.max_schur_residual, local.schur_residual);
    diagnostics.min_schur_reciprocal_condition =
        std::min(diagnostics.min_schur_reciprocal_condition, local.schur_rcond);
    diagnostics.max_vcycle_relative_residual = std::max(
        diagnostics.max_vcycle_relative_residual,
        local.vcycle_relative_residual);
    diagnostics.max_vcycle_levels =
        std::max(diagnostics.max_vcycle_levels, local.vcycle_levels);
    diagnostics.max_vcycle_coarse_dofs = std::max(
        diagnostics.max_vcycle_coarse_dofs, local.vcycle_coarse_dofs);
    diagnostics.max_vcycle_finest_dofs = std::max(
        diagnostics.max_vcycle_finest_dofs, local.vcycle_finest_dofs);
    diagnostics.gmres_right_hand_sides += local.gmres_right_hand_sides;
    diagnostics.gmres_iterations += local.gmres_total_iterations;
    diagnostics.gmres_max_iterations =
        std::max(diagnostics.gmres_max_iterations, local.gmres_max_iterations);
    diagnostics.gmres_restarts += local.gmres_restarts;
    if (local.direct_fallback) ++diagnostics.direct_fallbacks;
    if (local.symbolic_reused)
        ++diagnostics.symbolic_reuses;
    else
        ++diagnostics.symbolic_analyses;
    if (local.factorization_reused)
        ++diagnostics.factorization_reuses;
    if (patch.touches_physical_boundary)
        ++diagnostics.patches_touching_physical_boundary;
}

} // namespace

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
    const HelmholtzPatchSolverConfig &solver_config) {
    HelmholtzPatchAssembler assembler(
        coarse, fine, fine_element_prolongation, fine_dg_prolongation,
        quasi_interpolation, patches, hierarchy_meshes,
        node_level_prolongations, element_level_prolongations, operators);
    const int patch_count = assembler.patch_count();
    std::vector<int> target_order(patch_count);
    std::iota(target_order.begin(), target_order.end(), 0);
    std::stable_sort(target_order.begin(), target_order.end(), [&](int lhs, int rhs) {
        const std::size_t lhs_cost = assembler.patch_cost(lhs);
        const std::size_t rhs_cost = assembler.patch_cost(rhs);
        return lhs_cost == rhs_cost ? lhs < rhs : lhs_cost > rhs_cost;
    });

    std::vector<PatchResult> patch_results(patch_count);
    std::atomic<bool> failed{false};
    std::exception_ptr first_exception;
    std::mutex exception_mutex;
    int used_threads = 1;
    const auto run_target = [&](int ordered_index) {
        if (failed.load(std::memory_order_relaxed)) return;
        const int target = target_order[ordered_index];
        try {
            patch_results[target] =
                solve_and_pack_patch(assembler, target, solver_config);
        } catch (...) {
            failed.store(true, std::memory_order_relaxed);
            std::lock_guard<std::mutex> lock(exception_mutex);
            if (!first_exception) first_exception = std::current_exception();
        }
    };

#ifdef _OPENMP
    #pragma omp parallel
    {
        #pragma omp single
        used_threads = omp_get_num_threads();
        #pragma omp for schedule(dynamic, 1)
        for (int ordered_index = 0; ordered_index < patch_count; ++ordered_index)
            run_target(ordered_index);
    }
#else
    for (int ordered_index = 0; ordered_index < patch_count; ++ordered_index)
        run_target(ordered_index);
#endif
    if (first_exception) std::rethrow_exception(first_exception);

    HelmholtzCorrectorResult result;
    result.primal.resize(patch_count);
    result.diagnostics.patch_count = patch_count;
    result.diagnostics.parallel_threads = used_threads;
    for (int target = 0; target < patch_count; ++target) {
        result.primal[target] = std::move(patch_results[target].primal);
        accumulate_diagnostics(patch_results[target], result.diagnostics);
    }
    return result;
}

} // namespace lod2d::helmholtz

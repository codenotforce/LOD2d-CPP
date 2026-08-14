#include "helmholtz/corrector.h"

#include "helmholtz/patch_system.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <exception>
#include <memory>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
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
    bool cache_hit = false;
    double assembly_seconds = 0.0;
    double solve_seconds = 0.0;
    double pack_seconds = 0.0;
    int patch_dofs = 0;
    int patch_constraints = 0;
    int patch_rhs = 0;
};

PatchResult pack_patch_solution(
    const HelmholtzPatchSystem &system,
    const HelmholtzPatchSolveResult &solved,
    const bool cache_hit) {
    PatchResult result;
    result.diagnostics = solved.diagnostics;
    result.touches_physical_boundary = system.touches_physical_boundary;
    result.cache_hit = cache_hit;
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
    if (patch.cache_hit) {
        ++diagnostics.patch_cache_hits;
    } else {
        ++diagnostics.patch_cache_misses;
        diagnostics.gmres_right_hand_sides += local.gmres_right_hand_sides;
        diagnostics.gmres_iterations += local.gmres_total_iterations;
        diagnostics.gmres_max_iterations = std::max(
            diagnostics.gmres_max_iterations, local.gmres_max_iterations);
        diagnostics.gmres_restarts += local.gmres_restarts;
        if (local.direct_fallback) ++diagnostics.direct_fallbacks;
        if (local.symbolic_reused)
            ++diagnostics.symbolic_reuses;
        else
            ++diagnostics.symbolic_analyses;
        if (local.factorization_reused)
            ++diagnostics.factorization_reuses;
    }
    if (patch.touches_physical_boundary)
        ++diagnostics.patches_touching_physical_boundary;
    diagnostics.patch_assembly_work_seconds += patch.assembly_seconds;
    diagnostics.patch_solve_work_seconds += patch.solve_seconds;
    diagnostics.patch_pack_work_seconds += patch.pack_seconds;
    diagnostics.maximum_patch_dofs = std::max(
        diagnostics.maximum_patch_dofs, patch.patch_dofs);
    diagnostics.maximum_patch_constraints = std::max(
        diagnostics.maximum_patch_constraints, patch.patch_constraints);
    diagnostics.maximum_patch_rhs = std::max(
        diagnostics.maximum_patch_rhs, patch.patch_rhs);
}

constexpr std::uint64_t fnv_offset = 14695981039346656037ull;
constexpr std::uint64_t fnv_prime = 1099511628211ull;

void hash_bytes(std::uint64_t &hash, const void *data, const std::size_t size) {
    const auto *bytes = static_cast<const unsigned char *>(data);
    for (std::size_t i = 0; i < size; ++i) {
        hash ^= bytes[i];
        hash *= fnv_prime;
    }
}

template <class Scalar>
void hash_scalar(std::uint64_t &hash, const Scalar &value) {
    static_assert(std::is_trivially_copyable_v<Scalar>);
    hash_bytes(hash, std::addressof(value), sizeof(value));
}

template <class Scalar>
void hash_vector(std::uint64_t &hash, const std::vector<Scalar> &values) {
    hash_scalar(hash, values.size());
    for (const Scalar &value : values) hash_scalar(hash, value);
}

template <class Scalar>
void hash_sparse(
    std::uint64_t &hash,
    const Eigen::SparseMatrix<Scalar> &matrix) {
    hash_scalar(hash, matrix.rows());
    hash_scalar(hash, matrix.cols());
    hash_scalar(hash, matrix.nonZeros());
    for (int outer = 0; outer < matrix.outerSize(); ++outer) {
        for (typename Eigen::SparseMatrix<Scalar>::InnerIterator it(matrix, outer);
             it; ++it) {
            hash_scalar(hash, it.row());
            hash_scalar(hash, it.col());
            hash_scalar(hash, it.value());
        }
    }
}

template <class Derived>
void hash_dense(std::uint64_t &hash, const Eigen::MatrixBase<Derived> &matrix) {
    hash_scalar(hash, matrix.rows());
    hash_scalar(hash, matrix.cols());
    for (Eigen::Index column = 0; column < matrix.cols(); ++column)
        for (Eigen::Index row = 0; row < matrix.rows(); ++row)
            hash_scalar(hash, matrix(row, column));
}

void hash_solver_config(
    std::uint64_t &hash,
    const HelmholtzPatchSolverConfig &config) {
    hash_scalar(hash, config.kind);
    hash_scalar(hash, config.symbolic_cache_slots);
    hash_scalar(hash, config.reuse_identical_factorization);
    hash_scalar(hash, config.gmres.restart);
    hash_scalar(hash, config.gmres.max_iterations);
    hash_scalar(hash, config.gmres.relative_tolerance);
    hash_scalar(hash, config.gmres.absolute_tolerance);
    hash_scalar(hash, config.gmres.reorthogonalize);
    hash_scalar(hash, config.shifted.rule);
    hash_scalar(hash, config.shifted.alpha);
    hash_scalar(hash, config.shifted.absolute_epsilon);
    hash_scalar(hash, config.shifted.inverse);
    hash_scalar(hash, config.shifted.pre_smooth);
    hash_scalar(hash, config.shifted.post_smooth);
    hash_scalar(hash, config.shifted.coarse_max_dofs);
    hash_scalar(hash, config.shifted.jacobi_weight);
    hash_scalar(hash, config.fallback_to_direct);
}

std::uint64_t patch_cache_hash(
    const HelmholtzPatchSystem &system,
    const HelmholtzPatchSolverConfig &config) {
    std::uint64_t hash = fnv_offset;
    // target_element is intentionally omitted: rhs and local/global mappings
    // already encode the mathematical source element.
    hash_vector(hash, system.local_vertices);
    hash_vector(hash, system.patch_elements);
    hash_sparse(hash, system.stiffness);
    hash_sparse(hash, system.mass);
    hash_sparse(hash, system.robin);
    hash_sparse(hash, system.helmholtz);
    hash_dense(hash, system.constraints);
    hash_dense(hash, system.rhs);
    hash_scalar(hash, system.wavenumber);
    hash_scalar(hash, system.diameter);
    hash_scalar(hash, system.touches_physical_boundary);
    hash_scalar(hash, system.geometric_prolongations.size());
    for (const auto &prolongation : system.geometric_prolongations)
        hash_sparse(hash, prolongation);
    hash_solver_config(hash, config);
    return hash;
}

template <class Scalar>
bool sparse_equal(
    const Eigen::SparseMatrix<Scalar> &lhs,
    const Eigen::SparseMatrix<Scalar> &rhs) {
    if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()
        || lhs.nonZeros() != rhs.nonZeros())
        return false;
    for (int outer = 0; outer < lhs.outerSize(); ++outer) {
        typename Eigen::SparseMatrix<Scalar>::InnerIterator a(lhs, outer);
        typename Eigen::SparseMatrix<Scalar>::InnerIterator b(rhs, outer);
        for (; a && b; ++a, ++b) {
            if (a.row() != b.row() || a.col() != b.col()
                || a.value() != b.value())
                return false;
        }
        if (a || b) return false;
    }
    return true;
}

template <class Left, class Right>
bool dense_equal(
    const Eigen::MatrixBase<Left> &lhs,
    const Eigen::MatrixBase<Right> &rhs) {
    return lhs.rows() == rhs.rows() && lhs.cols() == rhs.cols()
        && (lhs.array() == rhs.array()).all();
}

bool solver_config_equal(
    const HelmholtzPatchSolverConfig &lhs,
    const HelmholtzPatchSolverConfig &rhs) {
    return lhs.kind == rhs.kind
        && lhs.symbolic_cache_slots == rhs.symbolic_cache_slots
        && lhs.reuse_identical_factorization == rhs.reuse_identical_factorization
        && lhs.gmres.restart == rhs.gmres.restart
        && lhs.gmres.max_iterations == rhs.gmres.max_iterations
        && lhs.gmres.relative_tolerance == rhs.gmres.relative_tolerance
        && lhs.gmres.absolute_tolerance == rhs.gmres.absolute_tolerance
        && lhs.gmres.reorthogonalize == rhs.gmres.reorthogonalize
        && lhs.shifted.rule == rhs.shifted.rule
        && lhs.shifted.alpha == rhs.shifted.alpha
        && lhs.shifted.absolute_epsilon == rhs.shifted.absolute_epsilon
        && lhs.shifted.inverse == rhs.shifted.inverse
        && lhs.shifted.pre_smooth == rhs.shifted.pre_smooth
        && lhs.shifted.post_smooth == rhs.shifted.post_smooth
        && lhs.shifted.coarse_max_dofs == rhs.shifted.coarse_max_dofs
        && lhs.shifted.jacobi_weight == rhs.shifted.jacobi_weight
        && lhs.fallback_to_direct == rhs.fallback_to_direct;
}

bool patch_system_equal(
    const HelmholtzPatchSystem &lhs,
    const HelmholtzPatchSystem &rhs) {
    if (lhs.local_vertices != rhs.local_vertices
        || lhs.patch_elements != rhs.patch_elements
        || lhs.wavenumber != rhs.wavenumber
        || lhs.diameter != rhs.diameter
        || lhs.touches_physical_boundary != rhs.touches_physical_boundary
        || lhs.geometric_prolongations.size()
            != rhs.geometric_prolongations.size()
        || !sparse_equal(lhs.stiffness, rhs.stiffness)
        || !sparse_equal(lhs.mass, rhs.mass)
        || !sparse_equal(lhs.robin, rhs.robin)
        || !sparse_equal(lhs.helmholtz, rhs.helmholtz)
        || !dense_equal(lhs.constraints, rhs.constraints)
        || !dense_equal(lhs.rhs, rhs.rhs))
        return false;
    for (std::size_t i = 0; i < lhs.geometric_prolongations.size(); ++i)
        if (!sparse_equal(
                lhs.geometric_prolongations[i],
                rhs.geometric_prolongations[i]))
            return false;
    return true;
}

} // namespace

struct HelmholtzCorrectorPatchCache::Impl {
    struct Entry {
        HelmholtzPatchSystem system;
        HelmholtzPatchSolverConfig solver_config;
        HelmholtzPatchSolveResult solved;
    };

    explicit Impl(const std::size_t maximum) : maximum_entries(maximum) {
        if (maximum_entries == 0)
            throw std::invalid_argument(
                "Helmholtz corrector patch cache capacity must be positive");
    }

    bool lookup(
        const HelmholtzPatchSystem &system,
        const HelmholtzPatchSolverConfig &config,
        HelmholtzPatchSolveResult &solved) {
        const std::uint64_t hash = patch_cache_hash(system, config);
        std::lock_guard<std::mutex> lock(mutex);
        const auto range = entries.equal_range(hash);
        for (auto it = range.first; it != range.second; ++it) {
            if (solver_config_equal(it->second.solver_config, config)
                && patch_system_equal(it->second.system, system)) {
                solved = it->second.solved;
                ++stats.hits;
                return true;
            }
        }
        ++stats.misses;
        return false;
    }

    void store(
        HelmholtzPatchSystem system,
        const HelmholtzPatchSolverConfig &config,
        HelmholtzPatchSolveResult solved) {
        const std::uint64_t hash = patch_cache_hash(system, config);
        std::lock_guard<std::mutex> lock(mutex);
        if (entries.size() >= maximum_entries) {
            stats.evictions += entries.size();
            entries.clear();
        }
        entries.emplace(hash, Entry{std::move(system), config, std::move(solved)});
        ++stats.stores;
    }

    const std::size_t maximum_entries;
    mutable std::mutex mutex;
    std::unordered_multimap<std::uint64_t, Entry> entries;
    Statistics stats;
};

HelmholtzCorrectorPatchCache::HelmholtzCorrectorPatchCache(
    const std::size_t maximum_entries)
    : impl_(std::make_unique<Impl>(maximum_entries)) {}

HelmholtzCorrectorPatchCache::~HelmholtzCorrectorPatchCache() = default;
HelmholtzCorrectorPatchCache::HelmholtzCorrectorPatchCache(
    HelmholtzCorrectorPatchCache &&) noexcept = default;
HelmholtzCorrectorPatchCache &HelmholtzCorrectorPatchCache::operator=(
    HelmholtzCorrectorPatchCache &&) noexcept = default;

void HelmholtzCorrectorPatchCache::clear() {
    std::lock_guard<std::mutex> lock(impl_->mutex);
    impl_->entries.clear();
    impl_->stats.entries = 0;
}

HelmholtzCorrectorPatchCache::Statistics
HelmholtzCorrectorPatchCache::statistics() const {
    std::lock_guard<std::mutex> lock(impl_->mutex);
    Statistics result = impl_->stats;
    result.entries = impl_->entries.size();
    return result;
}

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
    const HelmholtzPatchSolverConfig &solver_config,
    HelmholtzCorrectorPatchCache *cache) {
    HelmholtzPatchAssembler assembler(
        coarse, fine, fine_element_prolongation, fine_dg_prolongation,
        quasi_interpolation, patches, hierarchy_meshes,
        node_level_prolongations, element_level_prolongations, operators);
    const int patch_count = assembler.patch_count();
    if (solver_config.maximum_parallel_solves < 0) {
        throw std::invalid_argument(
            "maximum parallel patch solves must be nonnegative");
    }
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
            const auto assembly_begin = std::chrono::steady_clock::now();
            HelmholtzPatchSystem system = assembler.assemble(target);
            const double assembly_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - assembly_begin).count();
            const auto solve_begin = std::chrono::steady_clock::now();
            HelmholtzPatchSolveResult solved;
            const bool cache_hit = cache != nullptr
                && cache->impl_->lookup(system, solver_config, solved);
            if (!cache_hit) {
                solved = solve_helmholtz_patch(system, solver_config);
                if (cache != nullptr)
                    cache->impl_->store(system, solver_config, solved);
            }
            const double solve_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - solve_begin).count();
            const auto pack_begin = std::chrono::steady_clock::now();
            PatchResult packed = pack_patch_solution(system, solved, cache_hit);
            packed.pack_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - pack_begin).count();
            packed.assembly_seconds = assembly_seconds;
            packed.solve_seconds = solve_seconds;
            packed.patch_dofs = system.helmholtz.rows();
            packed.patch_constraints = system.constraints.rows();
            packed.patch_rhs = system.rhs.cols();
            patch_results[target] = std::move(packed);
        } catch (...) {
            failed.store(true, std::memory_order_relaxed);
            std::lock_guard<std::mutex> lock(exception_mutex);
            if (!first_exception) first_exception = std::current_exception();
        }
    };

#ifdef _OPENMP
    const int requested_threads = std::max(
        1, solver_config.maximum_parallel_solves > 0
            ? std::min(solver_config.maximum_parallel_solves, patch_count)
            : std::min(omp_get_max_threads(), patch_count));
    #pragma omp parallel num_threads(requested_threads)
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

#include "helmholtz/corrector.h"

#include "helmholtz/patch_system.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <exception>
#include <list>
#include <map>
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
    bool cache_oversized = false;
    bool cache_budget_rejected = false;
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
        if (patch.cache_oversized)
            ++diagnostics.patch_cache_oversized_misses;
        if (patch.cache_budget_rejected)
            ++diagnostics.patch_cache_budget_rejections;
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
        std::uint64_t hash = 0;
        HelmholtzPatchSystem system;
        HelmholtzPatchSolverConfig solver_config;
        HelmholtzPatchSolveResult solved;
        std::size_t bytes = 0;
    };

    using EntryList = std::list<Entry>;
    using EntryIterator = EntryList::iterator;

    Impl(
        const std::size_t maximum,
        const std::size_t maximum_dofs,
        const std::size_t maximum_memory)
        : maximum_entries(maximum), maximum_patch_dofs(maximum_dofs),
          maximum_bytes(maximum_memory) {
        if (maximum_entries == 0 || maximum_patch_dofs == 0
            || maximum_bytes == 0)
            throw std::invalid_argument(
                "Helmholtz corrector patch cache bounds must be positive");
    }

    template <class Scalar, int Options, class Index>
    static std::size_t sparse_bytes(
        const Eigen::SparseMatrix<Scalar, Options, Index> &matrix) {
        return sizeof(Scalar) * static_cast<std::size_t>(matrix.nonZeros())
            + sizeof(Index) * static_cast<std::size_t>(matrix.nonZeros())
            + sizeof(Index) * static_cast<std::size_t>(matrix.outerSize() + 1);
    }

    template <class Matrix>
    static std::size_t dense_bytes(const Matrix &matrix) {
        return sizeof(typename Matrix::Scalar)
            * static_cast<std::size_t>(matrix.size());
    }

    static std::size_t entry_bytes(
        const HelmholtzPatchSystem &system,
        const HelmholtzPatchSolveResult &solved) {
        std::size_t bytes = sizeof(Entry)
            + sizeof(int) * (system.local_vertices.capacity()
                             + system.patch_elements.capacity())
            + sparse_bytes(system.stiffness)
            + sparse_bytes(system.mass)
            + sparse_bytes(system.robin)
            + sparse_bytes(system.helmholtz)
            + dense_bytes(system.constraints)
            + dense_bytes(system.rhs)
            + dense_bytes(solved.corrector)
            + dense_bytes(solved.multipliers);
        for (const auto &prolongation : system.geometric_prolongations)
            bytes += sparse_bytes(prolongation);
        return bytes;
    }

    bool lookup(
        const HelmholtzPatchSystem &system,
        const HelmholtzPatchSolverConfig &config,
        HelmholtzPatchSolveResult &solved) {
        if (static_cast<std::size_t>(system.helmholtz.rows())
            > maximum_patch_dofs) {
            std::lock_guard<std::mutex> lock(mutex);
            ++stats.misses;
            ++stats.oversized_misses;
            return false;
        }
        const std::uint64_t hash = patch_cache_hash(system, config);
        std::lock_guard<std::mutex> lock(mutex);
        const auto range = index.equal_range(hash);
        for (auto it = range.first; it != range.second; ++it) {
            const EntryIterator entry = it->second;
            if (solver_config_equal(entry->solver_config, config)
                && patch_system_equal(entry->system, system)) {
                solved = entry->solved;
                entries.splice(entries.end(), entries, entry);
                ++stats.hits;
                return true;
            }
        }
        ++stats.misses;
        return false;
    }

    bool store(
        HelmholtzPatchSystem system,
        const HelmholtzPatchSolverConfig &config,
        HelmholtzPatchSolveResult solved) {
        if (static_cast<std::size_t>(system.helmholtz.rows())
            > maximum_patch_dofs)
            return false;
        const std::size_t bytes = entry_bytes(system, solved);
        const std::uint64_t hash = patch_cache_hash(system, config);
        std::lock_guard<std::mutex> lock(mutex);
        if (bytes > maximum_bytes) {
            ++stats.budget_rejections;
            return false;
        }
        while (!entries.empty()
               && (entries.size() >= maximum_entries
                   || bytes > maximum_bytes - stats.current_bytes)) {
            const EntryIterator oldest = entries.begin();
            const auto range = index.equal_range(oldest->hash);
            for (auto it = range.first; it != range.second; ++it) {
                if (it->second == oldest) {
                    index.erase(it);
                    break;
                }
            }
            stats.current_bytes -= oldest->bytes;
            entries.erase(oldest);
            ++stats.evictions;
        }
        if (bytes > maximum_bytes - stats.current_bytes) {
            ++stats.budget_rejections;
            return false;
        }
        entries.push_back(Entry{
            hash, std::move(system), config, std::move(solved), bytes});
        index.emplace(hash, std::prev(entries.end()));
        ++stats.stores;
        stats.current_bytes += bytes;
        stats.peak_bytes = std::max(stats.peak_bytes, stats.current_bytes);
        return true;
    }

    const std::size_t maximum_entries;
    const std::size_t maximum_patch_dofs;
    const std::size_t maximum_bytes;
    mutable std::mutex mutex;
    EntryList entries;
    std::unordered_multimap<std::uint64_t, EntryIterator> index;
    Statistics stats;
};

HelmholtzCorrectorPatchCache::HelmholtzCorrectorPatchCache(
    const std::size_t maximum_entries,
    const std::size_t maximum_patch_dofs,
    const std::size_t maximum_bytes)
    : impl_(std::make_unique<Impl>(
          maximum_entries, maximum_patch_dofs, maximum_bytes)) {}

HelmholtzCorrectorPatchCache::~HelmholtzCorrectorPatchCache() = default;
HelmholtzCorrectorPatchCache::HelmholtzCorrectorPatchCache(
    HelmholtzCorrectorPatchCache &&) noexcept = default;
HelmholtzCorrectorPatchCache &HelmholtzCorrectorPatchCache::operator=(
    HelmholtzCorrectorPatchCache &&) noexcept = default;

void HelmholtzCorrectorPatchCache::clear() {
    std::lock_guard<std::mutex> lock(impl_->mutex);
    impl_->entries.clear();
    impl_->index.clear();
    impl_->stats.entries = 0;
    impl_->stats.current_bytes = 0;
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
    HelmholtzCorrectorPatchCache *cache,
    const std::vector<int> &skipped_coarse_elements) {
    HelmholtzPatchAssembler assembler(
        coarse, fine, fine_element_prolongation, fine_dg_prolongation,
        quasi_interpolation, patches, hierarchy_meshes,
        node_level_prolongations, element_level_prolongations, operators);
    const int patch_count = assembler.patch_count();
    if (solver_config.maximum_parallel_solves < 0) {
        throw std::invalid_argument(
            "maximum parallel patch solves must be nonnegative");
    }
    std::vector<char> skipped(patch_count, false);
    for (int target : skipped_coarse_elements) {
        if (target < 0 || target >= patch_count)
            throw std::out_of_range("skipped corrector element is out of range");
        if (skipped[target])
            throw std::invalid_argument("skipped corrector elements contain a duplicate");
        skipped[target] = true;
    }
    std::vector<int> target_order;
    target_order.reserve(patch_count - skipped_coarse_elements.size());
    for (int target = 0; target < patch_count; ++target)
        if (!skipped[target]) target_order.push_back(target);
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
    const auto run_target = [&](
        const int target,
        const bool shared_direct_schur_block) {
        if (failed.load(std::memory_order_relaxed)) return;
        try {
            const auto assembly_begin = std::chrono::steady_clock::now();
            HelmholtzPatchSystem system = assembler.assemble(target);
            const double assembly_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - assembly_begin).count();
            const auto solve_begin = std::chrono::steady_clock::now();
            HelmholtzPatchSolveResult solved;
            HelmholtzPatchSolverConfig effective_solver = solver_config;
            if (solver_config.kind == HelmholtzPatchSolverKind::DirectSchur
                && solver_config.reuse_identical_factorization
                && system.helmholtz.rows() <= 8192
                && (!shared_direct_schur_block
                    || system.constraints.rows() > 32)) {
                // Schur reduction pays off when one A factorization serves
                // several small constraint/RHS blocks.  A singleton support
                // or a large dense Schur complement on a modest patch is
                // faster with the mathematically equivalent saddle solve.
                // Large patches never take this fallback: sparse indefinite
                // saddle fill-in dominates the Schur overhead there.
                effective_solver.kind = HelmholtzPatchSolverKind::DirectSaddle;
                effective_solver.reuse_identical_factorization = false;
            }
            const bool cache_hit = cache != nullptr
                && cache->impl_->lookup(system, effective_solver, solved);
            const bool cache_oversized = cache != nullptr
                && static_cast<std::size_t>(system.helmholtz.rows())
                    > cache->impl_->maximum_patch_dofs;
            if (!cache_hit) {
                solved = solve_helmholtz_patch(system, effective_solver);
            }
            const double solve_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - solve_begin).count();
            const auto pack_begin = std::chrono::steady_clock::now();
            PatchResult packed = pack_patch_solution(system, solved, cache_hit);
            packed.cache_oversized = cache_oversized;
            packed.assembly_seconds = assembly_seconds;
            packed.solve_seconds = solve_seconds;
            packed.patch_dofs = system.helmholtz.rows();
            packed.patch_constraints = system.constraints.rows();
            packed.patch_rhs = system.rhs.cols();
            if (!cache_hit && cache != nullptr && !cache_oversized) {
                packed.cache_budget_rejected = !cache->impl_->store(
                    std::move(system), effective_solver, std::move(solved));
            }
            packed.pack_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - pack_begin).count();
            patch_results[target] = std::move(packed);
        } catch (...) {
            failed.store(true, std::memory_order_relaxed);
            std::lock_guard<std::mutex> lock(exception_mutex);
            if (!first_exception) first_exception = std::current_exception();
        }
    };

    std::vector<std::vector<int>> work_groups;
    if (solver_config.kind == HelmholtzPatchSolverKind::DirectSchur
        && solver_config.reuse_identical_factorization) {
        // DirectSchur factorizes only the local energy block.  All targets
        // with exactly the same fine-element support have the same block but
        // different constraints/RHS.  Keep such targets on one worker so the
        // thread-local SparseLU cache factorizes the block once.
        std::map<std::vector<int>, std::vector<int>> grouped;
        for (const int target : target_order)
            grouped[assembler.patch_fine_elements(target)].push_back(target);
        work_groups.reserve(grouped.size());
        for (auto &[support, targets] : grouped)
            work_groups.push_back(std::move(targets));
        std::stable_sort(
            work_groups.begin(), work_groups.end(), [&](const auto &lhs, const auto &rhs) {
                const auto work = [&](const std::vector<int> &targets) {
                    return std::accumulate(
                        targets.begin(), targets.end(), std::size_t{0},
                        [&](const std::size_t sum, const int target) {
                            return sum + assembler.patch_cost(target);
                        });
                };
                return work(lhs) > work(rhs);
            });
    } else {
        work_groups.reserve(target_order.size());
        for (const int target : target_order) work_groups.push_back({target});
    }
    const int group_count = static_cast<int>(work_groups.size());
#ifdef _OPENMP
    const int requested_threads = std::max(
        1, solver_config.maximum_parallel_solves > 0
            ? std::min(solver_config.maximum_parallel_solves,
                       std::max(1, group_count))
            : std::min(omp_get_max_threads(),
                       std::max(1, group_count)));
    #pragma omp parallel num_threads(requested_threads)
    {
        #pragma omp single
        used_threads = omp_get_num_threads();
        #pragma omp for schedule(dynamic, 1)
        for (int group = 0; group < group_count; ++group)
            for (const int target : work_groups[group])
                run_target(target, work_groups[group].size() > 1);
    }
#else
    for (const auto &group : work_groups)
        for (const int target : group)
            run_target(target, group.size() > 1);
#endif
    if (first_exception) std::rethrow_exception(first_exception);

    HelmholtzCorrectorResult result;
    result.primal.resize(patch_count);
    result.diagnostics.patch_count = static_cast<int>(target_order.size());
    result.diagnostics.skipped_patch_count =
        static_cast<int>(skipped_coarse_elements.size());
    for (int target : skipped_coarse_elements)
        result.diagnostics.skipped_patch_work_units += assembler.patch_cost(target);
    result.diagnostics.parallel_threads = used_threads;
    for (int target : target_order) {
        result.primal[target] = std::move(patch_results[target].primal);
        accumulate_diagnostics(patch_results[target], result.diagnostics);
    }
    return result;
}

} // namespace lod2d::helmholtz

#include "helmholtz/corrector.h"

#include <Eigen/QR>
#include <Eigen/SparseLU>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <exception>
#include <limits>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <unordered_map>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lod2d::helmholtz {
namespace {

using ComplexSparseLu = Eigen::SparseLU<ComplexSparseMatrix>;

struct PatchSolveResult {
    HelmholtzElementCorrector primal;
    HelmholtzElementCorrector adjoint;
    double primal_residual = 0.0;
    double adjoint_residual = 0.0;
    double constraint_residual = 0.0;
    bool touches_physical_boundary = false;
    bool symbolic_reused = false;
};

struct PatchWorkspace {
    std::vector<int> node_count;
    std::vector<int> node_seen;
    std::vector<int> local_index;
    std::vector<int> local_seen;
    std::vector<int> coarse_row_index;
    std::vector<int> coarse_row_seen;
    int stamp = 0;

    void begin(size_t fine_nodes, size_t coarse_nodes) {
        if (node_count.size() != fine_nodes) {
            node_count.assign(fine_nodes, 0);
            node_seen.assign(fine_nodes, 0);
            local_index.assign(fine_nodes, -1);
            local_seen.assign(fine_nodes, 0);
            stamp = 0;
        }
        if (coarse_row_index.size() != coarse_nodes) {
            coarse_row_index.assign(coarse_nodes, -1);
            coarse_row_seen.assign(coarse_nodes, 0);
            stamp = 0;
            std::fill(node_seen.begin(), node_seen.end(), 0);
            std::fill(local_seen.begin(), local_seen.end(), 0);
        }
        if (stamp == std::numeric_limits<int>::max()) {
            std::fill(node_seen.begin(), node_seen.end(), 0);
            std::fill(local_seen.begin(), local_seen.end(), 0);
            std::fill(coarse_row_seen.begin(), coarse_row_seen.end(), 0);
            stamp = 0;
        }
        ++stamp;
    }
};

struct SparseLuPatternCache {
    ComplexSparseLu solver;
    std::vector<int> outer_indices;
    std::vector<int> inner_indices;
    int rows = -1;
    int cols = -1;

    bool matches(const ComplexSparseMatrix &matrix) const {
        if (rows != matrix.rows() || cols != matrix.cols()) return false;
        if (outer_indices.size() != static_cast<size_t>(matrix.outerSize() + 1)
            || inner_indices.size() != static_cast<size_t>(matrix.nonZeros()))
            return false;
        return std::equal(
                   outer_indices.begin(), outer_indices.end(), matrix.outerIndexPtr())
            && std::equal(
                   inner_indices.begin(), inner_indices.end(), matrix.innerIndexPtr());
    }

    ComplexSparseLu &factorize(const ComplexSparseMatrix &matrix, bool &reused) {
        reused = matches(matrix);
        if (!reused) {
            solver.analyzePattern(matrix);
            if (solver.info() != Eigen::Success)
                throw std::runtime_error("Helmholtz corrector symbolic analysis failed");
            rows = matrix.rows();
            cols = matrix.cols();
            outer_indices.assign(
                matrix.outerIndexPtr(), matrix.outerIndexPtr() + matrix.outerSize() + 1);
            inner_indices.assign(
                matrix.innerIndexPtr(), matrix.innerIndexPtr() + matrix.nonZeros());
        }
        solver.factorize(matrix);
        if (solver.info() != Eigen::Success)
            throw std::runtime_error("Helmholtz corrector saddle factorization failed");
        return solver;
    }
};

thread_local PatchWorkspace patch_workspace;
thread_local SparseLuPatternCache sparse_lu_cache;

std::uint64_t edge_key(int a, int b) {
    const auto lo = static_cast<std::uint32_t>(std::min(a, b));
    const auto hi = static_cast<std::uint32_t>(std::max(a, b));
    return (static_cast<std::uint64_t>(lo) << 32U) | hi;
}

std::vector<int> node_incidence(const TriMesh &mesh) {
    std::vector<int> incidence(mesh.nodes.size(), 0);
    for (const Triangle &tri : mesh.elems)
        for (int vertex : tri) ++incidence[vertex];
    return incidence;
}

std::vector<std::vector<int>> fine_element_children(
    const Eigen::SparseMatrix<double> &prolongation,
    int coarse_element_count) {
    if (prolongation.cols() != coarse_element_count)
        throw std::invalid_argument("fine element prolongation has the wrong column count");
    std::vector<std::vector<int>> children(coarse_element_count);
    for (int coarse_element = 0; coarse_element < prolongation.outerSize(); ++coarse_element) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(prolongation, coarse_element); it; ++it) {
            if (it.value() != 0.0) children[coarse_element].push_back(it.row());
        }
    }
    return children;
}

std::unordered_map<std::uint64_t, int> edge_counts(const TriMesh &mesh) {
    std::unordered_map<std::uint64_t, int> counts;
    counts.reserve(static_cast<size_t>(3 * mesh.elems.size()));
    for (const Triangle &tri : mesh.elems) {
        ++counts[edge_key(tri[0], tri[1])];
        ++counts[edge_key(tri[1], tri[2])];
        ++counts[edge_key(tri[2], tri[0])];
    }
    return counts;
}

PatchSolveResult solve_patch(
    int target,
    const TriMesh &fine,
    int coarse_node_count,
    const std::vector<std::vector<int>> &children,
    const std::vector<int> &fine_incidence,
    const std::unordered_map<std::uint64_t, int> &fine_edge_counts,
    const Eigen::SparseMatrix<double> &fine_dg_prolongation,
    const Eigen::SparseMatrix<double> &quasi_interpolation,
    const Eigen::SparseMatrix<double> &patches,
    const HelmholtzElementBlocks &element_blocks) {
    std::vector<int> patch_elements;
    for (Eigen::SparseMatrix<double>::InnerIterator it(patches, target); it; ++it) {
        if (it.value() == 0.0) continue;
        const auto &element_children = children[it.row()];
        patch_elements.insert(patch_elements.end(), element_children.begin(), element_children.end());
    }
    if (patch_elements.empty()) throw std::runtime_error("Helmholtz corrector patch is empty");

    patch_workspace.begin(fine.nodes.size(), coarse_node_count);
    const int stamp = patch_workspace.stamp;
    std::vector<int> touched_vertices;
    touched_vertices.reserve(3 * patch_elements.size());
    for (int element : patch_elements) {
        for (int vertex : fine.elems[element]) {
            if (patch_workspace.node_seen[vertex] != stamp) {
                patch_workspace.node_seen[vertex] = stamp;
                patch_workspace.node_count[vertex] = 0;
                touched_vertices.push_back(vertex);
            }
            ++patch_workspace.node_count[vertex];
        }
    }

    std::vector<int> local_vertices;
    local_vertices.reserve(touched_vertices.size());
    for (int vertex : touched_vertices) {
        if (patch_workspace.node_count[vertex] != fine_incidence[vertex]) continue;
        patch_workspace.local_seen[vertex] = stamp;
        patch_workspace.local_index[vertex] = static_cast<int>(local_vertices.size());
        local_vertices.push_back(vertex);
    }
    if (local_vertices.empty())
        throw std::runtime_error("Helmholtz corrector patch has no unconstrained fine vertices");
    auto local_index_of = [&](int vertex) {
        return patch_workspace.local_seen[vertex] == stamp
            ? patch_workspace.local_index[vertex]
            : -1;
    };

    PatchSolveResult result;
    static constexpr int local_edges[3][2] = {{0, 1}, {1, 2}, {2, 0}};
    for (int element : patch_elements) {
        const Triangle &tri = fine.elems[element];
        for (const auto &edge : local_edges) {
            if (fine_edge_counts.at(edge_key(tri[edge[0]], tri[edge[1]])) == 1) {
                result.touches_physical_boundary = true;
                break;
            }
        }
        if (result.touches_physical_boundary) break;
    }

    const int local_size = static_cast<int>(local_vertices.size());
    std::vector<ComplexTriplet> local_triplets;
    local_triplets.reserve(9 * patch_elements.size());
    for (int element : patch_elements) {
        const Triangle &tri = fine.elems[element];
        const Eigen::Matrix3cd &block = element_blocks[element];
        for (int i = 0; i < 3; ++i) {
            const int row = local_index_of(tri[i]);
            if (row < 0) continue;
            for (int j = 0; j < 3; ++j) {
                const int col = local_index_of(tri[j]);
                if (col >= 0 && std::abs(block(i, j)) > 0.0)
                    local_triplets.emplace_back(row, col, block(i, j));
            }
        }
    }
    ComplexSparseMatrix local_operator(local_size, local_size);
    local_operator.setFromTriplets(local_triplets.begin(), local_triplets.end());

    ComplexMatrix rhs = ComplexMatrix::Zero(local_size, 3);
    for (int element : children[target]) {
        const Triangle &tri = fine.elems[element];
        const Eigen::Matrix3cd &block = element_blocks[element];
        for (int i = 0; i < 3; ++i) {
            const int row = local_index_of(tri[i]);
            if (row < 0) continue;
            for (int coarse_local = 0; coarse_local < 3; ++coarse_local) {
                Complex value = 0.0;
                for (int trial_local = 0; trial_local < 3; ++trial_local) {
                    value += block(i, trial_local)
                           * fine_dg_prolongation.coeff(
                               3 * element + trial_local,
                               3 * target + coarse_local);
                }
                rhs(row, coarse_local) += value;
            }
        }
    }

    std::vector<int> active_constraint_rows;
    for (int local_col = 0; local_col < local_size; ++local_col) {
        const int global_col = local_vertices[local_col];
        for (Eigen::SparseMatrix<double>::InnerIterator it(quasi_interpolation, global_col); it; ++it) {
            const int coarse_row = it.row();
            if (patch_workspace.coarse_row_seen[coarse_row] == stamp) continue;
            patch_workspace.coarse_row_seen[coarse_row] = stamp;
            patch_workspace.coarse_row_index[coarse_row]
                = static_cast<int>(active_constraint_rows.size());
            active_constraint_rows.push_back(coarse_row);
        }
    }

    Eigen::MatrixXd active_constraints = Eigen::MatrixXd::Zero(
        active_constraint_rows.size(), local_size);
    for (int local_col = 0; local_col < local_size; ++local_col) {
        const int global_col = local_vertices[local_col];
        for (Eigen::SparseMatrix<double>::InnerIterator it(quasi_interpolation, global_col); it; ++it) {
            const int row = patch_workspace.coarse_row_index[it.row()];
            active_constraints(row, local_col) = it.value();
        }
    }

    Eigen::MatrixXd constraints;
    if (active_constraints.rows() > 0) {
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(active_constraints.transpose());
        qr.setThreshold(1e-11);
        const int rank = qr.rank();
        constraints.resize(rank, local_size);
        const auto permutation = qr.colsPermutation().indices();
        for (int row = 0; row < rank; ++row)
            constraints.row(row) = active_constraints.row(permutation(row));
    } else {
        constraints.resize(0, local_size);
    }

    const int constraint_count = static_cast<int>(constraints.rows());
    std::vector<ComplexTriplet> saddle_triplets;
    saddle_triplets.reserve(local_operator.nonZeros()
        + 2 * static_cast<size_t>(constraint_count * local_size));
    for (int col = 0; col < local_operator.outerSize(); ++col) {
        for (ComplexSparseMatrix::InnerIterator it(local_operator, col); it; ++it)
            saddle_triplets.emplace_back(it.row(), it.col(), it.value());
    }
    for (int row = 0; row < constraint_count; ++row) {
        for (int col = 0; col < local_size; ++col) {
            const double value = constraints(row, col);
            if (value == 0.0) continue;
            saddle_triplets.emplace_back(col, local_size + row, value);
            saddle_triplets.emplace_back(local_size + row, col, value);
        }
    }

    ComplexSparseMatrix saddle(local_size + constraint_count, local_size + constraint_count);
    saddle.setFromTriplets(saddle_triplets.begin(), saddle_triplets.end());
    saddle.makeCompressed();
    ComplexMatrix saddle_rhs = ComplexMatrix::Zero(local_size + constraint_count, 3);
    saddle_rhs.topRows(local_size) = rhs;

    bool symbolic_reused = false;
    ComplexSparseLu &solver = sparse_lu_cache.factorize(saddle, symbolic_reused);
    const ComplexMatrix solution = solver.solve(saddle_rhs);
    if (solver.info() != Eigen::Success || !solution.allFinite())
        throw std::runtime_error("Helmholtz corrector saddle solve failed");
    result.symbolic_reused = symbolic_reused;

    const ComplexMatrix corrector = solution.topRows(local_size);
    const ComplexMatrix multipliers = solution.bottomRows(constraint_count);
    ComplexMatrix equation_residual = local_operator * corrector - rhs;
    if (constraint_count > 0)
        equation_residual.noalias() += constraints.transpose().cast<Complex>() * multipliers;
    const double rhs_scale = std::max(1.0, rhs.norm());
    result.primal_residual = equation_residual.norm() / rhs_scale;
    result.constraint_residual = constraint_count > 0
        ? (constraints.cast<Complex>() * corrector).norm() / rhs_scale
        : 0.0;

    const ComplexMatrix adjoint_corrector = corrector.conjugate();
    const ComplexMatrix adjoint_multipliers = multipliers.conjugate();
    ComplexMatrix adjoint_residual = local_operator.adjoint() * adjoint_corrector
                                   - rhs.conjugate();
    if (constraint_count > 0) {
        adjoint_residual.noalias() += constraints.transpose().cast<Complex>()
                                   * adjoint_multipliers;
    }
    result.adjoint_residual = adjoint_residual.norm() / rhs_scale;

    result.primal.reserve(static_cast<size_t>(local_size) * 3);
    result.adjoint.reserve(static_cast<size_t>(local_size) * 3);
    for (int row = 0; row < local_size; ++row) {
        for (int col = 0; col < 3; ++col) {
            const Complex value = corrector(row, col);
            if (std::abs(value) <= 1e-14) continue;
            result.primal.push_back({local_vertices[row], col, value});
            result.adjoint.push_back({local_vertices[row], col, std::conj(value)});
        }
    }
    return result;
}

} // namespace

HelmholtzCorrectorResult build_helmholtz_correctors(
    const TriMesh &coarse,
    const TriMesh &fine,
    const Eigen::SparseMatrix<double> &fine_element_prolongation,
    const Eigen::SparseMatrix<double> &fine_dg_prolongation,
    const Eigen::SparseMatrix<double> &quasi_interpolation,
    const Eigen::SparseMatrix<double> &patches,
    const HelmholtzElementBlocks &element_blocks) {
    const int coarse_element_count = static_cast<int>(coarse.elems.size());
    if (patches.rows() != coarse_element_count || patches.cols() != coarse_element_count)
        throw std::invalid_argument("Helmholtz patch matrix has the wrong dimensions");
    if (fine_dg_prolongation.rows() != 3 * static_cast<int>(fine.elems.size())
        || fine_dg_prolongation.cols() != 3 * coarse_element_count)
        throw std::invalid_argument("Helmholtz DG prolongation has the wrong dimensions");
    if (quasi_interpolation.rows() != static_cast<int>(coarse.nodes.size())
        || quasi_interpolation.cols() != static_cast<int>(fine.nodes.size()))
        throw std::invalid_argument("Helmholtz quasi-interpolation has the wrong dimensions");
    if (element_blocks.size() != fine.elems.size())
        throw std::invalid_argument("Helmholtz element block count must match fine elements");

    const auto children = fine_element_children(fine_element_prolongation, coarse_element_count);
    const auto incidence = node_incidence(fine);
    const auto boundary_edges = edge_counts(fine);
    std::vector<PatchSolveResult> patch_results(coarse_element_count);
    std::vector<int> target_order(coarse_element_count);
    std::iota(target_order.begin(), target_order.end(), 0);
    std::vector<size_t> patch_cost(coarse_element_count, 0);
    for (int target = 0; target < coarse_element_count; ++target) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(patches, target); it; ++it) {
            if (it.value() != 0.0) patch_cost[target] += children[it.row()].size();
        }
    }
    std::stable_sort(target_order.begin(), target_order.end(), [&](int lhs, int rhs) {
        if (patch_cost[lhs] != patch_cost[rhs]) return patch_cost[lhs] > patch_cost[rhs];
        return lhs < rhs;
    });

    std::atomic<bool> failed{false};
    std::exception_ptr first_exception;
    std::mutex exception_mutex;
    int used_threads = 1;
    auto run_target = [&](int ordered_index) {
        if (failed.load(std::memory_order_relaxed)) return;
        const int target = target_order[ordered_index];
        try {
            patch_results[target] = solve_patch(
                target, fine, static_cast<int>(coarse.nodes.size()),
                children, incidence, boundary_edges,
                fine_dg_prolongation, quasi_interpolation, patches, element_blocks);
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
        for (int ordered_index = 0; ordered_index < coarse_element_count; ++ordered_index)
            run_target(ordered_index);
    }
#else
    for (int ordered_index = 0; ordered_index < coarse_element_count; ++ordered_index)
        run_target(ordered_index);
#endif
    if (first_exception) std::rethrow_exception(first_exception);

    HelmholtzCorrectorResult result;
    result.primal.resize(coarse_element_count);
    result.adjoint.resize(coarse_element_count);
    result.diagnostics.patch_count = coarse_element_count;
    result.diagnostics.parallel_threads = used_threads;
    for (int target = 0; target < coarse_element_count; ++target) {
        result.primal[target] = std::move(patch_results[target].primal);
        result.adjoint[target] = std::move(patch_results[target].adjoint);
        result.diagnostics.max_primal_residual = std::max(
            result.diagnostics.max_primal_residual, patch_results[target].primal_residual);
        result.diagnostics.max_adjoint_residual = std::max(
            result.diagnostics.max_adjoint_residual, patch_results[target].adjoint_residual);
        result.diagnostics.max_constraint_residual = std::max(
            result.diagnostics.max_constraint_residual, patch_results[target].constraint_residual);
        if (patch_results[target].touches_physical_boundary)
            ++result.diagnostics.patches_touching_physical_boundary;
        if (patch_results[target].symbolic_reused)
            ++result.diagnostics.symbolic_reuses;
        else
            ++result.diagnostics.symbolic_analyses;
    }
    return result;
}

ComplexSparseMatrix build_helmholtz_corrector_matrix(
    const TriMesh &coarse,
    int fine_node_count,
    const std::vector<HelmholtzElementCorrector> &correctors) {
    if (correctors.size() != coarse.elems.size())
        throw std::invalid_argument("Helmholtz corrector count must match coarse elements");
    std::vector<ComplexTriplet> triplets;
    size_t nonzeros = 0;
    for (const auto &element : correctors) nonzeros += element.size();
    triplets.reserve(nonzeros);
    for (int element = 0; element < static_cast<int>(coarse.elems.size()); ++element) {
        for (const auto &entry : correctors[element]) {
            triplets.emplace_back(
                entry.row,
                coarse.elems[element][entry.local_coarse_vertex],
                entry.value);
        }
    }
    ComplexSparseMatrix matrix(fine_node_count, static_cast<int>(coarse.nodes.size()));
    matrix.setFromTriplets(triplets.begin(), triplets.end());
    return matrix;
}

ComplexSparseMatrix build_helmholtz_corrected_basis(
    const Eigen::SparseMatrix<double> &coarse_to_fine,
    const TriMesh &coarse,
    int fine_node_count,
    const std::vector<HelmholtzElementCorrector> &correctors) {
    ComplexSparseMatrix basis = coarse_to_fine.cast<Complex>();
    basis -= build_helmholtz_corrector_matrix(coarse, fine_node_count, correctors);
    basis.prune(Complex(0.0, 0.0), 1e-14);
    basis.makeCompressed();
    return basis;
}

} // namespace lod2d::helmholtz

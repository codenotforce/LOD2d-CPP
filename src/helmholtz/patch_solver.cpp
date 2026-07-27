#include "helmholtz/patch_solver.h"
#include "helmholtz/patch_multigrid.h"
#include "helmholtz/shifted_laplacian.h"


#include <Eigen/LU>
#include <Eigen/SparseLU>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <memory>
#include <string>
#include <cstdint>
#include <vector>

namespace lod2d::helmholtz {
namespace {

using ComplexSparseLu = Eigen::SparseLU<ComplexSparseMatrix>;

struct SparseLuPatternCache {
    ComplexSparseLu solver;
    std::vector<int> outer_indices;
    std::vector<int> inner_indices;
    int rows = -1;
    int cols = -1;
    std::vector<Complex> values;
    std::uint64_t last_use = 0;

    bool matches(const ComplexSparseMatrix &matrix) const {
        if (rows != matrix.rows() || cols != matrix.cols()) return false;
        if (outer_indices.size() != static_cast<std::size_t>(matrix.outerSize() + 1)
            || inner_indices.size() != static_cast<std::size_t>(matrix.nonZeros()))
            return false;
        return std::equal(
                   outer_indices.begin(), outer_indices.end(), matrix.outerIndexPtr())
            && std::equal(
                   inner_indices.begin(), inner_indices.end(), matrix.innerIndexPtr());
    }

    bool values_match(const ComplexSparseMatrix &matrix) const {
        return matches(matrix)
            && values.size() == static_cast<std::size_t>(matrix.nonZeros())
            && std::equal(values.begin(), values.end(), matrix.valuePtr());
    }

    ComplexSparseLu &factorize(
        const ComplexSparseMatrix &matrix,
        bool reuse_identical,
        bool &symbolic_reused,
        bool &factorization_reused,
        const char *description) {
        symbolic_reused = matches(matrix);
        factorization_reused = reuse_identical && values_match(matrix);
        if (!symbolic_reused) {
            solver.analyzePattern(matrix);
            rows = matrix.rows();
            cols = matrix.cols();
            outer_indices.assign(
                matrix.outerIndexPtr(), matrix.outerIndexPtr() + matrix.outerSize() + 1);
            inner_indices.assign(
                matrix.innerIndexPtr(), matrix.innerIndexPtr() + matrix.nonZeros());
        }
        if (factorization_reused) return solver;
        solver.factorize(matrix);
        if (solver.info() != Eigen::Success)
            throw std::runtime_error(std::string(description) + " factorization failed");
        values.assign(matrix.valuePtr(), matrix.valuePtr() + matrix.nonZeros());
        return solver;
    }
};

struct SparseLuCachePool {
    std::vector<std::unique_ptr<SparseLuPatternCache>> entries;
    std::uint64_t clock = 0;

    ComplexSparseLu &factorize(
        const ComplexSparseMatrix &matrix,
        int requested_slots,
        bool reuse_identical,
        bool &symbolic_reused,
        bool &factorization_reused,
        const char *description) {
        const int slots = std::max(1, requested_slots);
        if (static_cast<int>(entries.size()) > slots)
            entries.resize(slots);

        SparseLuPatternCache *selected = nullptr;
        if (reuse_identical) {
            for (const auto &entry : entries) {
                if (entry->values_match(matrix)) {
                    selected = entry.get();
                    break;
                }
            }
        }
        if (!selected && !reuse_identical) {
            for (const auto &entry : entries) {
                if (entry->matches(matrix)) {
                    selected = entry.get();
                    break;
                }
            }
        }
        if (!selected && static_cast<int>(entries.size()) < slots) {
            entries.push_back(std::make_unique<SparseLuPatternCache>());
            selected = entries.back().get();
        }
        if (!selected) {
            selected = entries.front().get();
            for (const auto &entry : entries) {
                if (entry->last_use < selected->last_use)
                    selected = entry.get();
            }
        }
        selected->last_use = ++clock;
        return selected->factorize(
            matrix, reuse_identical, symbolic_reused,
            factorization_reused, description);
    }
};

thread_local SparseLuCachePool saddle_cache;
thread_local SparseLuCachePool helmholtz_cache;
thread_local SparseLuCachePool shifted_cache;

ComplexSparseMatrix build_saddle(const HelmholtzPatchSystem &system) {
    const int n = system.helmholtz.rows();
    const int m = static_cast<int>(system.constraints.rows());
    std::vector<ComplexTriplet> triplets;
    triplets.reserve(system.helmholtz.nonZeros()
        + 2 * static_cast<std::size_t>(m * n));
    for (int col = 0; col < system.helmholtz.outerSize(); ++col) {
        for (ComplexSparseMatrix::InnerIterator it(system.helmholtz, col); it; ++it)
            triplets.emplace_back(it.row(), it.col(), it.value());
    }
    for (int row = 0; row < m; ++row) {
        for (int col = 0; col < n; ++col) {
            const double value = system.constraints(row, col);
            if (value == 0.0) continue;
            triplets.emplace_back(col, n + row, value);
            triplets.emplace_back(n + row, col, value);
        }
    }
    ComplexSparseMatrix saddle(n + m, n + m);
    saddle.setFromTriplets(triplets.begin(), triplets.end());
    saddle.makeCompressed();
    return saddle;
}

HelmholtzPatchSolveResult solve_direct_saddle(
    const HelmholtzPatchSystem &system,
    const HelmholtzPatchSolverConfig &config) {
    const int n = system.helmholtz.rows();
    const int m = static_cast<int>(system.constraints.rows());
    const ComplexSparseMatrix saddle = build_saddle(system);
    ComplexMatrix saddle_rhs = ComplexMatrix::Zero(n + m, system.rhs.cols());
    saddle_rhs.topRows(n) = system.rhs;

    bool symbolic_reused = false;
    bool factorization_reused = false;
    ComplexSparseLu &solver = saddle_cache.factorize(
        saddle, config.symbolic_cache_slots,
        config.reuse_identical_factorization,
        symbolic_reused, factorization_reused,
        "Helmholtz corrector saddle");
    const ComplexMatrix solution = solver.solve(saddle_rhs);
    if (solver.info() != Eigen::Success || !solution.allFinite())
        throw std::runtime_error("Helmholtz corrector saddle solve failed");

    HelmholtzPatchSolveResult result;
    result.corrector = solution.topRows(n);
    result.multipliers = solution.bottomRows(m);
    result.diagnostics.symbolic_reused = symbolic_reused;
    result.diagnostics.factorization_reused = factorization_reused;
    return result;
}

struct SchurInputs {
    ComplexMatrix z;
    ComplexMatrix y;
    HelmholtzPatchSolveDiagnostics diagnostics;
};

SchurInputs solve_direct_helmholtz_blocks(
    const HelmholtzPatchSystem &system,
    const HelmholtzPatchSolverConfig &config) {
    const int n = system.helmholtz.rows();
    const int m = static_cast<int>(system.constraints.rows());
    ComplexMatrix combined(n, m + system.rhs.cols());
    if (m > 0)
        combined.leftCols(m) = system.constraints.transpose().cast<Complex>();
    combined.rightCols(system.rhs.cols()) = system.rhs;

    bool symbolic_reused = false;
    bool factorization_reused = false;
    ComplexSparseLu &solver = helmholtz_cache.factorize(
        system.helmholtz, config.symbolic_cache_slots,
        config.reuse_identical_factorization,
        symbolic_reused, factorization_reused,
        "Helmholtz patch A block");
    const ComplexMatrix solved = solver.solve(combined);
    if (solver.info() != Eigen::Success || !solved.allFinite())
        throw std::runtime_error("Helmholtz patch A-block solve failed");

    SchurInputs result;
    result.z = solved.leftCols(m);
    result.y = solved.rightCols(system.rhs.cols());
    result.diagnostics.symbolic_reused = symbolic_reused;
    result.diagnostics.factorization_reused = factorization_reused;
    return result;
}


SchurInputs solve_shifted_gmres_blocks(
    const HelmholtzPatchSystem &system,
    const HelmholtzPatchSolverConfig &config) {

    const int n = system.helmholtz.rows();
    const int m = static_cast<int>(system.constraints.rows());
    ComplexMatrix combined(n, m + system.rhs.cols());
    if (m > 0)
        combined.leftCols(m) = system.constraints.transpose().cast<Complex>();
    combined.rightCols(system.rhs.cols()) = system.rhs;

    const double epsilon = helmholtz_shift_epsilon(system, config.shifted);
    SchurInputs result;
    ComplexSparseLu *sparse_preconditioner = nullptr;
    std::unique_ptr<HelmholtzPatchVcycle> vcycle;
    if (config.shifted.inverse == HelmholtzShiftedInverseKind::SparseLu) {
        const ComplexSparseMatrix shifted =
            build_shifted_helmholtz_operator(system, epsilon);
        bool symbolic_reused = false;
        bool factorization_reused = false;
        sparse_preconditioner = &shifted_cache.factorize(
            shifted, config.symbolic_cache_slots,
            config.reuse_identical_factorization,
            symbolic_reused, factorization_reused,
            "shifted Helmholtz patch");
        result.diagnostics.symbolic_reused = symbolic_reused;
        result.diagnostics.factorization_reused = factorization_reused;
    } else if (config.shifted.inverse
               == HelmholtzShiftedInverseKind::GeometricVcycle) {
        if (system.geometric_prolongations.empty())
            throw std::runtime_error(
                "geometric V-cycle requires a nontrivial uniform NVB hierarchy");
        vcycle = std::make_unique<HelmholtzPatchVcycle>(
            system, epsilon, config.shifted.pre_smooth,
            config.shifted.post_smooth, config.shifted.coarse_max_dofs,
            config.shifted.jacobi_weight);
        result.diagnostics.vcycle_levels = vcycle->levels();
        result.diagnostics.vcycle_coarse_dofs = vcycle->coarse_dofs();
        result.diagnostics.vcycle_finest_dofs = vcycle->finest_dofs();
        result.diagnostics.vcycle_relative_residual =
            vcycle->relative_residual(combined.col(0));
    }
    ComplexMatrix solved(n, combined.cols());
    for (int column = 0; column < combined.cols(); ++column) {
        const auto apply_operator = [&](const solver::ComplexVector &vector) {
            return solver::ComplexVector(system.helmholtz * vector);
        };
        const auto apply_preconditioner = [&](const solver::ComplexVector &vector) {
            if (vcycle) return vcycle->apply(vector);
            if (sparse_preconditioner) {
                solver::ComplexVector value = sparse_preconditioner->solve(vector);
                if (sparse_preconditioner->info() != Eigen::Success)
                    throw std::runtime_error(
                        "shifted Helmholtz preconditioner solve failed");
                return solver::ComplexVector(value);
            }
            return solver::ComplexVector(vector);
        };
        const solver::RightGmresResult gmres =
            solver::solve_right_preconditioned_gmres(
                n, apply_operator, apply_preconditioner,
                combined.col(column), config.gmres);
        if (!gmres.converged) {
            throw std::runtime_error(
                "right-preconditioned Helmholtz GMRES failed: " + gmres.message);
        }
        solved.col(column) = gmres.solution;
        ++result.diagnostics.gmres_right_hand_sides;
        result.diagnostics.gmres_total_iterations += gmres.iterations;
        result.diagnostics.gmres_max_iterations = std::max(
            result.diagnostics.gmres_max_iterations, gmres.iterations);
        result.diagnostics.gmres_restarts += gmres.restarts;
        result.diagnostics.max_gmres_residual = std::max(
            result.diagnostics.max_gmres_residual, gmres.relative_residual);
    }
    result.z = solved.leftCols(m);
    result.y = solved.rightCols(system.rhs.cols());
    return result;
}

HelmholtzPatchSolveResult recover_schur_solution(
    const HelmholtzPatchSystem &system,
    SchurInputs inputs) {
    const int m = static_cast<int>(system.constraints.rows());
    HelmholtzPatchSolveResult result;
    result.diagnostics = inputs.diagnostics;
    if (m == 0) {
        result.corrector = std::move(inputs.y);
        result.multipliers.resize(0, system.rhs.cols());
        return result;
    }

    const ComplexMatrix constraints = system.constraints.cast<Complex>();
    const ComplexMatrix schur = constraints * inputs.z;
    const ComplexMatrix schur_rhs = constraints * inputs.y;
    Eigen::FullPivLU<ComplexMatrix> factorization(schur);
    result.diagnostics.schur_rcond = factorization.rcond();
    if (!factorization.isInvertible())
        throw std::runtime_error("Helmholtz patch Schur complement is singular");
    result.multipliers = factorization.solve(schur_rhs);
    if (!result.multipliers.allFinite())
        throw std::runtime_error("Helmholtz patch Schur solve returned non-finite values");
    result.corrector = inputs.y - inputs.z * result.multipliers;
    result.diagnostics.schur_residual =
        (schur * result.multipliers - schur_rhs).norm()
        / std::max(1.0, schur_rhs.norm());
    return result;
}

void compute_final_residuals(
    const HelmholtzPatchSystem &system,
    HelmholtzPatchSolveResult &result) {
    ComplexMatrix equation_residual =
        system.helmholtz * result.corrector - system.rhs;
    if (system.constraints.rows() > 0) {
        equation_residual.noalias() +=
            system.constraints.transpose().cast<Complex>() * result.multipliers;
    }
    const double rhs_scale = std::max(1.0, system.rhs.norm());
    result.diagnostics.primal_residual = equation_residual.norm() / rhs_scale;
    result.diagnostics.constraint_residual = system.constraints.rows() > 0
        ? (system.constraints.cast<Complex>() * result.corrector).norm() / rhs_scale
        : 0.0;

    ComplexMatrix adjoint_residual = system.helmholtz.adjoint()
        * result.corrector.conjugate() - system.rhs.conjugate();
    if (system.constraints.rows() > 0) {
        adjoint_residual.noalias() +=
            system.constraints.transpose().cast<Complex>()
            * result.multipliers.conjugate();
    }
    result.diagnostics.adjoint_residual = adjoint_residual.norm() / rhs_scale;
}

} // namespace

double helmholtz_shift_epsilon(
    const HelmholtzPatchSystem &system,
    const HelmholtzShiftedLaplacianConfig &config) {
    if (!std::isfinite(config.alpha)
        || !std::isfinite(config.absolute_epsilon))
        throw std::invalid_argument("Helmholtz shift parameters must be finite");
    double epsilon = 0.0;
    switch (config.rule) {
    case HelmholtzShiftRule::KappaSquared:
        epsilon = config.alpha * system.wavenumber * system.wavenumber;
        break;
    case HelmholtzShiftRule::PatchScaled:
        if (!(system.diameter > 0.0))
            throw std::invalid_argument("Helmholtz patch diameter must be positive");
        epsilon = config.alpha * system.wavenumber / system.diameter;
        break;
    case HelmholtzShiftRule::Absolute:
        epsilon = config.absolute_epsilon;
        break;
    }
    if (!(epsilon >= 0.0) || !std::isfinite(epsilon))
        throw std::invalid_argument("Helmholtz shift epsilon must be finite and nonnegative");
    return epsilon;
}

HelmholtzPatchSolveResult solve_helmholtz_patch(
    const HelmholtzPatchSystem &system,
    const HelmholtzPatchSolverConfig &config) {
    if (system.helmholtz.rows() <= 0
        || system.helmholtz.rows() != system.helmholtz.cols())
        throw std::invalid_argument("Helmholtz patch operator must be square and nonempty");
    if (system.rhs.rows() != system.helmholtz.rows())
        throw std::invalid_argument("Helmholtz patch right-hand side has the wrong size");
    if (system.constraints.cols() != system.helmholtz.cols())
        throw std::invalid_argument("Helmholtz patch constraints have the wrong size");
    if (config.symbolic_cache_slots <= 0)
        throw std::invalid_argument(
            "Helmholtz symbolic cache slots must be positive");

    try {
        HelmholtzPatchSolveResult result;
        switch (config.kind) {
        case HelmholtzPatchSolverKind::DirectSaddle:
            result = solve_direct_saddle(system, config);
            break;
        case HelmholtzPatchSolverKind::DirectSchur:
            result = recover_schur_solution(
                system, solve_direct_helmholtz_blocks(system, config));
            break;
        case HelmholtzPatchSolverKind::ShiftedGmres:
            result = recover_schur_solution(
                system, solve_shifted_gmres_blocks(system, config));
            break;
        }
        compute_final_residuals(system, result);
        return result;
    } catch (...) {
        if (!config.fallback_to_direct
            || config.kind == HelmholtzPatchSolverKind::DirectSaddle)
            throw;
        HelmholtzPatchSolveResult result = solve_direct_saddle(system, config);
        result.diagnostics.direct_fallback = true;
        compute_final_residuals(system, result);
        return result;
    }
}

} // namespace lod2d::helmholtz

#include "helmholtz/hp_model.h"

#include <Eigen/SparseLU>
#include <algorithm>
#include <atomic>
#include <exception>
#include <iterator>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lod2d::helmholtz {

struct HelmholtzHpLodModel::Implementation {
    HelmholtzHpProblemConfig config;
    HelmholtzProblemData problem;
    std::unique_ptr<HpTriSpace> fine_space;
    HelmholtzHpInterpolation interpolation;
    HelmholtzHpOperators operators;
    HelmholtzHpCorrectorDiagnostics diagnostics;
    ComplexSparseMatrix corrector;
    ComplexSparseMatrix adjoint_corrector;
    ComplexSparseMatrix trial_basis;
    ComplexSparseMatrix test_basis;
    ComplexSparseMatrix coarse_operator;
    Eigen::SparseLU<ComplexSparseMatrix> coarse_solver;
    Eigen::SparseLU<ComplexSparseMatrix> fine_solver;
    bool fine_factorized = false;
};

HelmholtzHpLodModel::HelmholtzHpLodModel() = default;
HelmholtzHpLodModel::~HelmholtzHpLodModel() = default;
HelmholtzHpLodModel::HelmholtzHpLodModel(
    HelmholtzHpLodModel &&) noexcept = default;
HelmholtzHpLodModel &HelmholtzHpLodModel::operator=(
    HelmholtzHpLodModel &&) noexcept = default;

HelmholtzHpLodModel HelmholtzHpLodModel::build(
    const HelmholtzHpProblemConfig &config) {
    if (config.degree < 1 || config.degree > 3)
        throw std::invalid_argument("hp LOD degree must be 1, 2, or 3");
    if (config.h < config.H)
        throw std::invalid_argument("hp LOD requires h >= H");
    if (!(config.wavenumber > 0.0))
        throw std::invalid_argument("hp LOD wavenumber must be positive");

    HelmholtzHpLodModel model;
    model.implementation_ = std::make_unique<Implementation>();
    Implementation &data = *model.implementation_;
    data.config = config;
    if (data.config.initial_mesh.nodes.empty())
        data.config.initial_mesh = make_helmholtz_unit_square_mesh();
    data.problem = build_helmholtz_problem_data(
        data.config.initial_mesh, data.config.H, data.config.h,
        data.config.ell);
    data.fine_space =
        std::make_unique<HpTriSpace>(data.problem.fine, data.config.degree);
    data.interpolation = build_helmholtz_hp_interpolation(
        data.problem.coarse, *data.fine_space,
        data.problem.fine_element_prolongation);
    data.operators = assemble_helmholtz_hp_operators(
        *data.fine_space, data.config.wavenumber, {}, {},
        data.config.boundary_beta);

    HelmholtzHpPatchAssembler assembler(
        data.problem.coarse, *data.fine_space,
        data.problem.fine_element_prolongation, data.problem.patches,
        data.interpolation, data.operators);
    const int patch_count = assembler.patch_count();
    std::vector<int> target_order(patch_count);
    std::iota(target_order.begin(), target_order.end(), 0);
    std::stable_sort(
        target_order.begin(), target_order.end(), [&](int lhs, int rhs) {
            const std::size_t lhs_cost = assembler.patch_cost(lhs);
            const std::size_t rhs_cost = assembler.patch_cost(rhs);
            return lhs_cost == rhs_cost ? lhs < rhs : lhs_cost > rhs_cost;
        });

    std::vector<ComplexTriplet> primal_triplets;
    std::vector<ComplexTriplet> adjoint_triplets;
    std::atomic<bool> failed{false};
    std::atomic<int> completed{0};
    std::exception_ptr first_exception;
    std::mutex exception_mutex;
    std::mutex progress_mutex;
    int used_threads = 1;

    const auto solve_target = [&] (
        int target,
        std::vector<ComplexTriplet> &local_primal,
        std::vector<ComplexTriplet> &local_adjoint,
        HelmholtzHpCorrectorDiagnostics &local_diagnostics) {
        if (failed.load(std::memory_order_relaxed)) return;
        try {
            const HelmholtzHpPatchSolveResult patch =
                assembler.solve(target, data.config.patch_solver);
            const Triangle &coarse_triangle =
                data.problem.coarse.elems[target];
            for (int row = 0;
                 row < static_cast<int>(patch.system.local_vertices.size());
                 ++row) {
                const int global_dof = patch.system.local_vertices[row];
                for (int local = 0; local < 3; ++local) {
                    const Complex primal =
                        patch.primal.corrector(row, local);
                    const Complex adjoint =
                        patch.adjoint_corrector(row, local);
                    if (std::abs(primal) > 1e-14)
                        local_primal.emplace_back(
                            global_dof, coarse_triangle[local], primal);
                    if (std::abs(adjoint) > 1e-14)
                        local_adjoint.emplace_back(
                            global_dof, coarse_triangle[local], adjoint);
                }
            }
            local_diagnostics.max_primal_residual = std::max(
                local_diagnostics.max_primal_residual,
                patch.primal.diagnostics.primal_residual);
            local_diagnostics.max_adjoint_residual = std::max(
                local_diagnostics.max_adjoint_residual,
                patch.primal.diagnostics.adjoint_residual);
            local_diagnostics.max_constraint_residual = std::max(
                local_diagnostics.max_constraint_residual,
                patch.primal.diagnostics.constraint_residual);
            local_diagnostics.max_schur_residual = std::max(
                local_diagnostics.max_schur_residual,
                patch.primal.diagnostics.schur_residual);
            local_diagnostics.min_schur_rcond = std::min(
                local_diagnostics.min_schur_rcond,
                patch.primal.diagnostics.schur_rcond);
            local_diagnostics.direct_fallback_count +=
                patch.primal.diagnostics.direct_fallback ? 1 : 0;
            ++local_diagnostics.patch_count;
            local_diagnostics.boundary_patch_count +=
                patch.system.touches_physical_boundary ? 1 : 0;

            const int done =
                completed.fetch_add(1, std::memory_order_relaxed) + 1;
            const bool report = data.config.corrector_progress
                && (done == patch_count
                    || (data.config.progress_interval > 0
                        && done % data.config.progress_interval == 0));
            if (report) {
                std::lock_guard<std::mutex> lock(progress_mutex);
                data.config.corrector_progress(done, patch_count);
            }
        } catch (...) {
            failed.store(true, std::memory_order_relaxed);
            std::lock_guard<std::mutex> lock(exception_mutex);
            if (!first_exception) first_exception = std::current_exception();
        }
    };

    const auto merge_thread_result = [&] (
        std::vector<ComplexTriplet> &local_primal,
        std::vector<ComplexTriplet> &local_adjoint,
        const HelmholtzHpCorrectorDiagnostics &local_diagnostics) {
        primal_triplets.insert(
            primal_triplets.end(),
            std::make_move_iterator(local_primal.begin()),
            std::make_move_iterator(local_primal.end()));
        adjoint_triplets.insert(
            adjoint_triplets.end(),
            std::make_move_iterator(local_adjoint.begin()),
            std::make_move_iterator(local_adjoint.end()));
        data.diagnostics.max_primal_residual = std::max(
            data.diagnostics.max_primal_residual,
            local_diagnostics.max_primal_residual);
        data.diagnostics.max_adjoint_residual = std::max(
            data.diagnostics.max_adjoint_residual,
            local_diagnostics.max_adjoint_residual);
        data.diagnostics.max_constraint_residual = std::max(
            data.diagnostics.max_constraint_residual,
            local_diagnostics.max_constraint_residual);
        data.diagnostics.max_schur_residual = std::max(
            data.diagnostics.max_schur_residual,
            local_diagnostics.max_schur_residual);
        data.diagnostics.min_schur_rcond = std::min(
            data.diagnostics.min_schur_rcond,
            local_diagnostics.min_schur_rcond);
        data.diagnostics.direct_fallback_count +=
            local_diagnostics.direct_fallback_count;
        data.diagnostics.patch_count += local_diagnostics.patch_count;
        data.diagnostics.boundary_patch_count +=
            local_diagnostics.boundary_patch_count;
    };

#ifdef _OPENMP
    const int requested_threads = std::max(1, data.config.corrector_threads);
    #pragma omp parallel num_threads(requested_threads)
    {
        std::vector<ComplexTriplet> local_primal;
        std::vector<ComplexTriplet> local_adjoint;
        HelmholtzHpCorrectorDiagnostics local_diagnostics;
        #pragma omp single
        used_threads = omp_get_num_threads();
        #pragma omp for schedule(dynamic, 1)
        for (int ordered = 0; ordered < patch_count; ++ordered)
            solve_target(
                target_order[ordered], local_primal, local_adjoint,
                local_diagnostics);
        #pragma omp critical(hp_corrector_merge)
        merge_thread_result(
            local_primal, local_adjoint, local_diagnostics);
    }
#else
    std::vector<ComplexTriplet> local_primal;
    std::vector<ComplexTriplet> local_adjoint;
    HelmholtzHpCorrectorDiagnostics local_diagnostics;
    for (int target : target_order)
        solve_target(
            target, local_primal, local_adjoint, local_diagnostics);
    merge_thread_result(local_primal, local_adjoint, local_diagnostics);
#endif
    data.diagnostics.parallel_threads = used_threads;
    if (first_exception) std::rethrow_exception(first_exception);
    if (data.diagnostics.max_primal_residual > 1e-8
        || data.diagnostics.max_adjoint_residual > 1e-8
        || data.diagnostics.max_constraint_residual > 1e-8
        || data.diagnostics.max_schur_residual > 1e-8)
        throw std::runtime_error(
            "hp LOD corrector residual exceeded the correctness threshold");

    const int hp_dofs = data.fine_space->dof_count();
    const int coarse_dofs =
        static_cast<int>(data.problem.coarse.nodes.size());
    data.corrector.resize(hp_dofs, coarse_dofs);
    data.corrector.setFromTriplets(
        primal_triplets.begin(), primal_triplets.end());
    data.adjoint_corrector.resize(hp_dofs, coarse_dofs);
    data.adjoint_corrector.setFromTriplets(
        adjoint_triplets.begin(), adjoint_triplets.end());
    const ComplexSparseMatrix injection =
        data.interpolation.coarse_injection.cast<Complex>();
    const ComplexSparseMatrix corrected_trial =
        injection - data.corrector;
    const ComplexSparseMatrix corrected_test =
        injection - data.adjoint_corrector;
    data.test_basis = corrected_test;
    data.trial_basis = data.config.mode == HelmholtzPetrovMode::TwoSided
        ? corrected_trial : injection;
    data.trial_basis.prune(Complex(0.0), 1e-14);
    data.test_basis.prune(Complex(0.0), 1e-14);
    data.trial_basis.makeCompressed();
    data.test_basis.makeCompressed();

    data.coarse_operator =
        data.test_basis.adjoint() * data.operators.system
        * data.trial_basis;
    data.coarse_operator.prune(Complex(0.0), 1e-14);
    data.coarse_operator.makeCompressed();
    data.coarse_solver.analyzePattern(data.coarse_operator);
    data.coarse_solver.factorize(data.coarse_operator);
    if (data.coarse_solver.info() != Eigen::Success)
        throw std::runtime_error("hp LOD coarse factorization failed");
    return model;
}

HelmholtzLodSolution HelmholtzHpLodModel::solve_load(
    const ComplexVector &load) const {
    if (!implementation_)
        throw std::logic_error("hp LOD model has not been built");
    Implementation &data = *implementation_;
    if (load.size() != data.fine_space->dof_count())
        throw std::invalid_argument("hp LOD load has the wrong size");
    const ComplexVector coarse_rhs = data.test_basis.adjoint() * load;
    HelmholtzLodSolution result;
    result.coarse_coefficients = data.coarse_solver.solve(coarse_rhs);
    if (data.coarse_solver.info() != Eigen::Success
        || !result.coarse_coefficients.allFinite())
        throw std::runtime_error("hp LOD coarse solve failed");
    result.fine_values =
        data.trial_basis * result.coarse_coefficients;
    const ComplexVector residual = data.test_basis.adjoint()
        * (data.operators.system * result.fine_values - load);
    result.petrov_residual =
        residual.norm() / std::max(1.0, coarse_rhs.norm());
    return result;
}

HelmholtzLodSolution HelmholtzHpLodModel::solve_source(
    const ComplexFunction &source) const {
    return solve_load(assemble_helmholtz_hp_load(fine_space(), source));
}

ComplexVector HelmholtzHpLodModel::solve_fine_reference(
    const ComplexVector &load) const {
    if (!implementation_)
        throw std::logic_error("hp LOD model has not been built");
    Implementation &data = *implementation_;
    if (!data.fine_factorized) {
        data.fine_solver.analyzePattern(data.operators.system);
        data.fine_solver.factorize(data.operators.system);
        if (data.fine_solver.info() != Eigen::Success)
            throw std::runtime_error("hp LOD fine factorization failed");
        data.fine_factorized = true;
    }
    ComplexVector solution = data.fine_solver.solve(load);
    if (data.fine_solver.info() != Eigen::Success
        || !solution.allFinite())
        throw std::runtime_error("hp LOD fine solve failed");
    return solution;
}

const HelmholtzHpProblemConfig &HelmholtzHpLodModel::config() const {
    return implementation_->config;
}
const HelmholtzProblemData &HelmholtzHpLodModel::problem() const {
    return implementation_->problem;
}
const HpTriSpace &HelmholtzHpLodModel::fine_space() const {
    return *implementation_->fine_space;
}
const HelmholtzHpInterpolation &
HelmholtzHpLodModel::interpolation() const {
    return implementation_->interpolation;
}
const HelmholtzHpOperators &HelmholtzHpLodModel::operators() const {
    return implementation_->operators;
}
const HelmholtzHpCorrectorDiagnostics &
HelmholtzHpLodModel::corrector_diagnostics() const {
    return implementation_->diagnostics;
}
const ComplexSparseMatrix &HelmholtzHpLodModel::corrector_matrix() const {
    return implementation_->corrector;
}
const ComplexSparseMatrix &
HelmholtzHpLodModel::adjoint_corrector_matrix() const {
    return implementation_->adjoint_corrector;
}
const ComplexSparseMatrix &HelmholtzHpLodModel::trial_basis() const {
    return implementation_->trial_basis;
}
const ComplexSparseMatrix &HelmholtzHpLodModel::test_basis() const {
    return implementation_->test_basis;
}
const ComplexSparseMatrix &HelmholtzHpLodModel::coarse_operator() const {
    return implementation_->coarse_operator;
}

} // namespace lod2d::helmholtz

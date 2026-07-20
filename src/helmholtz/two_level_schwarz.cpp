#include "helmholtz/two_level_schwarz.h"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_map>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lod2d::helmholtz {

namespace {

void hash_word(std::uint64_t word, std::uint64_t &hash) {
    constexpr std::uint64_t prime = 1099511628211ULL;
    for (int byte = 0; byte < 8; ++byte) {
        hash ^= (word >> (8 * byte)) & 0xffULL;
        hash *= prime;
    }
}

std::uint64_t sparse_matrix_hash(const ComplexSparseMatrix &matrix) {
    std::uint64_t hash = 14695981039346656037ULL;
    hash_word(static_cast<std::uint64_t>(matrix.rows()), hash);
    hash_word(static_cast<std::uint64_t>(matrix.cols()), hash);
    hash_word(static_cast<std::uint64_t>(matrix.nonZeros()), hash);
    for (int index = 0; index <= matrix.outerSize(); ++index)
        hash_word(static_cast<std::uint64_t>(matrix.outerIndexPtr()[index]), hash);
    for (int index = 0; index < matrix.nonZeros(); ++index) {
        hash_word(static_cast<std::uint64_t>(matrix.innerIndexPtr()[index]), hash);
        std::uint64_t real_bits = 0;
        std::uint64_t imag_bits = 0;
        const Complex value = matrix.valuePtr()[index];
        const double real = value.real();
        const double imag = value.imag();
        std::memcpy(&real_bits, &real, sizeof(double));
        std::memcpy(&imag_bits, &imag, sizeof(double));
        hash_word(real_bits, hash);
        hash_word(imag_bits, hash);
    }
    return hash;
}

bool sparse_matrices_equal(
    const ComplexSparseMatrix &left,
    const ComplexSparseMatrix &right) {
    if (left.rows() != right.rows() || left.cols() != right.cols()
        || left.nonZeros() != right.nonZeros()
        || left.outerSize() != right.outerSize())
        return false;
    for (int index = 0; index <= left.outerSize(); ++index) {
        if (left.outerIndexPtr()[index] != right.outerIndexPtr()[index])
            return false;
    }
    for (int index = 0; index < left.nonZeros(); ++index) {
        if (left.innerIndexPtr()[index] != right.innerIndexPtr()[index]
            || left.valuePtr()[index] != right.valuePtr()[index])
            return false;
    }
    return true;
}

} // namespace

HelmholtzTwoLevelSchwarzPreconditioner::
HelmholtzTwoLevelSchwarzPreconditioner(
    const HelmholtzLodModel &model,
    HelmholtzTwoLevelSchwarzConfig config)
    : config_(config),
      model_(&model),
      fine_operator_(&model.operators().system),
      dimension_(model.operators().system.rows()) {
    if (dimension_ <= 0 || fine_operator_->cols() != dimension_)
        throw std::invalid_argument("fine Helmholtz Schwarz operator must be square");

    const HelmholtzProblemData &problem = model.problem();
    if (!(config_.artificial_impedance_beta >= 0.0)
        || !std::isfinite(config_.artificial_impedance_beta))
        throw std::invalid_argument(
            "Schwarz artificial impedance beta must be finite and nonnegative");

    HelmholtzSchwarzPatchAssembler assembler(
        problem.fine,
        problem.fine_element_prolongation,
        problem.patches,
        static_cast<int>(problem.coarse.elems.size()),
        problem.fine_hierarchy_meshes,
        problem.fine_node_level_prolongations,
        problem.fine_element_level_prolongations,
        model.operators());

    const bool use_vcycle =
        config_.local_solver.kind
            == HelmholtzSchwarzLocalSolverKind::ShiftedGmres
        && config_.local_solver.shifted_inverse
            == HelmholtzSchwarzShiftedInverseKind::GeometricVcycle;
    if (use_vcycle
        && config_.artificial_boundary
            != HelmholtzSchwarzArtificialBoundary::HomogeneousDirichlet)
        throw std::invalid_argument(
            "Schwarz geometric V-cycle requires Dirichlet artificial boundaries");

    const int patch_count = assembler.patch_count();
    std::vector<ComplexSparseMatrix> local_matrices(patch_count);
    std::vector<Eigen::SparseMatrix<double>> local_masses(patch_count);
    std::vector<std::vector<Eigen::SparseMatrix<double>>>
        local_prolongations(patch_count);
    std::vector<std::vector<int>> core_dofs(patch_count);
    Eigen::VectorXi multiplicity = Eigen::VectorXi::Zero(dimension_);
    diagnostics_.min_local_dofs = std::numeric_limits<int>::max();
    subdomains_.reserve(patch_count);
    for (int target = 0; target < patch_count; ++target) {
        HelmholtzSchwarzLocalSystem system = assembler.assemble(
            target,
            config_.artificial_boundary,
            config_.artificial_impedance_beta,
            config_.local_solver.kind
                == HelmholtzSchwarzLocalSolverKind::ShiftedGmres,
            use_vcycle);
        Subdomain subdomain;
        subdomain.global_dofs = std::move(system.global_dofs);
        core_dofs[target] = std::move(system.core_global_dofs);
        if (subdomain.global_dofs.empty())
            throw std::runtime_error("fine Helmholtz Schwarz patch has no DOFs");
        for (int dof : subdomain.global_dofs) ++multiplicity(dof);

        const int local_dofs = static_cast<int>(subdomain.global_dofs.size());
        diagnostics_.min_local_dofs =
            std::min(diagnostics_.min_local_dofs, local_dofs);
        diagnostics_.max_local_dofs =
            std::max(diagnostics_.max_local_dofs, local_dofs);
        diagnostics_.artificial_boundary_edges +=
            system.artificial_boundary_edges;
        diagnostics_.physical_boundary_edges += system.physical_boundary_edges;
        local_matrices[target] = std::move(system.matrix);
        local_masses[target] = std::move(system.mass);
        local_prolongations[target] =
            std::move(system.geometric_prolongations);
        subdomains_.push_back(std::move(subdomain));
    }

    diagnostics_.uncovered_dofs = 0;
    Eigen::VectorXd weighted_injection = Eigen::VectorXd::Zero(dimension_);
    for (int dof = 0; dof < dimension_; ++dof) {
        if (multiplicity(dof) <= 0) {
            ++diagnostics_.uncovered_dofs;
        } else {
            weighted_injection(dof) =
                1.0 / static_cast<double>(multiplicity(dof));
        }
    }
    if (diagnostics_.uncovered_dofs != 0)
        throw std::runtime_error(
            "fine Helmholtz Schwarz patches do not cover every fine DOF");

    diagnostics_.min_owned_dofs = std::numeric_limits<int>::max();
    diagnostics_.max_owned_dofs = 0;
    if (config_.extension == HelmholtzSchwarzExtension::WeightedOverlap) {
        for (Subdomain &subdomain : subdomains_) {
            subdomain.injection_weights.resize(subdomain.global_dofs.size());
            for (int local = 0;
                 local < static_cast<int>(subdomain.global_dofs.size());
                 ++local) {
                subdomain.injection_weights(local) =
                    weighted_injection(subdomain.global_dofs[local]);
            }
            const int owned = static_cast<int>(subdomain.global_dofs.size());
            diagnostics_.min_owned_dofs =
                std::min(diagnostics_.min_owned_dofs, owned);
            diagnostics_.max_owned_dofs =
                std::max(diagnostics_.max_owned_dofs, owned);
        }
    } else if (config_.extension
               == HelmholtzSchwarzExtension::RestrictedCore) {
        Eigen::VectorXi owner = Eigen::VectorXi::Constant(dimension_, -1);
        for (int target = 0; target < patch_count; ++target) {
            for (int dof : core_dofs[target]) {
                if (owner(dof) < 0) owner(dof) = target;
            }
        }
        for (int dof = 0; dof < dimension_; ++dof) {
            if (owner(dof) < 0)
                throw std::runtime_error(
                    "restricted Schwarz cores do not own every fine DOF");
        }
        for (int target = 0; target < patch_count; ++target) {
            Subdomain &subdomain = subdomains_[target];
            subdomain.injection_weights =
                Eigen::VectorXd::Zero(subdomain.global_dofs.size());
            int owned = 0;
            for (int local = 0;
                 local < static_cast<int>(subdomain.global_dofs.size());
                 ++local) {
                if (owner(subdomain.global_dofs[local]) == target) {
                    subdomain.injection_weights(local) = 1.0;
                    ++owned;
                }
            }
            diagnostics_.min_owned_dofs =
                std::min(diagnostics_.min_owned_dofs, owned);
            diagnostics_.max_owned_dofs =
                std::max(diagnostics_.max_owned_dofs, owned);
        }
    } else {
        throw std::invalid_argument("unknown Schwarz extension");
    }
    if (diagnostics_.min_owned_dofs <= 0)
        throw std::runtime_error("fine Helmholtz Schwarz core is empty");
    if (config_.factorization_reuse
            == HelmholtzSchwarzFactorizationReuse::IdenticalMatrix
        && config_.local_solver.kind
            != HelmholtzSchwarzLocalSolverKind::SparseLu)
        throw std::invalid_argument(
            "Schwarz factorization reuse requires the SparseLU local solver");
    if (config_.factorization_reuse
            != HelmholtzSchwarzFactorizationReuse::None
        && config_.factorization_reuse
            != HelmholtzSchwarzFactorizationReuse::IdenticalMatrix)
        throw std::invalid_argument(
            "unknown Schwarz factorization reuse policy");

    std::unordered_map<std::uint64_t, std::vector<int>> hash_groups;
    solver_groups_.reserve(patch_count);
    for (int target = 0; target < patch_count; ++target) {
        local_matrices[target].makeCompressed();
        int group_index = -1;
        if (config_.factorization_reuse
            == HelmholtzSchwarzFactorizationReuse::IdenticalMatrix) {
            const std::uint64_t hash =
                sparse_matrix_hash(local_matrices[target]);
            const auto found = hash_groups.find(hash);
            if (found != hash_groups.end()) {
                for (const int candidate : found->second) {
                    const int representative =
                        solver_groups_[candidate].subdomains.front();
                    if (sparse_matrices_equal(
                            local_matrices[target],
                            local_matrices[representative])) {
                        group_index = candidate;
                        break;
                    }
                }
            }
            if (group_index < 0) {
                group_index = static_cast<int>(solver_groups_.size());
                solver_groups_.push_back(SolverGroup{});
                hash_groups[hash].push_back(group_index);
            }
        } else {
            group_index = static_cast<int>(solver_groups_.size());
            solver_groups_.push_back(SolverGroup{});
        }
        solver_groups_[group_index].subdomains.push_back(target);
    }

    diagnostics_.local_solver_groups =
        static_cast<int>(solver_groups_.size());
    diagnostics_.reused_factorizations =
        patch_count - diagnostics_.local_solver_groups;
    for (const SolverGroup &group : solver_groups_) {
        diagnostics_.max_reuse_group = std::max(
            diagnostics_.max_reuse_group,
            static_cast<int>(group.subdomains.size()));
    }

    std::vector<char> factorization_failed(solver_groups_.size(), 0);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int group_index = 0;
         group_index < static_cast<int>(solver_groups_.size());
         ++group_index) {
        const int target = solver_groups_[group_index].subdomains.front();
        try {
            solver_groups_[group_index].solver =
                std::make_unique<HelmholtzSchwarzLocalSolver>(
                    local_matrices[target],
                    local_masses[target],
                    local_prolongations[target],
                    model.operators().wavenumber,
                    config_.local_solver);
        } catch (...) {
            factorization_failed[group_index] = 1;
        }
    }
    for (int group_index = 0;
         group_index < static_cast<int>(solver_groups_.size());
         ++group_index) {
        if (factorization_failed[group_index])
            throw std::runtime_error(
                "fine Helmholtz Schwarz local solver setup failed for group "
                + std::to_string(group_index));
    }
    if (use_vcycle) {
        local_solver_diagnostics_.min_vcycle_coarse_dofs =
            std::numeric_limits<int>::max();
        for (const SolverGroup &group : solver_groups_) {
            local_solver_diagnostics_.max_vcycle_levels = std::max(
                local_solver_diagnostics_.max_vcycle_levels,
                group.solver->vcycle_levels());
            local_solver_diagnostics_.min_vcycle_coarse_dofs = std::min(
                local_solver_diagnostics_.min_vcycle_coarse_dofs,
                group.solver->vcycle_coarse_dofs());
            local_solver_diagnostics_.max_vcycle_finest_dofs = std::max(
                local_solver_diagnostics_.max_vcycle_finest_dofs,
                group.solver->vcycle_finest_dofs());
        }
    }

    Eigen::VectorXd partition_sum = Eigen::VectorXd::Zero(dimension_);
    for (const Subdomain &subdomain : subdomains_) {
        for (int local = 0;
             local < static_cast<int>(subdomain.global_dofs.size());
             ++local) {
            partition_sum(subdomain.global_dofs[local]) +=
                subdomain.injection_weights(local);
        }
    }
    diagnostics_.subdomains = static_cast<int>(subdomains_.size());
    diagnostics_.partition_unity_error =
        (partition_sum - Eigen::VectorXd::Ones(dimension_))
            .lpNorm<Eigen::Infinity>();
}

void HelmholtzTwoLevelSchwarzPreconditioner::validate_right_hand_side(
    const ComplexVector &right_hand_side) const {
    if (right_hand_side.size() != dimension_)
        throw std::invalid_argument(
            "fine Helmholtz Schwarz right-hand side size mismatch");
    if (!right_hand_side.allFinite())
        throw std::invalid_argument(
            "fine Helmholtz Schwarz right-hand side is non-finite");
}

ComplexVector HelmholtzTwoLevelSchwarzPreconditioner::apply_coarse(
    const ComplexVector &right_hand_side) const {
    validate_right_hand_side(right_hand_side);
    return model_->solve_load(right_hand_side).fine_values;
}

ComplexVector HelmholtzTwoLevelSchwarzPreconditioner::apply_local(
    const ComplexVector &right_hand_side) const {
    validate_right_hand_side(right_hand_side);

    int worker_count = 1;
#ifdef _OPENMP
    worker_count = omp_get_max_threads();
#endif
    std::vector<ComplexVector> partial_results(
        worker_count, ComplexVector::Zero(dimension_));
    std::vector<char> solve_failed(solver_groups_.size(), 0);
    std::vector<int> solve_iterations(subdomains_.size(), 0);
    std::vector<int> solve_restarts(subdomains_.size(), 0);
    std::vector<double> solve_residuals(subdomains_.size(), 0.0);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int group_index = 0;
         group_index < static_cast<int>(solver_groups_.size());
         ++group_index) {
        const SolverGroup &group = solver_groups_[group_index];
        int worker = 0;
#ifdef _OPENMP
        worker = omp_get_thread_num();
#endif
        try {
            if (config_.local_solver.kind
                    == HelmholtzSchwarzLocalSolverKind::SparseLu
                && group.subdomains.size() > 1) {
                const int local_size = static_cast<int>(
                    subdomains_[group.subdomains.front()].global_dofs.size());
                ComplexMatrix local_rhs(
                    local_size,
                    static_cast<int>(group.subdomains.size()));
                for (int column = 0;
                     column < static_cast<int>(group.subdomains.size());
                     ++column) {
                    const Subdomain &subdomain =
                        subdomains_[group.subdomains[column]];
                    for (int local = 0; local < local_size; ++local) {
                        local_rhs(local, column) =
                            right_hand_side(subdomain.global_dofs[local]);
                    }
                }
                const ComplexMatrix local_solutions =
                    group.solver->solve_block(local_rhs);
                for (int column = 0;
                     column < static_cast<int>(group.subdomains.size());
                     ++column) {
                    const Subdomain &subdomain =
                        subdomains_[group.subdomains[column]];
                    for (int local = 0; local < local_size; ++local) {
                        const int global = subdomain.global_dofs[local];
                        partial_results[worker](global) +=
                            subdomain.injection_weights(local)
                            * local_solutions(local, column);
                    }
                }
            } else {
                if (group.subdomains.size() != 1)
                    throw std::logic_error(
                        "scalar Schwarz solver group must contain one patch");
                const int subdomain_index = group.subdomains.front();
                const Subdomain &subdomain = subdomains_[subdomain_index];
                ComplexVector local_rhs(subdomain.global_dofs.size());
                for (int local = 0; local < local_rhs.size(); ++local) {
                    local_rhs(local) =
                        right_hand_side(subdomain.global_dofs[local]);
                }
                const HelmholtzSchwarzLocalSolveResult solve =
                    group.solver->solve(local_rhs);
                solve_iterations[subdomain_index] = solve.iterations;
                solve_restarts[subdomain_index] = solve.restarts;
                solve_residuals[subdomain_index] = solve.relative_residual;
                for (int local = 0; local < solve.solution.size(); ++local) {
                    const int global = subdomain.global_dofs[local];
                    partial_results[worker](global) +=
                        subdomain.injection_weights(local)
                        * solve.solution(local);
                }
            }
        } catch (...) {
            solve_failed[group_index] = 1;
        }
    }
    if (std::find(solve_failed.begin(), solve_failed.end(), 1)
        != solve_failed.end())
        throw std::runtime_error("fine Helmholtz Schwarz local solve failed");

    local_solver_diagnostics_.solve_calls += subdomains_.size();
    for (int index = 0; index < static_cast<int>(subdomains_.size()); ++index) {
        local_solver_diagnostics_.total_iterations += solve_iterations[index];
        local_solver_diagnostics_.total_restarts += solve_restarts[index];
        local_solver_diagnostics_.max_iterations = std::max(
            local_solver_diagnostics_.max_iterations,
            solve_iterations[index]);
        local_solver_diagnostics_.max_relative_residual = std::max(
            local_solver_diagnostics_.max_relative_residual,
            solve_residuals[index]);
    }

    ComplexVector result = ComplexVector::Zero(dimension_);
    for (const ComplexVector &partial : partial_results) result += partial;
    return result;
}

ComplexVector HelmholtzTwoLevelSchwarzPreconditioner::apply(
    const ComplexVector &right_hand_side,
    HelmholtzTwoLevelSchwarzMode mode) const {
    validate_right_hand_side(right_hand_side);
    const ComplexVector coarse = apply_coarse(right_hand_side);
    if (mode == HelmholtzTwoLevelSchwarzMode::Additive)
        return coarse + apply_local(right_hand_side);
    if (mode != HelmholtzTwoLevelSchwarzMode::Hybrid)
        throw std::invalid_argument("unknown two-level Schwarz mode");

    const ComplexVector right_projected =
        right_hand_side - *fine_operator_ * coarse;
    const ComplexVector local = apply_local(right_projected);
    return coarse + local - apply_coarse(*fine_operator_ * local);
}

} // namespace lod2d::helmholtz

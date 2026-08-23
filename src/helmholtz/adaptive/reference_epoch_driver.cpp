#include "helmholtz/adaptive/reference_epoch_driver.h"

#include <cmath>
#include <optional>
#include <stdexcept>
#include <utility>

namespace lod2d::helmholtz::adaptive {
namespace {

bool terminal(ReferenceEpochDriverState state) {
    return state == ReferenceEpochDriverState::Converged
        || state == ReferenceEpochDriverState::WorkLimitReached
        || state == ReferenceEpochDriverState::Failed;
}

} // namespace

const char *reference_epoch_driver_state_name(ReferenceEpochDriverState state) {
    switch (state) {
    case ReferenceEpochDriverState::EpochInit: return "EpochInit";
    case ReferenceEpochDriverState::CorrectorCheck: return "CorrectorCheck";
    case ReferenceEpochDriverState::SolveEstimate: return "SolveEstimate";
    case ReferenceEpochDriverState::ProposeCoarse: return "ProposeCoarse";
    case ReferenceEpochDriverState::CandidateEnrich: return "CandidateEnrich";
    case ReferenceEpochDriverState::LazyDualDecision: return "LazyDualDecision";
    case ReferenceEpochDriverState::CandidateDualCheck: return "CandidateDualCheck";
    case ReferenceEpochDriverState::CommitCoarse: return "CommitCoarse";
    case ReferenceEpochDriverState::ReferenceRefresh: return "ReferenceRefresh";
    case ReferenceEpochDriverState::Converged: return "Converged";
    case ReferenceEpochDriverState::WorkLimitReached: return "WorkLimitReached";
    case ReferenceEpochDriverState::Failed: return "Failed";
    }
    throw std::invalid_argument("unknown reference-epoch driver state");
}

const char *reference_epoch_driver_action_name(ReferenceEpochDriverAction action) {
    switch (action) {
    case ReferenceEpochDriverAction::BeginEpoch: return "BeginEpoch";
    case ReferenceEpochDriverAction::AcceptCorrector: return "AcceptCorrector";
    case ReferenceEpochDriverAction::IncreaseGlobalEll: return "IncreaseGlobalEll";
    case ReferenceEpochDriverAction::SolveAndEstimate: return "SolveAndEstimate";
    case ReferenceEpochDriverAction::ProposeCoarseRefinement: return "ProposeCoarseRefinement";
    case ReferenceEpochDriverAction::EnrichCandidate: return "EnrichCandidate";
    case ReferenceEpochDriverAction::SkipCandidateDual: return "SkipCandidateDual";
    case ReferenceEpochDriverAction::RequestCandidateDual: return "RequestCandidateDual";
    case ReferenceEpochDriverAction::ComputeCandidateDual: return "ComputeCandidateDual";
    case ReferenceEpochDriverAction::CommitCoarseRefinement: return "CommitCoarseRefinement";
    case ReferenceEpochDriverAction::RefreshReference: return "RefreshReference";
    case ReferenceEpochDriverAction::Complete: return "Complete";
    case ReferenceEpochDriverAction::StopWorkLimit: return "StopWorkLimit";
    case ReferenceEpochDriverAction::Fail: return "Fail";
    }
    throw std::invalid_argument("unknown reference-epoch driver action");
}

ReferenceEpochPracticalDriver::ReferenceEpochPracticalDriver(
    ReferenceEpochDriverBackend &backend,
    ReferenceEpochDriverConfig config)
    : backend_(backend), config_(std::move(config)) {
    if (config_.ell0 < 0 || config_.ell_max < config_.ell0
        || !(config_.tau_loc >= 0.0)
        || !(config_.q_dual > 0.0 && config_.q_dual < 1.0)
        || config_.m_dual == 0 || !(config_.tau_ep > 0.0)
        || config_.reference_refresh_level_gap < 0
        || config_.reference_refresh_target_gap < 0
        || (config_.reference_refresh_target_gap > 0
            && config_.reference_refresh_target_gap
                <= config_.reference_refresh_level_gap)
        || !(config_.tolerance_reference >= 0.0)
        || config_.limits.maximum_H_steps == 0
        || config_.minimum_solved_points_per_new_epoch
            > config_.limits.maximum_H_steps
        || config_.limits.maximum_epochs == 0
        || config_.limits.maximum_dual_checks == 0
        || !(config_.limits.maximum_wall_seconds > 0.0)) {
        throw std::invalid_argument("reference-epoch driver config is invalid");
    }
}

ReferenceEpochDriverResult ReferenceEpochPracticalDriver::run() {
    ReferenceEpochDriverResult result;
    ReferenceEpochDriverState state = ReferenceEpochDriverState::EpochInit;
    std::size_t epoch = 0;
    std::size_t H_steps = 0;
    std::size_t dual_checks = 0;
    int ell = config_.ell0;
    double U_practical = std::numeric_limits<double>::infinity();
    std::vector<int> marked_H;
    std::optional<double> last_dual_U;
    std::size_t last_dual_H_step = 0;
    std::size_t epoch_start_H_step = 0;
    bool structural_trigger = false;
    bool level_gap_trigger = false;
    bool termination_trigger = false;
    bool refresh_commit_pending = false;
    const auto start = std::chrono::steady_clock::now();

    const auto append = [&](ReferenceEpochDriverRecord record) {
        record.sequence = result.journal.size();
        record.epoch = epoch;
        record.H_step = H_steps;
        record.ell = ell;
        const ReferenceEpochResourceSnapshot snapshot = backend_.resources();
        record.coarse_unknowns = snapshot.coarse_unknowns;
        record.reference_unknowns = snapshot.reference_unknowns;
        record.candidate_unknowns = snapshot.candidate_unknowns;
        record.kappa_H_max = snapshot.kappa_H_max;
        record.time_total_cumulative = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        result.journal.push_back(std::move(record));
    };
    const auto stop_limit = [&](const std::string &reason) {
        ReferenceEpochDriverRecord record;
        record.state_before = state;
        state = ReferenceEpochDriverState::WorkLimitReached;
        record.state_after = state;
        record.action = ReferenceEpochDriverAction::StopWorkLimit;
        record.detail = reason;
        append(std::move(record));
        result.stop_reason = reason;
    };
    const auto check_limits = [&]() -> std::optional<std::string> {
        const double wall = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        if (wall >= config_.limits.maximum_wall_seconds)
            return "maximum_wall_seconds reached";
        const ReferenceEpochResourceSnapshot resources = backend_.resources();
        if (resources.reference_unknowns
            > config_.limits.maximum_reference_unknowns)
            return "maximum_reference_unknowns reached";
        if (resources.candidate_unknowns
            > config_.limits.maximum_candidate_unknowns)
            return "maximum_candidate_unknowns reached";
        return std::nullopt;
    };

    try {
        while (!terminal(state)) {
            if (const auto limit = check_limits()) {
                stop_limit(*limit);
                break;
            }
            ReferenceEpochDriverRecord record;
            record.state_before = state;
            switch (state) {
            case ReferenceEpochDriverState::EpochInit:
                if (epoch >= config_.limits.maximum_epochs) {
                    stop_limit("maximum_epochs reached");
                    continue;
                }
                {
                    const auto mesh_start = std::chrono::steady_clock::now();
                    backend_.begin_epoch();
                    record.time_mesh = std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - mesh_start).count();
                }
                epoch_start_H_step = H_steps;
                last_dual_U.reset();
                last_dual_H_step = H_steps;
                state = ReferenceEpochDriverState::CorrectorCheck;
                record.state_after = state;
                record.action = ReferenceEpochDriverAction::BeginEpoch;
                record.minimum_reference_level_gap =
                    backend_.minimum_reference_level_gap();
                record.detail =
                    epoch == 0
                    ? "candidate reset to the fixed reference mesh; initialized ell to ell0"
                    : "candidate reset to the fixed reference mesh; inherited ell across epoch";
                append(std::move(record));
                continue;

            case ReferenceEpochDriverState::CorrectorCheck: {
                const ReferenceEpochCorrectorObservation observation =
                    backend_.corrector_check(ell);
                record.theta_loc = observation.theta_loc;
                record.delta_loc_hat = observation.delta_loc_hat;
                record.rebuilt_correctors = observation.rebuilt_correctors;
                record.skipped_correctors = observation.skipped_correctors;
                record.skipped_corrector_work_units =
                    observation.skipped_corrector_work_units;
                record.hybrid_l_s = observation.hybrid_l_s;
                record.hybrid_minimum_physical_radius =
                    observation.hybrid_minimum_physical_radius;
                record.hybrid_covered_physical_radius =
                    observation.hybrid_covered_physical_radius;
                record.time_corrector = observation.time_corrector;
                record.time_theta = observation.time_theta;
                record.time_gram_prepare_structure =
                    observation.time_gram_prepare_structure;
                record.time_gram_prepare_factorization =
                    observation.time_gram_prepare_factorization;
                record.time_gram_action_rhs = observation.time_gram_action_rhs;
                record.time_gram_action_patch_solve =
                    observation.time_gram_action_patch_solve;
                record.time_gram_action_scatter =
                    observation.time_gram_action_scatter;
                record.gram_action_calls = observation.gram_action_calls;
                record.gram_patch_factorizations =
                    observation.gram_patch_factorizations;
                record.gram_factor_cache_hits =
                    observation.gram_factor_cache_hits;
                record.gram_factor_cache_misses =
                    observation.gram_factor_cache_misses;
                record.gram_parallel_threads = observation.gram_parallel_threads;
                record.localization_iterations =
                    observation.localization_iterations;
                record.localization_relative_residual =
                    observation.localization_relative_residual;
                record.localization_used_warm_start =
                    observation.localization_used_warm_start;
                if (!std::isfinite(observation.delta_loc_hat)
                    || observation.delta_loc_hat < 0.0)
                    throw std::runtime_error("corrector check returned an invalid bound");
                if (observation.delta_loc_hat > config_.tau_loc) {
                    if (ell >= config_.ell_max) {
                        stop_limit("corrector threshold failed at ell_max");
                        continue;
                    }
                    ++ell;
                    record.state_after = state;
                    record.action = ReferenceEpochDriverAction::IncreaseGlobalEll;
                    record.detail = "corrector failure changed only global ell";
                } else {
                    state = ReferenceEpochDriverState::SolveEstimate;
                    record.state_after = state;
                    record.action = ReferenceEpochDriverAction::AcceptCorrector;
                }
                append(std::move(record));
                continue;
            }

            case ReferenceEpochDriverState::SolveEstimate: {
                if (H_steps >= config_.limits.maximum_H_steps) {
                    stop_limit("maximum_H_steps reached");
                    continue;
                }
                const ReferenceEpochSolveObservation observation =
                    backend_.solve_and_estimate();
                if (!(observation.eta_H >= 0.0)
                    || !(observation.U_practical >= 0.0)
                    || !std::isfinite(observation.U_practical))
                    throw std::runtime_error("solve/estimate returned invalid values");
                U_practical = observation.U_practical;
                marked_H = observation.marked_H;
                if (marked_H.empty() && U_practical > config_.tolerance_reference)
                    throw std::runtime_error("positive error bound produced empty H marking");
                state = ReferenceEpochDriverState::ProposeCoarse;
                record.state_after = state;
                record.action = ReferenceEpochDriverAction::SolveAndEstimate;
                record.eta_H = observation.eta_H;
                record.U_practical = U_practical;
                record.marked_H = observation.marked_H.size();
                record.relative_reference_energy =
                    observation.relative_reference_energy;
                record.relative_exact_energy = observation.relative_exact_energy;
                record.relative_exact_L2 = observation.relative_exact_L2;
                record.time_lod_solve = observation.time_lod_solve;
                record.time_reference_riesz =
                    observation.time_reference_riesz;
                append(std::move(record));
                continue;
            }

            case ReferenceEpochDriverState::ProposeCoarse: {
                const auto mesh_start = std::chrono::steady_clock::now();
                backend_.propose_coarse_refinement(marked_H);
                record.time_mesh = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - mesh_start).count();
                state = ReferenceEpochDriverState::CandidateEnrich;
                record.state_after = state;
                record.action = ReferenceEpochDriverAction::ProposeCoarseRefinement;
                append(std::move(record));
                continue;
            }

            case ReferenceEpochDriverState::CandidateEnrich: {
                const ReferenceEpochCandidateObservation observation =
                    backend_.enrich_candidate();
                if (!(observation.eta_eq_c >= 0.0)
                    || !std::isfinite(observation.eta_eq_c))
                    throw std::runtime_error("candidate enrichment returned invalid eta_eq_c");
                state = ReferenceEpochDriverState::LazyDualDecision;
                record.state_after = state;
                record.action = ReferenceEpochDriverAction::EnrichCandidate;
                record.eta_eq_c = observation.eta_eq_c;
                record.marked_c = observation.marked_c.size();
                record.time_candidate_flux = observation.time_candidate_flux;
                record.time_candidate_close = observation.time_candidate_close;
                record.time_candidate_operator_assembly =
                    observation.time_candidate_operator_assembly;
                record.time_candidate_prolongation =
                    observation.time_candidate_prolongation;
                record.time_candidate_flux_reconstruction =
                    observation.time_candidate_flux_reconstruction;
                record.time_candidate_flux_prepare =
                    observation.time_candidate_flux_prepare;
                record.time_candidate_flux_patch_solve =
                    observation.time_candidate_flux_patch_solve;
                record.time_candidate_flux_merge =
                    observation.time_candidate_flux_merge;
                record.time_candidate_flux_audit =
                    observation.time_candidate_flux_audit;
                record.candidate_flux_parallel_threads =
                    observation.candidate_flux_parallel_threads;
                record.time_candidate_enrich =
                    observation.time_candidate_enrich;
                record.time_candidate_nvb_refine =
                    observation.time_candidate_nvb_refine;
                record.time_candidate_embedding_composition =
                    observation.time_candidate_embedding_composition;
                record.time_candidate_parent_map_update =
                    observation.time_candidate_parent_map_update;
                record.time_candidate_quasi_interpolation =
                    observation.time_candidate_quasi_interpolation;
                record.time_candidate_embedding_validation =
                    observation.time_candidate_embedding_validation;
                append(std::move(record));
                continue;
            }

            case ReferenceEpochDriverState::LazyDualDecision: {
                structural_trigger = !backend_.proposal_contained_in_reference();
                const int local_level_gap =
                    backend_.minimum_reference_level_gap();
                level_gap_trigger =
                    config_.reference_refresh_level_gap > 0
                    && H_steps - epoch_start_H_step
                        >= config_.minimum_H_steps_per_epoch
                    && local_level_gap
                        <= config_.reference_refresh_level_gap;
                termination_trigger = U_practical <= config_.tolerance_reference;
                const bool decrease_trigger = last_dual_U.has_value()
                    && U_practical <= config_.q_dual * *last_dual_U;
                const bool interval_trigger =
                    H_steps - last_dual_H_step >= config_.m_dual;
                const bool forced_refresh = structural_trigger
                    || level_gap_trigger;
                const std::size_t remaining_new_epoch_solves =
                    H_steps + 1 >= config_.limits.maximum_H_steps
                    ? 0
                    : config_.limits.maximum_H_steps - (H_steps + 1);
                const bool refresh_budget_available =
                    config_.minimum_solved_points_per_new_epoch == 0
                    || remaining_new_epoch_solves
                        >= config_.minimum_solved_points_per_new_epoch;
                const bool numerical_trigger = termination_trigger
                    || decrease_trigger || interval_trigger;
                record.structural_dual_trigger = structural_trigger;
                record.level_gap_dual_trigger = level_gap_trigger;
                record.minimum_reference_level_gap = local_level_gap;
                record.tolerance_dual_trigger = termination_trigger;
                record.decrease_dual_trigger = decrease_trigger;
                record.interval_dual_trigger = interval_trigger;
                state = forced_refresh
                    ? ReferenceEpochDriverState::ReferenceRefresh
                    : numerical_trigger
                    ? ReferenceEpochDriverState::CandidateDualCheck
                    : ReferenceEpochDriverState::CommitCoarse;
                record.state_after = state;
                record.action = numerical_trigger && !forced_refresh
                    ? ReferenceEpochDriverAction::RequestCandidateDual
                    : ReferenceEpochDriverAction::SkipCandidateDual;
                if (forced_refresh) {
                    record.detail = structural_trigger
                        ? "structural refresh decided without redundant candidate dual"
                        : "level-gap refresh decided without redundant candidate dual";
                }
                append(std::move(record));
                if (forced_refresh && !refresh_budget_available) {
                    stop_limit(
                        "insufficient remaining H-step budget for a new reference epoch");
                    continue;
                }
                continue;
            }

            case ReferenceEpochDriverState::CandidateDualCheck: {
                if (dual_checks >= config_.limits.maximum_dual_checks) {
                    stop_limit("maximum_dual_checks reached");
                    continue;
                }
                const ReferenceEpochDualObservation observation =
                    backend_.candidate_dual_check(U_practical);
                if (!(observation.eta_dual_c >= 0.0)
                    || !(observation.L_gap_c >= 0.0)
                    || !std::isfinite(observation.L_gap_c))
                    throw std::runtime_error("candidate dual check returned invalid values");
                ++dual_checks;
                last_dual_U = U_practical;
                last_dual_H_step = H_steps;
                const bool numerical_refresh = observation.L_gap_c
                    >= config_.tau_ep * U_practical;
                const std::size_t remaining_new_epoch_solves =
                    H_steps + 1 >= config_.limits.maximum_H_steps
                    ? 0
                    : config_.limits.maximum_H_steps - (H_steps + 1);
                const bool refresh_budget_available =
                    config_.minimum_solved_points_per_new_epoch == 0
                    || remaining_new_epoch_solves
                        >= config_.minimum_solved_points_per_new_epoch;
                if (numerical_refresh && refresh_budget_available) {
                    state = ReferenceEpochDriverState::ReferenceRefresh;
                } else if (numerical_refresh) {
                    state = ReferenceEpochDriverState::WorkLimitReached;
                    result.stop_reason =
                        "insufficient remaining H-step budget for a new reference epoch";
                } else if (termination_trigger) {
                    state = ReferenceEpochDriverState::Converged;
                } else {
                    state = ReferenceEpochDriverState::CommitCoarse;
                }
                record.state_after = state;
                record.action = ReferenceEpochDriverAction::ComputeCandidateDual;
                record.eta_dual_c = observation.eta_dual_c;
                record.L_gap_c = observation.L_gap_c;
                record.time_candidate_dual = observation.time_candidate_dual;
                record.time_candidate_dual_operator_assembly =
                    observation.time_candidate_dual_operator_assembly;
                record.time_candidate_dual_load_assembly =
                    observation.time_candidate_dual_load_assembly;
                record.time_candidate_dual_prolongation =
                    observation.time_candidate_dual_prolongation;
                record.time_candidate_dual_solve =
                    observation.time_candidate_dual_solve;
                record.time_candidate_dual_prepare =
                    observation.time_candidate_dual_prepare;
                record.time_candidate_dual_patch_solve =
                    observation.time_candidate_dual_patch_solve;
                record.time_candidate_dual_reduction =
                    observation.time_candidate_dual_reduction;
                record.candidate_dual_patch_factorizations =
                    observation.candidate_dual_patch_factorizations;
                record.candidate_dual_parallel_threads =
                    observation.candidate_dual_parallel_threads;
                record.structural_dual_trigger = structural_trigger;
                record.level_gap_dual_trigger = level_gap_trigger;
                record.minimum_reference_level_gap =
                    backend_.minimum_reference_level_gap();
                record.tolerance_dual_trigger = termination_trigger;
                append(std::move(record));
                if (state == ReferenceEpochDriverState::WorkLimitReached)
                    continue;
                if (state == ReferenceEpochDriverState::Converged) {
                    ReferenceEpochDriverRecord complete;
                    complete.state_before = ReferenceEpochDriverState::Converged;
                    complete.state_after = ReferenceEpochDriverState::Converged;
                    complete.action = ReferenceEpochDriverAction::Complete;
                    complete.detail = "reference tolerance reached after forced dual check";
                    append(std::move(complete));
                }
                continue;
            }

            case ReferenceEpochDriverState::ReferenceRefresh:
                if (epoch + 1 >= config_.limits.maximum_epochs) {
                    stop_limit("maximum_epochs reached before reference refresh");
                    continue;
                }
                {
                    const auto mesh_start = std::chrono::steady_clock::now();
                    backend_.refresh_reference(
                        config_.reference_refresh_target_gap);
                    record.time_mesh = std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - mesh_start).count();
                }
                record.minimum_reference_level_gap =
                    backend_.minimum_reference_level_gap();
                refresh_commit_pending = true;
                state = ReferenceEpochDriverState::CommitCoarse;
                record.state_after = state;
                record.action = ReferenceEpochDriverAction::RefreshReference;
                record.detail = structural_trigger
                    ? "structural hierarchy trigger took precedence"
                    : level_gap_trigger
                    ? "local reference/coarse level-gap guard requested reference refresh"
                    : "candidate dual gap requested reference refresh";
                append(std::move(record));
                continue;

            case ReferenceEpochDriverState::CommitCoarse: {
                const auto mesh_start = std::chrono::steady_clock::now();
                backend_.commit_coarse_refinement();
                record.time_mesh = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - mesh_start).count();
                ++H_steps;
                if (refresh_commit_pending) {
                    refresh_commit_pending = false;
                    ++epoch;
                    state = ReferenceEpochDriverState::EpochInit;
                    record.detail = "committed prospective coarse mesh in refreshed reference";
                } else {
                    state = ReferenceEpochDriverState::CorrectorCheck;
                    record.detail = "committed prospective coarse mesh in current epoch";
                }
                record.state_after = state;
                record.action = ReferenceEpochDriverAction::CommitCoarseRefinement;
                append(std::move(record));
                continue;
            }

            case ReferenceEpochDriverState::Converged:
            case ReferenceEpochDriverState::WorkLimitReached:
            case ReferenceEpochDriverState::Failed:
                break;
            }
        }
    } catch (const std::exception &error) {
        ReferenceEpochDriverRecord record;
        record.state_before = state;
        state = ReferenceEpochDriverState::Failed;
        record.state_after = state;
        record.action = ReferenceEpochDriverAction::Fail;
        record.detail = error.what();
        append(std::move(record));
        result.stop_reason = error.what();
    }
    result.state = state;
    result.epochs = epoch + (result.journal.empty() ? 0 : 1);
    result.H_steps = H_steps;
    result.dual_checks = dual_checks;
    result.ell = ell;
    if (state == ReferenceEpochDriverState::Converged)
        result.stop_reason = "reference tolerance reached";
    return result;
}

} // namespace lod2d::helmholtz::adaptive

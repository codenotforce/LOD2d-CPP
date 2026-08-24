#include "helmholtz/adaptive/reference_epoch_driver.h"

#include <algorithm>
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
    std::size_t solved_points = 0;
    std::size_t dual_checks = 0;
    int ell = config_.ell0;
    double U_practical = std::numeric_limits<double>::infinity();
    std::vector<int> marked_H;
    std::optional<double> last_dual_U;
    std::size_t last_dual_H_step = 0;
    bool structural_trigger = false;
    bool level_gap_trigger = false;
    bool termination_trigger = false;
    bool refresh_commit_pending = false;
    double validation_time_cumulative = 0.0;
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
        record.time_validation_cumulative = validation_time_cumulative;
        record.time_artifact_capture_cumulative =
            snapshot.artifact_capture_seconds;
        record.time_method_cumulative = std::max(
            0.0, record.time_total_cumulative - validation_time_cumulative
                - snapshot.artifact_capture_seconds);
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
                record.active_correctors = observation.active_correctors;
                record.rebuilt_correctors = observation.rebuilt_correctors;
                record.reused_correctors = observation.reused_correctors;
                record.corrector_cache_oversized_misses =
                    observation.corrector_cache_oversized_misses;
                record.corrector_cache_budget_rejections =
                    observation.corrector_cache_budget_rejections;
                record.corrector_cache_entries =
                    observation.corrector_cache_entries;
                record.corrector_cache_current_bytes =
                    observation.corrector_cache_current_bytes;
                record.corrector_cache_peak_bytes =
                    observation.corrector_cache_peak_bytes;
                record.skipped_correctors = observation.skipped_correctors;
                record.skipped_corrector_work_units =
                    observation.skipped_corrector_work_units;
                record.corrector_parallel_threads =
                    observation.corrector_parallel_threads;
                record.corrector_patch_assembly_work_seconds =
                    observation.corrector_patch_assembly_work_seconds;
                record.corrector_patch_solve_work_seconds =
                    observation.corrector_patch_solve_work_seconds;
                record.corrector_patch_pack_work_seconds =
                    observation.corrector_patch_pack_work_seconds;
                record.corrector_maximum_patch_dofs =
                    observation.corrector_maximum_patch_dofs;
                record.corrector_maximum_patch_constraints =
                    observation.corrector_maximum_patch_constraints;
                record.corrector_maximum_patch_rhs =
                    observation.corrector_maximum_patch_rhs;
                record.hybrid_l_s = observation.hybrid_l_s;
                record.hybrid_minimum_physical_radius =
                    observation.hybrid_minimum_physical_radius;
                record.hybrid_covered_physical_radius =
                    observation.hybrid_covered_physical_radius;
                record.hybrid_omega_s_elements =
                    observation.hybrid_omega_s_elements;
                record.hybrid_omega_f_elements =
                    observation.hybrid_omega_f_elements;
                record.time_reference_stability =
                    observation.time_reference_stability;
                record.time_corrector_check_total =
                    observation.time_corrector_check_total;
                record.time_lod_build_total = observation.time_lod_build_total;
                record.time_lod_mesh_and_interpolation =
                    observation.time_lod_mesh_and_interpolation;
                record.time_lod_operators = observation.time_lod_operators;
                record.time_corrector = observation.time_corrector;
                record.time_lod_corrected_basis =
                    observation.time_lod_corrected_basis;
                record.time_lod_coarse_operator =
                    observation.time_lod_coarse_operator;
                record.time_lod_coarse_factorization =
                    observation.time_lod_coarse_factorization;
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
                record.gram_structure_parallel_threads =
                    observation.gram_structure_parallel_threads;
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
                        const std::string reason =
                            "corrector threshold failed at ell_max";
                        record.state_after =
                            ReferenceEpochDriverState::WorkLimitReached;
                        record.action =
                            ReferenceEpochDriverAction::StopWorkLimit;
                        record.detail = reason;
                        state = ReferenceEpochDriverState::WorkLimitReached;
                        append(std::move(record));
                        result.stop_reason = reason;
                        continue;
                    }
                    record.state_after = state;
                    record.action = ReferenceEpochDriverAction::IncreaseGlobalEll;
                    record.detail = "corrector failure changed global ell from "
                        + std::to_string(ell) + " to "
                        + std::to_string(ell + 1);
                    // The observation belongs to the ell that was actually
                    // checked.  Advance only after journaling it so Theta,
                    // Gram timings, and cache counters are not mislabeled.
                    append(std::move(record));
                    ++ell;
                    continue;
                } else {
                    state = ReferenceEpochDriverState::SolveEstimate;
                    record.state_after = state;
                    record.action = ReferenceEpochDriverAction::AcceptCorrector;
                }
                append(std::move(record));
                continue;
            }

            case ReferenceEpochDriverState::SolveEstimate: {
                if (solved_points >= config_.limits.maximum_H_steps) {
                    stop_limit("maximum_H_steps reached");
                    continue;
                }
                const ReferenceEpochSolveObservation observation =
                    backend_.solve_and_estimate();
                if (!(observation.eta_H >= 0.0)
                    || !(observation.U_practical >= 0.0)
                    || !std::isfinite(observation.U_practical))
                    throw std::runtime_error("solve/estimate returned invalid values");
                ++solved_points;
                U_practical = observation.U_practical;
                marked_H = observation.marked_H;
                // Algorithm 1 initializes the lazy-dual reference values at
                // the first solved point of every epoch.  Without this, the
                // decrease trigger is disabled until an interval-triggered
                // dual solve happens, potentially delaying a numerical epoch
                // switch by an entire m_dual block.
                if (!last_dual_U) {
                    last_dual_U = U_practical;
                    last_dual_H_step = H_steps;
                }
                // Algorithm 2 permits an empty regular-region marking when
                // its regional indicator mass vanishes.  The backend opens
                // an identity proposal so candidate enrichment and the epoch
                // decision are still executed transactionally.
                state = ReferenceEpochDriverState::ProposeCoarse;
                record.state_after = state;
                record.action = ReferenceEpochDriverAction::SolveAndEstimate;
                record.eta_H = observation.eta_H;
                record.U_practical = U_practical;
                record.marked_H = observation.marked_H.size();
                record.hybrid_regular_indicator_mass =
                    observation.hybrid_regular_indicator_mass;
                record.hybrid_admissible_indicator_mass =
                    observation.hybrid_admissible_indicator_mass;
                record.hybrid_marked_H_indicator_mass =
                    observation.hybrid_marked_H_indicator_mass;
                record.hybrid_coarse_conformity_collar =
                    observation.hybrid_coarse_conformity_collar;
                record.hybrid_coarse_marking_closure_safe =
                    observation.hybrid_coarse_marking_closure_safe;
                record.hybrid_full_regular_doerfler =
                    observation.hybrid_full_regular_doerfler;
                record.hybrid_coarse_preview_attempts =
                    observation.hybrid_coarse_preview_attempts;
                record.hybrid_coarse_preview_cached =
                    observation.hybrid_coarse_preview_cached;
                record.relative_reference_energy =
                    observation.relative_reference_energy;
                record.relative_exact_energy = observation.relative_exact_energy;
                record.relative_exact_L2 = observation.relative_exact_L2;
                record.relative_exact_reference_energy =
                    observation.relative_exact_reference_energy;
                record.time_lod_solve = observation.time_lod_solve;
                record.time_reference_riesz =
                    observation.time_reference_riesz;
                record.time_hybrid_coarse_marking =
                    observation.time_hybrid_coarse_marking;
                record.time_hybrid_coarse_preview =
                    observation.time_hybrid_coarse_preview;
                record.time_hybrid_coarse_preview_nvb =
                    observation.time_hybrid_coarse_preview_nvb;
                record.time_hybrid_coarse_preview_reference_embedding =
                    observation.time_hybrid_coarse_preview_reference_embedding;
                record.time_reference_validation =
                    observation.time_reference_validation;
                record.time_exact_validation =
                    observation.time_exact_validation;
                validation_time_cumulative +=
                    observation.time_reference_validation
                    + observation.time_exact_validation;
                append(std::move(record));
                // maximum_H_steps is the number of solved coarse states,
                // not the number of committed refinements.  The final
                // solved point is already a valid trajectory datum; do not
                // spend candidate, proposal, commit, or corrector work after
                // recording it.
                if (solved_points >= config_.limits.maximum_H_steps) {
                    stop_limit("maximum_H_steps reached");
                }
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
                record.eta_eq_c_f = observation.eta_eq_c_f;
                record.eta_eq_c_r = observation.eta_eq_c_r;
                record.indicator_mass_c_f = observation.indicator_mass_c_f;
                record.indicator_mass_c_r = observation.indicator_mass_c_r;
                record.marked_mass_c_f = observation.marked_mass_c_f;
                record.marked_mass_c_r = observation.marked_mass_c_r;
                record.marked_c = observation.marked_c.size();
                record.marked_c_f = observation.marked_c_f;
                record.marked_c_r = observation.marked_c_r;
                record.candidate_cells_f = observation.candidate_cells_f;
                record.candidate_cells_r = observation.candidate_cells_r;
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
                    solved_points >= config_.limits.maximum_H_steps
                    ? 0
                    : config_.limits.maximum_H_steps - solved_points;
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
                    solved_points >= config_.limits.maximum_H_steps
                    ? 0
                    : config_.limits.maximum_H_steps - solved_points;
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
                record.time_candidate_wellposedness =
                    observation.time_candidate_wellposedness;
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
                if (state == ReferenceEpochDriverState::WorkLimitReached) {
                    ReferenceEpochDriverRecord stop;
                    stop.state_before =
                        ReferenceEpochDriverState::WorkLimitReached;
                    stop.state_after =
                        ReferenceEpochDriverState::WorkLimitReached;
                    stop.action = ReferenceEpochDriverAction::StopWorkLimit;
                    stop.detail = result.stop_reason;
                    append(std::move(stop));
                    continue;
                }
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
                const bool committed_refresh = refresh_commit_pending;
                if (committed_refresh) {
                    refresh_commit_pending = false;
                    state = ReferenceEpochDriverState::EpochInit;
                    record.detail = "committed prospective coarse mesh in refreshed reference";
                } else {
                    state = ReferenceEpochDriverState::CorrectorCheck;
                    record.detail = "committed prospective coarse mesh in current epoch";
                }
                record.state_after = state;
                record.action = ReferenceEpochDriverAction::CommitCoarseRefinement;
                // The commit closes the old epoch.  Journal it before
                // advancing the epoch counter so its mesh snapshot and CSV
                // row share the same epoch; the following BeginEpoch row is
                // the first row owned by the new epoch.
                append(std::move(record));
                if (committed_refresh) ++epoch;
                continue;
            }

            case ReferenceEpochDriverState::Converged:
            case ReferenceEpochDriverState::WorkLimitReached:
            case ReferenceEpochDriverState::Failed:
                break;
            }
        }
    } catch (const ReferenceEpochWorkLimitExceeded &error) {
        stop_limit(error.what());
    } catch (const ReferenceEpochReserveUnavailable &error) {
        stop_limit(error.what());
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
    result.H_steps = solved_points;
    result.dual_checks = dual_checks;
    result.ell = ell;
    if (state == ReferenceEpochDriverState::Converged)
        result.stop_reason = "reference tolerance reached";
    return result;
}

} // namespace lod2d::helmholtz::adaptive

#pragma once

#include <chrono>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace lod2d::helmholtz::adaptive {

// Backends use this exception only for resource guards that fire inside one
// otherwise atomic driver action (for example, a multi-step matching or
// reserve closure).  The driver converts it to WorkLimitReached instead of
// treating expected resource exhaustion as a numerical failure.
class ReferenceEpochWorkLimitExceeded final : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// A requested post-refresh reserve may be geometrically unavailable while
// preserving an exact hybrid matching region.  This is an expected bounded
// stop, not a numerical failure or an invitation to chase the interface.
class ReferenceEpochReserveUnavailable final : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

enum class ReferenceEpochDriverState {
    EpochInit,
    CorrectorCheck,
    SolveEstimate,
    ProposeCoarse,
    CandidateEnrich,
    LazyDualDecision,
    CandidateDualCheck,
    CommitCoarse,
    ReferenceRefresh,
    Converged,
    WorkLimitReached,
    Failed,
};

enum class ReferenceEpochDriverAction {
    BeginEpoch,
    AcceptCorrector,
    IncreaseGlobalEll,
    SolveAndEstimate,
    ProposeCoarseRefinement,
    EnrichCandidate,
    SkipCandidateDual,
    RequestCandidateDual,
    ComputeCandidateDual,
    CommitCoarseRefinement,
    RefreshReference,
    Complete,
    StopWorkLimit,
    Fail,
};

struct ReferenceEpochCorrectorObservation {
    double theta_loc = 0.0;
    double delta_loc_hat = 0.0;
    std::size_t active_correctors = 0;
    std::size_t rebuilt_correctors = 0;
    std::size_t reused_correctors = 0;
    std::size_t corrector_cache_oversized_misses = 0;
    std::size_t corrector_cache_budget_rejections = 0;
    std::size_t corrector_cache_entries = 0;
    std::size_t corrector_cache_current_bytes = 0;
    std::size_t corrector_cache_peak_bytes = 0;
    std::size_t skipped_correctors = 0;
    std::size_t skipped_corrector_work_units = 0;
    int corrector_parallel_threads = 1;
    double corrector_patch_assembly_work_seconds = 0.0;
    double corrector_patch_solve_work_seconds = 0.0;
    double corrector_patch_pack_work_seconds = 0.0;
    int corrector_maximum_patch_dofs = 0;
    int corrector_maximum_patch_constraints = 0;
    int corrector_maximum_patch_rhs = 0;
    int hybrid_l_s = -1;
    double hybrid_minimum_physical_radius =
        std::numeric_limits<double>::quiet_NaN();
    double hybrid_covered_physical_radius =
        std::numeric_limits<double>::quiet_NaN();
    std::size_t hybrid_omega_s_elements = 0;
    std::size_t hybrid_omega_f_elements = 0;
    // One-time solve/factorization of the fixed reference problem at the
    // first corrector check of an epoch.  This is the practical
    // well-posedness check used by the implementation study.
    double time_reference_stability = 0.0;
    double time_corrector_check_total = 0.0;
    double time_lod_build_total = 0.0;
    double time_lod_mesh_and_interpolation = 0.0;
    double time_lod_operators = 0.0;
    double time_corrector = 0.0;
    double time_lod_corrected_basis = 0.0;
    double time_lod_coarse_operator = 0.0;
    double time_lod_coarse_factorization = 0.0;
    double time_theta = 0.0;
    double time_gram_prepare_structure = 0.0;
    double time_gram_prepare_factorization = 0.0;
    double time_gram_action_rhs = 0.0;
    double time_gram_action_patch_solve = 0.0;
    double time_gram_action_scatter = 0.0;
    std::size_t gram_action_calls = 0;
    std::size_t gram_patch_factorizations = 0;
    std::size_t gram_factor_cache_hits = 0;
    std::size_t gram_factor_cache_misses = 0;
    int gram_structure_parallel_threads = 1;
    int gram_parallel_threads = 1;
    int localization_iterations = 0;
    double localization_relative_residual = 0.0;
    bool localization_used_warm_start = false;
};

struct ReferenceEpochSolveObservation {
    double eta_H = 0.0;
    double U_practical = 0.0;
    std::vector<int> marked_H;
    double hybrid_regular_indicator_mass =
        std::numeric_limits<double>::quiet_NaN();
    double hybrid_admissible_indicator_mass =
        std::numeric_limits<double>::quiet_NaN();
    double hybrid_marked_H_indicator_mass =
        std::numeric_limits<double>::quiet_NaN();
    int hybrid_coarse_conformity_collar = -1;
    bool hybrid_coarse_marking_closure_safe = false;
    bool hybrid_full_regular_doerfler = false;
    std::size_t hybrid_coarse_preview_attempts = 0;
    bool hybrid_coarse_preview_cached = false;
    double relative_reference_energy =
        std::numeric_limits<double>::quiet_NaN();
    double relative_exact_energy =
        std::numeric_limits<double>::quiet_NaN();
    double relative_exact_L2 =
        std::numeric_limits<double>::quiet_NaN();
    double relative_exact_reference_energy =
        std::numeric_limits<double>::quiet_NaN();
    double time_lod_solve = 0.0;
    double time_reference_riesz = 0.0;
    double time_hybrid_coarse_marking = 0.0;
    double time_hybrid_coarse_preview = 0.0;
    double time_hybrid_coarse_preview_nvb = 0.0;
    double time_hybrid_coarse_preview_reference_embedding = 0.0;
    // Experimental diagnostics, excluded from method time.  They compare the
    // computed state with the fixed reference/exact solutions and do not
    // participate in any adaptive decision.
    double time_reference_validation = 0.0;
    double time_exact_validation = 0.0;
};

struct ReferenceEpochCandidateObservation {
    bool accuracy_sweep_performed = true;
    double eta_eq_c = 0.0;
    double eta_eq_c_f = std::numeric_limits<double>::quiet_NaN();
    double eta_eq_c_r = std::numeric_limits<double>::quiet_NaN();
    double indicator_mass_c_f = std::numeric_limits<double>::quiet_NaN();
    double indicator_mass_c_r = std::numeric_limits<double>::quiet_NaN();
    double marked_mass_c_f = std::numeric_limits<double>::quiet_NaN();
    double marked_mass_c_r = std::numeric_limits<double>::quiet_NaN();
    std::size_t marked_c_f = 0;
    std::size_t marked_c_r = 0;
    std::size_t candidate_cells_f = 0;
    std::size_t candidate_cells_r = 0;
    std::vector<int> marked_c;
    std::size_t candidate_elements_before = 0;
    std::size_t containment_requested_marks = 0;
    std::size_t containment_added_elements = 0;
    std::size_t containment_closure_added_elements = 0;
    std::size_t accuracy_requested_marks = 0;
    std::size_t accuracy_added_elements = 0;
    std::size_t accuracy_closure_added_elements = 0;
    std::size_t candidate_elements_after = 0;
    double time_candidate_marking = 0.0;
    std::size_t candidate_marking_pool = 0;
    std::size_t candidate_estimated_selected_closure_cost = 0;
    double time_candidate_flux = 0.0;
    double time_candidate_close = 0.0;
    double time_candidate_operator_assembly = 0.0;
    double time_candidate_prolongation = 0.0;
    double time_candidate_flux_reconstruction = 0.0;
    double time_candidate_flux_prepare = 0.0;
    double time_candidate_flux_patch_solve = 0.0;
    double time_candidate_flux_merge = 0.0;
    double time_candidate_flux_audit = 0.0;
    int candidate_flux_parallel_threads = 1;
    double time_candidate_enrich = 0.0;
    double time_candidate_nvb_refine = 0.0;
    double time_candidate_embedding_composition = 0.0;
    double time_candidate_parent_map_update = 0.0;
    double time_candidate_quasi_interpolation = 0.0;
    double time_candidate_embedding_validation = 0.0;
};

struct ReferenceEpochDualObservation {
    double eta_dual_c = 0.0;
    double L_gap_c = 0.0;
    double time_candidate_dual = 0.0;
    double time_candidate_wellposedness = 0.0;
    double time_candidate_dual_operator_assembly = 0.0;
    double time_candidate_dual_load_assembly = 0.0;
    double time_candidate_dual_prolongation = 0.0;
    double time_candidate_dual_solve = 0.0;
    double time_candidate_dual_prepare = 0.0;
    double time_candidate_dual_patch_solve = 0.0;
    double time_candidate_dual_reduction = 0.0;
    std::size_t candidate_dual_patch_factorizations = 0;
    int candidate_dual_parallel_threads = 1;
};

struct ReferenceEpochResourceSnapshot {
    std::size_t coarse_unknowns = 0;
    std::size_t reference_unknowns = 0;
    std::size_t candidate_unknowns = 0;
    double kappa_H_max = std::numeric_limits<double>::quiet_NaN();
    double artifact_capture_seconds = 0.0;
};

// The production contract intentionally has no candidate Helmholtz solve.
// Candidate enrichment and candidate dual-Riesz consume the current U kept by
// the backend, but neither may construct u_c.
class ReferenceEpochDriverBackend {
public:
    virtual ~ReferenceEpochDriverBackend() = default;
    virtual void begin_epoch() = 0;
    virtual ReferenceEpochCorrectorObservation corrector_check(int ell) = 0;
    virtual ReferenceEpochSolveObservation solve_and_estimate() = 0;
    virtual void propose_coarse_refinement(
        const std::vector<int> &marked_H) = 0;
    // A false accuracy flag still requires hierarchy containment closure, but
    // skips the expensive RT2 estimator/reconstruction and its persistent
    // accuracy refinement.  The driver forces a sweep before every numerical
    // or structural refresh decision.
    virtual ReferenceEpochCandidateObservation enrich_candidate(
        bool perform_accuracy_sweep) = 0;
    virtual bool proposal_contained_in_reference() const = 0;
    virtual int minimum_reference_level_gap() const = 0;
    virtual ReferenceEpochDualObservation candidate_dual_check(
        double reference_upper_bound) = 0;
    virtual void commit_coarse_refinement() = 0;
    virtual void refresh_reference(int minimum_post_refresh_level_gap) = 0;
    virtual ReferenceEpochResourceSnapshot resources() const = 0;
};

struct ReferenceEpochDriverLimits {
    std::size_t maximum_H_steps = 100;
    std::size_t maximum_epochs = 10;
    std::size_t maximum_reference_unknowns = 1000000;
    std::size_t maximum_candidate_unknowns = 2000000;
    std::size_t maximum_dual_checks = 100;
    double maximum_wall_seconds = 3600.0;
};

struct ReferenceEpochDriverConfig {
    // Algorithm 2 uses a moving reference: after every nonterminal
    // candidate enrichment the frozen candidate is promoted immediately.
    // This bypasses the lazy dual/level-reserve state machine used by the
    // standard reference-epoch algorithm.
    bool moving_reference = false;
    int ell0 = 2;
    int ell_max = 6;
    // Threshold for the computable localization bound delta_loc_hat.  The
    // raw Theta_loc remains separately logged in each observation.
    double tau_loc = 0.5;
    double q_dual = 0.5;
    std::size_t m_dual = 3;
    double tau_ep = 0.5;
    // A positive value forces an epoch refresh once the smallest prospective
    // reference/coarse NVB-level gap reaches this threshold.  Since this is a
    // structural decision, the candidate dual solve is deliberately skipped.
    // Zero disables the guard.
    int reference_refresh_level_gap = 0;
    // Before promoting candidate to reference, locally deepen it until the
    // prospective coarse/candidate gap reaches this value.  This prevents a
    // new epoch from starting with an already exhausted reference.
    int reference_refresh_target_gap = 0;
    // Perform the candidate accuracy/RT2 sweep only every N committed H
    // refinements.  Containment closure still runs every step.  One preserves
    // the historical production trajectory.
    std::size_t candidate_update_stride = 1;
    // Force an accuracy sweep once the prospective reference reserve reaches
    // this value.  Zero disables the additional force condition.
    int candidate_force_level_gap = 0;
    // Deprecated experimental provenance field.  The refinement-reserve
    // trigger is mandatory in Algorithms 1--2 and is therefore never delayed
    // by this value.  Production configurations must set it to zero.
    std::size_t minimum_H_steps_per_epoch = 0;
    // Do not open a refreshed epoch unless the remaining fixed-work budget
    // can produce at least this many SolveAndEstimate observations.
    std::size_t minimum_solved_points_per_new_epoch = 0;
    double tolerance_reference = 1e-2;
    ReferenceEpochDriverLimits limits;
};

struct ReferenceEpochDriverRecord {
    std::size_t sequence = 0;
    std::size_t epoch = 0;
    std::size_t H_step = 0;
    int ell = 0;
    ReferenceEpochDriverState state_before =
        ReferenceEpochDriverState::EpochInit;
    ReferenceEpochDriverState state_after =
        ReferenceEpochDriverState::EpochInit;
    ReferenceEpochDriverAction action =
        ReferenceEpochDriverAction::BeginEpoch;
    double theta_loc = std::numeric_limits<double>::quiet_NaN();
    double delta_loc_hat = std::numeric_limits<double>::quiet_NaN();
    double eta_H = std::numeric_limits<double>::quiet_NaN();
    double U_practical = std::numeric_limits<double>::quiet_NaN();
    double eta_eq_c = std::numeric_limits<double>::quiet_NaN();
    double eta_eq_c_f = std::numeric_limits<double>::quiet_NaN();
    double eta_eq_c_r = std::numeric_limits<double>::quiet_NaN();
    double indicator_mass_c_f = std::numeric_limits<double>::quiet_NaN();
    double indicator_mass_c_r = std::numeric_limits<double>::quiet_NaN();
    double marked_mass_c_f = std::numeric_limits<double>::quiet_NaN();
    double marked_mass_c_r = std::numeric_limits<double>::quiet_NaN();
    double eta_dual_c = std::numeric_limits<double>::quiet_NaN();
    double L_gap_c = std::numeric_limits<double>::quiet_NaN();
    double relative_reference_energy = std::numeric_limits<double>::quiet_NaN();
    double relative_exact_energy = std::numeric_limits<double>::quiet_NaN();
    double relative_exact_L2 = std::numeric_limits<double>::quiet_NaN();
    double relative_exact_reference_energy =
        std::numeric_limits<double>::quiet_NaN();
    std::size_t coarse_unknowns = 0;
    std::size_t reference_unknowns = 0;
    std::size_t candidate_unknowns = 0;
    double kappa_H_max = std::numeric_limits<double>::quiet_NaN();
    std::size_t marked_H = 0;
    double hybrid_regular_indicator_mass =
        std::numeric_limits<double>::quiet_NaN();
    double hybrid_admissible_indicator_mass =
        std::numeric_limits<double>::quiet_NaN();
    double hybrid_marked_H_indicator_mass =
        std::numeric_limits<double>::quiet_NaN();
    int hybrid_coarse_conformity_collar = -1;
    bool hybrid_coarse_marking_closure_safe = false;
    bool hybrid_full_regular_doerfler = false;
    std::size_t hybrid_coarse_preview_attempts = 0;
    bool hybrid_coarse_preview_cached = false;
    std::size_t marked_c = 0;
    std::size_t marked_c_f = 0;
    std::size_t marked_c_r = 0;
    bool candidate_accuracy_sweep_performed = false;
    std::size_t candidate_elements_before = 0;
    std::size_t candidate_containment_requested_marks = 0;
    std::size_t candidate_containment_added_elements = 0;
    std::size_t candidate_containment_closure_added_elements = 0;
    std::size_t candidate_accuracy_requested_marks = 0;
    std::size_t candidate_accuracy_added_elements = 0;
    std::size_t candidate_accuracy_closure_added_elements = 0;
    std::size_t candidate_elements_after = 0;
    double time_candidate_marking = 0.0;
    std::size_t candidate_marking_pool = 0;
    std::size_t candidate_estimated_selected_closure_cost = 0;
    std::size_t active_correctors = 0;
    std::size_t rebuilt_correctors = 0;
    std::size_t reused_correctors = 0;
    std::size_t corrector_cache_oversized_misses = 0;
    std::size_t corrector_cache_budget_rejections = 0;
    std::size_t corrector_cache_entries = 0;
    std::size_t corrector_cache_current_bytes = 0;
    std::size_t corrector_cache_peak_bytes = 0;
    std::size_t skipped_correctors = 0;
    std::size_t skipped_corrector_work_units = 0;
    int corrector_parallel_threads = 1;
    double corrector_patch_assembly_work_seconds = 0.0;
    double corrector_patch_solve_work_seconds = 0.0;
    double corrector_patch_pack_work_seconds = 0.0;
    int corrector_maximum_patch_dofs = 0;
    int corrector_maximum_patch_constraints = 0;
    int corrector_maximum_patch_rhs = 0;
    int hybrid_l_s = -1;
    double hybrid_minimum_physical_radius =
        std::numeric_limits<double>::quiet_NaN();
    double hybrid_covered_physical_radius =
        std::numeric_limits<double>::quiet_NaN();
    std::size_t hybrid_omega_s_elements = 0;
    std::size_t hybrid_omega_f_elements = 0;
    std::size_t candidate_cells_f = 0;
    std::size_t candidate_cells_r = 0;
    double time_reference_stability = 0.0;
    double time_corrector_check_total = 0.0;
    double time_lod_build_total = 0.0;
    double time_lod_mesh_and_interpolation = 0.0;
    double time_lod_operators = 0.0;
    double time_corrector = 0.0;
    double time_lod_corrected_basis = 0.0;
    double time_lod_coarse_operator = 0.0;
    double time_lod_coarse_factorization = 0.0;
    double time_theta = 0.0;
    double time_gram_prepare_structure = 0.0;
    double time_gram_prepare_factorization = 0.0;
    double time_gram_action_rhs = 0.0;
    double time_gram_action_patch_solve = 0.0;
    double time_gram_action_scatter = 0.0;
    std::size_t gram_action_calls = 0;
    std::size_t gram_patch_factorizations = 0;
    std::size_t gram_factor_cache_hits = 0;
    std::size_t gram_factor_cache_misses = 0;
    int gram_structure_parallel_threads = 1;
    int gram_parallel_threads = 1;
    int localization_iterations = 0;
    double localization_relative_residual = 0.0;
    bool localization_used_warm_start = false;
    double time_lod_solve = 0.0;
    double time_reference_riesz = 0.0;
    double time_hybrid_coarse_marking = 0.0;
    double time_hybrid_coarse_preview = 0.0;
    double time_hybrid_coarse_preview_nvb = 0.0;
    double time_hybrid_coarse_preview_reference_embedding = 0.0;
    double time_reference_validation = 0.0;
    double time_exact_validation = 0.0;
    double time_candidate_flux = 0.0;
    double time_candidate_close = 0.0;
    double time_candidate_operator_assembly = 0.0;
    double time_candidate_prolongation = 0.0;
    double time_candidate_flux_reconstruction = 0.0;
    double time_candidate_flux_prepare = 0.0;
    double time_candidate_flux_patch_solve = 0.0;
    double time_candidate_flux_merge = 0.0;
    double time_candidate_flux_audit = 0.0;
    int candidate_flux_parallel_threads = 1;
    double time_candidate_enrich = 0.0;
    double time_candidate_nvb_refine = 0.0;
    double time_candidate_embedding_composition = 0.0;
    double time_candidate_parent_map_update = 0.0;
    double time_candidate_quasi_interpolation = 0.0;
    double time_candidate_embedding_validation = 0.0;
    double time_candidate_dual = 0.0;
    double time_candidate_wellposedness = 0.0;
    double time_candidate_dual_operator_assembly = 0.0;
    double time_candidate_dual_load_assembly = 0.0;
    double time_candidate_dual_prolongation = 0.0;
    double time_candidate_dual_solve = 0.0;
    double time_candidate_dual_prepare = 0.0;
    double time_candidate_dual_patch_solve = 0.0;
    double time_candidate_dual_reduction = 0.0;
    std::size_t candidate_dual_patch_factorizations = 0;
    int candidate_dual_parallel_threads = 1;
    double time_mesh = 0.0;
    double time_validation_cumulative = 0.0;
    double time_artifact_capture_cumulative = 0.0;
    double time_method_cumulative = 0.0;
    double time_total_cumulative = 0.0;
    bool structural_dual_trigger = false;
    bool level_gap_dual_trigger = false;
    int minimum_reference_level_gap =
        std::numeric_limits<int>::max();
    bool tolerance_dual_trigger = false;
    bool decrease_dual_trigger = false;
    bool interval_dual_trigger = false;
    std::string detail;
};

struct ReferenceEpochDriverResult {
    ReferenceEpochDriverState state = ReferenceEpochDriverState::Failed;
    std::string stop_reason;
    std::size_t epochs = 0;
    std::size_t H_steps = 0;
    std::size_t dual_checks = 0;
    int ell = 0;
    std::vector<ReferenceEpochDriverRecord> journal;
};

class ReferenceEpochPracticalDriver {
public:
    ReferenceEpochPracticalDriver(
        ReferenceEpochDriverBackend &backend,
        ReferenceEpochDriverConfig config);

    ReferenceEpochDriverResult run();

private:
    ReferenceEpochDriverBackend &backend_;
    ReferenceEpochDriverConfig config_;
};

const char *reference_epoch_driver_state_name(ReferenceEpochDriverState state);
const char *reference_epoch_driver_action_name(ReferenceEpochDriverAction action);

} // namespace lod2d::helmholtz::adaptive

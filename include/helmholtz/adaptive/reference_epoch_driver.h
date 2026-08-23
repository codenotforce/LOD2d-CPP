#pragma once

#include <chrono>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace lod2d::helmholtz::adaptive {

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
    std::size_t rebuilt_correctors = 0;
    std::size_t skipped_correctors = 0;
    std::size_t skipped_corrector_work_units = 0;
    double time_corrector = 0.0;
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
    int gram_parallel_threads = 1;
    int localization_iterations = 0;
    double localization_relative_residual = 0.0;
    bool localization_used_warm_start = false;
};

struct ReferenceEpochSolveObservation {
    double eta_H = 0.0;
    double U_practical = 0.0;
    std::vector<int> marked_H;
    double relative_reference_energy =
        std::numeric_limits<double>::quiet_NaN();
    double relative_exact_energy =
        std::numeric_limits<double>::quiet_NaN();
    double relative_exact_L2 =
        std::numeric_limits<double>::quiet_NaN();
    double time_lod_solve = 0.0;
    double time_reference_riesz = 0.0;
};

struct ReferenceEpochCandidateObservation {
    double eta_eq_c = 0.0;
    std::vector<int> marked_c;
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
    virtual ReferenceEpochCandidateObservation enrich_candidate() = 0;
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
    // Suppress only the proactive level-gap guard until this many coarse
    // commits have occurred in the current epoch. Structural non-containment
    // and a positive dual-gap certificate remain immediate.
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
    double eta_dual_c = std::numeric_limits<double>::quiet_NaN();
    double L_gap_c = std::numeric_limits<double>::quiet_NaN();
    double relative_reference_energy = std::numeric_limits<double>::quiet_NaN();
    double relative_exact_energy = std::numeric_limits<double>::quiet_NaN();
    double relative_exact_L2 = std::numeric_limits<double>::quiet_NaN();
    std::size_t coarse_unknowns = 0;
    std::size_t reference_unknowns = 0;
    std::size_t candidate_unknowns = 0;
    double kappa_H_max = std::numeric_limits<double>::quiet_NaN();
    std::size_t marked_H = 0;
    std::size_t marked_c = 0;
    std::size_t rebuilt_correctors = 0;
    std::size_t skipped_correctors = 0;
    std::size_t skipped_corrector_work_units = 0;
    double time_corrector = 0.0;
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
    int gram_parallel_threads = 1;
    int localization_iterations = 0;
    double localization_relative_residual = 0.0;
    bool localization_used_warm_start = false;
    double time_lod_solve = 0.0;
    double time_reference_riesz = 0.0;
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

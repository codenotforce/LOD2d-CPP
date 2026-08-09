#pragma once

#include "helmholtz/adaptive/certificates.h"
#include "helmholtz/adaptive/kernel_residual.h"

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace lod2d::helmholtz::adaptive {

inline constexpr int certified_driver_checkpoint_version = 1;

enum class CertifiedMethod {
    Calod,
    Hlod
};

enum class CertifiedErrorTarget {
    AuditSpace,
    Continuous
};

enum class CertifiedEvidencePolicy {
    RequireVerified,
    AllowConditional
};

enum class CertifiedEvidenceLevel {
    Verified,
    Conditional
};

enum class CertifiedDriverState {
    CoarseAdmissibility,
    CorrectorCertification,
    CoarseErrorControl,
    AuditControl,
    Done,
    WorkLimit,
    Failure
};

enum class CertifiedDriverAction {
    RefineCoarseAdmissibility,
    AcceptCoarseAdmissibility,
    RefineCorrectorFine,
    IncreaseOversampling,
    AcceptCorrectorCertificate,
    FormPendingCoarseMarking,
    RefineAudit,
    ApplyPendingCoarseMarking,
    CompleteAuditTolerance,
    CompleteContinuousTolerance,
    StopAtWorkLimit,
    Fail
};

enum class CertifiedMutationKind {
    RefineCoarse,
    RefineCorrectorFine,
    IncreaseOversampling,
    RefineAudit
};

// Terminal codes are deliberately one-to-one with terminal causes.  A run has
// exactly one code and never falls back from CALOD to the legacy H-only proxy.
enum class CertifiedStopCode {
    AuditToleranceReached,
    ContinuousToleranceReached,
    StateTransitionLimit,
    CoarseRefinementLimit,
    CorrectorRefinementLimit,
    OversamplingLimit,
    AuditRefinementLimit,
    CoarseDofLimit,
    FineDofLimit,
    AuditDofLimit,
    BackendWorkLimit,
    MemoryLimit,
    TimeLimit,
    UnverifiedEvidence,
    InvalidObservation,
    FrozenHlodPriorMismatch,
    FrozenHlodCorrectorFailure,
    CoarseRefinementUnavailable,
    CorrectorRefinementUnavailable,
    OversamplingUnavailable,
    AuditRefinementUnavailable,
    EmptyMarking,
    BackendFailure,
    CheckpointIncompatible
};

struct CertifiedWorkLimits {
    // Zero means unlimited, except max_state_transitions which must be positive.
    std::uint64_t max_state_transitions = 1000;
    std::uint64_t max_coarse_refinements = 0;
    std::uint64_t max_corrector_refinements = 0;
    std::uint64_t max_oversampling_increments = 0;
    std::uint64_t max_audit_refinements = 0;
    std::uint64_t max_coarse_dofs = 0;
    std::uint64_t max_fine_dofs = 0;
    std::uint64_t max_audit_dofs = 0;
    std::uint64_t max_backend_work_units = 0;
    std::uint64_t max_peak_memory_bytes = 0;
    double max_elapsed_seconds = 0.0;
    int max_oversampling = -1;
};

struct CertifiedDriverConfig {
    CertifiedMethod method = CertifiedMethod::Calod;
    CertifiedErrorTarget error_target = CertifiedErrorTarget::AuditSpace;
    CertifiedEvidencePolicy evidence_policy =
        CertifiedEvidencePolicy::RequireVerified;
    double mu0 = 0.5;
    double q0 = 0.25;
    double tau = 0.5;
    double theta_H = 0.5;
    double theta_h = 0.5;
    double rho_aud = 0.05;
    double tolerance = 0.1;
    CertifiedWorkLimits limits;

    // Required for HLOD.  These identifiers come from the frozen
    // baseline_parameters.csv entry and are checked before every decision.
    std::string hlod_prior_corrector_space_id;
    int hlod_prior_oversampling = -1;
};

struct CertifiedObservationEvidence {
    bool valid = true;
    CertifiedEvidenceLevel level = CertifiedEvidenceLevel::Conditional;
    std::string source;
    std::string hash;
    std::string invalid_reason;
};

struct CoarseAdmissibilityObservation {
    CertifiedObservationEvidence evidence;
    std::vector<double> mu_by_element;
    bool refinement_available = true;
};

struct CorrectorCertificationObservation {
    CertifiedObservationEvidence evidence;
    bool stability_condition_holds = false;
    double q_total_upper = 0.0;
    double delta_total_lower = 0.0;
    double delta_h_upper = 0.0;
    std::vector<double> eta_h_element_squared;
    bool corrector_refinement_available = true;
    bool oversampling_increment_available = true;
};

struct CoarseErrorObservation {
    CertifiedObservationEvidence evidence;
    double eta_H = 0.0;
    double lod_error_lower = 0.0;
    double lod_error_upper = 0.0;
    std::vector<double> eta_H_element_squared;
    bool coarse_refinement_available = true;
};

enum class AuditIntervalKind {
    VerifiedExternalEstimator,
    EmpiricalSaturation
};

struct AuditControlObservation {
    CertifiedObservationEvidence evidence;
    AuditIntervalKind interval_kind =
        AuditIntervalKind::VerifiedExternalEstimator;
    double audit_error_lower = 0.0;
    double audit_error_upper = 0.0;
    std::vector<int> marked_fine_elements;
    bool refinement_available = true;
};

struct CertifiedWorkSnapshot {
    std::uint64_t coarse_dofs = 0;
    std::uint64_t fine_dofs = 0;
    std::uint64_t audit_dofs = 0;
    std::uint64_t backend_work_units = 0;
    std::uint64_t peak_memory_bytes = 0;
    double elapsed_seconds = 0.0;
    std::string initial_problem_id;
    std::string state_fingerprint;
    std::string corrector_space_id;
    int oversampling = -1;
};

// This is the entire algorithm-facing surface.  It intentionally has no exact
// solution, evaluation-reference, or ErrorReference accessor.
class CertifiedDriverBackend {
public:
    virtual ~CertifiedDriverBackend() = default;

    virtual CertifiedWorkSnapshot work_snapshot() const = 0;
    virtual CoarseAdmissibilityObservation inspect_coarse_admissibility() = 0;
    virtual void refine_coarse(const std::vector<int> &marked_elements) = 0;
    virtual CorrectorCertificationObservation inspect_corrector_certification() = 0;
    virtual void refine_corrector_patches(
        const std::vector<int> &marked_coarse_sources) = 0;
    virtual void increase_oversampling() = 0;
    virtual CoarseErrorObservation solve_and_estimate_coarse_error() = 0;
    virtual AuditControlObservation inspect_audit_control() = 0;
    virtual void refine_audit(const std::vector<int> &marked_fine_elements) = 0;
};

struct CertifiedMutation {
    CertifiedMutationKind kind = CertifiedMutationKind::RefineCoarse;
    std::vector<int> marked_elements;
};

struct CertifiedTransitionRecord {
    std::uint64_t sequence = 0;
    CertifiedDriverState before = CertifiedDriverState::CoarseAdmissibility;
    CertifiedDriverState after = CertifiedDriverState::CoarseAdmissibility;
    CertifiedDriverAction action =
        CertifiedDriverAction::AcceptCoarseAdmissibility;
    CertifiedEvidenceLevel claim = CertifiedEvidenceLevel::Verified;
    std::vector<int> marked_elements;
    std::optional<double> mu_max;
    std::optional<double> q_total_upper;
    std::optional<double> delta_total_lower;
    std::optional<double> delta_h_upper;
    std::optional<double> lod_error_lower;
    std::optional<double> lod_error_upper;
    std::optional<double> audit_error_lower;
    std::optional<double> audit_error_upper;
    std::optional<double> true_error_lower;
    std::optional<double> true_error_upper;
    std::string detail;
};

struct CertifiedTermination {
    CertifiedDriverState state = CertifiedDriverState::Failure;
    CertifiedStopCode code = CertifiedStopCode::InvalidObservation;
    CertifiedEvidenceLevel claim = CertifiedEvidenceLevel::Verified;
    std::string detail;
};

struct CertifiedDriverCheckpoint {
    int schema_version = certified_driver_checkpoint_version;
    std::string config_fingerprint;
    CertifiedDriverState state = CertifiedDriverState::CoarseAdmissibility;
    CertifiedEvidenceLevel claim = CertifiedEvidenceLevel::Verified;
    std::uint64_t transition_count = 0;
    std::vector<int> pending_coarse_marking;
    bool pending_coarse_refinement_available = true;
    std::optional<double> lod_error_lower;
    std::optional<double> lod_error_upper;
    std::optional<double> audit_error_lower;
    std::optional<double> audit_error_upper;
    std::optional<double> true_error_lower;
    std::optional<double> true_error_upper;
    std::string initial_problem_id;
    std::string backend_fingerprint;
    std::string frozen_corrector_space_id;
    int frozen_oversampling = -1;
    std::vector<CertifiedMutation> mutations;
    std::vector<CertifiedTransitionRecord> history;
    std::optional<CertifiedTermination> termination;
};

class CertifiedAdaptiveDriver {
public:
    explicit CertifiedAdaptiveDriver(CertifiedDriverConfig config);

    // Resume reconstructs the numerical backend by replaying only recorded
    // mesh/oversampling mutations, then checks its deterministic fingerprint.
    static CertifiedAdaptiveDriver resume(
        CertifiedDriverConfig config,
        const CertifiedDriverCheckpoint &checkpoint,
        CertifiedDriverBackend &fresh_backend);
    static CertifiedAdaptiveDriver resume(
        CertifiedDriverConfig config,
        std::string_view serialized_checkpoint,
        CertifiedDriverBackend &fresh_backend);

    // Executes one paper-state decision.  Returns false after the unique
    // terminal outcome has been recorded.
    bool step(CertifiedDriverBackend &backend);
    void run(CertifiedDriverBackend &backend);

    bool terminal() const;
    std::string_view output_namespace() const;
    const CertifiedDriverConfig &config() const { return config_; }
    const CertifiedDriverCheckpoint &checkpoint() const { return checkpoint_; }
    const std::optional<CertifiedTermination> &termination() const {
        return checkpoint_.termination;
    }

private:
    CertifiedAdaptiveDriver(
        CertifiedDriverConfig config,
        CertifiedDriverCheckpoint checkpoint);

    CertifiedDriverConfig config_;
    CertifiedDriverCheckpoint checkpoint_;
};

std::string serialize_certified_checkpoint(
    const CertifiedDriverCheckpoint &checkpoint);
CertifiedDriverCheckpoint deserialize_certified_checkpoint(
    std::string_view serialized);

// Canonical adapters from WP3/WP4 results.  They do not accept evaluation
// reference data and therefore cannot leak exact/reference errors into MARK/STOP.
CoarseAdmissibilityObservation make_coarse_admissibility_observation(
    const TriMesh &coarse_mesh,
    double wavenumber,
    const CertificateConstantRegistry &constants);
CorrectorCertificationObservation make_corrector_certification_observation(
    const CorrectorCertificateResult &certificate);
CoarseErrorObservation make_coarse_error_observation(
    const AuditKernelResidualEstimate &estimate,
    const CorrectorCertificateResult &certificate);

std::string_view to_string(CertifiedMethod value);
std::string_view to_string(CertifiedErrorTarget value);
std::string_view to_string(CertifiedEvidencePolicy value);
std::string_view to_string(CertifiedEvidenceLevel value);
std::string_view to_string(CertifiedDriverState value);
std::string_view to_string(CertifiedDriverAction value);
std::string_view to_string(CertifiedMutationKind value);
std::string_view to_string(CertifiedStopCode value);
std::string_view to_string(AuditIntervalKind value);

} // namespace lod2d::helmholtz::adaptive

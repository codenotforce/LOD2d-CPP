#include "helmholtz/adaptive/certified_driver.h"

#include "helmholtz/adaptive/estimator.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace lod2d::helmholtz::adaptive {

std::string_view to_string(CertifiedMethod value) {
    switch (value) {
    case CertifiedMethod::Calod: return "calod";
    case CertifiedMethod::Hlod: return "hlod";
    }
    throw std::invalid_argument("unknown certified method");
}

std::string_view to_string(CertifiedErrorTarget value) {
    switch (value) {
    case CertifiedErrorTarget::AuditSpace: return "audit_space";
    case CertifiedErrorTarget::Continuous: return "continuous";
    }
    throw std::invalid_argument("unknown certified error target");
}

std::string_view to_string(CertifiedEvidencePolicy value) {
    switch (value) {
    case CertifiedEvidencePolicy::RequireVerified: return "require_verified";
    case CertifiedEvidencePolicy::AllowConditional: return "allow_conditional";
    }
    throw std::invalid_argument("unknown certified evidence policy");
}

std::string_view to_string(CertifiedEvidenceLevel value) {
    switch (value) {
    case CertifiedEvidenceLevel::Verified: return "verified";
    case CertifiedEvidenceLevel::Conditional: return "conditional";
    }
    throw std::invalid_argument("unknown certified evidence level");
}

std::string_view to_string(CertifiedDriverState value) {
    switch (value) {
    case CertifiedDriverState::CoarseAdmissibility: return "coarse_admissibility";
    case CertifiedDriverState::CorrectorCertification: return "corrector_certification";
    case CertifiedDriverState::CoarseErrorControl: return "coarse_error_control";
    case CertifiedDriverState::AuditControl: return "audit_control";
    case CertifiedDriverState::Done: return "done";
    case CertifiedDriverState::WorkLimit: return "work_limit";
    case CertifiedDriverState::Failure: return "failure";
    }
    throw std::invalid_argument("unknown certified driver state");
}

std::string_view to_string(CertifiedDriverAction value) {
    switch (value) {
    case CertifiedDriverAction::RefineCoarseAdmissibility:
        return "refine_coarse_admissibility";
    case CertifiedDriverAction::AcceptCoarseAdmissibility:
        return "accept_coarse_admissibility";
    case CertifiedDriverAction::RefineCorrectorFine:
        return "refine_corrector_fine";
    case CertifiedDriverAction::IncreaseOversampling:
        return "increase_oversampling";
    case CertifiedDriverAction::AcceptCorrectorCertificate:
        return "accept_corrector_certificate";
    case CertifiedDriverAction::FormPendingCoarseMarking:
        return "form_pending_coarse_marking";
    case CertifiedDriverAction::RefineAudit: return "refine_audit";
    case CertifiedDriverAction::ApplyPendingCoarseMarking:
        return "apply_pending_coarse_marking";
    case CertifiedDriverAction::CompleteAuditTolerance:
        return "complete_audit_tolerance";
    case CertifiedDriverAction::CompleteContinuousTolerance:
        return "complete_continuous_tolerance";
    case CertifiedDriverAction::StopAtWorkLimit: return "stop_at_work_limit";
    case CertifiedDriverAction::Fail: return "fail";
    }
    throw std::invalid_argument("unknown certified driver action");
}

std::string_view to_string(CertifiedMutationKind value) {
    switch (value) {
    case CertifiedMutationKind::RefineCoarse: return "refine_coarse";
    case CertifiedMutationKind::RefineCorrectorFine:
        return "refine_corrector_fine";
    case CertifiedMutationKind::IncreaseOversampling:
        return "increase_oversampling";
    case CertifiedMutationKind::RefineAudit: return "refine_audit";
    }
    throw std::invalid_argument("unknown certified mutation kind");
}

std::string_view to_string(CertifiedStopCode value) {
    switch (value) {
    case CertifiedStopCode::AuditToleranceReached:
        return "audit_tolerance_reached";
    case CertifiedStopCode::ContinuousToleranceReached:
        return "continuous_tolerance_reached";
    case CertifiedStopCode::StateTransitionLimit:
        return "state_transition_limit";
    case CertifiedStopCode::CoarseRefinementLimit:
        return "coarse_refinement_limit";
    case CertifiedStopCode::CorrectorRefinementLimit:
        return "corrector_refinement_limit";
    case CertifiedStopCode::OversamplingLimit: return "oversampling_limit";
    case CertifiedStopCode::AuditRefinementLimit:
        return "audit_refinement_limit";
    case CertifiedStopCode::CoarseDofLimit: return "coarse_dof_limit";
    case CertifiedStopCode::FineDofLimit: return "fine_dof_limit";
    case CertifiedStopCode::AuditDofLimit: return "audit_dof_limit";
    case CertifiedStopCode::BackendWorkLimit: return "backend_work_limit";
    case CertifiedStopCode::MemoryLimit: return "memory_limit";
    case CertifiedStopCode::TimeLimit: return "time_limit";
    case CertifiedStopCode::UnverifiedEvidence: return "unverified_evidence";
    case CertifiedStopCode::InvalidObservation: return "invalid_observation";
    case CertifiedStopCode::FrozenHlodPriorMismatch:
        return "frozen_hlod_prior_mismatch";
    case CertifiedStopCode::FrozenHlodCorrectorFailure:
        return "frozen_hlod_corrector_failure";
    case CertifiedStopCode::CoarseRefinementUnavailable:
        return "coarse_refinement_unavailable";
    case CertifiedStopCode::CorrectorRefinementUnavailable:
        return "corrector_refinement_unavailable";
    case CertifiedStopCode::OversamplingUnavailable:
        return "oversampling_unavailable";
    case CertifiedStopCode::AuditRefinementUnavailable:
        return "audit_refinement_unavailable";
    case CertifiedStopCode::EmptyMarking: return "empty_marking";
    case CertifiedStopCode::BackendFailure: return "backend_failure";
    case CertifiedStopCode::CheckpointIncompatible:
        return "checkpoint_incompatible";
    }
    throw std::invalid_argument("unknown certified stop code");
}

std::string_view to_string(AuditIntervalKind value) {
    switch (value) {
    case AuditIntervalKind::VerifiedExternalEstimator:
        return "verified_external_estimator";
    case AuditIntervalKind::EmpiricalSaturation:
        return "empirical_saturation";
    }
    throw std::invalid_argument("unknown audit interval kind");
}

namespace {

template <class Enum>
Enum parse_enum(std::string_view text);

template <>
CertifiedEvidenceLevel parse_enum(std::string_view text) {
    for (CertifiedEvidenceLevel value : {
             CertifiedEvidenceLevel::Verified,
             CertifiedEvidenceLevel::Conditional}) {
        if (to_string(value) == text) return value;
    }
    throw std::invalid_argument("unknown certified evidence level: "
                                + std::string(text));
}

template <>
CertifiedDriverState parse_enum(std::string_view text) {
    for (CertifiedDriverState value : {
             CertifiedDriverState::CoarseAdmissibility,
             CertifiedDriverState::CorrectorCertification,
             CertifiedDriverState::CoarseErrorControl,
             CertifiedDriverState::AuditControl,
             CertifiedDriverState::Done,
             CertifiedDriverState::WorkLimit,
             CertifiedDriverState::Failure}) {
        if (to_string(value) == text) return value;
    }
    throw std::invalid_argument("unknown certified driver state: "
                                + std::string(text));
}

template <>
CertifiedDriverAction parse_enum(std::string_view text) {
    for (CertifiedDriverAction value : {
             CertifiedDriverAction::RefineCoarseAdmissibility,
             CertifiedDriverAction::AcceptCoarseAdmissibility,
             CertifiedDriverAction::RefineCorrectorFine,
             CertifiedDriverAction::IncreaseOversampling,
             CertifiedDriverAction::AcceptCorrectorCertificate,
             CertifiedDriverAction::FormPendingCoarseMarking,
             CertifiedDriverAction::RefineAudit,
             CertifiedDriverAction::ApplyPendingCoarseMarking,
             CertifiedDriverAction::CompleteAuditTolerance,
             CertifiedDriverAction::CompleteContinuousTolerance,
             CertifiedDriverAction::StopAtWorkLimit,
             CertifiedDriverAction::Fail}) {
        if (to_string(value) == text) return value;
    }
    throw std::invalid_argument("unknown certified driver action: "
                                + std::string(text));
}

template <>
CertifiedMutationKind parse_enum(std::string_view text) {
    for (CertifiedMutationKind value : {
             CertifiedMutationKind::RefineCoarse,
             CertifiedMutationKind::RefineCorrectorFine,
             CertifiedMutationKind::IncreaseOversampling,
             CertifiedMutationKind::RefineAudit}) {
        if (to_string(value) == text) return value;
    }
    throw std::invalid_argument("unknown certified mutation kind: "
                                + std::string(text));
}

template <>
CertifiedStopCode parse_enum(std::string_view text) {
    for (CertifiedStopCode value : {
             CertifiedStopCode::AuditToleranceReached,
             CertifiedStopCode::ContinuousToleranceReached,
             CertifiedStopCode::StateTransitionLimit,
             CertifiedStopCode::CoarseRefinementLimit,
             CertifiedStopCode::CorrectorRefinementLimit,
             CertifiedStopCode::OversamplingLimit,
             CertifiedStopCode::AuditRefinementLimit,
             CertifiedStopCode::CoarseDofLimit,
             CertifiedStopCode::FineDofLimit,
             CertifiedStopCode::AuditDofLimit,
             CertifiedStopCode::BackendWorkLimit,
             CertifiedStopCode::MemoryLimit,
             CertifiedStopCode::TimeLimit,
             CertifiedStopCode::UnverifiedEvidence,
             CertifiedStopCode::InvalidObservation,
             CertifiedStopCode::FrozenHlodPriorMismatch,
             CertifiedStopCode::FrozenHlodCorrectorFailure,
             CertifiedStopCode::CoarseRefinementUnavailable,
             CertifiedStopCode::CorrectorRefinementUnavailable,
             CertifiedStopCode::OversamplingUnavailable,
             CertifiedStopCode::AuditRefinementUnavailable,
             CertifiedStopCode::EmptyMarking,
             CertifiedStopCode::BackendFailure,
             CertifiedStopCode::CheckpointIncompatible}) {
        if (to_string(value) == text) return value;
    }
    throw std::invalid_argument("unknown certified stop code: "
                                + std::string(text));
}

bool terminal_state(CertifiedDriverState state) {
    return state == CertifiedDriverState::Done
        || state == CertifiedDriverState::WorkLimit
        || state == CertifiedDriverState::Failure;
}

void validate_config(const CertifiedDriverConfig &config) {
    const auto finite = [](double value) { return std::isfinite(value); };
    if (!(finite(config.mu0) && config.mu0 > 0.0 && config.mu0 < 1.0))
        throw std::invalid_argument("mu0 must lie in (0,1)");
    if (!(finite(config.q0) && config.q0 > 0.0))
        throw std::invalid_argument("q0 must be positive and finite");
    if (!(finite(config.tau) && config.tau > 0.0 && config.tau < 1.0))
        throw std::invalid_argument("tau must lie in (0,1)");
    if (!(finite(config.theta_H) && config.theta_H > 0.0
          && config.theta_H <= 1.0))
        throw std::invalid_argument("theta_H must lie in (0,1]");
    if (!(finite(config.theta_h) && config.theta_h > 0.0
          && config.theta_h <= 1.0))
        throw std::invalid_argument("theta_h must lie in (0,1]");
    if (!(finite(config.rho_aud) && config.rho_aud > 0.0
          && config.rho_aud < 1.0))
        throw std::invalid_argument("rho_aud must lie in (0,1)");
    if (!(finite(config.tolerance) && config.tolerance > 0.0))
        throw std::invalid_argument("tolerance must be positive and finite");
    if (config.limits.max_state_transitions == 0)
        throw std::invalid_argument("max_state_transitions must be positive");
    if (!finite(config.limits.max_elapsed_seconds)
        || config.limits.max_elapsed_seconds < 0.0)
        throw std::invalid_argument("max_elapsed_seconds is invalid");
    if (config.limits.max_oversampling < -1)
        throw std::invalid_argument("max_oversampling is invalid");
    if (config.method == CertifiedMethod::Hlod
        && (config.hlod_prior_corrector_space_id.empty()
            || config.hlod_prior_oversampling < 0)) {
        throw std::invalid_argument(
            "HLOD requires a frozen corrector-space id and oversampling value");
    }
}

std::string config_fingerprint(const CertifiedDriverConfig &config) {
    std::ostringstream canonical;
    canonical << std::setprecision(std::numeric_limits<double>::max_digits10)
              << to_string(config.method) << '|'
              << to_string(config.error_target) << '|'
              << to_string(config.evidence_policy) << '|'
              << config.mu0 << '|' << config.q0 << '|' << config.tau << '|'
              << config.theta_H << '|' << config.theta_h << '|'
              << config.rho_aud << '|' << config.tolerance << '|'
              << config.limits.max_state_transitions << '|'
              << config.limits.max_coarse_refinements << '|'
              << config.limits.max_corrector_refinements << '|'
              << config.limits.max_oversampling_increments << '|'
              << config.limits.max_audit_refinements << '|'
              << config.limits.max_coarse_dofs << '|'
              << config.limits.max_fine_dofs << '|'
              << config.limits.max_audit_dofs << '|'
              << config.limits.max_backend_work_units << '|'
              << config.limits.max_peak_memory_bytes << '|'
              << config.limits.max_elapsed_seconds << '|'
              << config.limits.max_oversampling << '|'
              << std::quoted(config.hlod_prior_corrector_space_id) << '|'
              << config.hlod_prior_oversampling;
    std::uint64_t hash = 14695981039346656037ULL;
    for (unsigned char value : canonical.str()) {
        hash ^= value;
        hash *= 1099511628211ULL;
    }
    std::ostringstream encoded;
    encoded << std::hex << std::setfill('0') << std::setw(16) << hash;
    return encoded.str();
}

std::optional<std::string> validate_snapshot(
    const CertifiedWorkSnapshot &snapshot) {
    if (snapshot.initial_problem_id.empty())
        return "backend initial_problem_id is empty";
    if (snapshot.state_fingerprint.empty())
        return "backend state_fingerprint is empty";
    if (snapshot.corrector_space_id.empty())
        return "backend corrector_space_id is empty";
    if (snapshot.oversampling < 0)
        return "backend oversampling is negative";
    if (!std::isfinite(snapshot.elapsed_seconds)
        || snapshot.elapsed_seconds < 0.0)
        return "backend elapsed_seconds is invalid";
    return std::nullopt;
}

std::optional<std::pair<CertifiedStopCode, std::string>> resource_limit(
    const CertifiedWorkSnapshot &snapshot,
    const CertifiedWorkLimits &limits) {
    const auto exceeded = [](std::uint64_t value, std::uint64_t limit) {
        return limit != 0 && value > limit;
    };
    if (exceeded(snapshot.coarse_dofs, limits.max_coarse_dofs))
        return std::pair{CertifiedStopCode::CoarseDofLimit,
                         "coarse degree-of-freedom limit exceeded"};
    if (exceeded(snapshot.fine_dofs, limits.max_fine_dofs))
        return std::pair{CertifiedStopCode::FineDofLimit,
                         "corrector-fine degree-of-freedom limit exceeded"};
    if (exceeded(snapshot.audit_dofs, limits.max_audit_dofs))
        return std::pair{CertifiedStopCode::AuditDofLimit,
                         "audit degree-of-freedom limit exceeded"};
    if (exceeded(snapshot.backend_work_units, limits.max_backend_work_units))
        return std::pair{CertifiedStopCode::BackendWorkLimit,
                         "backend work-unit limit exceeded"};
    if (exceeded(snapshot.peak_memory_bytes, limits.max_peak_memory_bytes))
        return std::pair{CertifiedStopCode::MemoryLimit,
                         "peak-memory limit exceeded"};
    if (limits.max_elapsed_seconds > 0.0
        && snapshot.elapsed_seconds > limits.max_elapsed_seconds)
        return std::pair{CertifiedStopCode::TimeLimit,
                         "elapsed-time limit exceeded"};
    return std::nullopt;
}

std::uint64_t mutation_count(
    const CertifiedDriverCheckpoint &checkpoint,
    CertifiedMutationKind kind) {
    return static_cast<std::uint64_t>(std::count_if(
        checkpoint.mutations.begin(), checkpoint.mutations.end(),
        [kind](const CertifiedMutation &mutation) {
            return mutation.kind == kind;
        }));
}

bool finite_nonnegative(const std::vector<double> &values) {
    return std::all_of(values.begin(), values.end(), [](double value) {
        return std::isfinite(value) && value >= 0.0;
    });
}

std::vector<int> normalized_marking(const std::vector<int> &marked) {
    std::vector<int> result = marked;
    if (std::any_of(result.begin(), result.end(),
                    [](int value) { return value < 0; }))
        throw std::invalid_argument("marking contains a negative element index");
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

double element_diameter(const TriMesh &mesh, const Triangle &triangle) {
    return std::max({
        (mesh.nodes[triangle[0]] - mesh.nodes[triangle[1]]).norm(),
        (mesh.nodes[triangle[1]] - mesh.nodes[triangle[2]]).norm(),
        (mesh.nodes[triangle[2]] - mesh.nodes[triangle[0]]).norm()});
}

std::string join_reasons(const std::vector<std::string> &reasons) {
    std::ostringstream joined;
    for (std::size_t index = 0; index < reasons.size(); ++index) {
        if (index != 0) joined << "; ";
        joined << reasons[index];
    }
    return joined.str();
}

void expect_tag(std::istream &input, std::string_view expected) {
    std::string actual;
    if (!(input >> actual) || actual != expected)
        throw std::invalid_argument("checkpoint expected tag "
                                    + std::string(expected));
}

template <class T>
T read_value(std::istream &input, std::string_view name) {
    T result{};
    if (!(input >> result))
        throw std::invalid_argument("checkpoint cannot read "
                                    + std::string(name));
    return result;
}

std::string read_quoted(std::istream &input, std::string_view name) {
    std::string result;
    if (!(input >> std::quoted(result)))
        throw std::invalid_argument("checkpoint cannot read "
                                    + std::string(name));
    return result;
}

std::vector<int> read_int_vector(std::istream &input, std::string_view name) {
    const std::uint64_t size = read_value<std::uint64_t>(input, name);
    if (size > 100000000ULL)
        throw std::invalid_argument("checkpoint vector size is unreasonable");
    std::vector<int> result(static_cast<std::size_t>(size));
    for (int &value : result) value = read_value<int>(input, name);
    return result;
}

void write_optional(std::ostream &output, const std::optional<double> &value) {
    output << ' ' << (value.has_value() ? 1 : 0);
    if (value) output << ' ' << *value;
}

std::optional<double> read_optional(std::istream &input, std::string_view name) {
    const int present = read_value<int>(input, name);
    if (present == 0) return std::nullopt;
    if (present != 1)
        throw std::invalid_argument("checkpoint optional presence is not 0 or 1");
    const double value = read_value<double>(input, name);
    if (!std::isfinite(value))
        throw std::invalid_argument("checkpoint optional value is not finite");
    return value;
}

bool success_stop_code(CertifiedStopCode code) {
    return code == CertifiedStopCode::AuditToleranceReached
        || code == CertifiedStopCode::ContinuousToleranceReached;
}

bool work_limit_stop_code(CertifiedStopCode code) {
    switch (code) {
    case CertifiedStopCode::StateTransitionLimit:
    case CertifiedStopCode::CoarseRefinementLimit:
    case CertifiedStopCode::CorrectorRefinementLimit:
    case CertifiedStopCode::OversamplingLimit:
    case CertifiedStopCode::AuditRefinementLimit:
    case CertifiedStopCode::CoarseDofLimit:
    case CertifiedStopCode::FineDofLimit:
    case CertifiedStopCode::AuditDofLimit:
    case CertifiedStopCode::BackendWorkLimit:
    case CertifiedStopCode::MemoryLimit:
    case CertifiedStopCode::TimeLimit:
        return true;
    default:
        return false;
    }
}

void validate_checkpoint_shape(const CertifiedDriverCheckpoint &checkpoint) {
    if (checkpoint.schema_version != certified_driver_checkpoint_version)
        throw std::invalid_argument("unsupported certified checkpoint version");
    if (checkpoint.config_fingerprint.empty())
        throw std::invalid_argument("checkpoint config fingerprint is empty");
    if (!std::isfinite(checkpoint.cumulative_elapsed_seconds)
        || checkpoint.cumulative_elapsed_seconds < 0.0) {
        throw std::invalid_argument(
            "checkpoint cumulative elapsed time is invalid");
    }
    if (terminal_state(checkpoint.state) != checkpoint.termination.has_value())
        throw std::invalid_argument(
            "checkpoint terminal state and termination payload disagree");
    if (checkpoint.termination
        && checkpoint.termination->state != checkpoint.state)
        throw std::invalid_argument(
            "checkpoint termination state does not match driver state");
    if (checkpoint.termination
        && checkpoint.termination->claim != checkpoint.claim)
        throw std::invalid_argument(
            "checkpoint termination claim does not match driver claim");
    if (checkpoint.termination) {
        const bool code_matches_state =
            (checkpoint.state == CertifiedDriverState::Done
             && success_stop_code(checkpoint.termination->code))
            || (checkpoint.state == CertifiedDriverState::WorkLimit
                && work_limit_stop_code(checkpoint.termination->code))
            || (checkpoint.state == CertifiedDriverState::Failure
                && !success_stop_code(checkpoint.termination->code)
                && !work_limit_stop_code(checkpoint.termination->code));
        if (!code_matches_state)
            throw std::invalid_argument(
                "checkpoint stop code does not match its terminal state");
    }
    if (checkpoint.history.size() != checkpoint.transition_count)
        throw std::invalid_argument(
            "checkpoint transition count does not match history size");
    CertifiedDriverState replayed_state =
        CertifiedDriverState::CoarseAdmissibility;
    for (std::size_t index = 0; index < checkpoint.history.size(); ++index) {
        if (checkpoint.history[index].sequence != index)
            throw std::invalid_argument("checkpoint history sequence is not contiguous");
        if (checkpoint.history[index].before != replayed_state)
            throw std::invalid_argument("checkpoint history state chain is discontinuous");
        replayed_state = checkpoint.history[index].after;
    }
    if (replayed_state != checkpoint.state)
        throw std::invalid_argument(
            "checkpoint history does not end at the recorded driver state");
    for (int element : checkpoint.pending_coarse_marking)
        if (element < 0)
            throw std::invalid_argument("checkpoint pending marking is invalid");
    for (const CertifiedMutation &mutation : checkpoint.mutations) {
        for (int element : mutation.marked_elements)
            if (element < 0)
                throw std::invalid_argument("checkpoint mutation marking is invalid");
        if (mutation.kind == CertifiedMutationKind::IncreaseOversampling) {
            if (!mutation.marked_elements.empty())
                throw std::invalid_argument(
                    "checkpoint oversampling mutation carries element indices");
        } else if (mutation.marked_elements.empty()) {
            throw std::invalid_argument(
                "checkpoint refinement mutation has an empty marking");
        }
    }
    const auto validate_interval = [](
        const std::optional<double> &lower,
        const std::optional<double> &upper,
        const char *name) {
        if (lower.has_value() != upper.has_value()) {
            throw std::invalid_argument(
                std::string("checkpoint ") + name
                + " interval has only one endpoint");
        }
        if (!lower) return;
        if (!std::isfinite(*lower) || !std::isfinite(*upper)
            || *lower < 0.0 || *upper < *lower) {
            throw std::invalid_argument(
                std::string("checkpoint ") + name
                + " interval is invalid");
        }
    };
    validate_interval(
        checkpoint.lod_error_lower, checkpoint.lod_error_upper, "LOD");
    validate_interval(
        checkpoint.audit_error_lower, checkpoint.audit_error_upper, "audit");
    validate_interval(
        checkpoint.true_error_lower, checkpoint.true_error_upper, "true");
    if (checkpoint.audit_error_lower && !checkpoint.lod_error_lower)
        throw std::invalid_argument(
            "checkpoint audit interval has no LOD interval");
    if (checkpoint.true_error_lower
        && (!checkpoint.lod_error_lower || !checkpoint.audit_error_lower)) {
        throw std::invalid_argument(
            "checkpoint true interval lacks a component interval");
    }
    if (checkpoint.state == CertifiedDriverState::AuditControl
        && !checkpoint.lod_error_lower) {
        throw std::invalid_argument(
            "AuditControl checkpoint has no stored LOD interval");
    }
}

} // namespace

std::string serialize_certified_checkpoint(
    const CertifiedDriverCheckpoint &checkpoint) {
    validate_checkpoint_shape(checkpoint);
    std::ostringstream output;
    output << std::setprecision(std::numeric_limits<double>::max_digits10);
    output << "LOD2D_CERTIFIED_DRIVER_CHECKPOINT "
           << checkpoint.schema_version << '\n';
    output << "config_fingerprint "
           << std::quoted(checkpoint.config_fingerprint) << '\n';
    output << "state " << to_string(checkpoint.state) << '\n';
    output << "claim " << to_string(checkpoint.claim) << '\n';
    output << "transition_count " << checkpoint.transition_count << '\n';
    output << "pending " << checkpoint.pending_coarse_marking.size();
    for (int element : checkpoint.pending_coarse_marking)
        output << ' ' << element;
    output << '\n';
    output << "pending_available "
           << (checkpoint.pending_coarse_refinement_available ? 1 : 0) << '\n';
    output << "intervals";
    write_optional(output, checkpoint.lod_error_lower);
    write_optional(output, checkpoint.lod_error_upper);
    write_optional(output, checkpoint.audit_error_lower);
    write_optional(output, checkpoint.audit_error_upper);
    write_optional(output, checkpoint.true_error_lower);
    write_optional(output, checkpoint.true_error_upper);
    output << '\n';
    output << "initial_problem_id "
           << std::quoted(checkpoint.initial_problem_id) << '\n';
    output << "backend_fingerprint "
           << std::quoted(checkpoint.backend_fingerprint) << '\n';
    output << "work_accounting "
           << checkpoint.cumulative_backend_work_units << ' '
           << checkpoint.cumulative_peak_memory_bytes << ' '
           << checkpoint.cumulative_elapsed_seconds << '\n';
    output << "frozen_corrector_space_id "
           << std::quoted(checkpoint.frozen_corrector_space_id) << '\n';
    output << "frozen_oversampling "
           << checkpoint.frozen_oversampling << '\n';
    output << "termination " << (checkpoint.termination ? 1 : 0);
    if (checkpoint.termination) {
        output << ' ' << to_string(checkpoint.termination->state)
               << ' ' << to_string(checkpoint.termination->code)
               << ' ' << to_string(checkpoint.termination->claim)
               << ' ' << std::quoted(checkpoint.termination->detail);
    }
    output << '\n';
    output << "mutations " << checkpoint.mutations.size() << '\n';
    for (const CertifiedMutation &mutation : checkpoint.mutations) {
        output << "mutation " << to_string(mutation.kind)
               << ' ' << mutation.marked_elements.size();
        for (int element : mutation.marked_elements) output << ' ' << element;
        output << '\n';
    }
    output << "history " << checkpoint.history.size() << '\n';
    for (const CertifiedTransitionRecord &record : checkpoint.history) {
        output << "record " << record.sequence
               << ' ' << to_string(record.before)
               << ' ' << to_string(record.after)
               << ' ' << to_string(record.action)
               << ' ' << to_string(record.claim)
               << ' ' << record.marked_elements.size();
        for (int element : record.marked_elements) output << ' ' << element;
        write_optional(output, record.mu_max);
        write_optional(output, record.q_total_upper);
        write_optional(output, record.delta_total_lower);
        write_optional(output, record.delta_h_upper);
        write_optional(output, record.lod_error_lower);
        write_optional(output, record.lod_error_upper);
        write_optional(output, record.audit_error_lower);
        write_optional(output, record.audit_error_upper);
        write_optional(output, record.true_error_lower);
        write_optional(output, record.true_error_upper);
        output << ' ' << std::quoted(record.detail) << '\n';
    }
    output << "end\n";
    return output.str();
}

CertifiedDriverCheckpoint deserialize_certified_checkpoint(
    std::string_view serialized) {
    std::istringstream input{std::string(serialized)};
    expect_tag(input, "LOD2D_CERTIFIED_DRIVER_CHECKPOINT");
    CertifiedDriverCheckpoint checkpoint;
    checkpoint.schema_version = read_value<int>(input, "schema version");
    if (checkpoint.schema_version != certified_driver_checkpoint_version) {
        throw std::invalid_argument(
            "unsupported certified checkpoint version; version 1 cannot be resumed because it has no cumulative resource accounting");
    }
    expect_tag(input, "config_fingerprint");
    checkpoint.config_fingerprint = read_quoted(input, "config fingerprint");
    expect_tag(input, "state");
    checkpoint.state = parse_enum<CertifiedDriverState>(
        read_value<std::string>(input, "state"));
    expect_tag(input, "claim");
    checkpoint.claim = parse_enum<CertifiedEvidenceLevel>(
        read_value<std::string>(input, "claim"));
    expect_tag(input, "transition_count");
    checkpoint.transition_count =
        read_value<std::uint64_t>(input, "transition count");
    expect_tag(input, "pending");
    checkpoint.pending_coarse_marking = read_int_vector(input, "pending marking");
    expect_tag(input, "pending_available");
    const int pending_available = read_value<int>(input, "pending availability");
    if (pending_available != 0 && pending_available != 1)
        throw std::invalid_argument("checkpoint pending availability is invalid");
    checkpoint.pending_coarse_refinement_available = pending_available != 0;
    expect_tag(input, "intervals");
    checkpoint.lod_error_lower = read_optional(input, "LOD lower interval");
    checkpoint.lod_error_upper = read_optional(input, "LOD upper interval");
    checkpoint.audit_error_lower = read_optional(input, "audit lower interval");
    checkpoint.audit_error_upper = read_optional(input, "audit upper interval");
    checkpoint.true_error_lower = read_optional(input, "true lower interval");
    checkpoint.true_error_upper = read_optional(input, "true upper interval");
    expect_tag(input, "initial_problem_id");
    checkpoint.initial_problem_id = read_quoted(input, "initial problem id");
    expect_tag(input, "backend_fingerprint");
    checkpoint.backend_fingerprint = read_quoted(input, "backend fingerprint");
    expect_tag(input, "work_accounting");
    checkpoint.cumulative_backend_work_units =
        read_value<std::uint64_t>(input, "cumulative backend work");
    checkpoint.cumulative_peak_memory_bytes =
        read_value<std::uint64_t>(input, "cumulative peak memory");
    checkpoint.cumulative_elapsed_seconds =
        read_value<double>(input, "cumulative elapsed seconds");
    expect_tag(input, "frozen_corrector_space_id");
    checkpoint.frozen_corrector_space_id =
        read_quoted(input, "frozen corrector-space id");
    expect_tag(input, "frozen_oversampling");
    checkpoint.frozen_oversampling =
        read_value<int>(input, "frozen oversampling");
    expect_tag(input, "termination");
    const int has_termination = read_value<int>(input, "termination presence");
    if (has_termination == 1) {
        CertifiedTermination termination;
        termination.state = parse_enum<CertifiedDriverState>(
            read_value<std::string>(input, "termination state"));
        termination.code = parse_enum<CertifiedStopCode>(
            read_value<std::string>(input, "termination code"));
        termination.claim = parse_enum<CertifiedEvidenceLevel>(
            read_value<std::string>(input, "termination claim"));
        termination.detail = read_quoted(input, "termination detail");
        checkpoint.termination = std::move(termination);
    } else if (has_termination != 0) {
        throw std::invalid_argument("checkpoint termination presence is invalid");
    }
    expect_tag(input, "mutations");
    const std::uint64_t mutation_size =
        read_value<std::uint64_t>(input, "mutation count");
    if (mutation_size > 100000000ULL)
        throw std::invalid_argument("checkpoint mutation count is unreasonable");
    checkpoint.mutations.reserve(static_cast<std::size_t>(mutation_size));
    for (std::uint64_t index = 0; index < mutation_size; ++index) {
        expect_tag(input, "mutation");
        CertifiedMutation mutation;
        mutation.kind = parse_enum<CertifiedMutationKind>(
            read_value<std::string>(input, "mutation kind"));
        mutation.marked_elements = read_int_vector(input, "mutation marking");
        checkpoint.mutations.push_back(std::move(mutation));
    }
    expect_tag(input, "history");
    const std::uint64_t history_size =
        read_value<std::uint64_t>(input, "history count");
    if (history_size > 100000000ULL)
        throw std::invalid_argument("checkpoint history count is unreasonable");
    checkpoint.history.reserve(static_cast<std::size_t>(history_size));
    for (std::uint64_t index = 0; index < history_size; ++index) {
        expect_tag(input, "record");
        CertifiedTransitionRecord record;
        record.sequence = read_value<std::uint64_t>(input, "record sequence");
        record.before = parse_enum<CertifiedDriverState>(
            read_value<std::string>(input, "record before state"));
        record.after = parse_enum<CertifiedDriverState>(
            read_value<std::string>(input, "record after state"));
        record.action = parse_enum<CertifiedDriverAction>(
            read_value<std::string>(input, "record action"));
        record.claim = parse_enum<CertifiedEvidenceLevel>(
            read_value<std::string>(input, "record claim"));
        record.marked_elements = read_int_vector(input, "record marking");
        record.mu_max = read_optional(input, "record mu max");
        record.q_total_upper = read_optional(input, "record q total");
        record.delta_total_lower = read_optional(input, "record delta total");
        record.delta_h_upper = read_optional(input, "record delta h");
        record.lod_error_lower = read_optional(input, "record LOD lower");
        record.lod_error_upper = read_optional(input, "record LOD upper");
        record.audit_error_lower = read_optional(input, "record audit lower");
        record.audit_error_upper = read_optional(input, "record audit upper");
        record.true_error_lower = read_optional(input, "record true lower");
        record.true_error_upper = read_optional(input, "record true upper");
        record.detail = read_quoted(input, "record detail");
        checkpoint.history.push_back(std::move(record));
    }
    expect_tag(input, "end");
    std::string trailing;
    if (input >> trailing)
        throw std::invalid_argument("checkpoint contains trailing data");
    validate_checkpoint_shape(checkpoint);
    return checkpoint;
}

namespace {

void record_checkpoint_failure(
    CertifiedDriverCheckpoint &checkpoint,
    std::string detail) {
    CertifiedTransitionRecord record;
    record.sequence = checkpoint.transition_count;
    record.before = checkpoint.state;
    record.after = CertifiedDriverState::Failure;
    record.action = CertifiedDriverAction::Fail;
    record.claim = checkpoint.claim;
    record.detail = detail;
    checkpoint.state = CertifiedDriverState::Failure;
    checkpoint.termination = CertifiedTermination{
        CertifiedDriverState::Failure,
        CertifiedStopCode::CheckpointIncompatible,
        checkpoint.claim,
        std::move(detail)};
    checkpoint.history.push_back(std::move(record));
    ++checkpoint.transition_count;
}

void replay_mutation(
    CertifiedDriverBackend &backend,
    const CertifiedMutation &mutation) {
    switch (mutation.kind) {
    case CertifiedMutationKind::RefineCoarse:
        backend.refine_coarse(mutation.marked_elements);
        return;
    case CertifiedMutationKind::RefineCorrectorFine:
        backend.refine_corrector_patches(mutation.marked_elements);
        return;
    case CertifiedMutationKind::IncreaseOversampling:
        backend.increase_oversampling();
        return;
    case CertifiedMutationKind::RefineAudit:
        backend.refine_audit(mutation.marked_elements);
        return;
    }
    throw std::invalid_argument("checkpoint contains an unknown mutation");
}

} // namespace

CertifiedAdaptiveDriver::CertifiedAdaptiveDriver(CertifiedDriverConfig config)
    : config_(std::move(config)) {
    validate_config(config_);
    checkpoint_.config_fingerprint = config_fingerprint(config_);
    if (config_.method == CertifiedMethod::Hlod) {
        checkpoint_.frozen_corrector_space_id =
            config_.hlod_prior_corrector_space_id;
        checkpoint_.frozen_oversampling =
            config_.hlod_prior_oversampling;
    }
}

CertifiedAdaptiveDriver::CertifiedAdaptiveDriver(
    CertifiedDriverConfig config,
    CertifiedDriverCheckpoint checkpoint)
    : config_(std::move(config)), checkpoint_(std::move(checkpoint)) {
    validate_config(config_);
    try {
        validate_checkpoint_shape(checkpoint_);
    } catch (const std::exception &error) {
        // A malformed in-memory checkpoint cannot preserve its old history
        // invariants.  Report it through the constructor exception; serialized
        // checkpoints are validated before reaching this path.
        throw std::invalid_argument(
            std::string("invalid certified checkpoint: ") + error.what());
    }
    if (checkpoint_.config_fingerprint != config_fingerprint(config_)) {
        record_checkpoint_failure(
            checkpoint_, "checkpoint configuration fingerprint mismatch");
        return;
    }
    if (config_.method == CertifiedMethod::Hlod
        && (checkpoint_.frozen_corrector_space_id
                != config_.hlod_prior_corrector_space_id
            || checkpoint_.frozen_oversampling
                != config_.hlod_prior_oversampling)) {
        record_checkpoint_failure(
            checkpoint_, "checkpoint HLOD prior differs from frozen configuration");
    }
}

CertifiedAdaptiveDriver CertifiedAdaptiveDriver::resume(
    CertifiedDriverConfig config,
    const CertifiedDriverCheckpoint &checkpoint,
    CertifiedDriverBackend &fresh_backend) {
    CertifiedAdaptiveDriver driver(std::move(config), checkpoint);
    if (driver.terminal()) return driver;
    try {
        CertifiedWorkSnapshot initial = fresh_backend.work_snapshot();
        if (const auto issue = validate_snapshot(initial)) {
            record_checkpoint_failure(
                driver.checkpoint_, "fresh backend snapshot is invalid: " + *issue);
            return driver;
        }
        if (!driver.checkpoint_.initial_problem_id.empty()
            && initial.initial_problem_id
                != driver.checkpoint_.initial_problem_id) {
            record_checkpoint_failure(
                driver.checkpoint_, "fresh backend initial problem id mismatch");
            return driver;
        }
        if (driver.config_.method == CertifiedMethod::Hlod
            && (initial.corrector_space_id
                    != driver.config_.hlod_prior_corrector_space_id
                || initial.oversampling
                    != driver.config_.hlod_prior_oversampling)) {
            record_checkpoint_failure(
                driver.checkpoint_, "fresh backend violates the frozen HLOD prior");
            return driver;
        }
        for (const CertifiedMutation &mutation : driver.checkpoint_.mutations)
            replay_mutation(fresh_backend, mutation);
        CertifiedWorkSnapshot restored = fresh_backend.work_snapshot();
        if (const auto issue = validate_snapshot(restored)) {
            record_checkpoint_failure(
                driver.checkpoint_, "restored backend snapshot is invalid: " + *issue);
            return driver;
        }
        if (!driver.checkpoint_.initial_problem_id.empty()
            && restored.initial_problem_id
                != driver.checkpoint_.initial_problem_id) {
            record_checkpoint_failure(
                driver.checkpoint_, "restored backend initial problem id mismatch");
            return driver;
        }
        if (!driver.checkpoint_.backend_fingerprint.empty()
            && restored.state_fingerprint
                != driver.checkpoint_.backend_fingerprint) {
            record_checkpoint_failure(
                driver.checkpoint_, "mutation replay produced a different backend fingerprint");
            return driver;
        }
        if (driver.config_.method == CertifiedMethod::Hlod
            && (restored.corrector_space_id
                    != driver.config_.hlod_prior_corrector_space_id
                || restored.oversampling
                    != driver.config_.hlod_prior_oversampling)) {
            record_checkpoint_failure(
                driver.checkpoint_, "mutation replay changed the frozen HLOD prior");
            return driver;
        }
        if (driver.checkpoint_.initial_problem_id.empty())
            driver.checkpoint_.initial_problem_id = restored.initial_problem_id;
        driver.checkpoint_.backend_fingerprint = restored.state_fingerprint;
        // Replay establishes a fresh backend-local accounting epoch.  The
        // replay/setup cost itself is outside the resumed algorithm budget;
        // subsequent deltas continue from the persisted v2 cumulative values.
        // Peak memory remains conservative across both epochs.
        driver.backend_work_epoch_start_ = restored.backend_work_units;
        driver.backend_work_epoch_offset_ =
            driver.checkpoint_.cumulative_backend_work_units;
        driver.peak_memory_floor_ =
            driver.checkpoint_.cumulative_peak_memory_bytes;
        driver.elapsed_epoch_start_ = restored.elapsed_seconds;
        driver.elapsed_epoch_offset_ =
            driver.checkpoint_.cumulative_elapsed_seconds;
        if (const auto issue = driver.account_work_snapshot(restored)) {
            record_checkpoint_failure(
                driver.checkpoint_,
                "restored backend work accounting is invalid: " + *issue);
            return driver;
        }
    } catch (const std::exception &error) {
        record_checkpoint_failure(
            driver.checkpoint_,
            std::string("checkpoint mutation replay failed: ") + error.what());
    } catch (...) {
        record_checkpoint_failure(
            driver.checkpoint_, "checkpoint mutation replay failed with unknown error");
    }
    return driver;
}

CertifiedAdaptiveDriver CertifiedAdaptiveDriver::resume(
    CertifiedDriverConfig config,
    std::string_view serialized_checkpoint,
    CertifiedDriverBackend &fresh_backend) {
    return resume(
        std::move(config),
        deserialize_certified_checkpoint(serialized_checkpoint),
        fresh_backend);
}

bool CertifiedAdaptiveDriver::terminal() const {
    return terminal_state(checkpoint_.state);
}

std::string_view CertifiedAdaptiveDriver::output_namespace() const {
    return config_.method == CertifiedMethod::Calod
        ? std::string_view("helmholtz/calod")
        : std::string_view("helmholtz/hlod");
}

std::optional<std::string> CertifiedAdaptiveDriver::account_work_snapshot(
    CertifiedWorkSnapshot &snapshot) {
    if (snapshot.backend_work_units < backend_work_epoch_start_)
        return "backend work counter decreased within its accounting epoch";
    if (!std::isfinite(snapshot.elapsed_seconds)
        || snapshot.elapsed_seconds < elapsed_epoch_start_) {
        return "backend elapsed time decreased within its accounting epoch";
    }
    const std::uint64_t work_delta =
        snapshot.backend_work_units - backend_work_epoch_start_;
    if (work_delta > std::numeric_limits<std::uint64_t>::max()
                         - backend_work_epoch_offset_) {
        return "cumulative backend work counter overflowed";
    }
    const double elapsed_delta =
        snapshot.elapsed_seconds - elapsed_epoch_start_;
    const double cumulative_elapsed = elapsed_epoch_offset_ + elapsed_delta;
    if (!std::isfinite(cumulative_elapsed))
        return "cumulative backend elapsed time overflowed";

    const std::uint64_t cumulative_work =
        backend_work_epoch_offset_ + work_delta;
    if (cumulative_work < checkpoint_.cumulative_backend_work_units)
        return "cumulative backend work counter decreased";
    if (cumulative_elapsed < checkpoint_.cumulative_elapsed_seconds)
        return "cumulative backend elapsed time decreased";

    snapshot.backend_work_units = cumulative_work;
    snapshot.peak_memory_bytes = std::max(
        peak_memory_floor_, snapshot.peak_memory_bytes);
    snapshot.elapsed_seconds = cumulative_elapsed;
    checkpoint_.cumulative_backend_work_units = snapshot.backend_work_units;
    checkpoint_.cumulative_peak_memory_bytes = std::max(
        checkpoint_.cumulative_peak_memory_bytes,
        snapshot.peak_memory_bytes);
    checkpoint_.cumulative_elapsed_seconds = snapshot.elapsed_seconds;
    return std::nullopt;
}

bool CertifiedAdaptiveDriver::step(CertifiedDriverBackend &backend) {
    if (terminal()) return false;

    CertifiedTransitionRecord record;
    record.sequence = checkpoint_.transition_count;
    record.before = checkpoint_.state;
    record.after = checkpoint_.state;
    record.claim = checkpoint_.claim;

    const auto commit = [&](CertifiedDriverAction action,
                            CertifiedDriverState after,
                            std::string detail) {
        record.action = action;
        record.after = after;
        record.claim = checkpoint_.claim;
        record.detail = std::move(detail);
        checkpoint_.state = after;
        checkpoint_.history.push_back(record);
        ++checkpoint_.transition_count;
        return !terminal_state(after);
    };

    const auto terminate = [&](CertifiedDriverState terminal,
                               CertifiedStopCode code,
                               CertifiedDriverAction action,
                               std::string detail) {
        checkpoint_.termination = CertifiedTermination{
            terminal, code, checkpoint_.claim, detail};
        return commit(action, terminal, std::move(detail));
    };

    if (checkpoint_.transition_count >= config_.limits.max_state_transitions) {
        return terminate(
            CertifiedDriverState::WorkLimit,
            CertifiedStopCode::StateTransitionLimit,
            CertifiedDriverAction::StopAtWorkLimit,
            "maximum state-transition count reached");
    }

    try {
        CertifiedWorkSnapshot snapshot = backend.work_snapshot();
        if (const auto issue = validate_snapshot(snapshot)) {
            return terminate(
                CertifiedDriverState::Failure,
                CertifiedStopCode::InvalidObservation,
                CertifiedDriverAction::Fail,
                "invalid backend work snapshot: " + *issue);
        }
        if (const auto issue = account_work_snapshot(snapshot)) {
            return terminate(
                CertifiedDriverState::Failure,
                CertifiedStopCode::InvalidObservation,
                CertifiedDriverAction::Fail,
                "invalid backend work accounting: " + *issue);
        }
        if (checkpoint_.initial_problem_id.empty())
            checkpoint_.initial_problem_id = snapshot.initial_problem_id;
        else if (checkpoint_.initial_problem_id != snapshot.initial_problem_id) {
            return terminate(
                CertifiedDriverState::Failure,
                CertifiedStopCode::CheckpointIncompatible,
                CertifiedDriverAction::Fail,
                "backend initial problem id changed during the run");
        }
        if (!checkpoint_.backend_fingerprint.empty()
            && checkpoint_.backend_fingerprint != snapshot.state_fingerprint) {
            return terminate(
                CertifiedDriverState::Failure,
                CertifiedStopCode::CheckpointIncompatible,
                CertifiedDriverAction::Fail,
                "backend state changed outside the certified driver journal");
        }
        if (config_.method == CertifiedMethod::Hlod
            && (snapshot.corrector_space_id
                    != checkpoint_.frozen_corrector_space_id
                || snapshot.oversampling
                    != checkpoint_.frozen_oversampling)) {
            return terminate(
                CertifiedDriverState::Failure,
                CertifiedStopCode::FrozenHlodPriorMismatch,
                CertifiedDriverAction::Fail,
                "HLOD corrector fine space or oversampling differs from the frozen prior");
        }
        checkpoint_.backend_fingerprint = snapshot.state_fingerprint;
        if (const auto limit = resource_limit(snapshot, config_.limits)) {
            return terminate(
                CertifiedDriverState::WorkLimit,
                limit->first,
                CertifiedDriverAction::StopAtWorkLimit,
                limit->second);
        }

        const auto evidence_issue = [&](const CertifiedObservationEvidence &evidence,
                                        std::string_view context)
            -> std::optional<std::pair<CertifiedStopCode, std::string>> {
            if (!evidence.valid) {
                return std::pair{
                    CertifiedStopCode::InvalidObservation,
                    std::string(context) + " is invalid: "
                        + (evidence.invalid_reason.empty()
                               ? "backend supplied no reason"
                               : evidence.invalid_reason)};
            }
            if (evidence.source.empty()) {
                return std::pair{
                    CertifiedStopCode::InvalidObservation,
                    std::string(context) + " has no evidence source"};
            }
            if (evidence.level == CertifiedEvidenceLevel::Verified
                && evidence.hash.empty()) {
                return std::pair{
                    CertifiedStopCode::InvalidObservation,
                    std::string(context) + " has no verified evidence hash"};
            }
            if (evidence.level == CertifiedEvidenceLevel::Conditional)
                checkpoint_.claim = CertifiedEvidenceLevel::Conditional;
            if (evidence.level == CertifiedEvidenceLevel::Conditional
                && config_.evidence_policy
                    == CertifiedEvidencePolicy::RequireVerified) {
                return std::pair{
                    CertifiedStopCode::UnverifiedEvidence,
                    std::string(context)
                        + " is conditional while verified evidence is required"};
            }
            return std::nullopt;
        };

        const auto count_limit = [&](CertifiedMutationKind kind,
                                     std::uint64_t maximum,
                                     CertifiedStopCode code,
                                     std::string detail) -> std::optional<bool> {
            if (maximum != 0 && mutation_count(checkpoint_, kind) >= maximum) {
                return terminate(
                    CertifiedDriverState::WorkLimit,
                    code,
                    CertifiedDriverAction::StopAtWorkLimit,
                    std::move(detail));
            }
            return std::nullopt;
        };

        const auto finish_mutation = [&](CertifiedDriverAction action,
                                         CertifiedDriverState next,
                                         std::string detail) {
            CertifiedWorkSnapshot updated = backend.work_snapshot();
            if (const auto issue = validate_snapshot(updated)) {
                return terminate(
                    CertifiedDriverState::Failure,
                    CertifiedStopCode::InvalidObservation,
                    action,
                    "post-mutation backend snapshot is invalid: " + *issue);
            }
            if (const auto issue = account_work_snapshot(updated)) {
                return terminate(
                    CertifiedDriverState::Failure,
                    CertifiedStopCode::InvalidObservation,
                    action,
                    "post-mutation backend work accounting is invalid: "
                        + *issue);
            }
            if (updated.initial_problem_id != checkpoint_.initial_problem_id) {
                return terminate(
                    CertifiedDriverState::Failure,
                    CertifiedStopCode::CheckpointIncompatible,
                    action,
                    "backend initial problem id changed after a mutation");
            }
            if (config_.method == CertifiedMethod::Hlod
                && (updated.corrector_space_id
                        != checkpoint_.frozen_corrector_space_id
                    || updated.oversampling
                        != checkpoint_.frozen_oversampling)) {
                return terminate(
                    CertifiedDriverState::Failure,
                    CertifiedStopCode::FrozenHlodPriorMismatch,
                    action,
                    "a driver mutation changed the frozen HLOD prior");
            }
            checkpoint_.backend_fingerprint = updated.state_fingerprint;
            if (const auto limit = resource_limit(updated, config_.limits)) {
                return terminate(
                    CertifiedDriverState::WorkLimit,
                    limit->first,
                    action,
                    detail + "; " + limit->second);
            }
            return commit(action, next, std::move(detail));
        };

        // Observations may consume substantial work, memory, and elapsed
        // time, but they are not allowed to mutate the numerical state.
        // Re-snapshot after each observation so a backend cannot exceed a
        // configured resource limit and then report a successful terminal
        // state before the next driver step.
        const auto finish_observation = [&](std::string_view context)
            -> std::optional<bool> {
            CertifiedWorkSnapshot updated = backend.work_snapshot();
            if (const auto issue = validate_snapshot(updated)) {
                return terminate(
                    CertifiedDriverState::Failure,
                    CertifiedStopCode::InvalidObservation,
                    CertifiedDriverAction::Fail,
                    std::string(context)
                        + " produced an invalid backend snapshot: " + *issue);
            }
            if (const auto issue = account_work_snapshot(updated)) {
                return terminate(
                    CertifiedDriverState::Failure,
                    CertifiedStopCode::InvalidObservation,
                    CertifiedDriverAction::Fail,
                    std::string(context)
                        + " produced invalid backend work accounting: "
                        + *issue);
            }
            if (updated.initial_problem_id != checkpoint_.initial_problem_id) {
                return terminate(
                    CertifiedDriverState::Failure,
                    CertifiedStopCode::CheckpointIncompatible,
                    CertifiedDriverAction::Fail,
                    std::string(context)
                        + " changed the backend initial problem id");
            }
            if (updated.state_fingerprint != checkpoint_.backend_fingerprint) {
                return terminate(
                    CertifiedDriverState::Failure,
                    CertifiedStopCode::CheckpointIncompatible,
                    CertifiedDriverAction::Fail,
                    std::string(context)
                        + " mutated numerical state outside the driver journal");
            }
            if (config_.method == CertifiedMethod::Hlod
                && (updated.corrector_space_id
                        != checkpoint_.frozen_corrector_space_id
                    || updated.oversampling
                        != checkpoint_.frozen_oversampling)) {
                return terminate(
                    CertifiedDriverState::Failure,
                    CertifiedStopCode::FrozenHlodPriorMismatch,
                    CertifiedDriverAction::Fail,
                    std::string(context) + " changed the frozen HLOD prior");
            }
            if (const auto limit = resource_limit(updated, config_.limits)) {
                return terminate(
                    CertifiedDriverState::WorkLimit,
                    limit->first,
                    CertifiedDriverAction::StopAtWorkLimit,
                    std::string(context) + "; " + limit->second);
            }
            return std::nullopt;
        };

        switch (checkpoint_.state) {
        case CertifiedDriverState::CoarseAdmissibility: {
            const CoarseAdmissibilityObservation observation =
                backend.inspect_coarse_admissibility();
            if (const auto stopped = finish_observation(
                    "coarse admissibility inspection"))
                return *stopped;
            if (const auto issue = evidence_issue(
                    observation.evidence, "coarse admissibility observation")) {
                return terminate(
                    CertifiedDriverState::Failure,
                    issue->first,
                    CertifiedDriverAction::Fail,
                    issue->second);
            }
            if (observation.mu_by_element.empty()
                || !finite_nonnegative(observation.mu_by_element)) {
                return terminate(
                    CertifiedDriverState::Failure,
                    CertifiedStopCode::InvalidObservation,
                    CertifiedDriverAction::Fail,
                    "coarse admissibility values must be a nonempty finite nonnegative vector");
            }
            record.mu_max = *std::max_element(
                observation.mu_by_element.begin(),
                observation.mu_by_element.end());
            std::vector<int> violating;
            for (int element = 0;
                 element < static_cast<int>(observation.mu_by_element.size());
                 ++element) {
                if (observation.mu_by_element[element] > config_.mu0)
                    violating.push_back(element);
            }
            if (violating.empty()) {
                return commit(
                    CertifiedDriverAction::AcceptCoarseAdmissibility,
                    CertifiedDriverState::CorrectorCertification,
                    "all coarse elements satisfy mu_H <= mu0");
            }
            record.marked_elements = violating;
            if (!observation.refinement_available) {
                return terminate(
                    CertifiedDriverState::Failure,
                    CertifiedStopCode::CoarseRefinementUnavailable,
                    CertifiedDriverAction::Fail,
                    "coarse admissibility is violated and coarse refinement is unavailable");
            }
            if (const auto stopped = count_limit(
                    CertifiedMutationKind::RefineCoarse,
                    config_.limits.max_coarse_refinements,
                    CertifiedStopCode::CoarseRefinementLimit,
                    "coarse-refinement count limit reached"))
                return *stopped;
            backend.refine_coarse(violating);
            checkpoint_.mutations.push_back(
                {CertifiedMutationKind::RefineCoarse, violating});
            return finish_mutation(
                CertifiedDriverAction::RefineCoarseAdmissibility,
                CertifiedDriverState::CoarseAdmissibility,
                "refined every element violating coarse admissibility");
        }

        case CertifiedDriverState::CorrectorCertification: {
            const CorrectorCertificationObservation observation =
                backend.inspect_corrector_certification();
            if (const auto stopped = finish_observation(
                    "corrector certificate inspection"))
                return *stopped;
            if (const auto issue = evidence_issue(
                    observation.evidence, "corrector certificate")) {
                return terminate(
                    CertifiedDriverState::Failure,
                    issue->first,
                    CertifiedDriverAction::Fail,
                    issue->second);
            }
            if (!std::isfinite(observation.q_total_upper)
                || observation.q_total_upper < 0.0
                || !std::isfinite(observation.delta_total_lower)
                || observation.delta_total_lower < 0.0
                || !std::isfinite(observation.delta_h_upper)
                || observation.delta_h_upper < 0.0
                || !finite_nonnegative(observation.eta_h_element_squared)) {
                return terminate(
                    CertifiedDriverState::Failure,
                    CertifiedStopCode::InvalidObservation,
                    CertifiedDriverAction::Fail,
                    "corrector certificate contains invalid bounds or indicators");
            }
            record.q_total_upper = observation.q_total_upper;
            record.delta_total_lower = observation.delta_total_lower;
            record.delta_h_upper = observation.delta_h_upper;
            const bool accepted = observation.stability_condition_holds
                && observation.q_total_upper <= config_.q0;
            if (accepted) {
                return commit(
                    CertifiedDriverAction::AcceptCorrectorCertificate,
                    CertifiedDriverState::CoarseErrorControl,
                    "stability condition and q_total <= q0 are both satisfied");
            }
            if (config_.method == CertifiedMethod::Hlod) {
                return terminate(
                    CertifiedDriverState::Failure,
                    CertifiedStopCode::FrozenHlodCorrectorFailure,
                    CertifiedDriverAction::Fail,
                    "frozen HLOD prior does not satisfy corrector certification");
            }
            if (observation.delta_h_upper
                > config_.tau * observation.delta_total_lower) {
                if (!observation.corrector_refinement_available) {
                    return terminate(
                        CertifiedDriverState::Failure,
                        CertifiedStopCode::CorrectorRefinementUnavailable,
                        CertifiedDriverAction::Fail,
                        "corrector-h branch selected but local fine refinement is unavailable");
                }
                const double indicator_sum = std::accumulate(
                    observation.eta_h_element_squared.begin(),
                    observation.eta_h_element_squared.end(), 0.0);
                if (!(indicator_sum > 0.0)) {
                    return terminate(
                        CertifiedDriverState::Failure,
                        CertifiedStopCode::EmptyMarking,
                        CertifiedDriverAction::Fail,
                        "corrector-h branch has no positive eta_h indicator energy");
                }
                std::vector<int> marked = mark_doerfler(
                    observation.eta_h_element_squared, config_.theta_h);
                if (marked.empty()) {
                    return terminate(
                        CertifiedDriverState::Failure,
                        CertifiedStopCode::EmptyMarking,
                        CertifiedDriverAction::Fail,
                        "corrector-h Doerfler marking is empty");
                }
                record.marked_elements = marked;
                if (const auto stopped = count_limit(
                        CertifiedMutationKind::RefineCorrectorFine,
                        config_.limits.max_corrector_refinements,
                        CertifiedStopCode::CorrectorRefinementLimit,
                        "corrector-refinement count limit reached"))
                    return *stopped;
                backend.refine_corrector_patches(marked);
                checkpoint_.mutations.push_back(
                    {CertifiedMutationKind::RefineCorrectorFine, marked});
                return finish_mutation(
                    CertifiedDriverAction::RefineCorrectorFine,
                    CertifiedDriverState::CorrectorCertification,
                    "delta_h upper bound dominates tau times the total lower bound");
            }
            if (!observation.oversampling_increment_available) {
                return terminate(
                    CertifiedDriverState::Failure,
                    CertifiedStopCode::OversamplingUnavailable,
                    CertifiedDriverAction::Fail,
                    "reverse corrector branch selected but oversampling cannot be increased");
            }
            if (config_.limits.max_oversampling >= 0
                && snapshot.oversampling >= config_.limits.max_oversampling) {
                return terminate(
                    CertifiedDriverState::WorkLimit,
                    CertifiedStopCode::OversamplingLimit,
                    CertifiedDriverAction::StopAtWorkLimit,
                    "maximum oversampling value reached");
            }
            if (const auto stopped = count_limit(
                    CertifiedMutationKind::IncreaseOversampling,
                    config_.limits.max_oversampling_increments,
                    CertifiedStopCode::OversamplingLimit,
                    "oversampling-increment count limit reached"))
                return *stopped;
            backend.increase_oversampling();
            checkpoint_.mutations.push_back(
                {CertifiedMutationKind::IncreaseOversampling, {}});
            return finish_mutation(
                CertifiedDriverAction::IncreaseOversampling,
                CertifiedDriverState::CorrectorCertification,
                "reverse branch selected global ell <- ell + 1");
        }

        case CertifiedDriverState::CoarseErrorControl: {
            const CoarseErrorObservation observation =
                backend.solve_and_estimate_coarse_error();
            if (const auto stopped = finish_observation(
                    "coarse solve and residual estimate"))
                return *stopped;
            if (const auto issue = evidence_issue(
                    observation.evidence, "coarse error certificate")) {
                return terminate(
                    CertifiedDriverState::Failure,
                    issue->first,
                    CertifiedDriverAction::Fail,
                    issue->second);
            }
            if (!std::isfinite(observation.eta_H) || observation.eta_H < 0.0
                || !std::isfinite(observation.lod_error_lower)
                || !std::isfinite(observation.lod_error_upper)
                || observation.lod_error_lower < 0.0
                || observation.lod_error_upper < observation.lod_error_lower
                || !finite_nonnegative(observation.eta_H_element_squared)) {
                return terminate(
                    CertifiedDriverState::Failure,
                    CertifiedStopCode::InvalidObservation,
                    CertifiedDriverAction::Fail,
                    "coarse error observation contains an invalid interval or indicators");
            }
            const double allocated = std::accumulate(
                observation.eta_H_element_squared.begin(),
                observation.eta_H_element_squared.end(), 0.0);
            const double expected = observation.eta_H * observation.eta_H;
            if (std::abs(allocated - expected)
                > 1e-8 * std::max({1.0, allocated, expected})) {
                return terminate(
                    CertifiedDriverState::Failure,
                    CertifiedStopCode::InvalidObservation,
                    CertifiedDriverAction::Fail,
                    "eta_H element allocation does not sum to eta_H squared");
            }
            record.lod_error_lower = observation.lod_error_lower;
            record.lod_error_upper = observation.lod_error_upper;
            checkpoint_.lod_error_lower = observation.lod_error_lower;
            checkpoint_.lod_error_upper = observation.lod_error_upper;
            checkpoint_.audit_error_lower.reset();
            checkpoint_.audit_error_upper.reset();
            checkpoint_.true_error_lower.reset();
            checkpoint_.true_error_upper.reset();
            if (config_.error_target == CertifiedErrorTarget::AuditSpace
                && observation.lod_error_upper <= config_.tolerance) {
                return terminate(
                    CertifiedDriverState::Done,
                    CertifiedStopCode::AuditToleranceReached,
                    CertifiedDriverAction::CompleteAuditTolerance,
                    "audit-space upper bound reached the requested tolerance");
            }
            if (!(allocated > 0.0)
                && config_.error_target == CertifiedErrorTarget::AuditSpace) {
                return terminate(
                    CertifiedDriverState::Failure,
                    CertifiedStopCode::EmptyMarking,
                    CertifiedDriverAction::Fail,
                    "coarse error is above tolerance but eta_H marking energy is zero");
            }
            std::vector<int> marked;
            if (allocated > 0.0) {
                marked = mark_doerfler(
                    observation.eta_H_element_squared, config_.theta_H);
                if (marked.empty()) {
                    return terminate(
                        CertifiedDriverState::Failure,
                        CertifiedStopCode::EmptyMarking,
                        CertifiedDriverAction::Fail,
                        "eta_H Doerfler marking is empty");
                }
            }
            checkpoint_.pending_coarse_marking = marked;
            checkpoint_.pending_coarse_refinement_available =
                observation.coarse_refinement_available;
            record.marked_elements = marked;
            return commit(
                CertifiedDriverAction::FormPendingCoarseMarking,
                CertifiedDriverState::AuditControl,
                marked.empty()
                    ? "eta_H is zero; entered continuous audit without a coarse marking"
                    : "formed and postponed the eta_H coarse marking until audit control");
        }

        case CertifiedDriverState::AuditControl: {
            if (checkpoint_.pending_coarse_marking.empty()
                && config_.error_target == CertifiedErrorTarget::AuditSpace) {
                return terminate(
                    CertifiedDriverState::Failure,
                    CertifiedStopCode::InvalidObservation,
                    CertifiedDriverAction::Fail,
                    "audit control has no postponed coarse marking");
            }
            if (config_.error_target == CertifiedErrorTarget::AuditSpace) {
                if (!checkpoint_.pending_coarse_refinement_available) {
                    return terminate(
                        CertifiedDriverState::Failure,
                        CertifiedStopCode::CoarseRefinementUnavailable,
                        CertifiedDriverAction::Fail,
                        "audit-space tolerance is unmet and coarse refinement is unavailable");
                }
                if (const auto stopped = count_limit(
                        CertifiedMutationKind::RefineCoarse,
                        config_.limits.max_coarse_refinements,
                        CertifiedStopCode::CoarseRefinementLimit,
                        "coarse-refinement count limit reached"))
                    return *stopped;
                const std::vector<int> marked = checkpoint_.pending_coarse_marking;
                record.marked_elements = marked;
                backend.refine_coarse(marked);
                checkpoint_.mutations.push_back(
                    {CertifiedMutationKind::RefineCoarse, marked});
                checkpoint_.pending_coarse_marking.clear();
                return finish_mutation(
                    CertifiedDriverAction::ApplyPendingCoarseMarking,
                    CertifiedDriverState::CoarseAdmissibility,
                    "audit-space tolerance is unmet; applied postponed eta_H marking");
            }

            const AuditControlObservation observation =
                backend.inspect_audit_control();
            if (const auto stopped = finish_observation(
                    "continuous audit inspection"))
                return *stopped;
            if (observation.interval_kind == AuditIntervalKind::EmpiricalSaturation
                && observation.evidence.level
                    == CertifiedEvidenceLevel::Verified) {
                return terminate(
                    CertifiedDriverState::Failure,
                    CertifiedStopCode::InvalidObservation,
                    CertifiedDriverAction::Fail,
                    "an empirical saturation interval cannot carry verified evidence");
            }
            if (const auto issue = evidence_issue(
                    observation.evidence, "continuous audit interval")) {
                return terminate(
                    CertifiedDriverState::Failure,
                    issue->first,
                    CertifiedDriverAction::Fail,
                    issue->second);
            }
            if (!std::isfinite(observation.audit_error_lower)
                || !std::isfinite(observation.audit_error_upper)
                || observation.audit_error_lower < 0.0
                || observation.audit_error_upper
                    < observation.audit_error_lower
                || !checkpoint_.lod_error_lower
                || !checkpoint_.lod_error_upper) {
                return terminate(
                    CertifiedDriverState::Failure,
                    CertifiedStopCode::InvalidObservation,
                    CertifiedDriverAction::Fail,
                    "continuous audit observation or stored LOD interval is invalid");
            }
            record.audit_error_lower = observation.audit_error_lower;
            record.audit_error_upper = observation.audit_error_upper;
            checkpoint_.audit_error_lower = observation.audit_error_lower;
            checkpoint_.audit_error_upper = observation.audit_error_upper;
            const double audit_threshold_product =
                config_.rho_aud * *checkpoint_.lod_error_upper;
            const double audit_threshold_lower = audit_threshold_product > 0.0
                ? std::nextafter(
                      audit_threshold_product,
                      -std::numeric_limits<double>::infinity())
                : 0.0;
            if (observation.audit_error_upper > audit_threshold_lower) {
                if (!observation.refinement_available) {
                    return terminate(
                        CertifiedDriverState::Failure,
                        CertifiedStopCode::AuditRefinementUnavailable,
                        CertifiedDriverAction::Fail,
                        "audit fraction is violated and audit refinement is unavailable");
                }
                if (std::any_of(
                        observation.marked_fine_elements.begin(),
                        observation.marked_fine_elements.end(),
                        [](int element) { return element < 0; })) {
                    return terminate(
                        CertifiedDriverState::Failure,
                        CertifiedStopCode::InvalidObservation,
                        CertifiedDriverAction::Fail,
                        "audit marking contains a negative fine-element index");
                }
                std::vector<int> marked = normalized_marking(
                    observation.marked_fine_elements);
                if (marked.empty()) {
                    return terminate(
                        CertifiedDriverState::Failure,
                        CertifiedStopCode::EmptyMarking,
                        CertifiedDriverAction::Fail,
                        "audit fraction is violated but the audit marking is empty");
                }
                record.marked_elements = marked;
                if (const auto stopped = count_limit(
                        CertifiedMutationKind::RefineAudit,
                        config_.limits.max_audit_refinements,
                        CertifiedStopCode::AuditRefinementLimit,
                        "audit-refinement count limit reached"))
                    return *stopped;
                backend.refine_audit(marked);
                checkpoint_.mutations.push_back(
                    {CertifiedMutationKind::RefineAudit, marked});
                checkpoint_.pending_coarse_marking.clear();
                checkpoint_.true_error_lower.reset();
                checkpoint_.true_error_upper.reset();
                return finish_mutation(
                    CertifiedDriverAction::RefineAudit,
                    CertifiedDriverState::CorrectorCertification,
                    "audit error exceeds rho_aud times the LOD upper bound; returned to Step 2");
            }

            const double lower_from_lod = std::nextafter(
                *checkpoint_.lod_error_lower
                    - observation.audit_error_upper,
                -std::numeric_limits<double>::infinity());
            const double lower_from_audit = std::nextafter(
                observation.audit_error_lower
                    - *checkpoint_.lod_error_upper,
                -std::numeric_limits<double>::infinity());
            const double true_lower = std::max({
                0.0, lower_from_lod, lower_from_audit});
            const double upper_sum = *checkpoint_.lod_error_upper
                + observation.audit_error_upper;
            if (!std::isfinite(upper_sum)) {
                return terminate(
                    CertifiedDriverState::Failure,
                    CertifiedStopCode::InvalidObservation,
                    CertifiedDriverAction::Fail,
                    "continuous upper-bound addition overflowed");
            }
            const double true_upper = upper_sum > 0.0
                ? std::nextafter(
                      upper_sum,
                      std::numeric_limits<double>::infinity())
                : 0.0;
            if (!std::isfinite(true_upper)) {
                return terminate(
                    CertifiedDriverState::Failure,
                    CertifiedStopCode::InvalidObservation,
                    CertifiedDriverAction::Fail,
                    "continuous upper-bound outward rounding overflowed");
            }
            checkpoint_.true_error_lower = true_lower;
            checkpoint_.true_error_upper = true_upper;
            record.true_error_lower = true_lower;
            record.true_error_upper = true_upper;
            if (true_upper <= config_.tolerance) {
                return terminate(
                    CertifiedDriverState::Done,
                    CertifiedStopCode::ContinuousToleranceReached,
                    CertifiedDriverAction::CompleteContinuousTolerance,
                    "continuous-error upper bound reached the requested tolerance");
            }
            if (!checkpoint_.pending_coarse_refinement_available) {
                return terminate(
                    CertifiedDriverState::Failure,
                    CertifiedStopCode::CoarseRefinementUnavailable,
                    CertifiedDriverAction::Fail,
                    "continuous tolerance is unmet and coarse refinement is unavailable");
            }
            if (checkpoint_.pending_coarse_marking.empty()) {
                return terminate(
                    CertifiedDriverState::Failure,
                    CertifiedStopCode::EmptyMarking,
                    CertifiedDriverAction::Fail,
                    "continuous tolerance is unmet but eta_H supplied no coarse marking");
            }
            if (const auto stopped = count_limit(
                    CertifiedMutationKind::RefineCoarse,
                    config_.limits.max_coarse_refinements,
                    CertifiedStopCode::CoarseRefinementLimit,
                    "coarse-refinement count limit reached"))
                return *stopped;
            const std::vector<int> marked = checkpoint_.pending_coarse_marking;
            record.marked_elements = marked;
            backend.refine_coarse(marked);
            checkpoint_.mutations.push_back(
                {CertifiedMutationKind::RefineCoarse, marked});
            checkpoint_.pending_coarse_marking.clear();
            return finish_mutation(
                CertifiedDriverAction::ApplyPendingCoarseMarking,
                CertifiedDriverState::CoarseAdmissibility,
                "continuous tolerance is unmet; applied postponed eta_H marking");
        }

        case CertifiedDriverState::Done:
        case CertifiedDriverState::WorkLimit:
        case CertifiedDriverState::Failure:
            return false;
        }
    } catch (const std::exception &error) {
        return terminate(
            CertifiedDriverState::Failure,
            CertifiedStopCode::BackendFailure,
            CertifiedDriverAction::Fail,
            std::string("certified backend operation failed: ") + error.what());
    } catch (...) {
        return terminate(
            CertifiedDriverState::Failure,
            CertifiedStopCode::BackendFailure,
            CertifiedDriverAction::Fail,
            "certified backend operation failed with an unknown error");
    }
    return terminate(
        CertifiedDriverState::Failure,
        CertifiedStopCode::InvalidObservation,
        CertifiedDriverAction::Fail,
        "certified driver reached an unknown state");
}

void CertifiedAdaptiveDriver::run(CertifiedDriverBackend &backend) {
    while (step(backend)) {
    }
}

CoarseAdmissibilityObservation make_coarse_admissibility_observation(
    const TriMesh &coarse_mesh,
    double wavenumber,
    const CertificateConstantRegistry &constants) {
    CoarseAdmissibilityObservation observation;
    observation.evidence.source = "WP4 certificate constant registry:C_app";
    if (!(std::isfinite(wavenumber) && wavenumber > 0.0)) {
        observation.evidence.valid = false;
        observation.evidence.invalid_reason =
            "wavenumber must be positive and finite";
        return observation;
    }
    const CertificateConstant *constant = constants.find("C_app");
    if (constant == nullptr) {
        observation.evidence.valid = false;
        observation.evidence.invalid_reason =
            "C_app is absent from the certificate constant registry";
        return observation;
    }
    if (!(std::isfinite(constant->value) && constant->value >= 0.0)) {
        observation.evidence.valid = false;
        observation.evidence.invalid_reason = "C_app is invalid";
        return observation;
    }
    if (constant->direction != CertificateBoundDirection::Upper
        && constant->direction != CertificateBoundDirection::Exact) {
        observation.evidence.valid = false;
        observation.evidence.invalid_reason =
            "C_app is not an upper or exact bound";
        return observation;
    }
    if (coarse_mesh.elems.empty()) {
        observation.evidence.valid = false;
        observation.evidence.invalid_reason = "coarse mesh has no elements";
        return observation;
    }
    // Even a rigorously bounded C_app does not make the ordinary-double
    // products below directed-rounding enclosures.  This context-free adapter
    // therefore remains diagnostic and can never promote evidence.
    observation.evidence.level = CertifiedEvidenceLevel::Conditional;
    observation.evidence.source = constant->source.empty()
        ? observation.evidence.source
        : constant->source;
    observation.mu_by_element.reserve(coarse_mesh.elems.size());
    for (const Triangle &triangle : coarse_mesh.elems) {
        observation.mu_by_element.push_back(
            constant->value * wavenumber
            * element_diameter(coarse_mesh, triangle));
    }
    return observation;
}

CoarseAdmissibilityObservation make_coarse_admissibility_observation(
    const AdaptiveMeshHierarchy &hierarchy,
    const HelmholtzLodModel &model,
    const HelmholtzOperators &audit_operators,
    const CertificateConstantRegistry &constants) {
    const KernelPatchPolicy policy = audit_kernel_patch_policy(hierarchy);
    const CertificateContextFingerprint context =
        certificate_context_fingerprint(
            hierarchy, model, audit_operators, policy);
    CoarseAdmissibilityObservation observation =
        make_coarse_admissibility_observation(
            hierarchy.coarse_mesh(), model.config().wavenumber, constants);
    if (!observation.evidence.valid) return observation;

    observation.evidence.hash = context.mesh + '|' + context.pde + '|'
        + context.patch_policy + '|' + context.operators;
    observation.evidence.level = CertifiedEvidenceLevel::Conditional;
    const bool bound_to_current_context = constants.has_verified(
        "C_app", CertificateBoundDirection::Upper, context);
    observation.evidence.invalid_reason = bound_to_current_context
        ? "C_app is context-verified, but mu_T uses ordinary floating-point multiplication without directed rounding"
        : "C_app is not verified for the internally recomputed numerical context";
    return observation;
}

namespace {

bool complete_verified_corrector_chain(
    const CorrectorCertificateResult &certificate) {
    return certificate.corrector_status == CorrectorCertificateStatus::Certified
        && certificate.verification_metadata.verified
        && certificate.context_fingerprint.complete()
        && certificate.assembly_evidence.valid_for(
            certificate.context_fingerprint)
        && certificate.matrix_enclosure_arithmetic_verified
        && certificate.scalar_formula_enclosures_verified
        && certificate.total_spectrum.verified_lambda.metadata.verified
        && certificate.fine_spectrum.verified_lambda.metadata.verified
        && certificate.audit_infsup.metadata.verified
        && certificate.constants.missing_or_unverified_required(
               certificate.context_fingerprint).empty()
        && certificate.conjugation_passed
        && certificate.stability_verified;
}

} // namespace

CorrectorCertificationObservation make_corrector_certification_observation(
    const CorrectorCertificateResult &certificate) {
    CorrectorCertificationObservation observation;
    observation.evidence.valid =
        certificate.corrector_status != CorrectorCertificateStatus::Invalid;
    const bool verified_chain = complete_verified_corrector_chain(certificate);
    observation.evidence.level =
        verified_chain
        ? CertifiedEvidenceLevel::Verified
        : CertifiedEvidenceLevel::Conditional;
    observation.evidence.source = "WP4 build_corrector_certificates";
    if (!certificate.verification_metadata.backend.empty())
        observation.evidence.source += ":"
            + certificate.verification_metadata.backend;
    observation.evidence.hash = !certificate.assembly_evidence.hash.empty()
        ? certificate.assembly_evidence.hash
        : certificate.patch_policy.hash;
    if (!observation.evidence.valid) {
        observation.evidence.invalid_reason =
            certificate.verification_metadata.failure_reason.empty()
            ? "WP4 corrector certificate is invalid"
            : certificate.verification_metadata.failure_reason;
    } else if (certificate.corrector_status
                   == CorrectorCertificateStatus::Conditional
               && !certificate.corrector_conditional_reasons.empty()) {
        observation.evidence.invalid_reason =
            join_reasons(certificate.corrector_conditional_reasons);
    } else if (!verified_chain
               && certificate.corrector_status
                    == CorrectorCertificateStatus::Certified) {
        observation.evidence.invalid_reason =
            "WP4 corrector status is Certified but one or more verified-chain gates are absent";
    }
    // Whether the computed stability inequality holds is distinct from the
    // evidence level used to justify it.  AllowConditional must be able to
    // consume a finite nonnegative diagnostic margin without promoting the
    // resulting run to Verified.
    observation.stability_condition_holds = observation.evidence.valid
        && certificate.conjugation_passed
        && std::isfinite(certificate.stability_margin)
        && certificate.stability_margin >= 0.0;
    observation.q_total_upper = certificate.q_total;
    observation.delta_total_lower = certificate.delta_total_lower;
    observation.delta_h_upper = certificate.delta_h_upper;
    observation.eta_h_element_squared = certificate.eta_h_element_squared;
    return observation;
}

CoarseErrorObservation make_coarse_error_observation(
    const AuditKernelResidualEstimate &estimate,
    const CorrectorCertificateResult &certificate) {
    CoarseErrorObservation observation;
    observation.evidence.valid =
        certificate.status != CorrectorCertificateStatus::Invalid;
    const bool verified_chain = complete_verified_corrector_chain(certificate)
        && certificate.status == CorrectorCertificateStatus::Certified
        && certificate.verification_metadata.verified
        && estimate.evidence.verified()
        && certificate.eta_H_evidence.verified();
    observation.evidence.level = verified_chain
        ? CertifiedEvidenceLevel::Verified
        : CertifiedEvidenceLevel::Conditional;
    observation.evidence.source =
        "WP3 audit kernel residual + WP4 error enclosure";
    if (!estimate.evidence.backend().empty())
        observation.evidence.source += ":" + estimate.evidence.backend();
    observation.evidence.hash = !certificate.assembly_evidence.hash.empty()
        ? certificate.assembly_evidence.hash
        : certificate.patch_policy.hash;
    if (!estimate.evidence.diagnostic_fingerprint().empty())
        observation.evidence.hash += '|'
            + estimate.evidence.diagnostic_fingerprint();
    observation.eta_H = estimate.eta;
    observation.lod_error_lower = certificate.lod_error_lower;
    observation.lod_error_upper = certificate.lod_error_upper;
    observation.eta_H_element_squared = estimate.element_eta_squared;
    const double scale = std::max({1.0, std::abs(estimate.eta),
                                   std::abs(certificate.eta_H)});
    const bool evidence_matches =
        !estimate.evidence.context_fingerprint().empty()
        && !estimate.evidence.diagnostic_fingerprint().empty()
        && estimate.evidence.matches_result(estimate)
        && certificate.eta_H_evidence.matches_result(estimate)
        && certificate.eta_H_evidence.matches_eta(certificate.eta_H)
        && estimate.evidence.context_fingerprint()
            == certificate.eta_H_evidence.context_fingerprint()
        && estimate.evidence.diagnostic_fingerprint()
            == certificate.eta_H_evidence.diagnostic_fingerprint();
    if (!evidence_matches) {
        observation.evidence.valid = false;
        observation.evidence.invalid_reason =
            "WP3 eta_H evidence token does not match the WP4 enclosure";
    } else if (std::abs(estimate.eta - certificate.eta_H) > 1e-10 * scale) {
        observation.evidence.valid = false;
        observation.evidence.invalid_reason =
            "WP3 eta_H does not match the eta_H used to build the WP4 enclosure";
    } else if (!observation.evidence.valid) {
        observation.evidence.invalid_reason =
            certificate.verification_metadata.failure_reason.empty()
            ? "WP4 coarse error certificate is invalid"
            : certificate.verification_metadata.failure_reason;
    } else if (certificate.status == CorrectorCertificateStatus::Conditional
               && !certificate.conditional_reasons.empty()) {
        observation.evidence.invalid_reason =
            join_reasons(certificate.conditional_reasons);
    }
    return observation;
}

} // namespace lod2d::helmholtz::adaptive

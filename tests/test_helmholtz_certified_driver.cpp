#include "helmholtz/adaptive/certified_driver.h"
#include "helmholtz/adaptive/driver.h"
#include "helmholtz/adaptive/error_control.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using namespace lod2d;
using namespace lod2d::helmholtz;
using namespace lod2d::helmholtz::adaptive;

namespace {

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

void require_close(double first, double second, double tolerance,
                   const char *message) {
    const double scale = std::max({1.0, std::abs(first), std::abs(second)});
    if (std::abs(first - second) > tolerance * scale)
        throw std::runtime_error(message);
}

template <class T>
concept HasEvaluationReferenceAccessor = requires(T &value) {
    value.evaluation_reference();
};

static_assert(!HasEvaluationReferenceAccessor<CertifiedDriverBackend>);
static_assert(!std::is_constructible_v<
              CoarseErrorObservation, EvaluationReferenceError>);
static_assert(!std::is_constructible_v<
              AuditControlObservation, ErrorReference>);

CertifiedObservationEvidence evidence(bool conditional = false) {
    CertifiedObservationEvidence result;
    result.level = conditional
        ? CertifiedEvidenceLevel::Conditional
        : CertifiedEvidenceLevel::Verified;
    result.source = conditional ? "scripted conditional evidence"
                                : "scripted verified evidence";
    result.hash = conditional ? "" : "scripted-evidence-v1";
    return result;
}

enum class ObservationWorkStage {
    None,
    CoarseSolve,
    Audit
};

class ScriptedBackend final : public CertifiedDriverBackend {
public:
    explicit ScriptedBackend(bool always_stable = false,
                             bool freeze_corrector = false,
                             bool conditional = false,
                             bool zero_indicator = false,
                             ObservationWorkStage work_stage =
                                 ObservationWorkStage::None)
        : always_stable_(always_stable),
          freeze_corrector_(freeze_corrector),
          conditional_(conditional),
          zero_indicator_(zero_indicator),
          work_stage_(work_stage) {}

    CertifiedWorkSnapshot work_snapshot() const override {
        CertifiedWorkSnapshot snapshot;
        snapshot.coarse_dofs = 10 + 5 * coarse_refinements_;
        snapshot.fine_dofs = 40 + 10 * corrector_refinements_;
        snapshot.audit_dofs = 80 + 20 * audit_refinements_;
        snapshot.backend_work_units = static_cast<std::uint64_t>(
            coarse_refinements_ + corrector_refinements_
            + audit_refinements_ + (ell_ - 1)
            + observation_work_units_);
        if (forced_work_units_) snapshot.backend_work_units = *forced_work_units_;
        snapshot.initial_problem_id = "scripted-problem-v1";
        snapshot.state_fingerprint = fingerprint();
        snapshot.corrector_space_id = freeze_corrector_
            ? "prior-fine-v1"
            : "adaptive-fine-" + std::to_string(corrector_refinements_);
        snapshot.oversampling = ell_;
        return snapshot;
    }

    CoarseAdmissibilityObservation inspect_coarse_admissibility() override {
        CoarseAdmissibilityObservation result;
        result.evidence = evidence(conditional_);
        result.mu_by_element = coarse_refinements_ == 0
            ? std::vector<double>{0.8, 0.2, 0.9}
            : std::vector<double>{0.2, 0.3, 0.4};
        return result;
    }

    void refine_coarse(const std::vector<int> &marked_elements) override {
        require(!marked_elements.empty(), "scripted coarse marking is empty");
        ++coarse_refinements_;
        mutation_log_.push_back("H:" + encode(marked_elements));
    }

    CorrectorCertificationObservation
    inspect_corrector_certification() override {
        CorrectorCertificationObservation result;
        result.evidence = evidence(conditional_);
        result.q_total_upper = 0.1;
        result.delta_total_lower = 1.0;
        result.eta_h_element_squared = {9.0, 1.0, 0.0};
        if (always_stable_ || (corrector_refinements_ > 0 && ell_ > 1)) {
            result.stability_condition_holds = true;
            result.delta_h_upper = 0.1;
        } else if (corrector_refinements_ == 0) {
            result.delta_h_upper = 0.75;
        } else {
            result.delta_h_upper = 0.25;
        }
        return result;
    }

    void refine_corrector_patches(
        const std::vector<int> &marked_coarse_sources) override {
        require(!marked_coarse_sources.empty(),
                "scripted corrector marking is empty");
        ++corrector_refinements_;
        mutation_log_.push_back("h:" + encode(marked_coarse_sources));
    }

    void increase_oversampling() override {
        ++ell_;
        mutation_log_.push_back("ell:" + std::to_string(ell_));
    }

    CoarseErrorObservation solve_and_estimate_coarse_error() override {
        if (work_stage_ == ObservationWorkStage::CoarseSolve)
            ++observation_work_units_;
        CoarseErrorObservation result;
        result.evidence = evidence(conditional_);
        if (zero_indicator_) {
            result.eta_H = 0.0;
            result.lod_error_lower = 0.0;
            result.lod_error_upper = 0.0;
            result.eta_H_element_squared = {0.0, 0.0, 0.0};
        } else if (coarse_refinements_ < 2) {
            result.eta_H = 1.0;
            result.lod_error_lower = 0.2;
            result.lod_error_upper = 0.4;
            result.eta_H_element_squared = {0.64, 0.36, 0.0};
        } else {
            result.eta_H = 0.5;
            result.lod_error_lower = 0.07;
            result.lod_error_upper = 0.15;
            result.eta_H_element_squared = {0.16, 0.09, 0.0};
        }
        return result;
    }

    AuditControlObservation inspect_audit_control() override {
        if (work_stage_ == ObservationWorkStage::Audit)
            ++observation_work_units_;
        AuditControlObservation result;
        result.evidence = evidence(conditional_);
        result.audit_error_lower = 0.0;
        if (zero_indicator_) {
            result.audit_error_upper = 0.0;
        } else if (audit_refinements_ == 0) {
            result.audit_error_upper = 0.1;
            result.marked_fine_elements = {2, 0, 2};
        } else if (coarse_refinements_ < 2) {
            // The first refined value lies exactly on the unrounded audit
            // threshold.  Directed downward rounding must reject equality;
            // the next refinement is strictly below it.
            result.audit_error_upper = audit_refinements_ == 1
                ? 0.05 * 0.4 : 0.01;
            if (audit_refinements_ == 1)
                result.marked_fine_elements = {1};
        } else {
            result.audit_error_upper = 0.005;
        }
        return result;
    }

    void refine_audit(const std::vector<int> &marked_fine_elements) override {
        require(!marked_fine_elements.empty(),
                "scripted audit marking is empty");
        ++audit_refinements_;
        mutation_log_.push_back("audit:" + encode(marked_fine_elements));
    }

    int coarse_refinements() const { return coarse_refinements_; }
    int corrector_refinements() const { return corrector_refinements_; }
    int audit_refinements() const { return audit_refinements_; }
    int ell() const { return ell_; }
    void force_work_units(std::optional<std::uint64_t> value) {
        forced_work_units_ = value;
    }
    const std::vector<std::string> &mutation_log() const {
        return mutation_log_;
    }

private:
    static std::string encode(const std::vector<int> &values) {
        std::string result;
        for (int value : values) {
            if (!result.empty()) result += ',';
            result += std::to_string(value);
        }
        return result;
    }

    std::string fingerprint() const {
        return "H=" + std::to_string(coarse_refinements_)
            + ";h=" + std::to_string(corrector_refinements_)
            + ";ell=" + std::to_string(ell_)
            + ";audit=" + std::to_string(audit_refinements_);
    }

    bool always_stable_ = false;
    bool freeze_corrector_ = false;
    bool conditional_ = false;
    bool zero_indicator_ = false;
    ObservationWorkStage work_stage_ = ObservationWorkStage::None;
    int coarse_refinements_ = 0;
    int corrector_refinements_ = 0;
    int audit_refinements_ = 0;
    int observation_work_units_ = 0;
    int ell_ = 1;
    std::optional<std::uint64_t> forced_work_units_;
    std::vector<std::string> mutation_log_;
};

CertifiedDriverConfig calod_config() {
    CertifiedDriverConfig config;
    config.method = CertifiedMethod::Calod;
    config.error_target = CertifiedErrorTarget::Continuous;
    config.mu0 = 0.5;
    config.q0 = 0.25;
    config.tau = 0.5;
    config.theta_H = 0.5;
    config.theta_h = 0.5;
    config.rho_aud = 0.05;
    config.tolerance = 0.2;
    config.limits.max_state_transitions = 100;
    return config;
}

std::vector<CertifiedDriverAction> actions(
    const CertifiedDriverCheckpoint &checkpoint) {
    std::vector<CertifiedDriverAction> result;
    for (const CertifiedTransitionRecord &record : checkpoint.history)
        result.push_back(record.action);
    return result;
}

void verify_full_calod_order_and_branches() {
    ScriptedBackend backend;
    CertifiedAdaptiveDriver driver(calod_config());
    require(driver.output_namespace() == "helmholtz/calod",
            "CALOD output namespace is wrong");
    driver.run(backend);
    require(driver.terminal(), "CALOD scripted run did not terminate");
    if (driver.checkpoint().state != CertifiedDriverState::Done) {
        throw std::runtime_error(
            "CALOD scripted run did not succeed: "
            + driver.termination()->detail);
    }
    require(driver.termination()->code
                == CertifiedStopCode::ContinuousToleranceReached,
            "CALOD used the wrong success stop code");
    require(driver.termination()->claim == CertifiedEvidenceLevel::Verified,
            "verified CALOD run lost its verified claim");
    require(backend.coarse_refinements() == 2,
            "CALOD did not execute admissibility-H and eta_H refinements");
    require(backend.corrector_refinements() == 1,
            "CALOD did not execute exactly one corrector-h refinement");
    require(backend.ell() == 2,
            "CALOD did not execute the global reverse ell branch");
    require(backend.audit_refinements() == 2,
            "CALOD did not conservatively refine the threshold-equality audit");

    const std::vector<CertifiedDriverAction> expected{
        CertifiedDriverAction::RefineCoarseAdmissibility,
        CertifiedDriverAction::AcceptCoarseAdmissibility,
        CertifiedDriverAction::RefineCorrectorFine,
        CertifiedDriverAction::IncreaseOversampling,
        CertifiedDriverAction::AcceptCorrectorCertificate,
        CertifiedDriverAction::FormPendingCoarseMarking,
        CertifiedDriverAction::RefineAudit,
        CertifiedDriverAction::AcceptCorrectorCertificate,
        CertifiedDriverAction::FormPendingCoarseMarking,
        CertifiedDriverAction::RefineAudit,
        CertifiedDriverAction::AcceptCorrectorCertificate,
        CertifiedDriverAction::FormPendingCoarseMarking,
        CertifiedDriverAction::ApplyPendingCoarseMarking,
        CertifiedDriverAction::AcceptCoarseAdmissibility,
        CertifiedDriverAction::AcceptCorrectorCertificate,
        CertifiedDriverAction::FormPendingCoarseMarking,
        CertifiedDriverAction::CompleteContinuousTolerance};
    require(actions(driver.checkpoint()) == expected,
            "CALOD transition order differs from the paper algorithm");
    require_close(*driver.checkpoint().true_error_upper, 0.155, 1e-14,
                  "continuous upper interval was formed incorrectly");
    require_close(*driver.checkpoint().true_error_lower, 0.065, 1e-14,
                  "continuous lower interval was formed incorrectly");
    require(*driver.checkpoint().true_error_upper > 0.155,
            "asymmetric continuous upper endpoint was not rounded outward");
    require(*driver.checkpoint().true_error_lower < 0.065,
            "asymmetric continuous lower endpoint was not rounded outward");
}

void verify_checkpoint_resume_is_identical() {
    const CertifiedDriverConfig config = calod_config();
    ScriptedBackend continuous_backend;
    CertifiedAdaptiveDriver continuous(config);
    continuous.run(continuous_backend);
    const std::string continuous_checkpoint =
        serialize_certified_checkpoint(continuous.checkpoint());

    ScriptedBackend interrupted_backend;
    CertifiedAdaptiveDriver interrupted(config);
    for (int step = 0; step < 7; ++step)
        require(interrupted.step(interrupted_backend),
                "scripted run terminated before checkpoint");
    require(interrupted.checkpoint().state
                == CertifiedDriverState::CorrectorCertification,
            "checkpoint was not taken after the audit branch returned to Step 2");
    const std::string serialized =
        serialize_certified_checkpoint(interrupted.checkpoint());
    require(interrupted.checkpoint().cumulative_backend_work_units > 0,
            "checkpoint did not persist cumulative backend work");
    const CertifiedDriverCheckpoint parsed =
        deserialize_certified_checkpoint(serialized);
    require(serialize_certified_checkpoint(parsed) == serialized,
            "checkpoint serialization is not canonical on round trip");

    ScriptedBackend resumed_backend;
    CertifiedAdaptiveDriver resumed = CertifiedAdaptiveDriver::resume(
        config, serialized, resumed_backend);
    require(!resumed.terminal(), "valid checkpoint replay failed");
    require(resumed.checkpoint().cumulative_backend_work_units
                == interrupted.checkpoint().cumulative_backend_work_units,
            "checkpoint replay reset or double-counted backend work");
    resumed.run(resumed_backend);
    require(serialize_certified_checkpoint(resumed.checkpoint())
                == continuous_checkpoint,
            "resumed and continuous state histories differ");
    require(resumed_backend.mutation_log()
                == continuous_backend.mutation_log(),
            "checkpoint mutation replay differs from continuous execution");
    require(resumed.checkpoint().cumulative_backend_work_units
                == continuous.checkpoint().cumulative_backend_work_units,
            "resumed backend work differs from uninterrupted execution");

    std::string legacy = serialized;
    const std::string current_header =
        "LOD2D_CERTIFIED_DRIVER_CHECKPOINT 2";
    const std::size_t header = legacy.find(current_header);
    require(header != std::string::npos,
            "v2 checkpoint header is missing");
    legacy.replace(header, current_header.size(),
                   "LOD2D_CERTIFIED_DRIVER_CHECKPOINT 1");
    bool rejected_v1 = false;
    try {
        (void)deserialize_certified_checkpoint(legacy);
    } catch (const std::invalid_argument &) {
        rejected_v1 = true;
    }
    require(rejected_v1,
            "v1 checkpoint without resource accounting was accepted");

    const auto replace_intervals = [&](const std::string &line) {
        std::string tampered = serialized;
        const std::size_t begin = tampered.find("intervals");
        require(begin != std::string::npos,
                "checkpoint intervals line is missing");
        const std::size_t end = tampered.find('\n', begin);
        require(end != std::string::npos,
                "checkpoint intervals line is unterminated");
        tampered.replace(begin, end - begin, line);
        return tampered;
    };
    for (const std::string &tampered : {
             replace_intervals("intervals 1 2 1 1 0 0 0 0"),
             replace_intervals("intervals 1 0.2 0 0 0 0 0")}) {
        bool rejected_interval = false;
        try {
            (void)deserialize_certified_checkpoint(tampered);
        } catch (const std::invalid_argument &) {
            rejected_interval = true;
        }
        require(rejected_interval,
                "malformed checkpoint interval was accepted");
    }

    CertifiedDriverConfig incompatible = config;
    incompatible.tolerance = 0.19;
    ScriptedBackend incompatible_backend;
    CertifiedAdaptiveDriver rejected = CertifiedAdaptiveDriver::resume(
        incompatible, serialized, incompatible_backend);
    require(rejected.terminal()
                && rejected.checkpoint().state == CertifiedDriverState::Failure,
            "incompatible checkpoint was not rejected structurally");
    require(rejected.termination()->code
                == CertifiedStopCode::CheckpointIncompatible,
            "incompatible checkpoint used the wrong stop code");
}

void verify_hlod_freeze_and_shared_eta_marking() {
    CertifiedDriverConfig config = calod_config();
    config.method = CertifiedMethod::Hlod;
    config.error_target = CertifiedErrorTarget::AuditSpace;
    config.hlod_prior_corrector_space_id = "prior-fine-v1";
    config.hlod_prior_oversampling = 1;
    ScriptedBackend backend(/*always_stable=*/true,
                            /*freeze_corrector=*/true);
    CertifiedAdaptiveDriver driver(config);
    require(driver.output_namespace() == "helmholtz/hlod",
            "HLOD output namespace is wrong");
    driver.run(backend);
    require(driver.checkpoint().state == CertifiedDriverState::Done,
            "frozen HLOD baseline did not complete");
    require(driver.termination()->code
                == CertifiedStopCode::AuditToleranceReached,
            "HLOD used the wrong audit-space success code");
    require(backend.corrector_refinements() == 0 && backend.ell() == 1,
            "HLOD changed its frozen h or ell prior");
    require(backend.coarse_refinements() == 2,
            "HLOD did not reuse admissibility and eta_H coarse marking");
    const auto run_actions = actions(driver.checkpoint());
    require(std::find(run_actions.begin(), run_actions.end(),
                      CertifiedDriverAction::RefineCorrectorFine)
                == run_actions.end(),
            "HLOD entered the corrector-h branch");
    require(std::find(run_actions.begin(), run_actions.end(),
                      CertifiedDriverAction::IncreaseOversampling)
                == run_actions.end(),
            "HLOD entered the ell branch");

    ScriptedBackend unstable_backend(/*always_stable=*/false,
                                     /*freeze_corrector=*/true);
    CertifiedAdaptiveDriver unstable(config);
    unstable.run(unstable_backend);
    require(unstable.checkpoint().state == CertifiedDriverState::Failure,
            "unstable frozen HLOD silently continued");
    require(unstable.termination()->code
                == CertifiedStopCode::FrozenHlodCorrectorFailure,
            "unstable frozen HLOD used the wrong failure code");
    require(unstable_backend.corrector_refinements() == 0
                && unstable_backend.ell() == 1,
            "failed HLOD mutated frozen corrector parameters");
}

void verify_conditional_and_work_limit_outcomes() {
    ScriptedBackend conditional_backend(
        /*always_stable=*/false,
        /*freeze_corrector=*/false,
        /*conditional=*/true);
    CertifiedAdaptiveDriver strict(calod_config());
    strict.run(conditional_backend);
    require(strict.checkpoint().state == CertifiedDriverState::Failure,
            "strict CALOD accepted conditional evidence");
    require(strict.termination()->code
                == CertifiedStopCode::UnverifiedEvidence,
            "strict CALOD used the wrong conditional-evidence failure code");
    require(strict.termination()->claim
                == CertifiedEvidenceLevel::Conditional,
            "strict CALOD failure retained a verified claim after observing conditional evidence");

    CertifiedDriverConfig conditional_config = calod_config();
    conditional_config.evidence_policy =
        CertifiedEvidencePolicy::AllowConditional;
    ScriptedBackend allowed_backend(
        /*always_stable=*/false,
        /*freeze_corrector=*/false,
        /*conditional=*/true);
    CertifiedAdaptiveDriver allowed(conditional_config);
    allowed.run(allowed_backend);
    require(allowed.checkpoint().state == CertifiedDriverState::Done,
            "explicit conditional CALOD mode did not complete");
    require(allowed.termination()->claim
                == CertifiedEvidenceLevel::Conditional,
            "conditional CALOD result was promoted to verified");

    CertifiedDriverConfig limited_config = calod_config();
    limited_config.limits.max_state_transitions = 2;
    ScriptedBackend limited_backend;
    CertifiedAdaptiveDriver limited(limited_config);
    limited.run(limited_backend);
    require(limited.checkpoint().state == CertifiedDriverState::WorkLimit,
            "state transition limit did not produce WORK_LIMIT");
    require(limited.termination()->code
                == CertifiedStopCode::StateTransitionLimit,
            "state transition limit used the wrong stop code");
    const std::string serialized =
        serialize_certified_checkpoint(limited.checkpoint());
    const CertifiedDriverCheckpoint restored =
        deserialize_certified_checkpoint(serialized);
    require(restored.termination.has_value()
                && restored.termination->code
                    == CertifiedStopCode::StateTransitionLimit,
            "structured work-limit reason did not survive serialization");
}

void verify_backend_work_counter_cannot_regress() {
    ScriptedBackend backend;
    CertifiedAdaptiveDriver driver(calod_config());
    require(driver.step(backend),
            "scripted run stopped before work-regression fixture");
    require(driver.checkpoint().cumulative_backend_work_units == 1,
            "work-regression fixture has the wrong initial counter");
    backend.force_work_units(0);
    driver.run(backend);
    require(driver.checkpoint().state == CertifiedDriverState::Failure,
            "decreasing backend work counter was accepted");
    require(driver.termination()->code
                == CertifiedStopCode::InvalidObservation,
            "decreasing backend work counter used the wrong failure code");
}

void verify_post_observation_limits_and_zero_continuous_case() {
    CertifiedDriverConfig zero_config = calod_config();
    ScriptedBackend zero_backend(
        /*always_stable=*/true,
        /*freeze_corrector=*/false,
        /*conditional=*/false,
        /*zero_indicator=*/true);
    CertifiedAdaptiveDriver zero_driver(zero_config);
    zero_driver.run(zero_backend);
    require(zero_driver.checkpoint().state == CertifiedDriverState::Done,
            "continuous zero-indicator case did not reach audit completion");
    require(zero_driver.termination()->code
                == CertifiedStopCode::ContinuousToleranceReached,
            "continuous zero-indicator case used the wrong stop code");

    for (ObservationWorkStage stage : {
             ObservationWorkStage::CoarseSolve,
             ObservationWorkStage::Audit}) {
        CertifiedDriverConfig limited = calod_config();
        limited.limits.max_backend_work_units = 1;
        ScriptedBackend backend(
            /*always_stable=*/true,
            /*freeze_corrector=*/false,
            /*conditional=*/false,
            /*zero_indicator=*/stage == ObservationWorkStage::Audit,
            stage);
        CertifiedAdaptiveDriver driver(limited);
        driver.run(backend);
        require(driver.checkpoint().state == CertifiedDriverState::WorkLimit,
                "post-observation backend work was not limited");
        require(driver.termination()->code
                    == CertifiedStopCode::BackendWorkLimit,
                "post-observation work used the wrong stop code");
    }
}

void verify_conditional_wp4_adapter_preserves_stability_condition() {
    CorrectorCertificateResult certificate;
    certificate.status = CorrectorCertificateStatus::Conditional;
    certificate.conditional_reasons = {"diagnostic floating-point evidence"};
    certificate.conjugation_passed = true;
    certificate.stability_margin = 0.125;
    certificate.q_total = 0.1;
    certificate.delta_total_lower = 1.0;
    certificate.delta_h_upper = 0.1;
    certificate.eta_h_element_squared = {1.0};
    const CorrectorCertificationObservation observation =
        make_corrector_certification_observation(certificate);
    require(observation.evidence.valid
                && observation.evidence.level
                    == CertifiedEvidenceLevel::Conditional,
            "conditional WP4 certificate used the wrong evidence level");
    require(observation.stability_condition_holds,
            "conditional WP4 adapter discarded a satisfied stability inequality");

    certificate.conjugation_passed = false;
    const CorrectorCertificationObservation missing_adjoint =
        make_corrector_certification_observation(certificate);
    require(!missing_adjoint.stability_condition_holds,
            "failed conjugation without an independent adjoint was accepted");

    CorrectorCertificateResult forged;
    forged.corrector_status = CorrectorCertificateStatus::Certified;
    forged.status = CorrectorCertificateStatus::Certified;
    forged.conjugation_passed = true;
    forged.stability_verified = true;
    forged.stability_margin = 0.125;
    const CorrectorCertificationObservation rejected_forgery =
        make_corrector_certification_observation(forged);
    require(rejected_forgery.evidence.level
                == CertifiedEvidenceLevel::Conditional,
            "public Certified status bypassed the WP4 verified-chain gates");
}

void verify_wp4_coarse_adapter_and_proxy_identity() {
    TriMesh mesh;
    mesh.nodes = {Point2(0.0, 0.0), Point2(1.0, 0.0), Point2(0.0, 1.0)};
    mesh.elems = {Triangle{0, 1, 2}};
    CertificateConstantRegistry constants;
    constants.set({"C_app", 0.25, CertificateBoundDirection::Upper,
                   "verified adapter test", "direct input", "unit triangle",
                   "adapter-test-v1", true,
                   "adapter-mesh-v1", "adapter-pde-v1",
                   "adapter-operators-v1"});
    const CoarseAdmissibilityObservation observation =
        make_coarse_admissibility_observation(mesh, 2.0, constants);
    require(observation.evidence.valid
                && observation.evidence.level
                    == CertifiedEvidenceLevel::Conditional,
            "ordinary-double mu_T was promoted to verified evidence");
    require_close(observation.mu_by_element.at(0), std::sqrt(2.0) / 2.0,
                  1e-14, "coarse admissibility adapter computed the wrong mu_T");

    AdaptiveHelmholtzResult proxy;
    require(proxy.output_namespace == "helmholtz/hlod_proxy",
            "legacy H-only result is not isolated in the proxy namespace");
    require(proxy.implementation_status == "diagnostic_h_only_proxy",
            "legacy H-only result lost its diagnostic proxy status");
}

} // namespace

int main() {
    try {
        verify_full_calod_order_and_branches();
        verify_checkpoint_resume_is_identical();
        verify_backend_work_counter_cannot_regress();
        verify_hlod_freeze_and_shared_eta_marking();
        verify_conditional_and_work_limit_outcomes();
        verify_post_observation_limits_and_zero_continuous_case();
        verify_conditional_wp4_adapter_preserves_stability_condition();
        verify_wp4_coarse_adapter_and_proxy_identity();
        std::cout << "Certified CALOD/HLOD states, branches, isolation, and replay passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_certified_driver failed: "
                  << error.what() << '\n';
        return 1;
    }
}

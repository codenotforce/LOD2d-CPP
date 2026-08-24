#include "helmholtz/adaptive/numerical_backend.h"
#include "helmholtz/benchmarks/paper_cases.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

using namespace lod2d;
using namespace lod2d::helmholtz;
using namespace lod2d::helmholtz::adaptive;
using namespace lod2d::helmholtz::benchmarks;
using namespace lod2d::helmholtz::experiments;

namespace {

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

template <class T>
concept HasExactMember = requires(T value) { value.exact; };

template <class T>
concept HasReferenceMember = requires(T value) { value.reference; };

static_assert(!HasExactMember<NumericalCertifiedBackendConfig>);
static_assert(!HasReferenceMember<NumericalCertifiedBackendConfig>);

CertificateConstantRegistry conditional_constants() {
    CertificateConstantRegistry registry;
    const auto upper = [&](const std::string &name, double value) {
        CertificateConstant constant;
        constant.name = name;
        constant.value = value;
        constant.direction = CertificateBoundDirection::Upper;
        constant.source = "R2a integration-test floating constant";
        constant.derivation = "diagnostic value, not an interval proof";
        constant.mesh_class = "runtime R2a mesh";
        constant.verified = false;
        registry.set(std::move(constant));
    };
    upper("C_app", 0.1);
    upper("C_st", 2.0);
    upper("C_sd", 2.0);
    upper("C_ov", 2.0);
    upper("C_a", 2.0);
    upper("C_loc", 2.0);
    upper("beta", 0.5);
    CertificateConstant shift;
    shift.name = "s";
    shift.value = 1.0;
    shift.direction = CertificateBoundDirection::Exact;
    shift.source = "R2a integration-test floating constant";
    shift.derivation = "diagnostic shift, not an interval proof";
    shift.mesh_class = "runtime R2a mesh";
    shift.verified = false;
    registry.set(std::move(shift));
    return registry;
}

NumericalCertifiedBackendConfig r1_backend_config() {
    const PaperCaseData data = make_paper_case(PaperCase::R2a, 1.0);
    NumericalCertifiedBackendConfig config;
    config.problem_id = "paper-R2a-k1-small";
    config.source_id = "R2a-localized-gaussian-source-v1";
    config.initial_mesh = data.initial_mesh;
    config.source = data.source;
    config.initial_coarse_level = 0;
    config.initial_fine_level = 1;
    config.initial_oversampling = 1;
    config.wavenumber = data.wavenumber;
    config.quadrature_context = data.quadrature_context;
    config.constants = conditional_constants();
    config.certificate.q0 = 1e6;
    config.coarse_doerfler_theta = 0.5;
    config.audit_doerfler_theta = 0.5;
    return config;
}

CertifiedDriverConfig driver_config(CertifiedEvidencePolicy policy) {
    CertifiedDriverConfig config;
    config.method = CertifiedMethod::Calod;
    config.error_target = CertifiedErrorTarget::AuditSpace;
    config.evidence_policy = policy;
    config.mu0 = 0.5;
    config.q0 = 1e6;
    config.tolerance = 1e6;
    config.limits.max_state_transitions = 20;
    return config;
}

void verify_real_wp3_wp4_and_audit_chain_is_conditional() {
    NumericalCertifiedBackend backend(r1_backend_config());
    require(backend.full_rebuild_count() == 1,
            "production backend did not perform its initial full rebuild");
    require(backend.certificate_context().complete(),
            "production backend did not bind a complete WP4 context");

    const CoarseAdmissibilityObservation admissibility =
        backend.inspect_coarse_admissibility();
    require(admissibility.evidence.valid,
            "R1 coarse admissibility observation is invalid");
    require(admissibility.evidence.level
                == CertifiedEvidenceLevel::Conditional,
            "floating C_app was promoted to Verified");
    require(!admissibility.mu_by_element.empty(),
            "R1 coarse admissibility produced no local values");

    const CorrectorCertificationObservation corrector =
        backend.inspect_corrector_certification();
    require(corrector.evidence.valid,
            "R1 WP4 corrector diagnostics are invalid");
    require(corrector.evidence.level == CertifiedEvidenceLevel::Conditional,
            "floating WP4 corrector diagnostics were promoted to Verified");
    require(std::isfinite(corrector.q_total_upper),
            "R1 WP4 corrector q bound is not finite");
    require(!corrector.eta_h_element_squared.empty(),
            "R1 WP4 corrector marking data are empty");

    const CoarseErrorObservation coarse =
        backend.solve_and_estimate_coarse_error();
    require(coarse.evidence.valid,
            "R1 WP3/WP4 coarse error observation is invalid");
    require(coarse.evidence.level == CertifiedEvidenceLevel::Conditional,
            "floating WP3/WP4 error observation was promoted to Verified");
    require(std::isfinite(coarse.eta_H) && coarse.eta_H >= 0.0,
            "R1 WP3 eta_H is invalid");
    require(!coarse.evidence.hash.empty(),
            "R1 WP3 evidence fingerprint was not propagated");
    require(coarse.eta_H_element_squared.size()
                == backend.hierarchy().coarse_mesh().elems.size(),
            "R1 WP3 local allocation does not match the coarse mesh");

    const AuditControlObservation audit = backend.inspect_audit_control();
    require(audit.evidence.valid,
            "R1 empirical audit interval is invalid");
    require(audit.evidence.level == CertifiedEvidenceLevel::Conditional,
            "empirical saturation audit was promoted to Verified");
    require(audit.interval_kind == AuditIntervalKind::EmpiricalSaturation,
            "R1 audit did not identify the empirical interval kind");
    require(audit.audit_error_upper >= audit.audit_error_lower,
            "R1 empirical audit interval has reversed endpoints");
    require(backend.observation_count() == 4,
            "not every production observation path was exercised");

    const std::string before = backend.work_snapshot().state_fingerprint;
    backend.increase_oversampling();
    require(backend.full_rebuild_count() == 2,
            "oversampling mutation did not trigger a full rebuild");
    require(backend.work_snapshot().state_fingerprint != before,
            "oversampling mutation did not change the deterministic state fingerprint");
}

void verify_driver_calls_backend_and_strict_policy_rejects() {
    NumericalCertifiedBackend backend(r1_backend_config());
    CertifiedAdaptiveDriver driver(
        driver_config(CertifiedEvidencePolicy::RequireVerified));
    driver.run(backend);
    require(driver.terminal(), "strict R1 driver did not terminate");
    require(backend.observation_count() == 1,
            "driver did not call the production coarse observation");
    require(driver.checkpoint().state == CertifiedDriverState::Failure,
            "strict floating backend did not fail closed");
    require(driver.termination()->code
                == CertifiedStopCode::UnverifiedEvidence,
            "strict floating backend used the wrong structured stop code");
    require(std::none_of(
                driver.checkpoint().history.begin(),
                driver.checkpoint().history.end(),
                [](const CertifiedTransitionRecord &record) {
                    return record.action
                            == CertifiedDriverAction::CompleteAuditTolerance
                        || record.action
                            == CertifiedDriverAction::CompleteContinuousTolerance;
                }),
            "strict floating backend reported a certified completion");
}

template <class Mutation>
void require_full_rebuild_after(Mutation mutate, const char *message) {
    NumericalCertifiedBackend backend(r1_backend_config());
    const std::string before = backend.work_snapshot().state_fingerprint;
    mutate(backend);
    require(backend.full_rebuild_count() == 2, message);
    require(backend.work_snapshot().state_fingerprint != before,
            "mutation did not change the deterministic backend fingerprint");
}

void verify_every_mutation_uses_full_rebuild() {
    require_full_rebuild_after(
        [](NumericalCertifiedBackend &backend) {
            backend.refine_coarse({0});
        },
        "coarse mutation did not trigger a full rebuild");
    require_full_rebuild_after(
        [](NumericalCertifiedBackend &backend) {
            backend.refine_corrector_patches({0});
        },
        "corrector-fine mutation did not trigger a full rebuild");
    require_full_rebuild_after(
        [](NumericalCertifiedBackend &backend) {
            backend.increase_oversampling();
        },
        "oversampling mutation did not trigger a full rebuild");
    require_full_rebuild_after(
        [](NumericalCertifiedBackend &backend) {
            backend.refine_audit({0});
        },
        "audit mutation did not trigger a full rebuild");
}

void verify_conditional_checkpoint_replay() {
    CertifiedDriverConfig config =
        driver_config(CertifiedEvidencePolicy::AllowConditional);
    // Stop immediately after replay at an honest, journalled work limit.  The
    // direct integration test above exercises the expensive WP3/WP4/audit
    // observations independently of this deterministic replay assertion.
    config.limits.max_state_transitions = 1;
    NumericalCertifiedBackend interrupted_backend(r1_backend_config());
    CertifiedAdaptiveDriver interrupted(config);
    require(interrupted.step(interrupted_backend),
            "conditional R1 driver stopped before its first transition");
    require(interrupted.checkpoint().state
                == CertifiedDriverState::CorrectorCertification,
            "conditional R1 driver did not accept coarse admissibility");
    require(interrupted.checkpoint().claim
                == CertifiedEvidenceLevel::Conditional,
            "conditional R1 evidence did not downgrade the run claim");

    const std::string checkpoint =
        serialize_certified_checkpoint(interrupted.checkpoint());
    NumericalCertifiedBackend fresh_backend(r1_backend_config());
    CertifiedAdaptiveDriver resumed = CertifiedAdaptiveDriver::resume(
        config, checkpoint, fresh_backend);
    require(!resumed.terminal(),
            "deterministic production backend checkpoint replay failed");
    require(resumed.checkpoint().backend_fingerprint
                == fresh_backend.work_snapshot().state_fingerprint,
            "replayed backend fingerprint differs from the checkpoint");
    resumed.run(fresh_backend);
    require(resumed.checkpoint().state == CertifiedDriverState::WorkLimit,
            "replayed conditional backend did not stop at its work limit");
    require(resumed.termination()->code
                == CertifiedStopCode::StateTransitionLimit,
            "replayed conditional backend used the wrong work-limit code");
    require(resumed.termination()->claim
                == CertifiedEvidenceLevel::Conditional,
            "replayed conditional backend lost its evidence level");
}

void verify_audit_control_checkpoint_replay() {
    CertifiedDriverConfig config =
        driver_config(CertifiedEvidencePolicy::AllowConditional);
    config.error_target = CertifiedErrorTarget::Continuous;
    config.rho_aud = 0.99;
    NumericalCertifiedBackend interrupted_backend(r1_backend_config());
    CertifiedAdaptiveDriver interrupted(config);
    require(interrupted.step(interrupted_backend),
            "R1 run stopped before corrector certification");
    require(interrupted.step(interrupted_backend),
            "R1 run stopped before coarse error control");
    require(interrupted.checkpoint().state
                == CertifiedDriverState::CoarseErrorControl,
            "conditional R1 corrector diagnostics were not accepted");
    if (!interrupted.step(interrupted_backend)) {
        throw std::runtime_error(
            "R1 run stopped before audit control: "
            + interrupted.termination()->detail);
    }
    require(interrupted.checkpoint().state
                == CertifiedDriverState::AuditControl,
            "R1 checkpoint was not captured at AuditControl");
    const std::uint64_t saved_work =
        interrupted.checkpoint().cumulative_backend_work_units;
    const double saved_elapsed =
        interrupted.checkpoint().cumulative_elapsed_seconds;
    const std::uint64_t saved_peak =
        interrupted.checkpoint().cumulative_peak_memory_bytes;

    NumericalCertifiedBackend fresh_backend(r1_backend_config());
    CertifiedAdaptiveDriver resumed = CertifiedAdaptiveDriver::resume(
        config, serialize_certified_checkpoint(interrupted.checkpoint()),
        fresh_backend);
    require(!resumed.terminal(),
            "AuditControl checkpoint replay failed before its observation");
    require(resumed.checkpoint().cumulative_backend_work_units == saved_work,
            "production checkpoint replay reset or double-counted work");
    require(resumed.checkpoint().cumulative_elapsed_seconds == saved_elapsed,
            "production checkpoint replay reset cumulative elapsed time");
    require(resumed.checkpoint().cumulative_peak_memory_bytes >= saved_peak,
            "production checkpoint replay reset peak memory");
    resumed.step(fresh_backend);
    require(fresh_backend.observation_count() == 1,
            "resumed AuditControl did not reproduce the audit observation");
    require(resumed.checkpoint().state != CertifiedDriverState::Failure,
            "resumed AuditControl failed after deterministic reconstruction");
    require(resumed.checkpoint().claim
                == CertifiedEvidenceLevel::Conditional,
            "resumed empirical audit was promoted to Verified");
    require(resumed.checkpoint().cumulative_backend_work_units > saved_work,
            "resumed audit work was not added to the cumulative budget");
    require(resumed.checkpoint().cumulative_elapsed_seconds >= saved_elapsed,
            "resumed audit elapsed time regressed");
}

} // namespace

int main() {
    try {
        verify_real_wp3_wp4_and_audit_chain_is_conditional();
        verify_driver_calls_backend_and_strict_policy_rejects();
        verify_every_mutation_uses_full_rebuild();
        verify_conditional_checkpoint_replay();
        verify_audit_control_checkpoint_replay();
        std::cout << "Conditional numerical certified backend integration passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_numerical_backend failed: "
                  << error.what() << '\n';
        return 1;
    }
}

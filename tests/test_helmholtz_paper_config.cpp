#include "helmholtz/experiments/paper_config.h"

#include <array>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

using namespace lod2d::helmholtz::experiments;

namespace {

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

template <class Function>
void require_invalid(Function &&function, const char *message) {
    try {
        function();
    } catch (const std::invalid_argument &) {
        return;
    }
    throw std::runtime_error(message);
}

PaperConfig sample_config() {
    PaperConfig config;
    config.case_id = PaperCase::R2b;
    config.method_id = PaperMethod::Calod;
    config.wavenumber = 16.0;
    config.theta_H = 0.7;
    config.error_target = lod2d::helmholtz::adaptive::CertifiedErrorTarget::Continuous;
    config.evidence_policy =
        lod2d::helmholtz::adaptive::CertifiedEvidencePolicy::AllowConditional;
    config.mu0 = 0.45;
    config.q0 = 0.2;
    config.tau = 0.4;
    config.theta_h = 0.6;
    config.rho_aud = 0.04;
    config.tolerance = 0.02;
    config.numerical_backend.initial_coarse_level = 1;
    config.numerical_backend.initial_fine_level = 3;
    config.numerical_backend.initial_oversampling = 2;
    config.numerical_backend.boundary_beta = 1.25;
    config.numerical_backend.petrov_mode =
        lod2d::helmholtz::HelmholtzPetrovMode::CorrectedTestOnly;
    config.numerical_backend.patch_solver.kind =
        lod2d::helmholtz::HelmholtzPatchSolverKind::ShiftedGmres;
    config.numerical_backend.patch_solver.symbolic_cache_slots = 3;
    config.numerical_backend.patch_solver.maximum_parallel_solves = 2;
    config.numerical_backend.patch_solver.reuse_identical_factorization = true;
    config.numerical_backend.patch_solver.gmres.restart = 17;
    config.numerical_backend.patch_solver.gmres.max_iterations = 91;
    config.numerical_backend.patch_solver.gmres.relative_tolerance = 2e-11;
    config.numerical_backend.patch_solver.gmres.absolute_tolerance = 3e-14;
    config.numerical_backend.patch_solver.gmres.reorthogonalize = false;
    config.numerical_backend.patch_solver.shifted.rule =
        lod2d::helmholtz::HelmholtzShiftRule::PatchScaled;
    config.numerical_backend.patch_solver.shifted.alpha = 0.3;
    config.numerical_backend.patch_solver.shifted.absolute_epsilon = 0.01;
    config.numerical_backend.patch_solver.shifted.inverse =
        lod2d::helmholtz::HelmholtzShiftedInverseKind::GeometricVcycle;
    config.numerical_backend.patch_solver.shifted.pre_smooth = 3;
    config.numerical_backend.patch_solver.shifted.post_smooth = 4;
    config.numerical_backend.patch_solver.shifted.coarse_max_dofs = 123;
    config.numerical_backend.patch_solver.shifted.jacobi_weight = 0.7;
    config.numerical_backend.patch_solver.fallback_to_direct = true;
    config.numerical_backend.certificate.precision_bits = 192;
    config.numerical_backend.certificate.cluster_relative_gap = 2e-8;
    config.numerical_backend.certificate.cluster_absolute_gap = 3e-12;
    config.numerical_backend.certificate.conjugation_tolerance = 4e-10;
    config.numerical_backend.certificate.q0 = config.q0;
    config.numerical_backend.kernel_riesz_solver =
        lod2d::helmholtz::adaptive::KernelRieszSolver::KernelBasisReference;
    config.numerical_backend.audit_doerfler_theta = 0.65;
    config.numerical_backend.audit_saturation_factor = 0.04;
    config.numerical_backend.certificate_constant_set_hash =
        "sha256:77d9b551711f339992234ae21b2c5dd9809f2440b7ae37251626b04d733ad9da";
    config.work_limits.max_state_transitions = 80;
    config.work_limits.max_coarse_refinements = 11;
    config.work_limits.max_corrector_refinements = 12;
    config.work_limits.max_oversampling_increments = 13;
    config.work_limits.max_audit_refinements = 14;
    config.work_limits.max_coarse_dofs = 1000;
    config.work_limits.max_fine_dofs = 2000;
    config.work_limits.max_audit_dofs = 3000;
    config.work_limits.max_backend_work_units = 4000;
    config.work_limits.max_peak_memory_bytes = 5000;
    config.work_limits.max_elapsed_seconds = 60.0;
    config.work_limits.max_oversampling = 9;
    config.hlod_prior_corrector_space_id = "prior-space-v1";
    config.hlod_prior_oversampling = 4;
    config.repeat_index = 3;
    config.git_commit = "c8b5520d4f2b1e0f225a64e837dce5220daa1a2d";
    config.build_hash = "gcc-13-release-abc123";
    return config;
}

template <class Mutator>
void require_identity_change(
    const PaperConfig &original,
    Mutator &&mutator,
    const char *field) {
    PaperConfig changed = original;
    mutator(changed);
    require(canonical_config_hash(changed) != canonical_config_hash(original), field);
    require(make_run_id(changed) != make_run_id(original), field);
}

void verify_registries() {
    require(paper_case_registry().size() == 4, "paper case registry is incomplete");
    require(paper_method_registry().size() == 7, "paper method registry is incomplete");
    require(case_definition(PaperCase::R2a).gaussian_sigma == 1.0 / 32.0,
            "R2a sigma is not frozen to 2^-5");
    require(case_definition(PaperCase::R2b).gaussian_sigma == 1.0 / 64.0,
            "R2b sigma is not frozen to 2^-6");
    require(case_definition(PaperCase::S).has_mixed_boundary,
            "S must declare mixed boundary data");
    require(method_definition(PaperMethod::HlodProxy).diagnostic_only,
            "HLOD-proxy must remain diagnostic");
    require(!method_definition(PaperMethod::HlodProxy).paper_comparator,
            "HLOD-proxy must not enter the paper comparator matrix");
}

void verify_round_trip_and_hash() {
    const PaperConfig original = sample_config();
    const std::string encoded = canonical_json(original);
    const PaperConfig decoded = parse_paper_config(encoded);
    require(decoded == original, "paper config JSON round trip lost fields");
    require(canonical_json(decoded) == encoded, "canonical JSON changed after round trip");
    require(canonical_config_hash(decoded) == canonical_config_hash(original),
            "canonical hash changed after round trip");
    require(make_run_id(original) == make_run_id(decoded),
            "run ID is not deterministic");

    PaperConfig changed = original;
    changed.repeat_index = 4;
    require(canonical_config_hash(changed) != canonical_config_hash(original),
            "canonical hash ignores repeat index");
    require(make_run_id(changed) != make_run_id(original),
            "run ID ignores repeat index");
}

void verify_driver_config_is_part_of_identity() {
    const PaperConfig original = sample_config();
    const auto driver = make_certified_driver_config(original);
    require(driver.method == lod2d::helmholtz::adaptive::CertifiedMethod::Calod,
            "paper method was not mapped into driver config");
    require(driver.error_target == original.error_target
                && driver.evidence_policy == original.evidence_policy
                && driver.mu0 == original.mu0 && driver.q0 == original.q0
                && driver.tau == original.tau && driver.theta_H == original.theta_H
                && driver.theta_h == original.theta_h
                && driver.rho_aud == original.rho_aud
                && driver.tolerance == original.tolerance,
            "paper decision parameters drifted while making driver config");
    require(driver.limits.max_state_transitions
                    == original.work_limits.max_state_transitions
                && driver.limits.max_coarse_refinements
                    == original.work_limits.max_coarse_refinements
                && driver.limits.max_corrector_refinements
                    == original.work_limits.max_corrector_refinements
                && driver.limits.max_oversampling_increments
                    == original.work_limits.max_oversampling_increments
                && driver.limits.max_audit_refinements
                    == original.work_limits.max_audit_refinements
                && driver.limits.max_coarse_dofs
                    == original.work_limits.max_coarse_dofs
                && driver.limits.max_fine_dofs
                    == original.work_limits.max_fine_dofs
                && driver.limits.max_audit_dofs
                    == original.work_limits.max_audit_dofs
                && driver.limits.max_backend_work_units
                    == original.work_limits.max_backend_work_units
                && driver.limits.max_peak_memory_bytes
                    == original.work_limits.max_peak_memory_bytes
                && driver.limits.max_elapsed_seconds
                    == original.work_limits.max_elapsed_seconds
                && driver.limits.max_oversampling
                    == original.work_limits.max_oversampling,
            "paper work limits drifted while making driver config");
    require(driver.hlod_prior_corrector_space_id
                    == original.hlod_prior_corrector_space_id
                && driver.hlod_prior_oversampling
                    == original.hlod_prior_oversampling,
            "paper HLOD prior drifted while making driver config");

    require_identity_change(original, [](PaperConfig &value) {
        value.error_target = lod2d::helmholtz::adaptive::CertifiedErrorTarget::AuditSpace;
    }, "run identity ignores error_target");
    require_identity_change(original, [](PaperConfig &value) {
        value.method_id = PaperMethod::Hlod;
    }, "run identity ignores certified method");
    require_identity_change(original, [](PaperConfig &value) {
        value.evidence_policy =
            lod2d::helmholtz::adaptive::CertifiedEvidencePolicy::RequireVerified;
    }, "run identity ignores evidence_policy");
    require_identity_change(original, [](PaperConfig &value) { value.mu0 += 0.01; },
                            "run identity ignores mu0");
    require_identity_change(original, [](PaperConfig &value) {
        value.q0 += 0.01;
        value.numerical_backend.certificate.q0 = value.q0;
    },
                            "run identity ignores q0");
    require_identity_change(original, [](PaperConfig &value) { value.tau += 0.01; },
                            "run identity ignores tau");
    require_identity_change(original, [](PaperConfig &value) { value.theta_H = 0.5; },
                            "run identity ignores theta_H");
    require_identity_change(original, [](PaperConfig &value) { value.theta_h += 0.01; },
                            "run identity ignores theta_h");
    require_identity_change(original, [](PaperConfig &value) { value.rho_aud += 0.01; },
                            "run identity ignores rho_aud");
    require_identity_change(original, [](PaperConfig &value) { value.tolerance += 0.01; },
                            "run identity ignores tolerance");
    require_identity_change(original, [](PaperConfig &value) {
        ++value.work_limits.max_state_transitions;
    }, "run identity ignores max_state_transitions");
    require_identity_change(original, [](PaperConfig &value) {
        ++value.work_limits.max_coarse_refinements;
    }, "run identity ignores max_coarse_refinements");
    require_identity_change(original, [](PaperConfig &value) {
        ++value.work_limits.max_corrector_refinements;
    }, "run identity ignores max_corrector_refinements");
    require_identity_change(original, [](PaperConfig &value) {
        ++value.work_limits.max_oversampling_increments;
    }, "run identity ignores max_oversampling_increments");
    require_identity_change(original, [](PaperConfig &value) {
        ++value.work_limits.max_audit_refinements;
    }, "run identity ignores max_audit_refinements");
    require_identity_change(original, [](PaperConfig &value) {
        ++value.work_limits.max_coarse_dofs;
    }, "run identity ignores max_coarse_dofs");
    require_identity_change(original, [](PaperConfig &value) {
        ++value.work_limits.max_fine_dofs;
    }, "run identity ignores max_fine_dofs");
    require_identity_change(original, [](PaperConfig &value) {
        ++value.work_limits.max_audit_dofs;
    }, "run identity ignores max_audit_dofs");
    require_identity_change(original, [](PaperConfig &value) {
        ++value.work_limits.max_backend_work_units;
    }, "run identity ignores max_backend_work_units");
    require_identity_change(original, [](PaperConfig &value) {
        ++value.work_limits.max_peak_memory_bytes;
    }, "run identity ignores max_peak_memory_bytes");
    require_identity_change(original, [](PaperConfig &value) {
        value.work_limits.max_elapsed_seconds += 1.0;
    }, "run identity ignores max_elapsed_seconds");
    require_identity_change(original, [](PaperConfig &value) {
        ++value.work_limits.max_oversampling;
    }, "run identity ignores max_oversampling");
    require_identity_change(original, [](PaperConfig &value) {
        value.hlod_prior_corrector_space_id += "-changed";
    }, "run identity ignores HLOD corrector-space prior");
    require_identity_change(original, [](PaperConfig &value) {
        ++value.hlod_prior_oversampling;
    }, "run identity ignores HLOD oversampling prior");

    PaperConfig hlod = original;
    hlod.method_id = PaperMethod::Hlod;
    require(make_certified_driver_config(hlod).method
                == lod2d::helmholtz::adaptive::CertifiedMethod::Hlod,
            "HLOD paper method was not mapped into driver config");
}

void verify_numerical_backend_config_is_part_of_identity() {
    const PaperConfig original = sample_config();
    require_identity_change(original, [](PaperConfig &v) {
        --v.numerical_backend.initial_coarse_level;
    }, "run identity ignores initial coarse level");
    require_identity_change(original, [](PaperConfig &v) {
        ++v.numerical_backend.initial_fine_level;
    }, "run identity ignores initial fine level");
    require_identity_change(original, [](PaperConfig &v) {
        ++v.numerical_backend.initial_oversampling;
    }, "run identity ignores initial oversampling");
    require_identity_change(original, [](PaperConfig &v) {
        v.numerical_backend.boundary_beta += 0.1;
    }, "run identity ignores boundary beta");
    require_identity_change(original, [](PaperConfig &v) {
        v.numerical_backend.petrov_mode =
            lod2d::helmholtz::HelmholtzPetrovMode::TwoSided;
    }, "run identity ignores Petrov mode");
    require_identity_change(original, [](PaperConfig &v) {
        v.numerical_backend.patch_solver.kind =
            lod2d::helmholtz::HelmholtzPatchSolverKind::DirectSchur;
    }, "run identity ignores patch solver kind");
    require_identity_change(original, [](PaperConfig &v) {
        ++v.numerical_backend.patch_solver.symbolic_cache_slots;
    }, "run identity ignores patch symbolic cache slots");
    require_identity_change(original, [](PaperConfig &v) {
        ++v.numerical_backend.patch_solver.maximum_parallel_solves;
    }, "run identity ignores maximum parallel patch solves");
    require_identity_change(original, [](PaperConfig &v) {
        v.numerical_backend.patch_solver.reuse_identical_factorization = false;
    }, "run identity ignores factorization reuse");
    require_identity_change(original, [](PaperConfig &v) {
        ++v.numerical_backend.patch_solver.gmres.restart;
    }, "run identity ignores GMRES restart");
    require_identity_change(original, [](PaperConfig &v) {
        ++v.numerical_backend.patch_solver.gmres.max_iterations;
    }, "run identity ignores GMRES iteration limit");
    require_identity_change(original, [](PaperConfig &v) {
        v.numerical_backend.patch_solver.gmres.relative_tolerance *= 2.0;
    }, "run identity ignores GMRES relative tolerance");
    require_identity_change(original, [](PaperConfig &v) {
        v.numerical_backend.patch_solver.gmres.absolute_tolerance *= 2.0;
    }, "run identity ignores GMRES absolute tolerance");
    require_identity_change(original, [](PaperConfig &v) {
        v.numerical_backend.patch_solver.gmres.reorthogonalize = true;
    }, "run identity ignores GMRES reorthogonalization");
    require_identity_change(original, [](PaperConfig &v) {
        v.numerical_backend.patch_solver.shifted.rule =
            lod2d::helmholtz::HelmholtzShiftRule::Absolute;
    }, "run identity ignores shift rule");
    require_identity_change(original, [](PaperConfig &v) {
        v.numerical_backend.patch_solver.shifted.alpha += 0.1;
    }, "run identity ignores shift alpha");
    require_identity_change(original, [](PaperConfig &v) {
        v.numerical_backend.patch_solver.shifted.absolute_epsilon += 0.01;
    }, "run identity ignores absolute shift");
    require_identity_change(original, [](PaperConfig &v) {
        v.numerical_backend.patch_solver.shifted.inverse =
            lod2d::helmholtz::HelmholtzShiftedInverseKind::SparseLu;
    }, "run identity ignores shifted inverse");
    require_identity_change(original, [](PaperConfig &v) {
        ++v.numerical_backend.patch_solver.shifted.pre_smooth;
    }, "run identity ignores pre-smoothing");
    require_identity_change(original, [](PaperConfig &v) {
        ++v.numerical_backend.patch_solver.shifted.post_smooth;
    }, "run identity ignores post-smoothing");
    require_identity_change(original, [](PaperConfig &v) {
        ++v.numerical_backend.patch_solver.shifted.coarse_max_dofs;
    }, "run identity ignores shifted coarse size");
    require_identity_change(original, [](PaperConfig &v) {
        v.numerical_backend.patch_solver.shifted.jacobi_weight += 0.01;
    }, "run identity ignores Jacobi weight");
    require_identity_change(original, [](PaperConfig &v) {
        v.numerical_backend.patch_solver.fallback_to_direct = false;
    }, "run identity ignores direct fallback");
    require_identity_change(original, [](PaperConfig &v) {
        ++v.numerical_backend.certificate.precision_bits;
    }, "run identity ignores certificate precision");
    require_identity_change(original, [](PaperConfig &v) {
        v.numerical_backend.certificate.cluster_relative_gap *= 2.0;
    }, "run identity ignores relative cluster gap");
    require_identity_change(original, [](PaperConfig &v) {
        v.numerical_backend.certificate.cluster_absolute_gap *= 2.0;
    }, "run identity ignores absolute cluster gap");
    require_identity_change(original, [](PaperConfig &v) {
        v.numerical_backend.certificate.conjugation_tolerance *= 2.0;
    }, "run identity ignores conjugation tolerance");
    require_identity_change(original, [](PaperConfig &v) {
        v.numerical_backend.kernel_riesz_solver =
            lod2d::helmholtz::adaptive::KernelRieszSolver::SaddlePoint;
    }, "run identity ignores kernel Riesz solver");
    require_identity_change(original, [](PaperConfig &v) {
        v.numerical_backend.audit_doerfler_theta += 0.01;
    }, "run identity ignores audit marking theta");
    require_identity_change(original, [](PaperConfig &v) {
        v.numerical_backend.audit_saturation_factor += 0.01;
    }, "run identity ignores audit saturation factor");
    require_identity_change(original, [](PaperConfig &v) {
        v.numerical_backend.certificate_constant_set_hash.back() = 'b';
    }, "run identity ignores certificate constant set");
}

void verify_strict_validation() {
    const std::string encoded = canonical_json(sample_config());
    std::string unknown = encoded;
    unknown.insert(unknown.size() - 1, ",\"future_field\":1");
    require_invalid([&] { (void)parse_paper_config(unknown); },
                    "unknown config field was accepted");

    std::string wrong_version = encoded;
    const std::string needle = "\"schema_version\":1";
    wrong_version.replace(wrong_version.find(needle), needle.size(), "\"schema_version\":2");
    require_invalid([&] { (void)parse_paper_config(wrong_version); },
                    "unknown schema version was accepted");

    PaperConfig invalid = sample_config();
    invalid.wavenumber = 7.0;
    require_invalid([&] { (void)canonical_json(invalid); },
                    "non-protocol wavenumber was accepted");

    const std::string number_needle = "\"wavenumber\":16";
    for (const std::string replacement : {
             "\"wavenumber\":+16", "\"wavenumber\":016",
             "\"wavenumber\":16.", "\"wavenumber\":16e"}) {
        std::string malformed = encoded;
        malformed.replace(
            malformed.find(number_needle), number_needle.size(), replacement);
        require_invalid([&] { (void)parse_paper_config(malformed); },
                        "non-RFC-8259 JSON number was accepted");
    }

    invalid = sample_config();
    invalid.tolerances.linear_relative_residual =
        std::numeric_limits<double>::infinity();
    require_invalid([&] { (void)canonical_json(invalid); },
                    "infinite tolerance was accepted");
    invalid = sample_config();
    invalid.build_hash.push_back('\x01');
    require_invalid([&] { (void)canonical_json(invalid); },
                    "control character in provenance was accepted");

    invalid = sample_config();
    invalid.work_limits.max_state_transitions = 0;
    require_invalid([&] { (void)canonical_json(invalid); },
                    "zero state-transition limit was accepted");

    invalid = sample_config();
    invalid.numerical_backend.initial_fine_level = 0;
    require_invalid([&] { (void)canonical_json(invalid); },
                    "backend fine level below coarse level was accepted");
    invalid = sample_config();
    invalid.numerical_backend.certificate.q0 += 0.01;
    require_invalid([&] { (void)canonical_json(invalid); },
                    "driver/backend q0 drift was accepted");
    invalid = sample_config();
    invalid.numerical_backend.certificate.precision_bits = 63;
    require_invalid([&] { (void)canonical_json(invalid); },
                    "certificate precision below 64 bits was accepted");
    invalid = sample_config();
    invalid.numerical_backend.certificate_constant_set_hash.clear();
    require_invalid([&] { (void)canonical_json(invalid); },
                    "missing certificate constant set digest was accepted");
    invalid = sample_config();
    invalid.numerical_backend.certificate_constant_set_hash =
        "constant-set-label";
    require_invalid([&] { (void)canonical_json(invalid); },
                    "non-digest certificate constant set label was accepted");

    invalid = sample_config();
    invalid.method_id = PaperMethod::Hlod;
    invalid.hlod_prior_corrector_space_id.clear();
    require_invalid([&] { (void)canonical_json(invalid); },
                    "HLOD without a frozen prior was accepted");
}

void verify_status_contract() {
    for (PaperRunStatus status : std::array{
             PaperRunStatus::Success, PaperRunStatus::Interrupted,
             PaperRunStatus::CensoredWorkLimit, PaperRunStatus::CensoredMemoryLimit,
             PaperRunStatus::CensoredTimeLimit, PaperRunStatus::CensoredIterationLimit,
             PaperRunStatus::LinearAlgebraFailure, PaperRunStatus::CertificateFailure,
             PaperRunStatus::Unavailable}) {
        require(parse_paper_run_status(to_string(status)) == status,
                "paper run status does not round trip");
    }
    for (PaperValueStatus status : std::array{
             PaperValueStatus::Valid, PaperValueStatus::NotApplicable,
             PaperValueStatus::NotComputed, PaperValueStatus::InvalidDenominator,
             PaperValueStatus::EnclosureFailed}) {
        require(parse_paper_value_status(to_string(status)) == status,
                "paper value status does not round trip");
    }
    for (PaperCertificateStatus status : std::array{
             PaperCertificateStatus::ImplementationStudy,
             PaperCertificateStatus::Conditional,
             PaperCertificateStatus::AuditCertified,
             PaperCertificateStatus::ContinuousCertified,
             PaperCertificateStatus::EmpiricalReference}) {
        require(parse_paper_certificate_status(to_string(status)) == status,
                "paper certificate status does not round trip");
    }
    require_invalid([&] { (void)parse_paper_run_status("timeout"); },
                    "unknown paper run status was accepted");
    require_invalid([&] { (void)parse_paper_value_status("nan"); },
                    "unknown paper value status was accepted");
    require_invalid([&] { (void)parse_paper_certificate_status("certified"); },
                    "ambiguous certified output status was accepted");

    require(paper_output_metric_registry().size() == 30,
            "paper output metric registry changed without a schema version bump");
    for (const std::string_view metric : paper_output_metric_registry())
        validate_paper_output_metric(metric);
    require_invalid([&] { validate_paper_output_metric("custom_metric"); },
                    "unknown paper output metric was accepted");
}

PracticalPaperConfig sample_practical_config() {
    PracticalPaperConfig config;
    config.case_id = PaperCase::S;
    config.method_id = PracticalPaperMethod::Palod;
    config.wavenumber = 16.0;
    config.initial_coarse_level = 2;
    config.reference_level = 6;
    config.reference_epoch = 1;
    config.minimum_reference_level_gap = 0;
    config.ell0 = 1;
    config.ell_max = 5;
    config.boundary_beta = 1.25;
    config.c_H = 0.45;
    config.theta_loc = 0.2;
    config.C0_usr = 1.1;
    config.C1_usr = 1.3;
    config.rho_star = 0.2;
    config.trajectory_policy =
        PracticalTrajectoryPolicy::FixedWorkHorizon;
    config.practical_stop_tolerance = 0.0123;
    config.plateau_diagnostic.minimum_geometric_mean_ratio = 0.92;
    config.plateau_diagnostic.maximum_relative_oscillation = 0.12;
    config.plateau_diagnostic.window_steps = 2;
    config.reference_adequacy.enabled = true;
    config.reference_adequacy.maximum_terminal_error_fraction = 0.2;
    config.petrov_mode = lod2d::helmholtz::HelmholtzPetrovMode::CorrectedTestOnly;
    config.patch_solver_kind =
        lod2d::helmholtz::HelmholtzPatchSolverKind::DirectSchur;
    config.maximum_patch_threads = 2;
    config.work_limits.maximum_iterations = 40;
    config.work_limits.maximum_H_steps = 8;
    config.work_limits.maximum_unknowns = 500000;
    config.work_limits.maximum_coarse_elements = 50000;
    config.work_limits.maximum_ambient_elements = 500000;
    config.work_limits.maximum_wall_seconds = 120.0;
    config.timing_repeats = 3;
    config.repeat_index = 2;
    config.git_commit = "4abdb3c4154579dbdc1887ae9ba8ff77d9e5a810";
    config.build_hash = "gcc-release-wp5";
    config.manuscript_sha256 =
        "03d83e0eb7128aa5ef00002c6dac110f548351e52e92e56d7adf709880854d20";
    return config;
}

void verify_practical_v4_contract() {
    const PracticalPaperConfig original = sample_practical_config();
    const std::string encoded = canonical_json(original);
    const PracticalPaperConfig decoded = parse_practical_paper_config(encoded);
    require(decoded == original, "practical v4 JSON round trip lost fields");
    require(canonical_json(decoded) == encoded,
            "practical v4 canonical JSON changed after round trip");
    require(make_run_id(decoded) == make_run_id(original),
            "practical v4 run ID is not deterministic");

    for (const double wavenumber : {2.0, 4.0, 8.0, 16.0, 32.0}) {
        PracticalPaperConfig supported = original;
        supported.wavenumber = wavenumber;
        validate_practical_paper_config(supported);
    }
    PracticalPaperConfig unsupported = original;
    unsupported.wavenumber = 3.0;
    require_invalid(
        [&] { validate_practical_paper_config(unsupported); },
        "practical v4 accepted a wave number outside the frozen experiment set");

    const auto driver = make_practical_driver_config(original);
    require(driver.initial_coarse_level == original.initial_coarse_level &&
                driver.reference_level == original.reference_level &&
                driver.reference_epoch == original.reference_epoch &&
                driver.minimum_reference_level_gap ==
                    original.minimum_reference_level_gap &&
                driver.ell0 == original.ell0 && driver.ell_max == original.ell_max &&
                driver.theta_loc == original.theta_loc &&
                driver.C0_usr == original.C0_usr &&
                driver.C1_usr == original.C1_usr &&
                driver.theta_H == original.theta_H &&
                driver.rho_star == original.rho_star &&
                driver.tolerance_reference ==
                    original.practical_stop_tolerance &&
                driver.patch_solver.maximum_parallel_solves ==
                    original.maximum_patch_threads &&
                driver.stop_policy ==
                    lod2d::helmholtz::adaptive::
                        PracticalStopPolicy::FixedWorkHorizon &&
                driver.limits.maximum_unknowns ==
                    original.work_limits.maximum_unknowns,
            "practical v4 fields drifted while making the driver config");

    PracticalPaperConfig changed = original;
    changed.reference_level += 1;
    require(canonical_config_hash(changed) != canonical_config_hash(original),
            "practical v4 identity ignores reference_mesh level");
    changed = original;
    changed.maximum_patch_threads = 1;
    require(canonical_config_hash(changed) != canonical_config_hash(original),
            "practical v4 identity ignores the patch concurrency cap");
    changed = original;
    changed.singular_oscillatory_fraction = 0.0;
    changed.singular_cutoff_outer_radius = 1.0;
    changed.singular_quintic_cutoff = true;
    changed.smooth_wave_amplitude = 0.1;
    require(parse_practical_paper_config(canonical_json(changed)) == changed
                && canonical_config_hash(changed)
                    != canonical_config_hash(original),
            "practical v4 identity ignores the S oscillatory fraction");
    changed.case_id = PaperCase::R1;
    require_invalid([&] { (void)canonical_json(changed); },
                    "a non-S practical case accepted an S exact-solution parameter");
    changed = original;
    changed.singular_oscillatory_fraction = 1.01;
    require_invalid([&] { (void)canonical_json(changed); },
                    "practical v4 accepted an invalid S oscillatory fraction");
    changed = original;
    changed.singular_cutoff_outer_radius = 1.01;
    require_invalid([&] { (void)canonical_json(changed); },
                    "practical v4 accepted an invalid S cut-off radius");
    changed = original;
    changed.smooth_wave_amplitude = 1.01;
    require_invalid([&] { (void)canonical_json(changed); },
                    "practical v4 accepted an invalid smooth wave amplitude");
    changed = original;
    changed.minimum_reference_level_gap = 2;
    require(canonical_config_hash(changed) != canonical_config_hash(original),
            "practical v4 identity ignores the minimum reference level gap");
    changed = original;
    changed.reference_epoch += 1;
    require(canonical_config_hash(changed) != canonical_config_hash(original),
            "practical v4 identity ignores the reference epoch");
    require(make_practical_driver_config(changed).reference_epoch
                == changed.reference_epoch,
            "practical v4 reference epoch did not reach the driver");
    changed = original;
    changed.rho_star = 0.3;
    require(canonical_config_hash(changed) != canonical_config_hash(original),
            "practical v4 identity ignores ambient ratio policy");
    changed = original;
    changed.C1_usr += 0.1;
    require(canonical_config_hash(changed) != canonical_config_hash(original),
            "practical v4 identity ignores user localization constant");
    changed = original;
    changed.work_limits.maximum_unknowns += 1;
    require(canonical_config_hash(changed) != canonical_config_hash(original),
            "practical v4 identity ignores work limits");
    changed = original;
    changed.practical_stop_tolerance *= 2.0;
    require(canonical_config_hash(changed) != canonical_config_hash(original),
            "practical v4 identity ignores the independent stop tolerance");
    changed = original;
    changed.plateau_diagnostic.minimum_geometric_mean_ratio = 0.9;
    require(canonical_config_hash(changed) != canonical_config_hash(original),
            "practical v4 identity ignores the plateau diagnostic policy");
    changed = original;
    changed.reference_adequacy.maximum_terminal_error_fraction = 0.3;
    require(canonical_config_hash(changed) != canonical_config_hash(original),
            "practical v4 identity ignores the reference adequacy policy");
    changed = original;
    changed.reference_refresh_H_steps = {2, 5};
    const std::string scheduled_encoded = canonical_json(changed);
    const PracticalPaperConfig scheduled_decoded =
        parse_practical_paper_config(scheduled_encoded);
    require(scheduled_decoded == changed
                && make_practical_driver_config(scheduled_decoded)
                       .reference_refresh_H_steps == changed.reference_refresh_H_steps
                && canonical_config_hash(changed)
                       != canonical_config_hash(original),
            "continuous reference-epoch schedule did not round trip into the driver identity");
    changed.reference_refresh_H_steps = {2, 2};
    require_invalid([&] { (void)canonical_json(changed); },
                    "practical v4 accepted a repeated reference refresh step");

    std::string with_legacy_theta_h = encoded;
    with_legacy_theta_h.insert(with_legacy_theta_h.size() - 1, ",\"theta_h\":0.5");
    require_invalid(
        [&] { (void)parse_practical_paper_config(with_legacy_theta_h); },
        "practical v4 accepted legacy theta_h");
    std::string with_legacy_q_h = encoded;
    with_legacy_q_h.insert(with_legacy_q_h.size() - 1, ",\"q_h\":0.25");
    require_invalid(
        [&] { (void)parse_practical_paper_config(with_legacy_q_h); },
        "practical v4 accepted legacy q_h");

    changed = original;
    changed.method_id = PracticalPaperMethod::HlodFixed;
    changed.ell_max = changed.ell0;
    validate_practical_paper_config(changed);
    const auto hlod_driver = make_practical_driver_config(changed);
    require(hlod_driver.localization_policy
                == lod2d::helmholtz::adaptive::
                    PracticalLocalizationPolicy::FixedGlobalEll
                && hlod_driver.ell0 == hlod_driver.ell_max,
            "HLOD-fixed did not map to the real fixed-ell backend");
    changed.method_id = PracticalPaperMethod::Ufem;
    changed.ell0 = 0;
    changed.ell_max = 0;
    validate_practical_paper_config(changed);
    require_invalid([&] { (void)make_practical_driver_config(changed); },
                    "UFEM was silently relabelled as a LOD backend");
    changed = original;
    changed.method_id = PracticalPaperMethod::Afem;
    changed.ell0 = 0;
    changed.ell_max = 0;
    validate_practical_paper_config(changed);
    require_invalid([&] { (void)make_practical_driver_config(changed); },
                    "AFEM was silently relabelled as a LOD backend");
    changed.ell_max = 1;
    require_invalid([&] { validate_practical_paper_config(changed); },
                    "AFEM accepted a nonzero localization radius");

    require(standard_lod_prior_ell(8.0) == 3
                && standard_lod_prior_ell(16.0) == 4
                && standard_lod_prior_ell(32.0) == 5,
            "frozen SLOD c_prior=1 policy drifted");
    changed = original;
    changed.method_id = PracticalPaperMethod::Slod;
    changed.ell0 = 2;
    changed.ell_max = changed.ell0;
    validate_practical_paper_config(changed);
    require_invalid([&] { (void)make_practical_driver_config(changed); },
                    "SLOD was silently relabelled as the adaptive LOD backend");
    changed.ell_max += 1;
    require_invalid([&] { validate_practical_paper_config(changed); },
                    "SLOD accepted a non-frozen empirical ell");
    changed.ell0 = 0;
    changed.ell_max = 0;
    require_invalid([&] { validate_practical_paper_config(changed); },
                    "SLOD accepted a nonpositive empirical ell");
    require_invalid([&] { (void)standard_lod_prior_ell(0.0); },
                    "SLOD accepted a nonpositive wavenumber");
}

} // namespace

int main() {
    try {
        verify_registries();
        verify_round_trip_and_hash();
        verify_driver_config_is_part_of_identity();
        verify_numerical_backend_config_is_part_of_identity();
        verify_strict_validation();
        verify_status_contract();
        verify_practical_v4_contract();
        std::cout << "Helmholtz paper configuration protocol passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_paper_config failed: " << error.what() << '\n';
        return 1;
    }
}

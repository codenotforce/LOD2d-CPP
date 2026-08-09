#include "helmholtz/adaptive/numerical_backend.h"
#include "helmholtz/model.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

using namespace lod2d;
using namespace lod2d::helmholtz;
using namespace lod2d::helmholtz::adaptive;

namespace {

std::string constant_set_content_hash() {
    constexpr std::string_view canonical =
        "C_a:upper:2;C_app:upper:0.1;C_loc:upper:2;C_ov:upper:2;"
        "C_sd:upper:2;C_st:upper:2;beta:upper:0.5;s:exact:1";
    std::uint64_t hash = 14695981039346656037ULL;
    for (const unsigned char byte : canonical) {
        hash ^= byte;
        hash *= 1099511628211ULL;
    }
    std::ostringstream out;
    out << "fnv1a64:" << std::hex << std::setfill('0') << std::setw(16)
        << hash;
    return out.str();
}

struct Options {
    CertifiedEvidencePolicy policy = CertifiedEvidencePolicy::RequireVerified;
    bool check = false;
};

Options parse_options(int argc, char **argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string_view argument(argv[i]);
        if (argument == "--help") {
            std::cout
                << "Usage: bench_helmholtz_certified "
                   "--evidence=strict|conditional [--check]\n";
            std::exit(0);
        }
        if (argument == "--check") {
            options.check = true;
            continue;
        }
        constexpr std::string_view prefix = "--evidence=";
        if (!argument.starts_with(prefix))
            throw std::invalid_argument("unknown argument: " + std::string(argument));
        const std::string_view value = argument.substr(prefix.size());
        if (value == "strict") {
            options.policy = CertifiedEvidencePolicy::RequireVerified;
        } else if (value == "conditional") {
            options.policy = CertifiedEvidencePolicy::AllowConditional;
        } else {
            throw std::invalid_argument(
                "--evidence must be strict or conditional");
        }
    }
    return options;
}

// Source-only definition of paper case R1.  This runner deliberately has no
// exact solution, exact gradient, or evaluation-reference object.
ComplexFunction r1_source(double wavenumber) {
    const auto phi = [](double t) {
        return 16.0 * t * t * (1.0 - t) * (1.0 - t);
    };
    const auto phi_prime = [](double t) {
        return 32.0 * t - 96.0 * t * t + 64.0 * t * t * t;
    };
    const auto phi_second = [](double t) {
        return 32.0 - 192.0 * t + 192.0 * t * t;
    };
    return [=](const Point2 &point) {
        const Complex phase =
            std::exp(Complex(0.0, wavenumber * point.x()));
        return -phase * (
            phi_second(point.x()) * phi(point.y())
            + phi(point.x()) * phi_second(point.y())
            + Complex(0.0, 2.0 * wavenumber)
                * phi_prime(point.x()) * phi(point.y()));
    };
}

CertificateConstantRegistry conditional_constants() {
    CertificateConstantRegistry registry;
    const std::string content_hash = constant_set_content_hash();
    const auto set_upper = [&](std::string name, double value) {
        CertificateConstant constant;
        constant.name = std::move(name);
        constant.value = value;
        constant.direction = CertificateBoundDirection::Upper;
        constant.source =
            "R1 implementation-study floating constant set " + content_hash;
        constant.derivation = "diagnostic value without interval proof";
        constant.mesh_class = "runtime R1 mesh";
        constant.verified = false;
        registry.set(std::move(constant));
    };
    set_upper("C_app", 0.1);
    set_upper("C_st", 2.0);
    set_upper("C_sd", 2.0);
    set_upper("C_ov", 2.0);
    set_upper("C_a", 2.0);
    set_upper("C_loc", 2.0);
    set_upper("beta", 0.5);

    CertificateConstant shift;
    shift.name = "s";
    shift.value = 1.0;
    shift.direction = CertificateBoundDirection::Exact;
    shift.source =
        "R1 implementation-study floating constant set " + content_hash;
    shift.derivation = "diagnostic value without interval proof";
    shift.mesh_class = "runtime R1 mesh";
    shift.verified = false;
    registry.set(std::move(shift));
    return registry;
}

NumericalCertifiedBackendConfig backend_config() {
    NumericalCertifiedBackendConfig config;
    config.problem_id = "paper-R1-k1-production-smoke";
    config.source_id = "R1-polynomial-plane-wave-source-v1";
    config.initial_mesh = make_helmholtz_unit_square_mesh();
    config.wavenumber = 1.0;
    config.source = r1_source(config.wavenumber);
    config.initial_coarse_level = 0;
    config.initial_fine_level = 1;
    config.initial_oversampling = 1;
    config.constants = conditional_constants();
    config.certificate.q0 = 1e6;
    config.coarse_doerfler_theta = 0.5;
    config.audit_doerfler_theta = 0.5;
    config.audit_saturation_factor = 0.05;
    return config;
}

CertifiedDriverConfig driver_config(CertifiedEvidencePolicy policy) {
    CertifiedDriverConfig config;
    config.method = CertifiedMethod::Calod;
    config.error_target = CertifiedErrorTarget::Continuous;
    config.evidence_policy = policy;
    config.mu0 = 0.5;
    config.q0 = 1e6;
    // The implementation smoke asks whether the inexpensive two-level audit
    // is subordinate to the LOD interval; it is not a formal paper threshold.
    config.rho_aud = 0.99;
    config.tolerance = 1e6;
    config.limits.max_state_transitions = 12;
    return config;
}

int run(const Options &options) {
    NumericalCertifiedBackend backend(backend_config());
    CertifiedAdaptiveDriver driver(driver_config(options.policy));
    driver.run(backend);

    const CertifiedDriverCheckpoint &checkpoint = driver.checkpoint();
    if (!driver.terminal() || !driver.termination())
        throw std::runtime_error("driver did not produce a structured terminal outcome");

    const CertifiedTermination &termination = *driver.termination();
    std::cout << "runner_scope=implementation-smoke\n"
              << "wp6_runner=false\n"
              << "claim=conditional/implementation-study\n"
              << "certified=false\n"
              << "case=R1\n"
              << "constant_set_content_hash=" << constant_set_content_hash() << '\n'
              << "evidence_policy=" << to_string(options.policy) << '\n'
              << "driver_claim=" << to_string(termination.claim) << '\n'
              << "terminal_state=" << to_string(termination.state) << '\n'
              << "stop_code=" << to_string(termination.code) << '\n'
              << "transitions=" << checkpoint.transition_count << '\n'
              << "backend_observations=" << backend.observation_count() << '\n'
              << "full_rebuilds=" << backend.full_rebuild_count() << '\n'
              << "detail=" << std::quoted(termination.detail) << '\n';
    for (const CertifiedTransitionRecord &record : checkpoint.history) {
        std::cout << "transition=" << record.sequence
                  << ",before=" << to_string(record.before)
                  << ",action=" << to_string(record.action)
                  << ",after=" << to_string(record.after)
                  << ",claim=" << to_string(record.claim) << '\n';
    }

    const auto has_action = [&](CertifiedDriverAction action) {
        for (const CertifiedTransitionRecord &record : checkpoint.history) {
            if (record.action == action) return true;
        }
        return false;
    };
    if (options.policy == CertifiedEvidencePolicy::RequireVerified) {
        if (termination.code != CertifiedStopCode::UnverifiedEvidence
            || termination.state != CertifiedDriverState::Failure
            || termination.claim != CertifiedEvidenceLevel::Conditional
            || backend.observation_count() != 1) {
            throw std::runtime_error(
                "strict floating-point run did not fail closed with UnverifiedEvidence");
        }
    } else if (termination.claim != CertifiedEvidenceLevel::Conditional) {
        throw std::runtime_error("conditional run was promoted above its evidence");
    }
    if (options.check
        && options.policy == CertifiedEvidencePolicy::AllowConditional) {
        if (termination.state != CertifiedDriverState::Done
            || termination.code
                != CertifiedStopCode::ContinuousToleranceReached
            || backend.observation_count() < 4
            || !has_action(CertifiedDriverAction::AcceptCoarseAdmissibility)
            || !has_action(CertifiedDriverAction::AcceptCorrectorCertificate)
            || !has_action(CertifiedDriverAction::FormPendingCoarseMarking)
            || !has_action(CertifiedDriverAction::CompleteContinuousTolerance)) {
            throw std::runtime_error(
                "conditional check did not complete the real coarse/corrector/error/audit path");
        }
    }
    return 0;
}

} // namespace

int main(int argc, char **argv) {
    try {
        return run(parse_options(argc, argv));
    } catch (const std::exception &error) {
        std::cerr << "bench_helmholtz_certified failed: " << error.what() << '\n';
        return 1;
    }
}

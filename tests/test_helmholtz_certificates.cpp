#include "helmholtz/adaptive/certificates.h"
#include "helmholtz/benchmarks/paper_cases.h"

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
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

void require_close(double first, double second, double tolerance,
                   const char *message) {
    const double scale = std::max({1.0, std::abs(first), std::abs(second)});
    if (std::abs(first - second) > tolerance * scale)
        throw std::runtime_error(message);
}

CertificateConstantRegistry diagnostic_constants() {
    CertificateConstantRegistry registry;
    const auto upper = [&](const std::string &name, double value) {
        registry.set({name, value, CertificateBoundDirection::Upper,
                      "diagnostic test input", "not a rigorous paper constant",
                      "test mesh", "test-policy", false});
    };
    upper("C_app", 0.1);
    upper("C_st", 2.0);
    upper("C_sd", 2.0);
    upper("C_ov", 2.0);
    upper("C_a", 2.0);
    upper("C_loc", 2.0);
    upper("beta", 0.5);
    registry.set({"s", 1.0, CertificateBoundDirection::Exact,
                  "diagnostic test input", "not a rigorous theorem shift",
                  "test mesh", "test-policy", false});
    return registry;
}

CertificateConstantRegistry verified_base_constants(
    const CertificateContextFingerprint &context,
    bool stale_operator_fingerprint = false) {
    CertificateConstantRegistry registry;
    const auto upper = [&](const std::string &name, double value) {
        CertificateConstant constant;
        constant.name = name;
        constant.value = value;
        constant.direction = CertificateBoundDirection::Upper;
        constant.source = "verified negative-test fixture";
        constant.derivation = "direct interval fixture";
        constant.mesh_class = "test hierarchy";
        constant.patch_policy_hash = context.patch_policy;
        constant.verified = true;
        constant.mesh_fingerprint = context.mesh;
        constant.pde_fingerprint = context.pde;
        constant.operator_fingerprint = stale_operator_fingerprint
            && name == "C_sd"
            ? "operators:fnv1a64:stale"
            : context.operators;
        registry.set(std::move(constant));
    };
    upper("C_app", 0.1);
    upper("C_st", 2.0);
    upper("C_sd", 2.0);
    upper("C_ov", 2.0);
    upper("C_a", 2.0);
    // C_loc, beta, and s are intentionally absent: the theorem permits the
    // delta_tot + delta_h localization fallback without those optional data.
    return registry;
}

bool contains_reason(const CorrectorCertificateResult &result,
                     const std::string &needle) {
    return std::any_of(
        result.conditional_reasons.begin(),
        result.conditional_reasons.end(),
        [&](const std::string &reason) {
            return reason.find(needle) != std::string::npos;
        });
}

CorrectorCertificateResult build_case(PaperCase id, int coarse_level,
                                      int fine_level, int ell,
                                      bool verified_algebraic_inputs = false,
                                       bool stale_evidence = false,
                                       bool verified_theorem_constants = false,
                                       bool stale_constant = false,
                                       bool break_conjugation = false) {
    const PaperCaseData data = make_paper_case(id, 2.0);
    AdaptiveMeshHierarchy hierarchy(
        data.initial_mesh, coarse_level, fine_level);
    std::vector<int> all_fine_elements(hierarchy.fine_mesh().elems.size());
    std::iota(all_fine_elements.begin(), all_fine_elements.end(), 0);
    hierarchy.refine_cert_audit_from_fine_elements(all_fine_elements);
    HelmholtzProblemConfig model_config;
    model_config.initial_mesh = data.initial_mesh;
    model_config.H = coarse_level;
    model_config.h = fine_level;
    model_config.ell = ell;
    model_config.wavenumber = data.wavenumber;
    HelmholtzLodModel model = HelmholtzLodModel::build(model_config);
    HelmholtzOperators audit_operators = assemble_helmholtz_operators(
        hierarchy.cert_audit_mesh(), data.wavenumber);
    if (break_conjugation) {
        require(audit_operators.system.rows() >= 2,
                "conjugation negative fixture has too few audit nodes");
        audit_operators.system.coeffRef(0, 1) += Complex(0.25, 0.125);
        audit_operators.system.makeCompressed();
    }
    const KernelPatchPolicy patch_policy = audit_kernel_patch_policy(hierarchy);
    const CertificateContextFingerprint context =
        certificate_context_fingerprint(
            hierarchy, model, audit_operators, patch_policy);
    CorrectorCertificateConfig certificate_config;
    certificate_config.precision_bits = 128;
    certificate_config.cluster_relative_gap = 1e-8;
    CertificateAssemblyEvidence evidence;
    if (verified_algebraic_inputs) {
        evidence.verified = true;
        evidence.energy_entry_radius = 2e-15;
        evidence.system_entry_radius = 2e-15;
        evidence.local_source_entry_radius = 2e-15;
        evidence.prolongation_entry_radius = 2e-15;
        evidence.corrector_entry_radius = 2e-15;
        evidence.constraint_entry_radius = 2e-15;
        evidence.source = "test algebraic input enclosure";
        evidence.hash = "test:wp4-algebraic-enclosure-v1";
        evidence.mesh_fingerprint = context.mesh;
        evidence.pde_fingerprint = context.pde;
        evidence.patch_policy_hash = context.patch_policy;
        evidence.operator_fingerprint = stale_evidence
            ? "operators:fnv1a64:stale"
            : context.operators;
    }
    CertificateConstantRegistry constants = verified_theorem_constants
        ? verified_base_constants(context, stale_constant)
        : diagnostic_constants();
    CorrectorCertificateResult result = build_corrector_certificates(
        hierarchy, model, audit_operators, std::move(constants),
        1.0, AuditKernelResidualEvidence{}, evidence, certificate_config);
    return result;
}

double explicit_generalized_largest(const ComplexMatrix &numerator,
                                    const ComplexMatrix &denominator) {
    Eigen::LLT<ComplexMatrix> factor(denominator);
    require(factor.info() == Eigen::Success,
            "explicit coarse energy Cholesky failed");
    const ComplexMatrix inverse_lower = factor.matrixL().solve(
        ComplexMatrix::Identity(denominator.rows(), denominator.cols()));
    Eigen::SelfAdjointEigenSolver<ComplexMatrix> eigensolver(
        inverse_lower * numerator * inverse_lower.adjoint());
    require(eigensolver.info() == Eigen::Success,
            "explicit generalized eigensolve failed");
    return eigensolver.eigenvalues().maxCoeff();
}

void verify_case(PaperCase id) {
    // Revised R1 has Dirichlet data on both horizontal sides, so its
    // four-vertex level-zero mesh has no free coarse degree of freedom.
    const int coarse_level = id == PaperCase::R1 ? 1 : 0;
    const CorrectorCertificateResult result = build_case(
        id, coarse_level, coarse_level + 1, 1);
    require(result.status == CorrectorCertificateStatus::Conditional,
            "unverified constants or assembly were promoted to certified");
    require(!result.conditional_reasons.empty(),
            "conditional certificate has no reason metadata");
    for (const char *name : {
             "C_app", "C_st", "C_sd", "C_ov", "C_a", "C_Fort",
             "c_W", "C_Pi", "C_loc", "C_ol(ell)",
             "C_ol(ell+s)", "beta", "s"}) {
        const CertificateConstant *constant = result.constants.find(name);
        require(constant != nullptr,
                "required certificate constant is absent from result metadata");
        require(!constant->source.empty() && !constant->mesh_class.empty(),
                "certificate constant lost source or mesh-class metadata");
    }
    require(result.matrices.total_primal.hermitian_relative_error <= 1e-14,
            "G_tot primal is not Hermitian");
    require(result.matrices.total_adjoint.hermitian_relative_error <= 1e-14,
            "G_tot adjoint is not Hermitian");
    require(result.matrices.fine_primal.hermitian_relative_error <= 1e-14,
            "G_h primal is not Hermitian");
    require(result.matrices.fine_adjoint.hermitian_relative_error <= 1e-14,
            "G_h adjoint is not Hermitian");
    require(result.matrices.total_primal.minimum_midpoint_eigenvalue >= -2e-10,
            "G_tot primal is not positive semidefinite within tolerance");
    require(result.matrices.total_adjoint.minimum_midpoint_eigenvalue >= -2e-10,
            "G_tot adjoint is not positive semidefinite within tolerance");
    require(result.matrices.fine_primal.minimum_midpoint_eigenvalue >= -2e-10,
            "G_h primal is not positive semidefinite within tolerance");
    require(result.matrices.fine_adjoint.minimum_midpoint_eigenvalue >= -2e-10,
            "G_h adjoint is not positive semidefinite within tolerance");
    require(result.total_conjugation_relative_error <= 1e-10,
            "formal case failed total primal/adjoint conjugation");
    require(result.fine_conjugation_relative_error <= 1e-10,
            "formal case failed fine primal/adjoint conjugation");
    require(result.conjugation_passed,
            "formal case did not pass the conjugation gate");
    require(result.total_primal_riesz.max_constraint_relative_residual <= 1e-10,
            "total primal Riesz solve violates kernel constraints");
    require(result.total_adjoint_riesz.max_stationarity_relative_residual <= 1e-10,
            "total adjoint Riesz stationarity residual is too large");
    require(result.fine_primal_riesz.max_energy_identity_relative_error <= 1e-10,
            "fine primal Riesz energy identity failed");
    require(result.fine_adjoint_riesz.max_constraint_relative_residual <= 1e-10,
            "fine adjoint Riesz solve violates kernel constraints");

    const double explicit_total = explicit_generalized_largest(
        result.matrices.total_primal.enclosure.midpoint,
        result.matrices.coarse_energy.enclosure.midpoint);
    const double explicit_fine = explicit_generalized_largest(
        result.matrices.fine_primal.enclosure.midpoint,
        result.matrices.coarse_energy.enclosure.midpoint);
    require_close(result.total_primal_spectrum.lambda_max_approximation,
                  explicit_total, 2e-12,
                  "G_tot generalized eigenvalue disagrees with explicit dense result");
    require_close(result.fine_primal_spectrum.lambda_max_approximation,
                  explicit_fine, 2e-12,
                  "G_h generalized eigenvalue disagrees with explicit dense result");
    const double allocated = std::accumulate(
        result.eta_h_element_squared.begin(),
        result.eta_h_element_squared.end(), 0.0);
    require_close(allocated,
                  result.fine_spectrum.lambda_max_approximation,
                  2e-12,
                  "cluster-aware eta_h,T allocation is not conservative");
    require(result.eta_h_allocation_relative_error <= 2e-12,
            "reported cluster allocation residual is too large");
}

void verify_resolution_sensitivity() {
    const CorrectorCertificateResult coarse_corrector =
        build_case(PaperCase::R2a, 0, 1, 1);
    const CorrectorCertificateResult fine_corrector =
        build_case(PaperCase::R2a, 0, 2, 1);
    std::cout << "Theta_h coarse/fine corrector="
              << coarse_corrector.fine_spectrum.theta_approximation << '/'
              << fine_corrector.fine_spectrum.theta_approximation << '\n';
    require(coarse_corrector.fine_spectrum.theta_approximation
                > fine_corrector.fine_spectrum.theta_approximation,
            "lower corrector fine resolution did not increase Theta_h");
}

void verify_localization_sensitivity() {
    const CorrectorCertificateResult short_patch =
        build_case(PaperCase::R2a, 1, 2, 0);
    const CorrectorCertificateResult long_patch =
        build_case(PaperCase::R2a, 1, 2, 1);
    std::cout << "Theta_tot ell0/ell1="
              << short_patch.total_spectrum.theta_approximation << '/'
              << long_patch.total_spectrum.theta_approximation << '\n';
    require(short_patch.total_spectrum.theta_approximation
                > long_patch.total_spectrum.theta_approximation,
            "smaller ell did not expose a larger total-corrector defect");
}

void verify_registry_direction_gate() {
    CertificateConstantRegistry registry = diagnostic_constants();
    registry.set({"C_sd", 2.0, CertificateBoundDirection::Lower,
                  "wrong-direction test", "intentional",
                  "test mesh", "test-policy", true});
    require(!registry.has_verified("C_sd", CertificateBoundDirection::Upper),
            "lower constant was accepted where an upper bound is required");
    require(!registry.missing_or_unverified_required().empty(),
            "incomplete constants registry passed the required-input gate");
}

void verify_plain_floating_inputs_fail_closed() {
    const CorrectorCertificateResult result = build_case(
        PaperCase::R2a, 0, 1, 1,
        false, false, true, false);
    require(result.status != CorrectorCertificateStatus::Certified,
            "ordinary floating-point matrices were promoted to Certified");
    require(!result.matrices.total_primal.enclosure.entries_verified,
            "ordinary floating-point corrector received a verified enclosure");
    require(!result.matrices.fine_primal.enclosure.entries_verified,
            "ordinary floating-point local corrector received a verified enclosure");
    require(contains_reason(result, "assembly enclosure is incomplete"),
            "fail-closed assembly reason is missing");
}

void verify_stale_context_rejected() {
    const CorrectorCertificateResult bound = build_case(
        PaperCase::R2a, 0, 1, 1,
        true, false, true, false);
    require(bound.assembly_evidence.valid_for(bound.context_fingerprint),
            "current-context evidence fixture did not bind");
    CertificateAssemblyEvidence mismatch = bound.assembly_evidence;
    mismatch.mesh_fingerprint += ":stale";
    require(!mismatch.valid_for(bound.context_fingerprint),
            "wrong mesh fingerprint matched the current context");
    mismatch = bound.assembly_evidence;
    mismatch.pde_fingerprint += ":stale";
    require(!mismatch.valid_for(bound.context_fingerprint),
            "wrong PDE fingerprint matched the current context");
    mismatch = bound.assembly_evidence;
    mismatch.patch_policy_hash += ":stale";
    require(!mismatch.valid_for(bound.context_fingerprint),
            "wrong patch-policy fingerprint matched the current context");
    mismatch = bound.assembly_evidence;
    mismatch.operator_fingerprint += ":stale";
    require(!mismatch.valid_for(bound.context_fingerprint),
            "wrong operator fingerprint matched the current context");

    const CorrectorCertificateResult stale_evidence = build_case(
        PaperCase::R2a, 0, 1, 1,
        true, true, true, false);
    require(stale_evidence.assembly_evidence.valid(),
            "stale evidence fixture is structurally invalid");
    require(!stale_evidence.assembly_evidence.valid_for(
                stale_evidence.context_fingerprint),
            "wrong operator fingerprint matched the current context");
    require(stale_evidence.status != CorrectorCertificateStatus::Certified,
            "wrong operator fingerprint was promoted to Certified");
    require(!stale_evidence.matrices.total_primal.enclosure.entries_verified,
            "stale evidence entered the verified matrix chain");
    require(contains_reason(stale_evidence, "context fingerprint mismatch"),
            "stale assembly evidence has no mismatch reason");

    const CorrectorCertificateResult stale_constant = build_case(
        PaperCase::R2a, 0, 1, 1,
        true, false, true, true);
    const CertificateConstant *c_sd = stale_constant.constants.find("C_sd");
    require(c_sd != nullptr && !c_sd->verified,
            "stale theorem constant retained its verified flag");
    require(stale_constant.status != CorrectorCertificateStatus::Certified,
            "stale theorem constant was promoted to Certified");
    require(contains_reason(stale_constant,
                            "constant context fingerprint mismatch: C_sd"),
            "stale theorem constant has no mismatch reason");
}

void verify_localization_triangle_fallback() {
    const CorrectorCertificateResult result = build_case(
        PaperCase::R2a, 0, 1, 1,
        true, false, true, false);
    require(!result.localization_decay_bound_used,
            "missing C_loc/beta/s unexpectedly enabled localization decay");
    require_close(result.delta_ell_upper,
                  result.delta_total_upper + result.delta_h_upper,
                  2e-12,
                  "localization triangle fallback used the wrong upper bound");
    for (const char *optional : {"C_loc", "beta", "s", "C_ol(ell+s)"}) {
        require(!contains_reason(
                    result,
                    std::string("missing or unverified constant: ") + optional),
                "optional localization data incorrectly blocked certification");
    }
}

void verify_inconsistent_operator_rejected() {
    bool rejected = false;
    try {
        (void)build_case(
            PaperCase::R2a, 0, 1, 1,
            true, false, true, false, true);
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    require(rejected,
            "same-size audit operator from another PDE context was accepted");
}

void verify_local_interval_correction_gate() {
    const CorrectorCertificateResult result =
        build_case(PaperCase::R2a, 0, 1, 1, true);
    require(!result.matrix_enclosure_arithmetic_verified,
            "ordinary-double matrix enclosure propagation was marked verified");
    require(result.total_primal_riesz.verified_solve_count == 0,
            "unverified matrix arithmetic entered a verified local Riesz solve");
    require(!result.total_spectrum.verified_lambda.metadata.verified,
            "unverified matrix arithmetic entered a verified spectrum");
    require(result.status == CorrectorCertificateStatus::Conditional,
            "unverified matrix/scalar arithmetic was promoted to Certified");
    require(contains_reason(result, "matrix enclosure propagation"),
            "matrix arithmetic fail-closed reason is missing");
    require(contains_reason(result, "scalar certificate formulas"),
            "scalar arithmetic fail-closed reason is missing");
}

} // namespace

int main() {
    try {
        verify_case(PaperCase::R1);
        verify_case(PaperCase::R2a);
        verify_case(PaperCase::R2b);
        verify_case(PaperCase::S);
        verify_resolution_sensitivity();
        verify_localization_sensitivity();
        verify_registry_direction_gate();
        verify_plain_floating_inputs_fail_closed();
        verify_stale_context_rejected();
        verify_localization_triangle_fallback();
        verify_inconsistent_operator_rejected();
        verify_local_interval_correction_gate();
        std::cout << "Corrector, spectrum, conjugation, and conditional gates passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_certificates failed: "
                  << error.what() << '\n';
        return 1;
    }
}

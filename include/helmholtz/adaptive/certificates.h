#pragma once

#include "helmholtz/adaptive/hierarchy.h"
#include "helmholtz/adaptive/kernel_residual.h"
#include "helmholtz/adaptive/verified_spectrum.h"
#include "helmholtz/model.h"

#include <Eigen/Dense>

#include <map>
#include <string>
#include <vector>

namespace lod2d::helmholtz::adaptive {

enum class CertificateBoundDirection {
    Lower,
    Upper,
    Exact
};

struct CertificateConstant {
    std::string name;
    double value = 0.0;
    CertificateBoundDirection direction = CertificateBoundDirection::Upper;
    std::string source;
    std::string derivation;
    std::string mesh_class;
    std::string patch_policy_hash;
    bool verified = false;
    // Verified constants are scoped to one complete certification context.
    // Human-readable mesh_class metadata is not a substitute for these
    // deterministic fingerprints.
    std::string mesh_fingerprint;
    std::string pde_fingerprint;
    std::string operator_fingerprint;
};

struct CertificateContextFingerprint {
    std::string mesh;
    std::string pde;
    std::string patch_policy;
    std::string operators;

    bool complete() const;
};

class CertificateConstantRegistry {
  public:
    void set(CertificateConstant constant);
    const CertificateConstant *find(const std::string &name) const;
    const std::map<std::string, CertificateConstant> &entries() const {
        return entries_;
    }

    bool has_verified(const std::string &name,
                      CertificateBoundDirection required_direction) const;
    bool has_verified(
        const std::string &name,
        CertificateBoundDirection required_direction,
        const CertificateContextFingerprint &context) const;
    std::vector<std::string> missing_or_unverified_required() const;
    std::vector<std::string> missing_or_unverified_required(
        const CertificateContextFingerprint &context) const;

  private:
    std::map<std::string, CertificateConstant> entries_;
};

// Adds c_W, C_Pi, C_Fort, and the exact combinatorial patch overlaps from
// verified primitive constants whenever their prerequisites are available.
// Existing entries are replaced only by entries with the same name.
void derive_certificate_constants(
    CertificateConstantRegistry &registry,
    const AdaptiveMeshHierarchy &hierarchy,
    double wavenumber,
    int ell,
    int localization_shift,
    const KernelPatchPolicy &patch_policy,
    const CertificateContextFingerprint &context);

// Fingerprints the actual coarse/fine/audit meshes, PDE data, patch policy,
// assembled fine/audit operators, and the corrector midpoint used below.
// Evidence from another hierarchy or operator must not be replayed.
CertificateContextFingerprint certificate_context_fingerprint(
    const AdaptiveMeshHierarchy &hierarchy,
    const HelmholtzLodModel &model,
    const HelmholtzOperators &audit_operators,
    const KernelPatchPolicy &patch_policy);

struct CertificateAssemblyEvidence {
    // A true value is an attestation about the supplied radii, not a request
    // to force certification. Every external floating-point input needs a
    // finite positive enclosure radius, provenance, and exact current-context
    // fingerprints before the evidence is accepted.
    bool verified = false;
    double energy_entry_radius = 0.0;
    double system_entry_radius = 0.0;
    double local_source_entry_radius = 0.0;
    double prolongation_entry_radius = 0.0;
    double corrector_entry_radius = 0.0;
    double constraint_entry_radius = 0.0;
    std::string source;
    std::string hash;
    std::string mesh_fingerprint;
    std::string pde_fingerprint;
    std::string patch_policy_hash;
    std::string operator_fingerprint;

    bool valid() const;
    bool valid_for(const CertificateContextFingerprint &context) const;
};

struct CorrectorCertificateConfig {
    int precision_bits = 128;
    double cluster_relative_gap = 1e-8;
    double cluster_absolute_gap = 1e-12;
    double conjugation_tolerance = 1e-10;
    double q0 = 0.5;
};

struct RieszCertificateDiagnostics {
    int solve_count = 0;
    int verified_solve_count = 0;
    double max_constraint_relative_residual = 0.0;
    double max_stationarity_relative_residual = 0.0;
    double max_energy_identity_relative_error = 0.0;
    double max_solution_error_bound = 0.0;
};

struct CertificateMatrixDiagnostics {
    solver::ComplexMatrixEnclosure enclosure;
    double hermitian_relative_error = 0.0;
    double minimum_midpoint_eigenvalue = 0.0;
};

struct GeneralizedSpectrumCertificate {
    double lambda_max_approximation = 0.0;
    double theta_approximation = 0.0;
    solver::ScalarInterval lambda_enclosure;
    solver::ScalarInterval theta_enclosure;
    solver::VerifiedScalarResult verified_lambda;
    Eigen::VectorXd eigenvalues;
    ComplexMatrix dominant_cluster_basis;
    int dominant_cluster_size = 0;
    double dominant_residual = 0.0;
};

struct CorrectorCertificateMatrices {
    CertificateMatrixDiagnostics coarse_energy;
    CertificateMatrixDiagnostics total_primal;
    CertificateMatrixDiagnostics total_adjoint;
    CertificateMatrixDiagnostics fine_primal;
    CertificateMatrixDiagnostics fine_adjoint;
    std::vector<solver::ComplexMatrixEnclosure> fine_element_primal;
    std::vector<solver::ComplexMatrixEnclosure> fine_element_adjoint;
};

enum class CorrectorCertificateStatus {
    Certified,
    Conditional,
    Invalid
};

struct CorrectorCertificateResult {
    // Step 2 depends only on the corrector/stability chain.  The overall
    // status below additionally includes eta_H evidence, which is not
    // available until Step 3.
    CorrectorCertificateStatus corrector_status =
        CorrectorCertificateStatus::Conditional;
    std::vector<std::string> corrector_conditional_reasons;
    CorrectorCertificateStatus status = CorrectorCertificateStatus::Conditional;
    std::vector<std::string> conditional_reasons;
    CertificateConstantRegistry constants;
    CertificateAssemblyEvidence assembly_evidence;
    CertificateContextFingerprint context_fingerprint;
    // Remains false until enclosure propagation itself uses directed
    // interval operations, not just caller-supplied input radii.
    bool matrix_enclosure_arithmetic_verified = false;
    double eta_H = 0.0;
    AuditKernelResidualEvidence eta_H_evidence;

    KernelPatchPolicy patch_policy;
    CorrectorCertificateMatrices matrices;
    RieszCertificateDiagnostics total_primal_riesz;
    RieszCertificateDiagnostics total_adjoint_riesz;
    RieszCertificateDiagnostics fine_primal_riesz;
    RieszCertificateDiagnostics fine_adjoint_riesz;

    GeneralizedSpectrumCertificate total_primal_spectrum;
    GeneralizedSpectrumCertificate total_adjoint_spectrum;
    GeneralizedSpectrumCertificate fine_primal_spectrum;
    GeneralizedSpectrumCertificate fine_adjoint_spectrum;
    GeneralizedSpectrumCertificate total_spectrum;
    GeneralizedSpectrumCertificate fine_spectrum;

    solver::VerifiedScalarResult audit_infsup;
    double gamma_audit_approximation = 0.0;
    double gamma_audit_lower = 0.0;

    double total_conjugation_relative_error = 0.0;
    double fine_conjugation_relative_error = 0.0;
    bool conjugation_passed = false;
    bool used_independent_worse_side = false;

    double mu = 0.0;
    double theta_total_lower = 0.0;
    double theta_total_upper = 0.0;
    double theta_h_upper = 0.0;
    double delta_total_lower = 0.0;
    double delta_total_upper = 0.0;
    double delta_h_upper = 0.0;
    double delta_ell_lower = 0.0;
    double delta_ell_upper = 0.0;
    bool localization_decay_bound_used = false;
    double q_total = 0.0;
    double q_h = 0.0;
    double q_ell = 0.0;
    bool q_acceptable = false;
    double stability_margin = 0.0;
    bool stability_verified = false;
    // False until every scalar combination below is evaluated with a
    // direction-preserving interval backend.  Matrix/spectrum enclosures
    // alone are insufficient for a Certified final bound.
    bool scalar_formula_enclosures_verified = false;
    double lod_error_lower = 0.0;
    double lod_error_upper = 0.0;

    std::vector<double> eta_h_element_squared;
    double eta_h_allocation_relative_error = 0.0;
    std::string fine_marker_side;

    solver::VerificationMetadata verification_metadata;
};

CorrectorCertificateResult build_corrector_certificates(
    const AdaptiveMeshHierarchy &hierarchy,
    const HelmholtzLodModel &model,
    const HelmholtzOperators &audit_operators,
    CertificateConstantRegistry constants,
    double eta_H,
    const AuditKernelResidualEvidence &eta_H_evidence,
    const CertificateAssemblyEvidence &assembly_evidence = {},
    const CorrectorCertificateConfig &config = {});

const char *certificate_bound_direction_name(CertificateBoundDirection direction);
const char *corrector_certificate_status_name(CorrectorCertificateStatus status);

} // namespace lod2d::helmholtz::adaptive

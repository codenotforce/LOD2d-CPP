#pragma once

#include "helmholtz/adaptive/certified_driver.h"
#include "helmholtz/adaptive/error_control.h"
#include "helmholtz/model.h"

#include <cstddef>
#include <memory>
#include <string>

namespace lod2d::helmholtz::adaptive {

// Configuration of the production, floating-point implementation of the
// paper driver backend.  The source_id is mandatory because std::function is
// opaque: it binds the supplied source callback to the deterministic problem
// identity without exposing an exact or evaluation-reference solution.
struct NumericalCertifiedBackendConfig {
    std::string problem_id;
    std::string source_id;
    TriMesh initial_mesh;
    ComplexFunction source;

    int initial_coarse_level = 0;
    int initial_fine_level = 1;
    int initial_oversampling = 1;
    double wavenumber = 2.0;
    double boundary_beta = 1.0;
    HelmholtzPetrovMode petrov_mode = HelmholtzPetrovMode::TwoSided;
    HelmholtzPatchSolverConfig patch_solver;
    QuadraturePolicy quadrature;
    QuadratureContext quadrature_context;

    CertificateConstantRegistry constants;
    CorrectorCertificateConfig certificate;
    KernelRieszSolver kernel_riesz_solver = KernelRieszSolver::SaddlePoint;

    double coarse_doerfler_theta = 0.5;
    double audit_doerfler_theta = 0.5;
    double audit_saturation_factor = 0.05;
};

// A real numerical backend for CertifiedAdaptiveDriver.  It wires the
// adaptive hierarchy, LOD model, WP3 residual estimator, WP4 certificate
// builder, and certification-audit solve together.  Since its finite-element
// inputs and solves are ordinary floating-point values, every observation is
// deliberately Conditional.  The driver therefore rejects this backend with
// UnverifiedEvidence under RequireVerified instead of promoting diagnostics
// to a Certified result.
class NumericalCertifiedBackend final : public CertifiedDriverBackend {
public:
    explicit NumericalCertifiedBackend(NumericalCertifiedBackendConfig config);
    ~NumericalCertifiedBackend() override;
    NumericalCertifiedBackend(NumericalCertifiedBackend &&) noexcept;
    NumericalCertifiedBackend &operator=(NumericalCertifiedBackend &&) noexcept;
    NumericalCertifiedBackend(const NumericalCertifiedBackend &) = delete;
    NumericalCertifiedBackend &operator=(const NumericalCertifiedBackend &) = delete;

    CertifiedWorkSnapshot work_snapshot() const override;
    CoarseAdmissibilityObservation inspect_coarse_admissibility() override;
    void refine_coarse(const std::vector<int> &marked_elements) override;
    CorrectorCertificationObservation
    inspect_corrector_certification() override;
    void refine_corrector_patches(
        const std::vector<int> &marked_coarse_sources) override;
    void increase_oversampling() override;
    CoarseErrorObservation solve_and_estimate_coarse_error() override;
    AuditControlObservation inspect_audit_control() override;
    void refine_audit(const std::vector<int> &marked_fine_elements) override;

    // Diagnostics useful to runners and integration tests.  They expose no
    // exact/reference data and do not affect checkpoint fingerprints.
    std::size_t full_rebuild_count() const;
    std::size_t observation_count() const;
    const AdaptiveMeshHierarchy &hierarchy() const;
    const CertificateContextFingerprint &certificate_context() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace lod2d::helmholtz::adaptive

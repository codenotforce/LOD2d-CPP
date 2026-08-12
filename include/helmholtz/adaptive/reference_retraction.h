#pragma once

#include "helmholtz/adaptive/kernel_residual.h"

#include <Eigen/Sparse>

#include <cstddef>
#include <vector>

namespace lod2d::helmholtz::adaptive {

struct ReferenceRetractionDiagnostics {
    double projector_identity_relative_error = 0.0;
    double reference_kernel_identity_relative_error = 0.0;
    double kernel_constraint_relative_error = 0.0;
    int maximum_column_nonzeros = 0;
};

// Matrix realization of
//   R_ref = (I_ref - P_H,ref I_H,ref) J_ref | W_amb,
// where J_ref is an element-local L2-dual (Scott--Zhang type) projector from
// the nested ambient mesh onto the fixed reference mesh.  The stored matrix
// maps ambient coefficients to reference coefficients and may safely be
// applied to arbitrary ambient vectors; its range is always contained in
// W_ref and has homogeneous Dirichlet trace.
class ReferenceRetraction {
public:
    const Eigen::SparseMatrix<double> &projector() const {
        return projector_;
    }
    const Eigen::SparseMatrix<double> &matrix() const {
        return retraction_;
    }
    const ReferenceRetractionDiagnostics &diagnostics() const {
        return diagnostics_;
    }

    ComplexVector apply(const ComplexVector &ambient_values) const;
    ComplexMatrix apply(const ComplexMatrix &ambient_values) const;

private:
    friend ReferenceRetraction build_reference_retraction(
        const ReferenceEpochHierarchy &hierarchy);

    Eigen::SparseMatrix<double> projector_;
    Eigen::SparseMatrix<double> retraction_;
    ReferenceRetractionDiagnostics diagnostics_;
};

ReferenceRetraction build_reference_retraction(
    const ReferenceEpochHierarchy &hierarchy);

struct LocalizationEigenConfig {
    // Large production spectra can have a small dominant eigengap.  The
    // iteration is matrix-vector only after whitening, so prefer additional
    // iterations to an O(n^3) dense fallback once the coarse dimension grows.
    int maximum_iterations = 1000;
    double relative_tolerance = 1e-11;
    int dense_cross_check_max_dimension = 64;
    int dense_fallback_max_dimension = 1024;
    ComplexVector warm_start;
};

struct LocalizationSpectrum {
    double lambda_max = 0.0;
    ComplexVector dominant_vector;
    int iterations = 0;
    double relative_residual = 0.0;
    bool converged = false;
    bool dense_cross_checked = false;
    bool used_dense_fallback = false;
    double dense_lambda_max = 0.0;
    double dense_relative_difference = 0.0;
};

enum class LocalizationCertificateStatus {
    ImplementationStudy,
};

// Practical ambient-to-reference localization certificate for one localized
// adjoint corrected basis.  localized_adjoint_basis has one column for every
// entry in coarse_basis_nodes (normally the free coarse nodes).
struct ReferenceLocalizationCertificate {
    LocalizationCertificateStatus status =
        LocalizationCertificateStatus::ImplementationStudy;
    ReferenceRetraction retraction;
    ComplexMatrix defect_rhs;
    AmbientDefectRiesz ambient_riesz;
    ComplexMatrix coarse_energy;
    LocalizationSpectrum spectrum;
    double theta_loc = 0.0;
};

ReferenceLocalizationCertificate compute_reference_localization_certificate(
    const ReferenceEpochHierarchy &hierarchy,
    const HelmholtzOperators &reference_operators,
    const HelmholtzOperators &ambient_operators,
    const ComplexSparseMatrix &localized_adjoint_basis,
    const std::vector<int> &coarse_basis_nodes,
    KernelRieszSolver riesz_solver = KernelRieszSolver::SaddlePoint,
    const LocalizationEigenConfig &eigen_config = {});

// Dense E0/G3 diagnostic only.  It compares the certificate with the direct
// reference-space ideal/localized adjoint-corrector difference and computes
// finite-dimensional ambient kernel constants for the theorem's upper
// direction.  Production paths must respect maximum_free_dofs.
struct SmallMatrixLocalizationValidation {
    int ambient_kernel_dimension = 0;
    int local_kernel_columns = 0;
    double projector_stability_constant = 0.0;
    double retraction_stability_constant = 0.0;
    double ambient_kernel_coercivity = 0.0;
    double stable_decomposition_constant = 0.0;
    double decomposition_relative_residual = 0.0;
    double direct_delta = 0.0;
    double upper_certificate = 0.0;
    bool one_sided_control_holds = false;
};

// Lightweight E0 diagnostic for an ell sweep.  Unlike the full validation
// below, this computes only the direct ideal/localized reference-space
// perturbation and therefore does not repeat hierarchy-only dense constants.
double compute_reference_localization_direct_delta(
    const ReferenceEpochHierarchy &hierarchy,
    const HelmholtzOperators &reference_operators,
    const ComplexSparseMatrix &localized_adjoint_basis,
    const ComplexSparseMatrix &ideal_adjoint_basis,
    const std::vector<int> &coarse_basis_nodes,
    const ReferenceLocalizationCertificate &certificate);

SmallMatrixLocalizationValidation
validate_reference_localization_certificate_small_matrix(
    const ReferenceEpochHierarchy &hierarchy,
    const HelmholtzOperators &reference_operators,
    const HelmholtzOperators &ambient_operators,
    const ComplexSparseMatrix &localized_adjoint_basis,
    const ComplexSparseMatrix &ideal_adjoint_basis,
    const std::vector<int> &coarse_basis_nodes,
    const ReferenceLocalizationCertificate &certificate,
    int maximum_free_dofs = 512);

} // namespace lod2d::helmholtz::adaptive

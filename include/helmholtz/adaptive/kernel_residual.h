#pragma once

#include "helmholtz/adaptive/hierarchy.h"
#include "helmholtz/operators.h"

#include <Eigen/Dense>

#include <cstddef>
#include <string>
#include <vector>

namespace lod2d::helmholtz::adaptive {

enum class KernelRieszSolver {
    SaddlePoint,
    KernelBasisReference
};

// The two practical-adaptive uses of the same constrained local Riesz
// machinery.  The distinction is part of the public result type: an ambient
// defect computation is never an eta_H estimator.
enum class KernelRieszSpace {
    ReferenceResidual,
    AmbientDefect,
};

enum class KernelResidualEvidenceLevel {
    Diagnostic,
    Verified
};

struct AuditKernelResidualEstimate;

// Provenance for eta_H.  The verification state is deliberately read-only:
// callers can copy evidence produced by the estimator, but cannot construct a
// "verified" token or promote an ordinary floating-point result by setting a
// bool.  A default-constructed object is an uninitialized diagnostic token.
class AuditKernelResidualEvidence {
public:
    AuditKernelResidualEvidence() = default;

    KernelResidualEvidenceLevel level() const noexcept { return level_; }
    bool verified() const noexcept {
        return level_ == KernelResidualEvidenceLevel::Verified;
    }
    const std::string &source() const noexcept { return source_; }
    const std::string &backend() const noexcept { return backend_; }
    const std::string &context_fingerprint() const noexcept {
        return context_fingerprint_;
    }
    const std::string &diagnostic_fingerprint() const noexcept {
        return diagnostic_fingerprint_;
    }
    const std::string &result_fingerprint() const noexcept {
        return result_fingerprint_;
    }
    const std::string &failure_reason() const noexcept {
        return failure_reason_;
    }
    // A token is valid only for the exact estimator output from which it was
    // produced.  This prevents a copied token from being attached to a
    // caller-modified eta or element allocation.
    bool matches_eta(double eta) const noexcept;
    bool matches_result(const AuditKernelResidualEstimate &estimate) const;

private:
    friend AuditKernelResidualEstimate estimate_audit_kernel_residual(
        const AdaptiveMeshHierarchy &hierarchy,
        const HelmholtzOperators &audit_operators,
        const ComplexVector &audit_load,
        const ComplexVector &candidate_on_audit,
        double doerfler_theta,
        KernelRieszSolver solver);

    KernelResidualEvidenceLevel level_ =
        KernelResidualEvidenceLevel::Diagnostic;
    std::string source_ = "uninitialized audit-kernel residual evidence";
    std::string backend_;
    std::string context_fingerprint_;
    std::string diagnostic_fingerprint_;
    std::string result_fingerprint_;
    double eta_ = 0.0;
    std::string failure_reason_ = "the audit-kernel estimator has not run";
};

struct KernelPatchPolicy {
    int interpolation_support_layers = 0;
    int patch_layers = 1;
    int enlargement_layers = 1;
    std::string adjacency;
    std::string audit_enlargement;
    std::string hash;
};

struct KernelRieszPatch {
    int coarse_node = -1;
    std::vector<int> coarse_hat_support;
    std::vector<int> interpolation_support;
    std::vector<int> coarse_elements;
    std::vector<int> enlarged_coarse_elements;
    std::vector<int> discrete_dofs;
    std::vector<int> enlarged_discrete_elements;
    std::vector<int> active_constraint_rows;
    std::vector<int> independent_constraint_rows;
    Eigen::MatrixXd constraints;
};

// Compatibility name for the historical certified/audit implementation.
using AuditKernelPatch = KernelRieszPatch;

struct LocalKernelRieszResult {
    double eta = 0.0;
    double eta_squared = 0.0;
    double constraint_relative_residual = 0.0;
    double riesz_relative_residual = 0.0;
    double energy_identity_relative_error = 0.0;
    int kernel_dimension = 0;
    int saddle_unknowns = 0;
    ComplexVector local_values;
    ComplexVector lagrange_multipliers;
};

struct AuditKernelResidualEstimate {
    AuditKernelResidualEvidence evidence;
    KernelPatchPolicy policy;
    KernelRieszSolver solver = KernelRieszSolver::SaddlePoint;
    ComplexVector global_residual;
    std::vector<AuditKernelPatch> patches;
    std::vector<LocalKernelRieszResult> local_results;
    std::vector<double> node_eta;
    std::vector<double> node_eta_squared;
    std::vector<double> element_eta_squared;
    std::vector<int> marked_elements;
    double eta = 0.0;
    double allocation_relative_error = 0.0;
};

// Reference-space residual estimator used for H marking in one fixed
// reference epoch.  All vectors and operators live on reference_mesh().
struct ReferenceResidualRiesz {
    KernelRieszSpace space = KernelRieszSpace::ReferenceResidual;
    KernelPatchPolicy policy;
    KernelRieszSolver solver = KernelRieszSolver::SaddlePoint;
    ComplexVector global_residual;
    std::vector<KernelRieszPatch> patches;
    std::vector<LocalKernelRieszResult> local_results;
    std::vector<double> node_eta;
    std::vector<double> node_eta_squared;
    std::vector<double> element_eta_squared;
    std::vector<int> marked_elements;
    double eta = 0.0;
    double local_square_sum_relative_error = 0.0;
    double allocation_relative_error = 0.0;
};

struct AmbientDefectLocalRiesz {
    // Global defect-RHS column indices represented by columns/gram.  The
    // entries are strictly increasing, so a compact local Gram can be
    // scattered into the global matrix deterministically.
    std::vector<int> active_columns;
    std::vector<LocalKernelRieszResult> columns;
    ComplexMatrix gram;
    double hermitian_relative_error = 0.0;
};

enum class AmbientDefectDetail {
    SummaryOnly,
    StoreLocalDetails,
};

// Ambient-space operator data needed by WP3 to form G_loc.  defect_rhs has
// one column per selected coarse input basis function (normally the free
// coarse nodes); the returned Gram matrix is the sum of the local b_kappa
// Gram contributions.  It is deliberately not an eta_H result and contains
// no coarse-element marking.
struct AmbientDefectRiesz {
    KernelRieszSpace space = KernelRieszSpace::AmbientDefect;
    KernelPatchPolicy policy;
    KernelRieszSolver solver = KernelRieszSolver::SaddlePoint;
    std::size_t patch_count = 0;
    bool local_details_stored = false;
    // Populated only for StoreLocalDetails.  Production localization keeps
    // neither patch geometry nor local solution vectors after the streamed
    // Gram reduction.
    std::vector<KernelRieszPatch> patches;
    std::vector<AmbientDefectLocalRiesz> local_results;
    ComplexMatrix gram;
    Eigen::VectorXd column_eta_squared;
    std::size_t patch_factorizations = 0;
    std::size_t right_hand_side_solves = 0;
    std::size_t maximum_active_columns = 0;
    int parallel_threads = 1;
    // Wall-clock stage diagnostics for production profiling.  patch_solve
    // includes local assembly/factorization/solves and compact Gram formation;
    // gram_reduction is only the deterministic global scatter.
    double patch_solve_seconds = 0.0;
    double gram_reduction_seconds = 0.0;
    double local_square_sum_relative_error = 0.0;
    // The streamed implementation accumulates every compact patch Gram once
    // in patch order.  This remains zero unless that invariant is violated.
    double gram_accumulation_relative_error = 0.0;
};

KernelPatchPolicy kernel_riesz_patch_policy(
    const ReferenceEpochHierarchy &hierarchy,
    KernelRieszSpace space);

std::vector<KernelRieszPatch> build_kernel_riesz_patches(
    const ReferenceEpochHierarchy &hierarchy,
    KernelRieszSpace space,
    const KernelPatchPolicy &policy);

ReferenceResidualRiesz compute_reference_residual_riesz(
    const ReferenceEpochHierarchy &hierarchy,
    const HelmholtzOperators &reference_operators,
    const ComplexVector &reference_load,
    const ComplexVector &candidate_on_reference,
    double doerfler_theta,
    KernelRieszSolver solver = KernelRieszSolver::SaddlePoint);

AmbientDefectRiesz compute_ambient_defect_riesz(
    const ReferenceEpochHierarchy &hierarchy,
    const HelmholtzOperators &ambient_operators,
    const ComplexSparseMatrix &defect_rhs,
    KernelRieszSolver solver = KernelRieszSolver::SaddlePoint,
    AmbientDefectDetail detail = AmbientDefectDetail::SummaryOnly,
    int maximum_parallel_patch_solves = 0);

// Determine the verified support-propagation radius of I_H and freeze the
// corresponding D_z/D_z^+ construction in a deterministic policy hash.
KernelPatchPolicy audit_kernel_patch_policy(
    const AdaptiveMeshHierarchy &hierarchy);

// Context-only fingerprint shared with WP4.  Unlike the diagnostic
// fingerprint, this excludes the load, candidate, solver, and marking theta,
// so a certificate builder can recompute it from its current hierarchy and
// audit operator and reject replayed eta_H evidence.
std::string audit_kernel_residual_context_fingerprint(
    const AdaptiveMeshHierarchy &hierarchy,
    const HelmholtzOperators &audit_operators,
    const KernelPatchPolicy &policy);

std::vector<AuditKernelPatch> build_audit_kernel_patches(
    const AdaptiveMeshHierarchy &hierarchy,
    const KernelPatchPolicy &policy);

// Every local right-hand side is a direct restriction of the same algebraic
// audit residual load - A*candidate. No broken-residual ownership is used.
AuditKernelResidualEstimate estimate_audit_kernel_residual(
    const AdaptiveMeshHierarchy &hierarchy,
    const HelmholtzOperators &audit_operators,
    const ComplexVector &audit_load,
    const ComplexVector &candidate_on_audit,
    double doerfler_theta,
    KernelRieszSolver solver = KernelRieszSolver::SaddlePoint);

struct LocalKernelEffectivity {
    int coarse_node = -1;
    double eta = 0.0;
    double enlarged_patch_error = 0.0;
    double ratio = 0.0;
    bool valid = false;
};

struct EffectivityDistribution {
    std::size_t count = 0;
    std::size_t excluded_zero_error = 0;
    double minimum = 0.0;
    double first_quartile = 0.0;
    double median = 0.0;
    double third_quartile = 0.0;
    double percentile_90 = 0.0;
    double maximum = 0.0;
};

struct LocalKernelEffectivityReport {
    std::vector<LocalKernelEffectivity> local;
    EffectivityDistribution distribution;
};

LocalKernelEffectivityReport evaluate_local_kernel_effectivity(
    const AdaptiveMeshHierarchy &hierarchy,
    const HelmholtzOperators &audit_operators,
    const AuditKernelResidualEstimate &estimate,
    const ComplexVector &certification_solution,
    const ComplexVector &candidate_on_audit);

// Dense, small-mesh diagnostic. It explicitly constructs the global audit
// kernel and finite-dimensional overlap/stable-decomposition constants, then
// checks eta/C_ov <= ||r||_{W'} <= C_sd*eta in the theorem's directions.
struct ExplicitGlobalKernelComparison {
    int global_kernel_dimension = 0;
    int local_kernel_columns = 0;
    double global_eta = 0.0;
    double localized_eta = 0.0;
    double overlap_constant = 1.0;
    double stable_decomposition_constant = 1.0;
    double lower_bound = 0.0;
    double upper_bound = 0.0;
    double span_relative_residual = 0.0;
    bool lower_direction_holds = false;
    bool upper_direction_holds = false;
};

ExplicitGlobalKernelComparison compare_with_explicit_global_kernel(
    const AdaptiveMeshHierarchy &hierarchy,
    const HelmholtzOperators &audit_operators,
    const AuditKernelResidualEstimate &estimate,
    int maximum_free_dofs = 512);

const char *kernel_riesz_solver_name(KernelRieszSolver solver);

} // namespace lod2d::helmholtz::adaptive

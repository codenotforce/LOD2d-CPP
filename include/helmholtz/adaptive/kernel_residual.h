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

struct KernelPatchPolicy {
    int interpolation_support_layers = 0;
    int patch_layers = 1;
    int enlargement_layers = 1;
    std::string adjacency;
    std::string audit_enlargement;
    std::string hash;
};

struct AuditKernelPatch {
    int coarse_node = -1;
    std::vector<int> coarse_hat_support;
    std::vector<int> interpolation_support;
    std::vector<int> coarse_elements;
    std::vector<int> enlarged_coarse_elements;
    std::vector<int> audit_dofs;
    std::vector<int> enlarged_audit_elements;
    std::vector<int> active_constraint_rows;
    std::vector<int> independent_constraint_rows;
    Eigen::MatrixXd constraints;
};

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

// Determine the verified support-propagation radius of I_H and freeze the
// corresponding D_z/D_z^+ construction in a deterministic policy hash.
KernelPatchPolicy audit_kernel_patch_policy(
    const AdaptiveMeshHierarchy &hierarchy);

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

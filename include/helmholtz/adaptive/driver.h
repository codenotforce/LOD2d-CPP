#pragma once

// Legacy strong-residual H-only implementation.  Its output identity is
// intentionally distinct from the paper HLOD and CALOD state-machine driver.
#include "helmholtz/adaptive/estimator.h"
#include "helmholtz/adaptive/reliability.h"

#include <string>
#include <vector>

namespace lod2d::helmholtz::adaptive {

struct AdaptiveHelmholtzConfig {
    HelmholtzProblemConfig problem;
    int initial_coarse_level = 1;
    int fine_level = 6;
    int max_iterations = 6;
    int max_coarse_dofs = 20000;
    int dual_patch_layers = 1;
    int stability_max_dofs = 512;
    double theta = 0.5;
    double tolerance = 0.0;
    double q_limit = 0.5;
    diagnostics::ResidualEstimatorKind estimator =
        diagnostics::ResidualEstimatorKind::Mixed;
    bool compute_dual_calibration = true;
};

struct AdaptiveIterationRecord {
    int iteration = 0;
    int coarse_nodes = 0;
    int coarse_elements = 0;
    int fine_nodes = 0;
    int fine_elements = 0;
    int min_coarse_level = 0;
    int max_coarse_level = 0;
    int marked_elements = 0;
    int closure_added_elements = 0;
    double H_max = 0.0;
    double q_max = 0.0;
    double q_effective = 0.0;
    double energy_error = 0.0;
    double relative_energy_error = 0.0;
    double l2_error = 0.0;
    double relative_l2_error = 0.0;
    double exact_energy_error = 0.0;
    double exact_l2_error = 0.0;
    double fine_exact_energy_error = 0.0;
    double fine_exact_l2_error = 0.0;
    double eta_fine = 0.0;
    double eta_mixed = 0.0;
    double eta_macro = 0.0;
    double selected_estimator = 0.0;
    double selected_effectivity = 0.0;
    double exact_effectivity = 0.0;
    double residual_identity_error = 0.0;
    double dual_spearman = 0.0;
    double dual_marking_overlap = 0.0;
    double inf_sup = 0.0;
    double petrov_residual = 0.0;
    double corrector_residual = 0.0;
    double constraint_residual = 0.0;
    double build_ms = 0.0;
    double solve_ms = 0.0;
    double reference_ms = 0.0;
    double estimate_ms = 0.0;
    double dual_ms = 0.0;
    double mark_refine_ms = 0.0;
    double total_ms = 0.0;
};

struct AdaptiveHelmholtzResult {
    std::string output_namespace = "helmholtz/hlod_proxy";
    std::string implementation_status = "diagnostic_h_only_proxy";
    std::vector<AdaptiveIterationRecord> history;
    TriMesh final_coarse_mesh;
    std::vector<int> final_coarse_levels;
    std::vector<std::uint64_t> final_coarse_element_ids;
    std::vector<double> final_indicator_squared;
    std::string stop_reason;
};

// This function remains the diagnostic HLOD-proxy entry point.  Paper HLOD
// and CALOD use CertifiedAdaptiveDriver from certified_driver.h.
AdaptiveHelmholtzResult run_adaptive_helmholtz(
    const AdaptiveHelmholtzConfig &config,
    const ComplexFunction &source,
    const ComplexFunction &exact = {},
    const ComplexGradientFunction &exact_gradient = {});

const char *residual_estimator_name(diagnostics::ResidualEstimatorKind kind);

} // namespace lod2d::helmholtz::adaptive

#pragma once

#include "helmholtz/operators.h"

#include <Eigen/Dense>

#include <vector>

namespace lod2d::helmholtz::adaptive {

struct CandidateFluxConfig {
    double doerfler_theta = 0.5;
    bool closure_cost_aware_marking = false;
    std::size_t closure_cost_candidate_pool_factor = 2;
    QuadraturePolicy quadrature;
    QuadratureContext quadrature_context;
    // Empty means the theorem-level global reconstruction.  A nonempty list
    // activates only vertex patches incident to these elements and is a
    // marking heuristic; global reliability is then deliberately not claimed.
    std::vector<int> active_elements;
    bool compute_discrete_residual_audit = false;
};

struct CandidateFluxPatchDiagnostics {
    int vertex = -1;
    std::vector<int> elements;
    bool touches_dirichlet = false;
    Complex compatibility_defect = Complex(0.0, 0.0);
    Complex compatibility_correction = Complex(0.0, 0.0);
    double compatibility_identity_error = 0.0;
    double corrected_compatibility_error = 0.0;
    double maximum_divergence_residual = 0.0;
    double maximum_constrained_flux_residual = 0.0;
};

struct CandidateFluxResult {
    // On element K, sigma_c(x) = (a_x + b*x, a_y + b*y).
    std::vector<Eigen::Vector3cd> element_rt0_coefficients;
    std::vector<Complex> projected_source;
    std::vector<Complex> delta_pg;
    std::vector<double> element_eta_squared;
    std::vector<int> marked_elements;
    std::vector<CandidateFluxPatchDiagnostics> patches;
    double eta_eq = 0.0;
    double maximum_patch_compatibility_error = 0.0;
    double maximum_element_divergence_residual = 0.0;
    double maximum_boundary_flux_residual = 0.0;
    double doerfler_relative_error = 0.0;
    double discrete_residual_dual_norm = 0.0;
    bool discrete_residual_audit_performed = false;
    bool global_reconstruction = true;
    double time_marking = 0.0;
    std::size_t marking_candidate_pool = 0;
    std::size_t estimated_selected_closure_cost = 0;
};

// Coefficients in the reference-triangle basis
//   [P2]^2 + (xi,eta) span{xi^2,xi eta,eta^2}.
// The physical field is obtained with the contravariant Piola transform.
using CandidateRT2Coefficients = Eigen::Matrix<Complex, 15, 1>;

struct CandidateFluxRT2Result {
    std::vector<CandidateRT2Coefficients> element_rt2_coefficients;
    // P1 L2 projection of f, stored as its three local nodal values.
    std::vector<Eigen::Vector3cd> projected_source_p1;
    std::vector<Complex> delta_pg;
    std::vector<double> element_eta_squared;
    std::vector<int> marked_elements;
    std::vector<CandidateFluxPatchDiagnostics> patches;
    double eta_eq = 0.0;
    double maximum_patch_compatibility_error = 0.0;
    double maximum_element_divergence_residual = 0.0;
    double maximum_normal_continuity_residual = 0.0;
    double maximum_boundary_flux_residual = 0.0;
    double doerfler_relative_error = 0.0;
    double discrete_residual_dual_norm = 0.0;
    bool discrete_residual_audit_performed = false;
    bool global_reconstruction = true;
    double time_prepare = 0.0;
    double time_patch_solve = 0.0;
    double time_deterministic_merge = 0.0;
    double time_estimator_and_audit = 0.0;
    double time_marking = 0.0;
    std::size_t marking_candidate_pool = 0;
    std::size_t estimated_selected_closure_cost = 0;
    int parallel_threads = 1;
};

CandidateFluxResult reconstruct_candidate_flux_rt0(
    const TriMesh &candidate_mesh,
    const HelmholtzOperators &candidate_operators,
    const ComplexFunction &source,
    const ComplexVector &candidate_values,
    const CandidateFluxConfig &config = {});

Eigen::Vector2cd evaluate_candidate_rt0_flux(
    const Eigen::Vector3cd &coefficients,
    const Point2 &point);

// Production candidate estimator for a P1 primal/source projection.  RT0 is
// retained above solely as an explicit projected-data smoke path.
CandidateFluxRT2Result reconstruct_candidate_flux_rt2(
    const TriMesh &candidate_mesh,
    const HelmholtzOperators &candidate_operators,
    const ComplexFunction &source,
    const ComplexVector &candidate_values,
    const CandidateFluxConfig &config = {});

Eigen::Vector2cd evaluate_candidate_rt2_flux(
    const TriMesh &mesh,
    int element,
    const CandidateRT2Coefficients &coefficients,
    const Point2 &point);

Complex evaluate_candidate_rt2_divergence(
    const TriMesh &mesh,
    int element,
    const CandidateRT2Coefficients &coefficients,
    const Point2 &point);

} // namespace lod2d::helmholtz::adaptive

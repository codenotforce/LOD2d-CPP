#pragma once

#include "helmholtz/adaptive/kernel_residual.h"

#include <string>

namespace lod2d::helmholtz::adaptive {

enum class CandidateGapEvidenceMode {
    Practical,
    Certified,
};

struct CandidateDualGapConfig {
    double continuity_constant = 1.0;
    double overlap_constant = 1.0;
    double reference_upper_bound = 0.0;
    double epoch_switch_ratio = 0.5;
    CandidateGapEvidenceMode evidence_mode =
        CandidateGapEvidenceMode::Practical;
    KernelRieszSolver solver = KernelRieszSolver::SaddlePoint;
};

struct CandidateDualGapResult {
    CandidateResidualRiesz riesz;
    CandidateGapEvidenceMode evidence_mode =
        CandidateGapEvidenceMode::Practical;
    std::string gap_field_name = "L_gap_practical_c";
    double eta_dual_c = 0.0;
    double L_c = 0.0;
    double reference_upper_bound = 0.0;
    double L_gap_c = 0.0;
    double epoch_switch_threshold = 0.0;
    bool epoch_switch_requested = false;
};

// Production entry point.  It only restricts the already available residual
// to W_c; it never forms or solves the candidate Helmholtz Galerkin problem.
CandidateDualGapResult build_candidate_dual_gap(
    const ReferenceEpochHierarchy &hierarchy,
    const HelmholtzOperators &candidate_operators,
    const ComplexVector &candidate_load,
    const ComplexVector &lod_on_candidate,
    const CandidateDualGapConfig &config);

struct CandidateDualGapValidation {
    double candidate_error = 0.0;
    double reference_error = 0.0;
    double reference_candidate_gap = 0.0;
    bool L_c_lower_bound_holds = false;
    bool L_gap_lower_bound_holds = false;
};

// E0-only audit.  The caller deliberately supplies u_c; this API is kept
// separate so production drivers cannot accidentally solve for it.
CandidateDualGapValidation validate_candidate_dual_gap_small_mesh(
    const HelmholtzOperators &candidate_operators,
    const ComplexVector &candidate_solution,
    const ComplexVector &reference_solution_on_candidate,
    const ComplexVector &lod_on_candidate,
    const CandidateDualGapResult &gap,
    double tolerance = 1e-10);

} // namespace lod2d::helmholtz::adaptive

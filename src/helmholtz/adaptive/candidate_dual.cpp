#include "helmholtz/adaptive/candidate_dual.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace lod2d::helmholtz::adaptive {
namespace {

double energy_norm(
    const HelmholtzOperators &operators,
    const ComplexVector &values) {
    const Eigen::SparseMatrix<double> energy = operators.stiffness
        + operators.wavenumber * operators.wavenumber * operators.mass;
    return std::sqrt(std::max(
        0.0, std::real(values.dot(energy.cast<Complex>() * values))));
}

} // namespace

CandidateDualGapResult build_candidate_dual_gap(
    const ReferenceEpochHierarchy &hierarchy,
    const HelmholtzOperators &candidate_operators,
    const ComplexVector &candidate_load,
    const ComplexVector &lod_on_candidate,
    const CandidateDualGapConfig &config) {
    if (!(config.continuity_constant > 0.0)
        || !(config.overlap_constant > 0.0)
        || !(config.reference_upper_bound >= 0.0)
        || !std::isfinite(config.reference_upper_bound)
        || !(config.epoch_switch_ratio > 0.0)) {
        throw std::invalid_argument("candidate dual-gap constants are invalid");
    }
    CandidateDualGapResult result;
    result.evidence_mode = config.evidence_mode;
    result.gap_field_name = config.evidence_mode
            == CandidateGapEvidenceMode::Certified
        ? "L_gap_c" : "L_gap_practical_c";
    result.riesz = compute_candidate_residual_riesz(
        hierarchy, candidate_operators, candidate_load,
        lod_on_candidate, config.solver);
    result.eta_dual_c = result.riesz.eta;
    result.L_c = result.eta_dual_c
        / (config.continuity_constant * config.overlap_constant);
    result.reference_upper_bound = config.reference_upper_bound;
    result.L_gap_c = std::max(
        0.0, result.L_c - result.reference_upper_bound);
    result.epoch_switch_threshold = config.epoch_switch_ratio
        * result.reference_upper_bound;
    result.epoch_switch_requested =
        result.L_gap_c >= result.epoch_switch_threshold;
    return result;
}

CandidateDualGapValidation validate_candidate_dual_gap_small_mesh(
    const HelmholtzOperators &candidate_operators,
    const ComplexVector &candidate_solution,
    const ComplexVector &reference_solution_on_candidate,
    const ComplexVector &lod_on_candidate,
    const CandidateDualGapResult &gap,
    double tolerance) {
    const int nodes = candidate_operators.system.rows();
    if (candidate_solution.size() != nodes
        || reference_solution_on_candidate.size() != nodes
        || lod_on_candidate.size() != nodes
        || !candidate_solution.allFinite()
        || !reference_solution_on_candidate.allFinite()
        || !lod_on_candidate.allFinite()
        || !(tolerance >= 0.0)) {
        throw std::invalid_argument(
            "candidate dual-gap validation inputs are invalid");
    }
    CandidateDualGapValidation result;
    result.candidate_error = energy_norm(
        candidate_operators, candidate_solution - lod_on_candidate);
    result.reference_error = energy_norm(
        candidate_operators,
        reference_solution_on_candidate - lod_on_candidate);
    result.reference_candidate_gap = energy_norm(
        candidate_operators,
        candidate_solution - reference_solution_on_candidate);
    const double scale = std::max({
        1.0, result.candidate_error, result.reference_candidate_gap,
        gap.L_c, gap.L_gap_c});
    result.L_c_lower_bound_holds =
        gap.L_c <= result.candidate_error + tolerance * scale;
    result.L_gap_lower_bound_holds =
        gap.L_gap_c <= result.reference_candidate_gap + tolerance * scale;
    return result;
}

} // namespace lod2d::helmholtz::adaptive

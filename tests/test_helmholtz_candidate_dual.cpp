#include "helmholtz/adaptive/candidate_dual.h"
#include "helmholtz/benchmarks/paper_cases.h"
#include "helmholtz/operators.h"

#include <iostream>
#include <stdexcept>

using namespace lod2d;
using namespace lod2d::helmholtz;
using namespace lod2d::helmholtz::adaptive;
using namespace lod2d::helmholtz::benchmarks;
using namespace lod2d::helmholtz::experiments;

namespace {

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

} // namespace

int main() {
    try {
        const PaperCaseData data = make_paper_case(PaperCase::R1, 8.0);
        ReferenceEpochHierarchy hierarchy(data.initial_mesh, 1, 2);
        hierarchy.begin_reference_epoch();
        const ReferenceEpochRefinementResult enriched =
            hierarchy.enrich_candidate({0, 1});
        require(enriched.changed(), "candidate enrichment failed");

        const HelmholtzOperators reference_operators =
            assemble_helmholtz_operators(
                hierarchy.reference_mesh(), data.wavenumber);
        const ComplexVector reference_load = assemble_helmholtz_load(
            hierarchy.reference_mesh(), data.source, {},
            data.quadrature_context);
        const ComplexVector reference_solution = solve_helmholtz_fem(
            reference_operators, reference_load);

        const HelmholtzOperators candidate_operators =
            assemble_helmholtz_operators(
                hierarchy.candidate_mesh(), data.wavenumber);
        const ComplexVector candidate_load = assemble_helmholtz_load(
            hierarchy.candidate_mesh(), data.source, {},
            data.quadrature_context);
        const ComplexVector candidate_solution = solve_helmholtz_fem(
            candidate_operators, candidate_load);
        const ComplexVector lod_on_candidate =
            hierarchy.reference_to_candidate().cast<Complex>()
            * reference_solution;

        CandidateDualGapConfig config;
        config.continuity_constant = 8.0;
        config.overlap_constant = 8.0;
        config.reference_upper_bound = 0.0;
        config.epoch_switch_ratio = 0.5;
        config.evidence_mode = CandidateGapEvidenceMode::Certified;
        const CandidateDualGapResult gap = build_candidate_dual_gap(
            hierarchy, candidate_operators, candidate_load,
            lod_on_candidate, config);
        require(gap.riesz.space == KernelRieszSpace::CandidateResidual,
                "candidate dual used the wrong kernel mesh");
        require(gap.eta_dual_c > 0.0 && gap.L_c > 0.0,
                "enriched candidate produced a zero dual gap");
        require(gap.gap_field_name == "L_gap_c",
                "certified gap field was mislabeled");
        require(gap.riesz.marked_elements.empty(),
                "candidate dual-Riesz incorrectly produced marking");
        require(gap.riesz.patch_factorizations
                    + gap.riesz.skipped_zero_kernel_patches
                    == gap.riesz.patches.size()
                    && gap.riesz.parallel_threads >= 1
                    && gap.riesz.prepare_seconds >= 0.0
                    && gap.riesz.patch_solve_seconds >= 0.0
                    && gap.riesz.reduction_seconds >= 0.0,
                "candidate dual phase diagnostics are inconsistent");

        const CandidateDualGapValidation validation =
            validate_candidate_dual_gap_small_mesh(
                candidate_operators, candidate_solution,
                lod_on_candidate, lod_on_candidate, gap);
        require(validation.L_c_lower_bound_holds,
                "L_c is not a lower bound for the candidate error");
        require(validation.L_gap_lower_bound_holds,
                "L_gap is not a lower bound for the reference-candidate gap");

        CandidateDualGapConfig practical = config;
        practical.evidence_mode = CandidateGapEvidenceMode::Practical;
        const CandidateDualGapResult practical_gap =
            build_candidate_dual_gap(
                hierarchy, candidate_operators, candidate_load,
                lod_on_candidate, practical);
        require(practical_gap.gap_field_name == "L_gap_practical_c",
                "practical gap was promoted to a certified field");

        std::cout << "eta_dual_c=" << gap.eta_dual_c
                  << " L_c=" << gap.L_c
                  << " candidate_error=" << validation.candidate_error
                  << " L_gap_c=" << gap.L_gap_c
                  << " reference_candidate_gap="
                  << validation.reference_candidate_gap << '\n';
        std::cout << "Candidate dual-Riesz gap certificate passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_candidate_dual failed: "
                  << error.what() << '\n';
        return 1;
    }
}

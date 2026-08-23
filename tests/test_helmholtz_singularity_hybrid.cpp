#include "helmholtz/adaptive/kernel_residual.h"
#include "helmholtz/adaptive/singularity_hybrid.h"
#include "helmholtz/benchmarks/paper_cases.h"
#include "helmholtz/model.h"

#include <algorithm>
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

HelmholtzProblemConfig model_config(const PaperCaseData &data, int ell) {
    HelmholtzProblemConfig config;
    config.ell = ell;
    config.wavenumber = data.wavenumber;
    config.initial_mesh = data.initial_mesh;
    config.quadrature_context = data.quadrature_context;
    return config;
}

} // namespace

int main() {
    try {
        const PaperCaseData data = make_paper_case(PaperCase::S, 4.0, 0.05);
        ReferenceEpochHierarchy hierarchy(data.initial_mesh, 3, 5);
        hierarchy.begin_reference_epoch();
        const std::vector<Point2> singular_set{Point2(0.0, 0.0)};
        const std::uint64_t reference_version = hierarchy.reference_mesh_version();
        const std::size_t reference_cells = hierarchy.reference_mesh().elems.size();

        const HybridMatchingResult ell1 = restore_hybrid_reference_matching(
            hierarchy, singular_set, 1);
        require(ell1.matching_holds && ell1.reference_unchanged,
                "ell=1 hybrid matching was not restored");
        require(!ell1.regions.omega_s_elements.empty()
                    && !ell1.regions.omega_f_elements.empty(),
                "hybrid singular regions are empty");
        require(!ell1.regions.regular_elements.empty(),
                "small hybrid audit has no regular marking region");

        const HelmholtzLodModel hybrid =
            HelmholtzLodModel::build_adaptive_hybrid(
                model_config(data, 1), hierarchy,
                ell1.regions.omega_s_elements);
        require(hybrid.correctors().diagnostics.skipped_patch_count
                    == static_cast<int>(ell1.regions.omega_s_elements.size()),
                "singular-core corrector skips were not recorded");
        require(hybrid.correctors().diagnostics.skipped_patch_work_units > 0,
                "skipped corrector work dimension was not recorded");
        for (int element : ell1.regions.omega_s_elements)
            require(hybrid.correctors().primal[element].empty(),
                    "a zero-kernel singular corrector was solved");
        const HelmholtzLodSolution solution = hybrid.solve_source(data.source);
        require(solution.fine_values.allFinite(),
                "hybrid Petrov-Galerkin solution is non-finite");

        const ComplexVector load = assemble_helmholtz_load(
            hierarchy.reference_mesh(), data.source, {},
            data.quadrature_context);
        const ReferenceResidualRiesz global = compute_reference_residual_riesz(
            hierarchy, hybrid.operators(), load, solution.fine_values, 0.5);
        require(global.eta >= 0.0 && std::isfinite(global.eta),
                "hybrid global reference certificate is invalid");
        require(global.local_results.size() == global.patches.size(),
                "hybrid reference Riesz patch accounting is incomplete");
        for (std::size_t patch = 0; patch < global.patches.size(); ++patch) {
            if (global.patches[patch].discrete_dofs.empty()) {
                require(global.local_results[patch].local_values.size() == 0,
                        "a zero-dimensional reference Riesz patch was solved");
            }
        }

        std::vector<double> indicators(
            hierarchy.coarse_mesh().elems.size(), 1.0);
        const std::vector<int> marked = mark_hybrid_regular_region(
            indicators, ell1.regions, 0.5);
        require(!marked.empty(), "regular-region Dörfler marking is empty");
        for (int element : marked)
            require(ell1.regions.in_regular[element],
                    "hybrid coarse marking entered Omega_F");

        const std::size_t omega_f_before = ell1.regions.omega_f_elements.size();
        HybridMatchingResult ell2 = restore_hybrid_reference_matching(
            hierarchy, singular_set, 2);
        require(ell2.matching_holds && ell2.reference_unchanged
                    && ell2.regions.omega_f_elements.size() >= omega_f_before,
                "ell increase did not expand and restore the matching region");
        require(hierarchy.reference_mesh_version() == reference_version
                    && hierarchy.reference_mesh().elems.size() == reference_cells,
                "hybrid matching mutated the fixed reference");

        require(!ell2.regions.regular_elements.empty(),
                "ell=2 hybrid audit has no regular cell for the gap guard");
        require(hierarchy.propose_coarse_refinement(
                    {ell2.regions.regular_elements.front()}).changed()
                    && hierarchy.reference_contains_proposed_coarse(),
                "regular hybrid proposal is not contained in the reference");
        const SingularRegionClassification proposed_regions =
            classify_singular_regions(
                hierarchy.proposed_coarse_mesh(), singular_set, 2);
        const int global_gap = hierarchy.minimum_proposed_reference_level_gap();
        const int regular_positive_gap =
            hierarchy.minimum_proposed_reference_level_gap(
                proposed_regions.in_regular);
        require(global_gap == 0
                    && regular_positive_gap > 0,
                "hybrid gap guard did not exclude deliberately matched zero-gap cells");
        require(hierarchy.commit_coarse_refinement().changed(),
                "regular hybrid gap-audit proposal was not committed");
        ell2 = restore_hybrid_reference_matching(hierarchy, singular_set, 2);

        const auto singular_candidate = std::find(
            hierarchy.candidate_parent_coarse_elements().begin(),
            hierarchy.candidate_parent_coarse_elements().end(),
            ell2.regions.seed_elements.front());
        require(singular_candidate
                    != hierarchy.candidate_parent_coarse_elements().end(),
                "singular candidate parent was not found");
        const int singular_candidate_element = static_cast<int>(
            singular_candidate
            - hierarchy.candidate_parent_coarse_elements().begin());
        const std::size_t candidate_before = hierarchy.candidate_mesh().elems.size();
        const ReferenceEpochRefinementResult candidate =
            hierarchy.enrich_candidate({singular_candidate_element});
        require(candidate.changed()
                    && hierarchy.candidate_mesh().elems.size() > candidate_before
                    && hierarchy.reference_mesh_version() == reference_version,
                "candidate could not enrich inside the singular region");

        std::cout << "omega_s=" << ell2.regions.omega_s_elements.size()
                  << " omega_f=" << ell2.regions.omega_f_elements.size()
                  << " skipped_correctors="
                  << hybrid.correctors().diagnostics.skipped_patch_count
                  << " skipped_work_units="
                  << hybrid.correctors().diagnostics.skipped_patch_work_units
                  << " skipped_riesz="
                  << global.skipped_zero_kernel_patches << '\n';
        std::cout << "Singularity-aware hybrid matching passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_singularity_hybrid failed: "
                  << error.what() << '\n';
        return 1;
    }
}

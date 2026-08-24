#include "helmholtz/adaptive/kernel_residual.h"
#include "helmholtz/adaptive/singularity_hybrid.h"
#include "helmholtz/benchmarks/paper_cases.h"
#include "helmholtz/model.h"

#include <algorithm>
#include <cmath>
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

TriMesh make_graph_chain_mesh() {
    TriMesh mesh;
    mesh.nodes = {
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
        {2.0, 0.0}, {2.0, 1.0}, {3.0, 0.0}, {3.0, 1.0}};
    // Consecutive triangles meet at one vertex, so their element-graph
    // distances from the first triangle are exactly 0, 1, and 2.
    mesh.elems = {{0, 1, 2}, {1, 3, 4}, {4, 5, 6}};
    return mesh;
}

std::size_t embedding_children(
    const Eigen::SparseMatrix<double> &embedding, const int parent) {
    std::size_t count = 0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(embedding, parent); it;
         ++it) {
        if (std::abs(it.value()) > 1e-14) ++count;
    }
    return count;
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

        const TriMesh graph_chain = make_graph_chain_mesh();
        const SingularRegionClassification chain_ell1 =
            classify_singular_regions_with_physical_radius(
                graph_chain, {Point2(0.0, 0.0)}, 1, 2.3);
        const SingularRegionClassification chain_ell7 =
            classify_singular_regions_with_physical_radius(
                graph_chain, {Point2(0.0, 0.0)}, 7, 2.3);
        // The last covered cell has zero-based graph distance two, hence it
        // belongs to the manuscript neighborhood N^3(S), not N^2(S).
        require(chain_ell1.l_s == 3
                    && chain_ell7.l_s == chain_ell1.l_s
                    && chain_ell7.l_s < chain_ell7.corrector_ell,
                "physical ell_S is not the minimal ell-independent graph radius");
        require(chain_ell1.omega_s_elements.size() == 3
                    && chain_ell7.omega_s_elements
                        == chain_ell1.omega_s_elements,
                "physical singular core changed with the corrector width");
        require(chain_ell1.covered_physical_radius + 1e-12 >= 2.3,
                "reported physical-core coverage is too small");

        const SingularRegionClassification chain_local =
            classify_singular_regions_with_physical_radius(
                graph_chain, {Point2(0.0, 0.0)}, 0, 0.1);
        const HybridGradedReserveProfile chain_no_collar =
            make_hybrid_graded_reserve_profile(
                graph_chain, chain_local, 4, 0);
        const HybridGradedReserveProfile chain_reserve =
            make_hybrid_graded_reserve_profile(
                graph_chain, chain_local, 4, 1);
        require(chain_no_collar.target_level_gaps
                    == std::vector<int>({0, 1, 2})
                    && chain_reserve.target_level_gaps
                        == std::vector<int>({0, 0, 1})
                    && std::none_of(
                        chain_reserve.at_full_target.begin(),
                        chain_reserve.at_full_target.end(),
                        [](const char selected) { return selected; }),
                "hybrid reserve does not include the conformity collar");

        const SingularRegionClassification physical =
            classify_singular_regions_with_physical_radius(
                hierarchy.coarse_mesh(), singular_set, 1, 0.25);
        require(physical.minimum_physical_radius == 0.25
                    && physical.covered_physical_radius + 1e-12 >= 0.25,
                "physical hybrid neighborhood does not cover the frozen radius");
        const SingularRegionClassification physical_ell3 =
            classify_singular_regions_with_physical_radius(
                hierarchy.coarse_mesh(), singular_set, 3, 0.25);
        require(physical_ell3.l_s == physical.l_s
                    && physical_ell3.omega_s_elements
                        == physical.omega_s_elements
                    && physical_ell3.omega_f_elements.size()
                        >= physical.omega_f_elements.size(),
                "physical core depends on ell or transition did not expand");
        ReferenceEpochHierarchy physical_hierarchy(data.initial_mesh, 2, 4);
        physical_hierarchy.begin_reference_epoch();
        const HybridMatchingResult physical_matching =
            restore_hybrid_reference_matching_with_physical_radius(
                physical_hierarchy, singular_set, 1, 0.25);
        require(physical_matching.matching_holds
                    && physical_matching.reference_unchanged
                    && physical_matching.regions.covered_physical_radius
                        + 1e-12 >= 0.25,
                "fixed-radius hybrid matching invariant was not restored");

        // H3/h7 refresh-closure regression: exact matching gives gap zero on
        // Omega_F, while the first regular graph layer cannot support a
        // uniform positive reserve.  The candidate-only dry run must predict
        // the same spill-free refinement as the production mutation.
        ReferenceEpochHierarchy reserve_hierarchy(data.initial_mesh, 3, 7);
        reserve_hierarchy.begin_reference_epoch();
        const HybridMatchingResult reserve_matching =
            restore_hybrid_reference_matching_with_physical_radius(
                reserve_hierarchy, singular_set, 1, 0.01);
        require(reserve_matching.matching_holds,
                "reserve MRE could not establish exact matching");
        reserve_hierarchy.propose_identity_coarse();
        const SingularRegionClassification reserve_regions =
            classify_singular_regions_with_physical_radius(
                reserve_hierarchy.proposed_coarse_mesh(), singular_set, 1,
                0.01);
        constexpr int reserve_target = 5;
        const HybridGradedReserveProfile reserve_profile =
            make_hybrid_graded_reserve_profile(
                reserve_hierarchy.proposed_coarse_mesh(), reserve_regions,
                reserve_target, 0);
        const ReferenceEpochCandidateDeepeningProbe reserve_probe =
            reserve_hierarchy.probe_candidate_deepening_over_proposed_coarse(
                reserve_profile.target_level_gaps,
                reserve_regions.in_omega_f,
                reserve_profile.at_full_target);
        const int uniform_regular_gap =
            reserve_hierarchy.minimum_proposed_candidate_level_gap(
                reserve_regions.in_regular);
        const int reserve_margin_before =
            reserve_hierarchy.minimum_proposed_candidate_level_gap_margin(
                reserve_profile.target_level_gaps);
        bool reserve_deepened = false;
        if (reserve_margin_before < 0) {
            reserve_deepened =
                reserve_hierarchy.deepen_candidate_over_proposed_coarse(
                    reserve_profile.target_level_gaps).changed();
        }
        const int reserve_margin_after =
            reserve_hierarchy.minimum_proposed_candidate_level_gap_margin(
                reserve_profile.target_level_gaps);
        const int far_regular_gap =
            reserve_hierarchy.minimum_proposed_candidate_level_gap(
                reserve_profile.at_full_target);
        std::size_t reserve_matching_spill = 0;
        for (const int element : reserve_regions.omega_f_elements) {
            if (embedding_children(
                    reserve_hierarchy.proposed_coarse_elements_to_candidate(),
                    element) != 1) {
                ++reserve_matching_spill;
            }
        }
        require(uniform_regular_gap < reserve_target,
                "reserve MRE did not expose the incompatible uniform target");
        require(reserve_profile.full_target_elements > 0
                    && reserve_probe.target_satisfied
                    && reserve_probe.protected_parent_spill == 0
                    && reserve_probe.minimum_gap_margin >= 0
                    && reserve_probe.minimum_full_target_gap >= reserve_target,
                "candidate-only reserve probe rejected the H3/h7 collar");
        require(reserve_margin_before < 0 && reserve_deepened,
                "reserve MRE did not exercise graded candidate deepening");
        require(reserve_hierarchy.candidate_mesh().elems.size()
                    == reserve_probe.final_element_count,
                "candidate-only reserve probe did not predict production NVB");
        require(reserve_margin_after >= 0,
                "graded reserve profile was not closed");
        require(far_regular_gap >= reserve_target,
                "far regular region did not retain the full reserve");
        require(reserve_matching_spill == 0,
                "graded reserve closure spilled into frozen Omega_F");

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

        const std::vector<double> split_indicators{
            1e-12, 4.0, 3.0, 2.0, 1.0};
        const std::vector<char> in_omega_f{true, false, false, false, false};
        const double split_theta = 0.6;
        const SplitRegionalDoerflerMarking split =
            mark_split_regional_doerfler(
                split_indicators, in_omega_f, split_theta);
        require(split.omega_f_mass == split_indicators.front()
                    && split.regular_mass == 10.0,
                "split regional indicator masses are incorrect");
        require(split.marked_omega_f_mass
                        >= split_theta * split.omega_f_mass
                    && split.marked_regular_mass
                        >= split_theta * split.regular_mass,
                "a split region did not meet its own Doerfler target");
        require(split.marked_omega_f_elements.size() == 1
                    && split.marked_omega_f_elements.front() == 0
                    && std::find(
                           split.marked_elements.begin(),
                           split.marked_elements.end(), 0)
                        != split.marked_elements.end(),
                "tiny nonzero Omega_F mass was masked by regular mass");
        for (int element : split.marked_regular_elements)
            require(!in_omega_f[element],
                    "regular split marking entered Omega_F");

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
        const int regular_gap =
            hierarchy.minimum_proposed_reference_level_gap(
                proposed_regions.in_regular);
        require(global_gap == 0
                    && regular_gap >= global_gap,
                "hybrid regular-region gap guard is inconsistent");
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
                  << " physical_l_s=" << physical.l_s
                  << " covered_radius=" << physical.covered_physical_radius
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

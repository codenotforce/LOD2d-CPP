#include "helmholtz/adaptive/hierarchy.h"
#include "helmholtz/model.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace lod2d;
using namespace lod2d::helmholtz;
using namespace lod2d::helmholtz::adaptive;

namespace {

void require(const bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

HelmholtzProblemConfig model_config(const int ell, const double wavenumber = 4.0) {
    HelmholtzProblemConfig config;
    config.H = 2;
    config.h = 5;
    config.ell = ell;
    config.wavenumber = wavenumber;
    config.boundary_beta = 1.0;
    config.mode = HelmholtzPetrovMode::TwoSided;
    config.initial_mesh = make_helmholtz_unit_square_mesh();
    config.patch_solver.kind = HelmholtzPatchSolverKind::DirectSaddle;
    return config;
}

HelmholtzLodModel build_model(
    const ReferenceEpochHierarchy &hierarchy,
    const HelmholtzProblemConfig &config,
    HelmholtzCorrectorPatchCache *cache) {
    return HelmholtzLodModel::build_adaptive(
        config,
        hierarchy.coarse_mesh(),
        hierarchy.coarse_levels(),
        hierarchy.reference_mesh(),
        hierarchy.reference_element_levels(),
        cache);
}

void require_equivalent(
    const HelmholtzLodModel &cached,
    const HelmholtzLodModel &rebuilt) {
    require(cached.correctors().primal.size() == rebuilt.correctors().primal.size(),
            "cached/full corrector counts differ");
    for (std::size_t element = 0;
         element < cached.correctors().primal.size(); ++element) {
        const auto &lhs = cached.correctors().primal[element];
        const auto &rhs = rebuilt.correctors().primal[element];
        require(lhs.size() == rhs.size(),
                "cached/full local corrector sparsity differs");
        for (std::size_t entry = 0; entry < lhs.size(); ++entry) {
            require(lhs[entry].row == rhs[entry].row
                        && lhs[entry].local_coarse_vertex
                            == rhs[entry].local_coarse_vertex
                        && std::abs(lhs[entry].value - rhs[entry].value) <= 1e-13,
                    "cached/full corrector entry differs");
        }
    }
    const double basis_scale = std::max(1.0, rebuilt.corrected_trial_basis().norm());
    require((cached.corrected_trial_basis()
                - rebuilt.corrected_trial_basis()).norm() <= 1e-13 * basis_scale,
            "cached/full corrected bases differ");
    require((cached.coarse_operator() - rebuilt.coarse_operator()).norm()
                <= 1e-12 * std::max(1.0, rebuilt.coarse_operator().norm()),
            "cached/full coarse operators differ");

    ComplexVector load = ComplexVector::Zero(
        static_cast<int>(cached.problem().fine.nodes.size()));
    for (int i = 0; i < load.size(); ++i)
        load[i] = Complex(1.0 + 0.01 * i, -0.25);
    for (const int node : cached.operators().dirichlet_nodes) load[node] = 0.0;
    const HelmholtzLodSolution lhs = cached.solve_load(load);
    const HelmholtzLodSolution rhs = rebuilt.solve_load(load);
    require((lhs.fine_values - rhs.fine_values).norm()
                <= 1e-11 * std::max(1.0, rhs.fine_values.norm()),
            "cached/full LOD solutions differ");
}

void verify_exact_hits_and_full_rebuild_equivalence() {
    ReferenceEpochHierarchy hierarchy(
        make_helmholtz_unit_square_mesh(), 2, 5);
    HelmholtzCorrectorPatchCache cache;
    const HelmholtzProblemConfig config = model_config(1);

    HelmholtzLodModel first = build_model(hierarchy, config, &cache);
    require(first.correctors().diagnostics.patch_cache_hits == 0,
            "an empty corrector cache reported a hit");
    require(first.correctors().diagnostics.patch_cache_misses
                == static_cast<int>(hierarchy.coarse_mesh().elems.size()),
            "an empty corrector cache did not solve every patch");

    HelmholtzLodModel second = build_model(hierarchy, config, &cache);
    require(second.correctors().diagnostics.patch_cache_hits
                == static_cast<int>(hierarchy.coarse_mesh().elems.size())
                && second.correctors().diagnostics.patch_cache_misses == 0,
            "an identical corrector rebuild did not hit every patch");
    HelmholtzLodModel full = build_model(hierarchy, config, nullptr);
    require_equivalent(second, full);

    const auto stats = cache.statistics();
    require(stats.hits >= hierarchy.coarse_mesh().elems.size()
                && stats.entries == hierarchy.coarse_mesh().elems.size(),
            "corrector cache statistics are inconsistent");
}

void verify_local_refinement_matches_full_rebuild() {
    ReferenceEpochHierarchy hierarchy(
        make_helmholtz_unit_square_mesh(), 2, 5);
    HelmholtzCorrectorPatchCache cache;
    const HelmholtzProblemConfig config = model_config(1);
    (void)build_model(hierarchy, config, &cache);

    const ReferenceEpochRefinementResult refined =
        hierarchy.refine_coarse_preserving_reference({0});
    require(refined.changed(), "local H refinement did not change the hierarchy");
    HelmholtzLodModel cached = build_model(hierarchy, config, &cache);
    HelmholtzLodModel full = build_model(hierarchy, config, nullptr);
    require(cached.correctors().diagnostics.patch_cache_misses > 0,
            "local H refinement did not rebuild any changed patch system");
    require(cached.correctors().diagnostics.patch_cache_hits
                + cached.correctors().diagnostics.patch_cache_misses
            == static_cast<int>(hierarchy.coarse_mesh().elems.size()),
            "local H cache accounting does not match the patch count");
    require_equivalent(cached, full);
}

void verify_key_changes_fail_closed() {
    ReferenceEpochHierarchy hierarchy(
        make_helmholtz_unit_square_mesh(), 2, 5);
    HelmholtzCorrectorPatchCache cache;
    (void)build_model(hierarchy, model_config(1), &cache);

    HelmholtzLodModel changed_ell = build_model(hierarchy, model_config(2), &cache);
    require(changed_ell.correctors().diagnostics.patch_cache_hits
                + changed_ell.correctors().diagnostics.patch_cache_misses
            == static_cast<int>(hierarchy.coarse_mesh().elems.size()),
            "ell change produced inconsistent cache accounting");
    require_equivalent(
        changed_ell, build_model(hierarchy, model_config(2), nullptr));

    HelmholtzLodModel changed_kappa =
        build_model(hierarchy, model_config(2, 5.0), &cache);
    require(changed_kappa.correctors().diagnostics.patch_cache_hits == 0
                && changed_kappa.correctors().diagnostics.patch_cache_misses
                    == static_cast<int>(hierarchy.coarse_mesh().elems.size()),
            "PDE change reused stale corrector patches");
    require_equivalent(
        changed_kappa, build_model(hierarchy, model_config(2, 5.0), nullptr));
}

void verify_oversized_guard_is_observable() {
    ReferenceEpochHierarchy hierarchy(
        make_helmholtz_unit_square_mesh(), 2, 5);
    HelmholtzCorrectorPatchCache cache(8, 1);
    const HelmholtzLodModel model = build_model(
        hierarchy, model_config(1), &cache);
    const int patches = static_cast<int>(hierarchy.coarse_mesh().elems.size());
    require(model.correctors().diagnostics.patch_cache_oversized_misses
                == patches
            && model.correctors().diagnostics.patch_cache_misses == patches,
            "oversized corrector cache misses were not reported separately");
    const auto stats = cache.statistics();
    require(stats.oversized_misses == static_cast<std::size_t>(patches)
                && stats.entries == 0,
            "oversized corrector cache guard stored an ineligible patch");
}

void verify_memory_budget_is_observable() {
    ReferenceEpochHierarchy hierarchy(
        make_helmholtz_unit_square_mesh(), 2, 5);
    HelmholtzCorrectorPatchCache cache(8, 4096, 1);
    const HelmholtzLodModel model = build_model(
        hierarchy, model_config(1), &cache);
    const int patches = static_cast<int>(hierarchy.coarse_mesh().elems.size());
    require(model.correctors().diagnostics.patch_cache_budget_rejections
                == patches
            && model.correctors().diagnostics.patch_cache_oversized_misses == 0,
            "corrector cache byte-budget rejections were not distinguished");
    const auto stats = cache.statistics();
    require(stats.budget_rejections == static_cast<std::size_t>(patches)
                && stats.current_bytes == 0 && stats.entries == 0,
            "corrector cache exceeded a one-byte memory budget");
}

} // namespace

int main() {
    try {
        verify_exact_hits_and_full_rebuild_equivalence();
        verify_local_refinement_matches_full_rebuild();
        verify_key_changes_fail_closed();
        verify_oversized_guard_is_observable();
        verify_memory_budget_is_observable();
        std::cout << "Helmholtz corrector patch cache tests passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_corrector_cache failed: "
                  << error.what() << '\n';
        return 1;
    }
}

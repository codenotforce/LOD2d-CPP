#include "helmholtz/adaptive/candidate_flux.h"
#include "helmholtz/benchmarks/paper_cases.h"
#include "helmholtz/model.h"
#include "mesh/refine.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
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

void verify_case(PaperCase id) {
    const PaperCaseData data = id == PaperCase::S
        ? make_paper_case(id, 16.0, 0.0, 0.5, false, 0.05)
        : make_paper_case(id, 16.0);
    const TriMesh mesh = refine_mesh_nvb(data.initial_mesh, 2).mesh;
    const HelmholtzOperators operators = assemble_helmholtz_operators(
        mesh, data.wavenumber);
    ComplexVector values(mesh.nodes.size());
    for (int node = 0; node < static_cast<int>(mesh.nodes.size()); ++node)
        values(node) = data.exact(mesh.nodes[node]);

    CandidateFluxConfig config;
    config.doerfler_theta = 0.5;
    config.quadrature_context = data.quadrature_context;
    config.compute_discrete_residual_audit = true;
    const CandidateFluxRT2Result flux = reconstruct_candidate_flux_rt2(
        mesh, operators, data.source, values, config);
    require(flux.global_reconstruction,
            "global candidate reconstruction was labeled active-region only");
    require(flux.maximum_patch_compatibility_error < 2e-9,
            "candidate patch compatibility was not corrected");
    require(flux.maximum_element_divergence_residual < 2e-9,
            "candidate RT2 divergence identity failed");
    require(flux.maximum_normal_continuity_residual < 2e-9,
            "candidate RT2 normal continuity identity failed");
    require(flux.maximum_boundary_flux_residual < 2e-9,
            "candidate RT2 boundary-flux identity failed");
    require(flux.discrete_residual_audit_performed
                && flux.eta_eq + 2e-9
                    >= flux.discrete_residual_dual_norm,
            "candidate indicator does not control the discrete residual dual norm");
    require(!flux.marked_elements.empty(),
            "candidate Doerfler marking is empty");
    require(flux.parallel_threads >= 1
                && flux.time_prepare >= 0.0
                && flux.time_patch_solve >= 0.0
                && flux.time_deterministic_merge >= 0.0
                && flux.time_estimator_and_audit >= 0.0,
            "candidate RT2 phase diagnostics are invalid");
    double marked = 0.0;
    for (int element : flux.marked_elements)
        marked += flux.element_eta_squared[element];
    const double total = std::accumulate(
        flux.element_eta_squared.begin(),
        flux.element_eta_squared.end(), 0.0);
    require(marked + 2e-12 >= config.doerfler_theta * total,
            "candidate Doerfler marking misses its target");

    CandidateFluxConfig closure_aware = config;
    closure_aware.compute_discrete_residual_audit = false;
    closure_aware.closure_cost_aware_marking = true;
    closure_aware.closure_cost_candidate_pool_factor = 2;
    const CandidateFluxRT2Result economical = reconstruct_candidate_flux_rt2(
        mesh, operators, data.source, values, closure_aware);
    double economical_marked = 0.0;
    for (const int element : economical.marked_elements)
        economical_marked += economical.element_eta_squared[element];
    require(economical_marked + 2e-12
                >= closure_aware.doerfler_theta * total,
            "closure-cost-aware candidate marking misses its raw Doerfler target");
    require(economical.marking_candidate_pool
                >= economical.marked_elements.size()
                && economical.estimated_selected_closure_cost > 0
                && economical.time_marking >= 0.0,
            "closure-cost-aware candidate diagnostics are invalid");
    std::cout << "case=" << static_cast<int>(id)
              << " eta_eq=" << flux.eta_eq
              << " residual_dual=" << flux.discrete_residual_dual_norm
              << " compatibility_error="
              << flux.maximum_patch_compatibility_error
              << " divergence_error="
              << flux.maximum_element_divergence_residual
              << " boundary_flux_error="
              << flux.maximum_boundary_flux_residual
              << " normal_jump="
              << flux.maximum_normal_continuity_residual
              << " patches=" << flux.patches.size()
              << " threads=" << flux.parallel_threads
              << " marked=" << flux.marked_elements.size() << '\n';

    CandidateFluxConfig active = config;
    active.compute_discrete_residual_audit = false;
    active.active_elements = {0};
    const CandidateFluxRT2Result restricted = reconstruct_candidate_flux_rt2(
        mesh, operators, data.source, values, active);
    require(!restricted.global_reconstruction
                && restricted.patches.size() < flux.patches.size(),
            "active-region candidate reconstruction was not restricted");
}

} // namespace

int main() {
    try {
        verify_case(PaperCase::R1);
        verify_case(PaperCase::S);
        std::cout << "Compatibility-corrected candidate RT2/P2 flux passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_candidate_flux failed: "
                  << error.what() << '\n';
        return 1;
    }
}

#include "helmholtz/adaptive/hierarchy.h"
#include "helmholtz/boundary.h"
#include "helmholtz/model.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

using namespace lod2d;
using namespace lod2d::helmholtz;
using namespace lod2d::helmholtz::adaptive;

namespace {

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

using TriangleSignature = std::array<long long, 6>;

std::vector<TriangleSignature> canonical_triangles(const TriMesh &mesh) {
    std::vector<TriangleSignature> result;
    result.reserve(mesh.elems.size());
    for (const Triangle &triangle : mesh.elems) {
        std::array<std::pair<long long, long long>, 3> vertices;
        for (int local = 0; local < 3; ++local) {
            const Point2 &point = mesh.nodes[triangle[local]];
            vertices[local] = {
                std::llround(point.x() * 1e12),
                std::llround(point.y() * 1e12)};
        }
        std::sort(vertices.begin(), vertices.end());
        TriangleSignature signature{};
        for (int local = 0; local < 3; ++local) {
            signature[2 * local] = vertices[local].first;
            signature[2 * local + 1] = vertices[local].second;
        }
        result.push_back(signature);
    }
    std::sort(result.begin(), result.end());
    return result;
}

Eigen::MatrixXd node_coordinates(const TriMesh &mesh) {
    Eigen::MatrixXd result(mesh.nodes.size(), 2);
    for (int node = 0; node < static_cast<int>(mesh.nodes.size()); ++node)
        result.row(node) = mesh.nodes[node].transpose();
    return result;
}

void verify_embeddings(const ReferenceEpochHierarchy &hierarchy) {
    require(
        (hierarchy.reference_to_ambient()
             * hierarchy.coarse_to_reference()
         - hierarchy.coarse_to_ambient()).norm() < 1e-11,
        "reference-epoch nodal prolongations do not compose");
    require(
        (hierarchy.reference_elements_to_ambient()
             * hierarchy.coarse_elements_to_reference()
         - hierarchy.coarse_elements_to_ambient()).norm() < 1e-11,
        "reference-epoch element prolongations do not compose");
    require(
        (hierarchy.reference_dg_to_ambient()
             * hierarchy.coarse_dg_to_reference()
         - hierarchy.coarse_dg_to_ambient()).norm() < 1e-11,
        "reference-epoch DG prolongations do not compose");

    require(
        (hierarchy.coarse_to_reference()
             * node_coordinates(hierarchy.coarse_mesh())
         - node_coordinates(hierarchy.reference_mesh())).norm() < 1e-10,
        "coarse-to-reference prolongation does not reproduce coordinates");
    require(
        (hierarchy.reference_to_ambient()
             * node_coordinates(hierarchy.reference_mesh())
         - node_coordinates(hierarchy.ambient_mesh())).norm() < 1e-10,
        "reference-to-ambient prolongation does not reproduce coordinates");

    Eigen::MatrixXd expected = Eigen::MatrixXd::Identity(
        hierarchy.coarse_mesh().nodes.size(),
        hierarchy.coarse_mesh().nodes.size());
    for (int node : dirichlet_nodes(hierarchy.coarse_mesh()))
        expected.row(node).setZero();
    const Eigen::MatrixXd reference_right_inverse =
        hierarchy.reference_quasi_interpolation()
        * hierarchy.coarse_to_reference();
    const Eigen::MatrixXd ambient_right_inverse =
        hierarchy.ambient_quasi_interpolation()
        * hierarchy.coarse_to_ambient();
    require((reference_right_inverse - expected).norm() < 1e-9,
            "reference I_H is not a right inverse");
    require((ambient_right_inverse - expected).norm() < 1e-9,
            "ambient I_H is not a right inverse");
    require(hierarchy.reference_parent_coarse_elements().size()
                == hierarchy.reference_mesh().elems.size(),
            "reference parent map has the wrong size");
    require(hierarchy.ambient_parent_coarse_elements().size()
                == hierarchy.ambient_mesh().elems.size(),
            "ambient/coarse parent map has the wrong size");
    require(hierarchy.ambient_parent_reference_elements().size()
                == hierarchy.ambient_mesh().elems.size(),
            "ambient/reference parent map has the wrong size");
}

int deepest_coarse_element(const ReferenceEpochHierarchy &hierarchy) {
    return static_cast<int>(std::distance(
        hierarchy.coarse_levels().begin(),
        std::max_element(
            hierarchy.coarse_levels().begin(),
            hierarchy.coarse_levels().end())));
}

void verify_incremental_embedding_matches_full_rebuild() {
    const TriMesh initial = make_helmholtz_unit_square_mesh();
    const RefineOutput old_parent = refine_mesh_nvb(initial, 2);
    const RefineOutput child = refine_mesh_nvb(initial, 5);
    const RefineOutput old_embedding = build_nested_mesh_embedding(
        old_parent.mesh, child.mesh);
    const std::vector<int> old_child_parents = fine_element_parents(
        old_embedding.P_elem,
        static_cast<int>(child.mesh.elems.size()),
        static_cast<int>(old_parent.mesh.elems.size()));
    const RefineOutput parent_step = bisect_newest_vertex(
        old_parent.mesh, {0, 3, 7});
    const RefineOutput incremental =
        update_nested_mesh_embedding_after_parent_refinement(
            old_parent.mesh, parent_step, child.mesh, old_child_parents);
    const RefineOutput rebuilt = build_nested_mesh_embedding(
        parent_step.mesh, child.mesh);
    require((incremental.P_node - rebuilt.P_node).norm() < 1e-12,
            "incremental nodal embedding differs from a full rebuild");
    require((incremental.P_elem - rebuilt.P_elem).norm() < 1e-12,
            "incremental element embedding differs from a full rebuild");
    require((incremental.P_dg - rebuilt.P_dg).norm() < 1e-12,
            "incremental DG embedding differs from a full rebuild");
}

void verify_fixed_reference_and_ambient_ratio() {
    constexpr double rho_star = 0.45;
    ReferenceEpochHierarchy hierarchy(
        make_helmholtz_unit_square_mesh(), 2, 7);
    const std::vector<TriangleSignature> fixed_reference =
        canonical_triangles(hierarchy.reference_mesh());
    const std::size_t fixed_reference_nodes =
        hierarchy.reference_mesh().nodes.size();
    const std::vector<int> fixed_reference_levels =
        hierarchy.reference_element_levels();
    const std::uint64_t fixed_reference_version =
        hierarchy.reference_mesh_version();
    const std::size_t initial_ambient_elements =
        hierarchy.ambient_mesh().elems.size();

    hierarchy.enforce_ambient_ratio(rho_star);
    for (int step = 0; step < 3; ++step) {
        const ReferenceEpochRefinementResult refinement =
            hierarchy.refine_coarse_preserving_reference(
                {deepest_coarse_element(hierarchy)});
        require(refinement.status == ReferenceEpochRefinementStatus::Refined,
                "a nested local coarse refinement was rejected");
        const AmbientRatioEnforcementResult ambient =
            hierarchy.enforce_ambient_ratio(rho_star);
        require(ambient.maximum_ratio <= rho_star + 1e-11,
                "ambient shadow violates the requested diameter ratio");
        require(hierarchy.reference_embedding_holds(),
                "accepted coarse refinement broke the reference embedding");
        require(canonical_triangles(hierarchy.reference_mesh())
                    == fixed_reference,
                "local coarse refinement changed the reference mesh");
        require(hierarchy.reference_mesh().nodes.size()
                    == fixed_reference_nodes,
                "local coarse refinement changed reference nodes");
        require(hierarchy.reference_element_levels()
                    == fixed_reference_levels,
                "local coarse refinement changed reference levels");
        require(hierarchy.reference_mesh_version()
                    == fixed_reference_version,
                "local coarse refinement advanced the reference version");
        verify_embeddings(hierarchy);
    }

    require(hierarchy.coarse_mesh_version() == 3,
            "three local H refinements did not advance the coarse version");
    require(hierarchy.reference_epoch() == 0,
            "local H refinement unexpectedly started a new reference epoch");
    require(hierarchy.ambient_mesh().elems.size() > initial_ambient_elements,
            "ambient mesh did not follow local H refinement");
    require(*std::min_element(
                hierarchy.ambient_element_levels().begin(),
                hierarchy.ambient_element_levels().end()) == 7,
            "ambient ratio enforcement refined the entire shadow mesh");
    require(*std::max_element(
                hierarchy.ambient_element_levels().begin(),
                hierarchy.ambient_element_levels().end()) > 7,
            "ambient ratio enforcement did not refine a violating region");
    require(hierarchy.ambient_ratio() <= rho_star + 1e-11,
            "reported ambient ratio exceeds the target");
    require(std::abs(boundary_measure(
                hierarchy.reference_mesh(), BoundaryTag::Robin) - 4.0)
                < 2e-12,
            "reference boundary tags changed during the epoch");
    require(std::abs(boundary_measure(
                hierarchy.ambient_mesh(), BoundaryTag::Robin) - 4.0)
                < 2e-12,
            "ambient refinement changed the physical boundary");
}

void verify_transactional_refresh_boundary() {
    constexpr double rho_star = 0.45;
    ReferenceEpochHierarchy hierarchy(
        make_helmholtz_unit_square_mesh(), 0, 1);
    std::vector<int> marked(hierarchy.coarse_mesh().elems.size());
    std::iota(marked.begin(), marked.end(), 0);
    require(hierarchy.refine_coarse_preserving_reference(marked).changed(),
            "coarse refinement up to the reference mesh was rejected");
    hierarchy.enforce_ambient_ratio(rho_star);

    const std::vector<TriangleSignature> coarse_before =
        canonical_triangles(hierarchy.coarse_mesh());
    const std::vector<TriangleSignature> reference_before =
        canonical_triangles(hierarchy.reference_mesh());
    const std::vector<TriangleSignature> ambient_before =
        canonical_triangles(hierarchy.ambient_mesh());
    const std::vector<int> coarse_levels_before = hierarchy.coarse_levels();
    const std::uint64_t coarse_version_before = hierarchy.coarse_mesh_version();
    const std::uint64_t reference_version_before =
        hierarchy.reference_mesh_version();
    const std::uint64_t ambient_version_before =
        hierarchy.ambient_mesh_version();

    marked.resize(hierarchy.coarse_mesh().elems.size());
    std::iota(marked.begin(), marked.end(), 0);
    const ReferenceEpochRefinementResult blocked =
        hierarchy.refine_coarse_preserving_reference(marked);
    require(blocked.status
                == ReferenceEpochRefinementStatus::ReferenceRefreshRequired,
            "coarse refinement beyond the reference did not request refresh");
    require(!blocked.detail.empty(),
            "reference-refresh request has no structured explanation");
    require(canonical_triangles(hierarchy.coarse_mesh()) == coarse_before,
            "rejected coarse refinement partially changed the coarse mesh");
    require(canonical_triangles(hierarchy.reference_mesh()) == reference_before,
            "rejected coarse refinement changed the reference mesh");
    require(canonical_triangles(hierarchy.ambient_mesh()) == ambient_before,
            "rejected coarse refinement changed the ambient mesh");
    require(hierarchy.coarse_levels() == coarse_levels_before,
            "rejected coarse refinement changed coarse levels");
    require(hierarchy.coarse_mesh_version() == coarse_version_before
                && hierarchy.reference_mesh_version()
                    == reference_version_before
                && hierarchy.ambient_mesh_version()
                    == ambient_version_before,
            "rejected coarse refinement advanced a mesh version");

    hierarchy.refresh_reference_from_ambient();
    require(hierarchy.reference_epoch() == 1,
            "explicit reference refresh did not advance the epoch");
    require(hierarchy.reference_mesh_version() == reference_version_before + 1,
            "explicit reference refresh did not advance its version");
    require(canonical_triangles(hierarchy.reference_mesh()) == ambient_before,
            "reference refresh did not promote the ambient mesh");
    require(canonical_triangles(hierarchy.ambient_mesh()) == ambient_before,
            "reference refresh did not reinitialize ambient from reference");
    require(hierarchy.refine_coarse_preserving_reference(marked).changed(),
            "coarse refinement remained blocked after explicit refresh");
    verify_embeddings(hierarchy);
}

void verify_mixed_boundary_contract() {
    ReferenceEpochHierarchy hierarchy(make_helmholtz_l_shape_mesh(), 0, 3);
    require(hierarchy.refine_coarse_preserving_reference({0}).changed(),
            "mixed-boundary local H refinement was rejected");
    hierarchy.enforce_ambient_ratio(0.5);
    require(std::abs(boundary_measure(
                hierarchy.reference_mesh(), BoundaryTag::Dirichlet) - 2.0)
                < 2e-12,
            "reference epoch changed the L-shape Dirichlet boundary");
    require(std::abs(boundary_measure(
                hierarchy.ambient_mesh(), BoundaryTag::Robin) - 6.0)
                < 2e-12,
            "ambient refinement changed the L-shape Robin boundary");
    verify_embeddings(hierarchy);
}

void verify_incremental_candidate_refinement_matches_full_rebuild() {
    ReferenceEpochHierarchy hierarchy(
        make_helmholtz_unit_square_mesh(), 2, 7);
    const Eigen::SparseMatrix<double> fixed_reference_interpolation =
        hierarchy.reference_quasi_interpolation();
    for (int step = 0; step < 2; ++step) {
        const int count = static_cast<int>(hierarchy.candidate_mesh().elems.size());
        const std::vector<int> marked{0, count / 3, 2 * count / 3};
        const ReferenceEpochRefinementResult refined =
            hierarchy.enrich_candidate(marked);
        require(refined.changed(),
                "incremental candidate refinement did not change the mesh");
        require(refined.time_nvb_refine >= 0.0
                    && refined.time_embedding_composition >= 0.0
                    && refined.time_parent_map_update >= 0.0
                    && refined.time_candidate_quasi_interpolation >= 0.0
                    && refined.time_embedding_validation >= 0.0,
                "incremental candidate phase timings are invalid");

        const RefineOutput reference_rebuild = build_nested_mesh_embedding(
            hierarchy.reference_mesh(), hierarchy.candidate_mesh());
        const RefineOutput coarse_rebuild = build_nested_mesh_embedding(
            hierarchy.coarse_mesh(), hierarchy.candidate_mesh());
        require((hierarchy.reference_to_candidate()
                    - reference_rebuild.P_node).norm() < 1e-12,
                "incremental reference/candidate nodal map differs from rebuild");
        require((hierarchy.reference_elements_to_candidate()
                    - reference_rebuild.P_elem).norm() < 1e-12,
                "incremental reference/candidate element map differs from rebuild");
        require((hierarchy.reference_dg_to_candidate()
                    - reference_rebuild.P_dg).norm() < 1e-12,
                "incremental reference/candidate DG map differs from rebuild");
        require((hierarchy.coarse_to_candidate()
                    - coarse_rebuild.P_node).norm() < 1e-12,
                "incremental coarse/candidate nodal map differs from rebuild");
        require((hierarchy.coarse_elements_to_candidate()
                    - coarse_rebuild.P_elem).norm() < 1e-12,
                "incremental coarse/candidate element map differs from rebuild");
        require((hierarchy.coarse_dg_to_candidate()
                    - coarse_rebuild.P_dg).norm() < 1e-12,
                "incremental coarse/candidate DG map differs from rebuild");
        require(hierarchy.candidate_parent_reference_elements()
                    == fine_element_parents(
                        reference_rebuild.P_elem,
                        static_cast<int>(hierarchy.candidate_mesh().elems.size()),
                        static_cast<int>(hierarchy.reference_mesh().elems.size())),
                "incremental candidate/reference parent map differs from rebuild");
        require(hierarchy.candidate_parent_coarse_elements()
                    == fine_element_parents(
                        coarse_rebuild.P_elem,
                        static_cast<int>(hierarchy.candidate_mesh().elems.size()),
                        static_cast<int>(hierarchy.coarse_mesh().elems.size())),
                "incremental candidate/coarse parent map differs from rebuild");
        require((hierarchy.reference_quasi_interpolation()
                    - fixed_reference_interpolation).norm() == 0.0,
                "candidate-only refinement rebuilt the fixed reference interpolation");
        verify_embeddings(hierarchy);
    }
}

void verify_revised_candidate_transaction() {
    ReferenceEpochHierarchy hierarchy(
        make_helmholtz_unit_square_mesh(), 0, 1);
    hierarchy.begin_reference_epoch();
    const std::vector<TriangleSignature> reference_epoch_zero =
        canonical_triangles(hierarchy.reference_mesh());
    require(canonical_triangles(hierarchy.candidate_mesh())
                == reference_epoch_zero,
            "candidate was not initialized from the reference mesh");

    std::vector<int> marked(hierarchy.coarse_mesh().elems.size());
    std::iota(marked.begin(), marked.end(), 0);
    require(hierarchy.propose_coarse_refinement(marked).changed(),
            "first coarse proposal was not created");
    require(hierarchy.has_proposed_coarse_refinement(),
            "pending coarse proposal was not retained");
    require(hierarchy.candidate_contains_proposed_coarse(),
            "reference-inherited candidate does not contain a representable proposal");
    require(hierarchy.commit_coarse_refinement().changed(),
            "representable coarse proposal was not committed");

    marked.resize(hierarchy.coarse_mesh().elems.size());
    std::iota(marked.begin(), marked.end(), 0);
    const std::size_t committed_elements = hierarchy.coarse_mesh().elems.size();
    require(hierarchy.propose_coarse_refinement(marked).changed(),
            "second coarse proposal was not created");
    require(hierarchy.coarse_mesh().elems.size() == committed_elements,
            "coarse proposal mutated the committed coarse mesh");
    require(!hierarchy.candidate_contains_proposed_coarse(),
            "under-resolved candidate unexpectedly contains the proposal");
    require(hierarchy.close_candidate_over_proposed_coarse().changed(),
            "candidate closure did not refine the proposed region");
    require(hierarchy.candidate_contains_proposed_coarse(),
            "candidate closure did not contain the proposed coarse mesh");
    require(hierarchy.minimum_proposed_candidate_level_gap() >= 0,
            "cached proposed/candidate embedding has an invalid level gap");
    require(canonical_triangles(hierarchy.reference_mesh())
                == reference_epoch_zero,
            "candidate enrichment changed the fixed reference mesh");

    const ReferenceEpochRefinementResult blocked =
        hierarchy.commit_coarse_refinement();
    require(blocked.status
                == ReferenceEpochRefinementStatus::ReferenceRefreshRequired,
            "coarse proposal crossed the reference without requesting refresh");
    require(hierarchy.coarse_mesh().elems.size() == committed_elements,
            "failed commit partially mutated the coarse mesh");

    require(hierarchy.deepen_candidate_over_proposed_coarse(2).changed(),
            "pre-refresh candidate deepening did not refine the candidate");
    require(hierarchy.minimum_proposed_candidate_level_gap() >= 2,
            "pre-refresh candidate deepening missed its target gap");
    const std::vector<TriangleSignature> enriched_candidate =
        canonical_triangles(hierarchy.candidate_mesh());
    hierarchy.refresh_reference_from_candidate();
    require(hierarchy.minimum_proposed_reference_level_gap() >= 2,
            "promoted reference did not preserve the target gap");
    require(hierarchy.reference_epoch() == 1,
            "candidate promotion did not advance the reference epoch");
    require(canonical_triangles(hierarchy.reference_mesh())
                == enriched_candidate,
            "reference refresh did not promote the candidate mesh");
    require(hierarchy.commit_coarse_refinement().changed(),
            "coarse proposal was not committed after reference refresh");
    require(hierarchy.reference_embedding_holds(),
            "committed coarse mesh is not contained in the refreshed reference");
    verify_embeddings(hierarchy);
}

} // namespace

int main() {
    try {
        verify_incremental_embedding_matches_full_rebuild();
        verify_incremental_candidate_refinement_matches_full_rebuild();
        verify_fixed_reference_and_ambient_ratio();
        verify_transactional_refresh_boundary();
        verify_mixed_boundary_contract();
        verify_revised_candidate_transaction();
        std::cout << "Reference-epoch hierarchy invariants passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_reference_epoch_hierarchy failed: "
                  << error.what() << '\n';
        return 1;
    }
}

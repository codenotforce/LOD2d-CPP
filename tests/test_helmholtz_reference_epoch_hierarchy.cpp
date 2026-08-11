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

} // namespace

int main() {
    try {
        verify_fixed_reference_and_ambient_ratio();
        verify_transactional_refresh_boundary();
        verify_mixed_boundary_contract();
        std::cout << "Reference-epoch hierarchy invariants passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_reference_epoch_hierarchy failed: "
                  << error.what() << '\n';
        return 1;
    }
}

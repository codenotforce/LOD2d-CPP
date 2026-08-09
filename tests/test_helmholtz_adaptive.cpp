#include "helmholtz/adaptive/estimator.h"
#include "helmholtz/adaptive/hierarchy.h"
#include "helmholtz/boundary.h"
#include "helmholtz/model.h"

#include <algorithm>
#include <cmath>
#include <array>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace lod2d;
using namespace lod2d::helmholtz;
using namespace lod2d::helmholtz::adaptive;

namespace {

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

ComplexFunction source_function() {
    return [](const Point2 &point) {
        const double dx = point.x() - 0.31;
        const double dy = point.y() - 0.63;
        return Complex(std::exp(-45.0 * (dx * dx + dy * dy)), 0.15 * point.x());
    };
}

using TriangleSignature = std::array<long long, 6>;

std::vector<TriangleSignature> canonical_triangles(const TriMesh &mesh) {
    std::vector<TriangleSignature> result;
    result.reserve(mesh.elems.size());
    for (const Triangle &triangle : mesh.elems) {
        std::array<std::pair<long long, long long>, 3> vertices;
        for (int i = 0; i < 3; ++i) {
            const Point2 &point = mesh.nodes[triangle[i]];
            vertices[i] = {
                std::llround(point.x() * 1e12),
                std::llround(point.y() * 1e12)};
        }
        std::sort(vertices.begin(), vertices.end());
        TriangleSignature signature{};
        for (int i = 0; i < 3; ++i) {
            signature[2 * i] = vertices[i].first;
            signature[2 * i + 1] = vertices[i].second;
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

void verify_three_level_prolongations(const AdaptiveMeshHierarchy &hierarchy) {
    const Eigen::SparseMatrix<double> composed =
        hierarchy.fine_to_cert_audit() * hierarchy.coarse_to_fine();
    require((composed - hierarchy.coarse_to_cert_audit()).norm() < 1e-12,
            "P_H_audit is inconsistent with P_h_audit P_H_h");
    const Eigen::SparseMatrix<double> composed_elements =
        hierarchy.fine_elements_to_cert_audit()
        * hierarchy.coarse_elements_to_fine();
    require((composed_elements - hierarchy.coarse_elements_to_cert_audit()).norm() < 1e-12,
            "element P_H_audit is inconsistent with P_h_audit P_H_h");
    const Eigen::SparseMatrix<double> composed_dg =
        hierarchy.fine_dg_to_cert_audit() * hierarchy.coarse_dg_to_fine();
    require((composed_dg - hierarchy.coarse_dg_to_cert_audit()).norm() < 1e-12,
            "DG P_H_audit is inconsistent with P_h_audit P_H_h");
    require(hierarchy.fine_parent_coarse_elements().size()
                == hierarchy.fine_mesh().elems.size(),
            "fine parent map has the wrong size");
    require(hierarchy.cert_audit_parent_fine_elements().size()
                == hierarchy.cert_audit_mesh().elems.size(),
            "cert-audit parent map has the wrong size");
    require((hierarchy.coarse_to_fine() * node_coordinates(hierarchy.coarse_mesh())
             - node_coordinates(hierarchy.fine_mesh())).norm() < 1e-11,
            "P_H_h does not reproduce fine node coordinates");
    require((hierarchy.fine_to_cert_audit() * node_coordinates(hierarchy.fine_mesh())
             - node_coordinates(hierarchy.cert_audit_mesh())).norm() < 1e-11,
            "P_h_audit does not reproduce audit node coordinates");
    require((hierarchy.coarse_to_cert_audit()
                 * node_coordinates(hierarchy.coarse_mesh())
             - node_coordinates(hierarchy.cert_audit_mesh())).norm() < 1e-11,
            "P_H_audit does not reproduce audit node coordinates");

    Eigen::MatrixXd expected = Eigen::MatrixXd::Identity(
        hierarchy.coarse_mesh().nodes.size(), hierarchy.coarse_mesh().nodes.size());
    for (int node : dirichlet_nodes(hierarchy.coarse_mesh()))
        expected.row(node).setZero();
    const Eigen::MatrixXd fine_right_inverse =
        hierarchy.fine_quasi_interpolation() * hierarchy.coarse_to_fine();
    const Eigen::MatrixXd audit_right_inverse =
        hierarchy.cert_audit_quasi_interpolation()
        * hierarchy.coarse_to_cert_audit();
    require((fine_right_inverse - expected).norm() < 1e-9,
            "fine I_H is not a right inverse of P_H_h");
    require((audit_right_inverse - expected).norm() < 1e-9,
            "cert-audit I_H is not a right inverse of P_H_audit");
}

void verify_hierarchy_and_estimator() {
    AdaptiveMeshHierarchy hierarchy(make_helmholtz_unit_square_mesh(), 1, 4);
    const NestedFineMesh initial_fine = hierarchy.build_nested_fine_mesh();
    require(std::all_of(
                initial_fine.element_levels.begin(),
                initial_fine.element_levels.end(),
                [](int level) { return level == 4; }),
            "fixed fine completion has the wrong level");
    const TriMesh uniform_fine =
        refine_mesh_nvb(make_helmholtz_unit_square_mesh(), 4).mesh;
    require(canonical_triangles(initial_fine.refinement.mesh) == canonical_triangles(uniform_fine),
            "initial fixed fine completion differs from the uniform master mesh");
    require(initial_fine.refinement.mesh.nodes.size() == uniform_fine.nodes.size(),
            "initial fixed fine completion contains duplicate node coordinates");
    verify_three_level_prolongations(hierarchy);

    const std::size_t initial_audit_elements = hierarchy.cert_audit_mesh().elems.size();
    hierarchy.refine_cert_audit_from_fine_elements({0});
    validate_boundary_tags(hierarchy.cert_audit_mesh());
    require(hierarchy.cert_audit_mesh().elems.size() > initial_audit_elements,
            "cert-audit refinement did not refine the selected fine element");
    require(*std::min_element(
                hierarchy.cert_audit_element_levels().begin(),
                hierarchy.cert_audit_element_levels().end()) == hierarchy.fine_level(),
            "local cert-audit refinement unexpectedly refined the entire fine mesh");
    require(*std::max_element(
                hierarchy.cert_audit_element_levels().begin(),
                hierarchy.cert_audit_element_levels().end()) > hierarchy.fine_level(),
            "cert-audit mesh is not a strict refinement of the fine mesh");
    require(std::abs(boundary_measure(
                hierarchy.cert_audit_mesh(), BoundaryTag::Robin) - 4.0) < 2e-12,
            "cert-audit refinement changed the physical boundary measure");
    verify_three_level_prolongations(hierarchy);
    const std::vector<TriangleSignature> audit_before_coarse_refinement =
        canonical_triangles(hierarchy.cert_audit_mesh());

    hierarchy.refine({0});
    require(*std::min_element(hierarchy.coarse_levels().begin(), hierarchy.coarse_levels().end()) == 1,
            "local NVB unexpectedly refined every coarse element");
    require(*std::max_element(hierarchy.coarse_levels().begin(), hierarchy.coarse_levels().end()) > 1,
            "local NVB did not refine the marked region");
    require(hierarchy.coarse_element_ids().size() == hierarchy.coarse_mesh().elems.size(),
            "stable coarse element IDs do not match the mesh");
    const NestedFineMesh adaptive_fine = hierarchy.build_nested_fine_mesh();
    require(canonical_triangles(adaptive_fine.refinement.mesh) == canonical_triangles(uniform_fine),
            "adaptive completion changed the fixed fine mesh geometry");

    HelmholtzProblemConfig config;
    config.H = 1;
    require(adaptive_fine.refinement.mesh.nodes.size() == uniform_fine.nodes.size(),
            "adaptive completion contains duplicate node coordinates");
    require(canonical_triangles(hierarchy.cert_audit_mesh())
                == audit_before_coarse_refinement,
            "coarse refinement changed the independent cert-audit mesh");
    verify_three_level_prolongations(hierarchy);
    for (int step = 0; step < 3; ++step) {
        std::vector<int> marked;
        for (int t = 0; t < static_cast<int>(hierarchy.coarse_levels().size()); ++t) {
            if (hierarchy.coarse_levels()[t] + 1 < hierarchy.fine_level()
                && (t + step) % 3 == 0) {
                marked.push_back(t);
            }
        }
        if (marked.empty()) break;
        hierarchy.refine(marked);
        const NestedFineMesh completed = hierarchy.build_nested_fine_mesh();
        require(
            canonical_triangles(completed.refinement.mesh) == canonical_triangles(uniform_fine),
            "repeated adaptive completion changed the fixed fine mesh geometry");
        require(completed.refinement.mesh.nodes.size() == uniform_fine.nodes.size(),
                "repeated adaptive completion contains duplicate node coordinates");
    }

    std::vector<int> minimum_audit_level(
        hierarchy.fine_mesh().elems.size(), std::numeric_limits<int>::max());
    std::vector<int> maximum_audit_level(
        hierarchy.fine_mesh().elems.size(), std::numeric_limits<int>::min());
    for (int audit = 0;
         audit < static_cast<int>(hierarchy.cert_audit_parent_fine_elements().size());
         ++audit) {
        const int fine_parent = hierarchy.cert_audit_parent_fine_elements()[audit];
        minimum_audit_level[fine_parent] = std::min(
            minimum_audit_level[fine_parent], hierarchy.cert_audit_element_levels()[audit]);
        maximum_audit_level[fine_parent] = std::max(
            maximum_audit_level[fine_parent], hierarchy.cert_audit_element_levels()[audit]);
    }
    int locally_unrefined_fine_element = -1;
    for (int element = 0; element < static_cast<int>(minimum_audit_level.size()); ++element) {
        if (minimum_audit_level[element] == hierarchy.fine_element_levels()[element]
            && maximum_audit_level[element] == hierarchy.fine_element_levels()[element]) {
            locally_unrefined_fine_element = element;
            break;
        }
    }
    require(locally_unrefined_fine_element >= 0,
            "test setup has no fine element without prior audit refinement");
    const std::size_t fine_elements_before = hierarchy.fine_mesh().elems.size();
    const std::uint64_t coarse_version_before = hierarchy.coarse_mesh_version();
    const std::uint64_t fine_version_before = hierarchy.fine_mesh_version();
    const std::uint64_t audit_version_before = hierarchy.cert_audit_mesh_version();
    const std::uint64_t interpolation_version_before = hierarchy.interpolation_version();
    const std::uint64_t boundary_version_before = hierarchy.boundary_version();
    const std::uint64_t corrector_version_before = hierarchy.corrector_space_version();
    const int selected_coarse_element =
        hierarchy.fine_parent_coarse_elements()[locally_unrefined_fine_element];
    const std::vector<int> covered_fine_elements =
        hierarchy.fine_elements_in_coarse_patch({selected_coarse_element});
    require(std::find(covered_fine_elements.begin(), covered_fine_elements.end(),
                      locally_unrefined_fine_element) != covered_fine_elements.end(),
            "coarse corrector patch omitted one of its fine descendants");
    hierarchy.refine_fine_in_coarse_patch({selected_coarse_element});
    require(hierarchy.fine_mesh().elems.size() > fine_elements_before,
            "local V_h refinement did not change the fine mesh");
    require(*std::min_element(
                hierarchy.fine_element_levels().begin(),
                hierarchy.fine_element_levels().end()) == hierarchy.fine_level(),
            "local V_h refinement unexpectedly refined the entire fine mesh");
    require(*std::max_element(
                hierarchy.fine_element_levels().begin(),
                hierarchy.fine_element_levels().end()) > hierarchy.fine_level(),
            "local V_h refinement did not create a finer level");
    require(hierarchy.coarse_mesh_version() == coarse_version_before,
            "local V_h refinement changed the coarse mesh version");
    require(hierarchy.fine_mesh_version() == fine_version_before + 1,
            "local V_h refinement did not advance its mesh version");
    require(hierarchy.cert_audit_mesh_version() > audit_version_before,
            "cert-audit mesh did not preserve nesting after local V_h refinement");
    require(hierarchy.interpolation_version() == interpolation_version_before + 1,
            "local V_h refinement did not invalidate interpolation caches");
    require(hierarchy.boundary_version() == boundary_version_before + 1,
            "local V_h refinement did not invalidate boundary caches");
    require(hierarchy.corrector_space_version() == corrector_version_before + 1,
            "local V_h refinement did not invalidate corrector-space caches");
    verify_three_level_prolongations(hierarchy);

    const std::vector<int> local_dofs = {
        0, static_cast<int>(hierarchy.fine_mesh().nodes.size() / 2)};
    const Eigen::SparseMatrix<double> local_constraints =
        hierarchy.fine_kernel_constraints(local_dofs);
    require(local_constraints.cols() == static_cast<int>(local_dofs.size()),
            "local kernel constraint restriction has the wrong column count");
    for (int column = 0; column < static_cast<int>(local_dofs.size()); ++column)
        require((local_constraints.col(column)
                 - hierarchy.fine_quasi_interpolation().col(local_dofs[column])).norm() < 1e-14,
                "local kernel constraint restriction selected the wrong column");

    std::vector<int> minimum_fine_level(
        hierarchy.coarse_mesh().elems.size(), std::numeric_limits<int>::max());
    for (int fine = 0; fine < static_cast<int>(hierarchy.fine_parent_coarse_elements().size());
         ++fine) {
        const int coarse_parent = hierarchy.fine_parent_coarse_elements()[fine];
        minimum_fine_level[coarse_parent] = std::min(
            minimum_fine_level[coarse_parent], hierarchy.fine_element_levels()[fine]);
    }
    int coarse_with_capacity = -1;
    for (int coarse = 0; coarse < static_cast<int>(minimum_fine_level.size()); ++coarse) {
        if (hierarchy.coarse_levels()[coarse] < minimum_fine_level[coarse]) {
            coarse_with_capacity = coarse;
            break;
        }
    }
    require(coarse_with_capacity >= 0,
            "test setup has no coarse element below its local V_h descendants");
    const std::vector<TriangleSignature> fine_before_second_coarse_refinement =
        canonical_triangles(hierarchy.fine_mesh());
    const std::vector<TriangleSignature> audit_before_second_coarse_refinement =
        canonical_triangles(hierarchy.cert_audit_mesh());
    hierarchy.refine({coarse_with_capacity});
    require(canonical_triangles(hierarchy.fine_mesh())
                == fine_before_second_coarse_refinement,
            "coarse refinement after local V_h refinement changed the fine mesh");
    require(canonical_triangles(hierarchy.cert_audit_mesh())
                == audit_before_second_coarse_refinement,
            "coarse refinement after local V_h refinement changed the audit mesh");
    verify_three_level_prolongations(hierarchy);

    const std::uint64_t second_audit_version = hierarchy.cert_audit_mesh_version();
    const std::uint64_t second_interpolation_version = hierarchy.interpolation_version();
    const std::uint64_t second_boundary_version = hierarchy.boundary_version();
    const std::uint64_t second_corrector_version = hierarchy.corrector_space_version();
    hierarchy.refine_cert_audit_from_fine_elements(
        {static_cast<int>(hierarchy.fine_mesh().elems.size() - 1)});
    require(hierarchy.cert_audit_mesh_version() == second_audit_version + 1,
            "audit refinement after H/h updates did not advance its mesh version");
    require(hierarchy.interpolation_version() == second_interpolation_version + 1,
            "audit refinement did not invalidate interpolation caches");
    require(hierarchy.boundary_version() == second_boundary_version + 1,
            "audit refinement did not invalidate boundary caches");
    require(hierarchy.corrector_space_version() == second_corrector_version,
            "audit-only refinement invalidated the corrector-space cache");
    verify_three_level_prolongations(hierarchy);

    config.h = 4;
    config.ell = 1;
    config.wavenumber = 3.0;
    HelmholtzLodModel model = HelmholtzLodModel::build_adaptive(
        config, hierarchy.coarse_mesh(), hierarchy.coarse_levels(),
        hierarchy.fine_mesh(), hierarchy.fine_element_levels());
    require(model.problem().coarse.elems.size() == hierarchy.coarse_mesh().elems.size(),
            "adaptive model changed the coarse mesh");
    require(model.problem().fine_level == 4, "adaptive model lost the fixed fine level");
    require(model.problem().max_fine_level > model.problem().fine_level,
            "adaptive model discarded local V_h levels");
    require(canonical_triangles(model.problem().fine)
                == canonical_triangles(hierarchy.fine_mesh()),
            "adaptive model reconstructed a different fine mesh");

    const ComplexFunction source = source_function();
    const ComplexVector load = assemble_helmholtz_load(model.problem().fine, source);
    const HelmholtzLodSolution solution = model.solve_load(load);
    const ComplexVector reference = model.solve_fine_reference(load);
    const ComplexVector error = reference - solution.fine_values;
    require(discrete_energy_norm(model.operators(), error) > 0.0,
            "adaptive test unexpectedly has zero LOD error");

    const diagnostics::HelmholtzResidualContributions residual =
        diagnostics::assemble_helmholtz_residual_contributions(
            model.problem(), model.operators(), solution.fine_values, load, source);
    std::cout << "Residual reconstruction difference: "
              << residual.algebraic_relative_difference << '\n';
    require(residual.algebraic_relative_difference < 1e-10,
            "integrated residual disagrees with b-Au");
    const diagnostics::HelmholtzIndicatorSet indicators = diagnostics::build_helmholtz_indicators(
        model.problem(), residual);
    require(indicators.fine > 0.0 && indicators.mixed > 0.0 && indicators.macro > 0.0,
            "adaptive residual indicators must be positive");
    require(indicators.fine_squared.size() == model.problem().coarse.elems.size(),
            "coarse indicator count is incorrect");

    const std::vector<double> dual = diagnostics::build_local_dual_indicators(
        model.problem(), model.operators(), residual, 1);
    require(dual.size() == model.problem().coarse.elems.size(),
            "local dual indicator count is incorrect");
    require(*std::max_element(dual.begin(), dual.end()) > 0.0,
            "local dual indicators are all zero");

    std::vector<char> eligible(hierarchy.coarse_levels().size(), true);
    const std::vector<int> marked = mark_doerfler(
        indicators.mixed_squared, 0.5, eligible);
    require(!marked.empty(), "Doerfler marking returned an empty set");
}

void verify_capacity_expansion() {
    AdaptiveMeshHierarchy hierarchy(make_helmholtz_unit_square_mesh(), 0, 1);
    verify_three_level_prolongations(hierarchy);

    for (int cycle = 0; cycle < 2; ++cycle) {
        std::vector<int> marked(hierarchy.coarse_mesh().elems.size());
        for (int element = 0; element < static_cast<int>(marked.size()); ++element)
            marked[element] = element;

        const std::size_t old_fine_elements = hierarchy.fine_mesh().elems.size();
        const std::size_t old_audit_elements = hierarchy.cert_audit_mesh().elems.size();
        const std::uint64_t old_fine_version = hierarchy.fine_mesh_version();
        const std::uint64_t old_audit_version = hierarchy.cert_audit_mesh_version();
        hierarchy.refine(marked);

        require(hierarchy.fine_mesh().elems.size() > old_fine_elements,
                "H capacity exhaustion did not expand the master fine mesh");
        require(hierarchy.cert_audit_mesh().elems.size() > old_audit_elements,
                "H capacity exhaustion did not expand the cert-audit mesh");
        require(hierarchy.fine_mesh_version() > old_fine_version,
                "H capacity expansion did not advance the fine mesh version");
        require(hierarchy.cert_audit_mesh_version() > old_audit_version,
                "H capacity expansion did not advance the audit mesh version");
        require(*std::max_element(
                    hierarchy.coarse_levels().begin(),
                    hierarchy.coarse_levels().end())
                    < *std::min_element(
                        hierarchy.fine_element_levels().begin(),
                        hierarchy.fine_element_levels().end()),
                "H capacity expansion did not restore a fine-level gap");
        require(std::abs(boundary_measure(
                    hierarchy.fine_mesh(), BoundaryTag::Robin) - 4.0) < 2e-12,
                "H capacity expansion changed the fine boundary contract");
        require(std::abs(boundary_measure(
                    hierarchy.cert_audit_mesh(), BoundaryTag::Robin) - 4.0) < 2e-12,
                "H capacity expansion changed the audit boundary contract");
        verify_three_level_prolongations(hierarchy);
    }

    AdaptiveMeshHierarchy mixed(make_helmholtz_l_shape_mesh(), 0, 1);
    std::vector<int> mixed_marked(mixed.coarse_mesh().elems.size());
    for (int element = 0; element < static_cast<int>(mixed_marked.size()); ++element)
        mixed_marked[element] = element;
    mixed.refine(mixed_marked);
    require(std::abs(boundary_measure(
                mixed.fine_mesh(), BoundaryTag::Dirichlet) - 2.0) < 2e-12,
            "mixed H capacity expansion changed the Dirichlet boundary");
    require(std::abs(boundary_measure(
                mixed.cert_audit_mesh(), BoundaryTag::Robin) - 6.0) < 2e-12,
            "mixed H capacity expansion changed the Robin boundary");
    verify_three_level_prolongations(mixed);
}

void verify_doerfler_minimal_set() {
    const std::vector<int> marked = mark_doerfler({4.0, 3.0, 2.0, 1.0}, 0.6);
    require(marked == std::vector<int>({0, 1}),
            "Doerfler marking is not deterministic or minimal");
}

} // namespace

int main() {
    try {
        verify_hierarchy_and_estimator();
        verify_capacity_expansion();
        verify_doerfler_minimal_set();
        std::cout << "Adaptive Helmholtz hierarchy and residual estimators passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_adaptive failed: " << error.what() << '\n';
        return 1;
    }
}

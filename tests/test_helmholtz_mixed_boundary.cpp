#include "helmholtz/boundary.h"
#include "helmholtz/adaptive/estimator.h"
#include "helmholtz/model.h"
#include "helmholtz/operators.h"
#include "mesh/refine.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace lod2d;
using namespace lod2d::helmholtz;

namespace {

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

void verify_l_shape_and_refinement() {
    const TriMesh initial = make_helmholtz_l_shape_mesh();
    double area = 0.0;
    for (double value : compute_area(initial)) area += value;
    require(std::abs(area - 3.0) < 1e-14, "L-shaped mesh area is not three");
    require(initial.dirichlet == std::vector<int>({1, 3, 4}),
            "L-shaped reentrant boundary has the wrong Dirichlet nodes");
    require(std::abs(boundary_measure(initial, BoundaryTag::Dirichlet) - 2.0) < 1e-14,
            "L-shaped Dirichlet boundary measure is wrong");
    require(std::abs(boundary_measure(initial, BoundaryTag::Robin) - 6.0) < 1e-14,
            "L-shaped Robin boundary measure is wrong");

    TriMesh mesh = initial;
    for (int level = 0; level < 4; ++level) {
        mesh = refine_nvb(mesh).mesh;
        validate_boundary_tags(mesh);
        require(std::abs(boundary_measure(mesh, BoundaryTag::Dirichlet) - 2.0) < 2e-12,
                "NVB changed the Dirichlet boundary measure");
        require(std::abs(boundary_measure(mesh, BoundaryTag::Robin) - 6.0) < 2e-12,
                "NVB changed the Robin boundary measure");
        for (const Edge &edge : boundary_edges_with_tag(mesh, BoundaryTag::Dirichlet)) {
            const Point2 midpoint = 0.5 * (mesh.nodes[edge[0]] + mesh.nodes[edge[1]]);
            require((std::abs(midpoint.x()) < 1e-12 && midpoint.y() <= 0.0)
                        || (std::abs(midpoint.y()) < 1e-12 && midpoint.x() >= 0.0),
                    "NVB moved a Dirichlet tag away from a reentrant edge");
        }
    }

    mesh = initial;
    for (int step = 0; step < 6; ++step) {
        const int marked = step % static_cast<int>(mesh.elems.size());
        mesh = bisect_newest_vertex(mesh, {marked}).mesh;
        validate_boundary_tags(mesh);
        require(std::abs(boundary_measure(mesh, BoundaryTag::Dirichlet) - 2.0) < 2e-12,
                "local NVB closure changed the Dirichlet boundary measure");
        require(std::abs(boundary_measure(mesh, BoundaryTag::Robin) - 6.0) < 2e-12,
                "local NVB closure changed the Robin boundary measure");
    }
}

TriMesh mixed_unit_square() {
    TriMesh mesh = make_helmholtz_unit_square_mesh();
    mesh.boundary_edges = {
        {{0, 1}, BoundaryTag::Dirichlet},
        {{1, 2}, BoundaryTag::Robin},
        {{2, 3}, BoundaryTag::Robin},
        {{0, 3}, BoundaryTag::Neumann}};
    synchronize_dirichlet_nodes(mesh);
    validate_boundary_tags(mesh);
    return mesh;
}

void verify_mixed_global_operator_and_solve() {
    const TriMesh mesh = mixed_unit_square();
    const HelmholtzOperators operators = assemble_helmholtz_operators(mesh, 2.0);
    require(operators.dirichlet_nodes == std::vector<int>({0, 1}),
            "operator lost the Dirichlet nodes");
    require(std::abs(boundary_measure(mesh, BoundaryTag::Neumann) - 1.0) < 1e-14,
            "mixed mesh lost its homogeneous Neumann edge");

    Eigen::Matrix4d expected = Eigen::Matrix4d::Zero();
    expected(1, 1) = 1.0 / 3.0;
    expected(2, 2) = 2.0 / 3.0;
    expected(3, 3) = 1.0 / 3.0;
    expected(1, 2) = expected(2, 1) = 1.0 / 6.0;
    expected(2, 3) = expected(3, 2) = 1.0 / 6.0;
    require((Eigen::MatrixXd(operators.boundary_mass) - expected).norm() < 1e-14,
            "mixed Robin matrix differs from the hand reference");

    ComplexVector load(4);
    load << Complex(4.0, 1.0), Complex(-3.0, 2.0),
            Complex(1.0, -0.5), Complex(2.0, 0.25);
    const ComplexVector solution = solve_helmholtz_fem(operators, load);
    require(std::abs(solution(0)) < 1e-14 && std::abs(solution(1)) < 1e-14,
            "homogeneous Dirichlet values are not zero");
    const ComplexVector residual = operators.system * solution - load;
    require(std::abs(residual(2)) < 1e-12 && std::abs(residual(3)) < 1e-12,
            "mixed FEM residual is not zero on free nodes");
}

void verify_mixed_lod_chain() {
    HelmholtzProblemConfig config;
    config.initial_mesh = mixed_unit_square();
    config.H = 1;
    config.h = 3;
    config.ell = 1;
    config.wavenumber = 4.0;
    HelmholtzLodModel model = HelmholtzLodModel::build(config);
    require(!model.problem().coarse.dirichlet.empty()
                && !model.problem().fine.dirichlet.empty(),
            "mixed LOD hierarchy lost Dirichlet nodes");
    const ComplexFunction source = [](const Point2 &point) {
        return Complex(1.0 + point.x(), -0.25 * point.y());
    };
    const HelmholtzLodSolution solution = model.solve_source(source);
    for (int node : model.problem().fine.dirichlet)
        require(std::abs(solution.fine_values(node)) < 1e-12,
                "mixed LOD solution violates homogeneous Dirichlet data");
    require(solution.petrov_residual < 1e-10,
            "mixed LOD Petrov residual is too large");
    const ComplexVector lod_load = assemble_helmholtz_load(model.problem().fine, source);
    const adaptive::diagnostics::HelmholtzResidualContributions residual =
        adaptive::diagnostics::assemble_helmholtz_residual_contributions(
            model.problem(), model.operators(), solution.fine_values, lod_load, source);
    require(residual.algebraic_relative_difference < 1e-10,
            "mixed residual reconstruction disagrees with the free-node algebraic residual");
    for (const auto &edge : residual.edges) {
        const BoundaryTag tag = boundary_tag(model.problem().fine, edge.nodes);
        if (tag == BoundaryTag::Dirichlet) {
            require(!edge.robin_boundary && !edge.neumann_boundary
                        && edge.residual_l2_squared == 0.0,
                    "Dirichlet edge was assembled as a natural residual edge");
        } else if (tag == BoundaryTag::Neumann) {
            require(edge.neumann_boundary && !edge.robin_boundary,
                    "homogeneous Neumann residual edge was not identified");
        }
    }
    const ComplexVector conforming = model.solve_fine_reference(lod_load);
    const adaptive::diagnostics::HelmholtzP1ResidualEstimate afem =
        adaptive::diagnostics::estimate_conforming_p1_residual(
            model.problem().fine, model.operators(), conforming,
            lod_load, source);
    require(afem.element_squared.size() == model.problem().fine.elems.size()
                && afem.eta > 0.0
                && afem.algebraic_relative_difference < 1e-10,
            "mixed-boundary conforming P1 AFEM estimate is invalid");
    require(*std::max_element(
                afem.body_squared.begin(), afem.body_squared.end()) > 0.0
                && *std::max_element(
                    afem.interior_jump_squared.begin(),
                    afem.interior_jump_squared.end()) > 0.0
                && *std::max_element(
                    afem.neumann_boundary_squared.begin(),
                    afem.neumann_boundary_squared.end()) > 0.0
                && *std::max_element(
                    afem.robin_boundary_squared.begin(),
                    afem.robin_boundary_squared.end()) > 0.0,
            "AFEM estimate omitted body, jump, Neumann, or impedance terms");
    require(!adaptive::mark_doerfler(afem.element_squared, 0.5).empty(),
            "AFEM residual did not produce a Doerfler marking");
    const ComplexVector load = assemble_helmholtz_load(
        model.problem().fine,
        [](const Point2 &) { return Complex(1.0, 0.0); });
    const ComplexVector reference = model.solve_fine_reference(load);
    for (int node : model.problem().fine.dirichlet)
        require(std::abs(reference(node)) < 1e-13,
                "mixed fine reference violates homogeneous Dirichlet data");
}

void verify_invalid_classification_is_rejected() {
    TriMesh mesh = mixed_unit_square();
    mesh.boundary_edges.pop_back();
    bool rejected = false;
    try {
        validate_boundary_tags(mesh);
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    require(rejected, "incomplete boundary classification was accepted");
}

} // namespace

int main() {
    try {
        verify_l_shape_and_refinement();
        verify_mixed_global_operator_and_solve();
        verify_mixed_lod_chain();
        verify_invalid_classification_is_rejected();
        std::cout << "Helmholtz mixed-boundary mesh and FEM passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_mixed_boundary failed: " << error.what() << '\n';
        return 1;
    }
}

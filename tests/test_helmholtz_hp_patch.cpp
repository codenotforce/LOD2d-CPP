#include "helmholtz/hp_interpolation.h"
#include "helmholtz/hp_operators.h"
#include "helmholtz/hp_patch.h"
#include "helmholtz/model.h"
#include "helmholtz/operators.h"
#include "mesh/refine.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace lod2d;
using namespace lod2d::helmholtz;

namespace {

void require(bool condition, const std::string &message) {
    if (!condition) throw std::runtime_error(message);
}

double relative_norm(
    const Eigen::MatrixXcd &residual,
    const Eigen::MatrixXcd &reference) {
    return residual.norm() / std::max(1.0, reference.norm());
}

Eigen::Vector3d barycentric(
    const TriMesh &mesh,
    const Triangle &triangle,
    const Point2 &point) {
    Eigen::Matrix2d jacobian;
    jacobian.col(0) =
        mesh.nodes[triangle[1]] - mesh.nodes[triangle[0]];
    jacobian.col(1) =
        mesh.nodes[triangle[2]] - mesh.nodes[triangle[0]];
    const Eigen::Vector2d reference =
        jacobian.inverse() * (point - mesh.nodes[triangle[0]]);
    return {1.0 - reference.x() - reference.y(),
            reference.x(), reference.y()};
}

std::vector<int> target_types(
    const TriMesh &coarse,
    const HelmholtzHpPatchAssembler &assembler) {
    auto is_corner = [&](int vertex) {
        const Point2 &point = coarse.nodes[vertex];
        return (std::abs(point.x()) < 1e-13
                || std::abs(point.x() - 1.0) < 1e-13)
            && (std::abs(point.y()) < 1e-13
                || std::abs(point.y() - 1.0) < 1e-13);
    };

    int interior = -1;
    int boundary_target = -1;
    int corner = -1;
    for (int element = 0;
         element < static_cast<int>(coarse.elems.size()); ++element) {
        bool has_corner = false;
        for (int vertex : coarse.elems[element]) {
            has_corner = has_corner || is_corner(vertex);
        }
        const bool touches_boundary =
            assembler.assemble(element).touches_physical_boundary;
        if (has_corner && corner < 0) corner = element;
        if (!has_corner && touches_boundary && boundary_target < 0)
            boundary_target = element;
        if (!touches_boundary && interior < 0) interior = element;
    }
    require(interior >= 0, "failed to find an interior coarse target");
    require(boundary_target >= 0, "failed to find a boundary coarse target");
    require(corner >= 0, "failed to find a corner coarse target");
    return {interior, boundary_target, corner};
}

void check_basis_and_space(const TriMesh &mesh) {
    const auto [edges, boundary] = compute_edges(mesh);
    (void)boundary;
    for (int degree = 1; degree <= 3; ++degree) {
        HpTriSpace space(mesh, degree);
        const int expected =
            static_cast<int>(mesh.nodes.size())
            + (degree - 1) * static_cast<int>(edges.size())
            + ((degree - 1) * (degree - 2) / 2)
                * static_cast<int>(mesh.elems.size());
        require(space.dof_count() == expected, "hp global DOF count is wrong");
        require(space.local_dof_count() == (degree + 1) * (degree + 2) / 2,
                "hp local DOF count is wrong");
        for (int node = 0; node < space.local_dof_count(); ++node) {
            const auto evaluation =
                space.basis().evaluate(space.basis().nodes()[node]);
            Eigen::VectorXd expected_values =
                Eigen::VectorXd::Zero(space.local_dof_count());
            expected_values(node) = 1.0;
            require((evaluation.values - expected_values).norm() < 2e-12,
                    "hp basis is not nodal");
        }
        const auto sample =
            space.basis().evaluate(Eigen::Vector3d(0.2, 0.3, 0.5));
        require(std::abs(sample.values.sum() - 1.0) < 2e-12,
                "hp basis does not form a partition of unity");
        require(sample.reference_gradients.colwise().sum().norm() < 2e-12,
                "hp basis gradients do not sum to zero");

        Eigen::VectorXd coefficients(space.dof_count());
        for (int dof = 0; dof < space.dof_count(); ++dof)
            coefficients(dof) = std::sin(0.37 * (dof + 1));
        bool checked_shared_edge = false;
        for (int edge = 0;
             edge < static_cast<int>(edges.size()) && !checked_shared_edge;
             ++edge) {
            if (boundary[edge]) continue;
            const Point2 point =
                0.37 * mesh.nodes[edges[edge][0]]
                + 0.63 * mesh.nodes[edges[edge][1]];
            std::vector<double> traces;
            for (int element = 0;
                 element < static_cast<int>(mesh.elems.size()); ++element) {
                const Triangle &triangle = mesh.elems[element];
                const bool has_first =
                    std::find(triangle.begin(), triangle.end(), edges[edge][0])
                    != triangle.end();
                const bool has_second =
                    std::find(triangle.begin(), triangle.end(), edges[edge][1])
                    != triangle.end();
                if (!has_first || !has_second) continue;
                const auto evaluation = space.basis().evaluate(
                    barycentric(mesh, triangle, point));
                double trace = 0.0;
                for (int local = 0;
                     local < space.local_dof_count(); ++local)
                    trace += coefficients(
                        space.element_dofs()[element][local])
                        * evaluation.values(local);
                traces.push_back(trace);
            }
            require(traces.size() == 2,
                    "interior edge does not have two traces");
            require(std::abs(traces[0] - traces[1]) < 2e-12,
                    "hp shared-edge traces are discontinuous");
            checked_shared_edge = true;
        }
        require(checked_shared_edge, "no interior edge was checked");
    }
}

void check_p1_compatibility(const HelmholtzProblemData &problem) {
    HpTriSpace space(problem.fine, 1);
    const auto hp_operators =
        assemble_helmholtz_hp_operators(space, 2.5);
    const auto p1_operators =
        assemble_helmholtz_operators(problem.fine, 2.5);
    require((hp_operators.stiffness - p1_operators.stiffness).norm() < 2e-11,
            "hp p=1 stiffness differs from the P1 implementation");
    require((hp_operators.mass - p1_operators.mass).norm() < 2e-11,
            "hp p=1 mass differs from the P1 implementation");
    require((hp_operators.boundary_mass
             - p1_operators.boundary_mass).norm() < 2e-11,
            "hp p=1 Robin matrix differs from the P1 implementation");

    const auto interpolation = build_helmholtz_hp_interpolation(
        problem.coarse, space, problem.fine_element_prolongation);
    require((interpolation.coarse_injection
             - problem.coarse_to_fine).norm() < 2e-11,
            "hp p=1 coarse injection differs from nodal prolongation");
    require((interpolation.quasi_interpolation
             - problem.quasi_interpolation).norm() < 2e-10,
            "hp p=1 quasi-interpolation differs from the P1 implementation");
}

void check_hp_interpolation(const HelmholtzProblemData &problem) {
    for (int degree : {2, 3}) {
        HpTriSpace space(problem.fine, degree);
        const auto interpolation = build_helmholtz_hp_interpolation(
            problem.coarse, space, problem.fine_element_prolongation);
        const Eigen::MatrixXd reproduction = Eigen::MatrixXd(
            interpolation.quasi_interpolation
            * interpolation.coarse_injection);
        const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(
            problem.coarse.nodes.size(), problem.coarse.nodes.size());
        require((reproduction - identity).norm() < 2e-10,
                "hp I_H does not reproduce coarse P1");
    }
}

void check_patch_solves(const HelmholtzProblemData &problem) {
    HpTriSpace space(problem.fine, 2);
    const auto interpolation = build_helmholtz_hp_interpolation(
        problem.coarse, space, problem.fine_element_prolongation);
    const auto operators =
        assemble_helmholtz_hp_operators(space, 2.5);
    HelmholtzHpPatchAssembler assembler(
        problem.coarse, space, problem.fine_element_prolongation,
        problem.patches, interpolation, operators);

    const auto targets = target_types(problem.coarse, assembler);
    for (int target : targets) {
        const auto result = assembler.solve_direct_saddle(target);
        const auto &system = result.system;
        const auto &primal = result.primal;
        const Eigen::MatrixXcd primal_residual =
            system.helmholtz * primal.corrector
            + system.constraints.transpose().cast<Complex>()
                * primal.multipliers
            - system.rhs;
        const Eigen::MatrixXcd constraint_residual =
            system.constraints.cast<Complex>() * primal.corrector;
        const Eigen::MatrixXcd adjoint_residual =
            system.helmholtz.adjoint() * result.adjoint_corrector
            + system.constraints.transpose().cast<Complex>()
                * result.adjoint_multipliers
            - system.rhs.conjugate();
        require(relative_norm(primal_residual, system.rhs) < 2e-10,
                "hp direct-saddle primal residual is too large");
        require(relative_norm(constraint_residual, primal.corrector) < 2e-10,
                "hp direct-saddle constraint residual is too large");
        require(relative_norm(adjoint_residual, system.rhs) < 2e-10,
                "hp direct-saddle adjoint residual is too large");

        std::vector<int> patch_incidence(space.dof_count(), 0);
        for (int element : system.patch_elements)
            for (int dof : space.element_dofs()[element])
                ++patch_incidence[dof];
        for (int dof : system.local_vertices)
            require(patch_incidence[dof] == space.dof_incidence()[dof],
                    "artificial-boundary hp DOF was not removed");
    }
    require(!assembler.assemble(targets[0]).touches_physical_boundary,
            "interior hp patch unexpectedly touches the physical boundary");
    require(assembler.assemble(targets[1]).touches_physical_boundary,
            "boundary hp patch lost its physical Robin boundary");
    require(assembler.assemble(targets[2]).touches_physical_boundary,
            "corner hp patch lost its physical Robin boundary");
}

} // namespace

int main() {
    try {
        const HelmholtzProblemData problem =
            build_helmholtz_problem_data(
                make_helmholtz_unit_square_mesh(), 5, 7, 1);
        check_basis_and_space(problem.fine);
        check_p1_compatibility(problem);
        check_hp_interpolation(problem);
        check_patch_solves(problem);
        std::cout << "Helmholtz hp patch: p=1,2,3 and DirectSaddle passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_hp_patch failed: "
                  << error.what() << '\n';
        return 1;
    }
}

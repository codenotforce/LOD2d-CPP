#include "helmholtz/hp_operators.h"

#include <Eigen/SparseLU>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace lod2d::helmholtz {

ComplexVector assemble_helmholtz_hp_load(
    const HpTriSpace &space,
    const ComplexFunction &source) {
    if (!source)
        throw std::invalid_argument("Helmholtz hp source function is empty");
    const TriMesh &mesh = space.mesh();
    ComplexVector load = ComplexVector::Zero(space.dof_count());
    const auto quadrature =
        triangle_gauss_quadrature(space.degree() + 5);
    for (int element = 0;
         element < static_cast<int>(mesh.elems.size()); ++element) {
        const Triangle &triangle = mesh.elems[element];
        Eigen::Matrix2d jacobian;
        jacobian.col(0) =
            mesh.nodes[triangle[1]] - mesh.nodes[triangle[0]];
        jacobian.col(1) =
            mesh.nodes[triangle[2]] - mesh.nodes[triangle[0]];
        const double determinant = std::abs(jacobian.determinant());
        const auto &dofs = space.element_dofs()[element];
        for (const auto &point : quadrature) {
            Point2 physical = Point2::Zero();
            for (int i = 0; i < 3; ++i)
                physical +=
                    point.barycentric(i) * mesh.nodes[triangle[i]];
            const auto evaluation =
                space.basis().evaluate(point.barycentric);
            const Complex value = source(physical);
            for (int local = 0; local < space.local_dof_count(); ++local)
                load(dofs[local]) += determinant * point.weight
                    * value * evaluation.values(local);
        }
    }
    return load;
}

ComplexVector solve_helmholtz_hp_fem(
    const HelmholtzHpOperators &operators,
    const ComplexVector &load) {
    if (operators.system.rows() != load.size())
        throw std::invalid_argument(
            "Helmholtz hp load size does not match the operator");
    Eigen::SparseLU<ComplexSparseMatrix> solver;
    solver.analyzePattern(operators.system);
    solver.factorize(operators.system);
    if (solver.info() != Eigen::Success)
        throw std::runtime_error(
            "Helmholtz hp sparse LU factorization failed");
    ComplexVector solution = solver.solve(load);
    if (solver.info() != Eigen::Success || !solution.allFinite())
        throw std::runtime_error("Helmholtz hp sparse LU solve failed");
    return solution;
}

HelmholtzError compute_helmholtz_hp_error(
    const HpTriSpace &space,
    const ComplexVector &solution,
    double wavenumber,
    const ComplexFunction &exact,
    const ComplexGradientFunction &exact_gradient) {
    if (solution.size() != space.dof_count())
        throw std::invalid_argument(
            "Helmholtz hp solution size does not match the space");
    if (!exact || !exact_gradient)
        throw std::invalid_argument(
            "Helmholtz hp exact solution callbacks must not be empty");
    const TriMesh &mesh = space.mesh();
    const auto quadrature =
        triangle_gauss_quadrature(space.degree() + 5);
    double l2_squared = 0.0;
    double gradient_squared = 0.0;
    for (int element = 0;
         element < static_cast<int>(mesh.elems.size()); ++element) {
        const Triangle &triangle = mesh.elems[element];
        Eigen::Matrix2d jacobian;
        jacobian.col(0) =
            mesh.nodes[triangle[1]] - mesh.nodes[triangle[0]];
        jacobian.col(1) =
            mesh.nodes[triangle[2]] - mesh.nodes[triangle[0]];
        const double determinant = std::abs(jacobian.determinant());
        const Eigen::Matrix2d inverse_transpose =
            jacobian.inverse().transpose();
        const auto &dofs = space.element_dofs()[element];
        for (const auto &point : quadrature) {
            Point2 physical = Point2::Zero();
            for (int i = 0; i < 3; ++i)
                physical +=
                    point.barycentric(i) * mesh.nodes[triangle[i]];
            const auto evaluation =
                space.basis().evaluate(point.barycentric);
            Complex discrete_value = 0.0;
            Eigen::Vector2cd discrete_gradient =
                Eigen::Vector2cd::Zero();
            for (int local = 0; local < space.local_dof_count(); ++local) {
                discrete_value +=
                    solution(dofs[local]) * evaluation.values(local);
                const Eigen::Vector2d gradient =
                    inverse_transpose
                    * evaluation.reference_gradients.row(local).transpose();
                discrete_gradient +=
                    solution(dofs[local]) * gradient.cast<Complex>();
            }
            const double weight = determinant * point.weight;
            l2_squared += weight
                * std::norm(exact(physical) - discrete_value);
            gradient_squared += weight
                * (exact_gradient(physical) - discrete_gradient).squaredNorm();
        }
    }
    HelmholtzError error;
    error.l2 = std::sqrt(std::max(0.0, l2_squared));
    error.energy = std::sqrt(std::max(
        0.0, gradient_squared + wavenumber * wavenumber * l2_squared));
    return error;
}

} // namespace lod2d::helmholtz

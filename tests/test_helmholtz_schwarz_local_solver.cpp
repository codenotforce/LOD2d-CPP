#include "helmholtz/model.h"
#include "helmholtz/patch_multigrid.h"
#include "helmholtz/schwarz_local_solver.h"
#include "helmholtz/schwarz_patch.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace lod2d::helmholtz;

namespace {

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

} // namespace

int main() {
    try {
        HelmholtzProblemConfig problem_config;
        problem_config.H = 3;
        problem_config.h = 6;
        problem_config.ell = 1;
        problem_config.wavenumber = 2.0;
        HelmholtzLodModel model = HelmholtzLodModel::build(problem_config);
        const HelmholtzProblemData &problem = model.problem();
        HelmholtzSchwarzPatchAssembler assembler(
            problem.fine,
            problem.fine_element_prolongation,
            problem.patches,
            static_cast<int>(problem.coarse.elems.size()),
            problem.fine_hierarchy_meshes,
            problem.fine_node_level_prolongations,
            problem.fine_element_level_prolongations,
            model.operators());
        const HelmholtzSchwarzLocalSystem system = assembler.assemble(
            0, HelmholtzSchwarzArtificialBoundary::HomogeneousDirichlet,
            1.0, true, true);
        require(!system.geometric_prolongations.empty(),
                "Dirichlet patch hierarchy is empty");
        require(system.geometric_prolongations.back().rows()
                    == system.matrix.rows(),
                "patch hierarchy does not end at local operator DOFs");

        ComplexVector right_hand_side(system.matrix.rows());
        for (int index = 0; index < right_hand_side.size(); ++index) {
            right_hand_side(index) = Complex(
                std::sin(0.31 * (index + 1)),
                std::cos(0.17 * (index + 1)));
        }

        HelmholtzSchwarzLocalSolver direct(
            system.matrix,
            system.mass,
            model.operators().wavenumber);
        const HelmholtzSchwarzLocalSolveResult reference =
            direct.solve(right_hand_side);

        ComplexMatrix block_rhs(system.matrix.rows(), 2);
        block_rhs.col(0) = right_hand_side;
        block_rhs.col(1) = right_hand_side.reverse();
        const ComplexMatrix block_solution = direct.solve_block(block_rhs);
        require(
            (block_solution.col(0) - reference.solution).norm()
                    / std::max(1.0, reference.solution.norm())
                < 1e-13,
            "SparseLU block solve first column differs from scalar solve");
        const HelmholtzSchwarzLocalSolveResult second_reference =
            direct.solve(block_rhs.col(1));
        require(
            (block_solution.col(1) - second_reference.solution).norm()
                    / std::max(1.0, second_reference.solution.norm())
                < 1e-13,
            "SparseLU block solve second column differs from scalar solve");

        HelmholtzSchwarzLocalSolverConfig exact_shift_config;
        exact_shift_config.kind =
            HelmholtzSchwarzLocalSolverKind::ShiftedGmres;
        exact_shift_config.shift_alpha = 0.0;
        exact_shift_config.gmres.restart = 10;
        exact_shift_config.gmres.max_iterations = 30;
        exact_shift_config.gmres.relative_tolerance = 1e-12;
        HelmholtzSchwarzLocalSolver exact_shift(
            system.matrix,
            system.mass,
            model.operators().wavenumber,
            exact_shift_config);
        const HelmholtzSchwarzLocalSolveResult one_step =
            exact_shift.solve(right_hand_side);
        require(one_step.iterations == 1,
                "zero-shift exact inverse did not converge in one iteration");
        require(one_step.relative_residual < 1e-12,
                "zero-shift exact inverse residual is too large");
        require(
            (one_step.solution - reference.solution).norm()
                    / std::max(1.0, reference.solution.norm())
                < 1e-11,
            "zero-shift GMRES differs from direct local solve");

        HelmholtzSchwarzLocalSolverConfig shifted_config =
            exact_shift_config;
        shifted_config.shift_alpha = 0.2;
        HelmholtzSchwarzLocalSolver shifted(
            system.matrix,
            system.mass,
            model.operators().wavenumber,
            shifted_config);
        const HelmholtzSchwarzLocalSolveResult shifted_result =
            shifted.solve(right_hand_side);
        require(shifted_result.iterations > 1,
                "nonzero shift did not exercise iterative correction");
        require(shifted_result.relative_residual < 1e-11,
                "shifted local GMRES residual is too large");
        require(
            (shifted_result.solution - reference.solution).norm()
                    / std::max(1.0, reference.solution.norm())
                < 1e-10,
            "shifted local GMRES differs from direct local solve");

        HelmholtzSchwarzLocalSolverConfig vcycle_config = shifted_config;
        vcycle_config.shifted_inverse =
            HelmholtzSchwarzShiftedInverseKind::GeometricVcycle;
        vcycle_config.gmres.restart = 30;
        vcycle_config.gmres.max_iterations = 200;
        HelmholtzSchwarzLocalSolver vcycle_solver(
            system.matrix,
            system.mass,
            system.geometric_prolongations,
            model.operators().wavenumber,
            vcycle_config);
        const HelmholtzSchwarzLocalSolveResult vcycle_result =
            vcycle_solver.solve(right_hand_side);
        require(vcycle_solver.vcycle_levels() >= 2,
                "Schwarz V-cycle did not retain a hierarchy");
        require(vcycle_solver.vcycle_finest_dofs() == system.matrix.rows(),
                "Schwarz V-cycle finest dimension is wrong");
        require(vcycle_result.relative_residual < 1e-11,
                "V-cycle-preconditioned local GMRES residual is too large");
        require(
            (vcycle_result.solution - reference.solution).norm()
                    / std::max(1.0, reference.solution.norm())
                < 1e-10,
            "V-cycle-preconditioned local GMRES differs from direct solve");

        const double epsilon = vcycle_config.shift_alpha
            * model.operators().wavenumber * model.operators().wavenumber;
        ComplexSparseMatrix shifted_operator = system.matrix
            - Complex(0.0, epsilon) * system.mass.cast<Complex>();
        HelmholtzPatchVcycle raw_vcycle(
            shifted_operator,
            system.geometric_prolongations,
            vcycle_config.vcycle_pre_smooth,
            vcycle_config.vcycle_post_smooth,
            vcycle_config.vcycle_coarse_max_dofs,
            vcycle_config.vcycle_jacobi_weight);
        const ComplexVector second = right_hand_side.reverse();
        const Complex a(0.37, -0.21);
        const Complex b(-0.14, 0.42);
        const ComplexVector combined =
            raw_vcycle.apply(a * right_hand_side + b * second);
        const ComplexVector separate =
            a * raw_vcycle.apply(right_hand_side)
            + b * raw_vcycle.apply(second);
        require((combined - separate).norm()
                    / std::max(1.0, combined.norm()) < 1e-12,
                "fixed geometric V-cycle is not linear");

        bool rejected_impedance_hierarchy = false;
        try {
            (void)assembler.assemble(
                0, HelmholtzSchwarzArtificialBoundary::Impedance,
                1.0, true, true);
        } catch (const std::invalid_argument &) {
            rejected_impedance_hierarchy = true;
        }
        require(rejected_impedance_hierarchy,
                "V-cycle hierarchy accepted impedance patch boundaries");

        const ComplexVector zero = ComplexVector::Zero(system.matrix.rows());
        const HelmholtzSchwarzLocalSolveResult zero_result =
            shifted.solve(zero);
        require(zero_result.iterations == 0 && zero_result.solution.norm() == 0.0,
                "shifted local GMRES mishandled a zero right-hand side");

        std::cout << "local_dofs=" << system.matrix.rows()
                  << " alpha0_it=" << one_step.iterations
                  << " alpha02_it=" << shifted_result.iterations
                  << " alpha02_residual="
                  << shifted_result.relative_residual
                  << " vcycle_levels=" << vcycle_solver.vcycle_levels()
                  << " vcycle_it=" << vcycle_result.iterations
                  << " vcycle_residual="
                  << vcycle_result.relative_residual << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_schwarz_local_solver failed: "
                  << error.what() << '\n';
        return 1;
    }
}

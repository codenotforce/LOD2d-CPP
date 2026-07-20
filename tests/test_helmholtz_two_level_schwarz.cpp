#include "helmholtz/two_level_schwarz.h"
#include "solver/right_gmres.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace lod2d::helmholtz;
namespace solver = lod2d::solver;

namespace {

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

double energy_norm(
    const HelmholtzOperators &operators,
    const ComplexVector &vector) {
    const ComplexVector weighted =
        operators.stiffness.cast<Complex>() * vector
        + operators.wavenumber * operators.wavenumber
            * (operators.mass.cast<Complex>() * vector);
    return std::sqrt(std::max(0.0, std::real(vector.dot(weighted))));
}

void check_mode(
    const HelmholtzLodModel &model,
    const HelmholtzTwoLevelSchwarzPreconditioner &preconditioner,
    const ComplexVector &right_hand_side,
    const ComplexVector &reference,
    HelmholtzTwoLevelSchwarzMode mode,
    const char *name) {
    const ComplexSparseMatrix &matrix = model.operators().system;
    solver::RightGmresConfig config;
    config.restart = matrix.rows();
    config.max_iterations = 500;
    config.relative_tolerance = 1e-10;

    const auto apply_operator = [&](const solver::ComplexVector &vector) {
        return solver::ComplexVector(matrix * vector);
    };
    const auto apply_preconditioner = [&](const solver::ComplexVector &vector) {
        return solver::ComplexVector(preconditioner.apply(vector, mode));
    };
    const bool flexible =
        preconditioner.config().local_solver.kind
        == HelmholtzSchwarzLocalSolverKind::ShiftedGmres;
    const solver::RightGmresResult result = flexible
        ? solver::solve_right_preconditioned_fgmres(
            matrix.rows(), apply_operator, apply_preconditioner,
            right_hand_side, config)
        : solver::solve_right_preconditioned_gmres(
            matrix.rows(), apply_operator, apply_preconditioner,
            right_hand_side, config);
    require(result.converged, "two-level Schwarz GMRES did not converge");
    require(result.relative_residual < 1e-9,
            "two-level Schwarz true residual is too large");

    const double reference_energy =
        std::max(1e-30, energy_norm(model.operators(), reference));
    const double relative_energy =
        energy_norm(model.operators(), result.solution - reference)
        / reference_energy;
    require(relative_energy < 1e-8,
            "two-level Schwarz disagrees with the fine SparseLU reference");
    std::cout << name << ": iterations=" << result.iterations
              << " residual=" << result.relative_residual
              << " energy_rel=" << relative_energy << '\n';
}

} // namespace

int main() {
    try {
        HelmholtzProblemConfig config;
        config.H = 3;
        config.h = 6;
        config.ell = 1;
        config.wavenumber = 2.0;
        HelmholtzLodModel model = HelmholtzLodModel::build(config);
        HelmholtzTwoLevelSchwarzPreconditioner preconditioner(model);

        const auto &diagnostics = preconditioner.diagnostics();
        require(
            diagnostics.subdomains
                == static_cast<int>(model.problem().coarse.elems.size()),
            "two-level Schwarz did not create one subdomain per element patch");
        require(diagnostics.min_local_dofs > 0,
                "two-level Schwarz contains an empty subdomain");
        require(diagnostics.max_local_dofs < model.operators().system.rows(),
                "two-level Schwarz test did not exercise proper subdomains");
        require(diagnostics.uncovered_dofs == 0,
                "two-level Schwarz patches do not cover the fine space");
        require(diagnostics.partition_unity_error < 1e-14,
                "two-level Schwarz weights do not form a partition of unity");

        ComplexVector right_hand_side(model.operators().system.rows());
        for (int index = 0; index < right_hand_side.size(); ++index) {
            right_hand_side(index) = Complex(
                std::sin(0.37 * (index + 1)),
                std::cos(0.19 * (index + 1)));
        }

        const ComplexVector coarse =
            preconditioner.apply_coarse(right_hand_side);
        const ComplexVector petrov_residual = model.test_basis().adjoint()
            * (model.operators().system * coarse - right_hand_side);
        require(
            petrov_residual.norm() / right_hand_side.norm() < 1e-10,
            "LOD coarse correction does not satisfy the Petrov equation");

        const ComplexVector right_projected =
            right_hand_side - model.operators().system * coarse;
        const ComplexVector local = preconditioner.apply_local(right_projected);
        HelmholtzTwoLevelSchwarzConfig reuse_config;
        reuse_config.factorization_reuse =
            HelmholtzSchwarzFactorizationReuse::IdenticalMatrix;
        HelmholtzTwoLevelSchwarzPreconditioner reused(
            model, reuse_config);
        const auto &reuse_diagnostics = reused.diagnostics();
        require(
            reuse_diagnostics.local_solver_groups
                < reuse_diagnostics.subdomains,
            "identical-matrix reuse did not reduce factorization count");
        require(reuse_diagnostics.reused_factorizations > 0
                    && reuse_diagnostics.max_reuse_group > 1,
                "identical-matrix reuse diagnostics are inconsistent");
        const ComplexVector reused_local =
            reused.apply_local(right_projected);
        require(
            (reused_local - local).norm()
                    / std::max(1.0, local.norm())
                < 1e-12,
            "grouped block local solve differs from per-patch SparseLU");
        const ComplexVector expected_hybrid = coarse + local
            - preconditioner.apply_coarse(model.operators().system * local);
        const ComplexVector actual_hybrid = preconditioner.apply(
            right_hand_side, HelmholtzTwoLevelSchwarzMode::Hybrid);
        require(
            (actual_hybrid - expected_hybrid).norm()
                    / std::max(1.0, expected_hybrid.norm())
                < 1e-12,
            "two-level hybrid action does not match its factorized formula");

        const ComplexVector reference =
            model.solve_fine_reference(right_hand_side);
        check_mode(
            model, preconditioner, right_hand_side, reference,
            HelmholtzTwoLevelSchwarzMode::Additive, "additive");
        check_mode(
            model, preconditioner, right_hand_side, reference,
            HelmholtzTwoLevelSchwarzMode::Hybrid, "hybrid");
        check_mode(
            model, reused, right_hand_side, reference,
            HelmholtzTwoLevelSchwarzMode::Hybrid, "reused_hybrid");

        HelmholtzTwoLevelSchwarzConfig impedance_config;
        impedance_config.artificial_boundary =
            HelmholtzSchwarzArtificialBoundary::Impedance;
        HelmholtzTwoLevelSchwarzPreconditioner impedance(
            model, impedance_config);
        require(impedance.diagnostics().artificial_boundary_edges > 0,
                "impedance Schwarz did not assemble artificial boundaries");
        require(impedance.diagnostics().physical_boundary_edges > 0,
                "impedance Schwarz did not retain physical boundaries");
        require(impedance.diagnostics().uncovered_dofs == 0,
                "impedance Schwarz patches do not cover the fine space");
        require(impedance.diagnostics().partition_unity_error < 1e-14,
                "impedance Schwarz weights do not form a partition of unity");
        check_mode(
            model, impedance, right_hand_side, reference,
            HelmholtzTwoLevelSchwarzMode::Additive, "impedance_additive");
        check_mode(
            model, impedance, right_hand_side, reference,
            HelmholtzTwoLevelSchwarzMode::Hybrid, "impedance_hybrid");

        HelmholtzTwoLevelSchwarzConfig restricted_config = impedance_config;
        restricted_config.extension =
            HelmholtzSchwarzExtension::RestrictedCore;
        HelmholtzTwoLevelSchwarzPreconditioner restricted(
            model, restricted_config);
        const auto &restricted_diagnostics = restricted.diagnostics();
        require(restricted_diagnostics.uncovered_dofs == 0,
                "restricted Schwarz cores do not cover the fine space");
        require(restricted_diagnostics.partition_unity_error < 1e-14,
                "restricted Schwarz injection is not a partition of unity");
        require(restricted_diagnostics.min_owned_dofs > 0,
                "restricted Schwarz contains an empty core");
        require(
            restricted_diagnostics.max_owned_dofs
                < restricted_diagnostics.max_local_dofs,
            "restricted Schwarz test did not remove overlap on injection");
        check_mode(
            model, restricted, right_hand_side, reference,
            HelmholtzTwoLevelSchwarzMode::Additive, "oras_additive");
        check_mode(
            model, restricted, right_hand_side, reference,
            HelmholtzTwoLevelSchwarzMode::Hybrid, "oras_hybrid");
        HelmholtzTwoLevelSchwarzConfig shifted_config;
        shifted_config.local_solver.kind =
            HelmholtzSchwarzLocalSolverKind::ShiftedGmres;
        shifted_config.local_solver.shift_alpha = 0.2;
        shifted_config.local_solver.gmres.restart = 10;
        shifted_config.local_solver.gmres.max_iterations = 50;
        shifted_config.local_solver.gmres.relative_tolerance = 1e-12;
        HelmholtzTwoLevelSchwarzPreconditioner shifted(
            model, shifted_config);
        check_mode(
            model, shifted, right_hand_side, reference,
            HelmholtzTwoLevelSchwarzMode::Hybrid, "shifted_local_hybrid");
        const HelmholtzSchwarzLocalSolverDiagnostics shifted_diagnostics =
            shifted.local_solver_diagnostics();
        require(shifted_diagnostics.solve_calls > 0,
                "shifted Schwarz did not record local solves");
        require(
            shifted_diagnostics.total_iterations
                > shifted_diagnostics.solve_calls,
            "shifted Schwarz did not exercise iterative local solves");
        require(shifted_diagnostics.max_relative_residual < 1e-11,
                "shifted Schwarz local true residual is too large");

        HelmholtzTwoLevelSchwarzConfig vcycle_config = shifted_config;
        vcycle_config.local_solver.shifted_inverse =
            HelmholtzSchwarzShiftedInverseKind::GeometricVcycle;
        vcycle_config.local_solver.gmres.restart = 30;
        vcycle_config.local_solver.gmres.max_iterations = 200;
        HelmholtzTwoLevelSchwarzPreconditioner vcycle(
            model, vcycle_config);
        require(
            vcycle.local_solver_diagnostics().max_vcycle_levels >= 2,
            "two-level Schwarz did not construct local V-cycle levels");
        require(
            vcycle.local_solver_diagnostics().min_vcycle_coarse_dofs > 0,
            "two-level Schwarz V-cycle has an empty coarse solve");
        check_mode(
            model, vcycle, right_hand_side, reference,
            HelmholtzTwoLevelSchwarzMode::Hybrid, "vcycle_local_hybrid");
        const HelmholtzSchwarzLocalSolverDiagnostics vcycle_diagnostics =
            vcycle.local_solver_diagnostics();
        require(vcycle_diagnostics.solve_calls > 0,
                "V-cycle Schwarz did not record local solves");
        require(vcycle_diagnostics.max_relative_residual < 1e-11,
                "V-cycle Schwarz local true residual is too large");
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_two_level_schwarz failed: "
                  << error.what() << '\n';
        return 1;
    }
}

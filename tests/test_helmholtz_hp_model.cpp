#include "helmholtz/hp_model.h"

#include <algorithm>
#include <iostream>
#include <stdexcept>

using namespace lod2d;
using namespace lod2d::helmholtz;

namespace {

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

template <typename MatrixA, typename MatrixB>
double relative_difference(const MatrixA &a, const MatrixB &b) {
    return (a - b).norm() / std::max(1.0, b.norm());
}

} // namespace

int main() {
    try {
        HelmholtzProblemConfig old_config;
        old_config.H = 3;
        old_config.h = 5;
        old_config.ell = 1;
        old_config.wavenumber = 4.0;
        old_config.mode = HelmholtzPetrovMode::TwoSided;
        old_config.patch_solver.kind =
            HelmholtzPatchSolverKind::DirectSaddle;
        HelmholtzLodModel old_model =
            HelmholtzLodModel::build(old_config);

        HelmholtzHpProblemConfig hp_config;
        hp_config.H = old_config.H;
        hp_config.h = old_config.h;
        hp_config.ell = old_config.ell;
        hp_config.degree = 1;
        hp_config.wavenumber = old_config.wavenumber;
        hp_config.mode = old_config.mode;
        HelmholtzHpLodModel hp_model =
            HelmholtzHpLodModel::build(hp_config);

        require(relative_difference(
                    hp_model.interpolation().coarse_injection,
                    old_model.problem().coarse_to_fine) < 2e-11,
                "hp model p=1 injection differs from the P1 model");
        require(relative_difference(
                    hp_model.interpolation().quasi_interpolation,
                    old_model.problem().quasi_interpolation) < 2e-10,
                "hp model p=1 I_H differs from the P1 model");
        require(relative_difference(
                    hp_model.operators().system,
                    old_model.operators().system) < 2e-11,
                "hp model p=1 fine operator differs from the P1 model");
        require(relative_difference(
                    hp_model.trial_basis(),
                    old_model.trial_basis()) < 2e-9,
                "hp model p=1 trial basis differs from the P1 model");
        require(relative_difference(
                    hp_model.test_basis(),
                    old_model.test_basis()) < 2e-9,
                "hp model p=1 test basis differs from the P1 model");
        require(relative_difference(
                    hp_model.coarse_operator(),
                    old_model.coarse_operator()) < 2e-9,
                "hp model p=1 coarse operator differs from the P1 model");

        HelmholtzHpProblemConfig parallel_config = hp_config;
        parallel_config.corrector_threads = 4;
        HelmholtzHpLodModel parallel_model =
            HelmholtzHpLodModel::build(parallel_config);
        require(relative_difference(
                    parallel_model.corrector_matrix(),
                    hp_model.corrector_matrix()) < 2e-11,
                "parallel hp corrector differs from serial assembly");
        require(relative_difference(
                    parallel_model.adjoint_corrector_matrix(),
                    hp_model.adjoint_corrector_matrix()) < 2e-11,
                "parallel hp adjoint corrector differs from serial assembly");
        require(relative_difference(
                    parallel_model.coarse_operator(),
                    hp_model.coarse_operator()) < 2e-11,
                "parallel hp coarse operator differs from serial assembly");
        require(parallel_model.corrector_diagnostics().patch_count
                    == hp_model.corrector_diagnostics().patch_count,
                "parallel hp assembly lost a patch");

        HelmholtzHpProblemConfig schur_config = hp_config;
        schur_config.patch_solver.kind =
            HelmholtzPatchSolverKind::DirectSchur;
        schur_config.patch_solver.fallback_to_direct = false;
        HelmholtzHpLodModel schur_model =
            HelmholtzHpLodModel::build(schur_config);
        require(relative_difference(
                    schur_model.corrector_matrix(),
                    hp_model.corrector_matrix()) < 2e-9,
                "DirectSchur hp corrector differs from DirectSaddle");
        require(relative_difference(
                    schur_model.coarse_operator(),
                    hp_model.coarse_operator()) < 2e-9,
                "DirectSchur hp coarse operator differs from DirectSaddle");
        require(schur_model.corrector_diagnostics().max_schur_residual < 1e-10,
                "DirectSchur hp Schur residual is too large");
        require(schur_model.corrector_diagnostics().min_schur_rcond > 1e-14,
                "DirectSchur hp Schur complement is numerically singular");
        require(schur_model.corrector_diagnostics().direct_fallback_count == 0,
                "DirectSchur hp model unexpectedly used direct fallback");

        ComplexVector load(old_model.problem().fine.nodes.size());
        for (int i = 0; i < load.size(); ++i)
            load(i) = Complex(
                std::sin(0.13 * (i + 1)),
                std::cos(0.19 * (i + 1)));
        const HelmholtzLodSolution old_solution =
            old_model.solve_load(load);
        const HelmholtzLodSolution hp_solution =
            hp_model.solve_load(load);
        require(relative_difference(
                    hp_solution.coarse_coefficients,
                    old_solution.coarse_coefficients) < 3e-9,
                "hp model p=1 coarse solution differs from the P1 model");
        require(relative_difference(
                    hp_solution.fine_values,
                    old_solution.fine_values) < 3e-9,
                "hp model p=1 fine solution differs from the P1 model");
        require(hp_solution.petrov_residual < 1e-10,
                "hp model Petrov residual is too large");
        const HelmholtzLodSolution parallel_solution =
            parallel_model.solve_load(load);
        require(relative_difference(
                    parallel_solution.coarse_coefficients,
                    hp_solution.coarse_coefficients) < 2e-11,
                "parallel hp coarse solution differs from serial solve");
        require(relative_difference(
                    parallel_solution.fine_values,
                    hp_solution.fine_values) < 2e-11,
                "parallel hp fine solution differs from serial solve");
        const HelmholtzLodSolution schur_solution =
            schur_model.solve_load(load);
        require(relative_difference(
                    schur_solution.fine_values,
                    hp_solution.fine_values) < 3e-9,
                "DirectSchur hp solution differs from DirectSaddle");
        require(schur_solution.petrov_residual < 1e-10,
                "DirectSchur hp Petrov residual is too large");

        const ComplexSparseMatrix corrector_before =
            hp_model.corrector_matrix();
        ComplexVector second_load(load.size());
        for (int i = 0; i < second_load.size(); ++i)
            second_load(i) = Complex(
                std::cos(0.07 * (i + 1)),
                std::sin(0.11 * (i + 1)));
        const HelmholtzLodSolution second_solution =
            hp_model.solve_load(second_load);
        require(second_solution.petrov_residual < 1e-10,
                "hp model repeated-RHS Petrov residual is too large");
        require(relative_difference(
                    hp_model.corrector_matrix(),
                    corrector_before) == 0.0,
                "hp model changed correctors during repeated-RHS solve");

        std::cout << "Helmholtz hp model p=1 equivalence passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_hp_model failed: "
                  << error.what() << '\n';
        return 1;
    }
}

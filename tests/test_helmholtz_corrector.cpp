#include "helmholtz/model.h"
#include "helmholtz/corrector.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace lod2d;
using namespace lod2d::helmholtz;

namespace {

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

bool on_boundary(const Point2 &point) {
    constexpr double tolerance = 1e-12;
    return std::abs(point.x()) < tolerance || std::abs(point.x() - 1.0) < tolerance
        || std::abs(point.y()) < tolerance || std::abs(point.y() - 1.0) < tolerance;
}

} // namespace

int main() {
    try {
        TriMesh assembly_mesh;
        assembly_mesh.nodes = {
            Point2(0.0, 0.0), Point2(1.0, 0.0),
            Point2(0.0, 1.0), Point2(1.0, 1.0)};
        assembly_mesh.elems = {{0, 1, 2}, {3, 2, 1}};
        std::vector<HelmholtzElementCorrector> assembly_correctors(2);
        assembly_correctors[0] = {
            {0, 0, Complex(1.0, 2.0)},
            {0, 0, Complex(3.0, -1.0)},
            {1, 1, Complex(2.0, 0.0)}};
        assembly_correctors[1] = {
            {0, 2, Complex(-1.0, 0.5)},
            {0, 2, Complex(0.25, -0.5)},
            {2, 1, Complex(4.0, 0.0)}};
        const ComplexSparseMatrix assembled = build_helmholtz_corrector_matrix(
            assembly_mesh, 3, assembly_correctors);
        require(std::abs(assembled.coeff(0, 0) - Complex(4.0, 1.0)) < 1e-14,
                "column-wise corrector assembly did not sum duplicate entries");
        require(std::abs(assembled.coeff(0, 1) - Complex(-0.75, 0.0)) < 1e-14,
                "column-wise corrector assembly mapped a local vertex incorrectly");
        require(std::abs(assembled.coeff(1, 1) - Complex(2.0, 0.0)) < 1e-14,
                "column-wise corrector assembly lost an element contribution");
        require(std::abs(assembled.coeff(2, 2) - Complex(4.0, 0.0)) < 1e-14,
                "column-wise corrector assembly mapped a shared vertex incorrectly");
        require(assembled.nonZeros() == 4,
                "column-wise corrector assembly retained unexpected entries");

        HelmholtzProblemConfig config;
        config.H = 1;
        config.h = 4;
        config.ell = 1;
        config.wavenumber = 1.5;
        config.patch_solver.maximum_parallel_solves = 1;
        HelmholtzLodModel model = HelmholtzLodModel::build(config);

        const auto &diagnostics = model.correctors().diagnostics;
        require(diagnostics.patch_count == static_cast<int>(model.problem().coarse.elems.size()),
                "not every coarse element produced a Helmholtz corrector");
        require(diagnostics.patches_touching_physical_boundary > 0,
                "no corrector patch retained the physical Robin boundary");
        require(diagnostics.max_primal_residual < 1e-10,
                "primal Helmholtz corrector residual is too large");
        require(diagnostics.max_adjoint_residual < 1e-10,
                "adjoint Helmholtz corrector residual is too large");
        require(diagnostics.max_constraint_residual < 1e-10,
                "Helmholtz corrector violates ker(I_H)");
        require(diagnostics.parallel_threads == 1,
                "Helmholtz corrector ignored the patch concurrency cap");
        require(diagnostics.symbolic_analyses + diagnostics.symbolic_reuses
                    == diagnostics.patch_count,
                "Helmholtz symbolic factorization statistics do not cover all patches");

        std::vector<int> global_incidence(model.problem().fine.nodes.size(), 0);
        for (const Triangle &triangle : model.problem().fine.elems)
            for (int vertex : triangle) ++global_incidence[vertex];
        for (int target = 0;
             target < static_cast<int>(model.problem().coarse.elems.size());
             ++target) {
            std::vector<int> patch_incidence(model.problem().fine.nodes.size(), 0);
            for (Eigen::SparseMatrix<double>::InnerIterator patch_it(
                     model.problem().patches, target);
                 patch_it;
                 ++patch_it) {
                if (patch_it.value() == 0.0) continue;
                const int coarse_element = patch_it.row();
                for (Eigen::SparseMatrix<double>::InnerIterator child_it(
                         model.problem().fine_element_prolongation, coarse_element);
                     child_it;
                     ++child_it) {
                    if (child_it.value() == 0.0) continue;
                    for (int vertex : model.problem().fine.elems[child_it.row()])
                        ++patch_incidence[vertex];
                }
            }
            for (const auto &entry : model.correctors().primal[target]) {
                require(patch_incidence[entry.row] == global_incidence[entry.row],
                        "artificial patch-boundary vertex entered a corrector");
            }
        }

        const ComplexSparseMatrix prolongation = model.problem().coarse_to_fine.cast<Complex>();
        const ComplexSparseMatrix primal_corrector = prolongation - model.corrected_trial_basis();
        const ComplexMatrix constraint_residual = model.problem().quasi_interpolation.cast<Complex>()
                                                * primal_corrector;
        require(constraint_residual.norm() < 1e-10,
                "assembled corrector matrix is not in ker(I_H)");

        const ComplexMatrix primal_dense(model.corrected_trial_basis());
        const ComplexMatrix adjoint_dense(model.corrected_test_basis());
        require((adjoint_dense - primal_dense.conjugate()).norm() < 1e-11,
                "adjoint corrected basis is not the conjugate primal basis");

        HelmholtzProblemConfig reuse_config = config;
        reuse_config.patch_solver.symbolic_cache_slots = 8;
        reuse_config.patch_solver.reuse_identical_factorization = true;
#ifdef _OPENMP
        const int original_thread_count = omp_get_max_threads();
        // The factorization cache is thread-local.  A parallel dynamic
        // schedule is not required to send two identical patch patterns to
        // the same worker, so use one worker for this cache-behavior test.
        omp_set_num_threads(1);
#endif
        HelmholtzLodModel reuse_model = HelmholtzLodModel::build(reuse_config);
#ifdef _OPENMP
        omp_set_num_threads(original_thread_count);
#endif
        const ComplexMatrix reuse_basis(reuse_model.corrected_trial_basis());
        require((reuse_basis - primal_dense).norm() < 1e-11,
                "identical patch factorization reuse changed the corrector basis");
        require(reuse_model.correctors().diagnostics.factorization_reuses > 0,
                "identical patch factorization reuse was not exercised");

        double boundary_corrector_max = 0.0;
        for (int col = 0; col < primal_corrector.outerSize(); ++col) {
            for (ComplexSparseMatrix::InnerIterator it(primal_corrector, col); it; ++it) {
                if (!on_boundary(model.problem().fine.nodes[it.row()])) continue;
                boundary_corrector_max = std::max(boundary_corrector_max, std::abs(it.value()));
            }
        }
        require(boundary_corrector_max > 1e-12,
                "physical Robin boundary DOFs were incorrectly eliminated from all correctors");

        HelmholtzProblemConfig schur_config = config;
        schur_config.patch_solver.kind = HelmholtzPatchSolverKind::DirectSchur;
        HelmholtzLodModel schur_model = HelmholtzLodModel::build(schur_config);
        const ComplexMatrix schur_basis(schur_model.corrected_trial_basis());
        const double direct_schur_difference =
            (schur_basis - primal_dense).norm() / std::max(1.0, primal_dense.norm());
        require(direct_schur_difference < 1e-10,
                "DirectSchur correctors disagree with DirectSaddle");
        require(schur_model.correctors().diagnostics.max_schur_residual < 1e-11,
                "DirectSchur multiplier residual is too large");
        require(schur_model.correctors().diagnostics.min_schur_reciprocal_condition > 0.0,
                "DirectSchur reported a singular Schur complement");

        HelmholtzProblemConfig gmres_config = config;
        gmres_config.patch_solver.kind = HelmholtzPatchSolverKind::ShiftedGmres;
        gmres_config.patch_solver.shifted.rule = HelmholtzShiftRule::KappaSquared;
        gmres_config.patch_solver.shifted.alpha = 0.2;
        gmres_config.patch_solver.gmres.restart = 30;
        gmres_config.patch_solver.gmres.max_iterations = 120;
        gmres_config.patch_solver.gmres.relative_tolerance = 1e-11;
        HelmholtzLodModel gmres_model = HelmholtzLodModel::build(gmres_config);
        const auto &gmres_diagnostics = gmres_model.correctors().diagnostics;
        const ComplexMatrix gmres_basis(gmres_model.corrected_trial_basis());
        const double gmres_difference =
            (gmres_basis - primal_dense).norm() / std::max(1.0, primal_dense.norm());
        require(gmres_difference < 1e-8,
                "shifted-Laplacian GMRES correctors disagree with DirectSaddle");
        require(gmres_diagnostics.gmres_right_hand_sides > 0,
                "shifted-Laplacian path did not execute GMRES");
        require(gmres_diagnostics.gmres_iterations
                    >= gmres_diagnostics.gmres_right_hand_sides,
                "GMRES iteration statistics are inconsistent");
        require(gmres_diagnostics.max_gmres_relative_residual < 1e-10,
                "shifted-Laplacian GMRES true residual is too large");
        require(gmres_diagnostics.max_primal_residual < 1e-9,
                "shifted-Laplacian corrector primal residual is too large");
        require(gmres_diagnostics.max_constraint_residual < 1e-9,
                "shifted-Laplacian corrector constraint residual is too large");
        require(gmres_diagnostics.direct_fallbacks == 0,
                "shifted-Laplacian correctness test unexpectedly used a fallback");

        HelmholtzProblemConfig exact_shift_config = gmres_config;
        exact_shift_config.patch_solver.shifted.rule = HelmholtzShiftRule::Absolute;
        exact_shift_config.patch_solver.shifted.absolute_epsilon = 0.0;
        HelmholtzLodModel exact_shift_model =
            HelmholtzLodModel::build(exact_shift_config);
        const auto &exact_shift_diagnostics =
            exact_shift_model.correctors().diagnostics;
        require(exact_shift_diagnostics.gmres_max_iterations == 1,
                "epsilon=0 shifted GMRES should converge in one iteration");
        const ComplexMatrix exact_shift_basis(
            exact_shift_model.corrected_trial_basis());
        require(
            (exact_shift_basis - primal_dense).norm()
                / std::max(1.0, primal_dense.norm()) < 1e-10,
            "epsilon=0 shifted GMRES disagrees with DirectSaddle");

        HelmholtzProblemConfig vcycle_config = gmres_config;
        vcycle_config.patch_solver.shifted.inverse =
            HelmholtzShiftedInverseKind::GeometricVcycle;
        vcycle_config.patch_solver.shifted.pre_smooth = 2;
        vcycle_config.patch_solver.shifted.post_smooth = 2;
        vcycle_config.patch_solver.shifted.coarse_max_dofs = 20;
        vcycle_config.patch_solver.shifted.jacobi_weight = 0.6;
        HelmholtzLodModel vcycle_model =
            HelmholtzLodModel::build(vcycle_config);
        const auto &vcycle_diagnostics =
            vcycle_model.correctors().diagnostics;
        const ComplexMatrix vcycle_basis(
            vcycle_model.corrected_trial_basis());
        require(
            (vcycle_basis - primal_dense).norm()
                / std::max(1.0, primal_dense.norm()) < 1e-8,
            "V-cycle shifted GMRES disagrees with DirectSaddle");
        require(vcycle_diagnostics.max_vcycle_levels >= 2,
                "V-cycle regression test did not use a multilevel hierarchy");
        require(vcycle_diagnostics.max_vcycle_coarse_dofs
                    < vcycle_diagnostics.max_vcycle_finest_dofs,
                "V-cycle coarse grid is not smaller than the fine patch grid");
        require(vcycle_diagnostics.max_vcycle_relative_residual < 0.5,
                "one V-cycle does not sufficiently reduce the shifted residual");
        require(vcycle_diagnostics.max_gmres_relative_residual < 1e-10,
                "V-cycle GMRES true residual is too large");
        require(vcycle_diagnostics.direct_fallbacks == 0,
                "V-cycle regression test unexpectedly used a fallback");

        HelmholtzProblemConfig identity_config = gmres_config;
        identity_config.patch_solver.shifted.inverse =
            HelmholtzShiftedInverseKind::Identity;
        HelmholtzLodModel identity_model =
            HelmholtzLodModel::build(identity_config);
        const auto &identity_diagnostics =
            identity_model.correctors().diagnostics;
        const ComplexMatrix identity_basis(
            identity_model.corrected_trial_basis());
        require(
            (identity_basis - primal_dense).norm()
                / std::max(1.0, primal_dense.norm()) < 1e-8,
            "unpreconditioned GMRES correctors disagree with DirectSaddle");
        require(identity_diagnostics.max_gmres_relative_residual < 1e-10,
                "unpreconditioned GMRES true residual is too large");
        require(identity_diagnostics.gmres_max_iterations
                    > gmres_diagnostics.gmres_max_iterations,
                "unpreconditioned GMRES unexpectedly used no extra iterations");
        require(identity_diagnostics.direct_fallbacks == 0,
                "unpreconditioned GMRES unexpectedly used a fallback");
        std::cout << "Helmholtz correctors: primal=" << diagnostics.max_primal_residual
                  << " adjoint=" << diagnostics.max_adjoint_residual
                  << " constraint=" << diagnostics.max_constraint_residual << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_corrector failed: " << error.what() << '\n';
        return 1;
    }
}

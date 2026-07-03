#include "helmholtz/model.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

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
        HelmholtzProblemConfig config;
        config.H = 1;
        config.h = 4;
        config.ell = 1;
        config.wavenumber = 1.5;
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
        require(diagnostics.parallel_threads >= 1,
                "Helmholtz corrector reported an invalid thread count");
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

        double boundary_corrector_max = 0.0;
        for (int col = 0; col < primal_corrector.outerSize(); ++col) {
            for (ComplexSparseMatrix::InnerIterator it(primal_corrector, col); it; ++it) {
                if (!on_boundary(model.problem().fine.nodes[it.row()])) continue;
                boundary_corrector_max = std::max(boundary_corrector_max, std::abs(it.value()));
            }
        }
        require(boundary_corrector_max > 1e-12,
                "physical Robin boundary DOFs were incorrectly eliminated from all correctors");

        std::cout << "Helmholtz correctors: primal=" << diagnostics.max_primal_residual
                  << " adjoint=" << diagnostics.max_adjoint_residual
                  << " constraint=" << diagnostics.max_constraint_residual << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_corrector failed: " << error.what() << '\n';
        return 1;
    }
}

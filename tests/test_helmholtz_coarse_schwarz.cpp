#include "helmholtz/coarse_schwarz.h"
#include "helmholtz/model.h"
#include "solver/right_gmres.h"

#include <Eigen/SparseLU>
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

} // namespace

int main() {
    try {
        HelmholtzProblemConfig config;
        config.H = 3;
        config.h = 6;
        config.ell = 1;
        config.wavenumber = 2.0;
        HelmholtzLodModel model = HelmholtzLodModel::build(config);

        const ComplexSparseMatrix &matrix = model.coarse_operator();
        HelmholtzCoarseRasPreconditioner ras(
            matrix, model.problem().coarse, model.problem().patches);
        const auto &diagnostics = ras.diagnostics();
        require(
            diagnostics.subdomains
                == static_cast<int>(model.problem().coarse.elems.size()),
            "coarse RAS did not create one subdomain per element patch");
        require(diagnostics.min_local_dofs > 0,
                "coarse RAS contains an empty subdomain");
        require(diagnostics.max_local_dofs < matrix.rows(),
                "coarse RAS test did not exercise proper local subdomains");
        require(diagnostics.partition_unity_error < 1e-14,
                "coarse RAS weights do not form a partition of unity");

        ComplexVector rhs(matrix.rows());
        for (int index = 0; index < rhs.size(); ++index) {
            rhs(index) = Complex(
                std::sin(0.37 * (index + 1)),
                std::cos(0.19 * (index + 1)));
        }

        Eigen::SparseLU<ComplexSparseMatrix> direct;
        direct.analyzePattern(matrix);
        direct.factorize(matrix);
        require(direct.info() == Eigen::Success,
                "coarse LOD reference factorization failed");
        const ComplexVector reference = direct.solve(rhs);
        require(direct.info() == Eigen::Success && reference.allFinite(),
                "coarse LOD reference solve failed");

        solver::RightGmresConfig gmres_config;
        gmres_config.restart = 50;
        gmres_config.max_iterations = 500;
        gmres_config.relative_tolerance = 1e-11;
        const auto apply_operator = [&](const solver::ComplexVector &vector) {
            return solver::ComplexVector(matrix * vector);
        };
        const auto apply_ras = [&](const solver::ComplexVector &vector) {
            return solver::ComplexVector(ras.apply(vector));
        };
        const solver::RightGmresResult result =
            solver::solve_right_preconditioned_gmres(
                matrix.rows(), apply_operator, apply_ras, rhs, gmres_config);
        require(result.converged, "coarse RAS GMRES did not converge");
        require(result.relative_residual < 1e-10,
                "coarse RAS GMRES true residual is too large");
        require(
            (result.solution - reference).norm()
                    / std::max(1.0, reference.norm())
                < 1e-9,
            "coarse RAS GMRES disagrees with SparseLU");

        std::cout << "coarse RAS: subdomains=" << diagnostics.subdomains
                  << " local_dofs=" << diagnostics.min_local_dofs
                  << ".." << diagnostics.max_local_dofs
                  << " iterations=" << result.iterations
                  << " residual=" << result.relative_residual << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_coarse_schwarz failed: "
                  << error.what() << '\n';
        return 1;
    }
}

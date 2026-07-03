#include "helmholtz/model.h"

#include <Eigen/LU>
#include <iostream>
#include <stdexcept>

using namespace lod2d;
using namespace lod2d::helmholtz;

namespace {

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

void verify_mode(HelmholtzPetrovMode mode) {
    HelmholtzProblemConfig config;
    config.H = 1;
    config.h = 4;
    config.ell = 1;
    config.wavenumber = 2.25;
    config.mode = mode;
    HelmholtzLodModel model = HelmholtzLodModel::build(config);

    auto source = [](const Point2 &point) {
        return Complex(1.0 + point.x(), 0.5 - point.y());
    };
    const ComplexVector load = assemble_helmholtz_load(model.problem().fine, source);
    const HelmholtzLodSolution solution = model.solve_load(load);
    require(solution.petrov_residual < 1e-11,
            "Petrov-Galerkin residual is too large");

    const ComplexVector coarse_rhs = model.test_basis().adjoint() * load;
    const ComplexMatrix dense_operator(model.coarse_operator());
    const ComplexVector dense_coefficients = dense_operator.fullPivLu().solve(coarse_rhs);
    require((dense_coefficients - solution.coarse_coefficients).norm()
                < 1e-10 * std::max(1.0, dense_coefficients.norm()),
            "sparse Petrov-Galerkin solve disagrees with dense reference");

    const ComplexVector reconstructed = model.trial_basis() * solution.coarse_coefficients;
    require((reconstructed - solution.fine_values).norm() < 1e-12,
            "Petrov-Galerkin reconstruction uses the wrong trial basis");

    const ComplexVector fine_reference = model.solve_fine_reference(load);
    const double fine_residual = (model.operators().system * fine_reference - load).norm()
                               / std::max(1.0, load.norm());
    require(fine_residual < 1e-11, "fine Helmholtz reference solve residual is too large");

    const ComplexVector second_load = assemble_helmholtz_load(
        model.problem().fine,
        [](const Point2 &point) { return Complex(point.x() * point.y(), -0.25); });
    const HelmholtzLodSolution second_solution = model.solve_load(second_load);
    require(second_solution.petrov_residual < 1e-11,
            "reused Petrov-Galerkin factorization failed for a second right-hand side");
}

} // namespace

int main() {
    try {
        verify_mode(HelmholtzPetrovMode::TwoSided);
        verify_mode(HelmholtzPetrovMode::CorrectedTestOnly);
        std::cout << "Helmholtz Petrov-Galerkin modes agree with dense references\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_model failed: " << error.what() << '\n';
        return 1;
    }
}

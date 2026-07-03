#include "helmholtz/model.h"

#include <Eigen/Dense>
#include <complex>
#include <iostream>
#include <stdexcept>

using namespace lod2d;
using namespace lod2d::helmholtz;

namespace {

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

} // namespace

int main() {
    try {
        const HelmholtzProblemData problem = build_helmholtz_problem_data(
            make_helmholtz_unit_square_mesh(), 1, 4, 1);
        require(problem.coarse.dirichlet.empty() && problem.fine.dirichlet.empty(),
                "pure Robin NVB hierarchy must keep every boundary node free");

        const Eigen::MatrixXd reproduction = Eigen::MatrixXd(
            problem.quasi_interpolation * problem.coarse_to_fine);
        const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(
            reproduction.rows(), reproduction.cols());
        require((reproduction - identity).norm() < 1e-11,
                "I_H does not reproduce the coarse P1 space");

        ComplexVector x(problem.fine.nodes.size());
        ComplexVector y(problem.fine.nodes.size());
        for (int i = 0; i < x.size(); ++i) {
            x(i) = Complex(0.25 * i, -0.1 * (i + 1));
            y(i) = Complex(std::sin(i + 1.0), std::cos(0.5 * i));
        }
        const Complex alpha(0.7, -1.3);
        const Complex beta(-0.2, 0.4);
        const ComplexSparseMatrix interpolation = problem.quasi_interpolation.cast<Complex>();
        const ComplexVector combined = interpolation * (alpha * x + beta * y);
        const ComplexVector separate = alpha * (interpolation * x) + beta * (interpolation * y);
        require((combined - separate).norm() < 1e-11,
                "complex quasi-interpolation is not complex linear");

        const ComplexVector coarse = ComplexVector::LinSpaced(
            problem.coarse.nodes.size(), Complex(0.0, 0.0), Complex(1.0, 1.0));
        const ComplexVector projected = interpolation
            * (problem.coarse_to_fine.cast<Complex>() * coarse);
        require((projected - coarse).norm() < 1e-11,
                "complex coarse function changed under I_H");

        std::cout << "Helmholtz interpolation: coarse=" << problem.coarse.nodes.size()
                  << " fine=" << problem.fine.nodes.size() << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_interp failed: " << error.what() << '\n';
        return 1;
    }
}

#include "helmholtz/hp_operators.h"
#include "helmholtz/model.h"
#include "mesh/refine.h"

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

double phi(double value) {
    return 16.0 * value * value
        * (1.0 - value) * (1.0 - value);
}

double phi_first(double value) {
    return 32.0 * value - 96.0 * value * value
        + 64.0 * value * value * value;
}

double phi_second(double value) {
    return 32.0 - 192.0 * value + 192.0 * value * value;
}

} // namespace

int main() {
    try {
        constexpr double kappa = 2.0;
        const ComplexFunction exact = [](const Point2 &point) {
            return phi(point.x()) * phi(point.y())
                * std::exp(Complex(0.0, kappa * point.x()));
        };
        const ComplexGradientFunction exact_gradient = [](const Point2 &point) {
            const Complex phase =
                std::exp(Complex(0.0, kappa * point.x()));
            Eigen::Vector2cd gradient;
            gradient(0) = (phi_first(point.x())
                + Complex(0.0, kappa) * phi(point.x()))
                * phi(point.y()) * phase;
            gradient(1) = phi(point.x())
                * phi_first(point.y()) * phase;
            return gradient;
        };
        const ComplexFunction source = [](const Point2 &point) {
            const Complex phase =
                std::exp(Complex(0.0, kappa * point.x()));
            return -phase * (
                phi_second(point.x()) * phi(point.y())
                + phi(point.x()) * phi_second(point.y())
                + Complex(0.0, 2.0 * kappa)
                    * phi_first(point.x()) * phi(point.y()));
        };

        for (int degree = 1; degree <= 3; ++degree) {
            std::vector<double> diameters;
            std::vector<double> errors;
            for (int level : {2, 4, 6}) {
                TriMesh mesh = refine_mesh_nvb(
                    make_helmholtz_unit_square_mesh(), level).mesh;
                HpTriSpace space(mesh, degree);
                const auto operators =
                    assemble_helmholtz_hp_operators(space, kappa);
                const ComplexVector load =
                    assemble_helmholtz_hp_load(space, source);
                const ComplexVector solution =
                    solve_helmholtz_hp_fem(operators, load);
                const HelmholtzError error =
                    compute_helmholtz_hp_error(
                        space, solution, kappa, exact, exact_gradient);
                diameters.push_back(max_element_diameter(mesh));
                errors.push_back(error.energy);
            }
            require(errors[2] < errors[1] && errors[1] < errors[0],
                    "fine-hp energy error does not decrease");
            const double rate = std::log(errors[1] / errors[2])
                / std::log(diameters[1] / diameters[2]);
            require(rate > degree - 0.45,
                    "fine-hp energy convergence rate is too low");
            std::cout << "p=" << degree
                      << " energy_rate=" << rate
                      << " final_error=" << errors.back() << '\n';
        }
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_hp_fem failed: "
                  << error.what() << '\n';
        return 1;
    }
}

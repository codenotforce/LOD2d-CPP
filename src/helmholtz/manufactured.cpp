#include "helmholtz/manufactured.h"

#include <cmath>
#include <stdexcept>

namespace lod2d::helmholtz {

HelmholtzManufacturedSolution make_polynomial_plane_wave_solution(
    double wavenumber) {
    if (!(wavenumber > 0.0))
        throw std::invalid_argument("manufactured Helmholtz wavenumber must be positive");

    auto phi = [](double t) {
        return 16.0 * t * t * (1.0 - t) * (1.0 - t);
    };
    auto phi_prime = [](double t) {
        return 32.0 * t - 96.0 * t * t + 64.0 * t * t * t;
    };
    auto phi_second = [](double t) {
        return 32.0 - 192.0 * t + 192.0 * t * t;
    };

    HelmholtzManufacturedSolution result;
    result.value = [=](const Point2 &point) {
        return phi(point.x()) * phi(point.y())
            * std::exp(Complex(0.0, wavenumber * point.x()));
    };
    result.gradient = [=](const Point2 &point) {
        const Complex phase =
            std::exp(Complex(0.0, wavenumber * point.x()));
        Eigen::Vector2cd gradient;
        gradient << (phi_prime(point.x())
                        + Complex(0.0, wavenumber) * phi(point.x()))
                            * phi(point.y()) * phase,
                    phi(point.x()) * phi_prime(point.y()) * phase;
        return gradient;
    };
    result.source = [=](const Point2 &point) {
        const Complex phase =
            std::exp(Complex(0.0, wavenumber * point.x()));
        return -phase * (
            phi_second(point.x()) * phi(point.y())
            + phi(point.x()) * phi_second(point.y())
            + Complex(0.0, 2.0 * wavenumber)
                * phi_prime(point.x()) * phi(point.y()));
    };
    return result;
}

} // namespace lod2d::helmholtz

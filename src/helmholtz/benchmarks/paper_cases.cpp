#include "helmholtz/benchmarks/paper_cases.h"

#include "helmholtz/model.h"

#include <cmath>
#include <stdexcept>

namespace lod2d::helmholtz::benchmarks {
namespace {

constexpr double kR0 = 0.25;
constexpr double kR1 = 0.5;
constexpr double kAlpha = 2.0 / 3.0;

struct SmoothStepValue {
    double value = 0.0;
    double first = 0.0;
    double second = 0.0;
};

SmoothStepValue psi(double argument) {
    if (!(argument > 0.0)) return {};
    const double value = std::exp(-1.0 / argument);
    const double inverse = 1.0 / argument;
    return {
        value,
        value * inverse * inverse,
        value * (std::pow(inverse, 4) - 2.0 * std::pow(inverse, 3))};
}

SmoothStepValue cutoff_jet(double radius) {
    if (radius <= kR0) return {1.0, 0.0, 0.0};
    if (radius >= kR1) return {};
    const SmoothStepValue left_psi = psi(kR1 - radius);
    const SmoothStepValue right_psi = psi(radius - kR0);
    const double a = left_psi.value;
    const double a_first = -left_psi.first;
    const double a_second = left_psi.second;
    const double sum = a + right_psi.value;
    const double sum_first = a_first + right_psi.first;
    const double sum_second = a_second + right_psi.second;
    const double value = a / sum;
    const double first = a_first / sum - a * sum_first / (sum * sum);
    const double second = a_second / sum
        - a * sum_second / (sum * sum)
        - 2.0 * a_first * sum_first / (sum * sum)
        + 2.0 * a * sum_first * sum_first / (sum * sum * sum);
    return {value, first, second};
}

struct SingularAmplitude {
    double value = 0.0;
    Eigen::Vector2d gradient = Eigen::Vector2d::Zero();
    double laplacian = 0.0;
};

SingularAmplitude singular_amplitude(const Point2 &point) {
    const double radius = point.norm();
    if (radius == 0.0) return {};
    double angle = std::atan2(point.y(), point.x());
    if (angle < 0.0) angle += 2.0 * std::acos(-1.0);
    if (angle > 1.5 * std::acos(-1.0) + 1e-12)
        throw std::invalid_argument("singular manufactured solution evaluated outside the L-shaped domain");

    const SmoothStepValue cutoff = cutoff_jet(radius);
    if (cutoff.value == 0.0 && cutoff.first == 0.0) return {};
    const double sine = std::sin(kAlpha * angle);
    const double cosine = std::cos(kAlpha * angle);
    const double radial_power = std::pow(radius, kAlpha);
    const double base = radial_power * sine;
    const double base_radial = kAlpha * std::pow(radius, kAlpha - 1.0) * sine;
    const double base_angular_over_radius =
        kAlpha * std::pow(radius, kAlpha - 1.0) * cosine;
    const Eigen::Vector2d radial(point.x() / radius, point.y() / radius);
    const Eigen::Vector2d angular(-radial.y(), radial.x());

    SingularAmplitude result;
    result.value = cutoff.value * base;
    result.gradient =
        (cutoff.first * base + cutoff.value * base_radial) * radial
        + cutoff.value * base_angular_over_radius * angular;
    result.laplacian = 2.0 * cutoff.first * base_radial
        + (cutoff.second + cutoff.first / radius) * base;
    return result;
}

PaperCaseData make_r1(double wavenumber) {
    const HelmholtzManufacturedSolution manufactured =
        make_polynomial_plane_wave_solution(wavenumber);
    PaperCaseData result;
    result.id = experiments::PaperCase::R1;
    result.wavenumber = wavenumber;
    result.initial_mesh = make_helmholtz_unit_square_mesh();
    result.source = manufactured.source;
    result.exact = manufactured.value;
    result.exact_gradient = manufactured.gradient;
    result.exact_laplacian = [=](const Point2 &point) {
        return -manufactured.source(point)
            - wavenumber * wavenumber * manufactured.value(point);
    };
    return result;
}

PaperCaseData make_r2(
    experiments::PaperCase id,
    double wavenumber,
    double sigma) {
    const Point2 center(0.35, 0.55);
    const double normalization = normalized_gaussian_constant(sigma, center);
    PaperCaseData result;
    result.id = id;
    result.wavenumber = wavenumber;
    result.initial_mesh = make_helmholtz_unit_square_mesh();
    result.gaussian_sigma = sigma;
    result.gaussian_normalization = normalization;
    result.quadrature_context.integrand_class = QuadratureClass::LocalizedGaussian;
    result.quadrature_context.feature_point = center;
    result.quadrature_context.feature_scale = sigma;
    result.source = [=](const Point2 &point) {
        return Complex(normalization * std::exp(
            -(point - center).squaredNorm() / (2.0 * sigma * sigma)), 0.0);
    };
    return result;
}

PaperCaseData make_s(double wavenumber) {
    PaperCaseData result;
    result.id = experiments::PaperCase::S;
    result.wavenumber = wavenumber;
    result.initial_mesh = make_helmholtz_l_shape_mesh();
    result.quadrature_context.integrand_class = QuadratureClass::ReentrantSingular;
    result.quadrature_context.feature_point = Point2::Zero();
    result.quadrature_context.feature_scale = kR0;
    result.exact = [=](const Point2 &point) {
        const SingularAmplitude amplitude = singular_amplitude(point);
        return amplitude.value
            * std::exp(Complex(0.0, wavenumber * point.x()));
    };
    result.exact_gradient = [=](const Point2 &point) {
        const SingularAmplitude amplitude = singular_amplitude(point);
        const Complex phase = std::exp(Complex(0.0, wavenumber * point.x()));
        Eigen::Vector2cd gradient = amplitude.gradient.cast<Complex>();
        gradient.x() += Complex(0.0, wavenumber * amplitude.value);
        return (phase * gradient).eval();
    };
    result.exact_laplacian = [=](const Point2 &point) {
        const SingularAmplitude amplitude = singular_amplitude(point);
        const Complex phase = std::exp(Complex(0.0, wavenumber * point.x()));
        return phase * (amplitude.laplacian
            + Complex(0.0, 2.0 * wavenumber * amplitude.gradient.x())
            - wavenumber * wavenumber * amplitude.value);
    };
    result.source = [=](const Point2 &point) {
        const SingularAmplitude amplitude = singular_amplitude(point);
        const Complex phase = std::exp(Complex(0.0, wavenumber * point.x()));
        // -Delta(a exp(ikx)) - k^2 a exp(ikx); the k^2 terms cancel.
        return -phase * (amplitude.laplacian
            + Complex(0.0, 2.0 * wavenumber * amplitude.gradient.x()));
    };
    return result;
}

} // namespace

double normalized_gaussian_constant(double sigma, const Point2 &center) {
    if (!(sigma > 0.0)) throw std::invalid_argument("Gaussian sigma must be positive");
    if (center.x() <= 0.0 || center.x() >= 1.0
        || center.y() <= 0.0 || center.y() >= 1.0) {
        throw std::invalid_argument("Gaussian center must lie inside the unit square");
    }
    const double sqrt_pi = std::sqrt(std::acos(-1.0));
    auto squared_integral = [&](double coordinate) {
        return 0.5 * sigma * sqrt_pi
            * (std::erf((1.0 - coordinate) / sigma)
               + std::erf(coordinate / sigma));
    };
    return 1.0 / std::sqrt(squared_integral(center.x()) * squared_integral(center.y()));
}

double singular_cutoff(double radius) { return cutoff_jet(radius).value; }
double singular_cutoff_prime(double radius) { return cutoff_jet(radius).first; }
double singular_cutoff_second(double radius) { return cutoff_jet(radius).second; }

PaperCaseData make_paper_case(experiments::PaperCase id, double wavenumber) {
    if (!(wavenumber > 0.0))
        throw std::invalid_argument("paper case wavenumber must be positive");
    switch (id) {
    case experiments::PaperCase::R1: return make_r1(wavenumber);
    case experiments::PaperCase::R2a: return make_r2(id, wavenumber, 1.0 / 32.0);
    case experiments::PaperCase::R2b: return make_r2(id, wavenumber, 1.0 / 64.0);
    case experiments::PaperCase::S: return make_s(wavenumber);
    }
    throw std::invalid_argument("unknown paper case");
}

} // namespace lod2d::helmholtz::benchmarks

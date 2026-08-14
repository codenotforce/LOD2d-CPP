#include "helmholtz/benchmarks/paper_cases.h"

#include "helmholtz/model.h"

#include <array>
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

SmoothStepValue cutoff_jet(
    double radius,
    double outer_radius = kR1,
    bool quintic = false) {
    if (radius <= kR0) return {1.0, 0.0, 0.0};
    if (radius >= outer_radius) return {};
    if (quintic) {
        const double width = outer_radius - kR0;
        const double t = (radius - kR0) / width;
        const double t2 = t * t;
        const double t3 = t2 * t;
        const double t4 = t3 * t;
        const double t5 = t4 * t;
        return {
            1.0 - 10.0 * t3 + 15.0 * t4 - 6.0 * t5,
            (-30.0 * t2 + 60.0 * t3 - 30.0 * t4) / width,
            (-60.0 * t + 180.0 * t2 - 120.0 * t3)
                / (width * width)};
    }
    const SmoothStepValue left_psi = psi(outer_radius - radius);
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

struct SmoothWaveEnvelope {
    double value = 0.0;
    Eigen::Vector2d gradient = Eigen::Vector2d::Zero();
    double laplacian = 0.0;
};

SmoothWaveEnvelope smooth_wave_envelope(const Point2 &point) {
    // g(t)=t^2(1-t^2)^2 vanishes together with g' at t=+-1 and
    // vanishes at t=0.  Thus psi=C g(x)g(y) satisfies the outer homogeneous
    // Robin contract and both reentrant Dirichlet rays.  C normalizes max psi
    // to one at |x|=|y|=1/sqrt(3).
    constexpr double normalization = 729.0 / 16.0;
    const auto jet = [](const double value) {
        const double v2 = value * value;
        const double v3 = v2 * value;
        const double v4 = v2 * v2;
        const double v5 = v4 * value;
        const double v6 = v3 * v3;
        return std::array<double, 3>{
            v2 - 2.0 * v4 + v6,
            2.0 * value - 8.0 * v3 + 6.0 * v5,
            2.0 - 24.0 * v2 + 30.0 * v4};
    };
    const auto x = jet(point.x());
    const auto y = jet(point.y());
    SmoothWaveEnvelope result;
    result.value = normalization * x[0] * y[0];
    result.gradient = normalization
        * Eigen::Vector2d(x[1] * y[0], x[0] * y[1]);
    result.laplacian = normalization
        * (x[2] * y[0] + x[0] * y[2]);
    return result;
}

SingularAmplitude singular_amplitude(
    const Point2 &point,
    double cutoff_outer_radius = kR1,
    bool quintic_cutoff = false) {
    const double radius = point.norm();
    if (radius == 0.0) return {};
    double angle = std::atan2(point.y(), point.x());
    if (angle < 0.0) angle += 2.0 * std::acos(-1.0);
    if (angle > 1.5 * std::acos(-1.0) + 1e-12)
        throw std::invalid_argument("singular manufactured solution evaluated outside the L-shaped domain");

    const SmoothStepValue cutoff = cutoff_jet(
        radius, cutoff_outer_radius, quintic_cutoff);
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

PaperCaseData make_s(
    double wavenumber,
    double oscillatory_fraction,
    double cutoff_outer_radius,
    bool quintic_cutoff,
    double smooth_wave_amplitude) {
    if (!std::isfinite(oscillatory_fraction)
        || oscillatory_fraction < 0.0 || oscillatory_fraction > 1.0) {
        throw std::invalid_argument(
            "S oscillatory fraction must lie in [0,1]");
    }
    if (!std::isfinite(cutoff_outer_radius)
        || !(cutoff_outer_radius > kR0)
        || cutoff_outer_radius > 1.0) {
        throw std::invalid_argument(
            "S cutoff outer radius must lie in (0.25,1]");
    }
    if (!std::isfinite(smooth_wave_amplitude)
        || smooth_wave_amplitude < 0.0 || smooth_wave_amplitude > 1.0) {
        throw std::invalid_argument(
            "S smooth wave amplitude must lie in [0,1]");
    }
    PaperCaseData result;
    result.id = experiments::PaperCase::S;
    result.wavenumber = wavenumber;
    result.initial_mesh = make_helmholtz_l_shape_mesh();
    result.quadrature_context.integrand_class = QuadratureClass::ReentrantSingular;
    result.quadrature_context.feature_point = Point2::Zero();
    result.quadrature_context.feature_scale = kR0;
    result.singular_oscillatory_fraction = oscillatory_fraction;
    result.singular_cutoff_outer_radius = cutoff_outer_radius;
    result.singular_quintic_cutoff = quintic_cutoff;
    result.smooth_wave_amplitude = smooth_wave_amplitude;
    result.exact = [=](const Point2 &point) {
        const SingularAmplitude amplitude = singular_amplitude(
            point, cutoff_outer_radius, quintic_cutoff);
        const Complex phase = std::exp(
            Complex(0.0, wavenumber * point.x()));
        const Complex multiplier = (1.0 - oscillatory_fraction)
            + oscillatory_fraction * phase;
        const SmoothWaveEnvelope wave = smooth_wave_envelope(point);
        return amplitude.value * multiplier
            + smooth_wave_amplitude * wave.value * phase;
    };
    result.exact_gradient = [=](const Point2 &point) {
        const SingularAmplitude amplitude = singular_amplitude(
            point, cutoff_outer_radius, quintic_cutoff);
        const Complex phase = std::exp(Complex(0.0, wavenumber * point.x()));
        const Complex multiplier = (1.0 - oscillatory_fraction)
            + oscillatory_fraction * phase;
        Eigen::Vector2cd gradient =
            multiplier * amplitude.gradient.cast<Complex>();
        gradient.x() += oscillatory_fraction * phase
            * Complex(0.0, wavenumber * amplitude.value);
        const SmoothWaveEnvelope wave = smooth_wave_envelope(point);
        gradient += smooth_wave_amplitude * phase
            * wave.gradient.cast<Complex>();
        gradient.x() += smooth_wave_amplitude * phase
            * Complex(0.0, wavenumber * wave.value);
        return gradient;
    };
    result.exact_laplacian = [=](const Point2 &point) {
        const SingularAmplitude amplitude = singular_amplitude(
            point, cutoff_outer_radius, quintic_cutoff);
        const Complex phase = std::exp(Complex(0.0, wavenumber * point.x()));
        const Complex multiplier = (1.0 - oscillatory_fraction)
            + oscillatory_fraction * phase;
        const SmoothWaveEnvelope wave = smooth_wave_envelope(point);
        return multiplier * amplitude.laplacian
            + oscillatory_fraction * phase
                * (Complex(0.0, 2.0 * wavenumber * amplitude.gradient.x())
                   - wavenumber * wavenumber * amplitude.value)
            + smooth_wave_amplitude * phase
                * (wave.laplacian
                   + Complex(0.0, 2.0 * wavenumber * wave.gradient.x())
                   - wavenumber * wavenumber * wave.value);
    };
    result.source = [=](const Point2 &point) {
        const SingularAmplitude amplitude = singular_amplitude(
            point, cutoff_outer_radius, quintic_cutoff);
        const Complex phase = std::exp(Complex(0.0, wavenumber * point.x()));
        const Complex multiplier = (1.0 - oscillatory_fraction)
            + oscillatory_fraction * phase;
        // The oscillatory component retains the original k^2 cancellation;
        // the nonoscillatory corner component contributes -k^2 a.
        const SmoothWaveEnvelope wave = smooth_wave_envelope(point);
        return -multiplier * amplitude.laplacian
            - oscillatory_fraction * phase
                * Complex(0.0, 2.0 * wavenumber * amplitude.gradient.x())
            - (1.0 - oscillatory_fraction)
                * wavenumber * wavenumber * amplitude.value
            - smooth_wave_amplitude * phase
                * (wave.laplacian
                   + Complex(0.0, 2.0 * wavenumber * wave.gradient.x()));
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

PaperCaseData make_paper_case(
    experiments::PaperCase id,
    double wavenumber,
    double singular_oscillatory_fraction,
    double singular_cutoff_outer_radius,
    bool singular_quintic_cutoff,
    double smooth_wave_amplitude) {
    if (!(wavenumber > 0.0))
        throw std::invalid_argument("paper case wavenumber must be positive");
    if (id != experiments::PaperCase::S
        && (singular_oscillatory_fraction != 1.0
            || singular_cutoff_outer_radius != kR1
            || singular_quintic_cutoff
            || smooth_wave_amplitude != 0.0)) {
        throw std::invalid_argument(
            "S manufactured-solution parameters were supplied to a non-S case");
    }
    switch (id) {
    case experiments::PaperCase::R1: return make_r1(wavenumber);
    case experiments::PaperCase::R2a: return make_r2(id, wavenumber, 1.0 / 32.0);
    case experiments::PaperCase::R2b: return make_r2(id, wavenumber, 1.0 / 64.0);
    case experiments::PaperCase::S:
        return make_s(
            wavenumber, singular_oscillatory_fraction,
            singular_cutoff_outer_radius, singular_quintic_cutoff,
            smooth_wave_amplitude);
    }
    throw std::invalid_argument("unknown paper case");
}

} // namespace lod2d::helmholtz::benchmarks

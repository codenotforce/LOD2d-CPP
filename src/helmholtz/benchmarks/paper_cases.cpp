#include "helmholtz/benchmarks/paper_cases.h"

#include "helmholtz/boundary.h"
#include "helmholtz/model.h"

#include <cmath>
#include <stdexcept>

namespace lod2d::helmholtz::benchmarks {
namespace {

constexpr double kR0 = 0.25;
constexpr double kR1 = 0.5;
constexpr double kAlpha = 2.0 / 3.0;
constexpr double kR1Localization = 80.0;
constexpr double kOscillatorySupportHalfWidth = 0.25;

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

TriMesh make_r1_mixed_boundary_mesh() {
    TriMesh mesh = make_helmholtz_unit_square_mesh();
    const auto [edges, boundary] = compute_edges(mesh);
    mesh.boundary_edges.clear();
    for (std::size_t index = 0; index < edges.size(); ++index) {
        if (!boundary[index]) continue;
        const Edge edge = edges[index];
        const Point2 midpoint = 0.5 * (
            mesh.nodes[edge[0]] + mesh.nodes[edge[1]]);
        BoundaryTag tag = BoundaryTag::Robin;
        if (std::abs(midpoint.y()) < 1e-14
            || std::abs(midpoint.y() - 1.0) < 1e-14) {
            tag = BoundaryTag::Dirichlet;
        } else if (std::abs(midpoint.x()) < 1e-14) {
            tag = BoundaryTag::Neumann;
        }
        mesh.boundary_edges.push_back({edge, tag});
    }
    synchronize_dirichlet_nodes(mesh);
    validate_boundary_tags(mesh);
    return mesh;
}

SmoothWaveEnvelope localized_r1_amplitude(const Point2 &point) {
    constexpr double pi = 3.141592653589793238462643383279502884;
    const double x = point.x();
    const double y = point.y();
    const double one_minus_x = 1.0 - x;
    const double polynomial = x * x * one_minus_x * one_minus_x;
    const double polynomial_first = 2.0 * x - 6.0 * x * x
        + 4.0 * x * x * x;
    const double polynomial_second = 2.0 - 12.0 * x + 12.0 * x * x;
    const double sine = std::sin(pi * y);
    const double sine_first = pi * std::cos(pi * y);
    const double sine_second = -pi * pi * sine;
    const double dx = x - 0.75;
    const double dy = y - 0.5;
    const double gaussian = std::exp(
        -kR1Localization * (dx * dx + dy * dy));
    const double gaussian_x = -2.0 * kR1Localization * dx * gaussian;
    const double gaussian_y = -2.0 * kR1Localization * dy * gaussian;
    const double gaussian_laplacian = (
        4.0 * kR1Localization * kR1Localization * (dx * dx + dy * dy)
        - 4.0 * kR1Localization) * gaussian;

    SmoothWaveEnvelope result;
    result.value = polynomial * sine * gaussian;
    result.gradient.x() = sine * (
        polynomial_first * gaussian + polynomial * gaussian_x);
    result.gradient.y() = polynomial * (
        sine_first * gaussian + sine * gaussian_y);
    result.laplacian = polynomial_second * sine * gaussian
        + polynomial * sine_second * gaussian
        + polynomial * sine * gaussian_laplacian
        + 2.0 * polynomial_first * sine * gaussian_x
        + 2.0 * polynomial * sine_first * gaussian_y;
    return result;
}

SmoothWaveEnvelope smooth_wave_envelope(const Point2 &point) {
    // Tensor-product C-infinity bump with support
    // (-0.75,-0.25)x(0.25,0.75). Its support is separated from the corner
    // and boundary, while its dyadic support edges align with the E2 mesh
    // hierarchy and avoid quadrature cells straddling the flat cut-off.
    const auto bump_jet = [](double coordinate, double center) {
        const double s = (coordinate - center)
            / kOscillatorySupportHalfWidth;
        if (!(std::abs(s) < 1.0)) return SmoothStepValue{};
        const double t = 1.0 - s * s;
        const double log_value = 1.0 - 1.0 / t;
        if (log_value < -700.0) return SmoothStepValue{};
        const double value = std::exp(log_value);
        const double log_first = -2.0 * s / (t * t);
        const double log_second = -2.0 / (t * t)
            - 8.0 * s * s / (t * t * t);
        return SmoothStepValue{
            value,
            value * log_first / kOscillatorySupportHalfWidth,
            value * (log_first * log_first + log_second)
                / (kOscillatorySupportHalfWidth
                   * kOscillatorySupportHalfWidth)};
    };
    const SmoothStepValue x = bump_jet(point.x(), -0.5);
    const SmoothStepValue y = bump_jet(point.y(), 0.5);
    SmoothWaveEnvelope result;
    result.value = x.value * y.value;
    result.gradient = Eigen::Vector2d(
        x.first * y.value, x.value * y.first);
    result.laplacian = x.second * y.value + x.value * y.second;
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
    PaperCaseData result;
    result.id = experiments::PaperCase::R1;
    result.wavenumber = wavenumber;
    result.initial_mesh = make_r1_mixed_boundary_mesh();
    result.quadrature_context.integrand_class = QuadratureClass::LocalizedGaussian;
    result.quadrature_context.feature_point = Point2(0.75, 0.5);
    result.quadrature_context.feature_scale = 1.0 / std::sqrt(160.0);
    result.exact = [=](const Point2 &point) {
        const SmoothWaveEnvelope amplitude = localized_r1_amplitude(point);
        return amplitude.value
            * std::exp(Complex(0.0, wavenumber * point.x()));
    };
    result.exact_gradient = [=](const Point2 &point) {
        const SmoothWaveEnvelope amplitude = localized_r1_amplitude(point);
        const Complex phase = std::exp(
            Complex(0.0, wavenumber * point.x()));
        Eigen::Vector2cd gradient = phase
            * amplitude.gradient.cast<Complex>();
        gradient.x() += phase * Complex(0.0, wavenumber * amplitude.value);
        return gradient;
    };
    result.exact_laplacian = [=](const Point2 &point) {
        const SmoothWaveEnvelope amplitude = localized_r1_amplitude(point);
        const Complex phase = std::exp(
            Complex(0.0, wavenumber * point.x()));
        return phase * (
            amplitude.laplacian
            + Complex(0.0, 2.0 * wavenumber * amplitude.gradient.x())
            - wavenumber * wavenumber * amplitude.value);
    };
    result.source = [=](const Point2 &point) {
        const SmoothWaveEnvelope amplitude = localized_r1_amplitude(point);
        const Complex phase = std::exp(
            Complex(0.0, wavenumber * point.x()));
        return -phase * (
            amplitude.laplacian
            + Complex(0.0, 2.0 * wavenumber * amplitude.gradient.x()));
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
    double wavenumber) {
    if (id == experiments::PaperCase::S) {
        return make_paper_case(
            id, wavenumber, 0.0, kR1, false, 0.05);
    }
    return make_paper_case(
        id, wavenumber, 1.0, kR1, false, 0.0);
}

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

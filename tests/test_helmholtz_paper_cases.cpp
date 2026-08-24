#include "helmholtz/benchmarks/paper_cases.h"
#include "helmholtz/boundary.h"
#include "mesh/refine.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace lod2d;
using namespace lod2d::helmholtz;
using namespace lod2d::helmholtz::benchmarks;
using lod2d::helmholtz::experiments::PaperCase;

namespace {

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

double simpson_integral(const std::function<double(double)> &function, int intervals) {
    if (intervals % 2 != 0) ++intervals;
    const double step = 1.0 / intervals;
    double sum = function(0.0) + function(1.0);
    for (int index = 1; index < intervals; ++index)
        sum += (index % 2 == 0 ? 2.0 : 4.0) * function(index * step);
    return step * sum / 3.0;
}

void verify_gaussian(PaperCase id, double expected_sigma) {
    const PaperCaseData data = make_paper_case(id, 8.0);
    require(data.gaussian_sigma == expected_sigma, "paper Gaussian has the wrong sigma");
    const double normalization = *data.gaussian_normalization;
    const int intervals = expected_sigma < 0.02 ? 32768 : 16384;
    const double integral_x = simpson_integral([&](double x) {
        const double value = std::real(data.source(Point2(x, 0.55)));
        return value * value / (normalization * normalization);
    }, intervals);
    const double integral_y = simpson_integral([&](double y) {
        const double value = std::real(data.source(Point2(0.35, y)));
        return value * value / (normalization * normalization);
    }, intervals);
    const double norm = normalization * std::sqrt(integral_x * integral_y);
    require(std::abs(norm - 1.0) < 1e-12,
            "paper Gaussian is not L2-normalized on the unit square");
}

void verify_impedance_unit_square(PaperCase id) {
    const PaperCaseData data = make_paper_case(id, 8.0);
    require(!data.initial_mesh.boundary_edges.empty(),
            "formal unit-square case has no explicit boundary-edge contract");
    validate_boundary_tags(data.initial_mesh);
    require(data.initial_mesh.boundary_edges.size() == 4,
            "formal unit-square case has an incomplete boundary-edge contract");
    require(boundary_edges_with_tag(
                data.initial_mesh, BoundaryTag::Dirichlet).empty(),
            "formal unit-square case unexpectedly has a Dirichlet edge");
    require(boundary_edges_with_tag(
                data.initial_mesh, BoundaryTag::Robin).size() == 4,
            "formal unit-square case is not pure impedance");
    require(std::abs(boundary_measure(
                data.initial_mesh, BoundaryTag::Robin) - 4.0) < 1e-14,
            "formal unit-square impedance measure is wrong");

    // Explicit edge tags are authoritative.  Poisoning the compatibility node
    // list must not switch a formal R2 edge to the legacy classification.
    TriMesh poisoned = data.initial_mesh;
    poisoned.dirichlet = {0, 1, 2, 3};
    for (const BoundaryEdge &entry : poisoned.boundary_edges) {
        require(boundary_tag(poisoned, entry.nodes) == BoundaryTag::Robin,
                "formal unit-square boundary depends on the legacy node fallback");
    }
}

Complex finite_difference_laplacian(
    const ComplexFunction &value,
    const Point2 &point,
    double step) {
    const Point2 ex(step, 0.0);
    const Point2 ey(0.0, step);
    return (value(point + ex) - 2.0 * value(point) + value(point - ex)) / (step * step)
        + (value(point + ey) - 2.0 * value(point) + value(point - ey)) / (step * step);
}

Eigen::Vector2cd finite_difference_gradient(
    const ComplexFunction &value,
    const Point2 &point,
    double step) {
    Eigen::Vector2cd result;
    result.x() = (value(point + Point2(step, 0.0))
                - value(point - Point2(step, 0.0))) / (2.0 * step);
    result.y() = (value(point + Point2(0.0, step))
                - value(point - Point2(0.0, step))) / (2.0 * step);
    return result;
}

void verify_localized_smooth_case() {
    constexpr double kappa = 8.0;
    constexpr double pi = 3.141592653589793238462643383279502884;
    const PaperCaseData data = make_paper_case(PaperCase::R1, kappa);
    validate_boundary_tags(data.initial_mesh);
    require(std::abs(boundary_measure(
                data.initial_mesh, BoundaryTag::Dirichlet) - 2.0) < 1e-14,
            "R1 top/bottom Dirichlet partition is wrong");
    require(std::abs(boundary_measure(
                data.initial_mesh, BoundaryTag::Neumann) - 1.0) < 1e-14,
            "R1 left Neumann partition is wrong");
    require(std::abs(boundary_measure(
                data.initial_mesh, BoundaryTag::Robin) - 1.0) < 1e-14,
            "R1 right Robin partition is wrong");
    TriMesh refined = data.initial_mesh;
    for (int level = 0; level < 3; ++level) {
        refined = refine_nvb(refined).mesh;
        validate_boundary_tags(refined);
        require(std::abs(boundary_measure(
                    refined, BoundaryTag::Dirichlet) - 2.0) < 2e-12
                    && std::abs(boundary_measure(
                        refined, BoundaryTag::Neumann) - 1.0) < 2e-12
                    && std::abs(boundary_measure(
                        refined, BoundaryTag::Robin) - 1.0) < 2e-12,
                "NVB failed to preserve the revised R1 D/N/R partition");
    }

    const Point2 sample(0.73, 0.52);
    const double polynomial = sample.x() * sample.x()
        * (1.0 - sample.x()) * (1.0 - sample.x());
    const double amplitude = polynomial * std::sin(pi * sample.y())
        * std::exp(-80.0 * (
            std::pow(sample.x() - 0.75, 2)
            + std::pow(sample.y() - 0.5, 2)));
    const Complex expected = amplitude
        * std::exp(Complex(0.0, kappa * sample.x()));
    require(std::abs(data.exact(sample) - expected) < 2e-16,
            "R1 exact solution does not match the revised localized formula");
    require((finite_difference_gradient(data.exact, sample, 2e-6)
             - data.exact_gradient(sample)).norm()
                / std::max(1.0, data.exact_gradient(sample).norm()) < 1e-7,
            "R1 analytic gradient disagrees with finite differences");
    require(std::abs(finite_difference_laplacian(data.exact, sample, 2e-5)
                     - data.exact_laplacian(sample))
                / std::max(1.0, std::abs(data.exact_laplacian(sample))) < 2e-6,
            "R1 analytic Laplacian disagrees with finite differences");
    require(std::abs(data.source(sample) + data.exact_laplacian(sample)
                     + kappa * kappa * data.exact(sample)) < 1e-12,
            "R1 source violates the manufactured PDE identity");
    require(std::abs(data.exact(Point2(0.75, 0.5)))
                > 1e5 * std::abs(data.exact(Point2(0.2, 0.5))),
            "R1 solution is not localized near (3/4,1/2)");

    for (double x : {0.2, 0.8}) {
        require(std::abs(data.exact(Point2(x, 0.0))) < 1e-15
                    && std::abs(data.exact(Point2(x, 1.0))) < 1e-15,
                "R1 violates a homogeneous Dirichlet edge");
    }
    const Point2 left(0.0, 0.43);
    require(std::abs(data.exact_gradient(left).x()) < 1e-14,
            "R1 violates the homogeneous left Neumann condition");
    const Point2 right(1.0, 0.43);
    require(std::abs(data.exact_gradient(right).x()
                     - Complex(0.0, kappa) * data.exact(right)) < 1e-14,
            "R1 violates the homogeneous right Robin condition");
}

void verify_singular_case() {
    constexpr double kappa = 8.0;
    const PaperCaseData data = make_paper_case(PaperCase::S, kappa);
    const PaperCaseData legacy_singular_only =
        make_paper_case(PaperCase::S, kappa, 0.0, 0.5, false, 0.0);
    const PaperCaseData revised_singular_only =
        make_paper_case(
            PaperCase::S, kappa, 0.0, 0.5, false, 0.0,
            "boundary-weight-gaussian");
    const PaperCaseData corner_dominant =
        make_paper_case(PaperCase::S, kappa, 0.0, 1.0, true);
    const PaperCaseData additive_wave =
        make_paper_case(PaperCase::S, kappa, 0.0, 1.0, true, 0.1);
    require(data.singular_solution_profile == "boundary-weight-gaussian"
                && legacy_singular_only.singular_solution_profile
                    == "radial-cutoff"
                && data.singular_oscillatory_fraction == 0.0
                && corner_dominant.singular_oscillatory_fraction == 0.0
                && data.singular_cutoff_outer_radius == 0.5
                && corner_dominant.singular_cutoff_outer_radius == 1.0
                && !data.singular_quintic_cutoff
                && corner_dominant.singular_quintic_cutoff
                && data.smooth_wave_amplitude == 0.25
                && additive_wave.smooth_wave_amplitude == 0.1,
            "S oscillatory fraction provenance is missing");
    require(std::abs(boundary_measure(data.initial_mesh, BoundaryTag::Dirichlet) - 2.0) < 1e-14,
            "S case lost its reentrant Dirichlet boundary");
    for (const Point2 point : {Point2(0.2, 0.0), Point2(0.0, -0.2)})
        require(std::abs(data.exact(point)) < 1e-14,
                "S exact solution violates the reentrant Dirichlet boundary");
    for (const Point2 point : {Point2(-1.0, 0.4), Point2(0.3, 1.0), Point2(1.0, 0.7)}) {
        require(std::abs(data.exact(point)) < 1e-14,
                "S cut-off reaches the outer Robin boundary");
        require(data.exact_gradient(point).norm() < 1e-14,
                "S gradient reaches the outer Robin boundary");
    }

    for (const PaperCaseData *variant : {
             &data, &corner_dominant, &additive_wave}) {
        for (const Point2 point : {
                 Point2(-0.18, 0.11), Point2(-0.31, 0.19),
                 Point2(0.16, 0.28), Point2(-0.46, 0.52)}) {
            const Complex numerical = finite_difference_laplacian(
                variant->exact, point, 2e-5);
            const Complex analytic = variant->exact_laplacian(point);
            require(std::abs(numerical - analytic)
                        / std::max(1.0, std::abs(analytic)) < 2e-5,
                    "S analytic Laplacian disagrees with an independent finite difference");
            require(std::abs(variant->source(point) + analytic
                             + kappa * kappa * variant->exact(point)) < 1e-11,
                    "S source and manufactured solution violate the PDE identity");
        }
    }
    require(std::abs(std::imag(corner_dominant.exact(Point2(-0.18, 0.11))))
                < 1e-14,
            "corner-dominant S exact solution is unexpectedly oscillatory");
    const Point2 oscillation_center(-0.5, 0.5);
    require(std::abs(data.exact(oscillation_center)
                     - revised_singular_only.exact(oscillation_center)
                     - Complex(0.25, 0.0)) < 1e-13,
            "revised S Gaussian wave lost its normalized center amplitude");
    require(std::abs(data.exact(Point2(-0.5, 0.5))
                     - revised_singular_only.exact(Point2(-0.5, 0.5)))
                > 50.0
                    * std::abs(data.exact(Point2(-0.9, 0.5))
                               - revised_singular_only.exact(Point2(-0.9, 0.5))),
            "revised S oscillatory wave is not localized near (-1/2,1/2)");
    const Point2 uncut_sample(-0.8, 0.2);
    require(std::abs(revised_singular_only.exact(uncut_sample)) > 1e-3
                && std::abs(legacy_singular_only.exact(uncut_sample)) < 1e-14,
            "revised S singularity still contains the legacy radial cut-off");
    require(std::abs(std::imag(additive_wave.exact(Point2(-0.45, 0.42))))
                > 1e-5,
            "additive S exact solution lost its smooth oscillatory component");
    for (const Point2 point : {
             Point2(-1.0, 0.4), Point2(0.3, 1.0), Point2(1.0, 0.7)}) {
        require(std::abs(additive_wave.exact(point)) < 1e-14
                    && additive_wave.exact_gradient(point).norm() < 1e-13,
                "additive smooth wave violates the homogeneous outer boundary");
    }

    bool invalid_fraction_rejected = false;
    try {
        (void)make_paper_case(PaperCase::S, kappa, 1.01);
    } catch (const std::invalid_argument &) {
        invalid_fraction_rejected = true;
    }
    require(invalid_fraction_rejected,
            "S accepted an invalid oscillatory fraction");
    bool invalid_radius_rejected = false;
    try {
        (void)make_paper_case(PaperCase::S, kappa, 0.0, 1.01);
    } catch (const std::invalid_argument &) {
        invalid_radius_rejected = true;
    }
    require(invalid_radius_rejected,
            "S accepted an invalid cut-off radius");
    bool invalid_wave_amplitude_rejected = false;
    try {
        (void)make_paper_case(PaperCase::S, kappa, 0.0, 1.0, true, 1.01);
    } catch (const std::invalid_argument &) {
        invalid_wave_amplitude_rejected = true;
    }
    require(invalid_wave_amplitude_rejected,
            "S accepted an invalid smooth wave amplitude");
    bool invalid_profile_rejected = false;
    try {
        (void)make_paper_case(
            PaperCase::S, kappa, 0.0, 0.5, false, 0.25,
            "unknown-profile");
    } catch (const std::invalid_argument &) {
        invalid_profile_rejected = true;
    }
    require(invalid_profile_rejected,
            "S accepted an unknown manufactured-solution profile");
    bool non_s_parameters_rejected = false;
    try {
        (void)make_paper_case(PaperCase::R1, kappa, 0.0, 1.0, true, 0.1);
    } catch (const std::invalid_argument &) {
        non_s_parameters_rejected = true;
    }
    require(non_s_parameters_rejected,
            "a non-S case accepted S exact-solution parameters");
}

void verify_cutoff_jet() {
    require(singular_cutoff(0.2) == 1.0 && singular_cutoff(0.6) == 0.0,
            "singular cut-off support is wrong");
    for (double radius : {0.3, 0.37, 0.45}) {
        const double step = 1e-6;
        const double first = (singular_cutoff(radius + step)
                            - singular_cutoff(radius - step)) / (2.0 * step);
        const double second = (singular_cutoff(radius + step)
            - 2.0 * singular_cutoff(radius)
            + singular_cutoff(radius - step)) / (step * step);
        require(std::abs(first - singular_cutoff_prime(radius)) < 5e-8,
                "singular cut-off first derivative is wrong");
        require(std::abs(second - singular_cutoff_second(radius)) < 2e-3,
                "singular cut-off second derivative is wrong");
    }
}

} // namespace

int main() {
    try {
        verify_localized_smooth_case();
        verify_impedance_unit_square(PaperCase::R2a);
        verify_impedance_unit_square(PaperCase::R2b);
        verify_gaussian(PaperCase::R2a, 1.0 / 32.0);
        verify_gaussian(PaperCase::R2b, 1.0 / 64.0);
        verify_cutoff_jet();
        verify_singular_case();
        std::cout << "Helmholtz paper case data passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_paper_cases failed: " << error.what() << '\n';
        return 1;
    }
}

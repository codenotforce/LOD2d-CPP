#include "helmholtz/benchmarks/paper_cases.h"
#include "helmholtz/boundary.h"

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
    // list must not switch a formal R1/R2 edge to the legacy classification.
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

void verify_singular_case() {
    constexpr double kappa = 8.0;
    const PaperCaseData data = make_paper_case(PaperCase::S, kappa);
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

    for (const Point2 point : {Point2(-0.18, 0.11), Point2(-0.31, 0.19), Point2(0.16, 0.28)}) {
        const Complex numerical = finite_difference_laplacian(data.exact, point, 2e-5);
        const Complex analytic = data.exact_laplacian(point);
        require(std::abs(numerical - analytic) / std::max(1.0, std::abs(analytic)) < 2e-5,
                "S analytic Laplacian disagrees with an independent finite difference");
        require(std::abs(data.source(point) + analytic
                         + kappa * kappa * data.exact(point)) < 1e-11,
                "S source and manufactured solution violate the PDE identity");
    }
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
        verify_impedance_unit_square(PaperCase::R1);
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

#include "helmholtz/benchmarks/paper_cases.h"
#include "helmholtz/operators.h"
#include "helmholtz/quadrature.h"
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

void verify_geometry_and_polynomials() {
    const TriMesh square = make_paper_case(PaperCase::R2a, 8.0).initial_mesh;
    const double area = integrate_scalar_function(
        square, [](const Point2 &) { return 1.0; });
    require(std::abs(area - 1.0) < 2e-14,
            "unified quadrature does not integrate unit-square area");
    const double polynomial = integrate_scalar_function(square, [](const Point2 &point) {
        return point.x() * point.x() * point.y() + 2.0 * point.x();
    });
    require(std::abs(polynomial - (1.0 / 6.0 + 1.0)) < 3e-14,
            "unified quadrature does not integrate a low polynomial exactly");

    const TriMesh l_shape = make_paper_case(PaperCase::S, 8.0).initial_mesh;
    const double l_area = integrate_scalar_function(
        l_shape, [](const Point2 &) { return 1.0; });
    require(std::abs(l_area - 3.0) < 5e-14,
            "unified quadrature does not integrate L-shaped area");
}

void verify_gaussian(PaperCase id) {
    const PaperCaseData data = make_paper_case(id, 8.0);
    QuadraturePolicy policy;
    const double squared_norm = integrate_scalar_function(
        data.initial_mesh,
        [&](const Point2 &point) { return std::norm(data.source(point)); },
        policy,
        data.quadrature_context);
    require(std::abs(std::sqrt(squared_norm) - 1.0) < 1e-10,
            "case-specific Gaussian quadrature misses L2 normalization");

    QuadraturePolicy stronger = policy;
    stronger.gaussian_triangle_order = 24;
    const double strengthened = integrate_scalar_function(
        data.initial_mesh,
        [&](const Point2 &point) { return std::norm(data.source(point)); },
        stronger,
        data.quadrature_context);
    require(std::abs(strengthened - squared_norm) < 2e-11,
            "Gaussian integral is not stable under quadrature strengthening");

    const ComplexVector load = assemble_helmholtz_load(
        data.initial_mesh, data.source, policy, data.quadrature_context);
    const ComplexVector strengthened_load = assemble_helmholtz_load(
        data.initial_mesh, data.source, stronger, data.quadrature_context);
    require((strengthened_load - load).norm() /
                std::max(1.0, strengthened_load.norm()) < 1e-10,
            "Gaussian load vector is not stable under quadrature strengthening");
}

void verify_singular_convergence() {
    const PaperCaseData data = make_paper_case(PaperCase::S, 8.0);
    const PaperCaseData singular_only =
        make_paper_case(PaperCase::S, 8.0, 0.0, 0.5, false, 0.0);
    // Production E2 assembles on a deeply refined reference/candidate mesh.
    // Resolve the radius-0.2 interior bump geometrically before comparing
    // quadrature orders; the six-element seed mesh is intentionally far
    // outside that production regime.
    const TriMesh mesh = refine_mesh_nvb(data.initial_mesh, 4).mesh;
    QuadraturePolicy policy;
    policy.singular_triangle_order = 32;
    const double baseline = integrate_scalar_function(
        mesh,
        [&](const Point2 &point) { return std::norm(data.source(point)); },
        policy,
        data.quadrature_context);
    QuadraturePolicy stronger = policy;
    stronger.singular_triangle_order = 48;
    const double strengthened = integrate_scalar_function(
        mesh,
        [&](const Point2 &point) { return std::norm(data.source(point)); },
        stronger,
        data.quadrature_context);
    require(baseline > 0.0 && std::isfinite(baseline),
            "singular source integral is not finite and positive");
    const double source_relative_change =
        std::abs(strengthened - baseline) / strengthened;
    std::cout << "S source squared integral on level4 mesh: order32=" << baseline
              << " order48=" << strengthened
              << " relative_change="
              << source_relative_change << '\n';

    const ComplexVector load = assemble_helmholtz_load(
        mesh, data.source, policy, data.quadrature_context);
    const ComplexVector strengthened_load = assemble_helmholtz_load(
        mesh, data.source, stronger, data.quadrature_context);
    const double load_relative_change = (strengthened_load - load).norm()
        / std::max(1.0, strengthened_load.norm());
    const ComplexVector singular_load = assemble_helmholtz_load(
        mesh, singular_only.source, policy, singular_only.quadrature_context);
    const ComplexVector strengthened_singular_load = assemble_helmholtz_load(
        mesh, singular_only.source, stronger, singular_only.quadrature_context);
    const double singular_load_relative_change =
        (strengthened_singular_load - singular_load).norm()
        / std::max(1.0, strengthened_singular_load.norm());
    std::cout << "S load relative_change: full=" << load_relative_change
              << " singular_only=" << singular_load_relative_change << '\n';

    const ComplexVector zero = ComplexVector::Zero(
        static_cast<Eigen::Index>(mesh.nodes.size()));
    const HelmholtzError error = compute_helmholtz_error(
        mesh, zero, data.wavenumber, data.exact,
        data.exact_gradient, policy, data.quadrature_context);
    const HelmholtzError strengthened_error = compute_helmholtz_error(
        mesh, zero, data.wavenumber, data.exact,
        data.exact_gradient, stronger, data.quadrature_context);
    const double energy_relative_change =
        std::abs(strengthened_error.energy - error.energy)
        / strengthened_error.energy;
    const double l2_relative_change =
        std::abs(strengthened_error.l2 - error.l2)
        / strengthened_error.l2;
    std::cout << "S exact norm relative_change: energy="
              << energy_relative_change << " L2=" << l2_relative_change << '\n';
    require(source_relative_change < 2e-6,
            "S source integral is not stable under Duffy order refinement");
    require(load_relative_change < 2e-6,
            "S load vector is not stable under Duffy order refinement");
    require(singular_load_relative_change < 2e-6,
            "S singular-only load is not stable under Duffy order refinement");
    require(energy_relative_change < 2e-6 && l2_relative_change < 2e-6,
            "S exact-error norms are not stable under Duffy order refinement");
}

} // namespace

int main() {
    try {
        verify_geometry_and_polynomials();
        verify_gaussian(PaperCase::R2a);
        verify_gaussian(PaperCase::R2b);
        verify_singular_convergence();
        std::cout << "Helmholtz unified quadrature policy passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_quadrature failed: " << error.what() << '\n';
        return 1;
    }
}

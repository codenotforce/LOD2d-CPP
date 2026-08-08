#pragma once

#include "helmholtz/experiments/paper_config.h"
#include "helmholtz/manufactured.h"
#include "mesh/types.h"

#include <optional>

namespace lod2d::helmholtz::benchmarks {

struct PaperCaseData {
    experiments::PaperCase id = experiments::PaperCase::R1;
    double wavenumber = 0.0;
    TriMesh initial_mesh;
    ComplexFunction source;
    ComplexFunction exact;
    ComplexGradientFunction exact_gradient;
    ComplexFunction exact_laplacian;
    QuadratureContext quadrature_context;
    std::optional<double> gaussian_sigma;
    std::optional<double> gaussian_normalization;
};

PaperCaseData make_paper_case(
    experiments::PaperCase id,
    double wavenumber);

double normalized_gaussian_constant(
    double sigma,
    const Point2 &center = Point2(0.35, 0.55));

// Frozen C-infinity cut-off and its first two radial derivatives.
double singular_cutoff(double radius);
double singular_cutoff_prime(double radius);
double singular_cutoff_second(double radius);

} // namespace lod2d::helmholtz::benchmarks

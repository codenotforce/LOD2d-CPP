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
    std::optional<double> singular_oscillatory_fraction;
    std::optional<double> singular_cutoff_outer_radius;
    bool singular_quintic_cutoff = false;
    std::optional<double> smooth_wave_amplitude;
};

PaperCaseData make_paper_case(
    experiments::PaperCase id,
    double wavenumber,
    double singular_oscillatory_fraction = 1.0,
    double singular_cutoff_outer_radius = 0.5,
    bool singular_quintic_cutoff = false,
    double smooth_wave_amplitude = 0.0);

double normalized_gaussian_constant(
    double sigma,
    const Point2 &center = Point2(0.35, 0.55));

// Frozen C-infinity cut-off and its first two radial derivatives.
double singular_cutoff(double radius);
double singular_cutoff_prime(double radius);
double singular_cutoff_second(double radius);

} // namespace lod2d::helmholtz::benchmarks

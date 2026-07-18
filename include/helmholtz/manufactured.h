#pragma once

#include "helmholtz/types.h"

namespace lod2d::helmholtz {

struct HelmholtzManufacturedSolution {
    ComplexFunction value;
    ComplexGradientFunction gradient;
    ComplexFunction source;
};

HelmholtzManufacturedSolution make_polynomial_plane_wave_solution(
    double wavenumber);

} // namespace lod2d::helmholtz

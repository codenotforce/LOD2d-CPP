#include "helmholtz/shifted_laplacian.h"

#include <cmath>
#include <stdexcept>

namespace lod2d::helmholtz {

ComplexSparseMatrix build_shifted_helmholtz_operator(
    const HelmholtzPatchSystem &system,
    double epsilon) {
    if (!(epsilon >= 0.0) || !std::isfinite(epsilon))
        throw std::invalid_argument(
            "shifted Helmholtz epsilon must be finite and nonnegative");
    ComplexSparseMatrix shifted = system.stiffness.cast<Complex>();
    shifted -= Complex(
        system.wavenumber * system.wavenumber, epsilon)
        * system.mass.cast<Complex>();
    shifted -= Complex(0.0, system.wavenumber)
        * system.robin.cast<Complex>();
    shifted.makeCompressed();
    return shifted;
}

} // namespace lod2d::helmholtz

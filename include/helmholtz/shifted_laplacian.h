#pragma once

#include "helmholtz/patch_system.h"

namespace lod2d::helmholtz {

ComplexSparseMatrix build_shifted_helmholtz_operator(
    const HelmholtzPatchSystem &system,
    double epsilon);

} // namespace lod2d::helmholtz

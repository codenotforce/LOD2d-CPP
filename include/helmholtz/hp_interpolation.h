#pragma once

#include "fem/hp_triangle.h"

#include <Eigen/Sparse>

namespace lod2d::helmholtz {

struct HelmholtzHpInterpolation {
    Eigen::SparseMatrix<double> coarse_injection;
    Eigen::SparseMatrix<double> quasi_interpolation;
};

HelmholtzHpInterpolation build_helmholtz_hp_interpolation(
    const TriMesh &coarse,
    const HpTriSpace &fine_space,
    const Eigen::SparseMatrix<double> &fine_element_prolongation);

} // namespace lod2d::helmholtz

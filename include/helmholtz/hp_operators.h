#pragma once

#include "fem/hp_triangle.h"
#include "helmholtz/operators.h"

#include <Eigen/Sparse>
#include <vector>

namespace lod2d::helmholtz {

struct HelmholtzHpOperators {
    double wavenumber = 0.0;
    double boundary_beta = 1.0;
    std::vector<double> diffusion;
    std::vector<double> refractive_index;
    std::vector<ComplexMatrix> element_blocks;
    Eigen::SparseMatrix<double> stiffness;
    Eigen::SparseMatrix<double> mass;
    Eigen::SparseMatrix<double> boundary_mass;
    ComplexSparseMatrix system;
};

HelmholtzHpOperators assemble_helmholtz_hp_operators(
    const HpTriSpace &space,
    double wavenumber,
    const std::vector<double> &diffusion = {},
    const std::vector<double> &refractive_index = {},
    double boundary_beta = 1.0);

ComplexVector assemble_helmholtz_hp_load(
    const HpTriSpace &space,
    const ComplexFunction &source);

ComplexVector solve_helmholtz_hp_fem(
    const HelmholtzHpOperators &operators,
    const ComplexVector &load);

HelmholtzError compute_helmholtz_hp_error(
    const HpTriSpace &space,
    const ComplexVector &solution,
    double wavenumber,
    const ComplexFunction &exact,
    const ComplexGradientFunction &exact_gradient);

} // namespace lod2d::helmholtz

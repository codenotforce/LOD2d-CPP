#pragma once

#include "helmholtz/types.h"
#include <Eigen/Dense>
#include <functional>
#include <vector>

namespace lod2d::helmholtz {

using HelmholtzElementBlocks = std::vector<Eigen::Matrix3cd>;

struct HelmholtzOperators {
    double wavenumber = 0.0;
    HelmholtzElementBlocks element_blocks;
    Eigen::SparseMatrix<double> stiffness;
    Eigen::SparseMatrix<double> mass;
    Eigen::SparseMatrix<double> boundary_mass;
    ComplexSparseMatrix system;
};

struct HelmholtzError {
    double energy = 0.0;
    double l2 = 0.0;
};

HelmholtzOperators assemble_helmholtz_operators(
    const TriMesh &mesh,
    double wavenumber,
    const std::vector<double> &diffusion = {},
    const std::vector<double> &refractive_index = {},
    double boundary_beta = 1.0);

ComplexVector assemble_helmholtz_load(
    const TriMesh &mesh,
    const ComplexFunction &source);

ComplexVector solve_helmholtz_fem(
    const HelmholtzOperators &operators,
    const ComplexVector &load);

HelmholtzError compute_helmholtz_error(
    const TriMesh &mesh,
    const ComplexVector &solution,
    double wavenumber,
    const ComplexFunction &exact,
    const ComplexGradientFunction &exact_gradient);

double max_element_diameter(const TriMesh &mesh);

} // namespace lod2d::helmholtz

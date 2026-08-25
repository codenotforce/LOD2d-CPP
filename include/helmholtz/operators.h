#pragma once

#include "helmholtz/types.h"
#include "helmholtz/quadrature.h"
#include <Eigen/Dense>
#include <functional>
#include <vector>

namespace lod2d::helmholtz {

using HelmholtzElementBlocks = std::vector<Eigen::Matrix3cd>;

struct HelmholtzOperators {
    double wavenumber = 0.0;
    HelmholtzElementBlocks element_blocks;
    std::vector<double> diffusion;
    std::vector<double> refractive_index;
    double boundary_beta = 1.0;
    Eigen::SparseMatrix<double> stiffness;
    Eigen::SparseMatrix<double> mass;
    Eigen::SparseMatrix<double> boundary_mass;
    ComplexSparseMatrix system;
    std::vector<int> dirichlet_nodes;
};

struct HelmholtzError {
    double energy = 0.0;
    double l2 = 0.0;
};

enum class HelmholtzFemSolverKind {
    SparseLu,
    Umfpack,
};

struct HelmholtzFemSolveTimings {
    double reduction_seconds = 0.0;
    double analysis_seconds = 0.0;
    double factorization_seconds = 0.0;
    double solve_seconds = 0.0;
    double total_seconds = 0.0;
};

const char *helmholtz_fem_solver_kind_name(HelmholtzFemSolverKind kind);
bool helmholtz_fem_solver_available(HelmholtzFemSolverKind kind);

HelmholtzOperators assemble_helmholtz_operators(
    const TriMesh &mesh,
    double wavenumber,
    const std::vector<double> &diffusion = {},
    const std::vector<double> &refractive_index = {},
    double boundary_beta = 1.0);

ComplexVector assemble_helmholtz_load(
    const TriMesh &mesh,
    const ComplexFunction &source,
    const QuadraturePolicy &quadrature = {},
    const QuadratureContext &quadrature_context = {});

ComplexVector solve_helmholtz_fem(
    const HelmholtzOperators &operators,
    const ComplexVector &load,
    HelmholtzFemSolverKind solver_kind = HelmholtzFemSolverKind::SparseLu,
    HelmholtzFemSolveTimings *timings = nullptr);

HelmholtzError compute_helmholtz_error(
    const TriMesh &mesh,
    const ComplexVector &solution,
    double wavenumber,
    const ComplexFunction &exact,
    const ComplexGradientFunction &exact_gradient,
    const QuadraturePolicy &quadrature = {},
    const QuadratureContext &quadrature_context = {});

double max_element_diameter(const TriMesh &mesh);

} // namespace lod2d::helmholtz

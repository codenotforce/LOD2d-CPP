#pragma once

#include "helmholtz/adaptive/hierarchy.h"
#include "helmholtz/model.h"

#include <array>
#include <vector>

namespace lod2d::helmholtz::adaptive {

// Strong/broken residual quantities are retained strictly as diagnostics.
// The paper estimator eta_H is declared in kernel_residual.h.
namespace diagnostics {

enum class ResidualEstimatorKind {
    Fine,
    Mixed,
    Macro
};

struct ResidualEdgeContribution {
    Edge nodes{-1, -1};
    int left_element = -1;
    int right_element = -1;
    int left_parent = -1;
    int right_parent = -1;
    bool robin_boundary = false;
    bool neumann_boundary = false;
    double length = 0.0;
    double residual_l2_squared = 0.0;
    std::array<Complex, 2> residual_nodal{Complex(0.0), Complex(0.0)};
};

struct HelmholtzResidualContributions {
    std::vector<int> fine_element_parent;
    std::vector<double> body_l2_squared;
    std::vector<std::array<Complex, 3>> body_residual_nodal;
    std::vector<ResidualEdgeContribution> edges;
    ComplexVector reconstructed_residual;
    double algebraic_relative_difference = 0.0;
};

struct HelmholtzIndicatorSet {
    std::vector<double> fine_squared;
    std::vector<double> mixed_squared;
    std::vector<double> macro_squared;
    double fine = 0.0;
    double mixed = 0.0;
    double macro = 0.0;

    const std::vector<double> &squared(ResidualEstimatorKind kind) const;
    double global(ResidualEstimatorKind kind) const;
};

// Standard conforming P1 AFEM estimator on one mesh:
//   h_T^2 ||f + kappa^2 n u_H||_T^2
// + h_E ||[A grad u_H . n]||_E^2
// + h_E ||A grad u_H . n||_E^2 on homogeneous Neumann edges
// + h_E ||A grad u_H . n - i kappa beta u_H||_E^2 on Robin edges.
// Interior-edge contributions are split equally between their neighbours.
struct HelmholtzP1ResidualEstimate {
    std::vector<double> body_squared;
    std::vector<double> interior_jump_squared;
    std::vector<double> neumann_boundary_squared;
    std::vector<double> robin_boundary_squared;
    std::vector<double> element_squared;
    double eta = 0.0;
    double algebraic_relative_difference = 0.0;
};

HelmholtzResidualContributions assemble_helmholtz_residual_contributions(
    const HelmholtzProblemData &problem,
    const HelmholtzOperators &operators,
    const ComplexVector &solution,
    const ComplexVector &load,
    const ComplexFunction &source,
    const QuadraturePolicy &quadrature = {},
    const QuadratureContext &quadrature_context = {});

HelmholtzIndicatorSet build_helmholtz_indicators(
    const HelmholtzProblemData &problem,
    const HelmholtzResidualContributions &contributions);

HelmholtzP1ResidualEstimate estimate_conforming_p1_residual(
    const TriMesh &mesh,
    const HelmholtzOperators &operators,
    const ComplexVector &solution,
    const ComplexVector &load,
    const ComplexFunction &source,
    const QuadraturePolicy &quadrature = {},
    const QuadratureContext &quadrature_context = {});

// Historical calibration diagnostic. Unlike eta_H, this partitions a broken
// residual by coarse-element ownership and must never drive paper MARK/STOP.
std::vector<double> build_local_dual_indicators(
    const HelmholtzProblemData &problem,
    const HelmholtzOperators &operators,
    const HelmholtzResidualContributions &contributions,
    int patch_layers);

} // namespace diagnostics

std::vector<int> mark_doerfler(
    const std::vector<double> &indicator_squared,
    double theta,
    const std::vector<char> &eligible = {});

double discrete_energy_norm(
    const HelmholtzOperators &operators,
    const ComplexVector &values);

double discrete_l2_norm(
    const HelmholtzOperators &operators,
    const ComplexVector &values);

} // namespace lod2d::helmholtz::adaptive

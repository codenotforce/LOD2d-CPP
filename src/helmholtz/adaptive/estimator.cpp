#include "helmholtz/adaptive/estimator.h"
#include "helmholtz/boundary.h"

#include "lod/patches.h"

#include <Eigen/SparseCholesky>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <unordered_map>

namespace lod2d::helmholtz::adaptive {
namespace {

struct ElementGeometry {
    double area = 0.0;
    double diameter = 0.0;
    std::array<Eigen::Vector2d, 3> gradients;
};

struct EdgeIncident {
    int element = -1;
    int opposite_node = -1;
};

std::uint64_t edge_key(int a, int b) {
    const auto lo = static_cast<std::uint32_t>(std::min(a, b));
    const auto hi = static_cast<std::uint32_t>(std::max(a, b));
    return (static_cast<std::uint64_t>(lo) << 32U) | hi;
}

Edge edge_from_key(std::uint64_t key) {
    return Edge{
        static_cast<int>(key >> 32U),
        static_cast<int>(key & 0xffffffffU)};
}

ElementGeometry element_geometry(const TriMesh &mesh, const Triangle &tri) {
    const Point2 &a = mesh.nodes[tri[0]];
    const Point2 &b = mesh.nodes[tri[1]];
    const Point2 &c = mesh.nodes[tri[2]];
    const double det = (b.x() - a.x()) * (c.y() - a.y())
                     - (b.y() - a.y()) * (c.x() - a.x());
    if (std::abs(det) <= 1e-15)
        throw std::invalid_argument("residual estimator encountered a degenerate triangle");
    ElementGeometry result;
    result.area = 0.5 * std::abs(det);
    result.gradients[0] = Eigen::Vector2d(b.y() - c.y(), c.x() - b.x()) / det;
    result.gradients[1] = Eigen::Vector2d(c.y() - a.y(), a.x() - c.x()) / det;
    result.gradients[2] = Eigen::Vector2d(a.y() - b.y(), b.x() - a.x()) / det;
    result.diameter = std::max({(a - b).norm(), (b - c).norm(), (c - a).norm()});
    return result;
}

Eigen::Vector2d outward_normal(
    const TriMesh &mesh,
    const Edge &edge,
    int opposite_node) {
    const Point2 &a = mesh.nodes[edge[0]];
    const Point2 &b = mesh.nodes[edge[1]];
    const Point2 tangent = b - a;
    const double length = tangent.norm();
    if (length <= 0.0) throw std::runtime_error("residual estimator found a zero-length edge");
    Eigen::Vector2d normal(tangent.y(), -tangent.x());
    normal /= length;
    const Point2 midpoint = 0.5 * (a + b);
    if (normal.dot(mesh.nodes[opposite_node] - midpoint) > 0.0) normal = -normal;
    return normal;
}

double quadratic_form(
    const Eigen::SparseMatrix<double> &matrix,
    const ComplexVector &values) {
    if (matrix.rows() != values.size())
        throw std::invalid_argument("quadratic form dimensions do not match");
    return std::max(0.0, std::real(values.dot(matrix.cast<Complex>() * values)));
}

double indicator_global(const std::vector<double> &squared) {
    return std::sqrt(std::max(0.0, std::accumulate(squared.begin(), squared.end(), 0.0)));
}

ComplexVector coarse_partition_residual(
    int coarse_element,
    int fine_node_count,
    const TriMesh &fine,
    const diagnostics::HelmholtzResidualContributions &contributions) {
    ComplexVector rhs = ComplexVector::Zero(fine_node_count);
    for (int element = 0; element < static_cast<int>(fine.elems.size()); ++element) {
        if (contributions.fine_element_parent[element] != coarse_element) continue;
        for (int local = 0; local < 3; ++local)
            rhs(fine.elems[element][local]) += contributions.body_residual_nodal[element][local];
    }
    for (const diagnostics::ResidualEdgeContribution &edge : contributions.edges) {
        double share = 0.0;
        if (edge.right_element < 0) {
            if (edge.left_parent == coarse_element) share = 1.0;
        } else if (edge.left_parent == edge.right_parent) {
            if (edge.left_parent == coarse_element) share = 1.0;
        } else if (edge.left_parent == coarse_element || edge.right_parent == coarse_element) {
            share = 0.5;
        }
        if (share == 0.0) continue;
        rhs(edge.nodes[0]) += share * edge.residual_nodal[0];
        rhs(edge.nodes[1]) += share * edge.residual_nodal[1];
    }
    return rhs;
}

} // namespace

namespace diagnostics {

const std::vector<double> &HelmholtzIndicatorSet::squared(ResidualEstimatorKind kind) const {
    if (kind == ResidualEstimatorKind::Fine) return fine_squared;
    if (kind == ResidualEstimatorKind::Mixed) return mixed_squared;
    return macro_squared;
}

double HelmholtzIndicatorSet::global(ResidualEstimatorKind kind) const {
    if (kind == ResidualEstimatorKind::Fine) return fine;
    if (kind == ResidualEstimatorKind::Mixed) return mixed;
    return macro;
}

HelmholtzResidualContributions assemble_helmholtz_residual_contributions(
    const HelmholtzProblemData &problem,
    const HelmholtzOperators &operators,
    const ComplexVector &solution,
    const ComplexVector &load,
    const ComplexFunction &source,
    const QuadraturePolicy &quadrature,
    const QuadratureContext &quadrature_context) {
    const TriMesh &fine = problem.fine;
    const int element_count = static_cast<int>(fine.elems.size());
    const int node_count = static_cast<int>(fine.nodes.size());
    if (!source) throw std::invalid_argument("residual estimator source is empty");
    if (solution.size() != node_count || load.size() != node_count)
        throw std::invalid_argument("residual estimator vector dimensions do not match the fine mesh");
    if (operators.diffusion.size() != fine.elems.size()
        || operators.refractive_index.size() != fine.elems.size()) {
        throw std::invalid_argument("residual estimator coefficient counts do not match fine elements");
    }

    HelmholtzResidualContributions result;
    result.fine_element_parent = fine_element_parents(
        problem.fine_element_prolongation,
        element_count,
        static_cast<int>(problem.coarse.elems.size()));
    result.body_l2_squared.assign(element_count, 0.0);
    result.body_residual_nodal.assign(
        element_count,
        std::array<Complex, 3>{Complex(0.0), Complex(0.0), Complex(0.0)});
    result.reconstructed_residual = ComplexVector::Zero(node_count);

    std::vector<ElementGeometry> geometry(element_count);
    std::vector<Eigen::Vector2cd> gradients(element_count, Eigen::Vector2cd::Zero());
    std::unordered_map<std::uint64_t, std::vector<EdgeIncident>> edge_incidence;
    edge_incidence.reserve(3 * fine.elems.size());
    static constexpr int local_edges[3][3] = {{0, 1, 2}, {1, 2, 0}, {2, 0, 1}};

    for (int element = 0; element < element_count; ++element) {
        const Triangle &tri = fine.elems[element];
        geometry[element] = element_geometry(fine, tri);
        for (int local = 0; local < 3; ++local)
            gradients[element] += solution(tri[local]) * geometry[element].gradients[local].cast<Complex>();
        for (const auto &local_edge : local_edges) {
            edge_incidence[edge_key(tri[local_edge[0]], tri[local_edge[1]])].push_back(
                {element, tri[local_edge[2]]});
        }

        for (const auto &point : triangle_quadrature_points(
                 fine, element, quadrature, quadrature_context)) {
            Complex value = 0.0;
            for (int local = 0; local < 3; ++local) {
                value += point.barycentric[local] * solution(tri[local]);
            }
            const Complex residual = source(point.point)
                + operators.wavenumber * operators.wavenumber
                * operators.refractive_index[element] * value;
            result.body_l2_squared[element] += point.weight * std::norm(residual);
            for (int local = 0; local < 3; ++local) {
                result.body_residual_nodal[element][local]
                    += point.weight * residual * point.barycentric[local];
            }
        }
        for (int local = 0; local < 3; ++local)
            result.reconstructed_residual(tri[local]) += result.body_residual_nodal[element][local];
    }

    result.edges.reserve(edge_incidence.size());
    for (const auto &[key, incidents] : edge_incidence) {
        if (incidents.empty() || incidents.size() > 2)
            throw std::runtime_error("residual estimator found invalid edge incidence");
        ResidualEdgeContribution edge;
        edge.nodes = edge_from_key(key);
        edge.left_element = incidents[0].element;
        edge.left_parent = result.fine_element_parent[edge.left_element];
        edge.length = (fine.nodes[edge.nodes[1]] - fine.nodes[edge.nodes[0]]).norm();
        const Eigen::Vector2d left_normal = outward_normal(
            fine, edge.nodes, incidents[0].opposite_node);
        const Complex left_flux = operators.diffusion[edge.left_element]
                                * left_normal.cast<Complex>().dot(gradients[edge.left_element]);

        if (incidents.size() == 2) {
            edge.right_element = incidents[1].element;
            edge.right_parent = result.fine_element_parent[edge.right_element];
            const Eigen::Vector2d right_normal = outward_normal(
                fine, edge.nodes, incidents[1].opposite_node);
            const Complex right_flux = operators.diffusion[edge.right_element]
                                     * right_normal.cast<Complex>().dot(gradients[edge.right_element]);
            const Complex jump = left_flux + right_flux;
            edge.residual_l2_squared = edge.length * std::norm(jump);
            edge.residual_nodal[0] = -0.5 * edge.length * jump;
            edge.residual_nodal[1] = -0.5 * edge.length * jump;
        } else if (boundary_tag(fine, canonical_edge(edge.nodes[0], edge.nodes[1]))
                   == BoundaryTag::Robin) {
            edge.robin_boundary = true;
            const Complex impedance(0.0, operators.wavenumber * operators.boundary_beta);
            const Complex boundary0 = left_flux - impedance * solution(edge.nodes[0]);
            const Complex boundary1 = left_flux - impedance * solution(edge.nodes[1]);
            edge.residual_l2_squared = edge.length / 3.0 * (
                std::norm(boundary0) + std::norm(boundary1)
                + std::real(boundary0 * std::conj(boundary1)));
            const Complex residual0 = -boundary0;
            const Complex residual1 = -boundary1;
            edge.residual_nodal[0] = edge.length * (2.0 * residual0 + residual1) / 6.0;
            edge.residual_nodal[1] = edge.length * (residual0 + 2.0 * residual1) / 6.0;
        }
        result.reconstructed_residual(edge.nodes[0]) += edge.residual_nodal[0];
        result.reconstructed_residual(edge.nodes[1]) += edge.residual_nodal[1];
        result.edges.push_back(edge);
    }

    ComplexVector algebraic = load - operators.system * solution;
    for (int node : dirichlet_nodes(fine)) {
        result.reconstructed_residual(node) = Complex(0.0, 0.0);
        algebraic(node) = Complex(0.0, 0.0);
    }
    result.algebraic_relative_difference = (result.reconstructed_residual - algebraic).norm()
        / std::max(1.0, algebraic.norm());
    return result;
}

HelmholtzIndicatorSet build_helmholtz_indicators(
    const HelmholtzProblemData &problem,
    const HelmholtzResidualContributions &contributions) {
    const int coarse_count = static_cast<int>(problem.coarse.elems.size());
    const int fine_count = static_cast<int>(problem.fine.elems.size());
    if (contributions.body_l2_squared.size() != problem.fine.elems.size()
        || contributions.fine_element_parent.size() != problem.fine.elems.size()) {
        throw std::invalid_argument("residual contribution count does not match the fine mesh");
    }

    std::vector<double> coarse_diameter(coarse_count);
    for (int element = 0; element < coarse_count; ++element)
        coarse_diameter[element] = element_geometry(problem.coarse, problem.coarse.elems[element]).diameter;

    HelmholtzIndicatorSet result;
    result.fine_squared.assign(coarse_count, 0.0);
    result.mixed_squared.assign(coarse_count, 0.0);
    result.macro_squared.assign(coarse_count, 0.0);
    for (int element = 0; element < fine_count; ++element) {
        const int parent = contributions.fine_element_parent[element];
        const double h = element_geometry(problem.fine, problem.fine.elems[element]).diameter;
        const double H = coarse_diameter[parent];
        const double body = contributions.body_l2_squared[element];
        result.fine_squared[parent] += h * h * body;
        result.mixed_squared[parent] += H * H * body;
        result.macro_squared[parent] += H * H * body;
    }

    for (const ResidualEdgeContribution &edge : contributions.edges) {
        auto add = [&](int parent, double share) {
            result.fine_squared[parent] += share * edge.length * edge.residual_l2_squared;
            result.mixed_squared[parent] += share * edge.length * edge.residual_l2_squared;
            result.macro_squared[parent] += share * coarse_diameter[parent]
                                          * edge.residual_l2_squared;
        };
        if (edge.right_element < 0 || edge.left_parent == edge.right_parent) {
            add(edge.left_parent, 1.0);
        } else {
            add(edge.left_parent, 0.5);
            add(edge.right_parent, 0.5);
        }
    }

    result.fine = indicator_global(result.fine_squared);
    result.mixed = indicator_global(result.mixed_squared);
    result.macro = indicator_global(result.macro_squared);
    return result;
}

HelmholtzP1ResidualEstimate estimate_conforming_p1_residual(
    const TriMesh &mesh,
    const HelmholtzOperators &operators,
    const ComplexVector &solution,
    const ComplexVector &load,
    const ComplexFunction &source,
    const QuadraturePolicy &quadrature,
    const QuadratureContext &quadrature_context) {
    HelmholtzProblemData identity_problem;
    identity_problem.coarse = mesh;
    identity_problem.fine = mesh;
    const int element_count = static_cast<int>(mesh.elems.size());
    identity_problem.fine_element_prolongation.resize(
        element_count, element_count);
    identity_problem.fine_element_prolongation.setIdentity();

    const HelmholtzResidualContributions contributions =
        assemble_helmholtz_residual_contributions(
            identity_problem, operators, solution, load, source,
            quadrature, quadrature_context);
    const HelmholtzIndicatorSet total =
        build_helmholtz_indicators(identity_problem, contributions);

    HelmholtzP1ResidualEstimate result;
    result.body_squared.assign(element_count, 0.0);
    result.interior_jump_squared.assign(element_count, 0.0);
    result.robin_boundary_squared.assign(element_count, 0.0);
    result.element_squared.assign(element_count, 0.0);
    for (int element = 0; element < element_count; ++element) {
        const double diameter =
            element_geometry(mesh, mesh.elems[element]).diameter;
        result.body_squared[element] = diameter * diameter
            * contributions.body_l2_squared[element];
    }
    for (const ResidualEdgeContribution &edge : contributions.edges) {
        const double weighted = edge.length * edge.residual_l2_squared;
        if (edge.right_element >= 0) {
            result.interior_jump_squared[edge.left_element] += 0.5 * weighted;
            result.interior_jump_squared[edge.right_element] += 0.5 * weighted;
        } else if (edge.robin_boundary) {
            result.robin_boundary_squared[edge.left_element] += weighted;
        }
    }
    for (int element = 0; element < element_count; ++element) {
        result.element_squared[element] = result.body_squared[element]
            + result.interior_jump_squared[element]
            + result.robin_boundary_squared[element];
    }
    const double component_sum = std::accumulate(
        result.element_squared.begin(), result.element_squared.end(), 0.0);
    const double total_sum = std::accumulate(
        total.fine_squared.begin(), total.fine_squared.end(), 0.0);
    const double consistency_error = std::abs(component_sum - total_sum)
        / std::max(1.0, total_sum);
    if (consistency_error > 5e-13) {
        throw std::runtime_error(
            "conforming P1 residual component allocation is inconsistent");
    }
    result.eta = std::sqrt(std::max(0.0, component_sum));
    result.algebraic_relative_difference =
        contributions.algebraic_relative_difference;
    return result;
}

std::vector<double> build_local_dual_indicators(
    const HelmholtzProblemData &problem,
    const HelmholtzOperators &operators,
    const HelmholtzResidualContributions &contributions,
    int patch_layers) {
    if (patch_layers < 0) throw std::invalid_argument("dual residual patch layers must be nonnegative");
    const int coarse_count = static_cast<int>(problem.coarse.elems.size());
    const int fine_node_count = static_cast<int>(problem.fine.nodes.size());
    const Eigen::SparseMatrix<double> patches = build_patches(problem.coarse, patch_layers);
    Eigen::SparseMatrix<double> energy = operators.stiffness;
    energy += operators.wavenumber * operators.wavenumber * operators.mass;

    std::vector<std::vector<int>> children(coarse_count);
    for (int fine_element = 0; fine_element < static_cast<int>(problem.fine.elems.size()); ++fine_element)
        children[contributions.fine_element_parent[fine_element]].push_back(fine_element);
    std::vector<int> global_incidence(fine_node_count, 0);
    for (const Triangle &tri : problem.fine.elems)
        for (int node : tri) ++global_incidence[node];

    std::vector<double> indicators(coarse_count, 0.0);
    for (int target = 0; target < coarse_count; ++target) {
        std::vector<int> patch_elements;
        for (Eigen::SparseMatrix<double>::InnerIterator it(patches, target); it; ++it) {
            patch_elements.insert(
                patch_elements.end(), children[it.row()].begin(), children[it.row()].end());
        }
        std::vector<int> patch_incidence(fine_node_count, 0);
        std::vector<int> touched;
        for (int element : patch_elements) {
            for (int node : problem.fine.elems[element]) {
                if (patch_incidence[node]++ == 0) touched.push_back(node);
            }
        }
        std::vector<int> local_nodes;
        std::vector<int> local_index(fine_node_count, -1);
        for (int node : touched) {
            if (patch_incidence[node] != global_incidence[node]) continue;
            local_index[node] = static_cast<int>(local_nodes.size());
            local_nodes.push_back(node);
        }
        if (local_nodes.empty()) continue;

        std::vector<Eigen::Triplet<double>> triplets;
        for (int global_col : local_nodes) {
            const int local_col = local_index[global_col];
            for (Eigen::SparseMatrix<double>::InnerIterator it(energy, global_col); it; ++it) {
                const int local_row = local_index[it.row()];
                if (local_row >= 0) triplets.emplace_back(local_row, local_col, it.value());
            }
        }
        Eigen::SparseMatrix<double> local_energy(local_nodes.size(), local_nodes.size());
        local_energy.setFromTriplets(triplets.begin(), triplets.end());
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(local_energy);
        if (solver.info() != Eigen::Success)
            throw std::runtime_error("local dual residual Riesz factorization failed");

        const ComplexVector global_rhs = coarse_partition_residual(
            target, fine_node_count, problem.fine, contributions);
        Eigen::VectorXd rhs_real(local_nodes.size());
        Eigen::VectorXd rhs_imag(local_nodes.size());
        for (int local = 0; local < static_cast<int>(local_nodes.size()); ++local) {
            rhs_real(local) = global_rhs(local_nodes[local]).real();
            rhs_imag(local) = global_rhs(local_nodes[local]).imag();
        }
        const Eigen::VectorXd solution_real = solver.solve(rhs_real);
        const Eigen::VectorXd solution_imag = solver.solve(rhs_imag);
        if (solver.info() != Eigen::Success)
            throw std::runtime_error("local dual residual Riesz solve failed");
        const double squared = rhs_real.dot(solution_real) + rhs_imag.dot(solution_imag);
        indicators[target] = std::sqrt(std::max(0.0, squared));
    }
    return indicators;
}

} // namespace diagnostics

std::vector<int> mark_doerfler(
    const std::vector<double> &indicator_squared,
    double theta,
    const std::vector<char> &eligible) {
    if (!(theta > 0.0 && theta <= 1.0))
        throw std::invalid_argument("Doerfler theta must lie in (0,1]");
    if (!eligible.empty() && eligible.size() != indicator_squared.size())
        throw std::invalid_argument("Doerfler eligibility count does not match indicators");
    for (double value : indicator_squared) {
        if (!(value >= 0.0) || !std::isfinite(value))
            throw std::invalid_argument("Doerfler indicators must be finite and nonnegative");
    }

    std::vector<int> order;
    double eligible_sum = 0.0;
    const double total = std::accumulate(indicator_squared.begin(), indicator_squared.end(), 0.0);
    for (int i = 0; i < static_cast<int>(indicator_squared.size()); ++i) {
        if (eligible.empty() || eligible[i]) {
            order.push_back(i);
            eligible_sum += indicator_squared[i];
        }
    }
    const double target = theta * total;
    if (eligible_sum + 1e-14 * std::max(1.0, total) < target)
        throw std::runtime_error("eligible coarse elements cannot satisfy Doerfler marking");
    std::sort(order.begin(), order.end(), [&](int left, int right) {
        if (indicator_squared[left] != indicator_squared[right])
            return indicator_squared[left] > indicator_squared[right];
        return left < right;
    });

    std::vector<int> marked;
    double accumulated = 0.0;
    for (int element : order) {
        if (accumulated >= target && !marked.empty()) break;
        marked.push_back(element);
        accumulated += indicator_squared[element];
    }
    return marked;
}

double discrete_energy_norm(
    const HelmholtzOperators &operators,
    const ComplexVector &values) {
    return std::sqrt(quadratic_form(operators.stiffness, values)
        + operators.wavenumber * operators.wavenumber
        * quadratic_form(operators.mass, values));
}

double discrete_l2_norm(
    const HelmholtzOperators &operators,
    const ComplexVector &values) {
    return std::sqrt(quadratic_form(operators.mass, values));
}

} // namespace lod2d::helmholtz::adaptive

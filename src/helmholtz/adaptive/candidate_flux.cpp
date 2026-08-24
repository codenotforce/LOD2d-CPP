#include "helmholtz/adaptive/candidate_flux.h"

#include "helmholtz/boundary.h"
#include "mesh/refine.h"

#include <Eigen/Cholesky>
#include <Eigen/LU>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <exception>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lod2d::helmholtz::adaptive {
namespace {

struct ElementGeometry {
    double area = 0.0;
    std::array<Eigen::Vector2d, 3> gradients;
};

ElementGeometry element_geometry(const TriMesh &mesh, int element) {
    const Triangle &triangle = mesh.elems[element];
    const Point2 &p0 = mesh.nodes[triangle[0]];
    const Point2 &p1 = mesh.nodes[triangle[1]];
    const Point2 &p2 = mesh.nodes[triangle[2]];
    const double determinant =
        (p1.x() - p0.x()) * (p2.y() - p0.y())
        - (p1.y() - p0.y()) * (p2.x() - p0.x());
    if (std::abs(determinant) <= 1e-15)
        throw std::invalid_argument("candidate mesh contains a degenerate triangle");
    ElementGeometry result;
    result.area = 0.5 * std::abs(determinant);
    result.gradients[0] = Eigen::Vector2d(
        p1.y() - p2.y(), p2.x() - p1.x()) / determinant;
    result.gradients[1] = Eigen::Vector2d(
        p2.y() - p0.y(), p0.x() - p2.x()) / determinant;
    result.gradients[2] = Eigen::Vector2d(
        p0.y() - p1.y(), p1.x() - p0.x()) / determinant;
    return result;
}

struct EdgeIncidence {
    int element = -1;
    int opposite_vertex = -1;
    double orientation = 0.0;
};

struct EdgeData {
    Edge nodes{};
    std::vector<EdgeIncidence> incidences;
};

struct MeshTopology {
    std::vector<EdgeData> edges;
    std::vector<std::array<int, 3>> element_edges;
    std::vector<std::array<double, 3>> element_orientations;
};

MeshTopology build_topology(const TriMesh &mesh) {
    MeshTopology result;
    result.element_edges.resize(mesh.elems.size());
    result.element_orientations.resize(mesh.elems.size());
    std::map<Edge, int> edge_index;
    for (int element = 0; element < static_cast<int>(mesh.elems.size());
         ++element) {
        const Triangle &triangle = mesh.elems[element];
        const Point2 centroid = (
            mesh.nodes[triangle[0]] + mesh.nodes[triangle[1]]
            + mesh.nodes[triangle[2]]) / 3.0;
        for (int opposite = 0; opposite < 3; ++opposite) {
            const int first = triangle[(opposite + 1) % 3];
            const int second = triangle[(opposite + 2) % 3];
            const Edge edge = canonical_edge(first, second);
            auto [found, inserted] = edge_index.emplace(
                edge, static_cast<int>(result.edges.size()));
            if (inserted) result.edges.push_back({edge, {}});
            const int index = found->second;
            const Point2 direction =
                mesh.nodes[edge[1]] - mesh.nodes[edge[0]];
            const double length = direction.norm();
            if (!(length > 0.0))
                throw std::invalid_argument("candidate mesh contains a zero edge");
            const Eigen::Vector2d canonical_normal(
                direction.y() / length, -direction.x() / length);
            const Point2 midpoint =
                0.5 * (mesh.nodes[edge[0]] + mesh.nodes[edge[1]]);
            const double orientation =
                canonical_normal.dot(midpoint - centroid) > 0.0 ? 1.0 : -1.0;
            result.element_edges[element][opposite] = index;
            result.element_orientations[element][opposite] = orientation;
            result.edges[index].incidences.push_back(
                {element, opposite, orientation});
        }
    }
    return result;
}

Eigen::Matrix3d rt0_flux_to_coefficients(
    const TriMesh &mesh,
    int element,
    const MeshTopology &topology) {
    Eigen::Matrix3d flux_matrix;
    for (int opposite = 0; opposite < 3; ++opposite) {
        const Edge &edge = topology.edges[
            topology.element_edges[element][opposite]].nodes;
        const Point2 direction =
            mesh.nodes[edge[1]] - mesh.nodes[edge[0]];
        const double length = direction.norm();
        Eigen::Vector2d normal(direction.y(), -direction.x());
        normal /= length;
        normal *= topology.element_orientations[element][opposite];
        const Point2 midpoint =
            0.5 * (mesh.nodes[edge[0]] + mesh.nodes[edge[1]]);
        flux_matrix.row(opposite) <<
            length * normal.x(), length * normal.y(),
            length * normal.dot(midpoint);
    }
    Eigen::FullPivLU<Eigen::Matrix3d> factorization(flux_matrix);
    if (!factorization.isInvertible())
        throw std::runtime_error("RT0 edge-flux map is singular");
    return factorization.inverse();
}

Eigen::Matrix3d rt0_mass_matrix(
    const TriMesh &mesh,
    int element,
    const Eigen::Matrix3d &flux_to_coefficients,
    const CandidateFluxConfig &config) {
    Eigen::Matrix3d result = Eigen::Matrix3d::Zero();
    for (const PhysicalTriangleQuadraturePoint &point :
         triangle_quadrature_points(
             mesh, element, config.quadrature,
             config.quadrature_context)) {
        Eigen::Matrix<double, 2, 3> basis;
        for (int edge = 0; edge < 3; ++edge) {
            const Eigen::Vector3d coefficient =
                flux_to_coefficients.col(edge);
            basis.col(edge) = Eigen::Vector2d(
                coefficient(0) + coefficient(2) * point.point.x(),
                coefficient(1) + coefficient(2) * point.point.y());
        }
        result += point.weight * basis.transpose() * basis;
    }
    return result;
}

Eigen::Vector2cd element_gradient(
    const TriMesh &mesh,
    int element,
    const ComplexVector &values,
    const ElementGeometry &geometry) {
    Eigen::Vector2cd result = Eigen::Vector2cd::Zero();
    for (int local = 0; local < 3; ++local)
        result += values(mesh.elems[element][local])
            * geometry.gradients[local].cast<Complex>();
    return result;
}

Complex edge_hat_solution_integral(
    const TriMesh &mesh,
    const Edge &edge,
    int patch_vertex,
    const ComplexVector &values) {
    static constexpr std::array<double, 3> points{
        0.1127016653792583, 0.5, 0.8872983346207417};
    static constexpr std::array<double, 3> weights{
        5.0 / 18.0, 8.0 / 18.0, 5.0 / 18.0};
    const double length =
        (mesh.nodes[edge[1]] - mesh.nodes[edge[0]]).norm();
    Complex result(0.0, 0.0);
    for (int q = 0; q < 3; ++q) {
        const double t = points[q];
        const double psi = patch_vertex == edge[0]
            ? 1.0 - t : (patch_vertex == edge[1] ? t : 0.0);
        const Complex value =
            (1.0 - t) * values(edge[0]) + t * values(edge[1]);
        result += length * weights[q] * psi * value;
    }
    return result;
}

std::vector<int> doerfler_mark(
    const std::vector<double> &indicators,
    double theta,
    double &relative_error) {
    if (!(theta > 0.0 && theta <= 1.0))
        throw std::invalid_argument("candidate Doerfler theta must lie in (0,1]");
    std::vector<int> order(indicators.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](int left, int right) {
        return indicators[left] > indicators[right];
    });
    const double total = std::accumulate(
        indicators.begin(), indicators.end(), 0.0);
    if (!(total > 0.0)) {
        relative_error = 0.0;
        return {};
    }
    const double target = theta * total;
    double accumulated = 0.0;
    std::vector<int> result;
    for (int element : order) {
        if (!(indicators[element] > 0.0)) break;
        result.push_back(element);
        accumulated += indicators[element];
        if (accumulated >= target) break;
    }
    relative_error = std::max(0.0, target - accumulated)
        / std::max(1.0, total);
    return result;
}

double discrete_residual_dual_norm(
    const TriMesh &mesh,
    const HelmholtzOperators &operators,
    const ComplexFunction &source,
    const ComplexVector &values,
    const CandidateFluxConfig &config) {
    ComplexVector residual = assemble_helmholtz_load(
        mesh, source, config.quadrature, config.quadrature_context)
        - operators.system * values;
    std::vector<char> constrained(mesh.nodes.size(), false);
    for (int node : operators.dirichlet_nodes) constrained[node] = true;
    std::vector<int> index(mesh.nodes.size(), -1);
    int free_count = 0;
    for (int node = 0; node < static_cast<int>(mesh.nodes.size()); ++node)
        if (!constrained[node]) index[node] = free_count++;
    Eigen::SparseMatrix<double> energy = operators.stiffness
        + operators.wavenumber * operators.wavenumber * operators.mass;
    std::vector<Eigen::Triplet<double>> triplets;
    for (int column = 0; column < energy.outerSize(); ++column) {
        if (index[column] < 0) continue;
        for (Eigen::SparseMatrix<double>::InnerIterator it(energy, column);
             it; ++it) {
            if (index[it.row()] >= 0) {
                triplets.emplace_back(
                    index[it.row()], index[column], it.value());
            }
        }
    }
    Eigen::SparseMatrix<double> reduced(free_count, free_count);
    reduced.setFromTriplets(triplets.begin(), triplets.end());
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> factorization(reduced);
    if (factorization.info() != Eigen::Success)
        throw std::runtime_error("candidate residual energy factorization failed");
    Eigen::VectorXd rhs_real(free_count), rhs_imag(free_count);
    for (int node = 0; node < static_cast<int>(mesh.nodes.size()); ++node) {
        if (index[node] < 0) continue;
        rhs_real(index[node]) = residual(node).real();
        rhs_imag(index[node]) = residual(node).imag();
    }
    const Eigen::VectorXd solution_real = factorization.solve(rhs_real);
    const Eigen::VectorXd solution_imag = factorization.solve(rhs_imag);
    if (factorization.info() != Eigen::Success)
        throw std::runtime_error("candidate residual energy solve failed");
    return std::sqrt(std::max(
        0.0, rhs_real.dot(solution_real) + rhs_imag.dot(solution_imag)));
}

} // namespace

Eigen::Vector2cd evaluate_candidate_rt0_flux(
    const Eigen::Vector3cd &coefficients,
    const Point2 &point) {
    return Eigen::Vector2cd(
        coefficients(0) + coefficients(2) * point.x(),
        coefficients(1) + coefficients(2) * point.y());
}

CandidateFluxResult reconstruct_candidate_flux_rt0(
    const TriMesh &mesh,
    const HelmholtzOperators &operators,
    const ComplexFunction &source,
    const ComplexVector &values,
    const CandidateFluxConfig &config) {
    validate_boundary_tags(mesh);
    validate_quadrature_policy(config.quadrature);
    const int element_count = static_cast<int>(mesh.elems.size());
    const int node_count = static_cast<int>(mesh.nodes.size());
    if (!source || values.size() != node_count || !values.allFinite()
        || operators.system.rows() != node_count
        || operators.diffusion.size() != mesh.elems.size()
        || operators.refractive_index.size() != mesh.elems.size()
        || !(operators.wavenumber > 0.0)
        || !(config.doerfler_theta > 0.0
             && config.doerfler_theta <= 1.0)) {
        throw std::invalid_argument(
            "candidate RT0 reconstruction inputs are inconsistent");
    }

    const MeshTopology topology = build_topology(mesh);
    std::vector<ElementGeometry> geometry(element_count);
    std::vector<Eigen::Matrix3d> flux_to_coefficients(element_count);
    std::vector<Eigen::Matrix3d> rt0_mass(element_count);
    std::vector<Eigen::Vector2cd> gradients(element_count);
    std::vector<double> areas(element_count);
    std::vector<Complex> source_projection(element_count);
    std::vector<Complex> reaction_projection(element_count);
    for (int element = 0; element < element_count; ++element) {
        geometry[element] = element_geometry(mesh, element);
        areas[element] = geometry[element].area;
        flux_to_coefficients[element] = rt0_flux_to_coefficients(
            mesh, element, topology);
        rt0_mass[element] = rt0_mass_matrix(
            mesh, element, flux_to_coefficients[element], config);
        gradients[element] = element_gradient(
            mesh, element, values, geometry[element]);
        Complex source_integral(0.0, 0.0);
        Complex solution_integral(0.0, 0.0);
        for (const PhysicalTriangleQuadraturePoint &point :
             triangle_quadrature_points(
                 mesh, element, config.quadrature,
                 config.quadrature_context)) {
            Complex value(0.0, 0.0);
            for (int local = 0; local < 3; ++local)
                value += point.barycentric[local]
                    * values(mesh.elems[element][local]);
            source_integral += point.weight * source(point.point);
            solution_integral += point.weight * value;
        }
        source_projection[element] = source_integral / areas[element];
        reaction_projection[element] = source_projection[element]
            + operators.wavenumber * operators.wavenumber
                * operators.refractive_index[element]
                * solution_integral / areas[element];
    }

    std::vector<std::vector<int>> vertex_elements(node_count);
    for (int element = 0; element < element_count; ++element)
        for (int node : mesh.elems[element])
            vertex_elements[node].push_back(element);
    std::vector<char> active_vertex(node_count, config.active_elements.empty());
    if (!config.active_elements.empty()) {
        for (int element : config.active_elements) {
            if (element < 0 || element >= element_count)
                throw std::out_of_range("active candidate element is out of range");
            for (int node : mesh.elems[element]) active_vertex[node] = true;
        }
    }

    CandidateFluxResult result;
    result.global_reconstruction = config.active_elements.empty();
    result.element_rt0_coefficients.assign(
        element_count, Eigen::Vector3cd::Zero());
    result.projected_source = source_projection;
    result.delta_pg.assign(element_count, Complex(0.0, 0.0));

    for (int vertex = 0; vertex < node_count; ++vertex) {
        if (!active_vertex[vertex] || vertex_elements[vertex].empty()) continue;
        CandidateFluxPatchDiagnostics diagnostic;
        diagnostic.vertex = vertex;
        diagnostic.elements = vertex_elements[vertex];
        std::vector<char> in_patch(element_count, false);
        for (int element : diagnostic.elements) in_patch[element] = true;

        std::set<int> patch_edges;
        for (int element : diagnostic.elements)
            for (int edge : topology.element_edges[element])
                patch_edges.insert(edge);
        std::map<int, Complex> fixed_flux;
        std::vector<int> free_edges;
        std::map<int, int> free_index;
        Complex physical_flux_integral(0.0, 0.0);
        for (int edge_index : patch_edges) {
            const EdgeData &edge_data = topology.edges[edge_index];
            int patch_incidence = 0;
            const EdgeIncidence *patch_side = nullptr;
            for (const EdgeIncidence &incidence : edge_data.incidences) {
                if (in_patch[incidence.element]) {
                    ++patch_incidence;
                    patch_side = &incidence;
                }
            }
            if (patch_incidence == 2) {
                free_index[edge_index] = static_cast<int>(free_edges.size());
                free_edges.push_back(edge_index);
                continue;
            }
            if (patch_incidence != 1 || patch_side == nullptr)
                throw std::runtime_error("candidate vertex patch has invalid edge incidence");
            if (edge_data.incidences.size() == 2) {
                fixed_flux[edge_index] = Complex(0.0, 0.0);
                continue;
            }
            const BoundaryTag tag = boundary_tag(mesh, edge_data.nodes);
            if (tag == BoundaryTag::Interior)
                throw std::runtime_error(
                    "candidate flux found an unclassified physical boundary edge");
            if (tag == BoundaryTag::Dirichlet) {
                diagnostic.touches_dirichlet = true;
                free_index[edge_index] = static_cast<int>(free_edges.size());
                free_edges.push_back(edge_index);
                continue;
            }
            Complex local_flux(0.0, 0.0);
            if (tag == BoundaryTag::Robin) {
                local_flux = -Complex(0.0, operators.wavenumber
                    * operators.boundary_beta)
                    * edge_hat_solution_integral(
                        mesh, edge_data.nodes, vertex, values);
            } else if (tag != BoundaryTag::Neumann) {
                throw std::runtime_error("candidate flux found an unknown boundary tag");
            }
            fixed_flux[edge_index] = patch_side->orientation * local_flux;
            physical_flux_integral += local_flux;
        }

        const int patch_elements = static_cast<int>(diagnostic.elements.size());
        std::vector<Complex> d_integrals(patch_elements);
        double patch_area = 0.0;
        Complex chi_volume(0.0, 0.0);
        Complex chi_boundary(0.0, 0.0);
        for (int local_element = 0; local_element < patch_elements;
             ++local_element) {
            const int element = diagnostic.elements[local_element];
            const Triangle &triangle = mesh.elems[element];
            const auto vertex_it = std::find(
                triangle.begin(), triangle.end(), vertex);
            if (vertex_it == triangle.end())
                throw std::runtime_error("candidate vertex patch is inconsistent");
            const int local_vertex = static_cast<int>(
                std::distance(triangle.begin(), vertex_it));
            Complex integral(0.0, 0.0);
            Complex f_hat(0.0, 0.0);
            Complex mass_hat(0.0, 0.0);
            for (const PhysicalTriangleQuadraturePoint &point :
                 triangle_quadrature_points(
                     mesh, element, config.quadrature,
                     config.quadrature_context)) {
                Complex solution(0.0, 0.0);
                for (int local = 0; local < 3; ++local)
                    solution += point.barycentric[local]
                        * values(triangle[local]);
                const double psi = point.barycentric[local_vertex];
                f_hat += point.weight * psi * source_projection[element];
                mass_hat += point.weight * psi * solution;
            }
            integral = f_hat
                + operators.wavenumber * operators.wavenumber
                    * operators.refractive_index[element] * mass_hat
                - areas[element] * operators.diffusion[element]
                    * geometry[element].gradients[local_vertex]
                        .cast<Complex>().dot(gradients[element]);
            d_integrals[local_element] = integral;
            patch_area += areas[element];
            chi_volume += integral;
        }
        for (int edge_index : patch_edges) {
            const EdgeData &edge_data = topology.edges[edge_index];
            if (edge_data.incidences.size() != 1
                || boundary_tag(mesh, edge_data.nodes) != BoundaryTag::Robin)
                continue;
            chi_boundary += Complex(0.0, operators.wavenumber
                * operators.boundary_beta)
                * edge_hat_solution_integral(
                    mesh, edge_data.nodes, vertex, values);
        }
        diagnostic.compatibility_defect = chi_volume + chi_boundary;
        diagnostic.compatibility_identity_error = std::abs(
            diagnostic.compatibility_defect
            - (std::accumulate(
                    d_integrals.begin(), d_integrals.end(),
                    Complex(0.0, 0.0))
               - physical_flux_integral));
        diagnostic.compatibility_correction = diagnostic.touches_dirichlet
            ? Complex(0.0, 0.0)
            : diagnostic.compatibility_defect / patch_area;
        for (int local_element = 0; local_element < patch_elements;
             ++local_element) {
            const int element = diagnostic.elements[local_element];
            d_integrals[local_element] -=
                diagnostic.compatibility_correction * areas[element];
            result.delta_pg[element] += diagnostic.compatibility_correction;
        }
        if (!diagnostic.touches_dirichlet) {
            diagnostic.corrected_compatibility_error = std::abs(
                std::accumulate(
                    d_integrals.begin(), d_integrals.end(),
                    Complex(0.0, 0.0))
                - physical_flux_integral);
        }

        const int free_count = static_cast<int>(free_edges.size());
        Eigen::MatrixXcd mass = Eigen::MatrixXcd::Zero(free_count, free_count);
        ComplexVector gradient = ComplexVector::Zero(free_count);
        Eigen::MatrixXcd divergence = Eigen::MatrixXcd::Zero(
            patch_elements, free_count);
        ComplexVector divergence_rhs(patch_elements);
        std::vector<Eigen::Vector3cd> fixed_local_flux(
            patch_elements, Eigen::Vector3cd::Zero());
        for (int local_element = 0; local_element < patch_elements;
             ++local_element) {
            const int element = diagnostic.elements[local_element];
            const Triangle &triangle = mesh.elems[element];
            const int local_vertex = static_cast<int>(std::distance(
                triangle.begin(),
                std::find(triangle.begin(), triangle.end(), vertex)));
            Eigen::Vector3cd target_integral = Eigen::Vector3cd::Zero();
            for (const PhysicalTriangleQuadraturePoint &point :
                 triangle_quadrature_points(
                     mesh, element, config.quadrature,
                     config.quadrature_context)) {
                Eigen::Matrix<double, 2, 3> basis;
                for (int edge = 0; edge < 3; ++edge) {
                    const Eigen::Vector3d coefficient =
                        flux_to_coefficients[element].col(edge);
                    basis.col(edge) = Eigen::Vector2d(
                        coefficient(0) + coefficient(2) * point.point.x(),
                        coefficient(1) + coefficient(2) * point.point.y());
                }
                target_integral += point.weight
                    * point.barycentric[local_vertex]
                    * basis.transpose().cast<Complex>()
                    * (operators.diffusion[element] * gradients[element]);
            }
            for (int local_edge = 0; local_edge < 3; ++local_edge) {
                const int edge_index = topology.element_edges[element][local_edge];
                const double sign = topology.element_orientations[element][local_edge];
                const auto free = free_index.find(edge_index);
                if (free == free_index.end()) {
                    fixed_local_flux[local_element](local_edge) = sign
                        * fixed_flux.at(edge_index);
                } else {
                    divergence(local_element, free->second) += sign;
                }
            }
            divergence_rhs(local_element) = d_integrals[local_element]
                - fixed_local_flux[local_element].sum();
            for (int local_i = 0; local_i < 3; ++local_i) {
                const int edge_i = topology.element_edges[element][local_i];
                const auto free_i = free_index.find(edge_i);
                if (free_i == free_index.end()) continue;
                const double sign_i = topology.element_orientations[element][local_i];
                Complex local_gradient = target_integral(local_i);
                for (int fixed_j = 0; fixed_j < 3; ++fixed_j) {
                    local_gradient += rt0_mass[element](local_i, fixed_j)
                        * fixed_local_flux[local_element](fixed_j);
                }
                gradient(free_i->second) += sign_i * local_gradient;
                for (int local_j = 0; local_j < 3; ++local_j) {
                    const int edge_j = topology.element_edges[element][local_j];
                    const auto free_j = free_index.find(edge_j);
                    if (free_j == free_index.end()) continue;
                    const double sign_j =
                        topology.element_orientations[element][local_j];
                    mass(free_i->second, free_j->second) +=
                        sign_i * sign_j * rt0_mass[element](local_i, local_j);
                }
            }
        }

        const int constraint_count = diagnostic.touches_dirichlet
            ? patch_elements : std::max(0, patch_elements - 1);
        ComplexVector free_solution = ComplexVector::Zero(free_count);
        if (free_count > 0) {
            Eigen::MatrixXcd saddle = Eigen::MatrixXcd::Zero(
                free_count + constraint_count,
                free_count + constraint_count);
            saddle.topLeftCorner(free_count, free_count) = mass;
            if (constraint_count > 0) {
                saddle.topRightCorner(free_count, constraint_count) =
                    divergence.topRows(constraint_count).adjoint();
                saddle.bottomLeftCorner(constraint_count, free_count) =
                    divergence.topRows(constraint_count);
            }
            ComplexVector rhs = ComplexVector::Zero(
                free_count + constraint_count);
            rhs.head(free_count) = -gradient;
            if (constraint_count > 0)
                rhs.tail(constraint_count) = divergence_rhs.head(constraint_count);
            Eigen::FullPivLU<Eigen::MatrixXcd> factorization(saddle);
            if (factorization.rank() != saddle.rows())
                throw std::runtime_error("candidate RT0 patch saddle system is singular");
            const ComplexVector solution = factorization.solve(rhs);
            if (!solution.allFinite())
                throw std::runtime_error("candidate RT0 patch solve is non-finite");
            free_solution = solution.head(free_count);
        } else if (constraint_count > 0
                   && divergence_rhs.head(constraint_count).norm() > 1e-10) {
            throw std::runtime_error(
                "candidate RT0 patch has constraints but no free flux");
        }

        for (int local_element = 0; local_element < patch_elements;
             ++local_element) {
            const int element = diagnostic.elements[local_element];
            Eigen::Vector3cd local_flux = fixed_local_flux[local_element];
            for (int local_edge = 0; local_edge < 3; ++local_edge) {
                const int edge_index = topology.element_edges[element][local_edge];
                const auto free = free_index.find(edge_index);
                if (free != free_index.end()) {
                    local_flux(local_edge) =
                        topology.element_orientations[element][local_edge]
                        * free_solution(free->second);
                }
            }
            diagnostic.maximum_divergence_residual = std::max(
                diagnostic.maximum_divergence_residual,
                std::abs(local_flux.sum() - d_integrals[local_element]));
            result.element_rt0_coefficients[element] +=
                flux_to_coefficients[element].cast<Complex>() * local_flux;
        }
        result.maximum_patch_compatibility_error = std::max({
            result.maximum_patch_compatibility_error,
            diagnostic.compatibility_identity_error,
            diagnostic.corrected_compatibility_error});
        result.patches.push_back(std::move(diagnostic));
    }

    result.element_eta_squared.assign(element_count, 0.0);
    for (int element = 0; element < element_count; ++element) {
        const Complex expected_divergence = reaction_projection[element]
            - result.delta_pg[element];
        const Complex actual_divergence =
            2.0 * result.element_rt0_coefficients[element](2);
        result.maximum_element_divergence_residual = std::max(
            result.maximum_element_divergence_residual,
            std::abs(actual_divergence - expected_divergence));
        double flux_error_squared = 0.0;
        double oscillation_squared = 0.0;
        for (const PhysicalTriangleQuadraturePoint &point :
             triangle_quadrature_points(
                 mesh, element, config.quadrature,
                 config.quadrature_context)) {
            Complex solution(0.0, 0.0);
            for (int local = 0; local < 3; ++local)
                solution += point.barycentric[local]
                    * values(mesh.elems[element][local]);
            const Eigen::Vector2cd mismatch = evaluate_candidate_rt0_flux(
                result.element_rt0_coefficients[element], point.point)
                + operators.diffusion[element] * gradients[element];
            flux_error_squared += point.weight * mismatch.squaredNorm();
            const Complex projection_remainder = source(point.point)
                + operators.wavenumber * operators.wavenumber
                    * operators.refractive_index[element] * solution
                - reaction_projection[element];
            oscillation_squared += point.weight
                * std::norm(projection_remainder)
                / (operators.wavenumber * operators.wavenumber);
        }
        const double correction_squared = areas[element]
            * std::norm(result.delta_pg[element])
            / (operators.wavenumber * operators.wavenumber);
        result.element_eta_squared[element] = std::max(
            0.0, flux_error_squared + correction_squared
                + oscillation_squared);
    }

    if (result.global_reconstruction) {
        for (const EdgeData &edge_data : topology.edges) {
            if (edge_data.incidences.size() != 1) continue;
            const BoundaryTag tag = boundary_tag(mesh, edge_data.nodes);
            if (tag == BoundaryTag::Dirichlet) continue;
            if (tag == BoundaryTag::Interior)
                throw std::runtime_error(
                    "candidate flux audit found an unclassified physical boundary edge");
            const EdgeIncidence &side = edge_data.incidences.front();
            const Point2 midpoint = 0.5 * (
                mesh.nodes[edge_data.nodes[0]]
                + mesh.nodes[edge_data.nodes[1]]);
            const Point2 direction =
                mesh.nodes[edge_data.nodes[1]]
                - mesh.nodes[edge_data.nodes[0]];
            const double length = direction.norm();
            Eigen::Vector2d normal(direction.y(), -direction.x());
            normal *= side.orientation / length;
            const Complex actual = length * normal.cast<Complex>().dot(
                evaluate_candidate_rt0_flux(
                    result.element_rt0_coefficients[side.element], midpoint));
            Complex expected(0.0, 0.0);
            if (tag == BoundaryTag::Robin) {
                expected = -Complex(0.0, operators.wavenumber
                    * operators.boundary_beta) * length * 0.5
                    * (values(edge_data.nodes[0]) + values(edge_data.nodes[1]));
            } else if (tag != BoundaryTag::Neumann) {
                throw std::runtime_error(
                    "candidate flux audit found an unknown boundary tag");
            }
            result.maximum_boundary_flux_residual = std::max(
                result.maximum_boundary_flux_residual,
                std::abs(actual - expected));
        }
    }

    const double total_squared = std::accumulate(
        result.element_eta_squared.begin(),
        result.element_eta_squared.end(), 0.0);
    result.eta_eq = std::sqrt(std::max(0.0, total_squared));
    result.marked_elements = doerfler_mark(
        result.element_eta_squared, config.doerfler_theta,
        result.doerfler_relative_error);
    if (config.compute_discrete_residual_audit) {
        result.discrete_residual_dual_norm = discrete_residual_dual_norm(
            mesh, operators, source, values, config);
        result.discrete_residual_audit_performed = true;
    }
    return result;
}

namespace {

struct RT2ElementMap {
    Point2 origin;
    Eigen::Matrix2d jacobian = Eigen::Matrix2d::Zero();
    Eigen::Matrix2d inverse = Eigen::Matrix2d::Zero();
    double determinant = 0.0;
    double area = 0.0;
};

RT2ElementMap rt2_element_map(const TriMesh &mesh, int element) {
    const Triangle &triangle = mesh.elems[element];
    RT2ElementMap map;
    map.origin = mesh.nodes[triangle[0]];
    map.jacobian.col(0) = mesh.nodes[triangle[1]] - map.origin;
    map.jacobian.col(1) = mesh.nodes[triangle[2]] - map.origin;
    map.determinant = map.jacobian.determinant();
    if (std::abs(map.determinant) <= 1e-15)
        throw std::invalid_argument("candidate RT2 mesh contains a degenerate triangle");
    map.inverse = map.jacobian.inverse();
    map.area = 0.5 * std::abs(map.determinant);
    return map;
}

Eigen::Vector2d rt2_reference_point(
    const RT2ElementMap &map, const Point2 &point) {
    return map.inverse * (point - map.origin);
}

Eigen::Matrix<double, 2, 15> rt2_reference_basis(
    const Eigen::Vector2d &reference) {
    const double x = reference.x();
    const double y = reference.y();
    const Eigen::Matrix<double, 6, 1> p =
        (Eigen::Matrix<double, 6, 1>() <<
             1.0, x, y, x * x, x * y, y * y).finished();
    Eigen::Matrix<double, 2, 15> basis =
        Eigen::Matrix<double, 2, 15>::Zero();
    basis.block<1, 6>(0, 0) = p.transpose();
    basis.block<1, 6>(1, 6) = p.transpose();
    basis.col(12) = Eigen::Vector2d(x * x * x, x * x * y);
    basis.col(13) = Eigen::Vector2d(x * x * y, x * y * y);
    basis.col(14) = Eigen::Vector2d(x * y * y, y * y * y);
    return basis;
}

Eigen::Matrix<double, 15, 1> rt2_reference_divergence(
    const Eigen::Vector2d &reference) {
    const double x = reference.x();
    const double y = reference.y();
    Eigen::Matrix<double, 15, 1> divergence =
        Eigen::Matrix<double, 15, 1>::Zero();
    divergence.segment<6>(0) << 0.0, 1.0, 0.0, 2.0 * x, y, 0.0;
    divergence.segment<6>(6) << 0.0, 0.0, 1.0, 0.0, x, 2.0 * y;
    divergence(12) = 4.0 * x * x;
    divergence(13) = 4.0 * x * y;
    divergence(14) = 4.0 * y * y;
    return divergence;
}

Eigen::Matrix<double, 2, 15> rt2_physical_basis(
    const RT2ElementMap &map, const Point2 &point) {
    return map.jacobian
        * rt2_reference_basis(rt2_reference_point(map, point))
        / map.determinant;
}

Eigen::Matrix<double, 15, 1> rt2_physical_divergence_basis(
    const RT2ElementMap &map, const Point2 &point) {
    return rt2_reference_divergence(rt2_reference_point(map, point))
        / map.determinant;
}

Eigen::Vector3d physical_barycentric(
    const RT2ElementMap &map, const Point2 &point) {
    const Eigen::Vector2d reference = rt2_reference_point(map, point);
    return Eigen::Vector3d(
        1.0 - reference.x() - reference.y(),
        reference.x(), reference.y());
}

std::array<Point2, 6> p2_interpolation_points(
    const TriMesh &mesh, int element) {
    const Triangle &triangle = mesh.elems[element];
    const Point2 &p0 = mesh.nodes[triangle[0]];
    const Point2 &p1 = mesh.nodes[triangle[1]];
    const Point2 &p2 = mesh.nodes[triangle[2]];
    return {p0, p1, p2, 0.5 * (p0 + p1),
            0.5 * (p1 + p2), 0.5 * (p2 + p0)};
}

Eigen::Vector3cd project_source_p1(
    const ComplexFunction &source,
    const std::vector<PhysicalTriangleQuadraturePoint> &quadrature,
    double area) {
    Eigen::Vector3cd rhs = Eigen::Vector3cd::Zero();
    for (const PhysicalTriangleQuadraturePoint &point : quadrature) {
        for (int local = 0; local < 3; ++local)
            rhs(local) += point.weight * point.barycentric[local]
                * source(point.point);
    }
    Eigen::Matrix3d mass = Eigen::Matrix3d::Constant(area / 12.0);
    mass.diagonal().array() = area / 6.0;
    return mass.ldlt().solve(rhs);
}

Eigen::Vector2d edge_outward_normal(
    const TriMesh &mesh,
    const EdgeData &edge,
    const EdgeIncidence &side) {
    const Point2 direction =
        mesh.nodes[edge.nodes[1]] - mesh.nodes[edge.nodes[0]];
    const double length = direction.norm();
    return side.orientation
        * Eigen::Vector2d(direction.y(), -direction.x()) / length;
}

constexpr std::array<double, 3> kEdgePoints{
    0.1127016653792583, 0.5, 0.8872983346207417};
constexpr std::array<double, 3> kEdgeWeights{
    5.0 / 18.0, 8.0 / 18.0, 5.0 / 18.0};

Complex interpolate_p1(
    const TriMesh &mesh,
    int element,
    const ComplexVector &values,
    const Eigen::Vector3d &barycentric) {
    Complex result(0.0, 0.0);
    for (int local = 0; local < 3; ++local)
        result += barycentric(local) * values(mesh.elems[element][local]);
    return result;
}

} // namespace

Eigen::Vector2cd evaluate_candidate_rt2_flux(
    const TriMesh &mesh,
    int element,
    const CandidateRT2Coefficients &coefficients,
    const Point2 &point) {
    const RT2ElementMap map = rt2_element_map(mesh, element);
    return rt2_physical_basis(map, point).cast<Complex>() * coefficients;
}

Complex evaluate_candidate_rt2_divergence(
    const TriMesh &mesh,
    int element,
    const CandidateRT2Coefficients &coefficients,
    const Point2 &point) {
    const RT2ElementMap map = rt2_element_map(mesh, element);
    return rt2_physical_divergence_basis(map, point)
        .cast<Complex>().dot(coefficients);
}

CandidateFluxRT2Result reconstruct_candidate_flux_rt2(
    const TriMesh &mesh,
    const HelmholtzOperators &operators,
    const ComplexFunction &source,
    const ComplexVector &values,
    const CandidateFluxConfig &config) {
    validate_boundary_tags(mesh);
    validate_quadrature_policy(config.quadrature);
    const int element_count = static_cast<int>(mesh.elems.size());
    const int node_count = static_cast<int>(mesh.nodes.size());
    if (!source || values.size() != node_count || !values.allFinite()
        || operators.system.rows() != node_count
        || operators.diffusion.size() != mesh.elems.size()
        || operators.refractive_index.size() != mesh.elems.size()
        || !(operators.wavenumber > 0.0)
        || !(config.doerfler_theta > 0.0
             && config.doerfler_theta <= 1.0)) {
        throw std::invalid_argument(
            "candidate RT2 reconstruction inputs are inconsistent");
    }

    const auto prepare_begin = std::chrono::steady_clock::now();
    const MeshTopology topology = build_topology(mesh);
    std::vector<RT2ElementMap> maps(element_count);
    std::vector<ElementGeometry> geometry(element_count);
    std::vector<Eigen::Vector2cd> gradients(element_count);
    std::vector<Eigen::Vector3cd> source_projection(element_count);
    std::vector<std::vector<PhysicalTriangleQuadraturePoint>> quadrature(
        element_count);
    std::vector<std::vector<Eigen::Matrix<double, 2, 15>>> rt2_basis(
        element_count);
    std::vector<Eigen::Matrix<double, 15, 15>> rt2_mass(element_count);
    std::vector<std::array<Point2, 6>> p2_points(element_count);
    std::vector<std::exception_ptr> prepare_errors(element_count);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if(element_count >= 64)
#endif
    for (int element = 0; element < element_count; ++element) {
        try {
        maps[element] = rt2_element_map(mesh, element);
        geometry[element] = element_geometry(mesh, element);
        gradients[element] = element_gradient(
            mesh, element, values, geometry[element]);
        quadrature[element] = triangle_quadrature_points(
            mesh, element, config.quadrature, config.quadrature_context);
        source_projection[element] = project_source_p1(
            source, quadrature[element], maps[element].area);
        rt2_mass[element].setZero();
        rt2_basis[element].reserve(quadrature[element].size());
        for (const PhysicalTriangleQuadraturePoint &point :
             quadrature[element]) {
            rt2_basis[element].push_back(
                rt2_physical_basis(maps[element], point.point));
            const auto &basis = rt2_basis[element].back();
            rt2_mass[element] += point.weight
                * basis.transpose() * basis;
        }
        p2_points[element] = p2_interpolation_points(mesh, element);
        } catch (...) {
            prepare_errors[element] = std::current_exception();
        }
    }
    for (const std::exception_ptr &error : prepare_errors)
        if (error) std::rethrow_exception(error);

    std::vector<std::vector<int>> vertex_elements(node_count);
    for (int element = 0; element < element_count; ++element)
        for (int node : mesh.elems[element])
            vertex_elements[node].push_back(element);
    std::vector<char> active_vertex(node_count, config.active_elements.empty());
    if (!config.active_elements.empty()) {
        for (int element : config.active_elements) {
            if (element < 0 || element >= element_count)
                throw std::out_of_range("active candidate element is out of range");
            for (int node : mesh.elems[element]) active_vertex[node] = true;
        }
    }

    CandidateFluxRT2Result result;
    result.global_reconstruction = config.active_elements.empty();
    result.element_rt2_coefficients.assign(
        element_count, CandidateRT2Coefficients::Zero());
    result.projected_source_p1 = source_projection;
    result.delta_pg.assign(element_count, Complex(0.0, 0.0));
    result.time_prepare = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - prepare_begin).count();
    const auto patch_begin = std::chrono::steady_clock::now();

    struct PatchContribution {
        std::vector<std::pair<int, Complex>> delta_pg;
        std::vector<std::pair<int, CandidateRT2Coefficients>> flux;
        CandidateFluxPatchDiagnostics diagnostic;
        std::exception_ptr error;
        bool active = false;
    };
    std::vector<PatchContribution> patch_contributions(node_count);
#if defined(_OPENMP)
#pragma omp parallel
    {
#pragma omp single
        result.parallel_threads = omp_get_num_threads();
#pragma omp for schedule(dynamic, 8)
#endif
    for (int vertex = 0; vertex < node_count; ++vertex) {
        if (!active_vertex[vertex] || vertex_elements[vertex].empty()) continue;
        PatchContribution &patch_output = patch_contributions[vertex];
        patch_output.active = true;
        try {
        CandidateFluxPatchDiagnostics diagnostic;
        diagnostic.vertex = vertex;
        diagnostic.elements = vertex_elements[vertex];
        const int patch_elements = static_cast<int>(diagnostic.elements.size());
        std::map<int, int> patch_local;
        std::vector<char> in_patch(element_count, false);
        double patch_area = 0.0;
        for (int local = 0; local < patch_elements; ++local) {
            const int element = diagnostic.elements[local];
            patch_local[element] = local;
            in_patch[element] = true;
            patch_area += maps[element].area;
        }

        std::set<int> patch_edges;
        for (int element : diagnostic.elements)
            for (int edge : topology.element_edges[element])
                patch_edges.insert(edge);

        Complex volume_integral(0.0, 0.0);
        Complex boundary_integral(0.0, 0.0);
        std::vector<int> local_vertex_index(patch_elements, -1);
        for (int local = 0; local < patch_elements; ++local) {
            const int element = diagnostic.elements[local];
            const Triangle &triangle = mesh.elems[element];
            const auto found = std::find(
                triangle.begin(), triangle.end(), vertex);
            local_vertex_index[local] = static_cast<int>(
                std::distance(triangle.begin(), found));
            for (const PhysicalTriangleQuadraturePoint &point :
                 quadrature[element]) {
                const Eigen::Vector3d bary(
                    point.barycentric[0], point.barycentric[1],
                    point.barycentric[2]);
                const Complex projected_f = bary.cast<Complex>().dot(
                    source_projection[element]);
                const Complex solution = interpolate_p1(
                    mesh, element, values, bary);
                const double psi = point.barycentric[
                    local_vertex_index[local]];
                const Complex d = psi * (
                    projected_f
                    + operators.wavenumber * operators.wavenumber
                        * operators.refractive_index[element] * solution)
                    - operators.diffusion[element]
                        * geometry[element].gradients[
                            local_vertex_index[local]].cast<Complex>()
                            .dot(gradients[element]);
                volume_integral += point.weight * d;
            }
        }

        for (int edge_index : patch_edges) {
            const EdgeData &edge = topology.edges[edge_index];
            std::vector<const EdgeIncidence *> sides;
            for (const EdgeIncidence &side : edge.incidences)
                if (in_patch[side.element]) sides.push_back(&side);
            if (sides.size() != 1 || edge.incidences.size() != 1) continue;
            const BoundaryTag tag = boundary_tag(mesh, edge.nodes);
            if (tag == BoundaryTag::Dirichlet) {
                diagnostic.touches_dirichlet = true;
                continue;
            }
            if (tag == BoundaryTag::Neumann) continue;
            if (tag != BoundaryTag::Robin)
                throw std::runtime_error(
                    "candidate RT2 flux found an unclassified physical boundary edge");
            const double length = (
                mesh.nodes[edge.nodes[1]] - mesh.nodes[edge.nodes[0]]).norm();
            const int element = sides.front()->element;
            const int local = patch_local.at(element);
            for (int q = 0; q < 3; ++q) {
                const double t = kEdgePoints[q];
                const Point2 point = (1.0 - t) * mesh.nodes[edge.nodes[0]]
                    + t * mesh.nodes[edge.nodes[1]];
                const Eigen::Vector3d bary = physical_barycentric(
                    maps[element], point);
                const double psi = bary(local_vertex_index[local]);
                const Complex solution = interpolate_p1(
                    mesh, element, values, bary);
                boundary_integral += length * kEdgeWeights[q]
                    * (-Complex(0.0, operators.wavenumber
                        * operators.boundary_beta) * psi * solution);
            }
        }

        const Complex chi = volume_integral - boundary_integral;
        const Complex correction = diagnostic.touches_dirichlet
            ? Complex(0.0, 0.0) : chi / patch_area;
        diagnostic.compatibility_defect = chi;
        diagnostic.compatibility_correction = correction;
        diagnostic.compatibility_identity_error = std::abs(
            chi - (volume_integral - boundary_integral));
        diagnostic.corrected_compatibility_error = diagnostic.touches_dirichlet
            ? 0.0 : std::abs(
                volume_integral - correction * patch_area
                - boundary_integral);
        for (int element : diagnostic.elements)
            patch_output.delta_pg.emplace_back(element, correction);

        const int unknowns = 15 * patch_elements;
        Eigen::MatrixXd mass = Eigen::MatrixXd::Zero(unknowns, unknowns);
        ComplexVector gradient = ComplexVector::Zero(unknowns);
        for (int local = 0; local < patch_elements; ++local) {
            const int element = diagnostic.elements[local];
            Eigen::Matrix<Complex, 15, 1> local_gradient =
                Eigen::Matrix<Complex, 15, 1>::Zero();
            for (int q = 0;
                 q < static_cast<int>(quadrature[element].size()); ++q) {
                const PhysicalTriangleQuadraturePoint &point =
                    quadrature[element][q];
                const auto &basis = rt2_basis[element][q];
                const double psi = point.barycentric[
                    local_vertex_index[local]];
                local_gradient += point.weight * psi
                    * basis.transpose().cast<Complex>()
                    * (operators.diffusion[element] * gradients[element]);
            }
            mass.block<15, 15>(15 * local, 15 * local) = rt2_mass[element];
            gradient.segment<15>(15 * local) = local_gradient;
        }

        std::vector<Eigen::VectorXd> constraint_rows;
        std::vector<Complex> constraint_rhs;
        for (int local = 0; local < patch_elements; ++local) {
            const int element = diagnostic.elements[local];
            for (const Point2 &point : p2_points[element]) {
                Eigen::VectorXd row = Eigen::VectorXd::Zero(unknowns);
                row.segment<15>(15 * local) =
                    rt2_physical_divergence_basis(maps[element], point);
                const Eigen::Vector3d bary = physical_barycentric(
                    maps[element], point);
                const Complex projected_f = bary.cast<Complex>().dot(
                    source_projection[element]);
                const Complex solution = interpolate_p1(
                    mesh, element, values, bary);
                const double psi = bary(local_vertex_index[local]);
                const Complex d = psi * (
                    projected_f
                    + operators.wavenumber * operators.wavenumber
                        * operators.refractive_index[element] * solution)
                    - operators.diffusion[element]
                        * geometry[element].gradients[
                            local_vertex_index[local]].cast<Complex>()
                            .dot(gradients[element])
                    - correction;
                constraint_rows.push_back(std::move(row));
                constraint_rhs.push_back(d);
            }
        }

        for (int edge_index : patch_edges) {
            const EdgeData &edge = topology.edges[edge_index];
            std::vector<const EdgeIncidence *> sides;
            for (const EdgeIncidence &side : edge.incidences)
                if (in_patch[side.element]) sides.push_back(&side);
            if (sides.size() == 2) {
                for (double t : kEdgePoints) {
                    const Point2 point = (1.0 - t) * mesh.nodes[edge.nodes[0]]
                        + t * mesh.nodes[edge.nodes[1]];
                    Eigen::VectorXd row = Eigen::VectorXd::Zero(unknowns);
                    for (const EdgeIncidence *side : sides) {
                        const int local = patch_local.at(side->element);
                        const Eigen::Vector2d normal = edge_outward_normal(
                            mesh, edge, *side);
                        row.segment<15>(15 * local) +=
                            rt2_physical_basis(
                                maps[side->element], point).transpose()
                            * normal;
                    }
                    constraint_rows.push_back(std::move(row));
                    constraint_rhs.emplace_back(0.0, 0.0);
                }
                continue;
            }
            if (sides.size() != 1)
                throw std::runtime_error("candidate RT2 patch edge incidence is invalid");
            const EdgeIncidence &side = *sides.front();
            bool free_dirichlet = false;
            BoundaryTag tag = BoundaryTag::Interior;
            if (edge.incidences.size() == 1) {
                tag = boundary_tag(mesh, edge.nodes);
                if (tag == BoundaryTag::Interior)
                    throw std::runtime_error(
                        "candidate RT2 flux found an unclassified physical boundary edge");
                free_dirichlet = tag == BoundaryTag::Dirichlet;
            }
            if (free_dirichlet) continue;
            const int local = patch_local.at(side.element);
            const Eigen::Vector2d normal = edge_outward_normal(mesh, edge, side);
            for (double t : kEdgePoints) {
                const Point2 point = (1.0 - t) * mesh.nodes[edge.nodes[0]]
                    + t * mesh.nodes[edge.nodes[1]];
                Eigen::VectorXd row = Eigen::VectorXd::Zero(unknowns);
                row.segment<15>(15 * local) =
                    rt2_physical_basis(maps[side.element], point).transpose()
                    * normal;
                Complex prescribed(0.0, 0.0);
                if (edge.incidences.size() == 1 && tag == BoundaryTag::Robin) {
                    const Eigen::Vector3d bary = physical_barycentric(
                        maps[side.element], point);
                    prescribed = -Complex(0.0, operators.wavenumber
                        * operators.boundary_beta)
                        * bary(local_vertex_index[local])
                        * interpolate_p1(mesh, side.element, values, bary);
                } else if (edge.incidences.size() == 1
                           && tag != BoundaryTag::Neumann) {
                    throw std::runtime_error(
                        "candidate RT2 flux found an unknown boundary tag");
                }
                constraint_rows.push_back(std::move(row));
                constraint_rhs.push_back(prescribed);
            }
        }

        Eigen::MatrixXd constraints(
            static_cast<int>(constraint_rows.size()), unknowns);
        ComplexVector rhs(static_cast<int>(constraint_rows.size()));
        for (int row = 0; row < constraints.rows(); ++row) {
            constraints.row(row) = constraint_rows[row].transpose();
            rhs(row) = constraint_rhs[row];
        }
        // One rank-revealing factorization supplies both a particular
        // constrained field and the constraint kernel.  The previous path
        // decomposed the same dense patch matrix twice (COD for solve,
        // FullPivLU for kernel), which dominated large RT2 reconstructions.
        Eigen::FullPivLU<Eigen::MatrixXd> constraint_lu(constraints);
        constraint_lu.setThreshold(1e-11);
        ComplexVector particular(unknowns);
        particular.real() = constraint_lu.solve(rhs.real());
        particular.imag() = constraint_lu.solve(rhs.imag());
        const double compatibility_residual = (
            constraints.cast<Complex>() * particular - rhs).norm()
            / std::max(1.0, rhs.norm());
        if (compatibility_residual > 2e-8)
            throw std::runtime_error(
                "candidate RT2 patch constraints are incompatible");

        const Eigen::MatrixXd kernel = constraint_lu.kernel();
        ComplexVector coefficients = particular;
        if (kernel.cols() > 0) {
            const Eigen::MatrixXd reduced =
                kernel.transpose() * mass * kernel;
            Eigen::LDLT<Eigen::MatrixXd> factorization(reduced);
            if (factorization.info() != Eigen::Success)
                throw std::runtime_error("candidate RT2 reduced mass factorization failed");
            const ComplexVector reduced_rhs = -kernel.transpose().cast<Complex>()
                * (mass.cast<Complex>() * particular + gradient);
            coefficients += kernel.cast<Complex>()
                * factorization.solve(reduced_rhs);
        }
        if (!coefficients.allFinite())
            throw std::runtime_error("candidate RT2 patch solve is non-finite");
        const ComplexVector constraint_error =
            constraints.cast<Complex>() * coefficients - rhs;
        diagnostic.maximum_constrained_flux_residual =
            constraint_error.size() == 0 ? 0.0
            : constraint_error.cwiseAbs().maxCoeff();

        for (int local = 0; local < patch_elements; ++local) {
            const int element = diagnostic.elements[local];
            const CandidateRT2Coefficients local_coefficients =
                coefficients.segment<15>(15 * local);
            patch_output.flux.emplace_back(element, local_coefficients);
            if (!config.compute_discrete_residual_audit) continue;
            for (const Point2 &point : p2_points[element]) {
                const Complex actual = evaluate_candidate_rt2_divergence(
                    mesh, element, local_coefficients, point);
                const Eigen::Vector3d bary = physical_barycentric(
                    maps[element], point);
                const Complex projected_f = bary.cast<Complex>().dot(
                    source_projection[element]);
                const Complex solution = interpolate_p1(
                    mesh, element, values, bary);
                const double psi = bary(local_vertex_index[local]);
                const Complex expected = psi * (
                    projected_f
                    + operators.wavenumber * operators.wavenumber
                        * operators.refractive_index[element] * solution)
                    - operators.diffusion[element]
                        * geometry[element].gradients[
                            local_vertex_index[local]].cast<Complex>()
                            .dot(gradients[element])
                    - correction;
                diagnostic.maximum_divergence_residual = std::max(
                    diagnostic.maximum_divergence_residual,
                    std::abs(actual - expected));
            }
        }
        patch_output.diagnostic = std::move(diagnostic);
        } catch (...) {
            patch_output.error = std::current_exception();
        }
    }
#if defined(_OPENMP)
    }
#endif
    result.time_patch_solve = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - patch_begin).count();

    // Merge in vertex order.  This avoids concurrent writes to overlapping
    // elements and keeps floating-point summation reproducible across thread
    // counts and OpenMP schedules.
    const auto merge_begin = std::chrono::steady_clock::now();
    for (int vertex = 0; vertex < node_count; ++vertex) {
        PatchContribution &patch_output = patch_contributions[vertex];
        if (!patch_output.active) continue;
        if (patch_output.error) std::rethrow_exception(patch_output.error);
        for (const auto &[element, correction] : patch_output.delta_pg)
            result.delta_pg[element] += correction;
        for (const auto &[element, coefficients] : patch_output.flux)
            result.element_rt2_coefficients[element] += coefficients;
        const CandidateFluxPatchDiagnostics &diagnostic =
            patch_output.diagnostic;
        result.maximum_patch_compatibility_error = std::max({
            result.maximum_patch_compatibility_error,
            diagnostic.compatibility_identity_error,
            diagnostic.corrected_compatibility_error});
        result.patches.push_back(std::move(patch_output.diagnostic));
    }
    result.time_deterministic_merge = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - merge_begin).count();

    const auto audit_begin = std::chrono::steady_clock::now();
    result.element_eta_squared.assign(element_count, 0.0);
    std::vector<double> element_divergence_residual(element_count, 0.0);
    std::vector<std::exception_ptr> audit_errors(element_count);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if(element_count >= 64)
#endif
    for (int element = 0; element < element_count; ++element) {
        try {
        if (config.compute_discrete_residual_audit) {
        for (const Point2 &point : p2_points[element]) {
            const Eigen::Vector3d bary = physical_barycentric(
                maps[element], point);
            const Complex expected = bary.cast<Complex>().dot(
                source_projection[element])
                + operators.wavenumber * operators.wavenumber
                    * operators.refractive_index[element]
                    * interpolate_p1(mesh, element, values, bary)
                - result.delta_pg[element];
            element_divergence_residual[element] = std::max(
                element_divergence_residual[element],
                std::abs(evaluate_candidate_rt2_divergence(
                    mesh, element, result.element_rt2_coefficients[element],
                point) - expected));
        }
        }
        double flux_error_squared = 0.0;
        double oscillation_squared = 0.0;
        for (const PhysicalTriangleQuadraturePoint &point :
             quadrature[element]) {
            const Eigen::Vector3d bary(
                point.barycentric[0], point.barycentric[1],
                point.barycentric[2]);
            const Complex projected_f = bary.cast<Complex>().dot(
                source_projection[element]);
            const Eigen::Vector2cd mismatch = evaluate_candidate_rt2_flux(
                mesh, element, result.element_rt2_coefficients[element],
                point.point) + operators.diffusion[element] * gradients[element];
            flux_error_squared += point.weight * mismatch.squaredNorm();
            oscillation_squared += point.weight
                * std::norm(source(point.point) - projected_f)
                / (operators.wavenumber * operators.wavenumber);
        }
        const double correction_squared = maps[element].area
            * std::norm(result.delta_pg[element])
            / (operators.wavenumber * operators.wavenumber);
        result.element_eta_squared[element] = std::max(
            0.0, flux_error_squared + correction_squared
                + oscillation_squared);
        } catch (...) {
            audit_errors[element] = std::current_exception();
        }
    }
    for (int element = 0; element < element_count; ++element) {
        if (audit_errors[element])
            std::rethrow_exception(audit_errors[element]);
        result.maximum_element_divergence_residual = std::max(
            result.maximum_element_divergence_residual,
            element_divergence_residual[element]);
    }

    if (result.global_reconstruction
        && config.compute_discrete_residual_audit) {
        for (const EdgeData &edge : topology.edges) {
            for (double t : kEdgePoints) {
                const Point2 point = (1.0 - t) * mesh.nodes[edge.nodes[0]]
                    + t * mesh.nodes[edge.nodes[1]];
                if (edge.incidences.size() == 2) {
                    Complex jump(0.0, 0.0);
                    for (const EdgeIncidence &side : edge.incidences) {
                        jump += edge_outward_normal(mesh, edge, side)
                            .cast<Complex>().dot(evaluate_candidate_rt2_flux(
                                mesh, side.element,
                                result.element_rt2_coefficients[side.element],
                                point));
                    }
                    result.maximum_normal_continuity_residual = std::max(
                        result.maximum_normal_continuity_residual,
                        std::abs(jump));
                    continue;
                }
                const BoundaryTag tag = boundary_tag(mesh, edge.nodes);
                if (tag == BoundaryTag::Dirichlet) continue;
                if (tag == BoundaryTag::Interior)
                    throw std::runtime_error(
                        "candidate RT2 flux audit found an unclassified physical boundary edge");
                const EdgeIncidence &side = edge.incidences.front();
                const Complex actual = edge_outward_normal(mesh, edge, side)
                    .cast<Complex>().dot(evaluate_candidate_rt2_flux(
                        mesh, side.element,
                        result.element_rt2_coefficients[side.element], point));
                Complex expected(0.0, 0.0);
                if (tag == BoundaryTag::Robin) {
                    const Eigen::Vector3d bary = physical_barycentric(
                        maps[side.element], point);
                    expected = -Complex(0.0, operators.wavenumber
                        * operators.boundary_beta)
                        * interpolate_p1(mesh, side.element, values, bary);
                } else if (tag != BoundaryTag::Neumann) {
                    throw std::runtime_error(
                        "candidate RT2 flux audit found an unknown boundary tag");
                }
                result.maximum_boundary_flux_residual = std::max(
                    result.maximum_boundary_flux_residual,
                    std::abs(actual - expected));
            }
        }
    }

    const double total_squared = std::accumulate(
        result.element_eta_squared.begin(),
        result.element_eta_squared.end(), 0.0);
    result.eta_eq = std::sqrt(std::max(0.0, total_squared));
    result.marked_elements = doerfler_mark(
        result.element_eta_squared, config.doerfler_theta,
        result.doerfler_relative_error);
    if (config.compute_discrete_residual_audit) {
        result.discrete_residual_dual_norm = discrete_residual_dual_norm(
            mesh, operators, source, values, config);
        result.discrete_residual_audit_performed = true;
    }
    result.time_estimator_and_audit = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - audit_begin).count();
    return result;
}

} // namespace lod2d::helmholtz::adaptive

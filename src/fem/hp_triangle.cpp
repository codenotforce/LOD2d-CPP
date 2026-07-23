#include "fem/hp_triangle.h"
#include "mesh/refine.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace lod2d {
namespace {

std::uint64_t edge_key(int a, int b) {
    const auto lo = static_cast<std::uint32_t>(std::min(a, b));
    const auto hi = static_cast<std::uint32_t>(std::max(a, b));
    return (static_cast<std::uint64_t>(lo) << 32U) | hi;
}

double integer_power(double value, int exponent) {
    double result = 1.0;
    for (int i = 0; i < exponent; ++i) result *= value;
    return result;
}

std::vector<EdgeQuadraturePoint> gauss_legendre_unit_interval(int order) {
    if (order < 1) throw std::invalid_argument("Gauss order must be positive");
    std::vector<EdgeQuadraturePoint> points(order);
    const int half = (order + 1) / 2;
    constexpr double pi = 3.141592653589793238462643383279502884;
    for (int i = 0; i < half; ++i) {
        double z = std::cos(pi * (i + 0.75) / (order + 0.5));
        double derivative = 0.0;
        for (int iteration = 0; iteration < 100; ++iteration) {
            double previous = 1.0;
            double polynomial = z;
            for (int degree = 2; degree <= order; ++degree) {
                const double next = ((2.0 * degree - 1.0) * z * polynomial
                    - (degree - 1.0) * previous) / degree;
                previous = polynomial;
                polynomial = next;
            }
            if (order == 1) {
                previous = 1.0;
                polynomial = z;
            }
            derivative =
                order * (z * polynomial - previous) / (z * z - 1.0);
            const double next_z = z - polynomial / derivative;
            if (std::abs(next_z - z) < 1e-15) {
                z = next_z;
                break;
            }
            z = next_z;
        }
        const double weight =
            1.0 / ((1.0 - z * z) * derivative * derivative);
        points[i] = {(1.0 - z) * 0.5, weight};
        points[order - 1 - i] = {(1.0 + z) * 0.5, weight};
    }
    return points;
}

} // namespace

std::vector<EdgeQuadraturePoint> edge_gauss_quadrature(int order) {
    return gauss_legendre_unit_interval(order);
}

std::vector<TriangleQuadraturePoint> triangle_gauss_quadrature(int order) {
    const auto line = gauss_legendre_unit_interval(order);
    std::vector<TriangleQuadraturePoint> result;
    result.reserve(line.size() * line.size());
    for (const auto &r : line) {
        for (const auto &s : line) {
            const double x = r.coordinate;
            const double y = (1.0 - r.coordinate) * s.coordinate;
            result.push_back({
                Eigen::Vector3d(1.0 - x - y, x, y),
                r.weight * s.weight * (1.0 - r.coordinate)
            });
        }
    }
    return result;
}

HpTriangleBasis::HpTriangleBasis(int degree) : degree_(degree) {
    if (degree < 1 || degree > 3)
        throw std::invalid_argument(
            "continuous hp triangle currently supports p=1,2,3");

    nodes_.push_back(Eigen::Vector3d(1.0, 0.0, 0.0));
    nodes_.push_back(Eigen::Vector3d(0.0, 1.0, 0.0));
    nodes_.push_back(Eigen::Vector3d(0.0, 0.0, 1.0));
    static constexpr int edges[3][2] = {{0, 1}, {1, 2}, {2, 0}};
    for (const auto &edge : edges) {
        for (int r = 1; r < degree; ++r) {
            Eigen::Vector3d node = Eigen::Vector3d::Zero();
            node(edge[0]) = static_cast<double>(degree - r) / degree;
            node(edge[1]) = static_cast<double>(r) / degree;
            nodes_.push_back(node);
        }
    }
    for (int i = 1; i < degree; ++i) {
        for (int j = 1; j < degree - i; ++j) {
            const int k = degree - i - j;
            nodes_.push_back(Eigen::Vector3d(i, j, k) / degree);
        }
    }

    for (int total = 0; total <= degree; ++total)
        for (int x_power = 0; x_power <= total; ++x_power)
            monomials_.push_back({x_power, total - x_power});
    if (monomials_.size() != nodes_.size())
        throw std::logic_error("hp triangle node and monomial counts differ");

    Eigen::MatrixXd vandermonde(nodes_.size(), monomials_.size());
    for (int row = 0; row < static_cast<int>(nodes_.size()); ++row) {
        const double x = nodes_[row](1);
        const double y = nodes_[row](2);
        for (int col = 0; col < static_cast<int>(monomials_.size()); ++col) {
            vandermonde(row, col) =
                integer_power(x, monomials_[col][0])
                * integer_power(y, monomials_[col][1]);
        }
    }
    inverse_vandermonde_ = vandermonde.inverse();
    if (!inverse_vandermonde_.allFinite())
        throw std::runtime_error("hp triangle Vandermonde inversion failed");
}

HpBasisEvaluation HpTriangleBasis::evaluate(
    const Eigen::Vector3d &barycentric) const {
    const double x = barycentric(1);
    const double y = barycentric(2);
    Eigen::VectorXd monomial(monomials_.size());
    Eigen::VectorXd derivative_x(monomials_.size());
    Eigen::VectorXd derivative_y(monomials_.size());
    for (int i = 0; i < static_cast<int>(monomials_.size()); ++i) {
        const int a = monomials_[i][0];
        const int b = monomials_[i][1];
        monomial(i) = integer_power(x, a) * integer_power(y, b);
        derivative_x(i) = a == 0 ? 0.0
            : a * integer_power(x, a - 1) * integer_power(y, b);
        derivative_y(i) = b == 0 ? 0.0
            : b * integer_power(x, a) * integer_power(y, b - 1);
    }
    HpBasisEvaluation result;
    result.values = inverse_vandermonde_.transpose() * monomial;
    result.reference_gradients.resize(size(), 2);
    result.reference_gradients.col(0) =
        inverse_vandermonde_.transpose() * derivative_x;
    result.reference_gradients.col(1) =
        inverse_vandermonde_.transpose() * derivative_y;
    return result;
}

HpTriSpace::HpTriSpace(const TriMesh &mesh, int degree)
    : mesh_(&mesh), basis_(degree) {
    if (mesh.nodes.empty() || mesh.elems.empty())
        throw std::invalid_argument("hp space requires a nonempty mesh");

    std::tie(edges_, boundary_edges_) = compute_edges(mesh);
    edge_indices_.reserve(edges_.size());
    for (int edge = 0; edge < static_cast<int>(edges_.size()); ++edge)
        edge_indices_.emplace(edge_key(edges_[edge][0], edges_[edge][1]), edge);

    dof_points_ = mesh.nodes;
    for (const Edge &edge : edges_) {
        for (int r = 1; r < degree; ++r) {
            dof_points_.push_back(
                ((degree - r) * mesh.nodes[edge[0]]
                 + r * mesh.nodes[edge[1]]) / degree);
        }
    }
    const int interior_per_element = (degree - 1) * (degree - 2) / 2;
    for (const Triangle &triangle : mesh.elems) {
        for (int i = 1; i < degree; ++i) {
            for (int j = 1; j < degree - i; ++j) {
                const int k = degree - i - j;
                dof_points_.push_back(
                    (i * mesh.nodes[triangle[0]]
                     + j * mesh.nodes[triangle[1]]
                     + k * mesh.nodes[triangle[2]]) / degree);
            }
        }
    }

    element_dofs_.resize(mesh.elems.size());
    const int edge_dof_offset = static_cast<int>(mesh.nodes.size());
    const int interior_offset =
        edge_dof_offset + static_cast<int>(edges_.size()) * (degree - 1);
    static constexpr int local_edges[3][2] = {{0, 1}, {1, 2}, {2, 0}};
    for (int element = 0; element < static_cast<int>(mesh.elems.size());
         ++element) {
        const Triangle &triangle = mesh.elems[element];
        auto &dofs = element_dofs_[element];
        dofs.reserve(basis_.size());
        dofs.insert(dofs.end(), triangle.begin(), triangle.end());
        for (const auto &local_edge : local_edges) {
            const int a = triangle[local_edge[0]];
            const int b = triangle[local_edge[1]];
            const int edge = edge_index(a, b);
            for (int r = 1; r < degree; ++r) {
                const int global_r = a < b ? r : degree - r;
                dofs.push_back(
                    edge_dof_offset + edge * (degree - 1) + global_r - 1);
            }
        }
        for (int local = 0; local < interior_per_element; ++local)
            dofs.push_back(
                interior_offset + element * interior_per_element + local);
        if (dofs.size() != static_cast<std::size_t>(basis_.size()))
            throw std::logic_error("hp element DOF count is inconsistent");
    }

    dof_incidence_.assign(dof_count(), 0);
    for (const auto &dofs : element_dofs_)
        for (int dof : dofs) ++dof_incidence_[dof];
}

int HpTriSpace::edge_index(int a, int b) const {
    const auto it = edge_indices_.find(edge_key(a, b));
    if (it == edge_indices_.end())
        throw std::out_of_range("mesh edge not found");
    return it->second;
}

} // namespace lod2d

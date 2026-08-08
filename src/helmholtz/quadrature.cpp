#include "helmholtz/quadrature.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>

namespace lod2d::helmholtz {
namespace {

struct GaussPoint {
    double coordinate = 0.0;
    double weight = 0.0;
};

struct Subtriangle {
    std::array<Point2, 3> points;
    std::array<std::array<double, 3>, 3> parent_barycentric;
};

std::vector<GaussPoint> compute_gauss_legendre_unit(int order) {
    if (order <= 0) throw std::invalid_argument("Gauss order must be positive");
    std::vector<GaussPoint> result(order);
    const double pi = std::acos(-1.0);
    const int half = (order + 1) / 2;
    for (int index = 0; index < half; ++index) {
        double root = std::cos(pi * (index + 0.75) / (order + 0.5));
        double derivative = 0.0;
        for (int iteration = 0; iteration < 100; ++iteration) {
            double previous = 1.0;
            double current = root;
            for (int degree = 2; degree <= order; ++degree) {
                const double next = ((2.0 * degree - 1.0) * root * current
                    - (degree - 1.0) * previous) / degree;
                previous = current;
                current = next;
            }
            derivative = order * (root * current - previous) / (root * root - 1.0);
            const double update = current / derivative;
            root -= update;
            if (std::abs(update) < 4.0e-16) break;
        }
        const double weight = 1.0 / ((1.0 - root * root) * derivative * derivative);
        result[index] = {0.5 * (1.0 - root), weight};
        result[order - 1 - index] = {0.5 * (1.0 + root), weight};
    }
    return result;
}

const std::vector<GaussPoint> &gauss_legendre_unit(int order) {
    thread_local std::map<int, std::vector<GaussPoint>> cache;
    const auto found = cache.find(order);
    if (found != cache.end()) return found->second;
    return cache.emplace(order, compute_gauss_legendre_unit(order)).first->second;
}

double signed_determinant(const std::array<Point2, 3> &points) {
    const Point2 first = points[1] - points[0];
    const Point2 second = points[2] - points[0];
    return first.x() * second.y() - first.y() * second.x();
}

double max_edge_length(const std::array<Point2, 3> &points) {
    return std::max({
        (points[0] - points[1]).norm(),
        (points[1] - points[2]).norm(),
        (points[2] - points[0]).norm()});
}

double box_distance(const std::array<Point2, 3> &points, const Point2 &point) {
    Point2 lower = points[0];
    Point2 upper = points[0];
    for (int vertex = 1; vertex < 3; ++vertex) {
        lower = lower.cwiseMin(points[vertex]);
        upper = upper.cwiseMax(points[vertex]);
    }
    const double dx = std::max({lower.x() - point.x(), 0.0, point.x() - upper.x()});
    const double dy = std::max({lower.y() - point.y(), 0.0, point.y() - upper.y()});
    return std::hypot(dx, dy);
}

std::array<double, 3> midpoint_barycentric(
    const std::array<double, 3> &first,
    const std::array<double, 3> &second) {
    return {
        0.5 * (first[0] + second[0]),
        0.5 * (first[1] + second[1]),
        0.5 * (first[2] + second[2])};
}

std::array<Subtriangle, 4> split_red(const Subtriangle &triangle) {
    const Point2 p01 = 0.5 * (triangle.points[0] + triangle.points[1]);
    const Point2 p12 = 0.5 * (triangle.points[1] + triangle.points[2]);
    const Point2 p20 = 0.5 * (triangle.points[2] + triangle.points[0]);
    const auto b01 = midpoint_barycentric(
        triangle.parent_barycentric[0], triangle.parent_barycentric[1]);
    const auto b12 = midpoint_barycentric(
        triangle.parent_barycentric[1], triangle.parent_barycentric[2]);
    const auto b20 = midpoint_barycentric(
        triangle.parent_barycentric[2], triangle.parent_barycentric[0]);
    std::array<Subtriangle, 4> children;
    children[0].points = {triangle.points[0], p01, p20};
    children[0].parent_barycentric = {
        triangle.parent_barycentric[0], b01, b20};
    children[1].points = {p01, triangle.points[1], p12};
    children[1].parent_barycentric = {
        b01, triangle.parent_barycentric[1], b12};
    children[2].points = {p20, p12, triangle.points[2]};
    children[2].parent_barycentric = {
        b20, b12, triangle.parent_barycentric[2]};
    children[3].points = {p01, p12, p20};
    children[3].parent_barycentric = {b01, b12, b20};
    return children;
}

void append_duffy_points(
    const Subtriangle &triangle,
    int order,
    int singular_vertex,
    std::vector<PhysicalTriangleQuadraturePoint> &result) {
    const double determinant = std::abs(signed_determinant(triangle.points));
    if (determinant <= 1e-30) throw std::invalid_argument("quadrature triangle is degenerate");
    const auto &gauss = gauss_legendre_unit(order);
    const int first = singular_vertex;
    const int second = (singular_vertex + 1) % 3;
    const int third = (singular_vertex + 2) % 3;
    result.reserve(result.size() + gauss.size() * gauss.size());
    for (const GaussPoint &radial : gauss) {
        for (const GaussPoint &angular : gauss) {
            const std::array<double, 3> local{{
                1.0 - radial.coordinate,
                radial.coordinate * (1.0 - angular.coordinate),
                radial.coordinate * angular.coordinate}};
            const std::array<int, 3> permutation{{first, second, third}};
            PhysicalTriangleQuadraturePoint point;
            for (int local_vertex = 0; local_vertex < 3; ++local_vertex) {
                const int vertex = permutation[local_vertex];
                point.point += local[local_vertex] * triangle.points[vertex];
                for (int parent = 0; parent < 3; ++parent) {
                    point.barycentric[parent] += local[local_vertex]
                        * triangle.parent_barycentric[vertex][parent];
                }
            }
            point.weight = determinant * radial.coordinate
                * radial.weight * angular.weight;
            result.push_back(point);
        }
    }
}

void append_gaussian_points(
    const Subtriangle &triangle,
    const QuadraturePolicy &policy,
    const QuadratureContext &context,
    int depth,
    std::vector<PhysicalTriangleQuadraturePoint> &result) {
    const bool near_feature = box_distance(triangle.points, context.feature_point)
        <= 8.0 * context.feature_scale;
    const bool unresolved = max_edge_length(triangle.points)
        > 2.0 * context.feature_scale;
    if (near_feature && unresolved && depth < policy.max_recursive_subdivisions) {
        for (const Subtriangle &child : split_red(triangle))
            append_gaussian_points(child, policy, context, depth + 1, result);
        return;
    }
    append_duffy_points(triangle, policy.gaussian_triangle_order, 0, result);
}

void append_singular_points(
    const Subtriangle &triangle,
    const QuadraturePolicy &policy,
    const QuadratureContext &context,
    int depth,
    std::vector<PhysicalTriangleQuadraturePoint> &result) {
    const bool near_corner = box_distance(triangle.points, context.feature_point)
        <= 2.0 * context.feature_scale;
    const bool unresolved = max_edge_length(triangle.points)
        > 0.5 * context.feature_scale;
    if (near_corner && unresolved && depth < policy.max_recursive_subdivisions) {
        for (const Subtriangle &child : split_red(triangle))
            append_singular_points(child, policy, context, depth + 1, result);
        return;
    }
    int singular_vertex = 0;
    double distance =
        (triangle.points[0] - context.feature_point).squaredNorm();
    for (int vertex = 1; vertex < 3; ++vertex) {
        const double candidate =
            (triangle.points[vertex] - context.feature_point).squaredNorm();
        if (candidate < distance) {
            distance = candidate;
            singular_vertex = vertex;
        }
    }
    append_duffy_points(
        triangle, policy.singular_triangle_order, singular_vertex, result);
}

Subtriangle initial_triangle(const TriMesh &mesh, int element) {
    if (element < 0 || element >= static_cast<int>(mesh.elems.size()))
        throw std::out_of_range("quadrature element index is out of range");
    const Triangle &indices = mesh.elems[element];
    Subtriangle result;
    result.points = {
        mesh.nodes[indices[0]], mesh.nodes[indices[1]], mesh.nodes[indices[2]]};
    result.parent_barycentric = {{
        {{1.0, 0.0, 0.0}},
        {{0.0, 1.0, 0.0}},
        {{0.0, 0.0, 1.0}}}};
    return result;
}

} // namespace

void validate_quadrature_policy(const QuadraturePolicy &policy) {
    if (policy.base_triangle_order <= 0 || policy.gaussian_triangle_order <= 0
        || policy.singular_triangle_order <= 0
        || policy.max_recursive_subdivisions < 0) {
        throw std::invalid_argument("quadrature policy contains invalid orders");
    }
}

std::vector<PhysicalTriangleQuadraturePoint> triangle_quadrature_points(
    const TriMesh &mesh,
    int element,
    const QuadraturePolicy &policy,
    const QuadratureContext &context) {
    validate_quadrature_policy(policy);
    Subtriangle triangle = initial_triangle(mesh, element);
    std::vector<PhysicalTriangleQuadraturePoint> result;
    if (context.integrand_class == QuadratureClass::LocalizedGaussian) {
        if (!(context.feature_scale > 0.0))
            throw std::invalid_argument("Gaussian quadrature requires a positive feature scale");
        append_gaussian_points(triangle, policy, context, 0, result);
    } else if (context.integrand_class == QuadratureClass::ReentrantSingular) {
        if (!(context.feature_scale > 0.0))
            throw std::invalid_argument("singular quadrature requires a positive feature scale");
        append_singular_points(triangle, policy, context, 0, result);
    } else {
        append_duffy_points(triangle, policy.base_triangle_order, 0, result);
    }
    return result;
}

double integrate_scalar_function(
    const TriMesh &mesh,
    const std::function<double(const Point2 &)> &function,
    const QuadraturePolicy &policy,
    const QuadratureContext &context) {
    if (!function) throw std::invalid_argument("quadrature scalar function is empty");
    double result = 0.0;
    for (int element = 0; element < static_cast<int>(mesh.elems.size()); ++element) {
        for (const auto &point : triangle_quadrature_points(
                 mesh, element, policy, context)) {
            result += point.weight * function(point.point);
        }
    }
    return result;
}

bool operator==(const QuadraturePolicy &lhs, const QuadraturePolicy &rhs) {
    return lhs.base_triangle_order == rhs.base_triangle_order
        && lhs.gaussian_triangle_order == rhs.gaussian_triangle_order
        && lhs.singular_triangle_order == rhs.singular_triangle_order
        && lhs.max_recursive_subdivisions == rhs.max_recursive_subdivisions;
}

} // namespace lod2d::helmholtz

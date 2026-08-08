#pragma once

#include "mesh/types.h"

#include <array>
#include <functional>
#include <vector>

namespace lod2d::helmholtz {

struct QuadraturePolicy {
    int base_triangle_order = 12;
    int gaussian_triangle_order = 20;
    int singular_triangle_order = 24;
    int max_recursive_subdivisions = 8;
};

enum class QuadratureClass {
    Regular,
    LocalizedGaussian,
    ReentrantSingular
};

struct QuadratureContext {
    QuadratureClass integrand_class = QuadratureClass::Regular;
    Point2 feature_point = Point2::Zero();
    double feature_scale = 0.0;
};

struct PhysicalTriangleQuadraturePoint {
    Point2 point = Point2::Zero();
    std::array<double, 3> barycentric{};
    double weight = 0.0;
};

void validate_quadrature_policy(const QuadraturePolicy &policy);

std::vector<PhysicalTriangleQuadraturePoint> triangle_quadrature_points(
    const TriMesh &mesh,
    int element,
    const QuadraturePolicy &policy = {},
    const QuadratureContext &context = {});

double integrate_scalar_function(
    const TriMesh &mesh,
    const std::function<double(const Point2 &)> &function,
    const QuadraturePolicy &policy = {},
    const QuadratureContext &context = {});

bool operator==(const QuadraturePolicy &lhs, const QuadraturePolicy &rhs);

} // namespace lod2d::helmholtz

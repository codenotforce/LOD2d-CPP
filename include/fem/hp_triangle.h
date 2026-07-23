#pragma once

#include "mesh/types.h"

#include <Eigen/Dense>
#include <array>
#include <unordered_map>
#include <vector>

namespace lod2d {

struct TriangleQuadraturePoint {
    Eigen::Vector3d barycentric;
    double weight = 0.0;
};

struct EdgeQuadraturePoint {
    double coordinate = 0.0;
    double weight = 0.0;
};

std::vector<TriangleQuadraturePoint> triangle_gauss_quadrature(int order);
std::vector<EdgeQuadraturePoint> edge_gauss_quadrature(int order);

struct HpBasisEvaluation {
    Eigen::VectorXd values;
    Eigen::MatrixXd reference_gradients;
};

class HpTriangleBasis {
public:
    explicit HpTriangleBasis(int degree);

    int degree() const { return degree_; }
    int size() const { return static_cast<int>(nodes_.size()); }
    const std::vector<Eigen::Vector3d> &nodes() const { return nodes_; }
    HpBasisEvaluation evaluate(const Eigen::Vector3d &barycentric) const;

private:
    int degree_ = 1;
    std::vector<Eigen::Vector3d> nodes_;
    std::vector<std::array<int, 2>> monomials_;
    Eigen::MatrixXd inverse_vandermonde_;
};

class HpTriSpace {
public:
    HpTriSpace(const TriMesh &mesh, int degree);

    const TriMesh &mesh() const { return *mesh_; }
    int degree() const { return basis_.degree(); }
    int dof_count() const { return static_cast<int>(dof_points_.size()); }
    int local_dof_count() const { return basis_.size(); }
    const HpTriangleBasis &basis() const { return basis_; }
    const std::vector<Point2> &dof_points() const { return dof_points_; }
    const std::vector<std::vector<int>> &element_dofs() const {
        return element_dofs_;
    }
    const std::vector<int> &dof_incidence() const { return dof_incidence_; }
    const std::vector<Edge> &edges() const { return edges_; }
    const std::vector<bool> &boundary_edges() const { return boundary_edges_; }
    int edge_index(int a, int b) const;

private:
    const TriMesh *mesh_ = nullptr;
    HpTriangleBasis basis_;
    std::vector<Edge> edges_;
    std::vector<bool> boundary_edges_;
    std::unordered_map<std::uint64_t, int> edge_indices_;
    std::vector<Point2> dof_points_;
    std::vector<std::vector<int>> element_dofs_;
    std::vector<int> dof_incidence_;
};

} // namespace lod2d

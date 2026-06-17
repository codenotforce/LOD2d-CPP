#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <array>
#include <cstdint>
#include <vector>

namespace lod2d {

/// Triangle element: 3 vertex indices (0‑based)
using Triangle = std::array<int, 3>;

/// Edge: sorted vertex pair (i < j)
using Edge = std::array<int, 2>;

/// 2D point
using Point2 = Eigen::Vector2d;

/// Triangulation
struct TriMesh {
    std::vector<Point2>     nodes;       // vertex coordinates
    std::vector<Triangle>   elems;       // element connectivity
    std::vector<int>        dirichlet;   // Dirichlet boundary node indices
};

/// Refinement output: refined mesh + three prolongation matrices
struct RefineOutput {
    TriMesh mesh;
    // Prolongation:  fine = P_node * coarse  (nodal)
    //               fine = P_elem * coarse  (element‑wise)
    //               fine = P_dg   * coarse  (DG dof)
    // Stored as CSC sparse (Eigen default for SparseMatrix<double>)
    Eigen::SparseMatrix<double> P_node;
    Eigen::SparseMatrix<double> P_elem;
    Eigen::SparseMatrix<double> P_dg;
};

} // namespace lod2d

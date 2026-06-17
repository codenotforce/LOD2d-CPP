#pragma once
#include <Eigen/Sparse>
namespace lod2d {
Eigen::VectorXd solve_coarse(const Eigen::SparseMatrix<double>&, const Eigen::VectorXd&, const std::vector<int>&);
}

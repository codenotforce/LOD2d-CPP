#pragma once
#include "mesh/types.h"
#include <Eigen/Sparse>
namespace lod2d {
Eigen::SparseMatrix<double> compute_corrector(int, const Eigen::SparseMatrix<double>&, const TriMesh&, int, const std::vector<int>&, const Eigen::SparseMatrix<double>&, const TriMesh&, int, const std::vector<int>&, const Eigen::SparseMatrix<double>&, const Eigen::SparseMatrix<double>&, const Eigen::SparseMatrix<double>&, const Eigen::SparseMatrix<double>&, const Eigen::SparseMatrix<double>&, int);
}

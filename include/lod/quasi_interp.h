#pragma once
#include "mesh/types.h"
#include <Eigen/Sparse>
namespace lod2d {
Eigen::SparseMatrix<double> build_quasi_interp(const TriMesh&, const TriMesh&, const Eigen::SparseMatrix<double>&, const Eigen::SparseMatrix<double>&, int, int);
}

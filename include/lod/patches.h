#pragma once
#include "mesh/types.h"
#include <Eigen/Sparse>
namespace lod2d {
Eigen::SparseMatrix<double> build_patches(const TriMesh&, int);
}

#pragma once
#include "mesh/types.h"
#include <Eigen/Sparse>

namespace lod2d {

Eigen::SparseMatrix<double> assemble_dg(const TriMesh &mesh,
                                         const std::vector<double> &coeff);

} // namespace lod2d

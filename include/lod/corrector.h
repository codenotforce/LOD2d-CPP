#pragma once
#include "mesh/types.h"
#include <Eigen/Sparse>
#include <vector>

namespace lod2d {

Eigen::SparseMatrix<double>
compute_corrector(int k,
    const Eigen::SparseMatrix<double> &patch,
    const TriMesh &coarse, int NH, const std::vector<int> &nngH,
    const Eigen::SparseMatrix<double> &P0,
    const TriMesh &fine,   int Nh, const std::vector<int> &nngh,
    const std::vector<std::array<int,3>> &dghidx,
    const Eigen::SparseMatrix<double> &cg2dgh,
    const Eigen::SparseMatrix<double> &Shdg,
    const Eigen::SparseMatrix<double> &P1dg,
    const std::vector<std::array<int,3>> &dgHidx,
    const Eigen::SparseMatrix<double> &IH,
    int d);

} // namespace lod2d

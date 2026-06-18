#pragma once
#include "fem/assemble_dg.h"
#include "mesh/types.h"
#include <Eigen/Sparse>
#include <vector>

namespace lod2d {

enum class CorrectorSolver {
    EigenLLT,
    Cholmod
};

using FineElementChildren = std::vector<std::vector<int>>;

FineElementChildren build_fine_element_children(
    const Eigen::SparseMatrix<double> &P0,
    int coarse_element_count);

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
    int d,
    CorrectorSolver solver = CorrectorSolver::EigenLLT,
    const ElementStiffnessBlocks *element_stiffness = nullptr,
    const FineElementChildren *fine_element_children = nullptr);

} // namespace lod2d

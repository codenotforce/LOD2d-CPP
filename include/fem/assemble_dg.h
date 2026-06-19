#pragma once
#include "mesh/types.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

namespace lod2d {

using ElementStiffnessBlocks = std::vector<Eigen::Matrix3d>;

ElementStiffnessBlocks assemble_element_stiffness(const TriMesh &mesh,
                                                   const std::vector<double> &coeff);

Eigen::SparseMatrix<double> assemble_dg_from_element_stiffness(
    const ElementStiffnessBlocks &element_stiffness);

Eigen::SparseMatrix<double> assemble_dg(const TriMesh &mesh,
                                         const std::vector<double> &coeff);

Eigen::SparseMatrix<double> assemble_cg_from_element_stiffness(
    const TriMesh &mesh,
    const ElementStiffnessBlocks &element_stiffness);

Eigen::SparseMatrix<double> assemble_cg_mass(
    const TriMesh &mesh,
    const std::vector<double> &areas);

} // namespace lod2d

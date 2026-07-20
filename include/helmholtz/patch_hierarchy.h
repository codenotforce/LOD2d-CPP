#pragma once

#include "mesh/types.h"

#include <Eigen/Sparse>
#include <vector>

namespace lod2d::helmholtz {

class HelmholtzDirichletPatchHierarchyBuilder {
public:
    HelmholtzDirichletPatchHierarchyBuilder(
        const std::vector<TriMesh> &meshes,
        const std::vector<Eigen::SparseMatrix<double>> &node_prolongations,
        const std::vector<Eigen::SparseMatrix<double>> &element_prolongations);

    std::vector<Eigen::SparseMatrix<double>> build(
        const std::vector<int> &coarse_patch_elements,
        const std::vector<int> &finest_patch_dofs) const;

    bool available() const { return meshes_.size() > 1; }

private:
    const std::vector<TriMesh> &meshes_;
    const std::vector<Eigen::SparseMatrix<double>> &node_prolongations_;
    const std::vector<Eigen::SparseMatrix<double>> &element_prolongations_;
    std::vector<std::vector<int>> node_incidence_;
};

} // namespace lod2d::helmholtz

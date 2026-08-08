#pragma once

#include "helmholtz/operators.h"

#include <Eigen/Sparse>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace lod2d::helmholtz {

struct HelmholtzPatchSystem {
    int target_element = -1;
    std::vector<int> local_vertices;
    std::vector<int> patch_elements;
    Eigen::SparseMatrix<double> stiffness;
    Eigen::SparseMatrix<double> mass;
    Eigen::SparseMatrix<double> robin;
    ComplexSparseMatrix helmholtz;
    Eigen::MatrixXd constraints;
    ComplexMatrix rhs;
    double wavenumber = 0.0;
    double diameter = 0.0;
    std::vector<Eigen::SparseMatrix<double>> geometric_prolongations;
    bool touches_physical_boundary = false;
};

class HelmholtzPatchAssembler {
public:
    HelmholtzPatchAssembler(
        const TriMesh &coarse,
        const TriMesh &fine,
        const Eigen::SparseMatrix<double> &fine_element_prolongation,
        const Eigen::SparseMatrix<double> &fine_dg_prolongation,
        const Eigen::SparseMatrix<double> &quasi_interpolation,
        const Eigen::SparseMatrix<double> &patches,
        const std::vector<TriMesh> &hierarchy_meshes,
        const std::vector<Eigen::SparseMatrix<double>> &node_level_prolongations,
        const std::vector<Eigen::SparseMatrix<double>> &element_level_prolongations,
        const HelmholtzOperators &operators);

    HelmholtzPatchSystem assemble(int target) const;
    std::size_t patch_cost(int target) const;
    int patch_count() const { return static_cast<int>(coarse_.elems.size()); }

private:
    const TriMesh &coarse_;
    const TriMesh &fine_;
    const Eigen::SparseMatrix<double> &fine_dg_prolongation_;
    const Eigen::SparseMatrix<double> &quasi_interpolation_;
    const Eigen::SparseMatrix<double> &patches_;
    const HelmholtzOperators &operators_;
    std::vector<std::vector<int>> children_;
    std::vector<int> fine_incidence_;
    std::vector<char> fine_dirichlet_;
    std::unordered_map<std::uint64_t, int> fine_edge_counts_;
    const std::vector<TriMesh> &hierarchy_meshes_;
    const std::vector<Eigen::SparseMatrix<double>> &node_level_prolongations_;
    const std::vector<Eigen::SparseMatrix<double>> &element_level_prolongations_;
    std::vector<std::vector<int>> hierarchy_incidence_;
};
} // namespace lod2d::helmholtz

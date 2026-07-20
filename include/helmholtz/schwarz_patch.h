#pragma once

#include "helmholtz/operators.h"
#include "helmholtz/patch_hierarchy.h"

#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

namespace lod2d::helmholtz {

enum class HelmholtzSchwarzArtificialBoundary {
    HomogeneousDirichlet,
    Impedance
};

struct HelmholtzSchwarzLocalSystem {
    std::vector<int> global_dofs;
    std::vector<int> core_global_dofs;
    ComplexSparseMatrix matrix;
    Eigen::SparseMatrix<double> mass;
    std::vector<Eigen::SparseMatrix<double>> geometric_prolongations;
    int artificial_boundary_edges = 0;
    int physical_boundary_edges = 0;
};

class HelmholtzSchwarzPatchAssembler {
public:
    HelmholtzSchwarzPatchAssembler(
        const TriMesh &fine_mesh,
        const Eigen::SparseMatrix<double> &fine_element_prolongation,
        const Eigen::SparseMatrix<double> &element_patches,
        int coarse_element_count,
        const HelmholtzOperators &operators);
    HelmholtzSchwarzPatchAssembler(
        const TriMesh &fine_mesh,
        const Eigen::SparseMatrix<double> &fine_element_prolongation,
        const Eigen::SparseMatrix<double> &element_patches,
        int coarse_element_count,
        const std::vector<TriMesh> &hierarchy_meshes,
        const std::vector<Eigen::SparseMatrix<double>> &node_prolongations,
        const std::vector<Eigen::SparseMatrix<double>> &element_prolongations,
        const HelmholtzOperators &operators);

    HelmholtzSchwarzLocalSystem assemble(
        int target,
        HelmholtzSchwarzArtificialBoundary boundary,
        double artificial_impedance_beta = 1.0,
        bool assemble_mass = false,
        bool assemble_hierarchy = false) const;

    int patch_count() const { return static_cast<int>(children_.size()); }

private:
    const TriMesh &fine_mesh_;
    const Eigen::SparseMatrix<double> &element_patches_;
    const HelmholtzOperators &operators_;
    std::vector<std::vector<int>> children_;
    std::vector<int> fine_incidence_;
    std::unordered_map<std::uint64_t, int> global_edge_counts_;
    std::unique_ptr<HelmholtzDirichletPatchHierarchyBuilder>
        hierarchy_builder_;
};

} // namespace lod2d::helmholtz

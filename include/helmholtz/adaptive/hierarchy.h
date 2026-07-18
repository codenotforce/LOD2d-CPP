#pragma once

#include "mesh/refine.h"

#include <cstdint>
#include <vector>

namespace lod2d::helmholtz::adaptive {

struct NestedFineMesh {
    RefineOutput refinement;
    std::vector<int> element_levels;
};

class AdaptiveMeshHierarchy {
public:
    AdaptiveMeshHierarchy(
        const TriMesh &initial_mesh,
        int initial_coarse_level,
        int fine_level);

    const TriMesh &initial_mesh() const { return initial_mesh_; }
    const TriMesh &coarse_mesh() const { return coarse_mesh_; }
    const std::vector<int> &coarse_levels() const { return coarse_levels_; }
    const std::vector<std::uint64_t> &coarse_element_ids() const { return coarse_element_ids_; }
    const std::vector<std::uint64_t> &coarse_parent_ids() const { return coarse_parent_ids_; }
    int fine_level() const { return fine_level_; }

    NestedFineMesh build_nested_fine_mesh() const;
    void refine(const std::vector<int> &marked_elements);

private:
    TriMesh initial_mesh_;
    TriMesh coarse_mesh_;
    std::vector<int> coarse_levels_;
    std::vector<std::uint64_t> coarse_element_ids_;
    std::vector<std::uint64_t> coarse_parent_ids_;
    int fine_level_ = 0;
    std::uint64_t next_element_id_ = 1;
};

NestedFineMesh complete_to_fine_level(
    const TriMesh &coarse_mesh,
    const std::vector<int> &coarse_levels,
    int fine_level);

std::vector<int> refinement_child_levels(
    const TriMesh &parent_mesh,
    const std::vector<int> &parent_levels,
    const RefineOutput &refinement);

std::vector<int> fine_element_parents(
    const Eigen::SparseMatrix<double> &fine_element_prolongation,
    int fine_element_count,
    int coarse_element_count);

} // namespace lod2d::helmholtz::adaptive

#pragma once

#include "mesh/types.h"

namespace lod2d {

/// Single-level conforming longest-edge bisection.
/// Each marked triangle contributes its longest edge to the split set. Every
/// triangle incident to a split edge is then subdivided consistently, so no
/// hanging nodes are introduced.
/// Returns refined mesh + per-level prolongation matrices.
RefineOutput refine(const TriMesh &coarse);

/// Local conforming longest-edge bisection.
/// The input marks coarse element indices. The closure is edge-based: if a
/// marked longest edge is shared, the neighboring triangle is also subdivided
/// across the same edge. Empty marks return identity prolongations.
RefineOutput refine_marked(const TriMesh &coarse, const std::vector<int> &marked_elements);

/// Single-level conforming newest-vertex bisection.
/// The first local vertex of every triangle is the newest vertex, hence the
/// reference edge is the local edge (1,2) in zero-based local numbering.
RefineOutput refine_nvb(const TriMesh &coarse);

/// Local conforming newest-vertex bisection with closure.
/// Marked elements request bisection of their reference edge. The closure adds
/// neighboring reference edges recursively until no hanging nodes remain.
RefineOutput bisect_newest_vertex(const TriMesh &coarse, const std::vector<int> &marked_elements);

/// Multi-level newest-vertex bisection.
RefineOutput refine_mesh_nvb(const TriMesh &initial, int nref);

/// Multi-level refinement: applies refine() nref times.
/// Returns final mesh + accumulated prolongation matrices.
RefineOutput refine_mesh(const TriMesh &initial, int nref);

/// Enumerate unique edges.  Returns edge list and a bool vector
/// where true = boundary edge (appears in exactly 1 triangle).
std::pair<std::vector<Edge>, std::vector<bool>>
compute_edges(const TriMesh &mesh);

/// Compute signed area of each triangle.
std::vector<double> compute_area(const TriMesh &mesh);

} // namespace lod2d

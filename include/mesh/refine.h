#pragma once

#include "mesh/types.h"

namespace lod2d {

/// Single‑level uniform red refinement.
/// Each triangle is split into 4 congruent sub‑triangles
/// by connecting edge midpoints.
/// Returns refined mesh + per‑level prolongation matrices.
RefineOutput refine(const TriMesh &coarse);

/// Multi‑level refinement: applies refine() nref times.
/// Returns final mesh + accumulated prolongation matrices.
RefineOutput refine_mesh(const TriMesh &initial, int nref);

/// Enumerate unique edges.  Returns edge list and a bool vector
/// where true = boundary edge (appears in exactly 1 triangle).
std::pair<std::vector<Edge>, std::vector<bool>>
compute_edges(const TriMesh &mesh);

/// Compute signed area of each triangle.
std::vector<double> compute_area(const TriMesh &mesh);

} // namespace lod2d

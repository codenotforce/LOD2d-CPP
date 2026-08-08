#pragma once

#include "mesh/types.h"

#include <vector>

namespace lod2d::helmholtz {

Edge canonical_edge(int first, int second);

// Validate an explicit edge classification. Empty boundary_edges is accepted
// for legacy meshes and interpreted from geometry plus mesh.dirichlet.
void validate_boundary_tags(const TriMesh &mesh);

// Interior is returned for non-boundary edges. Legacy meshes classify a
// boundary edge as Dirichlet iff both endpoints occur in mesh.dirichlet;
// every other boundary edge is Robin.
BoundaryTag boundary_tag(const TriMesh &mesh, const Edge &edge);

std::vector<Edge> boundary_edges_with_tag(
    const TriMesh &mesh,
    BoundaryTag tag);

std::vector<int> dirichlet_nodes(const TriMesh &mesh);

double boundary_measure(const TriMesh &mesh, BoundaryTag tag);

// Create an explicit complete classification for a legacy mesh.
void tag_all_boundary_edges(TriMesh &mesh, BoundaryTag tag);

// Recompute the legacy node list from explicit Dirichlet edge tags.
void synchronize_dirichlet_nodes(TriMesh &mesh);

// Transfer an explicit coarse classification to a conforming refinement.
void propagate_boundary_tags(const TriMesh &coarse, TriMesh &fine);

} // namespace lod2d::helmholtz

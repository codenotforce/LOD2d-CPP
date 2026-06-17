#include "mesh/refine.h"
#include <algorithm>
#include <map>
#include <set>

namespace lod2d {

std::pair<std::vector<Edge>, std::vector<bool>>
compute_edges(const TriMesh &mesh) {
    // Collect all edges: sort each edge (i<j), count occurrences
    std::map<Edge, int> edge_count;
    for (const auto &tri : mesh.elems) {
        Edge e0{std::min(tri[0], tri[1]), std::max(tri[0], tri[1])};
        Edge e1{std::min(tri[1], tri[2]), std::max(tri[1], tri[2])};
        Edge e2{std::min(tri[2], tri[0]), std::max(tri[2], tri[0])};
        ++edge_count[e0];  ++edge_count[e1];  ++edge_count[e2];
    }

    std::vector<Edge> edges;
    std::vector<bool> is_boundary;
    edges.reserve(edge_count.size());
    is_boundary.reserve(edge_count.size());

    for (const auto &[e, cnt] : edge_count) {
        edges.push_back(e);
        is_boundary.push_back(cnt == 1);  // appears in exactly 1 triangle
    }
    return {edges, is_boundary};
}

std::vector<double> compute_area(const TriMesh &mesh) {
    std::vector<double> areas;
    areas.reserve(mesh.elems.size());
    for (const auto &tri : mesh.elems) {
        const auto &a = mesh.nodes[tri[0]];
        const auto &b = mesh.nodes[tri[1]];
        const auto &c = mesh.nodes[tri[2]];
        double area = 0.5 * ((b.x()-a.x())*(c.y()-a.y()) - (b.y()-a.y())*(c.x()-a.x()));
        areas.push_back(std::abs(area));  // signed → absolute
    }
    return areas;
}

} // namespace lod2d

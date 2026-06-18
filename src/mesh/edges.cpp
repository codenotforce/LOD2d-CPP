#include "mesh/refine.h"
#include <algorithm>
#include <vector>
#include <map>

namespace lod2d {

std::pair<std::vector<Edge>, std::vector<bool>>
compute_edges(const TriMesh &mesh) {
    // Count edge occurrences (same as MATLAB's sparse edge counter)
    std::map<Edge, int> edge_count;
    for (const auto &tri : mesh.elems) {
        Edge e0{std::min(tri[0], tri[1]), std::max(tri[0], tri[1])};
        Edge e1{std::min(tri[1], tri[2]), std::max(tri[1], tri[2])};
        Edge e2{std::min(tri[2], tri[0]), std::max(tri[2], tri[0])};
        ++edge_count[e0];  ++edge_count[e1];  ++edge_count[e2];
    }

    // MATLAB ordering: find(sparse(row, col, 1)) returns entries
    // in COLUMN-MAJOR order.  Edge = {i,j} with i<j means i=row, j=col.
    // Column-major: sort by j (second element) first, then i.
    std::vector<std::pair<Edge, int>> sorted(edge_count.begin(), edge_count.end());
    std::sort(sorted.begin(), sorted.end(),
        [](const auto &a, const auto &b) {
            if (a.first[1] != b.first[1]) return a.first[1] < b.first[1];  // col first
            return a.first[0] < b.first[0];  // then row
        });

    std::vector<Edge> edges;
    std::vector<bool> is_boundary;
    edges.reserve(sorted.size());
    is_boundary.reserve(sorted.size());

    for (const auto &[e, cnt] : sorted) {
        edges.push_back(e);
        is_boundary.push_back(cnt == 1);
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
        areas.push_back(std::abs(area));
    }
    return areas;
}

} // namespace lod2d

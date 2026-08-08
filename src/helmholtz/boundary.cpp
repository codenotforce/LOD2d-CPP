#include "helmholtz/boundary.h"

#include "mesh/refine.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>

namespace lod2d::helmholtz {
namespace {

std::map<Edge, bool> geometric_boundary_map(const TriMesh &mesh) {
    const auto [edges, boundary] = compute_edges(mesh);
    std::map<Edge, bool> result;
    for (std::size_t i = 0; i < edges.size(); ++i)
        result.emplace(edges[i], boundary[i]);
    return result;
}

bool point_on_segment(
    const Point2 &point,
    const Point2 &first,
    const Point2 &second) {
    const Point2 direction = second - first;
    const double length = direction.norm();
    if (length <= 1e-15) return false;
    const Point2 relative = point - first;
    const double cross = direction.x() * relative.y()
                       - direction.y() * relative.x();
    const double scale = std::max(1.0, length);
    if (std::abs(cross) > 1e-11 * scale * scale) return false;
    const double projection = relative.dot(direction);
    return projection >= -1e-11 * scale * scale
        && projection <= direction.squaredNorm() + 1e-11 * scale * scale;
}

} // namespace

Edge canonical_edge(int first, int second) {
    return first < second ? Edge{first, second} : Edge{second, first};
}

void validate_boundary_tags(const TriMesh &mesh) {
    if (mesh.boundary_edges.empty()) return;
    const auto geometric = geometric_boundary_map(mesh);
    std::map<Edge, BoundaryTag> explicit_tags;
    for (const BoundaryEdge &entry : mesh.boundary_edges) {
        const Edge edge = canonical_edge(entry.nodes[0], entry.nodes[1]);
        if (edge[0] < 0 || edge[1] >= static_cast<int>(mesh.nodes.size())
            || edge[0] == edge[1]) {
            throw std::invalid_argument("boundary edge contains an invalid node index");
        }
        const auto found = geometric.find(edge);
        if (found == geometric.end() || !found->second)
            throw std::invalid_argument("tagged edge is not a geometric boundary edge");
        if (entry.tag == BoundaryTag::Interior)
            throw std::invalid_argument("physical boundary edge cannot be tagged Interior");
        if (!explicit_tags.emplace(edge, entry.tag).second)
            throw std::invalid_argument("physical boundary edge is tagged more than once");
    }
    for (const auto &[edge, is_boundary] : geometric) {
        if (is_boundary && !explicit_tags.contains(edge))
            throw std::invalid_argument("explicit boundary classification is incomplete");
    }
}

BoundaryTag boundary_tag(const TriMesh &mesh, const Edge &input) {
    const Edge edge = canonical_edge(input[0], input[1]);
    if (!mesh.boundary_edges.empty()) {
        for (const BoundaryEdge &entry : mesh.boundary_edges) {
            if (canonical_edge(entry.nodes[0], entry.nodes[1]) == edge)
                return entry.tag;
        }
        return BoundaryTag::Interior;
    }
    const auto geometric = geometric_boundary_map(mesh);
    const auto found = geometric.find(edge);
    if (found == geometric.end() || !found->second) return BoundaryTag::Interior;
    const std::set<int> dirichlet(mesh.dirichlet.begin(), mesh.dirichlet.end());
    return dirichlet.contains(edge[0]) && dirichlet.contains(edge[1])
        ? BoundaryTag::Dirichlet
        : BoundaryTag::Robin;
}

std::vector<Edge> boundary_edges_with_tag(const TriMesh &mesh, BoundaryTag tag) {
    const auto [edges, boundary] = compute_edges(mesh);
    std::vector<Edge> result;
    for (std::size_t i = 0; i < edges.size(); ++i) {
        if (boundary[i] && boundary_tag(mesh, edges[i]) == tag)
            result.push_back(edges[i]);
    }
    return result;
}

std::vector<int> dirichlet_nodes(const TriMesh &mesh) {
    if (mesh.boundary_edges.empty()) {
        std::vector<int> result = mesh.dirichlet;
        std::sort(result.begin(), result.end());
        result.erase(std::unique(result.begin(), result.end()), result.end());
        return result;
    }
    std::set<int> result;
    for (const BoundaryEdge &entry : mesh.boundary_edges) {
        if (entry.tag == BoundaryTag::Dirichlet) {
            result.insert(entry.nodes[0]);
            result.insert(entry.nodes[1]);
        }
    }
    return {result.begin(), result.end()};
}

double boundary_measure(const TriMesh &mesh, BoundaryTag tag) {
    double result = 0.0;
    for (const Edge &edge : boundary_edges_with_tag(mesh, tag))
        result += (mesh.nodes[edge[0]] - mesh.nodes[edge[1]]).norm();
    return result;
}

void synchronize_dirichlet_nodes(TriMesh &mesh) {
    if (mesh.boundary_edges.empty()) return;
    std::set<int> nodes;
    for (const BoundaryEdge &entry : mesh.boundary_edges) {
        if (entry.tag != BoundaryTag::Dirichlet) continue;
        nodes.insert(entry.nodes[0]);
        nodes.insert(entry.nodes[1]);
    }
    mesh.dirichlet.assign(nodes.begin(), nodes.end());
}

void tag_all_boundary_edges(TriMesh &mesh, BoundaryTag tag) {
    if (tag == BoundaryTag::Interior)
        throw std::invalid_argument("physical boundary cannot be tagged Interior");
    const auto [edges, boundary] = compute_edges(mesh);
    mesh.boundary_edges.clear();
    for (std::size_t i = 0; i < edges.size(); ++i) {
        if (boundary[i]) mesh.boundary_edges.push_back({edges[i], tag});
    }
    synchronize_dirichlet_nodes(mesh);
    validate_boundary_tags(mesh);
}

void propagate_boundary_tags(const TriMesh &coarse, TriMesh &fine) {
    if (coarse.boundary_edges.empty()) return;
    validate_boundary_tags(coarse);
    const auto [fine_edges, fine_boundary] = compute_edges(fine);
    fine.boundary_edges.clear();
    for (std::size_t index = 0; index < fine_edges.size(); ++index) {
        if (!fine_boundary[index]) continue;
        const Edge fine_edge = fine_edges[index];
        const Point2 &first = fine.nodes[fine_edge[0]];
        const Point2 &second = fine.nodes[fine_edge[1]];
        bool classified = false;
        for (const BoundaryEdge &coarse_entry : coarse.boundary_edges) {
            const Point2 &coarse_first = coarse.nodes[coarse_entry.nodes[0]];
            const Point2 &coarse_second = coarse.nodes[coarse_entry.nodes[1]];
            if (point_on_segment(first, coarse_first, coarse_second)
                && point_on_segment(second, coarse_first, coarse_second)) {
                fine.boundary_edges.push_back({fine_edge, coarse_entry.tag});
                classified = true;
                break;
            }
        }
        if (!classified)
            throw std::runtime_error("refined physical boundary edge has no coarse parent tag");
    }
    synchronize_dirichlet_nodes(fine);
    validate_boundary_tags(fine);
}

} // namespace lod2d::helmholtz

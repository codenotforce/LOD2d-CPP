#include "mesh/refine.h"
#include <Eigen/Sparse>
#include <algorithm>
#include <array>
#include <cmath>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace lod2d {

namespace {

struct EdgeHash {
    size_t operator()(const Edge &e) const {
        const auto a = static_cast<size_t>(e[0]);
        const auto b = static_cast<size_t>(e[1]);
        return a * 1000003ULL ^ (b + 0x9e3779b97f4a7c15ULL + (a << 6) + (a >> 2));
    }
};

Edge make_edge(int a, int b) {
    if (a < b) return Edge{a, b};
    return Edge{b, a};
}

std::array<Edge, 3> triangle_edges(const Triangle &tri) {
    return {make_edge(tri[0], tri[1]), make_edge(tri[1], tri[2]), make_edge(tri[2], tri[0])};
}

bool edge_less_matlab(const Edge &a, const Edge &b) {
    if (a[1] != b[1]) return a[1] < b[1];
    return a[0] < b[0];
}

double edge_length2(const TriMesh &mesh, const Edge &edge) {
    return (mesh.nodes[edge[0]] - mesh.nodes[edge[1]]).squaredNorm();
}

Edge longest_edge(const TriMesh &mesh, const Triangle &tri) {
    const auto edges = triangle_edges(tri);
    Edge best = edges[0];
    double best_len = edge_length2(mesh, best);
    for (int i = 1; i < 3; ++i) {
        const double len = edge_length2(mesh, edges[i]);
        if (len > best_len + 1e-14 || (std::abs(len - best_len) <= 1e-14 && edge_less_matlab(edges[i], best))) {
            best = edges[i];
            best_len = len;
        }
    }
    return best;
}

struct VertexRef {
    int node = -1;
    std::array<double, 3> w{0.0, 0.0, 0.0};
};

VertexRef original_vertex(int local, int node) {
    VertexRef v;
    v.node = node;
    v.w[local] = 1.0;
    return v;
}

VertexRef midpoint_vertex(int local_a, int local_b, int node) {
    VertexRef v;
    v.node = node;
    v.w[local_a] = 0.5;
    v.w[local_b] = 0.5;
    return v;
}

double signed_area(const TriMesh &mesh, const std::array<VertexRef, 3> &tri) {
    const auto &a = mesh.nodes[tri[0].node];
    const auto &b = mesh.nodes[tri[1].node];
    const auto &c = mesh.nodes[tri[2].node];
    return 0.5 * ((b.x() - a.x()) * (c.y() - a.y()) - (b.y() - a.y()) * (c.x() - a.x()));
}

void add_child(const TriMesh &mesh,
               int parent,
               std::array<VertexRef, 3> child,
               std::vector<Triangle> &fine_elems,
               std::vector<Eigen::Triplet<double>> &pe_triplets,
               std::vector<Eigen::Triplet<double>> &pdg_triplets) {
    if (signed_area(mesh, child) < 0.0) std::swap(child[1], child[2]);

    const int row = static_cast<int>(fine_elems.size());
    fine_elems.push_back({child[0].node, child[1].node, child[2].node});
    pe_triplets.emplace_back(row, parent, 1.0);

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (child[i].w[j] != 0.0) {
                pdg_triplets.emplace_back(3 * row + i, 3 * parent + j, child[i].w[j]);
            }
        }
    }
}

RefineOutput identity_refinement(const TriMesh &mesh) {
    const int np = static_cast<int>(mesh.nodes.size());
    const int nt = static_cast<int>(mesh.elems.size());
    Eigen::SparseMatrix<double> P_node(np, np);
    P_node.setIdentity();
    Eigen::SparseMatrix<double> P_elem(nt, nt);
    P_elem.setIdentity();
    Eigen::SparseMatrix<double> P_dg(3 * nt, 3 * nt);
    P_dg.setIdentity();
    return {mesh, std::move(P_node), std::move(P_elem), std::move(P_dg)};
}

} // namespace

RefineOutput refine_marked(const TriMesh &coarse, const std::vector<int> &marked_elements) {
    const int np = static_cast<int>(coarse.nodes.size());
    const int nt = static_cast<int>(coarse.elems.size());
    if (marked_elements.empty()) return identity_refinement(coarse);

    std::unordered_set<Edge, EdgeHash> split_edge_set;
    split_edge_set.reserve(marked_elements.size() * 2 + 8);
    for (int elem : marked_elements) {
        if (elem < 0 || elem >= nt) throw std::out_of_range("marked element index out of range");
        split_edge_set.insert(longest_edge(coarse, coarse.elems[elem]));
    }
    if (split_edge_set.empty()) return identity_refinement(coarse);

    std::vector<Edge> split_edges(split_edge_set.begin(), split_edge_set.end());
    std::sort(split_edges.begin(), split_edges.end(), edge_less_matlab);

    std::unordered_map<Edge, int, EdgeHash> edge_midpoint;
    edge_midpoint.reserve(split_edges.size() * 2 + 8);
    for (size_t i = 0; i < split_edges.size(); ++i) {
        edge_midpoint.emplace(split_edges[i], np + static_cast<int>(i));
    }

    std::unordered_map<Edge, int, EdgeHash> edge_count;
    edge_count.reserve(static_cast<size_t>(3 * nt));
    for (const auto &tri : coarse.elems) {
        for (const auto &edge : triangle_edges(tri)) ++edge_count[edge];
    }

    TriMesh fine;
    fine.nodes = coarse.nodes;
    fine.nodes.reserve(np + static_cast<int>(split_edges.size()));
    for (const auto &edge : split_edges) {
        fine.nodes.push_back(0.5 * (coarse.nodes[edge[0]] + coarse.nodes[edge[1]]));
    }

    std::vector<Eigen::Triplet<double>> pn_triplets;
    pn_triplets.reserve(np + 2 * split_edges.size());
    for (int i = 0; i < np; ++i) pn_triplets.emplace_back(i, i, 1.0);
    for (size_t i = 0; i < split_edges.size(); ++i) {
        const int row = np + static_cast<int>(i);
        pn_triplets.emplace_back(row, split_edges[i][0], 0.5);
        pn_triplets.emplace_back(row, split_edges[i][1], 0.5);
    }
    Eigen::SparseMatrix<double> P_node(np + static_cast<int>(split_edges.size()), np);
    P_node.setFromTriplets(pn_triplets.begin(), pn_triplets.end());

    fine.elems.reserve(static_cast<size_t>(nt) + split_edges.size() * 2);
    std::vector<Eigen::Triplet<double>> pe_triplets;
    pe_triplets.reserve(static_cast<size_t>(nt) + split_edges.size() * 2);
    std::vector<Eigen::Triplet<double>> pdg_triplets;
    pdg_triplets.reserve(6 * static_cast<size_t>(nt) + 18 * split_edges.size());

    for (int t = 0; t < nt; ++t) {
        const Triangle &tri = coarse.elems[t];
        const Edge e01 = make_edge(tri[0], tri[1]);
        const Edge e12 = make_edge(tri[1], tri[2]);
        const Edge e20 = make_edge(tri[2], tri[0]);
        const bool s01 = edge_midpoint.count(e01) != 0;
        const bool s12 = edge_midpoint.count(e12) != 0;
        const bool s20 = edge_midpoint.count(e20) != 0;

        const VertexRef v0 = original_vertex(0, tri[0]);
        const VertexRef v1 = original_vertex(1, tri[1]);
        const VertexRef v2 = original_vertex(2, tri[2]);
        const VertexRef m01 = s01 ? midpoint_vertex(0, 1, edge_midpoint.at(e01)) : VertexRef{};
        const VertexRef m12 = s12 ? midpoint_vertex(1, 2, edge_midpoint.at(e12)) : VertexRef{};
        const VertexRef m20 = s20 ? midpoint_vertex(2, 0, edge_midpoint.at(e20)) : VertexRef{};

        const int nsplit = static_cast<int>(s01) + static_cast<int>(s12) + static_cast<int>(s20);
        if (nsplit == 0) {
            add_child(fine, t, {v0, v1, v2}, fine.elems, pe_triplets, pdg_triplets);
        } else if (nsplit == 1) {
            if (s01) {
                add_child(fine, t, {v0, m01, v2}, fine.elems, pe_triplets, pdg_triplets);
                add_child(fine, t, {m01, v1, v2}, fine.elems, pe_triplets, pdg_triplets);
            } else if (s12) {
                add_child(fine, t, {v1, m12, v0}, fine.elems, pe_triplets, pdg_triplets);
                add_child(fine, t, {m12, v2, v0}, fine.elems, pe_triplets, pdg_triplets);
            } else {
                add_child(fine, t, {v2, m20, v1}, fine.elems, pe_triplets, pdg_triplets);
                add_child(fine, t, {m20, v0, v1}, fine.elems, pe_triplets, pdg_triplets);
            }
        } else if (nsplit == 2) {
            if (s01 && s12) {
                add_child(fine, t, {m01, v1, m12}, fine.elems, pe_triplets, pdg_triplets);
                add_child(fine, t, {v0, m01, m12}, fine.elems, pe_triplets, pdg_triplets);
                add_child(fine, t, {v0, m12, v2}, fine.elems, pe_triplets, pdg_triplets);
            } else if (s12 && s20) {
                add_child(fine, t, {m12, v2, m20}, fine.elems, pe_triplets, pdg_triplets);
                add_child(fine, t, {v1, m12, m20}, fine.elems, pe_triplets, pdg_triplets);
                add_child(fine, t, {v1, m20, v0}, fine.elems, pe_triplets, pdg_triplets);
            } else {
                add_child(fine, t, {m20, v0, m01}, fine.elems, pe_triplets, pdg_triplets);
                add_child(fine, t, {v2, m20, m01}, fine.elems, pe_triplets, pdg_triplets);
                add_child(fine, t, {v2, m01, v1}, fine.elems, pe_triplets, pdg_triplets);
            }
        } else {
            add_child(fine, t, {v0, m01, m20}, fine.elems, pe_triplets, pdg_triplets);
            add_child(fine, t, {m01, v1, m12}, fine.elems, pe_triplets, pdg_triplets);
            add_child(fine, t, {m20, m12, v2}, fine.elems, pe_triplets, pdg_triplets);
            add_child(fine, t, {m01, m12, m20}, fine.elems, pe_triplets, pdg_triplets);
        }
    }

    std::set<int> dir_set(coarse.dirichlet.begin(), coarse.dirichlet.end());
    for (const auto &edge : split_edges) {
        const bool is_boundary = edge_count[edge] == 1;
        if (is_boundary && dir_set.count(edge[0]) && dir_set.count(edge[1])) {
            dir_set.insert(edge_midpoint.at(edge));
        }
    }
    fine.dirichlet.assign(dir_set.begin(), dir_set.end());

    Eigen::SparseMatrix<double> P_elem(static_cast<int>(fine.elems.size()), nt);
    P_elem.setFromTriplets(pe_triplets.begin(), pe_triplets.end());
    Eigen::SparseMatrix<double> P_dg(3 * static_cast<int>(fine.elems.size()), 3 * nt);
    P_dg.setFromTriplets(pdg_triplets.begin(), pdg_triplets.end());

    return {std::move(fine), std::move(P_node), std::move(P_elem), std::move(P_dg)};
}


bool edge_less_lex(const Edge &a, const Edge &b) {
    if (a[0] != b[0]) return a[0] < b[0];
    return a[1] < b[1];
}

Edge reference_edge(const Triangle &tri) {
    return make_edge(tri[1], tri[2]);
}

std::vector<Edge> sorted_edges_from_set(const std::unordered_set<Edge, EdgeHash> &edges) {
    std::vector<Edge> out(edges.begin(), edges.end());
    std::sort(out.begin(), out.end(), edge_less_lex);
    return out;
}

std::unordered_set<Edge, EdgeHash> edge_set_from_vector(const std::vector<Edge> &edges) {
    std::unordered_set<Edge, EdgeHash> out;
    out.reserve(edges.size() * 2 + 8);
    for (const auto &edge : edges) out.insert(edge);
    return out;
}

std::vector<Edge> newest_vertex_closure(const TriMesh &mesh, std::vector<Edge> marked_edges) {
    std::unordered_set<Edge, EdgeHash> marked = edge_set_from_vector(marked_edges);
    bool changed = true;
    while (changed) {
        changed = false;
        const size_t before = marked.size();
        for (const auto &tri : mesh.elems) {
            bool touched = false;
            for (const auto &edge : triangle_edges(tri)) {
                if (marked.count(edge) != 0) {
                    touched = true;
                    break;
                }
            }
            if (touched) marked.insert(reference_edge(tri));
        }
        changed = marked.size() != before;
    }
    return sorted_edges_from_set(marked);
}

RefineOutput bisect_reference_edges(const TriMesh &coarse, const std::vector<char> &cut_elem) {
    const int np = static_cast<int>(coarse.nodes.size());
    const int nt = static_cast<int>(coarse.elems.size());
    if (static_cast<int>(cut_elem.size()) != nt) throw std::invalid_argument("cut_elem size mismatch");

    std::vector<int> keep_idx;
    std::vector<int> cut_idx;
    keep_idx.reserve(nt);
    cut_idx.reserve(nt);
    for (int t = 0; t < nt; ++t) {
        if (cut_elem[t]) cut_idx.push_back(t);
        else keep_idx.push_back(t);
    }
    if (cut_idx.empty()) return identity_refinement(coarse);

    std::unordered_map<Edge, int, EdgeHash> edge_count;
    edge_count.reserve(static_cast<size_t>(3 * nt));
    for (const auto &tri : coarse.elems) {
        for (const auto &edge : triangle_edges(tri)) ++edge_count[edge];
    }

    std::unordered_set<Edge, EdgeHash> new_edge_set;
    new_edge_set.reserve(cut_idx.size() * 2 + 8);
    for (int t : cut_idx) new_edge_set.insert(reference_edge(coarse.elems[t]));
    std::vector<Edge> new_edges = sorted_edges_from_set(new_edge_set);

    std::unordered_map<Edge, int, EdgeHash> edge_midpoint;
    edge_midpoint.reserve(new_edges.size() * 2 + 8);
    for (size_t i = 0; i < new_edges.size(); ++i) {
        edge_midpoint.emplace(new_edges[i], np + static_cast<int>(i));
    }

    TriMesh fine;
    fine.nodes = coarse.nodes;
    fine.nodes.reserve(np + static_cast<int>(new_edges.size()));
    for (const auto &edge : new_edges) {
        fine.nodes.push_back(0.5 * (coarse.nodes[edge[0]] + coarse.nodes[edge[1]]));
    }

    std::vector<Eigen::Triplet<double>> pn_triplets;
    pn_triplets.reserve(np + 2 * new_edges.size());
    for (int i = 0; i < np; ++i) pn_triplets.emplace_back(i, i, 1.0);
    for (size_t i = 0; i < new_edges.size(); ++i) {
        const int row = np + static_cast<int>(i);
        pn_triplets.emplace_back(row, new_edges[i][0], 0.5);
        pn_triplets.emplace_back(row, new_edges[i][1], 0.5);
    }
    Eigen::SparseMatrix<double> P_node(np + static_cast<int>(new_edges.size()), np);
    P_node.setFromTriplets(pn_triplets.begin(), pn_triplets.end());

    const int n_keep = static_cast<int>(keep_idx.size());
    const int n_cut = static_cast<int>(cut_idx.size());
    fine.elems.resize(n_keep + 2 * n_cut);
    std::vector<int> parent(fine.elems.size(), -1);

    std::vector<Eigen::Triplet<double>> pe_triplets;
    pe_triplets.reserve(fine.elems.size());
    std::vector<Eigen::Triplet<double>> pdg_triplets;
    pdg_triplets.reserve(3 * n_keep + 10 * n_cut);

    auto add_pdg_block = [&](int row_elem, int parent_elem, const double block[3][3]) {
        pe_triplets.emplace_back(row_elem, parent_elem, 1.0);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                if (block[i][j] != 0.0) {
                    pdg_triplets.emplace_back(3 * row_elem + i, 3 * parent_elem + j, block[i][j]);
                }
            }
        }
    };

    const double identity[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    const double child1[3][3] = {{0, 0.5, 0.5}, {1, 0, 0}, {0, 1, 0}};
    const double child2[3][3] = {{0, 0.5, 0.5}, {0, 0, 1}, {1, 0, 0}};

    for (int i = 0; i < n_keep; ++i) {
        const int t = keep_idx[i];
        fine.elems[i] = coarse.elems[t];
        parent[i] = t;
        add_pdg_block(i, t, identity);
    }

    for (int i = 0; i < n_cut; ++i) {
        const int t = cut_idx[i];
        const Triangle &tri = coarse.elems[t];
        const int a = tri[0];
        const int b = tri[1];
        const int c = tri[2];
        const int m = edge_midpoint.at(reference_edge(tri));
        const int row1 = n_keep + i;
        const int row2 = n_keep + n_cut + i;
        fine.elems[row1] = {m, a, b};
        fine.elems[row2] = {m, c, a};
        parent[row1] = t;
        parent[row2] = t;
        add_pdg_block(row1, t, child1);
        add_pdg_block(row2, t, child2);
    }

    std::set<int> dir_set(coarse.dirichlet.begin(), coarse.dirichlet.end());
    for (const auto &edge : new_edges) {
        const bool is_boundary = edge_count[edge] == 1;
        if (is_boundary && dir_set.count(edge[0]) && dir_set.count(edge[1])) {
            dir_set.insert(edge_midpoint.at(edge));
        }
    }
    fine.dirichlet.assign(dir_set.begin(), dir_set.end());

    Eigen::SparseMatrix<double> P_elem(static_cast<int>(fine.elems.size()), nt);
    P_elem.setFromTriplets(pe_triplets.begin(), pe_triplets.end());
    Eigen::SparseMatrix<double> P_dg(3 * static_cast<int>(fine.elems.size()), 3 * nt);
    P_dg.setFromTriplets(pdg_triplets.begin(), pdg_triplets.end());

    return {std::move(fine), std::move(P_node), std::move(P_elem), std::move(P_dg)};
}

RefineOutput bisect_newest_vertex(const TriMesh &coarse, const std::vector<int> &marked_elements) {
    const int nt = static_cast<int>(coarse.elems.size());
    if (marked_elements.empty()) return identity_refinement(coarse);

    std::vector<Edge> marked_edges;
    marked_edges.reserve(marked_elements.size());
    for (int elem : marked_elements) {
        if (elem < 0 || elem >= nt) throw std::out_of_range("marked element index out of range");
        marked_edges.push_back(reference_edge(coarse.elems[elem]));
    }

    RefineOutput total = identity_refinement(coarse);
    TriMesh current = coarse;

    while (!marked_edges.empty()) {
        marked_edges = newest_vertex_closure(current, marked_edges);
        const auto marked_set = edge_set_from_vector(marked_edges);

        std::vector<char> cut_elem(current.elems.size(), false);
        std::unordered_set<Edge, EdgeHash> cut_edges;
        cut_edges.reserve(marked_edges.size() * 2 + 8);
        for (int t = 0; t < static_cast<int>(current.elems.size()); ++t) {
            const Edge ref = reference_edge(current.elems[t]);
            if (marked_set.count(ref) != 0) {
                cut_elem[t] = true;
                cut_edges.insert(ref);
            }
        }
        if (cut_edges.empty()) break;

        RefineOutput step = bisect_reference_edges(current, cut_elem);
        total.P_node = step.P_node * total.P_node;
        total.P_elem = step.P_elem * total.P_elem;
        total.P_dg = step.P_dg * total.P_dg;
        current = std::move(step.mesh);
        total.mesh = current;

        std::unordered_set<Edge, EdgeHash> current_edges;
        current_edges.reserve(3 * current.elems.size());
        for (const auto &tri : current.elems) {
            for (const auto &edge : triangle_edges(tri)) current_edges.insert(edge);
        }
        std::vector<Edge> remaining;
        remaining.reserve(marked_edges.size());
        for (const auto &edge : marked_edges) {
            if (current_edges.count(edge) != 0) remaining.push_back(edge);
        }
        marked_edges = std::move(remaining);
    }

    return total;
}

RefineOutput refine_nvb(const TriMesh &coarse) {
    std::vector<int> marked(coarse.elems.size());
    for (int i = 0; i < static_cast<int>(marked.size()); ++i) marked[i] = i;
    return bisect_newest_vertex(coarse, marked);
}

RefineOutput refine_mesh_nvb(const TriMesh &initial, int nref) {
    if (nref < 0) throw std::invalid_argument("nref must be nonnegative");
    int np0 = static_cast<int>(initial.nodes.size());
    int nt0 = static_cast<int>(initial.elems.size());

    Eigen::SparseMatrix<double> P_node(np0, np0);
    P_node.setIdentity();
    Eigen::SparseMatrix<double> P_elem(nt0, nt0);
    P_elem.setIdentity();
    Eigen::SparseMatrix<double> P_dg(3 * nt0, 3 * nt0);
    P_dg.setIdentity();

    TriMesh current = initial;
    for (int lev = 0; lev < nref; ++lev) {
        RefineOutput out = refine_nvb(current);
        P_node = out.P_node * P_node;
        P_elem = out.P_elem * P_elem;
        P_dg = out.P_dg * P_dg;
        current = std::move(out.mesh);
    }

    return {std::move(current), std::move(P_node), std::move(P_elem), std::move(P_dg)};
}

RefineOutput refine(const TriMesh &coarse) {
    std::vector<int> marked(coarse.elems.size());
    for (int i = 0; i < static_cast<int>(marked.size()); ++i) marked[i] = i;
    return refine_marked(coarse, marked);
}

RefineOutput refine_mesh(const TriMesh &initial, int nref) {
    if (nref < 0) throw std::invalid_argument("nref must be nonnegative");
    int np0 = static_cast<int>(initial.nodes.size());
    int nt0 = static_cast<int>(initial.elems.size());

    Eigen::SparseMatrix<double> P_node(np0, np0);
    P_node.setIdentity();
    Eigen::SparseMatrix<double> P_elem(nt0, nt0);
    P_elem.setIdentity();
    Eigen::SparseMatrix<double> P_dg(3 * nt0, 3 * nt0);
    P_dg.setIdentity();

    TriMesh current = initial;

    for (int lev = 0; lev < nref; ++lev) {
        RefineOutput out = refine(current);
        P_node = out.P_node * P_node;
        P_elem = out.P_elem * P_elem;
        P_dg = out.P_dg * P_dg;
        current = std::move(out.mesh);
    }

    return {std::move(current), std::move(P_node), std::move(P_elem), std::move(P_dg)};
}

} // namespace lod2d

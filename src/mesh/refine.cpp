#include "mesh/refine.h"
#include <Eigen/Sparse>
#include <unordered_map>
#include <set>
#include <stdexcept>

namespace lod2d {

// ---- Internal helpers ----

/// Hash for Edge (sorted pair)
struct EdgeHash {
    size_t operator()(const Edge &e) const {
        return std::hash<int>{}(e[0]) ^ (std::hash<int>{}(e[1]) << 1);
    }
};

/// Edge → midpoint index map (0‑based lookup within current refinement level)
static std::unordered_map<Edge, int, EdgeHash>
build_edge_map(const std::vector<Edge> &edges) {
    std::unordered_map<Edge, int, EdgeHash> emap;
    emap.reserve(edges.size());
    for (size_t i = 0; i < edges.size(); ++i)
        emap[edges[i]] = static_cast<int>(i);
    return emap;
}

// ---- Single-level refinement ----

RefineOutput refine(const TriMesh &coarse) {
    int np = static_cast<int>(coarse.nodes.size());
    int nt = static_cast<int>(coarse.elems.size());

    // 1. Enumerate edges
    auto [edges, is_bdy] = compute_edges(coarse);
    int ne = static_cast<int>(edges.size());
    auto edge_map = build_edge_map(edges);

    // 2. New nodes at edge midpoints
    TriMesh fine;
    fine.nodes = coarse.nodes;
    fine.nodes.reserve(np + ne);
    for (const auto &e : edges) {
        fine.nodes.push_back(0.5 * (coarse.nodes[e[0]] + coarse.nodes[e[1]]));
    }

    // 3. Build nodal prolongation P_node  (fine_nodes × coarse_nodes)
    //    Identity for old nodes, 0.5 for midpoints
    std::vector<Eigen::Triplet<double>> pn_triplets;
    pn_triplets.reserve(np + 2 * ne);
    for (int i = 0; i < np; ++i)
        pn_triplets.emplace_back(i, i, 1.0);
    for (int e = 0; e < ne; ++e) {
        pn_triplets.emplace_back(np + e, edges[e][0], 0.5);
        pn_triplets.emplace_back(np + e, edges[e][1], 0.5);
    }
    Eigen::SparseMatrix<double> P_node(np + ne, np);
    P_node.setFromTriplets(pn_triplets.begin(), pn_triplets.end());

    // 4. Split triangles — MATLAB ordering: all sub‑1, then all sub‑2, …
    //    Pre‑compute vertex lookups for all coarse elements at once
    std::vector<int> a_(nt), b_(nt), c_(nt), m12_(nt), m23_(nt), m31_(nt);
    for (int t = 0; t < nt; ++t) {
        int a = coarse.elems[t][0], b = coarse.elems[t][1], c = coarse.elems[t][2];
        a_[t] = a;  b_[t] = b;  c_[t] = c;
        m12_[t] = np + edge_map.at(Edge{std::min(a,b), std::max(a,b)});
        m23_[t] = np + edge_map.at(Edge{std::min(b,c), std::max(b,c)});
        m31_[t] = np + edge_map.at(Edge{std::min(c,a), std::max(c,a)});
    }

    fine.elems.resize(4 * nt);
    int idx = 0;
    // sub 1: (v1, m12, m31) of all elements
    for (int t = 0; t < nt; ++t) fine.elems[idx++] = {a_[t],   m12_[t], m31_[t]};
    // sub 2: (m12, v2, m23)
    for (int t = 0; t < nt; ++t) fine.elems[idx++] = {m12_[t], b_[t],   m23_[t]};
    // sub 3: (m31, m23, v3)
    for (int t = 0; t < nt; ++t) fine.elems[idx++] = {m31_[t], m23_[t], c_[t]  };
    // sub 4: (m12, m23, m31)
    for (int t = 0; t < nt; ++t) fine.elems[idx++] = {m12_[t], m23_[t], m31_[t]};

    // 5. Dirichlet boundary: old Dirichlet + midpoints of boundary edges
    //    where both endpoints are Dirichlet
    std::set<int> dir_set(coarse.dirichlet.begin(), coarse.dirichlet.end());
    for (int e = 0; e < ne; ++e) {
        if (is_bdy[e] && dir_set.count(edges[e][0]) && dir_set.count(edges[e][1]))
            dir_set.insert(np + e);
    }
    fine.dirichlet.assign(dir_set.begin(), dir_set.end());

    // 6. Element prolongation P_elem  (4*nt × nt)
    //    Fine elems ordered by sub‑type: sub‑1(0..nt-1), sub‑2(nt..2nt-1), …
    std::vector<Eigen::Triplet<double>> pe_triplets;
    pe_triplets.reserve(4 * nt);
    for (int t = 0; t < nt; ++t) {
        pe_triplets.emplace_back(t,          t, 1.0);   // sub 1
        pe_triplets.emplace_back(nt + t,     t, 1.0);   // sub 2
        pe_triplets.emplace_back(2*nt + t,   t, 1.0);   // sub 3
        pe_triplets.emplace_back(3*nt + t,   t, 1.0);   // sub 4
    }
    Eigen::SparseMatrix<double> P_elem(4*nt, nt);
    P_elem.setFromTriplets(pe_triplets.begin(), pe_triplets.end());

    // 7. DG dof prolongation P_dg  (12*nt × 3*nt)
    //    MATLAB: [kron(E,sub1); kron(E,sub2); kron(E,sub3); kron(E,sub4)]
    //    Rows grouped by sub‑type, not by element.
    const double sub1[3][3] = {{1,0,0},{.5,.5,0},{.5,0,.5}};
    const double sub2[3][3] = {{.5,.5,0},{0,1,0},{0,.5,.5}};
    const double sub3[3][3] = {{.5,0,.5},{0,.5,.5},{0,0,1}};
    const double sub4[3][3] = {{.5,.5,0},{0,.5,.5},{.5,0,.5}};

    std::vector<Eigen::Triplet<double>> pdg_triplets;

    auto add_sub = [&](int row_off, const double sub[3][3]) {
        for (int t = 0; t < nt; ++t) {
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    if (sub[i][j] != 0.0) {
                        pdg_triplets.emplace_back(
                            row_off + 3*t + i,
                            3*t + j,
                            sub[i][j]);
                    }
                }
            }
        }
    };

    add_sub(0,      sub1);   // rows 0      .. 3*nt-1
    add_sub(3*nt,   sub2);   // rows 3*nt   .. 6*nt-1
    add_sub(6*nt,   sub3);   // rows 6*nt   .. 9*nt-1
    add_sub(9*nt,   sub4);   // rows 9*nt   ..12*nt-1

    Eigen::SparseMatrix<double> P_dg(12*nt, 3*nt);
    P_dg.setFromTriplets(pdg_triplets.begin(), pdg_triplets.end());

    return {std::move(fine), std::move(P_node), std::move(P_elem), std::move(P_dg)};
}

// ---- Multi‑level refinement ----

RefineOutput refine_mesh(const TriMesh &initial, int nref) {
    int np0 = static_cast<int>(initial.nodes.size());
    int nt0 = static_cast<int>(initial.elems.size());

    Eigen::SparseMatrix<double> P_node(np0, np0);
    P_node.setIdentity();
    Eigen::SparseMatrix<double> P_elem(nt0, nt0);
    P_elem.setIdentity();
    Eigen::SparseMatrix<double> P_dg(3*nt0, 3*nt0);
    P_dg.setIdentity();

    TriMesh current = initial;

    for (int lev = 0; lev < nref; ++lev) {
        RefineOutput out = refine(current);
        P_node = out.P_node * P_node;
        P_elem = out.P_elem * P_elem;
        P_dg   = out.P_dg   * P_dg;
        current = std::move(out.mesh);
    }

    return {std::move(current), std::move(P_node), std::move(P_elem), std::move(P_dg)};
}

} // namespace lod2d

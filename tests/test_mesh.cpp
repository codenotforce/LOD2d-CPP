#include "mesh/refine.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>

using namespace lod2d;

namespace {

TriMesh unit_square() {
    TriMesh mesh;
    mesh.nodes = {{0,0}, {1,0}, {1,1}, {0,1}};
    mesh.elems = {{0,1,3}, {1,2,3}};
    mesh.dirichlet = {0,1,2,3};
    return mesh;
}

bool near(double a, double b, double tol = 1e-12) {
    return std::abs(a - b) <= tol;
}

Edge edge_of(int a, int b) {
    return a < b ? Edge{a, b} : Edge{b, a};
}

bool check(bool condition, const std::string &name, int &failed) {
    if (condition) {
        std::cout << "  [PASS] " << name << "\n";
        return true;
    }
    std::cout << "  [FAIL] " << name << "\n";
    ++failed;
    return false;
}

double total_area(const TriMesh &mesh) {
    double total = 0.0;
    for (double a : compute_area(mesh)) total += a;
    return total;
}

bool all_positive_area(const TriMesh &mesh) {
    for (double a : compute_area(mesh)) {
        if (a <= 1e-14) return false;
    }
    return true;
}

bool is_conforming_tri_mesh(const TriMesh &mesh) {
    std::map<Edge, int> counts;
    for (const auto &tri : mesh.elems) {
        ++counts[edge_of(tri[0], tri[1])];
        ++counts[edge_of(tri[1], tri[2])];
        ++counts[edge_of(tri[2], tri[0])];
    }
    for (const auto &[edge, count] : counts) {
        if (count < 1 || count > 2) return false;
    }
    return true;
}

bool pnode_reproduces_coordinates(const TriMesh &coarse, const RefineOutput &out) {
    Eigen::MatrixXd coarse_xy(coarse.nodes.size(), 2);
    for (int i = 0; i < static_cast<int>(coarse.nodes.size()); ++i) {
        coarse_xy(i, 0) = coarse.nodes[i].x();
        coarse_xy(i, 1) = coarse.nodes[i].y();
    }
    Eigen::MatrixXd fine_xy = out.P_node * coarse_xy;
    if (fine_xy.rows() != static_cast<int>(out.mesh.nodes.size())) return false;
    for (int i = 0; i < fine_xy.rows(); ++i) {
        if (!near(fine_xy(i, 0), out.mesh.nodes[i].x()) || !near(fine_xy(i, 1), out.mesh.nodes[i].y())) {
            return false;
        }
    }
    return true;
}

bool pdg_reproduces_child_vertices(const TriMesh &coarse, const RefineOutput &out) {
    Eigen::MatrixXd coarse_dg(3 * coarse.elems.size(), 2);
    for (int e = 0; e < static_cast<int>(coarse.elems.size()); ++e) {
        for (int i = 0; i < 3; ++i) {
            const auto &p = coarse.nodes[coarse.elems[e][i]];
            coarse_dg(3 * e + i, 0) = p.x();
            coarse_dg(3 * e + i, 1) = p.y();
        }
    }
    Eigen::MatrixXd fine_dg = out.P_dg * coarse_dg;
    if (fine_dg.rows() != 3 * static_cast<int>(out.mesh.elems.size())) return false;
    for (int e = 0; e < static_cast<int>(out.mesh.elems.size()); ++e) {
        for (int i = 0; i < 3; ++i) {
            const auto &p = out.mesh.nodes[out.mesh.elems[e][i]];
            if (!near(fine_dg(3 * e + i, 0), p.x()) || !near(fine_dg(3 * e + i, 1), p.y())) {
                return false;
            }
        }
    }
    return true;
}

} // namespace

int main() {
    std::cout << "=== Mesh Refinement Tests ===\n";
    int failed = 0;

    {
        TriMesh mesh = unit_square();
        auto [edges, is_bdy] = compute_edges(mesh);
        const int bdy_count = static_cast<int>(std::count(is_bdy.begin(), is_bdy.end(), true));
        check(edges.size() == 5, "Edge enumeration count", failed);
        check(bdy_count == 4, "Boundary edge count", failed);
    }

    {
        TriMesh mesh = unit_square();
        auto areas = compute_area(mesh);
        check(areas.size() == 2 && near(areas[0], 0.5) && near(areas[1], 0.5), "Area calculation", failed);
    }

    {
        TriMesh mesh = unit_square();
        RefineOutput out = refine(mesh);
        check(out.mesh.nodes.size() == 5, "Uniform LEB node count", failed);
        check(out.mesh.elems.size() == 4, "Uniform LEB element count", failed);
        check(out.P_node.rows() == 5 && out.P_node.cols() == 4, "Uniform LEB P_node size", failed);
        check(out.P_elem.rows() == 4 && out.P_elem.cols() == 2, "Uniform LEB P_elem size", failed);
        check(out.P_dg.rows() == 12 && out.P_dg.cols() == 6, "Uniform LEB P_dg size", failed);
        check(is_conforming_tri_mesh(out.mesh), "Uniform LEB conforming edges", failed);
        check(all_positive_area(out.mesh) && near(total_area(out.mesh), 1.0), "Uniform LEB area", failed);
        check(pnode_reproduces_coordinates(mesh, out), "Uniform LEB P_node coordinates", failed);
        check(pdg_reproduces_child_vertices(mesh, out), "Uniform LEB P_dg coordinates", failed);
    }

    {
        TriMesh mesh = unit_square();
        RefineOutput out = refine_marked(mesh, {0});
        check(out.mesh.nodes.size() == 5, "Local LEB shared midpoint", failed);
        check(out.mesh.elems.size() == 4, "Local LEB neighbor closure", failed);
        check(is_conforming_tri_mesh(out.mesh), "Local LEB conforming edges", failed);
        check(all_positive_area(out.mesh) && near(total_area(out.mesh), 1.0), "Local LEB area", failed);
        check(pnode_reproduces_coordinates(mesh, out), "Local LEB P_node coordinates", failed);
        check(pdg_reproduces_child_vertices(mesh, out), "Local LEB P_dg coordinates", failed);
    }

    {
        TriMesh mesh = unit_square();
        for (int nref : {2, 4, 6}) {
            RefineOutput out = refine_mesh(mesh, nref);
            check(is_conforming_tri_mesh(out.mesh), "Multi-level LEB conforming nref=" + std::to_string(nref), failed);
            check(all_positive_area(out.mesh) && near(total_area(out.mesh), 1.0, 1e-10), "Multi-level LEB area nref=" + std::to_string(nref), failed);
            check(out.P_node.cols() == static_cast<int>(mesh.nodes.size()), "Multi-level P_node cols nref=" + std::to_string(nref), failed);
            check(out.P_elem.cols() == static_cast<int>(mesh.elems.size()), "Multi-level P_elem cols nref=" + std::to_string(nref), failed);
        }
    }

    if (failed == 0) {
        std::cout << "\nAll tests passed!\n";
        return 0;
    }
    std::cout << "\nMesh tests failed: " << failed << "\n";
    return 1;
}

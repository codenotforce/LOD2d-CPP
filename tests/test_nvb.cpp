#include "mesh/refine.h"
#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using namespace lod2d;

namespace {

struct SparseGolden {
    int rows = 0;
    int cols = 0;
    int nnz = 0;
    std::vector<int> i;
    std::vector<int> j;
    std::vector<double> v;
};

struct CaseGolden {
    std::string name;
    std::vector<Point2> nodes;
    std::vector<Triangle> elems;
    std::vector<int> dirichlet;
    SparseGolden P_node;
    SparseGolden P_elem;
    SparseGolden P_dg;
};

TriMesh matlab_seed_mesh() {
    TriMesh mesh;
    mesh.nodes = {{0,0}, {1,0}, {1,1}, {0,1}};
    mesh.elems = {{0,1,2}, {0,2,3}};
    mesh.dirichlet = {0,1,2,3};
    return mesh;
}

TriMesh lod_seed_mesh() {
    TriMesh mesh;
    mesh.nodes = {{0,0}, {1,0}, {1,1}, {0,1}};
    mesh.elems = {{0,1,3}, {1,2,3}};
    mesh.dirichlet = {0,1,2,3};
    return mesh;
}

SparseGolden read_sparse(std::istream &in, const std::string &expected_name) {
    std::string name;
    SparseGolden g;
    in >> name >> g.rows >> g.cols >> g.nnz;
    if (name != expected_name) throw std::runtime_error("expected " + expected_name + ", got " + name);
    g.i.resize(g.nnz);
    g.j.resize(g.nnz);
    g.v.resize(g.nnz);
    for (int k = 0; k < g.nnz; ++k) in >> g.i[k] >> g.j[k] >> g.v[k];
    return g;
}

std::vector<CaseGolden> read_golden(const std::string &path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::vector<CaseGolden> cases;
    std::string tok;
    while (in >> tok) {
        if (tok != "CASE") throw std::runtime_error("expected CASE");
        CaseGolden c;
        in >> c.name;
        int n_nodes = 0, n_elems = 0, n_dir = 0;
        in >> tok >> n_nodes;
        if (tok != "NODES") throw std::runtime_error("expected NODES");
        in >> tok >> n_elems;
        if (tok != "ELEMS") throw std::runtime_error("expected ELEMS");
        in >> tok >> n_dir;
        if (tok != "DNODES") throw std::runtime_error("expected DNODES");
        in >> tok;
        if (tok != "COORDS") throw std::runtime_error("expected COORDS");
        c.nodes.resize(n_nodes);
        for (auto &p : c.nodes) in >> p.x() >> p.y();
        in >> tok;
        if (tok != "CONN") throw std::runtime_error("expected CONN");
        c.elems.resize(n_elems);
        for (auto &tri : c.elems) in >> tri[0] >> tri[1] >> tri[2];
        in >> tok;
        if (tok != "DIRICHLET") throw std::runtime_error("expected DIRICHLET");
        c.dirichlet.resize(n_dir);
        for (int &d : c.dirichlet) in >> d;
        c.P_node = read_sparse(in, "P_NODE");
        c.P_elem = read_sparse(in, "P_ELEM");
        c.P_dg = read_sparse(in, "P_DG");
        in >> tok;
        if (tok != "END_CASE") throw std::runtime_error("expected END_CASE");
        cases.push_back(std::move(c));
    }
    return cases;
}

bool near(double a, double b, double tol = 1e-12) {
    return std::abs(a - b) <= tol;
}

bool compare_mesh(const TriMesh &mesh, const CaseGolden &gold, std::string &why) {
    if (mesh.nodes.size() != gold.nodes.size()) {
        why = "node count";
        return false;
    }
    if (mesh.elems.size() != gold.elems.size()) {
        why = "element count";
        return false;
    }
    if (mesh.dirichlet != gold.dirichlet) {
        why = "dirichlet nodes";
        return false;
    }
    for (size_t i = 0; i < mesh.nodes.size(); ++i) {
        if (!near(mesh.nodes[i].x(), gold.nodes[i].x()) || !near(mesh.nodes[i].y(), gold.nodes[i].y())) {
            why = "node coordinates at " + std::to_string(i);
            return false;
        }
    }
    for (size_t i = 0; i < mesh.elems.size(); ++i) {
        if (mesh.elems[i] != gold.elems[i]) {
            why = "connectivity at " + std::to_string(i);
            return false;
        }
    }
    return true;
}

bool compare_sparse(const Eigen::SparseMatrix<double> &M, const SparseGolden &gold, std::string &why) {
    if (M.rows() != gold.rows || M.cols() != gold.cols) {
        why = "sparse dimensions";
        return false;
    }
    if (M.nonZeros() != gold.nnz) {
        why = "sparse nnz got " + std::to_string(M.nonZeros()) + " expected " + std::to_string(gold.nnz);
        return false;
    }
    for (int k = 0; k < gold.nnz; ++k) {
        const double value = M.coeff(gold.i[k], gold.j[k]);
        if (!near(value, gold.v[k])) {
            why = "sparse value at (" + std::to_string(gold.i[k]) + "," + std::to_string(gold.j[k]) + ")";
            return false;
        }
    }
    return true;
}

Edge sorted_edge(int a, int b) {
    return a < b ? Edge{a, b} : Edge{b, a};
}

bool no_hanging_nodes(const TriMesh &mesh) {
    std::map<Edge, int> edges;
    for (const auto &tri : mesh.elems) {
        ++edges[sorted_edge(tri[0], tri[1])];
        ++edges[sorted_edge(tri[1], tri[2])];
        ++edges[sorted_edge(tri[2], tri[0])];
    }
    for (const auto &[edge, count] : edges) {
        if (count < 1 || count > 2) return false;
        const auto &a = mesh.nodes[edge[0]];
        const auto &b = mesh.nodes[edge[1]];
        const auto ab = b - a;
        const double len2 = ab.squaredNorm();
        for (int k = 0; k < static_cast<int>(mesh.nodes.size()); ++k) {
            if (k == edge[0] || k == edge[1]) continue;
            const auto ak = mesh.nodes[k] - a;
            const double cross = std::abs(ab.x() * ak.y() - ab.y() * ak.x());
            const double s = ak.dot(ab) / len2;
            if (cross < 1e-12 && s > 1e-12 && s < 1.0 - 1e-12) return false;
        }
    }
    return true;
}

bool check_case(const CaseGolden &gold) {
    TriMesh seed = matlab_seed_mesh();
    RefineOutput out;
    if (gold.name == "global1") out = refine_nvb(seed);
    else if (gold.name == "local_elem0") out = bisect_newest_vertex(seed, {0});
    else if (gold.name == "global4") out = refine_mesh_nvb(seed, 4);
    else throw std::runtime_error("unknown case " + gold.name);

    std::string why;
    if (!compare_mesh(out.mesh, gold, why)) {
        std::cout << "  [FAIL] " << gold.name << " mesh: " << why << "\n";
        return false;
    }
    if (!compare_sparse(out.P_node, gold.P_node, why)) {
        std::cout << "  [FAIL] " << gold.name << " P_node: " << why << "\n";
        return false;
    }
    if (!compare_sparse(out.P_elem, gold.P_elem, why)) {
        std::cout << "  [FAIL] " << gold.name << " P_elem: " << why << "\n";
        return false;
    }
    if (!compare_sparse(out.P_dg, gold.P_dg, why)) {
        std::cout << "  [FAIL] " << gold.name << " P_dg: " << why << "\n";
        return false;
    }
    if (!no_hanging_nodes(out.mesh)) {
        std::cout << "  [FAIL] " << gold.name << " hanging node check\n";
        return false;
    }
    std::cout << "  [PASS] " << gold.name << "\n";
    return true;
}

} // namespace

int main() {
    std::cout << "=== Newest-Vertex Bisection Golden Tests ===\n";
    int failed = 0;
    try {
        const auto cases = read_golden("tests/golden_nvb.txt");
        for (const auto &gold : cases) {
            if (!check_case(gold)) ++failed;
        }

        TriMesh incompatible = lod_seed_mesh();
        auto global = refine_nvb(incompatible);
        auto local = bisect_newest_vertex(incompatible, {0});
        if (no_hanging_nodes(global.mesh) && no_hanging_nodes(local.mesh)) {
            std::cout << "  [PASS] incompatible seed closure smoke test\n";
        } else {
            std::cout << "  [FAIL] incompatible seed closure smoke test\n";
            ++failed;
        }
    } catch (const std::exception &e) {
        std::cout << "  [FAIL] exception: " << e.what() << "\n";
        ++failed;
    }

    if (failed == 0) {
        std::cout << "\nAll NVB tests passed!\n";
        return 0;
    }
    std::cout << "\nNVB tests failed: " << failed << "\n";
    return 1;
}

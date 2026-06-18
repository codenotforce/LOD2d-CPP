/// Phase B: Golden‑data verification against MATLAB reference.
/// Reads `golden_mesh.txt` and compares C++ mesh output field‑by‑field.

#include "mesh/refine.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cassert>
#include <vector>

using namespace lod2d;

// ---- Simple golden‑file parser ----
struct GoldenLevel {
    int nref;
    int n_nodes, n_elems, n_dnodes;
    std::vector<Point2> coords;
    std::vector<Triangle> elems;
    std::vector<int> dirichlet;

    struct SparseMat {
        int rows, cols;
        int nnz;
        std::vector<int>    i, j;
        std::vector<double> v;
    };
    SparseMat P1, P0, P1dg;
};

static std::vector<GoldenLevel> parse_golden(const std::string &path) {
    std::ifstream f(path);
    if (!f) { std::cerr << "Cannot open " << path << "\n"; std::exit(1); }

    std::vector<GoldenLevel> levels;
    std::string line;
    GoldenLevel *L = nullptr;

    auto next_int = [&]() { int x; f >> x; f.ignore(); return x; };

    while (std::getline(f, line)) {
        if (line.rfind("LEVEL ", 0) == 0) {
            levels.emplace_back();
            L = &levels.back();
            L->nref = std::stoi(line.substr(6));
        }
        else if (line.rfind("NODES ", 0) == 0) L->n_nodes = std::stoi(line.substr(6));
        else if (line.rfind("ELEMS ", 0) == 0) L->n_elems = std::stoi(line.substr(6));
        else if (line.rfind("DNODES ", 0) == 0) L->n_dnodes = std::stoi(line.substr(7));
        else if (line == "COORDS") {
            L->coords.resize(L->n_nodes);
            for (int i = 0; i < L->n_nodes; ++i)
                f >> L->coords[i].x() >> L->coords[i].y();
            f.ignore();
        }
        else if (line == "CONN") {
            L->elems.resize(L->n_elems);
            for (int i = 0; i < L->n_elems; ++i) {
                int a, b, c; f >> a >> b >> c; f.ignore();
                L->elems[i] = {a, b, c};
            }
        }
        else if (line == "DIRICHLET") {
            L->dirichlet.resize(L->n_dnodes);
            for (int i = 0; i < L->n_dnodes; ++i) f >> L->dirichlet[i];
        }
        else if (line.rfind("P1_NNZ ", 0) == 0) {
            int nnz = std::stoi(line.substr(7));
            L->P1.nnz = nnz; L->P1.i.resize(nnz); L->P1.j.resize(nnz); L->P1.v.resize(nnz);
            for (int k = 0; k < nnz; ++k) f >> L->P1.i[k] >> L->P1.j[k] >> L->P1.v[k];
        }
        else if (line.rfind("P0_NNZ ", 0) == 0) {
            int nnz = std::stoi(line.substr(7));
            L->P0.nnz = nnz; L->P0.i.resize(nnz); L->P0.j.resize(nnz); L->P0.v.resize(nnz);
            for (int k = 0; k < nnz; ++k) f >> L->P0.i[k] >> L->P0.j[k] >> L->P0.v[k];
        }
        else if (line.rfind("P1DG_NNZ ", 0) == 0) {
            int nnz = std::stoi(line.substr(9));
            L->P1dg.nnz = nnz; L->P1dg.i.resize(nnz); L->P1dg.j.resize(nnz); L->P1dg.v.resize(nnz);
            for (int k = 0; k < nnz; ++k) f >> L->P1dg.i[k] >> L->P1dg.j[k] >> L->P1dg.v[k];
        }
    }
    return levels;
}

// ---- Helpers ----
static bool near(double a, double b, double tol = 1e-12) {
    return std::abs(a - b) < tol || std::abs(a - b) / std::max(1.0, std::abs(a)) < tol;
}

// ===================================================================
int main() {
    std::cout << "=== Phase B: Golden‑data Verification ===\n\n";

    std::string golden_path = "tests/golden_mesh.txt";
    auto levels = parse_golden(golden_path);
    std::cout << "Loaded " << levels.size() << " refinement levels\n";

    // Build C++ initial mesh (same as MATLAB: unit square, 2 triangles)
    TriMesh mesh;
    mesh.nodes = {{0,0}, {1,0}, {1,1}, {0,1}};
    mesh.elems = {{0,1,3}, {1,2,3}};
    mesh.dirichlet = {0,1,2,3};

    int passed = 0, failed = 0;

    for (const auto &gold : levels) {
        std::cout << "\n--- Level " << gold.nref << " ---\n";

        RefineOutput out = refine_mesh(mesh, gold.nref);
        const auto &fine = out.mesh;

        // ---- 1. Node count ----
        if (fine.nodes.size() == static_cast<size_t>(gold.n_nodes)) {
            std::cout << "  [PASS] node count: " << fine.nodes.size() << "\n"; ++passed;
        } else {
            std::cout << "  [FAIL] node count: got " << fine.nodes.size()
                      << " expected " << gold.n_nodes << "\n"; ++failed;
        }

        // ---- 2. Element count ----
        if (fine.elems.size() == static_cast<size_t>(gold.n_elems)) {
            std::cout << "  [PASS] elem count: " << fine.elems.size() << "\n"; ++passed;
        } else {
            std::cout << "  [FAIL] elem count: got " << fine.elems.size()
                      << " expected " << gold.n_elems << "\n"; ++failed;
        }

        // ---- 3. Node coordinates ----
        bool coords_ok = true;
        double max_err = 0;
        for (int i = 0; i < gold.n_nodes; ++i) {
            double dx = std::abs(fine.nodes[i].x() - gold.coords[i].x());
            double dy = std::abs(fine.nodes[i].y() - gold.coords[i].y());
            if (dx > max_err) max_err = dx;
            if (dy > max_err) max_err = dy;
            if (dx > 1e-12 || dy > 1e-12) { coords_ok = false; break; }
        }
        if (coords_ok) {
            std::cout << "  [PASS] coords match (max err " << max_err << ")\n"; ++passed;
        } else {
            std::cout << "  [FAIL] coords mismatch (max err " << max_err << ")\n"; ++failed;
        }

        // ---- 4. Element connectivity ----
        bool conn_ok = true;
        for (int i = 0; i < gold.n_elems; ++i) {
            if (fine.elems[i] != gold.elems[i]) { conn_ok = false; break; }
        }
        if (conn_ok) {
            std::cout << "  [PASS] connectivity exact match\n"; ++passed;
        } else {
            std::cout << "  [FAIL] connectivity differs\n"; ++failed;
        }

        // ---- 5. Dirichlet nodes ----
        bool dir_ok = true;
        for (int i = 0; i < gold.n_dnodes; ++i) {
            if (fine.dirichlet[i] != gold.dirichlet[i]) { dir_ok = false; break; }
        }
        if (dir_ok) {
            std::cout << "  [PASS] Dirichlet nodes match\n"; ++passed;
        } else {
            std::cout << "  [FAIL] Dirichlet nodes differ\n"; ++failed;
        }

        // ---- 6. Total area ----
        auto areas = compute_area(fine);
        double total = 0;
        for (double a : areas) total += a;
        if (std::abs(total - 1.0) < 1e-12) {
            std::cout << "  [PASS] total area = " << total << "\n"; ++passed;
        } else {
            std::cout << "  [FAIL] total area = " << total << " (expected 1.0)\n"; ++failed;
        }

        // ---- 7. Edge count ----
        auto [edges, is_bdy] = compute_edges(fine);
        int n_bdy = 0;
        for (bool b : is_bdy) if (b) ++n_bdy;
        std::cout << "  [INFO] edges: " << edges.size() << " (boundary: " << n_bdy << ")\n";

        // ---- 8. P_node prolongation ----
        int cpp_nnz = out.P_node.nonZeros();
        bool p1_ok = (cpp_nnz == gold.P1.nnz);
        if (!p1_ok) {
            std::cout << "  [FAIL] P_node nnz: got " << cpp_nnz
                      << " expected " << gold.P1.nnz << "\n"; ++failed;
        } else {
            double max_v = 0;
            for (int k = 0; k < gold.P1.nnz; ++k) {
                double v = out.P_node.coeff(gold.P1.i[k], gold.P1.j[k]);
                double d = std::abs(v - gold.P1.v[k]);
                if (d > max_v) max_v = d;
                if (d > 1e-12) { p1_ok = false; break; }
            }
            if (p1_ok) {
                std::cout << "  [PASS] P_node: " << cpp_nnz << " nnz (max diff " << max_v << ")\n"; ++passed;
            } else {
                std::cout << "  [FAIL] P_node: value mismatch\n"; ++failed;
            }
        }

        // ---- 9. P_elem prolongation ----
        int pe_nnz = out.P_elem.nonZeros();
        if (pe_nnz == gold.P0.nnz) {
            std::cout << "  [PASS] P_elem: " << pe_nnz << " nnz\n"; ++passed;
        } else {
            std::cout << "  [FAIL] P_elem nnz: got " << pe_nnz
                      << " expected " << gold.P0.nnz << "\n"; ++failed;
        }

        // ---- 10. P_dg prolongation ----
        int pdg_nnz = out.P_dg.nonZeros();
        bool pdg_ok = (pdg_nnz == gold.P1dg.nnz);
        if (!pdg_ok) {
            std::cout << "  [FAIL] P_dg nnz: got " << pdg_nnz
                      << " expected " << gold.P1dg.nnz << "\n"; ++failed;
        } else {
            double max_v = 0;
            for (int k = 0; k < gold.P1dg.nnz; ++k) {
                double v = out.P_dg.coeff(gold.P1dg.i[k], gold.P1dg.j[k]);
                double d = std::abs(v - gold.P1dg.v[k]);
                if (d > max_v) max_v = d;
                if (d > 1e-12) { pdg_ok = false; break; }
            }
            if (pdg_ok) {
                std::cout << "  [PASS] P_dg: " << pdg_nnz << " nnz (max diff " << max_v << ")\n"; ++passed;
            } else {
                std::cout << "  [FAIL] P_dg: value mismatch\n"; ++failed;
            }
        }
    }

    std::cout << "\n========================================\n";
    std::cout << "Phase B RESULTS: " << passed << " PASS, " << failed << " FAIL\n";
    std::cout << "========================================\n";
    return (failed > 0) ? 1 : 0;
}

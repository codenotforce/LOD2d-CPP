/// Phase C: Golden‑data verification for DG stiffness assembly.

#include "fem/assemble_dg.h"
#include "mesh/refine.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>

using namespace lod2d;

struct GoldenDG {
    int nref, N, NT, nnz;
    std::vector<double> coeff;
    std::vector<int>    i, j;
    std::vector<double> v;
};

static std::vector<GoldenDG> parse_dg_golden(const std::string &path) {
    std::ifstream f(path);
    if (!f) { std::cerr << "Cannot open " << path << "\n"; std::exit(1); }
    std::vector<GoldenDG> levels;
    std::string line;
    GoldenDG *L = nullptr;

    auto next_int = [&]() { int x; f >> x; f.ignore(); return x; };

    while (std::getline(f, line)) {
        if (line.rfind("LEVEL ", 0) == 0) {
            levels.emplace_back(); L = &levels.back();
            L->nref = std::stoi(line.substr(6));
        }
        else if (line.rfind("N ", 0) == 0)   L->N   = std::stoi(line.substr(2));
        else if (line.rfind("NT ", 0) == 0)  L->NT  = std::stoi(line.substr(3));
        else if (line.rfind("NNZ ", 0) == 0) L->nnz = std::stoi(line.substr(4));
        else if (line == "COEFF") {
            L->coeff.resize(L->NT);
            for (int i = 0; i < L->NT; ++i) f >> L->coeff[i];
            f.ignore();
        }
        else if (line == "MATRIX") {
            L->i.resize(L->nnz); L->j.resize(L->nnz); L->v.resize(L->nnz);
            for (int k = 0; k < L->nnz; ++k) f >> L->i[k] >> L->j[k] >> L->v[k];
        }
    }
    return levels;
}

int main() {
    std::cout << "=== Phase C: DG Assembly Golden‑data Verification ===\n\n";

    auto levels = parse_dg_golden("tests/golden_dg.txt");
    std::cout << "Loaded " << levels.size() << " refinement levels\n";

    TriMesh mesh;
    mesh.nodes = {{0,0}, {1,0}, {1,1}, {0,1}};
    mesh.elems = {{0,1,3}, {1,2,3}};
    mesh.dirichlet = {0,1,2,3};

    int passed = 0, failed = 0;

    for (const auto &gold : levels) {
        RefineOutput out = refine_mesh(mesh, gold.nref);
        const auto &fine = out.mesh;

        auto S = assemble_dg(fine, gold.coeff);

        std::cout << "\n--- Level " << gold.nref << " ---\n";
        std::cout << "  C++: " << S.rows() << "x" << S.cols() << " nnz=" << S.nonZeros() << "\n";
        std::cout << "  Gold: " << gold.N << "x" << gold.N << " nnz=" << gold.nnz << "\n";

        // 1. Size
        bool sz_ok = (S.rows() == gold.N && S.cols() == gold.N);
        std::cout << "  " << (sz_ok ? "[PASS]" : "[FAIL]") << " dimensions\n";
        if (sz_ok) ++passed; else ++failed;

        // 2. nnz
        bool nz_ok = (S.nonZeros() == gold.nnz);
        std::cout << "  " << (nz_ok ? "[PASS]" : "[FAIL]") << " nnz\n";
        if (nz_ok) ++passed; else ++failed;

        // 3. All entries
        bool val_ok = true;
        double max_err = 0;
        for (int k = 0; k < gold.nnz; ++k) {
            double v = S.coeff(gold.i[k], gold.j[k]);
            double d = std::abs(v - gold.v[k]);
            if (d > max_err) max_err = d;
            if (d > 1e-12) { val_ok = false; break; }
        }
        std::cout << "  " << (val_ok ? "[PASS]" : "[FAIL]") << " values (max err " << max_err << ")\n";
        if (val_ok) ++passed; else ++failed;
    }

    // ---- Deterministic A=1 test ----
    std::cout << "\n--- Deterministic A=1 test ---\n";
    {
        TriMesh simple;
        simple.nodes = {{0,0}, {1,0}, {0,1}};
        simple.elems = {{0,1,2}};
        std::vector<double> coeff = {1.0};
        auto S = assemble_dg(simple, coeff);

        double expected[9] = {1.0, -0.5, -0.5, -0.5, 0.5, 0.0, -0.5, 0.0, 0.5};
        bool ok = true;
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                if (std::abs(S.coeff(j, k) - expected[j*3 + k]) > 1e-12) ok = false;
        std::cout << "  " << (ok ? "[PASS]" : "[FAIL]") << " A=1 right triangle\n";
        if (ok) ++passed; else ++failed;
    }

    std::cout << "\n========================================\n";
    std::cout << "Phase C RESULTS: " << passed << " PASS, " << failed << " FAIL\n";
    std::cout << "========================================\n";
    return (failed > 0) ? 1 : 0;
}

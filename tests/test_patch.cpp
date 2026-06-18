/// Phase E: Golden‑data verification for element patches.

#include "lod/patches.h"
#include "mesh/refine.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace lod2d;

int main() {
    std::cout << "=== Phase E: Patch Construction Golden‑data ===\n\n";

    std::ifstream f("tests/golden_patch.txt");
    if (!f) { std::cerr << "Cannot open golden_patch.txt\n"; return 1; }

    TriMesh T0;
    T0.nodes = {{0,0},{1,0},{1,1},{0,1}};
    T0.elems = {{0,1,3},{1,2,3}};
    T0.dirichlet = {0,1,2,3};

    int passed = 0, failed = 0;
    std::string line;

    while (std::getline(f, line)) {
        int Hlevel, ell, NTH, nnz;
        if (sscanf(line.c_str(), "H %d L %d NTH %d NNZ %d", &Hlevel, &ell, &NTH, &nnz) != 4) continue;

        std::vector<int> gi(nnz), gj(nnz), gv(nnz);
        for (int k = 0; k < nnz; ++k) f >> gi[k] >> gj[k] >> gv[k];
        f.ignore();
        std::string end_tag; f >> end_tag;  // END

        auto coarse_out = refine_mesh(T0, Hlevel);
        auto patch = build_patches(coarse_out.mesh, ell);

        std::cout << "H=" << Hlevel << " L=" << ell << " NTH=" << NTH << "\n";
        std::cout << "  Gold nnz=" << nnz << "  C++ nnz=" << patch.nonZeros() << "\n";

        // Check all golden positions
        double max_err = 0;
        int matched = 0;
        for (int k = 0; k < nnz; ++k) {
            double v = patch.coeff(gi[k], gj[k]);
            double d = std::abs(v - gv[k]);
            if (d > max_err) max_err = d;
            if (d < 1e-12) matched++;
        }

        bool ok = (matched == nnz);
        std::cout << "  " << (ok ? "[PASS]" : "[FAIL]") << " " << matched << "/" << nnz
                  << " match (max err " << max_err << ")\n";
        ok ? ++passed : ++failed;
    }

    std::cout << "\n========================================\n";
    std::cout << "Phase E RESULTS: " << passed << " PASS, " << failed << " FAIL\n";
    std::cout << "========================================\n";
    return (failed > 0) ? 1 : 0;
}

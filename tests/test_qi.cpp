/// Phase D: Golden‑data verification for quasi‑interpolation.

#include "lod/quasi_interp.h"
#include "mesh/refine.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace lod2d;

int main() {
    std::cout << "=== Phase D: Quasi‑Interpolation Golden‑data ===\n\n";

    // Read golden data
    std::ifstream f("tests/golden_qi.txt");
    if (!f) { std::cerr << "Cannot open golden_qi.txt\n"; return 1; }

    int NH, NTH, Nh, NTh_f, nnz;
    std::string tok;
    f >> tok >> NH >> tok >> NTH >> tok >> Nh >> tok >> NTh_f >> tok >> nnz;

    std::vector<int> gi(nnz), gj(nnz);
    std::vector<double> gv(nnz);
    f >> tok;  // MATRIX
    for (int k = 0; k < nnz; ++k) f >> gi[k] >> gj[k] >> gv[k];

    std::cout << "Golden IH: " << NH << "x" << Nh << " nnz=" << nnz << "\n";

    // Build same meshes
    TriMesh T0;
    T0.nodes = {{0,0},{1,0},{1,1},{0,1}};
    T0.elems = {{0,1,3},{1,2,3}};
    T0.dirichlet = {0,1,2,3};

    int Hlevel = 3, hlevel = 5;
    auto coarse_out = refine_mesh(T0, Hlevel);
    auto fine_out   = refine_mesh(coarse_out.mesh, hlevel - Hlevel);
    const auto &coarse = coarse_out.mesh;
    const auto &fine   = fine_out.mesh;
    const auto &P1dg   = fine_out.P_dg;

    // Build cg2dgh (fine DG → CG)
    int Nhdg = static_cast<int>(fine.elems.size()) * 3;
    std::vector<Eigen::Triplet<double>> cg_t;
    for (size_t e = 0; e < fine.elems.size(); ++e) {
        for (int i = 0; i < 3; ++i) {
            int v = fine.elems[e][i];
            cg_t.emplace_back(3*static_cast<int>(e) + i, v, 1.0);
        }
    }
    Eigen::SparseMatrix<double> cg2dgh(Nhdg, static_cast<int>(fine.nodes.size()));
    cg2dgh.setFromTriplets(cg_t.begin(), cg_t.end());

    // Compute IH
    auto IH = build_quasi_interp(coarse, fine, P1dg, cg2dgh,
                                  static_cast<int>(fine.nodes.size()),
                                  static_cast<int>(coarse.nodes.size()));

    std::cout << "C++ IH: " << IH.rows() << "x" << IH.cols() << " nnz=" << IH.nonZeros() << "\n";

    int passed = 0, failed = 0;

    // 1. Dimensions
    bool sz_ok = (IH.rows() == NH && IH.cols() == Nh);
    std::cout << "  " << (sz_ok ? "[PASS]" : "[FAIL]") << " dimensions\n";
    sz_ok ? ++passed : ++failed;

    // 2. Non‑zero structure
    double max_err = 0;
    int matched = 0, missing = 0;
    for (int k = 0; k < nnz; ++k) {
        double v = IH.coeff(gi[k], gj[k]);
        double d = std::abs(v - gv[k]);
        if (d > max_err) max_err = d;
        if (d > 1e-12) missing++;
        else matched++;
    }
    std::cout << "  [INFO] golden positions matched: " << matched << "/" << nnz
              << " (max err " << max_err << ")\n";
    std::cout << "  [INFO] C++ nnz: " << IH.nonZeros() << "  Gold nnz: " << nnz << "\n";

    bool val_ok = (missing == 0);
    std::cout << "  " << (val_ok ? "[PASS]" : "[FAIL]") << " values ("
              << matched << "/" << nnz << " match, " << missing << " missing)\n";
    val_ok ? ++passed : ++failed;

    std::cout << "\n========================================\n";
    std::cout << "Phase D RESULTS: " << passed << " PASS, " << failed << " FAIL\n";
    std::cout << "========================================\n";
    return (failed > 0) ? 1 : 0;
}

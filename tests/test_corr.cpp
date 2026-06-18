/// Phase F: Golden‑data verification for element corrector.

#include "lod/corrector.h"
#include "lod/quasi_interp.h"
#include "lod/patches.h"
#include "fem/assemble_dg.h"
#include "mesh/refine.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace lod2d;

int main() {
    std::cout << "=== Phase F: Corrector Golden‑data ===\n\n";

    std::ifstream f("tests/golden_corr.txt");
    if (!f) { std::cerr << "Cannot open\n"; return 1; }

    int NH, NTH, Nh, NTh_f;
    std::string tok;
    f >> tok >> NH >> tok >> NTH >> tok >> Nh >> tok >> NTh_f;

    // Read Ah
    int nAh; f >> tok >> nAh;
    std::vector<double> Ah(nAh);
    for (int i = 0; i < nAh; ++i) f >> Ah[i];

    // Read nngH
    int nn; f >> tok >> nn;
    std::vector<int> nngH(nn);
    for (int i = 0; i < nn; ++i) f >> nngH[i];

    // Read nngh
    f >> tok >> nn;
    std::vector<int> nngh(nn);
    for (int i = 0; i < nn; ++i) f >> nngh[i];

    // Build mesh
    TriMesh T0;
    T0.nodes={{0,0},{1,0},{1,1},{0,1}};
    T0.elems={{0,1,3},{1,2,3}}; T0.dirichlet={0,1,2,3};

    int Hlevel=3, hlevel=5, ell=2, d=2;
    auto c_out = refine_mesh(T0, Hlevel);
    auto f_out = refine_mesh(c_out.mesh, hlevel-Hlevel);
    const auto &coarse = c_out.mesh;
    const auto &fine   = f_out.mesh;

    // dghidx, dgHidx
    std::vector<std::array<int,3>> dghidx(fine.elems.size());
    for (size_t e=0; e<fine.elems.size(); ++e) for (int i=0;i<3;++i) dghidx[e][i]=3*e+i;
    std::vector<std::array<int,3>> dgHidx(coarse.elems.size());
    for (size_t e=0; e<coarse.elems.size(); ++e) for (int i=0;i<3;++i) dgHidx[e][i]=3*e+i;

    // cg2dgh, Shdg, IH, patch, P0, P1dg
    int Nhdg=3*fine.elems.size();
    std::vector<Eigen::Triplet<double>> cg_t;
    for (size_t e=0; e<fine.elems.size(); ++e)
        for (int i=0; i<3; ++i)
            cg_t.emplace_back(3*(int)e+i, fine.elems[e][i], 1.0);
    Eigen::SparseMatrix<double> cg2dgh(Nhdg, Nh);
    cg2dgh.setFromTriplets(cg_t.begin(), cg_t.end());

    Eigen::SparseMatrix<double> Shdg = assemble_dg(fine, Ah);
    auto IH    = build_quasi_interp(coarse, fine, f_out.P_dg, cg2dgh, Nh, NH);
    auto patch = build_patches(coarse, ell);

    int passed=0, failed=0;
    std::string line;
    while (f >> tok && tok == "ELEM") {
        int k, nnz;
        f >> k >> tok >> nnz;
        std::vector<int> gi(nnz), gj(nnz);
        std::vector<double> gv(nnz);
        for (int n=0; n<nnz; ++n) f >> gi[n] >> gj[n] >> gv[n];
        f >> tok; // END_ELEM

        auto CTk = compute_corrector(k, patch, coarse, NH, nngH,
            f_out.P_elem, fine, Nh, nngh, dghidx, cg2dgh, Shdg,
            f_out.P_dg, dgHidx, IH, d);

        double max_err=0; int matched=0;
        for (int n=0; n<nnz; ++n) {
            double v = CTk.coeff(gi[n], gj[n]);
            double diff = std::abs(v - gv[n]);
            if (diff > max_err) max_err = diff;
            if (diff < 1e-10) matched++;
        }
        bool ok = (matched == nnz);
        std::cout << "Elem " << k << ": " << (ok?"[PASS]":"[FAIL]")
                  << " " << matched << "/" << nnz << " match (max err " << max_err << ")\n";
        ok ? ++passed : ++failed;
    }

    std::cout << "\n========================================\n";
    std::cout << "Phase F: " << passed << " PASS, " << failed << " FAIL\n";
    std::cout << "========================================\n";
    return (failed > 0) ? 1 : 0;
}

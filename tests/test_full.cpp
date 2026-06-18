/// Full end‑to‑end LOD test — compares C++ results & timing with MATLAB.
#include "lod/corrector.h"
#include "lod/quasi_interp.h"
#include "lod/patches.h"
#include "fem/assemble_dg.h"
#include "mesh/refine.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <chrono>
#include <Eigen/Dense>

using namespace lod2d;
namespace chr = std::chrono;

int main() {
    std::cout << "=== Full LOD Pipeline — C++ vs MATLAB ===\n\n";

    std::ifstream f("tests/golden_full.txt");
    if (!f) { std::cerr << "Cannot open\n"; return 1; }

    int NH, NTH, Nh, NTh_f;
    std::string tok;
    f >> tok >> NH >> tok >> NTH >> tok >> Nh >> tok >> NTh_f;

    int nAh; f >> tok >> nAh;
    std::vector<double> Ah(nAh);
    for (int i=0; i<nAh; ++i) f >> Ah[i];

    int nn; std::vector<int> nngH, nngh;
    f >> tok >> nn; nngH.resize(nn);
    for (int i=0; i<nn; ++i) f >> nngH[i];
    f >> tok >> nn; nngh.resize(nn);
    for (int i=0; i<nn; ++i) f >> nngh[i];

    // Golden solutions
    std::vector<double> guH(NH), guHms(Nh), guh(Nh);
    f >> tok >> nn; for (int i=0; i<nn; ++i) f >> guH[i];
    f >> tok >> nn; for (int i=0; i<nn; ++i) f >> guHms[i];
    f >> tok >> nn; for (int i=0; i<nn; ++i) f >> guh[i];
    double gE, gL2, gFE;
    f >> tok >> gE >> gL2 >> gFE;

    int Hlevel=3, hlevel=5, ell=2, d=2;
    TriMesh T0;
    T0.nodes={{0,0},{1,0},{1,1},{0,1}};
    T0.elems={{0,1,3},{1,2,3}}; T0.dirichlet={0,1,2,3};

    auto t0 = chr::high_resolution_clock::now();

    // Mesh
    auto c_out = refine_mesh(T0, Hlevel);
    auto f_out = refine_mesh(c_out.mesh, hlevel - Hlevel);
    const auto &coarse = c_out.mesh;
    const auto &fine   = f_out.mesh;

    // dghidx, dgHidx
    std::vector<std::array<int,3>> dghidx(NTh_f);
    for (int e=0; e<NTh_f; ++e) for (int i=0;i<3;++i) dghidx[e][i]=3*e+i;
    std::vector<std::array<int,3>> dgHidx(NTH);
    for (int e=0; e<NTH; ++e) for (int i=0;i<3;++i) dgHidx[e][i]=3*e+i;

    // cg2dgh
    int Nhdg = 3*NTh_f;
    std::vector<Eigen::Triplet<double>> cg_t;
    for (int e=0; e<NTh_f; ++e)
        for (int i=0; i<3; ++i)
            cg_t.emplace_back(3*e+i, fine.elems[e][i], 1.0);
    Eigen::SparseMatrix<double> cg2dgh(Nhdg, Nh);
    cg2dgh.setFromTriplets(cg_t.begin(), cg_t.end());

    // Operators
    auto element_stiffness = assemble_element_stiffness(fine, Ah);
    Eigen::SparseMatrix<double> Shdg = assemble_dg_from_element_stiffness(element_stiffness);
    auto IH    = build_quasi_interp(coarse, fine, f_out.P_dg, cg2dgh, Nh, NH);
    auto patch = build_patches(coarse, ell);
    auto fine_element_children = build_fine_element_children(f_out.P_elem, NTH);

    // Mass matrix Mhdg for coarse solve
    auto areas = compute_area(fine);
    double M3[3][3]={{2,1,1},{1,2,1},{1,1,2}};
    std::vector<Eigen::Triplet<double>> mh_t;
    for (int e=0; e<NTh_f; ++e) {
        double s=areas[e]/12.0;
        for (int i=0;i<3;++i) for (int j=0;j<3;++j)
            mh_t.emplace_back(3*e+i, 3*e+j, s*M3[i][j]);
    }
    Eigen::SparseMatrix<double> Mhdg(3*NTh_f, 3*NTh_f);
    Mhdg.setFromTriplets(mh_t.begin(), mh_t.end());

    // Correctors
    std::vector<Eigen::SparseMatrix<double>> CT(NTH);
    #pragma omp parallel for
    for (int k=0; k<NTH; ++k) {
        CT[k] = compute_corrector(k, patch, coarse, NH, nngH,
            f_out.P_elem, fine, Nh, nngh, dghidx, cg2dgh, Shdg,
            f_out.P_dg, dgHidx, IH, d, CorrectorSolver::EigenLLT,
            &element_stiffness, &fine_element_children);
    }

    // Global correction C_ell = cell2mat(CT) * cg2dgH
    // cg2dgH: NTH*3 × NH, coarse DG→CG
    int NHdg_c = 3*NTH;
    std::vector<Eigen::Triplet<double>> cg2dg_t;
    for (int e=0; e<NTH; ++e)
        for (int i=0; i<3; ++i)
            cg2dg_t.emplace_back(3*e+i, coarse.elems[e][i], 1.0);
    Eigen::SparseMatrix<double> cg2dgH(NHdg_c, NH);
    cg2dgH.setFromTriplets(cg2dg_t.begin(), cg2dg_t.end());

    // C_ell = [CT0 CT1 ... CT_{NTH-1}] * cg2dgH
    // Build horizontal concatenation
    std::vector<Eigen::Triplet<double>> cell_t;
    for (int k=0; k<NTH; ++k) {
        for (int c=0; c<CT[k].outerSize(); ++c)
            for (Eigen::SparseMatrix<double>::InnerIterator it(CT[k], c); it; ++it)
                cell_t.emplace_back(it.row(), k*(d+1) + static_cast<int>(it.col()), it.value());
    }
    Eigen::SparseMatrix<double> cell_mat(Nh, NTH*(d+1));
    cell_mat.setFromTriplets(cell_t.begin(), cell_t.end());
    Eigen::SparseMatrix<double> C_ell = cell_mat * cg2dgH;

    // G = P1 - C_ell
    Eigen::SparseMatrix<double> G = f_out.P_node - C_ell;

    // Coarse LOD system
    std::vector<int> dofH;
    for (int i=0; i<NH; ++i) {
        bool is_dir = false;
        for (int dv : coarse.dirichlet) if (dv == i) { is_dir=true; break; }
        if (!is_dir) dofH.push_back(i);
    }
    int nFree = static_cast<int>(dofH.size());

    // G0 = G(:, dofH)
    std::vector<Eigen::Triplet<double>> g0_t;
    for (int k=0; k<G.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(G, k); it; ++it)
            for (int j=0; j<nFree; ++j)
                if (static_cast<int>(it.col()) == dofH[j])
                    g0_t.emplace_back(it.row(), j, it.value());
    Eigen::SparseMatrix<double> G0(Nh, nFree);
    G0.setFromTriplets(g0_t.begin(), g0_t.end());

    // SHLOD0 = G0' * cg2dgh' * Shdg * cg2dgh * G0
    Eigen::SparseMatrix<double> T = cg2dgh.transpose() * Shdg * cg2dgh;
    Eigen::SparseMatrix<double> SHLOD0 = G0.transpose() * T * G0;

    // RHS: G0' * cg2dgh' * Mhdg * cg2dgh * (P1 * f)
    // f(TH.p) = ones(NH,1); P1 * f_coarse = fine-scale prolongation
    Eigen::VectorXd f_coarse = Eigen::VectorXd::Ones(NH);
    Eigen::VectorXd rhs = G0.transpose() * (cg2dgh.transpose() * (Mhdg * (cg2dgh * (f_out.P_node * f_coarse))));

    // Solve
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llts(SHLOD0);
    Eigen::VectorXd uH_free = llts.solve(rhs);
    Eigen::VectorXd uH = Eigen::VectorXd::Zero(NH);
    for (int j=0; j<nFree; ++j) uH(dofH[j]) = uH_free(j);

    // uHms = G * uH
    Eigen::VectorXd uHms = G * uH;

    // Reference
    std::vector<int> dofh;
    for (int i=0; i<Nh; ++i) {
        bool is_dir = false;
        for (int dv : fine.dirichlet) if (dv==i) { is_dir=true; break; }
        if (!is_dir) dofh.push_back(i);
    }
    int nFree_f = static_cast<int>(dofh.size());
    Eigen::SparseMatrix<double> Sh = cg2dgh.transpose() * Shdg * cg2dgh;
    Eigen::SparseMatrix<double> Mh = cg2dgh.transpose() * Mhdg * cg2dgh;

    // Sh_free = Sh(dofh, dofh), Mh_free = Mh(dofh, :)
    std::vector<Eigen::Triplet<double>> sh_t, mh_free;
    for (int k=0; k<Sh.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(Sh, k); it; ++it) {
            int ri=-1, ci=-1;
            for (int j=0; j<nFree_f; ++j) {
                if (dofh[j]==static_cast<int>(it.row())) ri=j;
                if (dofh[j]==static_cast<int>(it.col())) ci=j;
            }
            if (ri>=0 && ci>=0) sh_t.emplace_back(ri, ci, it.value());
        }
    Eigen::SparseMatrix<double> Sh_free(nFree_f, nFree_f);
    Sh_free.setFromTriplets(sh_t.begin(), sh_t.end());

    // RHS_ref = Mh(dofh, :) * f
    Eigen::VectorXd rhs_ref = Eigen::VectorXd::Zero(nFree_f);
    for (int j=0; j<nFree_f; ++j) {
        double sum=0;
        for (int k2=0; k2<Mh.outerSize(); ++k2)
            for (Eigen::SparseMatrix<double>::InnerIterator it(Mh, k2); it; ++it)
                if (static_cast<int>(it.row())==dofh[j]) sum += it.value();  // f=1 at all fine vertices
        rhs_ref(j) = sum;
    }
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt_ref(Sh_free);
    Eigen::VectorXd uh_free = llt_ref.solve(rhs_ref);
    Eigen::VectorXd uh = Eigen::VectorXd::Zero(Nh);
    for (int j=0; j<nFree_f; ++j) uh(dofh[j]) = uh_free(j);

    // Errors
    Eigen::VectorXd diff = uh - uHms;
    double errE  = std::sqrt(diff.dot(Sh * diff));
    double errL2 = std::sqrt(diff.dot(Mh * diff));
    Eigen::VectorXd diff_fe = uh - f_out.P_node * uH;
    double errFE = std::sqrt(diff_fe.dot(Mh * diff_fe));

    auto t1 = chr::high_resolution_clock::now();
    double ms = chr::duration<double,std::milli>(t1-t0).count();

    // Compare with golden
    double max_uH=0, max_uHms=0, max_uh=0;
    for (int i=0; i<NH; ++i) max_uH = std::max(max_uH, std::abs(uH(i)-guH[i]));
    for (int i=0; i<Nh; ++i) max_uHms = std::max(max_uHms, std::abs(uHms(i)-guHms[i]));
    for (int i=0; i<Nh; ++i) max_uh  = std::max(max_uh,  std::abs(uh(i)-guh[i]));

    int passed=0, failed=0;
    auto check = [&](const char* name, double v, double g, double tol=1e-10) {
        bool ok = std::abs(v-g)<tol || std::abs(v-g)/std::max(1.0,std::abs(g))<tol;
        std::cout << "  " << (ok?"[PASS]":"[FAIL]") << " " << name << ": " << v << " (gold " << g << ")\n";
        ok?++passed:++failed;
    };

    std::cout << "\nC++ total: " << ms << " ms (MATLAB: ~470 ms)\n\n";
    std::cout << "Solution comparison:\n";
    std::cout << "  max|uH diff|:   " << max_uH   << "\n";
    std::cout << "  max|uHms diff|: " << max_uHms << "\n";
    std::cout << "  max|uh diff|:   " << max_uh   << "\n\n";
    check("Energy err", errE,  gE);
    check("L2 err",     errL2, gL2);
    check("FE-L2 err",  errFE, gFE);

    std::cout << "\n========================================\n";
    std::cout << "Full LOD: " << passed << " PASS, " << failed << " FAIL\n";
    std::cout << "C++ time: " << ms << " ms  (MATLAB: ~470 ms)\n";
    std::cout << "========================================\n";
    return (failed>0) ? 1 : 0;
}

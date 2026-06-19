/// Phase F: Golden-data verification for element corrector.

#include "lod/corrector.h"
#include "lod/quasi_interp.h"
#include "lod/patches.h"
#include "fem/assemble_dg.h"
#include "mesh/refine.h"
#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace lod2d;

namespace {

struct GoldenElem {
    int k = 0;
    std::vector<int> gi;
    std::vector<int> gj;
    std::vector<double> gv;
};

const char *solver_name(CorrectorSolver solver) {
    if (solver == CorrectorSolver::Cholmod) return "cholmod";
    if (solver == CorrectorSolver::CholmodCached) return "cholmod_cached";
    return "eigen";
}

std::vector<CorrectorSolver> parse_solvers(int argc, char **argv) {
    std::string arg = "--solver=eigen";
    if (argc > 1) arg = argv[1];
    if (arg.rfind("--solver=", 0) != 0) {
        throw std::invalid_argument("usage: test_corr [--solver=eigen|cholmod|cholmod_cached|both]");
    }

    std::string value = arg.substr(std::string("--solver=").size());
    if (value == "eigen") return {CorrectorSolver::EigenLLT};
    if (value == "cholmod") return {CorrectorSolver::Cholmod};
    if (value == "cholmod_cached") return {CorrectorSolver::CholmodCached};
    if (value == "both") return {CorrectorSolver::EigenLLT, CorrectorSolver::Cholmod};
    throw std::invalid_argument("unknown solver: " + value);
}

} // namespace

int main(int argc, char **argv) {
    std::vector<CorrectorSolver> solvers;
    try {
        solvers = parse_solvers(argc, argv);
    } catch (const std::exception &e) {
        std::cerr << e.what() << "\n";
        return 2;
    }

    std::cout << "=== Phase F: Corrector Golden-data ===\n\n";

    std::ifstream f("tests/golden_corr.txt");
    if (!f) { std::cerr << "Cannot open\n"; return 1; }

    int NH, NTH, Nh, NTh_f;
    std::string tok;
    f >> tok >> NH >> tok >> NTH >> tok >> Nh >> tok >> NTh_f;

    int nAh; f >> tok >> nAh;
    std::vector<double> Ah(nAh);
    for (int i = 0; i < nAh; ++i) f >> Ah[i];

    int nn; f >> tok >> nn;
    std::vector<int> nngH(nn);
    for (int i = 0; i < nn; ++i) f >> nngH[i];

    f >> tok >> nn;
    std::vector<int> nngh(nn);
    for (int i = 0; i < nn; ++i) f >> nngh[i];

    std::vector<GoldenElem> golden;
    while (f >> tok && tok == "ELEM") {
        GoldenElem elem;
        int nnz;
        f >> elem.k >> tok >> nnz;
        elem.gi.resize(nnz);
        elem.gj.resize(nnz);
        elem.gv.resize(nnz);
        for (int n = 0; n < nnz; ++n) f >> elem.gi[n] >> elem.gj[n] >> elem.gv[n];
        f >> tok; // END_ELEM
        golden.push_back(std::move(elem));
    }

    TriMesh T0;
    T0.nodes={{0,0},{1,0},{1,1},{0,1}};
    T0.elems={{0,1,3},{1,2,3}}; T0.dirichlet={0,1,2,3};

    int Hlevel=3, hlevel=5, ell=2, d=2;
    auto c_out = refine_mesh(T0, Hlevel);
    auto f_out = refine_mesh(c_out.mesh, hlevel-Hlevel);
    const auto &coarse = c_out.mesh;
    const auto &fine   = f_out.mesh;

    std::vector<std::array<int,3>> dghidx(fine.elems.size());
    for (size_t e=0; e<fine.elems.size(); ++e) for (int i=0;i<3;++i) dghidx[e][i]=3*e+i;
    std::vector<std::array<int,3>> dgHidx(coarse.elems.size());
    for (size_t e=0; e<coarse.elems.size(); ++e) for (int i=0;i<3;++i) dgHidx[e][i]=3*e+i;

    int Nhdg=3*fine.elems.size();
    std::vector<Eigen::Triplet<double>> cg_t;
    for (size_t e=0; e<fine.elems.size(); ++e)
        for (int i=0; i<3; ++i)
            cg_t.emplace_back(3*(int)e+i, fine.elems[e][i], 1.0);
    Eigen::SparseMatrix<double> cg2dgh(Nhdg, Nh);
    cg2dgh.setFromTriplets(cg_t.begin(), cg_t.end());

    auto element_stiffness = assemble_element_stiffness(fine, Ah);
    Eigen::SparseMatrix<double> Shdg = assemble_dg_from_element_stiffness(element_stiffness);
    auto IH    = build_quasi_interp(coarse, fine, f_out.P_dg, cg2dgh, Nh, NH);
    auto patch = build_patches(coarse, ell);
    auto fine_element_children = build_fine_element_children(f_out.P_elem, NTH);

    int total_passed=0, total_failed=0;
    for (CorrectorSolver solver : solvers) {
        int passed=0, failed=0;
        std::cout << "--- Solver: " << solver_name(solver) << " ---\n";
        for (const auto &elem : golden) {
            auto CTk = compute_corrector(elem.k, patch, coarse, NH, nngH,
                f_out.P_elem, fine, Nh, nngh, dghidx, cg2dgh, Shdg,
                f_out.P_dg, dgHidx, IH, d, solver, &element_stiffness, &fine_element_children);

            double max_err=0; int matched=0;
            for (size_t n=0; n<elem.gv.size(); ++n) {
                double v = CTk.coeff(elem.gi[n], elem.gj[n]);
                double diff = std::abs(v - elem.gv[n]);
                if (diff > max_err) max_err = diff;
                if (diff < 1e-10) matched++;
            }
            bool ok = (matched == static_cast<int>(elem.gv.size()));
            std::cout << "Elem " << elem.k << ": " << (ok?"[PASS]":"[FAIL]")
                      << " " << matched << "/" << elem.gv.size()
                      << " match (max err " << max_err << ")\n";
            ok ? ++passed : ++failed;
        }
        std::cout << "Solver " << solver_name(solver) << ": "
                  << passed << " PASS, " << failed << " FAIL\n\n";
        total_passed += passed;
        total_failed += failed;
    }

    std::cout << "========================================\n";
    std::cout << "Phase F: " << total_passed << " PASS, " << total_failed << " FAIL\n";
    std::cout << "========================================\n";
    return (total_failed > 0) ? 1 : 0;
}


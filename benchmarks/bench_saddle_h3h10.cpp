#include "lod/corrector.h"
#include "lod/patches.h"
#include "lod/quasi_interp.h"
#include "fem/assemble_dg.h"
#include "mesh/refine.h"
#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace lod2d;
namespace chr = std::chrono;

namespace {

struct Options {
    int H = 3;
    int h = 10;
    int ell = 3;
    int threads = -1;
    bool skip_reference = false;
};

Options parse_options(int argc, char **argv) {
    Options opt;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.rfind("--threads=", 0) == 0) {
            std::string value = arg.substr(10);
            if (value == "auto") opt.threads = -1;
            else if (value == "env") opt.threads = 0;
            else opt.threads = std::stoi(value);
        } else if (arg == "--skip-reference") {
            opt.skip_reference = true;
        } else if (arg.rfind("--H=", 0) == 0) {
            opt.H = std::stoi(arg.substr(4));
        } else if (arg.rfind("--h=", 0) == 0) {
            opt.h = std::stoi(arg.substr(4));
        } else if (arg.rfind("--ell=", 0) == 0) {
            opt.ell = std::stoi(arg.substr(6));
        } else {
            throw std::invalid_argument("usage: bench_saddle_h3h10 [--H=N --h=N --ell=N] [--threads=auto|env|N] [--skip-reference]");
        }
    }
    if (opt.H < 0 || opt.h < opt.H || opt.ell < 0) throw std::invalid_argument("require 0 <= H <= h and ell >= 0");
    return opt;
}

#ifdef _OPENMP
void apply_threads(const Options &opt) {
    if (opt.threads > 0) omp_set_num_threads(opt.threads);
}
int thread_count() { return omp_get_max_threads(); }
#else
void apply_threads(const Options &) {}
int thread_count() { return 1; }
#endif

TriMesh unit_square_red_seed() {
    TriMesh mesh;
    mesh.nodes = {{0,0},{1,0},{1,1},{0,1}};
    mesh.elems = {{0,1,3},{1,2,3}};
    mesh.dirichlet = {0,1,2,3};
    return mesh;
}

std::vector<int> node_incidence(const TriMesh &mesh) {
    std::vector<int> nng(mesh.nodes.size(), 0);
    for (const auto &elem : mesh.elems) for (int v : elem) ++nng[v];
    for (int v : mesh.dirichlet) if (v >= 0 && v < static_cast<int>(nng.size())) nng[v] = 0;
    return nng;
}

Eigen::SparseMatrix<double> build_cg_to_dg_local(const TriMesh &mesh) {
    std::vector<Eigen::Triplet<double>> t;
    t.reserve(3 * mesh.elems.size());
    for (int e = 0; e < static_cast<int>(mesh.elems.size()); ++e)
        for (int i = 0; i < 3; ++i)
            t.emplace_back(3 * e + i, mesh.elems[e][i], 1.0);
    Eigen::SparseMatrix<double> M(3 * static_cast<int>(mesh.elems.size()), static_cast<int>(mesh.nodes.size()));
    M.setFromTriplets(t.begin(), t.end());
    return M;
}

struct SolveResult {
    Eigen::VectorXd uH;
    Eigen::VectorXd uHms;
    double coarse_ms = 0.0;
};

SolveResult solve_lod(
    const TriMesh &coarse,
    const TriMesh &fine,
    const Eigen::SparseMatrix<double> &P_node,
    const Eigen::SparseMatrix<double> &G,
    const Eigen::SparseMatrix<double> &Sh,
    const Eigen::SparseMatrix<double> &Mh) {
    auto t0 = chr::high_resolution_clock::now();
    const int NH = static_cast<int>(coarse.nodes.size());
    const int Nh = static_cast<int>(fine.nodes.size());
    std::vector<int> dofH;
    std::vector<int> dofH_map(NH, -1);
    std::vector<char> is_dirH(NH, false);
    for (int v : coarse.dirichlet) is_dirH[v] = true;
    for (int i = 0; i < NH; ++i) if (!is_dirH[i]) {
        dofH_map[i] = static_cast<int>(dofH.size());
        dofH.push_back(i);
    }
    std::vector<Eigen::Triplet<double>> g0_t;
    for (int col = 0; col < G.outerSize(); ++col)
        for (Eigen::SparseMatrix<double>::InnerIterator it(G, col); it; ++it)
            if (dofH_map[it.col()] >= 0) g0_t.emplace_back(it.row(), dofH_map[it.col()], it.value());
    Eigen::SparseMatrix<double> G0(Nh, static_cast<int>(dofH.size()));
    G0.setFromTriplets(g0_t.begin(), g0_t.end());
    Eigen::SparseMatrix<double> SHLOD0 = G0.transpose() * Sh * G0;
    Eigen::VectorXd f_coarse = Eigen::VectorXd::Ones(NH);
    Eigen::VectorXd rhs = G0.transpose() * (Mh * (P_node * f_coarse));
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt(SHLOD0);
    Eigen::VectorXd uf = llt.solve(rhs);
    SolveResult result;
    result.uH = Eigen::VectorXd::Zero(NH);
    for (int j = 0; j < static_cast<int>(dofH.size()); ++j) result.uH(dofH[j]) = uf(j);
    result.uHms = G * result.uH;
    auto t1 = chr::high_resolution_clock::now();
    result.coarse_ms = chr::duration<double, std::milli>(t1 - t0).count();
    return result;
}

struct Errors {
    double energy = 0.0;
    double l2 = 0.0;
    double fe_l2 = 0.0;
};

Eigen::VectorXd reference_solution(
    const TriMesh &fine,
    const Eigen::SparseMatrix<double> &Sh,
    const Eigen::SparseMatrix<double> &Mh) {
    const int Nh = static_cast<int>(fine.nodes.size());
    std::vector<int> dofh;
    std::vector<int> dofh_map(Nh, -1);
    std::vector<char> is_dir(Nh, false);
    for (int v : fine.dirichlet) is_dir[v] = true;
    for (int i = 0; i < Nh; ++i) if (!is_dir[i]) {
        dofh_map[i] = static_cast<int>(dofh.size());
        dofh.push_back(i);
    }
    std::vector<Eigen::Triplet<double>> sh_t;
    for (int col = 0; col < Sh.outerSize(); ++col)
        for (Eigen::SparseMatrix<double>::InnerIterator it(Sh, col); it; ++it) {
            int r = dofh_map[it.row()];
            int c = dofh_map[it.col()];
            if (r >= 0 && c >= 0) sh_t.emplace_back(r, c, it.value());
        }
    Eigen::SparseMatrix<double> Sh_free(static_cast<int>(dofh.size()), static_cast<int>(dofh.size()));
    Sh_free.setFromTriplets(sh_t.begin(), sh_t.end());
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(static_cast<int>(dofh.size()));
    for (int col = 0; col < Mh.outerSize(); ++col)
        for (Eigen::SparseMatrix<double>::InnerIterator it(Mh, col); it; ++it)
            if (dofh_map[it.row()] >= 0) rhs(dofh_map[it.row()]) += it.value();
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt(Sh_free);
    Eigen::VectorXd uf = llt.solve(rhs);
    Eigen::VectorXd uh = Eigen::VectorXd::Zero(Nh);
    for (int j = 0; j < static_cast<int>(dofh.size()); ++j) uh(dofh[j]) = uf(j);
    return uh;
}

Errors compute_errors(
    const Eigen::VectorXd &uh,
    const Eigen::VectorXd &uHms,
    const Eigen::VectorXd &uH_fe,
    const Eigen::SparseMatrix<double> &Sh,
    const Eigen::SparseMatrix<double> &Mh) {
    Eigen::VectorXd diff = uh - uHms;
    Eigen::VectorXd diff_fe = uh - uH_fe;
    return {std::sqrt(diff.dot(Sh * diff)), std::sqrt(diff.dot(Mh * diff)), std::sqrt(diff_fe.dot(Mh * diff_fe))};
}

uint64_t entry_key(int row, int col) {
    return (static_cast<uint64_t>(static_cast<uint32_t>(row)) << 32) | static_cast<uint32_t>(col);
}

double max_corrector_diff(const std::vector<CorrectorEntries> &a, const std::vector<CorrectorEntries> &b) {
    double max_diff = 0.0;
    for (size_t k = 0; k < a.size(); ++k) {
        std::unordered_map<uint64_t, double> values;
        values.reserve(a[k].size() * 2 + 8);
        for (const auto &e : a[k]) values[entry_key(e.row, e.col)] += e.value;
        for (const auto &e : b[k]) {
            uint64_t key = entry_key(e.row, e.col);
            auto it = values.find(key);
            double av = 0.0;
            if (it != values.end()) {
                av = it->second;
                values.erase(it);
            }
            max_diff = std::max(max_diff, std::abs(av - e.value));
        }
        for (const auto &[_, value] : values) max_diff = std::max(max_diff, std::abs(value));
    }
    return max_diff;
}

const char *solver_name(CorrectorSolver solver) {
    if (solver == CorrectorSolver::SaddleGmres) return "saddle_gmres";
    if (solver == CorrectorSolver::Cholmod) return "cholmod";
    if (solver == CorrectorSolver::CholmodCached) return "cholmod_cached";
    return "eigen_schur";
}

std::vector<CorrectorEntries> run_correctors(
    CorrectorSolver solver,
    const Eigen::SparseMatrix<double> &patch,
    const TriMesh &coarse,
    int NH,
    const std::vector<int> &nngH,
    const TriMesh &fine,
    int Nh,
    const std::vector<int> &nngh,
    const Eigen::SparseMatrix<double> &P_dg,
    const std::vector<std::array<int, 3>> &dgHidx,
    const ElementStiffnessBlocks &element_stiffness,
    const FineElementChildren &fine_element_children,
    const InterpolationRows &interpolation_rows,
    double &ms) {
    Eigen::SparseMatrix<double> unused;
    std::vector<std::array<int, 3>> dghidx;
    auto t0 = chr::high_resolution_clock::now();
    auto CT = compute_all_correctors(
        patch, coarse, NH, nngH, unused, fine, Nh, nngh,
        dghidx, unused, unused, P_dg, dgHidx, unused, 2, solver,
        &element_stiffness, &fine_element_children, &interpolation_rows);
    auto t1 = chr::high_resolution_clock::now();
    ms = chr::duration<double, std::milli>(t1 - t0).count();
    return CT;
}

} // namespace

int main(int argc, char **argv) {
    Options opt = parse_options(argc, argv);
    apply_threads(opt);
    auto total0 = chr::high_resolution_clock::now();

    std::cout << "=== Saddle corrector experiment: red mesh H=" << opt.H
              << " h=" << opt.h << " ell=" << opt.ell
              << " threads=" << thread_count()
              << (opt.skip_reference ? " skip_reference=1" : "") << " ===\n";

    TriMesh T0 = unit_square_red_seed();
    auto mesh0 = chr::high_resolution_clock::now();
    auto c_out = refine_mesh_red(T0, opt.H);
    auto f_out = refine_mesh_red(c_out.mesh, opt.h - opt.H);
    auto mesh1 = chr::high_resolution_clock::now();
    const auto &coarse = c_out.mesh;
    const auto &fine = f_out.mesh;
    const int NH = static_cast<int>(coarse.nodes.size());
    const int NTH = static_cast<int>(coarse.elems.size());
    const int Nh = static_cast<int>(fine.nodes.size());
    const int NTh = static_cast<int>(fine.elems.size());
    std::cout << "Coarse: " << NH << " v " << NTH << " t  Fine: " << Nh << " v " << NTh << " t\n";
    std::cout << "Mesh(red): " << chr::duration<double, std::milli>(mesh1 - mesh0).count() << " ms\n";

    auto op0 = chr::high_resolution_clock::now();
    std::vector<double> Ah(NTh, 1.0);
    auto element_stiffness = assemble_element_stiffness(fine, Ah);
    auto nngH = node_incidence(coarse);
    auto nngh = node_incidence(fine);
    std::vector<std::array<int, 3>> dgHidx(NTH);
    for (int e = 0; e < NTH; ++e) for (int i = 0; i < 3; ++i) dgHidx[e][i] = 3 * e + i;
    Eigen::SparseMatrix<double> cg2dgh = build_cg_to_dg_local(fine);
    Eigen::SparseMatrix<double> IH = build_quasi_interp(coarse, fine, f_out.P_dg, cg2dgh, Nh, NH);
    auto interpolation_rows = build_interpolation_rows(IH, NH);
    auto patch = build_patches(coarse, opt.ell);
    auto fine_element_children = build_fine_element_children(f_out.P_elem, NTH);
    auto areas = compute_area(fine);
    Eigen::SparseMatrix<double> Sh = assemble_cg_from_element_stiffness(fine, element_stiffness);
    Eigen::SparseMatrix<double> Mh = assemble_cg_mass(fine, areas);
    auto op1 = chr::high_resolution_clock::now();
    std::cout << "Operators: " << chr::duration<double, std::milli>(op1 - op0).count()
              << " ms  patch_nnz=" << patch.nonZeros() << "\n";

    double eigen_ms = 0.0;
    auto CT_eigen = run_correctors(CorrectorSolver::EigenLLT, patch, coarse, NH, nngH, fine, Nh, nngh,
        f_out.P_dg, dgHidx, element_stiffness, fine_element_children, interpolation_rows, eigen_ms);
    std::cout << "Correctors (" << solver_name(CorrectorSolver::EigenLLT) << "): " << eigen_ms << " ms\n";

    double saddle_ms = 0.0;
    auto CT_saddle = run_correctors(CorrectorSolver::SaddleGmres, patch, coarse, NH, nngH, fine, Nh, nngh,
        f_out.P_dg, dgHidx, element_stiffness, fine_element_children, interpolation_rows, saddle_ms);
    std::cout << "Correctors (" << solver_name(CorrectorSolver::SaddleGmres) << "): " << saddle_ms << " ms\n";

    const double max_ct_diff = max_corrector_diff(CT_eigen, CT_saddle);
    std::cout << "max |CT_eigen - CT_saddle_gmres| = " << max_ct_diff << "\n";

    auto G_eigen = build_multiscale_basis(f_out.P_node, coarse, Nh, CT_eigen);
    auto G_saddle = build_multiscale_basis(f_out.P_node, coarse, Nh, CT_saddle);
    auto sol_eigen = solve_lod(coarse, fine, f_out.P_node, G_eigen, Sh, Mh);
    auto sol_saddle = solve_lod(coarse, fine, f_out.P_node, G_saddle, Sh, Mh);
    std::cout << "Coarse solve eigen/saddle: " << sol_eigen.coarse_ms << " / " << sol_saddle.coarse_ms << " ms\n";
    std::cout << "max |uHms_eigen - uHms_saddle| = " << (sol_eigen.uHms - sol_saddle.uHms).cwiseAbs().maxCoeff() << "\n";

    if (!opt.skip_reference) {
        auto ref0 = chr::high_resolution_clock::now();
        Eigen::VectorXd uh = reference_solution(fine, Sh, Mh);
        auto ref1 = chr::high_resolution_clock::now();
        auto err_eigen = compute_errors(uh, sol_eigen.uHms, f_out.P_node * sol_eigen.uH, Sh, Mh);
        auto err_saddle = compute_errors(uh, sol_saddle.uHms, f_out.P_node * sol_saddle.uH, Sh, Mh);
        std::cout << "Reference: " << chr::duration<double, std::milli>(ref1 - ref0).count() << " ms\n";
        std::cout << "Errors eigen:  E=" << err_eigen.energy << " L2=" << err_eigen.l2 << " FE=" << err_eigen.fe_l2 << "\n";
        std::cout << "Errors saddle: E=" << err_saddle.energy << " L2=" << err_saddle.l2 << " FE=" << err_saddle.fe_l2 << "\n";
    }

    auto total1 = chr::high_resolution_clock::now();
    std::cout << "Total wall: " << chr::duration<double, std::milli>(total1 - total0).count() << " ms\n";
    return 0;
}

/// Benchmark repeated LOD solves for the same A/mesh/H/h/ell and different RHS.
#include "lod/corrector.h"
#include "lod/quasi_interp.h"
#include "lod/patches.h"
#include "fem/assemble_dg.h"
#include "mesh/refine.h"
#include "solver/lod_reuse.h"
#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace lod2d;
namespace chr = std::chrono;

namespace {
struct Options {
    CorrectorSolver solver = CorrectorSolver::Cholmod;
    int rhs_count = 8;
    int threads = -1;
};

const char *solver_name(CorrectorSolver s) {
    if (s == CorrectorSolver::Cholmod) return "cholmod";
    if (s == CorrectorSolver::CholmodCached) return "cholmod_cached";
    return "eigen";
}

Options parse_options(int argc, char **argv) {
    Options opt;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.rfind("--solver=", 0) == 0) {
            std::string value = arg.substr(std::string("--solver=").size());
            if (value == "eigen") opt.solver = CorrectorSolver::EigenLLT;
            else if (value == "cholmod") opt.solver = CorrectorSolver::Cholmod;
            else if (value == "cholmod_cached") opt.solver = CorrectorSolver::CholmodCached;
            else if (value == "auto") opt.solver = CorrectorSolver::Cholmod;
            else throw std::invalid_argument("unknown solver: " + value);
        } else if (arg.rfind("--rhs=", 0) == 0) {
            opt.rhs_count = std::stoi(arg.substr(std::string("--rhs=").size()));
        } else if (arg.rfind("--threads=", 0) == 0) {
            std::string value = arg.substr(std::string("--threads=").size());
            if (value == "auto") opt.threads = -1;
            else if (value == "env") opt.threads = 0;
            else opt.threads = std::stoi(value);
        } else {
            throw std::invalid_argument("usage: bench_reuse_rhs [--solver=eigen|cholmod|cholmod_cached|auto] [--rhs=N] [--threads=auto|env|N]");
        }
    }
    if (opt.rhs_count <= 0) throw std::invalid_argument("--rhs must be positive");
    return opt;
}

#ifdef _OPENMP
void apply_thread_option(const Options &opt) {
    if (opt.threads > 0) {
        omp_set_num_threads(opt.threads);
        return;
    }
    if (opt.threads == 0) return;
    if ((opt.solver == CorrectorSolver::Cholmod || opt.solver == CorrectorSolver::CholmodCached) && omp_get_max_threads() > 8)
        omp_set_num_threads(8);
}
#else
void apply_thread_option(const Options &) {}
#endif

Eigen::VectorXd make_rhs(const TriMesh &coarse, int m) {
    Eigen::VectorXd f(coarse.nodes.size());
    const double a = static_cast<double>(m + 1);
    const double b = static_cast<double>(m + 2);
    for (size_t i = 0; i < coarse.nodes.size(); ++i) {
        const double x = coarse.nodes[i][0];
        const double y = coarse.nodes[i][1];
        f(static_cast<int>(i)) = 1.0 + 0.2 * std::sin(a * M_PI * x) * std::cos(b * M_PI * y);
    }
    return f;
}

} // namespace

int main(int argc, char **argv) {
    Options opt;
    try {
        opt = parse_options(argc, argv);
    } catch (const std::exception &e) {
        std::cerr << e.what() << "\n";
        return 2;
    }
    apply_thread_option(opt);

    int H=4, h=10, ell=2, d=2;
    std::cout << "=== Reusable LOD RHS benchmark H=" << H << " h=" << h
              << " solver=" << solver_name(opt.solver)
              << " rhs=" << opt.rhs_count << " ===\n";

    std::ifstream af("benchmarks/data_H4h10.txt");
    int nAh; af >> nAh; std::vector<double> Ah(nAh);
    for (int i=0;i<nAh;++i) af >> Ah[i];

    auto t0 = chr::high_resolution_clock::now();
    TriMesh T0;
    T0.nodes={{0,0},{1,0},{1,1},{0,1}}; T0.elems={{0,1,3},{1,2,3}}; T0.dirichlet={0,1,2,3};
    auto c_out=refine_mesh(T0,H); auto f_out=refine_mesh(c_out.mesh,h-H);
    const auto &coarse=c_out.mesh, &fine=f_out.mesh;
    int NH=coarse.nodes.size(), NTH=coarse.elems.size();
    int Nh=fine.nodes.size(), NTh_f=fine.elems.size();

    std::vector<std::array<int,3>> dghidx, dgHidx(NTH);
    for(int e=0;e<NTH;++e)for(int i=0;i<3;++i)dgHidx[e][i]=3*e+i;

    std::vector<int> nngH(NH),nngh(Nh);
    for(auto&t:coarse.elems)for(int v:t)nngH[v]++; for(int v:coarse.dirichlet)nngH[v]=0;
    for(auto&t:fine.elems)for(int v:t)nngh[v]++; for(int v:fine.dirichlet)nngh[v]=0;

    auto element_stiffness=assemble_element_stiffness(fine,Ah);
    std::vector<double>().swap(Ah);
    Eigen::SparseMatrix<double> IH;
    {
        std::vector<Eigen::Triplet<double>> cg_t;
        cg_t.reserve(3 * NTh_f);
        for(int e=0;e<NTh_f;++e)for(int i=0;i<3;++i)cg_t.emplace_back(3*e+i,fine.elems[e][i],1.0);
        Eigen::SparseMatrix<double> cg2dgh(3*NTh_f,Nh); cg2dgh.setFromTriplets(cg_t.begin(),cg_t.end());
        IH=build_quasi_interp(coarse,fine,f_out.P_dg,cg2dgh,Nh,NH);
    }
    auto interpolation_rows=build_interpolation_rows(IH,NH);
    auto patch=build_patches(coarse,ell);
    auto fine_element_children=build_fine_element_children(f_out.P_elem,NTH);
    Eigen::SparseMatrix<double> empty_elem;
    f_out.P_elem.swap(empty_elem);
    auto areas=compute_area(fine);
    Eigen::SparseMatrix<double> Sh=assemble_cg_from_element_stiffness(fine,element_stiffness);
    Eigen::SparseMatrix<double> Mh=assemble_cg_mass(fine,areas);
    std::vector<double>().swap(areas);
    Eigen::SparseMatrix<double> empty_ih;
    IH.swap(empty_ih);

    auto tc0=chr::high_resolution_clock::now();
    std::vector<CorrectorEntries> CT(NTH);
    Eigen::SparseMatrix<double> unused_sparse;
    #pragma omp parallel for schedule(dynamic)
    for(int k=0;k<NTH;++k)
        CT[k]=compute_corrector_entries(k,patch,coarse,NH,nngH,unused_sparse,fine,Nh,nngh,dghidx,unused_sparse,unused_sparse,f_out.P_dg,dgHidx,unused_sparse,d,opt.solver,&element_stiffness,&fine_element_children,&interpolation_rows);
    auto tc1=chr::high_resolution_clock::now();

    ElementStiffnessBlocks().swap(element_stiffness);
    FineElementChildren().swap(fine_element_children);
    InterpolationRows().swap(interpolation_rows);
    Eigen::SparseMatrix<double> empty_pdg;
    f_out.P_dg.swap(empty_pdg);

    size_t corrector_nnz = 0;
    for (const auto &entries : CT) corrector_nnz += entries.size();
    std::vector<Eigen::Triplet<double>> g_t;
    g_t.reserve(static_cast<size_t>(f_out.P_node.nonZeros()) + corrector_nnz);
    for(int c=0;c<f_out.P_node.outerSize();++c)
        for(Eigen::SparseMatrix<double>::InnerIterator it(f_out.P_node,c);it;++it)
            g_t.emplace_back(it.row(),it.col(),it.value());
    for(int k=0;k<NTH;++k)
        for(const auto &entry:CT[k])
            g_t.emplace_back(entry.row,coarse.elems[k][entry.col],-entry.value);
    Eigen::SparseMatrix<double> G(Nh,NH); G.setFromTriplets(g_t.begin(),g_t.end());
    std::vector<CorrectorEntries>().swap(CT);

    LodReusableSystem lod_system(G, Sh, Mh, f_out.P_node, NH, coarse.dirichlet);
    auto t_setup = chr::high_resolution_clock::now();

    double checksum = 0.0;
    auto t_rhs0 = chr::high_resolution_clock::now();
    for (int r = 0; r < opt.rhs_count; ++r) {
        Eigen::VectorXd f = make_rhs(coarse, r);
        LodReuseSolution sol = lod_system.solve_from_coarse_values(f);
        checksum += sol.uHms.squaredNorm();
    }
    auto t_rhs1 = chr::high_resolution_clock::now();

    const double setup_ms = chr::duration<double,std::milli>(t_setup-t0).count();
    const double corr_ms = chr::duration<double,std::milli>(tc1-tc0).count();
    const double rhs_ms = chr::duration<double,std::milli>(t_rhs1-t_rhs0).count();

    std::cout << "Coarse:" << NH << "v " << NTH << "t  Fine:" << Nh << "v " << NTh_f << "t\n";
    std::cout << "Reusable setup: " << setup_ms << " ms\n";
    std::cout << "  Correctors: " << corr_ms << " ms";
#ifdef _OPENMP
    std::cout << " (" << omp_get_max_threads() << " threads)";
#endif
    std::cout << "\n";
    std::cout << "Repeated RHS total: " << rhs_ms << " ms\n";
    std::cout << "Repeated RHS avg: " << rhs_ms / opt.rhs_count << " ms\n";
    std::cout << "Checksum: " << checksum << "\n";
    return 0;
}

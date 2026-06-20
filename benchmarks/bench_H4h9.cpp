/// H=4 benchmark; this local copy is currently configured for h=10.
#include "lod/corrector.h"
#include "lod/quasi_interp.h"
#include "lod/patches.h"
#include "fem/assemble_dg.h"
#include "mesh/refine.h"
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
    CorrectorSolver solver = CorrectorSolver::EigenLLT;
    bool skip_reference = false;
    int threads = -1;  // -1: auto, 0: keep environment, >0: explicit
};

const char *solver_name(CorrectorSolver s) {
    if (s == CorrectorSolver::Cholmod) return "cholmod";
    if (s == CorrectorSolver::CholmodCached) return "cholmod_cached";
    return "eigen";
}

Options parse_options(int argc, char **argv, int h) {
    Options opt;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.rfind("--solver=", 0) == 0) {
            std::string value = arg.substr(std::string("--solver=").size());
            if (value == "eigen") opt.solver = CorrectorSolver::EigenLLT;
            else if (value == "cholmod") opt.solver = CorrectorSolver::Cholmod;
            else if (value == "cholmod_cached") opt.solver = CorrectorSolver::CholmodCached;
            else if (value == "auto") opt.solver = (h >= 10) ? CorrectorSolver::Cholmod : CorrectorSolver::EigenLLT;
            else throw std::invalid_argument("unknown solver: " + value);
        } else if (arg.rfind("--threads=", 0) == 0) {
            std::string value = arg.substr(std::string("--threads=").size());
            if (value == "auto") opt.threads = -1;
            else if (value == "env") opt.threads = 0;
            else opt.threads = std::stoi(value);
        } else if (arg == "--skip-reference") {
            opt.skip_reference = true;
        } else {
            throw std::invalid_argument("usage: bench_H4h9 [--solver=eigen|cholmod|cholmod_cached|auto] [--threads=auto|env|N] [--skip-reference]");
        }
    }
    return opt;
}
#ifdef _OPENMP
void apply_thread_option(const Options &opt, int h) {
    if (opt.threads > 0) {
        omp_set_num_threads(opt.threads);
        return;
    }
    if (opt.threads == 0) return;
    if (h >= 10 && (opt.solver == CorrectorSolver::Cholmod || opt.solver == CorrectorSolver::CholmodCached) && omp_get_max_threads() > 8)
        omp_set_num_threads(8);
}
#else
void apply_thread_option(const Options &, int) {}
#endif

} // namespace

int main(int argc, char **argv) {
    int H=4, h=10, ell=2, d=2;
    Options opt = parse_options(argc, argv, h);
    apply_thread_option(opt, h);
    CorrectorSolver solver = opt.solver;
    std::cout << "=== C++ LOD H=" << H << " h=" << h
              << " solver=" << solver_name(solver)
              << (opt.skip_reference ? " skip_reference=1" : "") << " ===\n";

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
    std::cout<<"Coarse:"<<NH<<"v "<<NTH<<"t  Fine:"<<Nh<<"v "<<NTh_f<<"t\n";

    std::vector<std::array<int,3>> dghidx, dgHidx(NTH);
    for(int e=0;e<NTH;++e)for(int i=0;i<3;++i)dgHidx[e][i]=3*e+i;

    std::vector<int> nngH(NH),nngh(Nh);
    for(auto&t:coarse.elems)for(int v:t)nngH[v]++; for(int v:coarse.dirichlet)nngH[v]=0;
    for(auto&t:fine.elems)for(int v:t)nngh[v]++; for(int v:fine.dirichlet)nngh[v]=0;

    auto tp0=chr::high_resolution_clock::now();
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
    auto tp1=chr::high_resolution_clock::now();
    std::cout<<"Setup: "<<chr::duration<double,std::milli>(tp1-tp0).count()<<" ms\n";

    auto tc0=chr::high_resolution_clock::now();
    Eigen::SparseMatrix<double> unused_sparse;
    std::vector<CorrectorEntries> CT = compute_all_correctors(
        patch, coarse, NH, nngH, unused_sparse, fine, Nh, nngh,
        dghidx, unused_sparse, unused_sparse, f_out.P_dg, dgHidx,
        unused_sparse, d, solver, &element_stiffness,
        &fine_element_children, &interpolation_rows);
    auto tc1=chr::high_resolution_clock::now();
    double tc=chr::duration<double,std::milli>(tc1-tc0).count();
    std::cout<<"Correctors ("<<solver_name(solver)<<"): "<<tc<<" ms";
#ifdef _OPENMP
    std::cout<<" ("<<omp_get_max_threads()<<" threads)";
#endif
    std::cout<<"\n";
    ElementStiffnessBlocks().swap(element_stiffness);
    FineElementChildren().swap(fine_element_children);
    InterpolationRows().swap(interpolation_rows);
    Eigen::SparseMatrix<double> empty_pdg;
    f_out.P_dg.swap(empty_pdg);

    Eigen::SparseMatrix<double> G = build_multiscale_basis(f_out.P_node, coarse, Nh, CT);
    std::vector<CorrectorEntries>().swap(CT);

    std::vector<int> dofH; std::vector<int> dofH_map(NH,-1);
    std::vector<char> is_dirH(NH,false);
    for(int dv:coarse.dirichlet)is_dirH[dv]=true;
    for(int i=0;i<NH;++i)if(!is_dirH[i]){dofH_map[i]=static_cast<int>(dofH.size());dofH.push_back(i);}
    std::vector<Eigen::Triplet<double>> g0_t;
    for(int k2=0;k2<G.outerSize();++k2)
        for(Eigen::SparseMatrix<double>::InnerIterator it(G,k2);it;++it)
            if(dofH_map[it.col()]>=0)g0_t.emplace_back(it.row(),dofH_map[it.col()],it.value());
    Eigen::SparseMatrix<double> G0(Nh,dofH.size()); G0.setFromTriplets(g0_t.begin(),g0_t.end());
    Eigen::SparseMatrix<double> SHLOD0=G0.transpose()*Sh*G0;
    Eigen::VectorXd f_coarse=Eigen::VectorXd::Ones(NH);
    Eigen::VectorXd rhs=G0.transpose()*(Mh*(f_out.P_node*f_coarse));
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llts(SHLOD0);
    Eigen::VectorXd uH=Eigen::VectorXd::Zero(NH), uf=llts.solve(rhs);
    for(size_t j=0;j<dofH.size();++j)uH(dofH[j])=uf(j);
    Eigen::VectorXd uHms=G*uH;

    if (opt.skip_reference) {
        auto t1=chr::high_resolution_clock::now();
        double ms=chr::duration<double,std::milli>(t1-t0).count();
        std::cout<<"Reference: skipped\n";
        std::cout<<"C++ LOD-only: "<<ms<<" ms\n";
        return 0;
    }

    std::vector<int> dofh; std::vector<int> dofh_map(Nh,-1);
    std::vector<char> is_dirh(Nh,false);
    for(int dv:fine.dirichlet)is_dirh[dv]=true;
    for(int i=0;i<Nh;++i)if(!is_dirh[i]){dofh_map[i]=static_cast<int>(dofh.size());dofh.push_back(i);}
    std::vector<Eigen::Triplet<double>> sh_t;
    for(int k2=0;k2<Sh.outerSize();++k2)
        for(Eigen::SparseMatrix<double>::InnerIterator it(Sh,k2);it;++it)
            if(dofh_map[it.row()]>=0&&dofh_map[it.col()]>=0)
                sh_t.emplace_back(dofh_map[it.row()],dofh_map[it.col()],it.value());
    Eigen::SparseMatrix<double> Sh_free(dofh.size(),dofh.size()); Sh_free.setFromTriplets(sh_t.begin(),sh_t.end());
    Eigen::VectorXd rhs_ref=Eigen::VectorXd::Zero(dofh.size());
    for(int k2=0;k2<Mh.outerSize();++k2)
        for(Eigen::SparseMatrix<double>::InnerIterator it(Mh,k2);it;++it)
            if(dofh_map[it.row()]>=0)rhs_ref(dofh_map[it.row()])+=it.value();
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt_ref(Sh_free);
    Eigen::VectorXd uh=Eigen::VectorXd::Zero(Nh), uf2=llt_ref.solve(rhs_ref);
    for(size_t j=0;j<dofh.size();++j)uh(dofh[j])=uf2(j);

    Eigen::VectorXd diff=uh-uHms, diff_fe=uh-f_out.P_node*uH;
    double errE=std::sqrt(diff.dot(Sh*diff)), errL2=std::sqrt(diff.dot(Mh*diff)), errFE=std::sqrt(diff_fe.dot(Mh*diff_fe));
    auto t1=chr::high_resolution_clock::now();
    double ms=chr::duration<double,std::milli>(t1-t0).count();

    std::cout<<"\nErrors: E="<<errE<<" L2="<<errL2<<" FE="<<errFE<<"\n";
    std::cout<<"C++: "<<ms<<" ms\n";
    return 0;
}

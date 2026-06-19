/// Phase timing breakdown — uses optimized C_ell path (same as bench_H4h8)
#include "lod/corrector.h"
#include "lod/quasi_interp.h"
#include "lod/patches.h"
#include "fem/assemble_dg.h"
#include "mesh/refine.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <unordered_map>
#include <vector>
#include <Eigen/Dense>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace lod2d;
namespace chr = std::chrono;

int main() {
    int H=4, h=10, ell=2, d=2;
    std::ifstream af("benchmarks/data_H4h10.txt");
    int nAh; af >> nAh; std::vector<double> Ah(nAh);
    for (int i=0;i<nAh;++i) af >> Ah[i];

    TriMesh T0;
    T0.nodes={{0,0},{1,0},{1,1},{0,1}}; T0.elems={{0,1,3},{1,2,3}}; T0.dirichlet={0,1,2,3};

    std::cout << "=== C++ LOD Phase Timing (H=" << H << ",h=" << h << ") ===\n\n";

    // 1. Mesh refinement
    auto t1a = chr::high_resolution_clock::now();
    auto c_out = refine_mesh(T0, H);
    auto f_out = refine_mesh(c_out.mesh, h-H);
    auto t1b = chr::high_resolution_clock::now();
    double t1 = chr::duration<double,std::milli>(t1b-t1a).count();
    const auto &coarse=c_out.mesh, &fine=f_out.mesh;
    int NH=coarse.nodes.size(), NTH=coarse.elems.size();
    int Nh=fine.nodes.size(), NTh_f=fine.elems.size();
    std::cout << "1. Mesh: " << t1 << " ms → Coarse:" << NH << "v " << NTH << "t  Fine:" << Nh << "v " << NTh_f << "t\n";

    std::vector<std::array<int,3>> dghidx(NTh_f), dgHidx(NTH);
    for(int e=0;e<NTh_f;++e)for(int i=0;i<3;++i)dghidx[e][i]=3*e+i;
    for(int e=0;e<NTH;++e)for(int i=0;i<3;++i)dgHidx[e][i]=3*e+i;

    std::vector<int> nngH(NH),nngh(Nh);
    for(auto&t:coarse.elems)for(int v:t)nngH[v]++; for(int v:coarse.dirichlet)nngH[v]=0;
    for(auto&t:fine.elems)for(int v:t)nngh[v]++; for(int v:fine.dirichlet)nngh[v]=0;

    // 2. Operators
    auto t2a = chr::high_resolution_clock::now();
    std::vector<Eigen::Triplet<double>> cg_t;
    for(int e=0;e<NTh_f;++e)for(int i=0;i<3;++i)cg_t.emplace_back(3*e+i,fine.elems[e][i],1.0);
    Eigen::SparseMatrix<double> cg2dgh(3*NTh_f,Nh); cg2dgh.setFromTriplets(cg_t.begin(),cg_t.end());
    auto element_stiffness=assemble_element_stiffness(fine,Ah);
    Eigen::SparseMatrix<double> Shdg=assemble_dg_from_element_stiffness(element_stiffness);
    auto IH=build_quasi_interp(coarse,fine,f_out.P_dg,cg2dgh,Nh,NH);
    auto interpolation_rows=build_interpolation_rows(IH,NH);
    auto patch=build_patches(coarse,ell);
    auto fine_element_children=build_fine_element_children(f_out.P_elem,NTH);
    auto areas=compute_area(fine);
    Eigen::SparseMatrix<double> Sh=assemble_cg_from_element_stiffness(fine,element_stiffness);
    Eigen::SparseMatrix<double> Mh=assemble_cg_mass(fine,areas);
    auto t2b = chr::high_resolution_clock::now();
    double t2 = chr::duration<double,std::milli>(t2b-t2a).count();
    std::cout << "2. Operators: " << t2 << " ms (DG:" << Shdg.nonZeros() << "nz IH:" << IH.nonZeros() << "nz patch:" << patch.nonZeros() << "nz)\n";

    // 3. Correctors
    auto t3a = chr::high_resolution_clock::now();
    std::vector<CorrectorEntries> CT(NTH);
    #pragma omp parallel for schedule(dynamic)
    for(int k=0;k<NTH;++k)
        CT[k]=compute_corrector_entries(k,patch,coarse,NH,nngH,f_out.P_elem,fine,Nh,nngh,dghidx,cg2dgh,Shdg,f_out.P_dg,dgHidx,IH,d,CorrectorSolver::EigenLLT,&element_stiffness,&fine_element_children,&interpolation_rows);
    auto t3b = chr::high_resolution_clock::now();
    double t3 = chr::duration<double,std::milli>(t3b-t3a).count();
#ifdef _OPENMP
    std::cout << "3. Correctors: " << t3 << " ms (" << omp_get_max_threads() << " threads, " << NTH/t3*1000 << " elem/s)\n";
#else
    std::cout << "3. Correctors: " << t3 << " ms (serial, " << NTH/t3*1000 << " elem/s)\n";
#endif

    // 4. G = P_node - C_ell (direct assembly, no cell_mat)
    auto t4a = chr::high_resolution_clock::now();
    std::vector<Eigen::Triplet<double>> g_t;
    g_t.reserve(f_out.P_node.nonZeros());
    for(int c=0;c<f_out.P_node.outerSize();++c)
        for(Eigen::SparseMatrix<double>::InnerIterator it(f_out.P_node,c);it;++it)
            g_t.emplace_back(it.row(),it.col(),it.value());
    for(int k=0;k<NTH;++k)
        for(const auto &entry:CT[k])
            g_t.emplace_back(entry.row,coarse.elems[k][entry.col],-entry.value);
    Eigen::SparseMatrix<double> G(Nh,NH); G.setFromTriplets(g_t.begin(),g_t.end());
    auto t4b = chr::high_resolution_clock::now();
    double t4 = chr::duration<double,std::milli>(t4b-t4a).count();
    std::cout << "4. C_ell + G: " << t4 << " ms\n";

    // 5. Coarse LOD solve
    auto t5a = chr::high_resolution_clock::now();
    std::vector<int> dofH; std::unordered_map<int,int> dofH_map;
    for(int i=0;i<NH;++i){bool dir=false;for(int dv:coarse.dirichlet)if(dv==i)dir=true;if(!dir){dofH_map[i]=dofH.size();dofH.push_back(i);}}
    std::vector<Eigen::Triplet<double>> g0_t;
    for(int k2=0;k2<G.outerSize();++k2)
        for(Eigen::SparseMatrix<double>::InnerIterator it(G,k2);it;++it)
            if(dofH_map.count(it.col()))g0_t.emplace_back(it.row(),dofH_map[it.col()],it.value());
    Eigen::SparseMatrix<double> G0(Nh,dofH.size()); G0.setFromTriplets(g0_t.begin(),g0_t.end());
    Eigen::SparseMatrix<double> SHLOD0=G0.transpose()*Sh*G0;
    Eigen::VectorXd f_coarse=Eigen::VectorXd::Ones(NH);
    Eigen::VectorXd rhs=G0.transpose()*(Mh*(f_out.P_node*f_coarse));
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llts(SHLOD0);
    Eigen::VectorXd uH=Eigen::VectorXd::Zero(NH), uf=llts.solve(rhs);
    for(size_t j=0;j<dofH.size();++j)uH(dofH[j])=uf(j);
    Eigen::VectorXd uHms=G*uH;
    auto t5b = chr::high_resolution_clock::now();
    double t5 = chr::duration<double,std::milli>(t5b-t5a).count();
    std::cout << "5. Coarse solve: " << t5 << " ms (SHLOD0: " << dofH.size() << "x" << dofH.size() << ")\n";

    // 6. Reference
    auto t6a = chr::high_resolution_clock::now();
    std::vector<int> dofh; std::unordered_map<int,int> dofh_map;
    for(int i=0;i<Nh;++i){bool dir=false;for(int dv:fine.dirichlet)if(dv==i)dir=true;if(!dir){dofh_map[i]=dofh.size();dofh.push_back(i);}}
    std::vector<Eigen::Triplet<double>> sh_t;
    for(int k2=0;k2<Sh.outerSize();++k2)
        for(Eigen::SparseMatrix<double>::InnerIterator it(Sh,k2);it;++it)
            if(dofh_map.count(it.row())&&dofh_map.count(it.col()))
                sh_t.emplace_back(dofh_map[it.row()],dofh_map[it.col()],it.value());
    Eigen::SparseMatrix<double> Sh_free(dofh.size(),dofh.size()); Sh_free.setFromTriplets(sh_t.begin(),sh_t.end());
    Eigen::VectorXd f_fine=Eigen::VectorXd::Ones(Nh);
    Eigen::VectorXd rhs_ref=Eigen::VectorXd::Zero(dofh.size());
    for(int k2=0;k2<Mh.outerSize();++k2)
        for(Eigen::SparseMatrix<double>::InnerIterator it(Mh,k2);it;++it)
            if(dofh_map.count(it.row()))rhs_ref(dofh_map[it.row()])+=it.value()*f_fine(it.col());
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt_ref(Sh_free);
    Eigen::VectorXd uh=Eigen::VectorXd::Zero(Nh), uf2=llt_ref.solve(rhs_ref);
    for(size_t j=0;j<dofh.size();++j)uh(dofh[j])=uf2(j);
    auto t6b = chr::high_resolution_clock::now();
    double t6 = chr::duration<double,std::milli>(t6b-t6a).count();
    std::cout << "6. Reference: " << t6 << " ms (DOFs: " << dofh.size() << ")\n";

    // Errors
    Eigen::VectorXd diff=uh-uHms, diff_fe=uh-f_out.P_node*uH;
    double errE=std::sqrt(diff.dot(Sh*diff)), errL2=std::sqrt(diff.dot(Mh*diff)), errFE=std::sqrt(diff_fe.dot(Mh*diff_fe));
    double total=t1+t2+t3+t4+t5+t6;
    std::cout << "\nErrors: E=" << errE << " L2=" << errL2 << " FE=" << errFE << "\n";
    std::cout << "Total: " << total << " ms\n";
    return 0;
}

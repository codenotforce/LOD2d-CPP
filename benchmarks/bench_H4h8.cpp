/// H=4,h=8 benchmark — same pipeline as test_full.cpp (verified correct)
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
#include <Eigen/Dense>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace lod2d;
namespace chr = std::chrono;

int main() {
    int H=4, h=8, ell=2, d=2;
    std::cout << "=== C++ LOD H=" << H << " h=" << h << " ===\n";

    std::ifstream af("benchmarks/data_H4h8.txt");
    int nAh; af >> nAh;
    std::vector<double> Ah(nAh);
    for (int i=0;i<nAh;++i) af >> Ah[i];

    auto t0 = chr::high_resolution_clock::now();

    TriMesh T0;
    T0.nodes={{0,0},{1,0},{1,1},{0,1}}; T0.elems={{0,1,3},{1,2,3}}; T0.dirichlet={0,1,2,3};
    auto c_out=refine_mesh(T0,H), f_out=refine_mesh(c_out.mesh,h-H);
    const auto &coarse=c_out.mesh, &fine=f_out.mesh;
    int NH=coarse.nodes.size(), NTH=coarse.elems.size();
    int Nh=fine.nodes.size(), NTh_f=fine.elems.size();
    std::cout<<"Coarse:"<<NH<<"v "<<NTH<<"t  Fine:"<<Nh<<"v "<<NTh_f<<"t\n";

    std::vector<std::array<int,3>> dghidx(NTh_f), dgHidx(NTH);
    for(int e=0;e<NTh_f;++e)for(int i=0;i<3;++i)dghidx[e][i]=3*e+i;
    for(int e=0;e<NTH;++e)for(int i=0;i<3;++i)dgHidx[e][i]=3*e+i;

    std::vector<int> nngH(NH),nngh(Nh);
    for(auto&t:coarse.elems)for(int v:t)nngH[v]++; for(int v:coarse.dirichlet)nngH[v]=0;
    for(auto&t:fine.elems)for(int v:t)nngh[v]++; for(int v:fine.dirichlet)nngh[v]=0;

    std::vector<Eigen::Triplet<double>> cg_t;
    for(int e=0;e<NTh_f;++e)for(int i=0;i<3;++i)cg_t.emplace_back(3*e+i,fine.elems[e][i],1.0);
    Eigen::SparseMatrix<double> cg2dgh(3*NTh_f,Nh); cg2dgh.setFromTriplets(cg_t.begin(),cg_t.end());

    auto tp0=chr::high_resolution_clock::now();
    Eigen::SparseMatrix<double> Shdg=assemble_dg(fine,Ah);
    auto IH=build_quasi_interp(coarse,fine,f_out.P_dg,cg2dgh,Nh,NH);
    auto patch=build_patches(coarse,ell);
    auto areas=compute_area(fine);
    double M3[3][3]={{2,1,1},{1,2,1},{1,1,2}};
    std::vector<Eigen::Triplet<double>> mh_t;
    for(int e=0;e<NTh_f;++e){double s=areas[e]/12.0;for(int i=0;i<3;++i)for(int j=0;j<3;++j)mh_t.emplace_back(3*e+i,3*e+j,s*M3[i][j]);}
    Eigen::SparseMatrix<double> Mhdg(3*NTh_f,3*NTh_f); Mhdg.setFromTriplets(mh_t.begin(),mh_t.end());
    auto tp1=chr::high_resolution_clock::now();
    std::cout<<"Setup: "<<chr::duration<double,std::milli>(tp1-tp0).count()<<" ms\n";

    // ---- Correctors (parallel) ----
    auto tc0=chr::high_resolution_clock::now();
    std::vector<Eigen::SparseMatrix<double>> CT(NTH);
    #pragma omp parallel for schedule(dynamic)
    for(int k=0;k<NTH;++k)
        CT[k]=compute_corrector(k,patch,coarse,NH,nngH,f_out.P_elem,fine,Nh,nngh,dghidx,cg2dgh,Shdg,f_out.P_dg,dgHidx,IH,d);
    auto tc1=chr::high_resolution_clock::now();
    double tc=chr::duration<double,std::milli>(tc1-tc0).count();
    std::cout<<"Correctors: "<<tc<<" ms";
#ifdef _OPENMP
    std::cout<<" ("<<omp_get_max_threads()<<" threads)";
#endif
    std::cout<<"\n";

    // ---- C_ell = cell2mat(CT) * cg2dgH ----
    // Build coarse DG→CG: cg2dgH(3*e+i, vertex) = 1
    std::vector<Eigen::Triplet<double>> cgH_t;
    for(int e=0;e<NTH;++e)for(int i=0;i<3;++i)cgH_t.emplace_back(3*e+i,coarse.elems[e][i],1.0);
    Eigen::SparseMatrix<double> cg2dgH(3*NTH,NH); cg2dgH.setFromTriplets(cgH_t.begin(),cgH_t.end());

    // Horizontal concat: [CT0 CT1 ... CT_{NTH-1}]
    std::vector<Eigen::Triplet<double>> cell_t;
    for(int k=0;k<NTH;++k)
        for(int c=0;c<CT[k].outerSize();++c)
            for(Eigen::SparseMatrix<double>::InnerIterator it(CT[k],c);it;++it)
                cell_t.emplace_back(it.row(),k*(d+1)+(int)it.col(),it.value());
    Eigen::SparseMatrix<double> cell_mat(Nh,NTH*(d+1)); cell_mat.setFromTriplets(cell_t.begin(),cell_t.end());
    Eigen::SparseMatrix<double> C_ell = cell_mat * cg2dgH;
    Eigen::SparseMatrix<double> G = f_out.P_node - C_ell;

    // ---- Coarse LOD ----
    std::vector<int> dofH;
    for(int i=0;i<NH;++i){
        bool dir=false; for(int dv:coarse.dirichlet)if(dv==i)dir=true;
        if(!dir)dofH.push_back(i);
    }
    int nFree=dofH.size();
    std::unordered_map<int,int> dofH_map;
    for(int j=0;j<nFree;++j)dofH_map[dofH[j]]=j;

    // G0 = G(:,dofH)
    std::vector<Eigen::Triplet<double>> g0_t;
    for(int k2=0;k2<G.outerSize();++k2)
        for(Eigen::SparseMatrix<double>::InnerIterator it(G,k2);it;++it)
            if(dofH_map.count(it.col()))g0_t.emplace_back(it.row(),dofH_map[it.col()],it.value());
    Eigen::SparseMatrix<double> G0(Nh,nFree); G0.setFromTriplets(g0_t.begin(),g0_t.end());

    auto T=cg2dgh.transpose()*Shdg*cg2dgh;
    Eigen::SparseMatrix<double> SHLOD0=G0.transpose()*T*G0;
    Eigen::VectorXd f_coarse=Eigen::VectorXd::Ones(NH);
    Eigen::VectorXd rhs=G0.transpose()*(cg2dgh.transpose()*(Mhdg*(cg2dgh*(f_out.P_node*f_coarse))));
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llts(SHLOD0);
    Eigen::VectorXd uH_free=llts.solve(rhs), uH=Eigen::VectorXd::Zero(NH);
    for(int j=0;j<nFree;++j)uH(dofH[j])=uH_free(j);
    Eigen::VectorXd uHms=G*uH;

    // ---- Reference (same pattern as test_full.cpp) ----
    std::vector<int> dofh;
    for(int i=0;i<Nh;++i){
        bool dir=false; for(int dv:fine.dirichlet)if(dv==i)dir=true;
        if(!dir)dofh.push_back(i);
    }
    int nFree_f=dofh.size();
    std::unordered_map<int,int> dofh_map;
    for(int j=0;j<nFree_f;++j)dofh_map[dofh[j]]=j;

    Eigen::SparseMatrix<double> Sh=cg2dgh.transpose()*Shdg*cg2dgh;
    Eigen::SparseMatrix<double> Mh=cg2dgh.transpose()*Mhdg*cg2dgh;

    // Sh_free = Sh(dofh, dofh)
    std::vector<Eigen::Triplet<double>> sh_t;
    for(int k2=0;k2<Sh.outerSize();++k2)
        for(Eigen::SparseMatrix<double>::InnerIterator it(Sh,k2);it;++it)
            if(dofh_map.count(it.row())&&dofh_map.count(it.col()))
                sh_t.emplace_back(dofh_map[it.row()],dofh_map[it.col()],it.value());
    Eigen::SparseMatrix<double> Sh_free(nFree_f,nFree_f); Sh_free.setFromTriplets(sh_t.begin(),sh_t.end());

    // RHS_ref = Mh(dofh,:)*f
    Eigen::VectorXd f_fine=Eigen::VectorXd::Ones(Nh);
    Eigen::VectorXd rhs_ref=Eigen::VectorXd::Zero(nFree_f);
    for(int k2=0;k2<Mh.outerSize();++k2)
        for(Eigen::SparseMatrix<double>::InnerIterator it(Mh,k2);it;++it)
            if(dofh_map.count(it.row()))rhs_ref(dofh_map[it.row()])+=it.value()*f_fine(it.col());
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt_ref(Sh_free);
    Eigen::VectorXd uh_free=llt_ref.solve(rhs_ref), uh=Eigen::VectorXd::Zero(Nh);
    for(int j=0;j<nFree_f;++j)uh(dofh[j])=uh_free(j);

    // ---- Errors ----
    Eigen::VectorXd diff=uh-uHms, diff_fe=uh-f_out.P_node*uH;
    double errE=std::sqrt(diff.dot(Sh*diff)), errL2=std::sqrt(diff.dot(Mh*diff)), errFE=std::sqrt(diff_fe.dot(Mh*diff_fe));
    auto t1=chr::high_resolution_clock::now();
    double ms=chr::duration<double,std::milli>(t1-t0).count();

    double mE=2.478159e-02, mL2=1.665442e-03, mFE=1.458422e-02;
    std::cout<<"\nErrors: E="<<errE<<" (MATLAB "<<mE<<")  L2="<<errL2<<" (MATLAB "<<mL2<<")  FE="<<errFE<<" (MATLAB "<<mFE<<")\n";
    bool ok=(std::abs(errE-mE)/mE<1e-5 && std::abs(errL2-mL2)/mL2<1e-5 && std::abs(errFE-mFE)/mFE<1e-5);
    std::cout<<(ok?"[PASS]":"[FAIL]")<<" All errors\n";
    std::cout<<"\nC++: "<<ms<<" ms  MATLAB: ~30400 ms  Speedup: "<<30400.0/ms<<"x\n";
    return ok?0:1;
}

/// CHOLMOD vs Eigen performance test
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <suitesparse/cholmod.h>
#include <iostream>
#include <chrono>
#include <cmath>

using namespace std::chrono;

// Convert Eigen sparse → CHOLMOD sparse (int API, safe copy)
cholmod_sparse* eigen_to_cholmod(const Eigen::SparseMatrix<double>& A, cholmod_common* c) {
    Eigen::SparseMatrix<double> Ac = A; Ac.makeCompressed();
    cholmod_triplet *T = cholmod_allocate_triplet(Ac.rows(), Ac.cols(), Ac.nonZeros(), 1, CHOLMOD_REAL, c);
    int *Ti=(int*)T->i, *Tj=(int*)T->j; double *Tx=(double*)T->x; int idx=0;
    for(int k=0;k<Ac.outerSize();++k)
        for(Eigen::SparseMatrix<double>::InnerIterator it(Ac,k);it;++it,++idx)
            {Ti[idx]=(int)it.row(); Tj[idx]=(int)it.col(); Tx[idx]=it.value();}
    T->nnz = Ac.nonZeros();
    cholmod_sparse *Asp = cholmod_triplet_to_sparse(T, T->nnz, c);
    Asp->stype = -1;
    cholmod_free_triplet(&T, c);
    return Asp;
}

int main() {
    std::cout << "=== CHOLMOD vs Eigen SimplicialLLT ===\n\n";

    for (int N : {500, 2000, 4000, 8000}) {
        int nrhs = (N <= 2000) ? 30 : 10;

        // Build sparse SPD (2D Poisson-like)
        Eigen::SparseMatrix<double> A(N,N);
        std::vector<Eigen::Triplet<double>> t;
        double h2 = N*N;  // scaled
        for (int i=0;i<N;++i) {
            t.emplace_back(i,i,4.0*h2);
            if(i+1<N){t.emplace_back(i,i+1,-1.0*h2);t.emplace_back(i+1,i,-1.0*h2);}
        }
        A.setFromTriplets(t.begin(),t.end());

        Eigen::MatrixXd B = Eigen::MatrixXd::Random(N, nrhs);

        // Eigen SimplicialLLT
        Eigen::MatrixXd Xe(N,nrhs);
        auto t0=high_resolution_clock::now();
        {
            Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt(A);
            for(int j=0;j<nrhs;++j) Xe.col(j)=llt.solve(B.col(j));
        }
        auto t1=high_resolution_clock::now();
        double ms_e = duration<double,std::milli>(t1-t0).count();

        // CHOLMOD
        Eigen::MatrixXd Xc(N,nrhs);
        auto t2=high_resolution_clock::now();
        {
            cholmod_common c; cholmod_start(&c);
            c.supernodal = CHOLMOD_SUPERNODAL; c.print = 0;

            cholmod_sparse *Asp = eigen_to_cholmod(A, &c);
            cholmod_factor *L = cholmod_analyze(Asp, &c);
            cholmod_factorize(Asp, L, &c);

            cholmod_dense *Bcd = cholmod_allocate_dense(N, nrhs, N, CHOLMOD_REAL, &c);
            double *bx=(double*)Bcd->x;
            for(int j=0;j<nrhs;++j) for(int i=0;i<N;++i) bx[j*N+i]=B(i,j);

            cholmod_dense *Xcd = cholmod_solve(CHOLMOD_A, L, Bcd, &c);
            double *xx=(double*)Xcd->x;
            for(int j=0;j<nrhs;++j) for(int i=0;i<N;++i) Xc(i,j)=xx[j*N+i];

            cholmod_free_sparse(&Asp,&c); cholmod_free_factor(&L,&c);
            cholmod_free_dense(&Bcd,&c); cholmod_free_dense(&Xcd,&c);
            cholmod_finish(&c);
        }
        auto t3=high_resolution_clock::now();
        double ms_c = duration<double,std::milli>(t3-t2).count();

        double err = (Xe-Xc).cwiseAbs().maxCoeff();
        std::cout << "N=" << N << " (" << nrhs << " RHS):  Eigen=" << ms_e << "ms  CHOLMOD=" << ms_c
                  << "ms  ratio=" << ms_e/ms_c << "x  err=" << err << "\n";
    }

    return 0;
}

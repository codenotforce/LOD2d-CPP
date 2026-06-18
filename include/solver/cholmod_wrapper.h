#pragma once
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <suitesparse/cholmod.h>
#include <vector>

namespace lod2d {

/// Solve A * X = B using CHOLMOD sparse Cholesky.
/// A must be SPD.  B is dense (multiple RHS columns).
/// Returns X with same dimensions as B.
inline Eigen::MatrixXd solve_cholmod(
    const Eigen::SparseMatrix<double> &A,
    const Eigen::MatrixXd &B)
{
    int n = static_cast<int>(A.rows());
    int nrhs = static_cast<int>(B.cols());

    cholmod_common c;
    cholmod_l_start(&c);
    c.supernodal = CHOLMOD_SUPERNODAL;  // faster for large systems
    c.print = 0;  // suppress output in parallel

    // Convert Eigen CSC → CHOLMOD sparse (make a compressed copy)
    Eigen::SparseMatrix<double> Ac = A;
    Ac.makeCompressed();
    cholmod_sparse *Asp = cholmod_l_allocate_sparse(
        n, n, Ac.nonZeros(), 1, 1, 0, CHOLMOD_REAL, &c);
    Asp->p = const_cast<int*>(Ac.outerIndexPtr());
    Asp->i = const_cast<int*>(Ac.innerIndexPtr());
    Asp->x = const_cast<double*>(Ac.valuePtr());
    Asp->stype = 0;  // full matrix (both triangles stored in Eigen CSC)

    // Factorize
    cholmod_factor *L = cholmod_l_analyze(Asp, &c);
    cholmod_l_factorize(Asp, L, &c);

    // Solve: allocate dense RHS, solve, copy back
    cholmod_dense *Bcd = cholmod_l_allocate_dense(n, nrhs, n, CHOLMOD_REAL, &c);
    double *Bx = static_cast<double*>(Bcd->x);
    for (int j = 0; j < nrhs; ++j)
        for (int i = 0; i < n; ++i)
            Bx[j * n + i] = B(i, j);

    cholmod_dense *Xcd = cholmod_l_solve(CHOLMOD_A, L, Bcd, &c);

    // Copy result back
    Eigen::MatrixXd X(n, nrhs);
    double *Xx = static_cast<double*>(Xcd->x);
    for (int j = 0; j < nrhs; ++j)
        for (int i = 0; i < n; ++i)
            X(i, j) = Xx[j * n + i];

    // Cleanup — don't free Asp->p/i/x (owned by Eigen)
    Asp->p = nullptr; Asp->i = nullptr; Asp->x = nullptr;
    cholmod_l_free_sparse(&Asp, &c);
    cholmod_l_free_factor(&L, &c);
    cholmod_l_free_dense(&Bcd, &c);
    cholmod_l_free_dense(&Xcd, &c);
    cholmod_l_finish(&c);

    return X;
}

} // namespace lod2d

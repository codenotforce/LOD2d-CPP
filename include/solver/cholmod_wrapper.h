#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <suitesparse/cholmod.h>
#include <stdexcept>
#include <string>

namespace lod2d {

/// Solve A * X = B using CHOLMOD sparse Cholesky.
/// A must be symmetric positive definite. Only A's lower triangle is copied.
inline Eigen::MatrixXd solve_cholmod(
    const Eigen::SparseMatrix<double> &A,
    const Eigen::MatrixXd &B)
{
    const int n = static_cast<int>(A.rows());
    const int nrhs = static_cast<int>(B.cols());
    if (A.rows() != A.cols())
        throw std::invalid_argument("solve_cholmod requires a square matrix");
    if (B.rows() != n)
        throw std::invalid_argument("solve_cholmod RHS row count does not match A");

    cholmod_common c;
    cholmod_start(&c);
    c.supernodal = CHOLMOD_SUPERNODAL;
    c.print = 0;

    auto fail = [&](const std::string &msg) -> void {
        cholmod_finish(&c);
        throw std::runtime_error("CHOLMOD failed: " + msg);
    };

    Eigen::SparseMatrix<double> Ac = A;
    Ac.makeCompressed();

    cholmod_triplet *T = cholmod_allocate_triplet(
        n, n, Ac.nonZeros(), -1, CHOLMOD_REAL, &c);
    if (!T) fail("allocate sparse triplet");

    int *Ti = static_cast<int*>(T->i);
    int *Tj = static_cast<int*>(T->j);
    double *Tx = static_cast<double*>(T->x);
    int idx = 0;
    for (int col = 0; col < Ac.outerSize(); ++col) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Ac, col); it; ++it) {
            if (it.row() < col) continue;
            Ti[idx] = static_cast<int>(it.row());
            Tj[idx] = col;
            Tx[idx] = it.value();
            ++idx;
        }
    }
    T->nnz = idx;

    cholmod_sparse *Asp = cholmod_triplet_to_sparse(T, T->nnz, &c);
    cholmod_free_triplet(&T, &c);
    if (!Asp) fail("convert triplet to sparse");
    Asp->stype = -1;

    cholmod_factor *L = cholmod_analyze(Asp, &c);
    if (!L) {
        cholmod_free_sparse(&Asp, &c);
        fail("analyze");
    }
    if (!cholmod_factorize(Asp, L, &c) || c.status < CHOLMOD_OK) {
        cholmod_free_sparse(&Asp, &c);
        cholmod_free_factor(&L, &c);
        fail("factorize");
    }

    cholmod_dense *Bcd = cholmod_allocate_dense(n, nrhs, n, CHOLMOD_REAL, &c);
    if (!Bcd) {
        cholmod_free_sparse(&Asp, &c);
        cholmod_free_factor(&L, &c);
        fail("allocate dense RHS");
    }

    double *Bx = static_cast<double*>(Bcd->x);
    for (int j = 0; j < nrhs; ++j)
        for (int i = 0; i < n; ++i)
            Bx[j * n + i] = B(i, j);

    cholmod_dense *Xcd = cholmod_solve(CHOLMOD_A, L, Bcd, &c);
    if (!Xcd) {
        cholmod_free_sparse(&Asp, &c);
        cholmod_free_factor(&L, &c);
        cholmod_free_dense(&Bcd, &c);
        fail("solve");
    }

    Eigen::MatrixXd X(n, nrhs);
    double *Xx = static_cast<double*>(Xcd->x);
    for (int j = 0; j < nrhs; ++j)
        for (int i = 0; i < n; ++i)
            X(i, j) = Xx[j * n + i];

    cholmod_free_sparse(&Asp, &c);
    cholmod_free_factor(&L, &c);
    cholmod_free_dense(&Bcd, &c);
    cholmod_free_dense(&Xcd, &c);
    cholmod_finish(&c);

    return X;
}

} // namespace lod2d

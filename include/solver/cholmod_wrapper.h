#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <suitesparse/cholmod.h>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace lod2d {

namespace detail {

constexpr std::size_t kMaxCachedCholmodPatterns = 1;

struct CholmodPatternKey {
    int n = 0;
    std::vector<int> outer;
    std::vector<int> inner;

    bool operator==(const CholmodPatternKey &other) const {
        return n == other.n && outer == other.outer && inner == other.inner;
    }
};

struct CholmodPatternHash {
    std::size_t operator()(const CholmodPatternKey &key) const {
        std::size_t h = static_cast<std::size_t>(key.n) * 1469598103934665603ull;
        auto mix = [&h](int value) {
            h ^= static_cast<std::size_t>(value) + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
        };
        for (int value : key.outer) mix(value);
        for (int value : key.inner) mix(value);
        return h;
    }
};

struct CholmodCachedContext {
    cholmod_common common;
    std::unordered_map<CholmodPatternKey, cholmod_factor*, CholmodPatternHash> factors;

    CholmodCachedContext() {
        cholmod_start(&common);
        common.supernodal = CHOLMOD_SUPERNODAL;
        common.print = 0;
    }

    CholmodCachedContext(const CholmodCachedContext&) = delete;
    CholmodCachedContext& operator=(const CholmodCachedContext&) = delete;

    ~CholmodCachedContext() {
        for (auto &entry : factors) {
            cholmod_factor *factor = entry.second;
            if (factor) cholmod_free_factor(&factor, &common);
        }
        cholmod_finish(&common);
    }
};

inline CholmodPatternKey pattern_key_from_sparse(const cholmod_sparse *A) {
    CholmodPatternKey key;
    key.n = static_cast<int>(A->nrow);

    const int ncol = static_cast<int>(A->ncol);
    const int *Ap = static_cast<const int*>(A->p);
    const int *Ai = static_cast<const int*>(A->i);
    const int nnz = Ap[ncol];

    key.outer.assign(Ap, Ap + ncol + 1);
    key.inner.assign(Ai, Ai + nnz);
    return key;
}

inline cholmod_sparse* copy_lower_triangle_to_cholmod(
    const Eigen::SparseMatrix<double> &A,
    cholmod_common *common)
{
    const int n = static_cast<int>(A.rows());
    Eigen::SparseMatrix<double> Ac = A;
    Ac.makeCompressed();

    cholmod_triplet *T = cholmod_allocate_triplet(
        n, n, Ac.nonZeros(), -1, CHOLMOD_REAL, common);
    if (!T) return nullptr;

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

    cholmod_sparse *Asp = cholmod_triplet_to_sparse(T, T->nnz, common);
    cholmod_free_triplet(&T, common);
    if (Asp) Asp->stype = -1;
    return Asp;
}

inline Eigen::MatrixXd solve_cholmod_impl(
    const Eigen::SparseMatrix<double> &A,
    const Eigen::MatrixXd &B,
    CholmodCachedContext *cached_context)
{
    const int n = static_cast<int>(A.rows());
    const int nrhs = static_cast<int>(B.cols());
    if (A.rows() != A.cols())
        throw std::invalid_argument("solve_cholmod requires a square matrix");
    if (B.rows() != n)
        throw std::invalid_argument("solve_cholmod RHS row count does not match A");

    cholmod_common local_common;
    cholmod_common *common = nullptr;
    if (cached_context) {
        common = &cached_context->common;
    } else {
        cholmod_start(&local_common);
        local_common.supernodal = CHOLMOD_SUPERNODAL;
        local_common.print = 0;
        common = &local_common;
    }
    common->status = CHOLMOD_OK;

    auto finish_local = [&]() {
        if (!cached_context) cholmod_finish(common);
    };
    auto fail = [&](const std::string &msg) -> void {
        finish_local();
        throw std::runtime_error("CHOLMOD failed: " + msg);
    };

    cholmod_sparse *Asp = copy_lower_triangle_to_cholmod(A, common);
    if (!Asp) fail("copy sparse matrix");

    cholmod_factor *L = nullptr;
    bool inserted_factor = false;
    typename std::unordered_map<CholmodPatternKey, cholmod_factor*, CholmodPatternHash>::iterator cached_it;
    if (cached_context) {
        CholmodPatternKey key = pattern_key_from_sparse(Asp);
        if (cached_context->factors.find(key) == cached_context->factors.end() &&
            cached_context->factors.size() >= kMaxCachedCholmodPatterns) {
            for (auto &entry : cached_context->factors) {
                cholmod_factor *factor = entry.second;
                if (factor) cholmod_free_factor(&factor, common);
            }
            cached_context->factors.clear();
        }
        auto result = cached_context->factors.emplace(std::move(key), nullptr);
        cached_it = result.first;
        inserted_factor = result.second;
        if (inserted_factor) {
            cached_it->second = cholmod_analyze(Asp, common);
            if (!cached_it->second || common->status < CHOLMOD_OK) {
                if (cached_it->second) cholmod_free_factor(&cached_it->second, common);
                cached_context->factors.erase(cached_it);
                cholmod_free_sparse(&Asp, common);
                fail("analyze cached pattern");
            }
        }
        L = cached_it->second;
    } else {
        L = cholmod_analyze(Asp, common);
        if (!L || common->status < CHOLMOD_OK) {
            cholmod_free_sparse(&Asp, common);
            fail("analyze");
        }
    }

    common->status = CHOLMOD_OK;
    if (!cholmod_factorize(Asp, L, common) || common->status < CHOLMOD_OK) {
        cholmod_free_sparse(&Asp, common);
        if (cached_context) {
            if (inserted_factor && L) {
                cholmod_free_factor(&L, common);
                cached_context->factors.erase(cached_it);
            }
        } else {
            cholmod_free_factor(&L, common);
        }
        fail("factorize");
    }

    cholmod_dense *Bcd = cholmod_allocate_dense(n, nrhs, n, CHOLMOD_REAL, common);
    if (!Bcd) {
        cholmod_free_sparse(&Asp, common);
        if (!cached_context) cholmod_free_factor(&L, common);
        fail("allocate dense RHS");
    }

    double *Bx = static_cast<double*>(Bcd->x);
    for (int j = 0; j < nrhs; ++j)
        for (int i = 0; i < n; ++i)
            Bx[j * n + i] = B(i, j);

    common->status = CHOLMOD_OK;
    cholmod_dense *Xcd = cholmod_solve(CHOLMOD_A, L, Bcd, common);
    if (!Xcd || common->status < CHOLMOD_OK) {
        cholmod_free_sparse(&Asp, common);
        if (!cached_context) cholmod_free_factor(&L, common);
        cholmod_free_dense(&Bcd, common);
        fail("solve");
    }

    Eigen::MatrixXd X(n, nrhs);
    double *Xx = static_cast<double*>(Xcd->x);
    for (int j = 0; j < nrhs; ++j)
        for (int i = 0; i < n; ++i)
            X(i, j) = Xx[j * n + i];

    cholmod_free_sparse(&Asp, common);
    if (!cached_context) cholmod_free_factor(&L, common);
    cholmod_free_dense(&Bcd, common);
    cholmod_free_dense(&Xcd, common);
    finish_local();

    return X;
}

} // namespace detail

/// Solve A * X = B using CHOLMOD sparse Cholesky.
/// A must be symmetric positive definite. Only A's lower triangle is copied.
inline Eigen::MatrixXd solve_cholmod(
    const Eigen::SparseMatrix<double> &A,
    const Eigen::MatrixXd &B)
{
    return detail::solve_cholmod_impl(A, B, nullptr);
}

/// Solve A * X = B using CHOLMOD with thread-local symbolic factor reuse.
/// The numeric factorization is refreshed every call; only matching sparsity
/// patterns reuse CHOLMOD's analyzed factor object.
inline Eigen::MatrixXd solve_cholmod_cached(
    const Eigen::SparseMatrix<double> &A,
    const Eigen::MatrixXd &B)
{
    thread_local detail::CholmodCachedContext context;
    return detail::solve_cholmod_impl(A, B, &context);
}

} // namespace lod2d


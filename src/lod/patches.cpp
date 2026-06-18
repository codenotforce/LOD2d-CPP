#include "lod/patches.h"
#include <Eigen/Sparse>
#include <vector>

namespace lod2d {

Eigen::SparseMatrix<double> build_patches(const TriMesh &coarse, int ell) {
    int NH  = static_cast<int>(coarse.nodes.size());
    int NTH = static_cast<int>(coarse.elems.size());
    int d   = 2;  // 2D

    // 1. Vertex-to-element incidence: Ivt(NH × NTH)
    std::vector<Eigen::Triplet<double>> ivt_t;
    for (int t = 0; t < NTH; ++t) {
        for (int i = 0; i < d + 1; ++i) {
            int v = coarse.elems[t][i];
            ivt_t.emplace_back(v, t, 1.0);
        }
    }
    Eigen::SparseMatrix<double> Ivt(NH, NTH);
    Ivt.setFromTriplets(ivt_t.begin(), ivt_t.end());

    // 2. Element adjacency: Itt = spones(Ivt' * Ivt)
    Eigen::SparseMatrix<double> Itt = (Ivt.transpose() * Ivt).pruned();
    // Threshold to binary
    for (int k = 0; k < Itt.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Itt, k); it; ++it) {
            it.valueRef() = 1.0;
        }
    }

    // 3. ℓ‑step neighborhood expansion
    Eigen::SparseMatrix<double> patch(NTH, NTH);
    patch.setIdentity();

    for (int lev = 0; lev < ell; ++lev) {
        patch = Itt * patch;
    }

    // 4. Binary threshold: convert to explicit triplets, then rebuild
    std::vector<Eigen::Triplet<double>> bin_t;
    for (int k = 0; k < patch.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(patch, k); it; ++it) {
            bin_t.emplace_back(it.row(), it.col(), 1.0);
        }
    }
    // Diagonal (MATLAB: patch = spones(patch) + speye(NTH))
    for (int i = 0; i < NTH; ++i) {
        bin_t.emplace_back(i, i, 1.0);
    }

    Eigen::SparseMatrix<double> result(NTH, NTH);
    result.setFromTriplets(bin_t.begin(), bin_t.end());
    return result;
}

} // namespace lod2d

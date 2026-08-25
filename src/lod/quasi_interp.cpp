#include "lod/quasi_interp.h"
#include "mesh/refine.h"
#include <Eigen/Sparse>
#include <algorithm>
#include <array>
#include <stdexcept>
#include <vector>

namespace lod2d {

namespace {

Eigen::SparseMatrix<double> build_quasi_interp_sparse_products(
    const TriMesh &coarse,
    const TriMesh &fine,
    const Eigen::SparseMatrix<double> &P1dg,
    const Eigen::SparseMatrix<double> &cg2dgh) {
    const int NTh = static_cast<int>(fine.elems.size());
    const int NTH = static_cast<int>(coarse.elems.size());
    const int NHdg = 3 * NTH;
    const auto fine_areas = compute_area(fine);
    const double M3[3][3] = {{2,1,1},{1,2,1},{1,1,2}};
    std::vector<Eigen::Triplet<double>> mhdg_t;
    mhdg_t.reserve(9 * static_cast<std::size_t>(NTh));
    for (int e = 0; e < NTh; ++e) {
        const double scale = fine_areas[e] / 12.0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                mhdg_t.emplace_back(3*e+i, 3*e+j, scale * M3[i][j]);
    }
    Eigen::SparseMatrix<double> Mhdg(3*NTh, 3*NTh);
    Mhdg.setFromTriplets(mhdg_t.begin(), mhdg_t.end());

    const auto coarse_areas = compute_area(coarse);
    const double B3[3][3] = {{9,-3,-3},{-3,9,-3},{-3,-3,9}};
    std::vector<Eigen::Triplet<double>> b_t;
    b_t.reserve(9 * static_cast<std::size_t>(NTH));
    for (int e = 0; e < NTH; ++e) {
        const double scale = 1.0 / coarse_areas[e];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                b_t.emplace_back(3*e+i, 3*e+j, scale * B3[i][j]);
    }
    Eigen::SparseMatrix<double> B(NHdg, NHdg);
    B.setFromTriplets(b_t.begin(), b_t.end());
    const Eigen::SparseMatrix<double> T1 = Mhdg * cg2dgh;
    const Eigen::SparseMatrix<double> T2 = P1dg.transpose() * T1;
    const Eigen::SparseMatrix<double> PiHdg = B * T2;

    const int coarse_nodes = static_cast<int>(coarse.nodes.size());
    std::vector<Eigen::Triplet<double>> incidence_triplets;
    incidence_triplets.reserve(3 * static_cast<std::size_t>(NTH));
    for (int e = 0; e < NTH; ++e)
        for (int i = 0; i < 3; ++i)
            incidence_triplets.emplace_back(coarse.elems[e][i], 3*e+i, 1.0);
    Eigen::SparseMatrix<double> incidence(coarse_nodes, NHdg);
    incidence.setFromTriplets(
        incidence_triplets.begin(), incidence_triplets.end());
    const Eigen::VectorXd degree =
        incidence * Eigen::VectorXd::Ones(NHdg);
    std::vector<Eigen::Triplet<double>> averaging_triplets;
    averaging_triplets.reserve(incidence.nonZeros());
    for (int column = 0; column < incidence.outerSize(); ++column) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(
                 incidence, column); it; ++it) {
            averaging_triplets.emplace_back(
                it.row(), it.col(), 1.0 / degree(it.row()));
        }
    }
    Eigen::SparseMatrix<double> averaging(coarse_nodes, NHdg);
    averaging.setFromTriplets(
        averaging_triplets.begin(), averaging_triplets.end());
    const Eigen::SparseMatrix<double> unrestricted = averaging * PiHdg;
    std::vector<char> is_dirichlet(coarse.nodes.size(), false);
    for (int node : coarse.dirichlet)
        if (node >= 0 && node < coarse_nodes) is_dirichlet[node] = true;
    std::vector<Eigen::Triplet<double>> restricted_triplets;
    restricted_triplets.reserve(unrestricted.nonZeros());
    for (int column = 0; column < unrestricted.outerSize(); ++column) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(
                 unrestricted, column); it; ++it) {
            if (!is_dirichlet[it.row()] && it.value() != 0.0)
                restricted_triplets.emplace_back(
                    it.row(), it.col(), it.value());
        }
    }
    Eigen::SparseMatrix<double> result(
        unrestricted.rows(), unrestricted.cols());
    result.setFromTriplets(
        restricted_triplets.begin(), restricted_triplets.end());
    return result;
}

} // namespace

Eigen::SparseMatrix<double>
build_quasi_interp(const TriMesh &coarse, const TriMesh &fine,
                   const Eigen::SparseMatrix<double> &P1dg,
                   const Eigen::SparseMatrix<double> &cg2dgh,
                   int /*Nh*/, int /*NH*/) {
    int NTh  = static_cast<int>(fine.elems.size());
    int NTH  = static_cast<int>(coarse.elems.size());
    int NHdg = 3 * NTH;
    if (P1dg.rows() != 3 * NTh || P1dg.cols() != NHdg
        || cg2dgh.rows() != 3 * NTh
        || cg2dgh.cols() != static_cast<int>(fine.nodes.size())) {
        throw std::invalid_argument(
            "quasi interpolation matrices have inconsistent dimensions");
    }
    // Sparse products are faster on small and medium hierarchies.  Their
    // temporary matrices become superlinear in the deep moving-reference
    // regime, where the element-local assembly below has the lower wall time
    // and memory growth.  The threshold is deliberately conservative and is
    // covered by identical-trajectory gates on both sides.
    constexpr int local_assembly_element_threshold = 80000;
    if (NTh < local_assembly_element_threshold) {
        return build_quasi_interp_sparse_products(
            coarse, fine, P1dg, cg2dgh);
    }

    // Algebraically this routine evaluates
    //
    //   I_H = E_H B_H P_{H,h}^{T} M_h C_h.
    //
    // All four factors are element-local before the final CG assembly.  The
    // historical implementation materialised the three large sparse products
    // on the right.  On moving-reference trajectories this became a dominant
    // cost (and memory peak) once the candidate reached O(10^5) elements.
    // Assemble the identical 3-by-3 element contributions directly instead.
    const auto fine_areas = compute_area(fine);
    const auto coarse_areas = compute_area(coarse);
    const Eigen::Matrix3d M3 = (Eigen::Matrix3d() <<
        2.0, 1.0, 1.0,
        1.0, 2.0, 1.0,
        1.0, 1.0, 2.0).finished();
    const Eigen::Matrix3d B3 = (Eigen::Matrix3d() <<
         9.0, -3.0, -3.0,
        -3.0,  9.0, -3.0,
        -3.0, -3.0,  9.0).finished();

    // E_H is the vertex-wise average of coarse DG values.
    int Nh_c = static_cast<int>(coarse.nodes.size());
    std::vector<int> vertex_degree(Nh_c, 0);
    for (int e = 0; e < NTH; ++e) {
        for (int i = 0; i < 3; ++i)
            ++vertex_degree[coarse.elems[e][i]];
    }
    std::vector<char> is_dirichlet(coarse.nodes.size(), false);
    for (int node : coarse.dirichlet) {
        if (node >= 0 && node < static_cast<int>(is_dirichlet.size()))
            is_dirichlet[node] = true;
    }

    using RowMajorSparse = Eigen::SparseMatrix<double, Eigen::RowMajor>;
    const RowMajorSparse prolongation(P1dg);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(9 * static_cast<std::size_t>(NTh));
    for (int fine_element = 0; fine_element < NTh; ++fine_element) {
        Eigen::Matrix3d local_prolongation = Eigen::Matrix3d::Zero();
        int coarse_parent = -1;
        for (int fine_local = 0; fine_local < 3; ++fine_local) {
            const int row = 3 * fine_element + fine_local;
            for (RowMajorSparse::InnerIterator it(prolongation, row); it; ++it) {
                const int parent = it.col() / 3;
                if (coarse_parent < 0) coarse_parent = parent;
                if (parent != coarse_parent) {
                    throw std::runtime_error(
                        "DG prolongation couples one fine element to multiple coarse elements");
                }
                local_prolongation(fine_local, it.col() % 3) = it.value();
            }
        }
        if (coarse_parent < 0 || coarse_parent >= NTH
            || !(coarse_areas[coarse_parent] > 0.0)) {
            throw std::runtime_error(
                "DG prolongation has no valid coarse parent element");
        }
        const Eigen::Matrix3d local =
            (B3 / coarse_areas[coarse_parent])
            * local_prolongation.transpose()
            * (fine_areas[fine_element] / 12.0 * M3);
        for (int coarse_local = 0; coarse_local < 3; ++coarse_local) {
            const int coarse_node = coarse.elems[coarse_parent][coarse_local];
            if (is_dirichlet[coarse_node]) continue;
            const double averaging = 1.0 / vertex_degree[coarse_node];
            for (int fine_local = 0; fine_local < 3; ++fine_local) {
                const double value = averaging * local(coarse_local, fine_local);
                if (value != 0.0) {
                    triplets.emplace_back(
                        coarse_node, fine.elems[fine_element][fine_local], value);
                }
            }
        }
    }
    Eigen::SparseMatrix<double> result(Nh_c, fine.nodes.size());
    result.setFromTriplets(triplets.begin(), triplets.end());
    result.makeCompressed();
    return result;
}

} // namespace lod2d

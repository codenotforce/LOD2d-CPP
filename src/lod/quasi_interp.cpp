#include "lod/quasi_interp.h"
#include "mesh/refine.h"
#include <Eigen/Sparse>
#include <vector>

namespace lod2d {

Eigen::SparseMatrix<double>
build_quasi_interp(const TriMesh &coarse, const TriMesh &fine,
                   const Eigen::SparseMatrix<double> &P1dg,
                   const Eigen::SparseMatrix<double> &cg2dgh,
                   int /*Nh*/, int /*NH*/) {
    int NTh  = static_cast<int>(fine.elems.size());
    int NTH  = static_cast<int>(coarse.elems.size());
    int NHdg = 3 * NTH;

    // ---- 1. Mhdg: fine DG mass  ----
    // block-diagonal: kron(diag(area/12), [2 1 1; 1 2 1; 1 1 2])
    auto fine_areas = compute_area(fine);
    const double M3[3][3] = {{2,1,1},{1,2,1},{1,1,2}};
    std::vector<Eigen::Triplet<double>> mhdg_t;
    for (int e = 0; e < NTh; ++e) {
        double s = fine_areas[e] / 12.0;
        int dg[3] = {3*e, 3*e+1, 3*e+2};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                mhdg_t.emplace_back(dg[i], dg[j], s * M3[i][j]);
    }
    Eigen::SparseMatrix<double> Mhdg(3*NTh, 3*NTh);
    Mhdg.setFromTriplets(mhdg_t.begin(), mhdg_t.end());

    // ---- 2. B: inverse coarse DG mass ----
    // block-diagonal: kron(diag(1/area), [9 -3 -3; -3 9 -3; -3 -3 9])
    auto coarse_areas = compute_area(coarse);
    const double B3[3][3] = {{9,-3,-3},{-3,9,-3},{-3,-3,9}};
    std::vector<Eigen::Triplet<double>> b_t;
    for (int e = 0; e < NTH; ++e) {
        double s = 1.0 / coarse_areas[e];
        int dg[3] = {3*e, 3*e+1, 3*e+2};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                b_t.emplace_back(dg[i], dg[j], s * B3[i][j]);
    }
    Eigen::SparseMatrix<double> B(NHdg, NHdg);
    B.setFromTriplets(b_t.begin(), b_t.end());

    // ---- 3. PiHdg = B * (P1dg' * Mhdg * cg2dgh) ----
    Eigen::SparseMatrix<double> T1 = Mhdg * cg2dgh;
    Eigen::SparseMatrix<double> T2 = P1dg.transpose() * T1;
    Eigen::SparseMatrix<double> PiHdg = B * T2;

    // ---- 4. EH: coarse DG → CG averaging ----
    int Nh_c = static_cast<int>(coarse.nodes.size());
    double inv3 = 1.0 / 3.0;  // each interior vertex has exactly 3 incident
                               // coarse triangles on a regular mesh;
                               // boundary vertices have fewer.
    // Build cg2dgH: sparse(NHdg, Nh_c) with 1 at (dg_dof, vertex)
    // Each coarse triangle has 3 DG dofs, each at one vertex.
    std::vector<Eigen::Triplet<double>> eh_t;
    for (int e = 0; e < NTH; ++e) {
        for (int i = 0; i < 3; ++i) {
            int v = coarse.elems[e][i];
            int dg = 3*e + i;
            eh_t.emplace_back(v, dg, 1.0);
        }
    }
    Eigen::SparseMatrix<double> cg2dgH_t(Nh_c, NHdg);  // transposed: CG×DG
    cg2dgH_t.setFromTriplets(eh_t.begin(), eh_t.end());

    // Row sums of cg2dgH_t = vertex degree (number of incident coarse triangles)
    Eigen::VectorXd vtx_degree = cg2dgH_t * Eigen::VectorXd::Ones(NHdg);

    // EH = diag(1./vtx_degree) * cg2dgH_t
    std::vector<Eigen::Triplet<double>> eh2_t;
    for (int k = 0; k < cg2dgH_t.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(cg2dgH_t, k); it; ++it) {
            int v = static_cast<int>(it.row());
            int dg = static_cast<int>(it.col());
            eh2_t.emplace_back(v, dg, 1.0 / vtx_degree(v));
        }
    }
    Eigen::SparseMatrix<double> EH(Nh_c, NHdg);
    EH.setFromTriplets(eh2_t.begin(), eh2_t.end());

    // ---- 5. IH = EH * PiHdg ----
    return EH * PiHdg;
}

} // namespace lod2d

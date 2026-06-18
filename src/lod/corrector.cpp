#include "lod/corrector.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <set>
#include <algorithm>
// #include "solver/cholmod_wrapper.h"  // disabled until debugged

namespace lod2d {

Eigen::SparseMatrix<double>
compute_corrector(int k,
    const Eigen::SparseMatrix<double> &patch,
    const TriMesh &coarse, int NH, const std::vector<int> &nngH,
    const Eigen::SparseMatrix<double> &P0,
    const TriMesh &fine,   int Nh, const std::vector<int> &nngh,
    const std::vector<std::array<int,3>> &dghidx,
    const Eigen::SparseMatrix<double> &cg2dgh,
    const Eigen::SparseMatrix<double> &Shdg,
    const Eigen::SparseMatrix<double> &P1dg,
    const std::vector<std::array<int,3>> &dgHidx,
    const Eigen::SparseMatrix<double> &IH,
    int d) {

    // ---- 1. Coarse patch DOFs ----
    Eigen::VectorXd pk = patch.col(k);
    // tpH: which coarse elements are in patch
    std::vector<int> tpH_list;
    for (int j = 0; j < pk.size(); ++j)
        if (pk[j] != 0.0) tpH_list.push_back(j);

    // dofpH: unique vertices of patch coarse elements, excluding boundary
    std::vector<int> vtx_count(NH, 0);
    for (int j : tpH_list)
        for (int vi : coarse.elems[j])
            vtx_count[vi]++;
    std::vector<int> dofpH;
    for (int v = 0; v < NH; ++v)
        if (vtx_count[v] > 0 && nngH[v] > 0)
            dofpH.push_back(v);

    // ---- 2. Fine patch DOFs ----
    // tph = P0 * patch(:,k) → which fine elements are in patch
    Eigen::VectorXd tph_vec = P0 * pk;
    std::vector<int> fine_patch_elems;
    std::vector<int> fine_target_elems;  // tTh = P0*pk > 1
    for (int e = 0; e < tph_vec.size(); ++e) {
        if (tph_vec[e] != 0.0) {
            fine_patch_elems.push_back(e);
            if (tph_vec[e] > 1.0) fine_target_elems.push_back(e);
        }
    }

    // dofph: interior fine vertices in patch
    std::vector<int> fine_vtx_count(Nh, 0);
    for (int e : fine_patch_elems)
        for (int vi : fine.elems[e])
            fine_vtx_count[vi]++;
    std::vector<int> dofph;
    for (int v = 0; v < Nh; ++v)
        if (fine_vtx_count[v] > 0 && fine_vtx_count[v] == nngh[v])
            dofph.push_back(v);

    // Build map: global vertex → local index in dofph
    std::vector<int> dofph_map(Nh, -1);
    for (size_t i = 0; i < dofph.size(); ++i) dofph_map[dofph[i]] = static_cast<int>(i);

    // ---- 3. Sph: direct CG assembly from patch elements (O(|patch|)) ----
    // Shdg is block-diagonal → Ke = 3×3 block at (3e, 3e+1, 3e+2).
    // Sph(v_i, v_j) = sum_{e in patch} Ke(local_i, local_j) for interior v_i, v_j.
    int Nph = static_cast<int>(dofph.size());
    std::vector<Eigen::Triplet<double>> sph_t;
    for (int e : fine_patch_elems) {
        int dg0 = 3*e;
        double Ke[3][3] = {{0}};
        for (int ci = 0; ci < 3; ++ci)
            for (Eigen::SparseMatrix<double>::InnerIterator it(Shdg, dg0+ci); it; ++it)
                if (it.row() >= dg0 && it.row() < dg0+3)
                    Ke[it.row()-dg0][ci] = it.value();

        for (int i = 0; i < 3; ++i) {
            int li = dofph_map[fine.elems[e][i]]; if (li < 0) continue;
            for (int j = 0; j < 3; ++j) {
                int lj = dofph_map[fine.elems[e][j]]; if (lj < 0) continue;
                if (Ke[i][j] != 0.0) sph_t.emplace_back(li, lj, Ke[i][j]);
            }
        }
    }
    Eigen::SparseMatrix<double> Sph(Nph, Nph);
    Sph.setFromTriplets(sph_t.begin(), sph_t.end());

    // ---- 4. rhsp: direct assembly from target elements (O(|target|)) ----
    std::vector<Eigen::Triplet<double>> rhs_t;
    for (int e : fine_target_elems) {
        int dg0 = 3*e;
        double Ke[3][3] = {{0}};
        for (int ci = 0; ci < 3; ++ci)
            for (Eigen::SparseMatrix<double>::InnerIterator it(Shdg, dg0+ci); it; ++it)
                if (it.row() >= dg0 && it.row() < dg0+3)
                    Ke[it.row()-dg0][ci] = it.value();

        for (int i = 0; i < 3; ++i) {
            int li = dofph_map[fine.elems[e][i]]; if (li < 0) continue;
            for (int j = 0; j < d+1; ++j) {
                double val = 0;
                for (int r = 0; r < 3; ++r)
                    val += Ke[i][r] * P1dg.coeff(dg0+r, dgHidx[k][j]);
                if (val != 0.0) rhs_t.emplace_back(li, j, val);
            }
        }
    }
    Eigen::SparseMatrix<double> rhsp(Nph, d+1);
    rhsp.setFromTriplets(rhs_t.begin(), rhs_t.end());

    // ---- 6. IHp = IH(dofpH, dofph) — simple double loop (IH small NH×Nh) ----
    int nd = static_cast<int>(dofpH.size());
    Eigen::SparseMatrix<double> IHp(nd, Nph);
    for (int i = 0; i < nd; ++i) {
        int row = dofpH[i];
        for (int j = 0; j < Nph; ++j) {
            double v = IH.coeff(row, dofph[j]);
            if (v != 0.0) IHp.insert(i, j) = v;
        }
    }

    // ---- 7. Solve X = Sph \ [IHp', rhsp] ----
    Eigen::MatrixXd RHS(Nph, nd + d + 1);
    RHS.setZero();
    // IHp': Nph × nd
    for (int k_ih = 0; k_ih < IHp.outerSize(); ++k_ih)
        for (Eigen::SparseMatrix<double>::InnerIterator it(IHp, k_ih); it; ++it)
            RHS(it.col(), it.row()) = it.value();
    // rhsp: Nph × (d+1)
    for (int k_rh = 0; k_rh < rhsp.outerSize(); ++k_rh)
        for (Eigen::SparseMatrix<double>::InnerIterator it(rhsp, k_rh); it; ++it)
            RHS(it.row(), nd + it.col()) = it.value();

    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt(Sph);
    Eigen::MatrixXd X(Nph, nd+d+1);
    for (int jj = 0; jj < nd+d+1; ++jj) X.col(jj) = llt.solve(RHS.col(jj));

    // ---- 8. mu = (IHp * X1) \ (IHp * X2) ----
    Eigen::MatrixXd IHp_dense = IHp;
    Eigen::MatrixXd X1 = X.leftCols(nd);
    Eigen::MatrixXd X2 = X.rightCols(d+1);
    Eigen::MatrixXd mu = (IHp_dense * X1).colPivHouseholderQr().solve(IHp_dense * X2);

    // ---- 9. Store ----
    Eigen::SparseMatrix<double> CTk(Nh, d+1);
    for (int i = 0; i < Nph; ++i) {
        int global_row = dofph[i];
        for (int j = 0; j < d+1; ++j) {
            double val = X2(i, j);
            for (int r = 0; r < nd; ++r) val -= X1(i, r) * mu(r, j);
            if (std::abs(val) > 1e-15) CTk.insert(global_row, j) = val;
        }
    }

    return CTk;
}

} // namespace lod2d

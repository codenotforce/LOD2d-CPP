#include "lod/corrector.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <set>
#include <algorithm>

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

    // ---- 3. dofphdg: DG dof indices for fine patch elements ----
    std::vector<int> dofphdg;
    for (int e : fine_patch_elems)
        for (int i = 0; i < 3; ++i)
            dofphdg.push_back(dghidx[e][i]);

    // ---- 4. Sph = cg2dgh(dofphdg, dofph)' * Shdg(dofphdg, dofphdg) * cg2dgh(dofphdg, dofph) ----
    int Nph = static_cast<int>(dofph.size());
    int Ndg_ph = static_cast<int>(dofphdg.size());

    // Build cg2dghk: Ndg_ph × Nph  (sparse, 1 at (dg_dof, local_cg_dof))
    std::vector<Eigen::Triplet<double>> cgk_t;
    for (size_t idx = 0; idx < fine_patch_elems.size(); ++idx) {
        int e = fine_patch_elems[idx];
        for (int i = 0; i < 3; ++i) {
            int v = fine.elems[e][i];
            int loc = dofph_map[v];
            if (loc >= 0) cgk_t.emplace_back(static_cast<int>(3*idx + i), loc, 1.0);
        }
    }
    Eigen::SparseMatrix<double> cg2dghk(Ndg_ph, Nph);
    cg2dghk.setFromTriplets(cgk_t.begin(), cgk_t.end());

    // Shdg_patch = Shdg(dofphdg, dofphdg)
    // Direct sparse submatrix extraction using triplet iteration
    std::vector<Eigen::Triplet<double>> s_sub_t;
    // Build index set for fast lookup
    std::set<int> dg_set(dofphdg.begin(), dofphdg.end());
    for (int k_shdg = 0; k_shdg < Shdg.outerSize(); ++k_shdg) {
        if (dg_set.count(k_shdg) == 0) continue;
        for (Eigen::SparseMatrix<double>::InnerIterator it(Shdg, k_shdg); it; ++it) {
            if (dg_set.count(static_cast<int>(it.row())) > 0)
                s_sub_t.emplace_back(it.row(), it.col(), it.value());
        }
    }
    // Map global DG indices → local (0..Ndg_ph-1)
    std::vector<int> dg_map(*std::max_element(dofphdg.begin(), dofphdg.end()) + 1, -1);
    for (size_t i = 0; i < dofphdg.size(); ++i) dg_map[dofphdg[i]] = static_cast<int>(i);
    std::vector<Eigen::Triplet<double>> s_loc_t;
    for (const auto &t : s_sub_t) {
        int li = dg_map[static_cast<int>(t.row())];
        int lj = dg_map[static_cast<int>(t.col())];
        if (li >= 0 && lj >= 0) s_loc_t.emplace_back(li, lj, t.value());
    }
    Eigen::SparseMatrix<double> S_sub(Ndg_ph, Ndg_ph);
    S_sub.setFromTriplets(s_loc_t.begin(), s_loc_t.end());

    Eigen::SparseMatrix<double> Sph = cg2dghk.transpose() * S_sub * cg2dghk;

    // ---- 5. rhsp ----
    // dofThdg: DG dofs of target elements
    std::vector<int> dofThdg;
    for (int e : fine_target_elems)
        for (int i = 0; i < 3; ++i)
            dofThdg.push_back(dghidx[e][i]);

    int Ndg_T = static_cast<int>(dofThdg.size());

    // cg2dghT: Ndg_T × Nph
    std::vector<Eigen::Triplet<double>> cgT_t;
    for (size_t idx = 0; idx < fine_target_elems.size(); ++idx) {
        int e = fine_target_elems[idx];
        for (int i = 0; i < 3; ++i) {
            int v = fine.elems[e][i];
            int loc = dofph_map[v];
            if (loc >= 0) cgT_t.emplace_back(static_cast<int>(3*idx + i), loc, 1.0);
        }
    }
    Eigen::SparseMatrix<double> cg2dghT(Ndg_T, Nph);
    cg2dghT.setFromTriplets(cgT_t.begin(), cgT_t.end());

    // S_th = Shdg(dofThdg, dofThdg)
    std::set<int> dgT_set(dofThdg.begin(), dofThdg.end());
    std::vector<Eigen::Triplet<double>> sT_t;
    for (int k_shdg = 0; k_shdg < Shdg.outerSize(); ++k_shdg) {
        if (dgT_set.count(k_shdg) == 0) continue;
        for (Eigen::SparseMatrix<double>::InnerIterator it(Shdg, k_shdg); it; ++it) {
            if (dgT_set.count(static_cast<int>(it.row())) > 0)
                sT_t.emplace_back(it.row(), it.col(), it.value());
        }
    }
    std::vector<int> dgT_map(*std::max_element(dofThdg.begin(), dofThdg.end()) + 1, -1);
    for (size_t i = 0; i < dofThdg.size(); ++i) dgT_map[dofThdg[i]] = static_cast<int>(i);
    std::vector<Eigen::Triplet<double>> sT_loc;
    for (const auto &t : sT_t) {
        int li = dgT_map[static_cast<int>(t.row())];
        int lj = dgT_map[static_cast<int>(t.col())];
        if (li >= 0 && lj >= 0) sT_loc.emplace_back(li, lj, t.value());
    }
    Eigen::SparseMatrix<double> S_T(Ndg_T, Ndg_T);
    S_T.setFromTriplets(sT_loc.begin(), sT_loc.end());

    // P1dg_block: P1dg(dofThdg, dgHidx(k,:))
    Eigen::SparseMatrix<double> P1dg_block(Ndg_T, d+1);
    for (size_t i = 0; i < dofThdg.size(); ++i) {
        for (int j = 0; j < d+1; ++j) {
            double val = P1dg.coeff(dofThdg[i], dgHidx[k][j]);
            if (val != 0.0) P1dg_block.insert(static_cast<int>(i), j) = val;
        }
    }

    Eigen::SparseMatrix<double> rhsp = cg2dghT.transpose() * S_T * P1dg_block;

    // ---- 6. IHp = IH(dofpH, dofph) ----
    int nd = static_cast<int>(dofpH.size());
    Eigen::SparseMatrix<double> IHp(nd, Nph);
    for (int i = 0; i < nd; ++i) {
        int row = dofpH[i];
        for (const auto &col : dofph) {
            double v = IH.coeff(row, col);
            if (v != 0.0) IHp.insert(i, dofph_map[col]) = v;
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

    Eigen::MatrixXd X = Sph.toDense().ldlt().solve(RHS);

    // ---- 8. mu = (IHp * X1) \ (IHp * X2) ----
    Eigen::MatrixXd IHp_dense = IHp.toDense();
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

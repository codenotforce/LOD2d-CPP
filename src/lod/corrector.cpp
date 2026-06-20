#include "lod/corrector.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <set>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include "solver/cholmod_wrapper.h"

namespace lod2d {

FineElementChildren build_fine_element_children(
    const Eigen::SparseMatrix<double> &P0,
    int coarse_element_count) {
    FineElementChildren children(coarse_element_count);
    for (int coarse_elem = 0; coarse_elem < P0.outerSize(); ++coarse_elem) {
        if (coarse_elem >= coarse_element_count) break;
        for (Eigen::SparseMatrix<double>::InnerIterator it(P0, coarse_elem); it; ++it) {
            if (it.value() != 0.0)
                children[coarse_elem].push_back(static_cast<int>(it.row()));
        }
    }
    return children;
}

InterpolationRows build_interpolation_rows(
    const Eigen::SparseMatrix<double> &IH,
    int coarse_vertex_count) {
    InterpolationRows rows(coarse_vertex_count);
    for (int col = 0; col < IH.outerSize(); ++col) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(IH, col); it; ++it) {
            const int row = static_cast<int>(it.row());
            if (row >= 0 && row < coarse_vertex_count && it.value() != 0.0)
                rows[row].emplace_back(static_cast<int>(it.col()), it.value());
        }
    }
    return rows;
}

CorrectorEntries
compute_corrector_entries(int k,
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
    int d,
    CorrectorSolver solver,
    const ElementStiffnessBlocks *element_stiffness,
    const FineElementChildren *fine_element_children,
    const InterpolationRows *interpolation_rows) {
    if (element_stiffness && element_stiffness->size() != fine.elems.size())
        throw std::invalid_argument("element_stiffness size must match fine element count");
    if (fine_element_children && fine_element_children->size() != coarse.elems.size())
        throw std::invalid_argument("fine_element_children size must match coarse element count");
    if (interpolation_rows && interpolation_rows->size() != coarse.nodes.size())
        throw std::invalid_argument("interpolation_rows size must match coarse vertex count");

    // ---- 1. Coarse patch DOFs and fine patch elements ----
    std::vector<int> tpH_list;
    std::vector<int> fine_patch_elems;
    std::vector<int> fine_target_elems;  // tTh = P0*patch(:,k) > 1

    if (fine_element_children) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(patch, k); it; ++it) {
            if (it.value() == 0.0) continue;
            const int coarse_elem = static_cast<int>(it.row());
            tpH_list.push_back(coarse_elem);

            const auto &children = (*fine_element_children)[coarse_elem];
            fine_patch_elems.insert(fine_patch_elems.end(), children.begin(), children.end());
            if (it.value() > 1.0)
                fine_target_elems.insert(fine_target_elems.end(), children.begin(), children.end());
        }
    } else {
        Eigen::VectorXd pk = patch.col(k);
        for (int j = 0; j < pk.size(); ++j)
            if (pk[j] != 0.0) tpH_list.push_back(j);

        Eigen::VectorXd tph_vec = P0 * pk;
        for (int e = 0; e < tph_vec.size(); ++e) {
            if (tph_vec[e] != 0.0) {
                fine_patch_elems.push_back(e);
                if (tph_vec[e] > 1.0) fine_target_elems.push_back(e);
            }
        }
    }

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
    // dofph: interior fine vertices in patch
    thread_local std::vector<int> fine_vtx_count;
    thread_local std::vector<int> fine_vtx_seen;
    thread_local int fine_vtx_stamp = 0;
    if (static_cast<int>(fine_vtx_count.size()) != Nh) {
        fine_vtx_count.assign(Nh, 0);
        fine_vtx_seen.assign(Nh, 0);
        fine_vtx_stamp = 0;
    }
    if (fine_vtx_stamp == std::numeric_limits<int>::max()) {
        std::fill(fine_vtx_seen.begin(), fine_vtx_seen.end(), 0);
        fine_vtx_stamp = 0;
    }
    ++fine_vtx_stamp;

    std::vector<int> touched_fine_vertices;
    touched_fine_vertices.reserve(3 * fine_patch_elems.size());
    for (int e : fine_patch_elems) {
        for (int vi : fine.elems[e]) {
            if (fine_vtx_seen[vi] != fine_vtx_stamp) {
                fine_vtx_seen[vi] = fine_vtx_stamp;
                fine_vtx_count[vi] = 0;
                touched_fine_vertices.push_back(vi);
            }
            fine_vtx_count[vi]++;
        }
    }
    std::vector<int> dofph;
    dofph.reserve(touched_fine_vertices.size());
    for (int v : touched_fine_vertices)
        if (fine_vtx_count[v] > 0 && fine_vtx_count[v] == nngh[v])
            dofph.push_back(v);

    // Build map: global vertex ??local index in dofph
    thread_local std::vector<int> dofph_map;
    thread_local std::vector<int> dofph_seen;
    thread_local int dofph_stamp = 0;
    if (static_cast<int>(dofph_map.size()) != Nh) {
        dofph_map.assign(Nh, -1);
        dofph_seen.assign(Nh, 0);
        dofph_stamp = 0;
    }
    if (dofph_stamp == std::numeric_limits<int>::max()) {
        std::fill(dofph_seen.begin(), dofph_seen.end(), 0);
        dofph_stamp = 0;
    }
    ++dofph_stamp;
    for (size_t i = 0; i < dofph.size(); ++i) {
        const int v = dofph[i];
        dofph_map[v] = static_cast<int>(i);
        dofph_seen[v] = dofph_stamp;
    }
    auto local_fine_dof = [&](int v) -> int {
        return dofph_seen[v] == dofph_stamp ? dofph_map[v] : -1;
    };

    // ---- 3. Sph: direct CG assembly from patch elements (O(|patch|)) ----
    // Shdg is block-diagonal ??Ke = 3?3 block at (3e, 3e+1, 3e+2).
    // Sph(v_i, v_j) = sum_{e in patch} Ke(local_i, local_j) for interior v_i, v_j.
    int Nph = static_cast<int>(dofph.size());
    thread_local std::vector<Eigen::Triplet<double>> sph_t;
    sph_t.clear();
    if (sph_t.capacity() < 9 * fine_patch_elems.size())
        sph_t.reserve(9 * fine_patch_elems.size());
    for (int e : fine_patch_elems) {
        Eigen::Matrix3d Ke;
        if (element_stiffness) {
            Ke = (*element_stiffness)[e];
        } else {
            int dg0 = 3*e;
            Ke.setZero();
            for (int ci = 0; ci < 3; ++ci)
                for (Eigen::SparseMatrix<double>::InnerIterator it(Shdg, dg0+ci); it; ++it)
                    if (it.row() >= dg0 && it.row() < dg0+3)
                        Ke(it.row()-dg0, ci) = it.value();
        }

        for (int i = 0; i < 3; ++i) {
            int li = local_fine_dof(fine.elems[e][i]); if (li < 0) continue;
            for (int j = 0; j < 3; ++j) {
                int lj = local_fine_dof(fine.elems[e][j]); if (lj < 0) continue;
                if (Ke(i, j) != 0.0) sph_t.emplace_back(li, lj, Ke(i, j));
            }
        }
    }
    Eigen::SparseMatrix<double> Sph(Nph, Nph);
    Sph.setFromTriplets(sph_t.begin(), sph_t.end());

    // ---- 4. rhsp: direct assembly from target elements (O(|target|)) ----
    thread_local std::vector<Eigen::Triplet<double>> rhs_t;
    rhs_t.clear();
    if (rhs_t.capacity() < 3 * (d + 1) * fine_target_elems.size())
        rhs_t.reserve(3 * (d + 1) * fine_target_elems.size());
    for (int e : fine_target_elems) {
        int dg0 = 3*e;
        Eigen::Matrix3d Ke;
        if (element_stiffness) {
            Ke = (*element_stiffness)[e];
        } else {
            Ke.setZero();
            for (int ci = 0; ci < 3; ++ci)
                for (Eigen::SparseMatrix<double>::InnerIterator it(Shdg, dg0+ci); it; ++it)
                    if (it.row() >= dg0 && it.row() < dg0+3)
                        Ke(it.row()-dg0, ci) = it.value();
        }

        for (int i = 0; i < 3; ++i) {
            int li = local_fine_dof(fine.elems[e][i]); if (li < 0) continue;
            for (int j = 0; j < d+1; ++j) {
                double val = 0;
                for (int r = 0; r < 3; ++r)
                    val += Ke(i, r) * P1dg.coeff(dg0+r, dgHidx[k][j]);
                if (val != 0.0) rhs_t.emplace_back(li, j, val);
            }
        }
    }
    Eigen::SparseMatrix<double> rhsp(Nph, d+1);
    rhsp.setFromTriplets(rhs_t.begin(), rhs_t.end());

    // ---- 6. IHp = IH(dofpH, dofph) ??dense local block ----
    int nd = static_cast<int>(dofpH.size());
    Eigen::MatrixXd IHp_dense = Eigen::MatrixXd::Zero(nd, Nph);
    for (int i = 0; i < nd; ++i) {
        int row = dofpH[i];
        if (interpolation_rows) {
            for (const auto &[fine_col, value] : (*interpolation_rows)[row]) {
                const int local_col = local_fine_dof(fine_col);
                if (local_col >= 0)
                    IHp_dense(i, local_col) = value;
            }
        } else {
            for (int j = 0; j < Nph; ++j) {
                double v = IH.coeff(row, dofph[j]);
                if (v != 0.0) IHp_dense(i, j) = v;
            }
        }
    }

    // ---- 7. Solve X = Sph \ [IHp', rhsp] ----
    Eigen::MatrixXd RHS(Nph, nd + d + 1);
    RHS.setZero();
    RHS.leftCols(nd) = IHp_dense.transpose();
    // rhsp: Nph ? (d+1)
    for (int k_rh = 0; k_rh < rhsp.outerSize(); ++k_rh)
        for (Eigen::SparseMatrix<double>::InnerIterator it(rhsp, k_rh); it; ++it)
            RHS(it.row(), nd + it.col()) = it.value();

    Eigen::MatrixXd X(Nph, nd+d+1);
    if (solver == CorrectorSolver::Cholmod) {
        X = solve_cholmod(Sph, RHS);
    } else if (solver == CorrectorSolver::CholmodCached) {
        X = solve_cholmod_cached(Sph, RHS);
    } else {
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt(Sph);
        X = llt.solve(RHS);
    }

    // ---- 8. mu = (IHp * X1) \ (IHp * X2) ----
    Eigen::MatrixXd X1 = X.leftCols(nd);
    Eigen::MatrixXd X2 = X.rightCols(d+1);
    Eigen::MatrixXd mu = (IHp_dense * X1).colPivHouseholderQr().solve(IHp_dense * X2);

    // ---- 9. Store ----
    CorrectorEntries entries;
    entries.reserve(static_cast<size_t>(Nph) * static_cast<size_t>(d + 1));
    for (int i = 0; i < Nph; ++i) {
        int global_row = dofph[i];
        for (int j = 0; j < d+1; ++j) {
            double val = X2(i, j);
            for (int r = 0; r < nd; ++r) val -= X1(i, r) * mu(r, j);
            if (std::abs(val) > 1e-15) entries.push_back({global_row, j, val});
        }
    }
    return entries;
}

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
    int d,
    CorrectorSolver solver,
    const ElementStiffnessBlocks *element_stiffness,
    const FineElementChildren *fine_element_children,
    const InterpolationRows *interpolation_rows) {
    CorrectorEntries entries = compute_corrector_entries(k, patch, coarse, NH, nngH,
        P0, fine, Nh, nngh, dghidx, cg2dgh, Shdg, P1dg, dgHidx, IH, d, solver,
        element_stiffness, fine_element_children, interpolation_rows);

    std::vector<Eigen::Triplet<double>> ctk_t;
    ctk_t.reserve(entries.size());
    for (const auto &entry : entries)
        ctk_t.emplace_back(entry.row, entry.col, entry.value);

    Eigen::SparseMatrix<double> CTk(Nh, d+1);
    CTk.setFromTriplets(ctk_t.begin(), ctk_t.end());

    return CTk;
}

} // namespace lod2d


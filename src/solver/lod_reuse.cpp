#include "solver/lod_reuse.h"
#include <stdexcept>
#include <utility>

namespace lod2d {

LodReusableSystem::LodReusableSystem(
    Eigen::SparseMatrix<double> G,
    Eigen::SparseMatrix<double> Sh,
    Eigen::SparseMatrix<double> Mh,
    Eigen::SparseMatrix<double> P_node,
    int coarse_node_count,
    const std::vector<int> &coarse_dirichlet)
    : NH_(coarse_node_count),
      Nh_(static_cast<int>(G.rows())),
      G_(std::move(G)),
      Mh_(std::move(Mh)),
      P_node_(std::move(P_node)) {
    if (G_.cols() != NH_)
        throw std::invalid_argument("G column count must match coarse_node_count");
    if (Mh_.rows() != Nh_ || Mh_.cols() != Nh_)
        throw std::invalid_argument("Mh dimensions must match fine node count");
    if (P_node_.rows() != Nh_ || P_node_.cols() != NH_)
        throw std::invalid_argument("P_node dimensions must be Nh x NH");

    std::vector<int> dofH_map(NH_, -1);
    std::vector<char> is_dirH(NH_, false);
    for (int dv : coarse_dirichlet) {
        if (dv >= 0 && dv < NH_) is_dirH[dv] = true;
    }
    for (int i = 0; i < NH_; ++i) {
        if (!is_dirH[i]) {
            dofH_map[i] = static_cast<int>(dofH_.size());
            dofH_.push_back(i);
        }
    }

    std::vector<Eigen::Triplet<double>> g0_t;
    g0_t.reserve(G_.nonZeros());
    for (int col = 0; col < G_.outerSize(); ++col) {
        const int local_col = dofH_map[col];
        if (local_col < 0) continue;
        for (Eigen::SparseMatrix<double>::InnerIterator it(G_, col); it; ++it)
            g0_t.emplace_back(it.row(), local_col, it.value());
    }
    G0_.resize(Nh_, static_cast<int>(dofH_.size()));
    G0_.setFromTriplets(g0_t.begin(), g0_t.end());

    Eigen::SparseMatrix<double> SHLOD0 = G0_.transpose() * Sh * G0_;
    coarse_factor_.compute(SHLOD0);
    if (coarse_factor_.info() != Eigen::Success)
        throw std::runtime_error("LOD coarse factorization failed");
}

LodReuseSolution LodReusableSystem::solve_from_coarse_values(const Eigen::VectorXd &f_coarse) const {
    if (f_coarse.size() != NH_)
        throw std::invalid_argument("coarse RHS size must match coarse node count");
    return solve_from_fine_values(P_node_ * f_coarse);
}

LodReuseSolution LodReusableSystem::solve_from_fine_values(const Eigen::VectorXd &f_fine) const {
    if (f_fine.size() != Nh_)
        throw std::invalid_argument("fine RHS size must match fine node count");

    Eigen::VectorXd rhs = G0_.transpose() * (Mh_ * f_fine);
    Eigen::VectorXd uf = coarse_factor_.solve(rhs);
    if (coarse_factor_.info() != Eigen::Success)
        throw std::runtime_error("LOD coarse solve failed");

    LodReuseSolution out;
    out.uH = Eigen::VectorXd::Zero(NH_);
    for (size_t j = 0; j < dofH_.size(); ++j)
        out.uH(dofH_[j]) = uf(static_cast<int>(j));
    out.uHms = G_ * out.uH;
    return out;
}

} // namespace lod2d

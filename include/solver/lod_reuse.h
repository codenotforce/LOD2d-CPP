#pragma once
#include "mesh/types.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

namespace lod2d {

struct LodReuseSolution {
    Eigen::VectorXd uH;
    Eigen::VectorXd uHms;
};

class LodReusableSystem {
public:
    LodReusableSystem() = default;

    LodReusableSystem(
        Eigen::SparseMatrix<double> G,
        Eigen::SparseMatrix<double> Sh,
        Eigen::SparseMatrix<double> Mh,
        Eigen::SparseMatrix<double> P_node,
        int coarse_node_count,
        const std::vector<int> &coarse_dirichlet);

    int coarse_node_count() const { return NH_; }
    int fine_node_count() const { return Nh_; }
    int free_coarse_dof_count() const { return static_cast<int>(dofH_.size()); }

    LodReuseSolution solve_from_coarse_values(const Eigen::VectorXd &f_coarse) const;
    LodReuseSolution solve_from_fine_values(const Eigen::VectorXd &f_fine) const;

private:
    int NH_ = 0;
    int Nh_ = 0;
    std::vector<int> dofH_;
    Eigen::SparseMatrix<double> G_;
    Eigen::SparseMatrix<double> G0_;
    Eigen::SparseMatrix<double> Mh_;
    Eigen::SparseMatrix<double> P_node_;
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> coarse_factor_;
};

} // namespace lod2d

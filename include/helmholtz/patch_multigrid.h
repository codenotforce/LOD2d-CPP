#pragma once

#include "helmholtz/patch_system.h"

#include <Eigen/SparseLU>
#include <vector>

namespace lod2d::helmholtz {

class HelmholtzPatchVcycle {
public:
    HelmholtzPatchVcycle(
        const HelmholtzPatchSystem &system,
        double epsilon,
        int pre_smooth,
        int post_smooth,
        int coarse_max_dofs,
        double jacobi_weight);
    HelmholtzPatchVcycle(
        const ComplexSparseMatrix &shifted_operator,
        const std::vector<Eigen::SparseMatrix<double>> &geometric_prolongations,
        int pre_smooth,
        int post_smooth,
        int coarse_max_dofs,
        double jacobi_weight);

    ComplexVector apply(const ComplexVector &right_hand_side) const;
    double relative_residual(const ComplexVector &right_hand_side) const;

    int levels() const { return static_cast<int>(operators_.size()); }
    int coarse_dofs() const { return operators_.front().rows(); }
    int finest_dofs() const { return operators_.back().rows(); }

private:
    ComplexVector cycle(int level, const ComplexVector &right_hand_side) const;
    void smooth(
        int level,
        const ComplexVector &right_hand_side,
        int steps,
        ComplexVector &iterate) const;

    std::vector<ComplexSparseMatrix> operators_;
    std::vector<ComplexSparseMatrix> prolongations_;
    std::vector<ComplexVector> inverse_diagonals_;
    int pre_smooth_ = 2;
    int post_smooth_ = 2;
    double jacobi_weight_ = 0.6;
    Eigen::SparseLU<ComplexSparseMatrix> coarse_solver_;
};

} // namespace lod2d::helmholtz

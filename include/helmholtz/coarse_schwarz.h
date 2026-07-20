#pragma once

#include "helmholtz/types.h"
#include "mesh/types.h"

#include <Eigen/SparseLU>
#include <memory>
#include <vector>

namespace lod2d::helmholtz {

struct HelmholtzCoarseSchwarzDiagnostics {
    int subdomains = 0;
    int min_local_dofs = 0;
    int max_local_dofs = 0;
    double partition_unity_error = 0.0;
};

class HelmholtzCoarseRasPreconditioner {
public:
    HelmholtzCoarseRasPreconditioner(
        const ComplexSparseMatrix &coarse_operator,
        const TriMesh &coarse_mesh,
        const Eigen::SparseMatrix<double> &element_patches);

    ComplexVector apply(const ComplexVector &right_hand_side) const;

    const HelmholtzCoarseSchwarzDiagnostics &diagnostics() const {
        return diagnostics_;
    }

private:
    struct Subdomain {
        std::vector<int> global_dofs;
        std::unique_ptr<Eigen::SparseLU<ComplexSparseMatrix>> solver;
    };

    int dimension_ = 0;
    std::vector<Subdomain> subdomains_;
    Eigen::VectorXd injection_weights_;
    HelmholtzCoarseSchwarzDiagnostics diagnostics_;
};

} // namespace lod2d::helmholtz

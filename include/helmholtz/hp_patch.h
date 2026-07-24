#pragma once

#include "helmholtz/hp_interpolation.h"
#include "helmholtz/hp_operators.h"
#include "helmholtz/patch_solver.h"

namespace lod2d::helmholtz {

struct HelmholtzHpPatchSolveResult {
    HelmholtzPatchSystem system;
    HelmholtzPatchSolveResult primal;
    ComplexMatrix adjoint_corrector;
    ComplexMatrix adjoint_multipliers;
};

class HelmholtzHpPatchAssembler {
public:
    HelmholtzHpPatchAssembler(
        const TriMesh &coarse,
        const HpTriSpace &fine_space,
        const Eigen::SparseMatrix<double> &fine_element_prolongation,
        const Eigen::SparseMatrix<double> &patches,
        const HelmholtzHpInterpolation &interpolation,
        const HelmholtzHpOperators &operators);

    HelmholtzPatchSystem assemble(int target) const;
    HelmholtzHpPatchSolveResult solve_direct_saddle(int target) const;
    std::size_t patch_cost(int target) const;
    int patch_count() const {
        return static_cast<int>(coarse_.elems.size());
    }

private:
    const TriMesh &coarse_;
    const HpTriSpace &fine_space_;
    const Eigen::SparseMatrix<double> &patches_;
    const HelmholtzHpInterpolation &interpolation_;
    const HelmholtzHpOperators &operators_;
    std::vector<std::vector<int>> children_;
};

} // namespace lod2d::helmholtz

#pragma once

#include "helmholtz/hp_patch.h"
#include "helmholtz/model.h"

#include <memory>

namespace lod2d::helmholtz {

struct HelmholtzHpProblemConfig {
    int H = 2;
    int h = 6;
    int ell = 2;
    int degree = 1;
    double wavenumber = 4.0;
    double boundary_beta = 1.0;
    HelmholtzPetrovMode mode = HelmholtzPetrovMode::TwoSided;
    TriMesh initial_mesh;
};

struct HelmholtzHpCorrectorDiagnostics {
    double max_primal_residual = 0.0;
    double max_adjoint_residual = 0.0;
    double max_constraint_residual = 0.0;
    int patch_count = 0;
    int boundary_patch_count = 0;
};

class HelmholtzHpLodModel {
public:
    HelmholtzHpLodModel();
    ~HelmholtzHpLodModel();
    HelmholtzHpLodModel(HelmholtzHpLodModel &&) noexcept;
    HelmholtzHpLodModel &operator=(HelmholtzHpLodModel &&) noexcept;
    HelmholtzHpLodModel(const HelmholtzHpLodModel &) = delete;
    HelmholtzHpLodModel &operator=(const HelmholtzHpLodModel &) = delete;

    static HelmholtzHpLodModel build(const HelmholtzHpProblemConfig &config);

    HelmholtzLodSolution solve_load(const ComplexVector &load) const;
    HelmholtzLodSolution solve_source(const ComplexFunction &source) const;
    ComplexVector solve_fine_reference(const ComplexVector &load) const;

    const HelmholtzHpProblemConfig &config() const;
    const HelmholtzProblemData &problem() const;
    const HpTriSpace &fine_space() const;
    const HelmholtzHpInterpolation &interpolation() const;
    const HelmholtzHpOperators &operators() const;
    const HelmholtzHpCorrectorDiagnostics &corrector_diagnostics() const;
    const ComplexSparseMatrix &corrector_matrix() const;
    const ComplexSparseMatrix &adjoint_corrector_matrix() const;
    const ComplexSparseMatrix &trial_basis() const;
    const ComplexSparseMatrix &test_basis() const;
    const ComplexSparseMatrix &coarse_operator() const;

private:
    struct Implementation;
    std::unique_ptr<Implementation> implementation_;
};

} // namespace lod2d::helmholtz

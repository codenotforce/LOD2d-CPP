#pragma once

#include "helmholtz/corrector.h"
#include "lod/patches.h"
#include "lod/quasi_interp.h"
#include "mesh/refine.h"
#include <memory>
#include <vector>

namespace lod2d::helmholtz {

enum class HelmholtzPetrovMode {
    TwoSided,
    CorrectedTestOnly
};

struct HelmholtzProblemConfig {
    int H = 1;
    int h = 3;
    int ell = 1;
    double wavenumber = 2.0;
    double boundary_beta = 1.0;
    HelmholtzPetrovMode mode = HelmholtzPetrovMode::TwoSided;
    bool progress = false;
    TriMesh initial_mesh;
    HelmholtzPatchSolverConfig patch_solver;
    std::vector<double> diffusion;
    std::vector<double> refractive_index;
    QuadraturePolicy quadrature;
    QuadratureContext quadrature_context;
};

struct HelmholtzProblemData {
    TriMesh coarse;
    TriMesh fine;
    Eigen::SparseMatrix<double> coarse_to_fine;
    Eigen::SparseMatrix<double> fine_element_prolongation;
    Eigen::SparseMatrix<double> fine_dg_prolongation;
    Eigen::SparseMatrix<double> quasi_interpolation;
    std::vector<TriMesh> fine_hierarchy_meshes;
    std::vector<Eigen::SparseMatrix<double>> fine_node_level_prolongations;
    std::vector<Eigen::SparseMatrix<double>> fine_element_level_prolongations;
    Eigen::SparseMatrix<double> patches;
    std::vector<int> coarse_element_levels;
    std::vector<int> fine_element_levels;
    int fine_level = -1;
    int max_fine_level = -1;
};

struct HelmholtzLodSolution {
    ComplexVector coarse_coefficients;
    ComplexVector fine_values;
    double petrov_residual = 0.0;
};

struct HelmholtzBuildTimings {
    double mesh_and_interpolation_ms = 0.0;
    double operators_ms = 0.0;
    double correctors_ms = 0.0;
    double basis_and_factorization_ms = 0.0;
    double total_ms = 0.0;
};

TriMesh make_helmholtz_unit_square_mesh();
TriMesh make_helmholtz_l_shape_mesh();

HelmholtzProblemData build_helmholtz_problem_data(
    const TriMesh &initial_mesh,
    int H,
    int h,
    int ell);
HelmholtzProblemData build_adaptive_helmholtz_problem_data(
    const TriMesh &coarse_mesh,
    const std::vector<int> &coarse_element_levels,
    int fine_level,
    int ell);
HelmholtzProblemData build_adaptive_helmholtz_problem_data(
    const TriMesh &coarse_mesh,
    const std::vector<int> &coarse_element_levels,
    const TriMesh &fine_mesh,
    const std::vector<int> &fine_element_levels,
    int ell);

class HelmholtzLodModel {
public:
    HelmholtzLodModel();
    ~HelmholtzLodModel();
    HelmholtzLodModel(HelmholtzLodModel &&) noexcept;
    HelmholtzLodModel &operator=(HelmholtzLodModel &&) noexcept;
    HelmholtzLodModel(const HelmholtzLodModel &) = delete;
    HelmholtzLodModel &operator=(const HelmholtzLodModel &) = delete;

    static HelmholtzLodModel build(const HelmholtzProblemConfig &config);
    static HelmholtzLodModel build_adaptive(
        const HelmholtzProblemConfig &config,
        const TriMesh &coarse_mesh,
        const std::vector<int> &coarse_element_levels);
    static HelmholtzLodModel build_adaptive(
        const HelmholtzProblemConfig &config,
        const TriMesh &coarse_mesh,
        const std::vector<int> &coarse_element_levels,
        const TriMesh &fine_mesh,
        const std::vector<int> &fine_element_levels);

    HelmholtzLodSolution solve_load(const ComplexVector &fine_load) const;
    HelmholtzLodSolution solve_source(const ComplexFunction &source) const;
    ComplexVector solve_fine_reference(const ComplexVector &fine_load) const;

    const HelmholtzProblemConfig &config() const { return config_; }
    const HelmholtzProblemData &problem() const { return problem_; }
    const HelmholtzOperators &operators() const { return operators_; }
    const HelmholtzCorrectorResult &correctors() const { return correctors_; }
    const ComplexSparseMatrix &corrected_trial_basis() const { return corrected_trial_basis_; }
    const ComplexSparseMatrix &corrected_test_basis() const { return corrected_test_basis_; }
    const ComplexSparseMatrix &trial_basis() const { return trial_basis_; }
    const ComplexSparseMatrix &test_basis() const { return test_basis_; }
    const ComplexSparseMatrix &coarse_operator() const { return coarse_operator_; }
    const HelmholtzBuildTimings &build_timings() const { return build_timings_; }

private:
    struct Factorization;

    static HelmholtzLodModel build_with_problem(
        HelmholtzProblemConfig config,
        HelmholtzProblemData problem,
        double mesh_and_interpolation_ms);

    HelmholtzProblemConfig config_;
    HelmholtzProblemData problem_;
    HelmholtzOperators operators_;
    HelmholtzCorrectorResult correctors_;
    ComplexSparseMatrix corrected_trial_basis_;
    ComplexSparseMatrix corrected_test_basis_;
    ComplexSparseMatrix trial_basis_;
    ComplexSparseMatrix test_basis_;
    ComplexSparseMatrix coarse_operator_;
    std::unique_ptr<Factorization> factorization_;
    HelmholtzBuildTimings build_timings_;
};

} // namespace lod2d::helmholtz

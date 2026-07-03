#include "helmholtz/model.h"

#include <Eigen/SparseLU>
#include <algorithm>
#include <chrono>
#include <stdexcept>
#include <utility>

namespace lod2d::helmholtz {
namespace {

Eigen::SparseMatrix<double> build_cg_to_dg(const TriMesh &mesh) {
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(3 * mesh.elems.size());
    for (int element = 0; element < static_cast<int>(mesh.elems.size()); ++element) {
        for (int local = 0; local < 3; ++local)
            triplets.emplace_back(3 * element + local, mesh.elems[element][local], 1.0);
    }
    Eigen::SparseMatrix<double> result(
        3 * static_cast<int>(mesh.elems.size()),
        static_cast<int>(mesh.nodes.size()));
    result.setFromTriplets(triplets.begin(), triplets.end());
    return result;
}

double elapsed_ms(const std::chrono::steady_clock::time_point &start) {
    return std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - start).count();
}

void validate_pure_robin_mesh(const TriMesh &mesh) {
    if (mesh.nodes.empty() || mesh.elems.empty())
        throw std::invalid_argument("Helmholtz initial mesh must not be empty");
    if (!mesh.dirichlet.empty())
        throw std::invalid_argument("pure Robin Helmholtz mesh must not contain Dirichlet nodes");
}

} // namespace

struct HelmholtzLodModel::Factorization {
    Eigen::SparseLU<ComplexSparseMatrix> solver;
};

HelmholtzLodModel::HelmholtzLodModel() = default;
HelmholtzLodModel::~HelmholtzLodModel() = default;
HelmholtzLodModel::HelmholtzLodModel(HelmholtzLodModel &&) noexcept = default;
HelmholtzLodModel &HelmholtzLodModel::operator=(HelmholtzLodModel &&) noexcept = default;

TriMesh make_helmholtz_unit_square_mesh() {
    TriMesh mesh;
    mesh.nodes = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    mesh.elems = {{0, 1, 3}, {2, 3, 1}};
    return mesh;
}

HelmholtzProblemData build_helmholtz_problem_data(
    const TriMesh &initial_mesh,
    int H,
    int h,
    int ell) {
    validate_pure_robin_mesh(initial_mesh);
    if (H < 0 || h < H)
        throw std::invalid_argument("Helmholtz refinement levels must satisfy 0 <= H <= h");
    if (ell < 0) throw std::invalid_argument("Helmholtz oversampling level must be nonnegative");

    RefineOutput coarse_output = refine_mesh_nvb(initial_mesh, H);
    RefineOutput fine_output = refine_mesh_nvb(coarse_output.mesh, h - H);

    HelmholtzProblemData problem;
    problem.coarse = std::move(coarse_output.mesh);
    problem.fine = std::move(fine_output.mesh);
    problem.coarse_to_fine = std::move(fine_output.P_node);
    problem.fine_element_prolongation = std::move(fine_output.P_elem);
    problem.fine_dg_prolongation = std::move(fine_output.P_dg);
    problem.patches = build_patches(problem.coarse, ell);

    const Eigen::SparseMatrix<double> cg_to_dg = build_cg_to_dg(problem.fine);
    problem.quasi_interpolation = build_quasi_interp(
        problem.coarse,
        problem.fine,
        problem.fine_dg_prolongation,
        cg_to_dg,
        static_cast<int>(problem.fine.nodes.size()),
        static_cast<int>(problem.coarse.nodes.size()));

    const Eigen::MatrixXd projection_check = Eigen::MatrixXd(
        problem.quasi_interpolation * problem.coarse_to_fine);
    const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(
        projection_check.rows(), projection_check.cols());
    if ((projection_check - identity).norm() > 1e-9)
        throw std::runtime_error("complex quasi-interpolation does not reproduce the coarse P1 space");
    return problem;
}

HelmholtzLodModel HelmholtzLodModel::build(const HelmholtzProblemConfig &config) {
    if (config.wavenumber <= 0.0)
        throw std::invalid_argument("Helmholtz wavenumber must be positive");
    HelmholtzProblemConfig resolved = config;
    if (resolved.initial_mesh.nodes.empty())
        resolved.initial_mesh = make_helmholtz_unit_square_mesh();

    const auto total_start = std::chrono::steady_clock::now();
    HelmholtzLodModel model;
    model.config_ = resolved;
    auto stage_start = std::chrono::steady_clock::now();
    model.problem_ = build_helmholtz_problem_data(
        resolved.initial_mesh, resolved.H, resolved.h, resolved.ell);
    model.build_timings_.mesh_and_interpolation_ms = elapsed_ms(stage_start);
    stage_start = std::chrono::steady_clock::now();
    model.operators_ = assemble_helmholtz_operators(
        model.problem_.fine,
        resolved.wavenumber,
        resolved.diffusion,
        resolved.refractive_index,
        resolved.boundary_beta);
    model.build_timings_.operators_ms = elapsed_ms(stage_start);
    stage_start = std::chrono::steady_clock::now();
    model.correctors_ = build_helmholtz_correctors(
        model.problem_.coarse,
        model.problem_.fine,
        model.problem_.fine_element_prolongation,
        model.problem_.fine_dg_prolongation,
        model.problem_.quasi_interpolation,
        model.problem_.patches,
        model.operators_.element_blocks);
    model.build_timings_.correctors_ms = elapsed_ms(stage_start);

    const auto &diagnostics = model.correctors_.diagnostics;
    if (diagnostics.max_primal_residual > 1e-8
        || diagnostics.max_adjoint_residual > 1e-8
        || diagnostics.max_constraint_residual > 1e-8) {
        throw std::runtime_error("Helmholtz corrector residual exceeded the correctness threshold");
    }

    stage_start = std::chrono::steady_clock::now();
    const int fine_node_count = static_cast<int>(model.problem_.fine.nodes.size());
    model.corrected_trial_basis_ = build_helmholtz_corrected_basis(
        model.problem_.coarse_to_fine,
        model.problem_.coarse,
        fine_node_count,
        model.correctors_.primal);
    model.corrected_test_basis_ = build_helmholtz_corrected_basis(
        model.problem_.coarse_to_fine,
        model.problem_.coarse,
        fine_node_count,
        model.correctors_.adjoint);
    model.test_basis_ = model.corrected_test_basis_;
    if (resolved.mode == HelmholtzPetrovMode::TwoSided)
        model.trial_basis_ = model.corrected_trial_basis_;
    else
        model.trial_basis_ = model.problem_.coarse_to_fine.cast<Complex>();

    model.coarse_operator_ = model.test_basis_.adjoint()
                           * model.operators_.system
                           * model.trial_basis_;
    model.coarse_operator_.prune(Complex(0.0, 0.0), 1e-14);
    model.coarse_operator_.makeCompressed();

    model.factorization_ = std::make_unique<Factorization>();
    model.factorization_->solver.analyzePattern(model.coarse_operator_);
    model.factorization_->solver.factorize(model.coarse_operator_);
    if (model.factorization_->solver.info() != Eigen::Success)
        throw std::runtime_error("Helmholtz Petrov-Galerkin coarse factorization failed");
    model.build_timings_.basis_and_factorization_ms = elapsed_ms(stage_start);
    model.build_timings_.total_ms = elapsed_ms(total_start);
    return model;
}

HelmholtzLodSolution HelmholtzLodModel::solve_load(const ComplexVector &fine_load) const {
    if (!factorization_)
        throw std::logic_error("Helmholtz LOD model has not been built");
    if (fine_load.size() != static_cast<int>(problem_.fine.nodes.size()))
        throw std::invalid_argument("Helmholtz fine load size does not match the model");

    const ComplexVector coarse_rhs = test_basis_.adjoint() * fine_load;
    HelmholtzLodSolution solution;
    solution.coarse_coefficients = factorization_->solver.solve(coarse_rhs);
    if (factorization_->solver.info() != Eigen::Success
        || !solution.coarse_coefficients.allFinite())
        throw std::runtime_error("Helmholtz Petrov-Galerkin coarse solve failed");
    solution.fine_values = trial_basis_ * solution.coarse_coefficients;
    const ComplexVector residual = test_basis_.adjoint()
        * (operators_.system * solution.fine_values - fine_load);
    solution.petrov_residual = residual.norm() / std::max(1.0, coarse_rhs.norm());
    return solution;
}

HelmholtzLodSolution HelmholtzLodModel::solve_source(const ComplexFunction &source) const {
    return solve_load(assemble_helmholtz_load(problem_.fine, source));
}

ComplexVector HelmholtzLodModel::solve_fine_reference(const ComplexVector &fine_load) const {
    return solve_helmholtz_fem(operators_, fine_load);
}

} // namespace lod2d::helmholtz

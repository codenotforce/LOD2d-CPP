#include "helmholtz/model.h"
#include "helmholtz/adaptive/hierarchy.h"

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

namespace {

HelmholtzProblemData finish_problem_data(
    TriMesh coarse,
    RefineOutput fine_output,
    std::vector<int> coarse_levels,
    int fine_level,
    int ell,
    std::vector<TriMesh> hierarchy_meshes,
    std::vector<Eigen::SparseMatrix<double>> node_level_prolongations,
    std::vector<Eigen::SparseMatrix<double>> element_level_prolongations) {
    if (ell < 0) throw std::invalid_argument("Helmholtz oversampling level must be nonnegative");

    HelmholtzProblemData problem;
    problem.coarse = std::move(coarse);
    problem.fine = std::move(fine_output.mesh);
    problem.coarse_to_fine = std::move(fine_output.P_node);
    problem.fine_element_prolongation = std::move(fine_output.P_elem);
    problem.fine_dg_prolongation = std::move(fine_output.P_dg);
    problem.patches = build_patches(problem.coarse, ell);
    problem.coarse_element_levels = std::move(coarse_levels);
    problem.fine_level = fine_level;
    problem.fine_hierarchy_meshes = std::move(hierarchy_meshes);
    problem.fine_node_level_prolongations = std::move(node_level_prolongations);
    problem.fine_element_level_prolongations = std::move(element_level_prolongations);

    const Eigen::SparseMatrix<double> coarse_cg_to_dg = build_cg_to_dg(problem.coarse);
    const Eigen::SparseMatrix<double> cg_to_dg = build_cg_to_dg(problem.fine);
    const Eigen::SparseMatrix<double> dg_from_coarse =
        problem.fine_dg_prolongation * coarse_cg_to_dg;
    const Eigen::SparseMatrix<double> dg_from_fine =
        cg_to_dg * problem.coarse_to_fine;
    if ((dg_from_coarse - dg_from_fine).norm() > 1e-10)
        throw std::runtime_error("DG and nodal prolongations represent different coarse P1 functions");
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
struct UniformFineHierarchy {
    RefineOutput accumulated;
    std::vector<TriMesh> meshes;
    std::vector<Eigen::SparseMatrix<double>> node_prolongations;
    std::vector<Eigen::SparseMatrix<double>> element_prolongations;
};

UniformFineHierarchy build_uniform_fine_hierarchy(
    const TriMesh &coarse,
    int levels) {
    UniformFineHierarchy result;
    result.meshes.push_back(coarse);
    const int node_count = static_cast<int>(coarse.nodes.size());
    const int element_count = static_cast<int>(coarse.elems.size());
    result.accumulated.P_node.resize(node_count, node_count);
    result.accumulated.P_node.setIdentity();
    result.accumulated.P_elem.resize(element_count, element_count);
    result.accumulated.P_elem.setIdentity();
    result.accumulated.P_dg.resize(3 * element_count, 3 * element_count);
    result.accumulated.P_dg.setIdentity();

    TriMesh current = coarse;
    for (int level = 0; level < levels; ++level) {
        RefineOutput refined = refine_nvb(current);
        result.accumulated.P_node =
            refined.P_node * result.accumulated.P_node;
        result.accumulated.P_elem =
            refined.P_elem * result.accumulated.P_elem;
        result.accumulated.P_dg =
            refined.P_dg * result.accumulated.P_dg;
        result.node_prolongations.push_back(refined.P_node);
        result.element_prolongations.push_back(refined.P_elem);
        current = std::move(refined.mesh);
        result.meshes.push_back(current);
    }
    result.accumulated.mesh = std::move(current);
    return result;
}

} // namespace

HelmholtzProblemData build_helmholtz_problem_data(
    const TriMesh &initial_mesh,
    int H,
    int h,
    int ell) {
    validate_pure_robin_mesh(initial_mesh);
    if (H < 0 || h < H)
        throw std::invalid_argument("Helmholtz refinement levels must satisfy 0 <= H <= h");

    RefineOutput coarse_output = refine_mesh_nvb(initial_mesh, H);
    TriMesh coarse = std::move(coarse_output.mesh);
    UniformFineHierarchy hierarchy =
        build_uniform_fine_hierarchy(coarse, h - H);
    const int coarse_elements = hierarchy.accumulated.P_elem.cols();
    return finish_problem_data(
        std::move(coarse),
        std::move(hierarchy.accumulated),
        std::vector<int>(coarse_elements, H),
        h,
        ell, std::move(hierarchy.meshes),
        std::move(hierarchy.node_prolongations),
        std::move(hierarchy.element_prolongations));
}

HelmholtzProblemData build_adaptive_helmholtz_problem_data(
    const TriMesh &coarse_mesh,
    const std::vector<int> &coarse_element_levels,
    int fine_level,
    int ell) {
    validate_pure_robin_mesh(coarse_mesh);
    if (coarse_element_levels.size() != coarse_mesh.elems.size())
        throw std::invalid_argument("adaptive coarse level count must match coarse elements");
    if (coarse_element_levels.empty()
        || *std::min_element(coarse_element_levels.begin(), coarse_element_levels.end()) < 0
        || *std::max_element(coarse_element_levels.begin(), coarse_element_levels.end()) > fine_level) {
        throw std::invalid_argument("adaptive coarse levels must lie between zero and the fine level");
    }
    adaptive::NestedFineMesh nested = adaptive::complete_to_fine_level(
        coarse_mesh, coarse_element_levels, fine_level);
    return finish_problem_data(
        coarse_mesh,
        std::move(nested.refinement),
        coarse_element_levels,
        fine_level,
        ell, {}, {}, {});
}

HelmholtzLodModel HelmholtzLodModel::build(const HelmholtzProblemConfig &config) {
    if (config.wavenumber <= 0.0)
        throw std::invalid_argument("Helmholtz wavenumber must be positive");
    HelmholtzProblemConfig resolved = config;
    if (resolved.initial_mesh.nodes.empty())
        resolved.initial_mesh = make_helmholtz_unit_square_mesh();

    const auto mesh_start = std::chrono::steady_clock::now();
    HelmholtzProblemData problem = build_helmholtz_problem_data(
        resolved.initial_mesh, resolved.H, resolved.h, resolved.ell);
    return build_with_problem(
        std::move(resolved), std::move(problem), elapsed_ms(mesh_start));
}

HelmholtzLodModel HelmholtzLodModel::build_adaptive(
    const HelmholtzProblemConfig &config,
    const TriMesh &coarse_mesh,
    const std::vector<int> &coarse_element_levels) {
    if (config.wavenumber <= 0.0)
        throw std::invalid_argument("Helmholtz wavenumber must be positive");
    HelmholtzProblemConfig resolved = config;
    if (!coarse_element_levels.empty())
        resolved.H = *std::min_element(coarse_element_levels.begin(), coarse_element_levels.end());

    const auto mesh_start = std::chrono::steady_clock::now();
    HelmholtzProblemData problem = build_adaptive_helmholtz_problem_data(
        coarse_mesh, coarse_element_levels, resolved.h, resolved.ell);
    return build_with_problem(
        std::move(resolved), std::move(problem), elapsed_ms(mesh_start));
}

HelmholtzLodModel HelmholtzLodModel::build_with_problem(
    HelmholtzProblemConfig config,
    HelmholtzProblemData problem,
    double mesh_and_interpolation_ms) {
    const auto total_start = std::chrono::steady_clock::now();
    HelmholtzLodModel model;
    model.config_ = std::move(config);
    model.problem_ = std::move(problem);
    model.build_timings_.mesh_and_interpolation_ms = mesh_and_interpolation_ms;

    auto stage_start = std::chrono::steady_clock::now();
    model.operators_ = assemble_helmholtz_operators(
        model.problem_.fine,
        model.config_.wavenumber,
        model.config_.diffusion,
        model.config_.refractive_index,
        model.config_.boundary_beta);
    model.build_timings_.operators_ms = elapsed_ms(stage_start);

    stage_start = std::chrono::steady_clock::now();
    model.correctors_ = build_helmholtz_correctors(
        model.problem_.coarse,
        model.problem_.fine,
        model.problem_.fine_element_prolongation,
        model.problem_.fine_dg_prolongation,
        model.problem_.quasi_interpolation,
        model.problem_.patches,
        model.problem_.fine_hierarchy_meshes,
        model.problem_.fine_node_level_prolongations,
        model.problem_.fine_element_level_prolongations,
        model.operators_,
        model.config_.patch_solver);
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
    if (model.config_.mode == HelmholtzPetrovMode::TwoSided)
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
    model.build_timings_.total_ms = mesh_and_interpolation_ms + elapsed_ms(total_start);
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

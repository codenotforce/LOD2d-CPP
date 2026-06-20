#include "lod/lod_model.h"
#include "fem/assemble_dg.h"
#include "lod/patches.h"
#include "lod/quasi_interp.h"
#include "mesh/refine.h"
#include <Eigen/Sparse>
#include <stdexcept>
#include <memory>
#include <utility>

namespace lod2d {
namespace {

std::vector<int> count_node_incidence(const TriMesh &mesh) {
    std::vector<int> counts(mesh.nodes.size(), 0);
    for (const auto &elem : mesh.elems) {
        for (int v : elem) ++counts[v];
    }
    for (int v : mesh.dirichlet) {
        if (v >= 0 && v < static_cast<int>(counts.size())) counts[v] = 0;
    }
    return counts;
}

Eigen::SparseMatrix<double> build_cg_to_dg(const TriMesh &mesh) {
    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(3 * mesh.elems.size());
    for (int e = 0; e < static_cast<int>(mesh.elems.size()); ++e) {
        for (int i = 0; i < 3; ++i)
            trips.emplace_back(3 * e + i, mesh.elems[e][i], 1.0);
    }
    Eigen::SparseMatrix<double> cg2dg(3 * static_cast<int>(mesh.elems.size()),
                                      static_cast<int>(mesh.nodes.size()));
    cg2dg.setFromTriplets(trips.begin(), trips.end());
    return cg2dg;
}

} // namespace

TriMesh make_unit_square_mesh() {
    TriMesh mesh;
    mesh.nodes = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    mesh.elems = {{0, 1, 3}, {1, 2, 3}};
    mesh.dirichlet = {0, 1, 2, 3};
    return mesh;
}

LodProblemData build_lod_problem_data(const TriMesh &initial_mesh, int H, int h) {
    if (H < 0 || h < H)
        throw std::invalid_argument("LOD refinement levels must satisfy 0 <= H <= h");
    if (initial_mesh.nodes.empty() || initial_mesh.elems.empty())
        throw std::invalid_argument("initial mesh must contain nodes and elements");

    auto coarse_out = refine_mesh(initial_mesh, H);
    auto fine_out = refine_mesh(coarse_out.mesh, h - H);

    LodProblemData problem;
    problem.coarse = std::move(coarse_out.mesh);
    problem.fine = std::move(fine_out.mesh);
    problem.NH = static_cast<int>(problem.coarse.nodes.size());
    problem.NTH = static_cast<int>(problem.coarse.elems.size());
    problem.Nh = static_cast<int>(problem.fine.nodes.size());
    problem.NTh = static_cast<int>(problem.fine.elems.size());
    problem.nngH = count_node_incidence(problem.coarse);
    problem.nngh = count_node_incidence(problem.fine);

    problem.dgHidx.resize(problem.NTH);
    for (int e = 0; e < problem.NTH; ++e) {
        for (int i = 0; i < 3; ++i) problem.dgHidx[e][i] = 3 * e + i;
    }

    problem.P_node = std::move(fine_out.P_node);
    problem.P_elem = std::move(fine_out.P_elem);
    problem.P_dg = std::move(fine_out.P_dg);
    return problem;
}

LodOperators build_lod_operators(
    const LodProblemData &problem,
    const std::vector<double> &Ah,
    int ell) {
    if (static_cast<int>(Ah.size()) != problem.NTh)
        throw std::invalid_argument("coefficient vector size must match fine element count");

    LodOperators operators;
    operators.element_stiffness = assemble_element_stiffness(problem.fine, Ah);
    operators.patch = build_patches(problem.coarse, ell);
    operators.fine_element_children = build_fine_element_children(problem.P_elem, problem.NTH);

    Eigen::SparseMatrix<double> cg2dgh = build_cg_to_dg(problem.fine);
    Eigen::SparseMatrix<double> IH = build_quasi_interp(
        problem.coarse, problem.fine, problem.P_dg, cg2dgh, problem.Nh, problem.NH);
    operators.interpolation_rows = build_interpolation_rows(IH, problem.NH);

    auto areas = compute_area(problem.fine);
    operators.Sh = assemble_cg_from_element_stiffness(problem.fine, operators.element_stiffness);
    operators.Mh = assemble_cg_mass(problem.fine, areas);
    return operators;
}

std::vector<CorrectorEntries> build_lod_correctors(
    const LodProblemData &problem,
    const LodOperators &operators,
    int d,
    CorrectorSolver solver) {
    Eigen::SparseMatrix<double> unused_sparse;
    return compute_all_correctors(
        operators.patch, problem.coarse, problem.NH, problem.nngH,
        unused_sparse, problem.fine, problem.Nh, problem.nngh,
        problem.dghidx, unused_sparse, unused_sparse, problem.P_dg,
        problem.dgHidx, unused_sparse, d, solver,
        &operators.element_stiffness,
        &operators.fine_element_children,
        &operators.interpolation_rows);
}

Eigen::SparseMatrix<double> build_lod_basis(
    const LodProblemData &problem,
    const std::vector<CorrectorEntries> &correctors) {
    return build_multiscale_basis(problem.P_node, problem.coarse, problem.Nh, correctors);
}

LodModel LodModel::build(const LodProblemConfig &config, const std::vector<double> &Ah) {
    LodProblemConfig resolved = config;
    if (resolved.initial_mesh.nodes.empty()) resolved.initial_mesh = make_unit_square_mesh();

    LodModel model;
    model.config_ = resolved;
    model.problem_ = build_lod_problem_data(resolved.initial_mesh, resolved.H, resolved.h);
    LodOperators operators = build_lod_operators(model.problem_, Ah, resolved.ell);
    std::vector<CorrectorEntries> correctors = build_lod_correctors(
        model.problem_, operators, resolved.d, resolved.solver);
    Eigen::SparseMatrix<double> G = build_lod_basis(model.problem_, correctors);
    if (!resolved.keep_setup_matrices) {
        Eigen::SparseMatrix<double>().swap(model.problem_.P_elem);
        Eigen::SparseMatrix<double>().swap(model.problem_.P_dg);
    }
    model.system_ = std::make_unique<LodReusableSystem>(
        std::move(G), std::move(operators.Sh), std::move(operators.Mh),
        model.problem_.P_node, model.problem_.NH, model.problem_.coarse.dirichlet);
    return model;
}

const LodReusableSystem &LodModel::reusable_system() const {
    if (!system_) throw std::logic_error("LOD model has not been built");
    return *system_;
}

LodReuseSolution LodModel::solve_from_coarse_values(const Eigen::VectorXd &f_coarse) const {
    return reusable_system().solve_from_coarse_values(f_coarse);
}

LodReuseSolution LodModel::solve_from_fine_values(const Eigen::VectorXd &f_fine) const {
    return reusable_system().solve_from_fine_values(f_fine);
}

} // namespace lod2d

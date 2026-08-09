#include "helmholtz/adaptive/driver.h"

#include <Eigen/Cholesky>
#include <Eigen/SVD>

#include <array>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace lod2d::helmholtz::adaptive {
namespace {

double elapsed_ms(const std::chrono::steady_clock::time_point &start) {
    return std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - start).count();
}

double element_diameter(const TriMesh &mesh, const Triangle &tri) {
    return std::max({
        (mesh.nodes[tri[0]] - mesh.nodes[tri[1]]).norm(),
        (mesh.nodes[tri[1]] - mesh.nodes[tri[2]]).norm(),
        (mesh.nodes[tri[2]] - mesh.nodes[tri[0]]).norm()});
}

double q_max_value(
    const HelmholtzProblemData &problem,
    const std::vector<int> &parents) {
    std::vector<double> coarse_h(problem.coarse.elems.size());
    for (int t = 0; t < static_cast<int>(problem.coarse.elems.size()); ++t)
        coarse_h[t] = element_diameter(problem.coarse, problem.coarse.elems[t]);
    double result = 0.0;
    for (int t = 0; t < static_cast<int>(problem.fine.elems.size()); ++t) {
        const double h = element_diameter(problem.fine, problem.fine.elems[t]);
        result = std::max(result, h / coarse_h[parents[t]]);
    }
    return result;
}

double q_effective_value(
    const HelmholtzProblemData &problem,
    const std::vector<int> &parents,
    const std::vector<double> &indicator_squared) {
    std::vector<double> coarse_h(problem.coarse.elems.size());
    std::vector<double> local_q(problem.coarse.elems.size(), 0.0);
    for (int t = 0; t < static_cast<int>(problem.coarse.elems.size()); ++t)
        coarse_h[t] = element_diameter(problem.coarse, problem.coarse.elems[t]);
    for (int t = 0; t < static_cast<int>(problem.fine.elems.size()); ++t) {
        const int parent = parents[t];
        const double h = element_diameter(problem.fine, problem.fine.elems[t]);
        local_q[parent] = std::max(local_q[parent], h / coarse_h[parent]);
    }
    const double denominator = std::accumulate(
        indicator_squared.begin(), indicator_squared.end(), 0.0);
    if (denominator <= 0.0) return 0.0;
    double numerator = 0.0;
    for (int t = 0; t < static_cast<int>(indicator_squared.size()); ++t)
        numerator += indicator_squared[t] * local_q[t] * local_q[t];
    return std::sqrt(numerator / denominator);
}

double relative_norm(double error, double reference) {
    return error / std::max(reference, 1e-30);
}

double energy_inf_sup(const HelmholtzLodModel &model, int max_dofs) {
    const int dofs = model.coarse_operator().rows();
    if (dofs == 0 || dofs > max_dofs)
        return std::numeric_limits<double>::quiet_NaN();

    Eigen::SparseMatrix<double> energy = model.operators().stiffness;
    energy += model.config().wavenumber * model.config().wavenumber
            * model.operators().mass;
    const ComplexSparseMatrix complex_energy = energy.cast<Complex>();
    const ComplexMatrix trial_gram(
        model.trial_basis().adjoint() * complex_energy * model.trial_basis());
    const ComplexMatrix test_gram(
        model.test_basis().adjoint() * complex_energy * model.test_basis());
    Eigen::LLT<ComplexMatrix> trial_llt(trial_gram);
    Eigen::LLT<ComplexMatrix> test_llt(test_gram);
    if (trial_llt.info() != Eigen::Success || test_llt.info() != Eigen::Success)
        return std::numeric_limits<double>::quiet_NaN();

    const ComplexMatrix coarse(model.coarse_operator());
    const ComplexMatrix left_scaled = test_llt.matrixL().solve(coarse);
    const ComplexMatrix scaled = trial_llt.matrixL()
        .solve(left_scaled.adjoint()).adjoint();
    Eigen::JacobiSVD<ComplexMatrix> svd(
        scaled, Eigen::ComputeThinU | Eigen::ComputeThinV);
    return svd.singularValues().minCoeff();
}

using TriangleGeometrySignature = std::array<long long, 6>;

std::vector<TriangleGeometrySignature> canonical_triangle_geometry(const TriMesh &mesh) {
    std::vector<TriangleGeometrySignature> result;
    result.reserve(mesh.elems.size());
    for (const Triangle &triangle : mesh.elems) {
        std::array<std::pair<long long, long long>, 3> vertices;
        for (int i = 0; i < 3; ++i) {
            const Point2 &point = mesh.nodes[triangle[i]];
            vertices[i] = {
                std::llround(point.x() * 1e12),
                std::llround(point.y() * 1e12)};
        }
        std::sort(vertices.begin(), vertices.end());
        TriangleGeometrySignature signature{};
        for (int i = 0; i < 3; ++i) {
            signature[2 * i] = vertices[i].first;
            signature[2 * i + 1] = vertices[i].second;
        }
        result.push_back(signature);
    }
    std::sort(result.begin(), result.end());
    return result;
}

} // namespace

const char *residual_estimator_name(diagnostics::ResidualEstimatorKind kind) {
    if (kind == diagnostics::ResidualEstimatorKind::Fine) return "fine";
    if (kind == diagnostics::ResidualEstimatorKind::Mixed) return "mixed";
    return "macro";
}

AdaptiveHelmholtzResult run_adaptive_helmholtz(
    const AdaptiveHelmholtzConfig &config,
    const ComplexFunction &source,
    const ComplexFunction &exact,
    const ComplexGradientFunction &exact_gradient) {
    if (!source) throw std::invalid_argument("adaptive Helmholtz source is empty");
    if (config.problem.wavenumber <= 0.0)
        throw std::invalid_argument("adaptive Helmholtz wavenumber must be positive");
    if (config.problem.ell < 0 || config.max_iterations <= 0 || config.max_coarse_dofs <= 0)
        throw std::invalid_argument("adaptive Helmholtz integer parameters are invalid");
    if (!(config.theta > 0.0 && config.theta <= 1.0)
        || !(config.q_limit > 0.0 && config.q_limit <= 1.0)) {
        throw std::invalid_argument("adaptive theta or q limit is invalid");
    }
    if (static_cast<bool>(exact) != static_cast<bool>(exact_gradient))
        throw std::invalid_argument("exact value and gradient callbacks must be provided together");

    TriMesh initial = config.problem.initial_mesh;
    if (initial.nodes.empty()) initial = make_helmholtz_unit_square_mesh();
    AdaptiveMeshHierarchy hierarchy(
        initial, config.initial_coarse_level, config.fine_level);
    AdaptiveHelmholtzResult result;
    std::vector<TriangleGeometrySignature> fixed_fine_geometry;
    if (exact) {
        fixed_fine_geometry = canonical_triangle_geometry(
            refine_mesh_nvb(initial, config.fine_level).mesh);
    }

    for (int iteration = 0; iteration < config.max_iterations; ++iteration) {
        const auto iteration_start = std::chrono::steady_clock::now();
        AdaptiveIterationRecord record;
        record.iteration = iteration;

        HelmholtzProblemConfig problem_config = config.problem;
        problem_config.H = config.initial_coarse_level;
        problem_config.h = config.fine_level;
        auto stage_start = std::chrono::steady_clock::now();
        HelmholtzLodModel model = HelmholtzLodModel::build_adaptive(
            problem_config, hierarchy.coarse_mesh(), hierarchy.coarse_levels());
        if (exact
            && canonical_triangle_geometry(model.problem().fine) != fixed_fine_geometry) {
            throw std::runtime_error("adaptive completion changed the fixed fine mesh");
        }
        record.build_ms = elapsed_ms(stage_start);
        record.coarse_nodes = static_cast<int>(model.problem().coarse.nodes.size());
        record.coarse_elements = static_cast<int>(model.problem().coarse.elems.size());
        record.fine_nodes = static_cast<int>(model.problem().fine.nodes.size());
        record.fine_elements = static_cast<int>(model.problem().fine.elems.size());
        record.min_coarse_level = *std::min_element(
            hierarchy.coarse_levels().begin(), hierarchy.coarse_levels().end());
        record.max_coarse_level = *std::max_element(
            hierarchy.coarse_levels().begin(), hierarchy.coarse_levels().end());
        record.H_max = max_element_diameter(model.problem().coarse);

        const ComplexVector load = assemble_helmholtz_load(
            model.problem().fine, source,
            model.config().quadrature, model.config().quadrature_context);
        stage_start = std::chrono::steady_clock::now();
        const HelmholtzLodSolution lod = model.solve_load(load);
        record.solve_ms = elapsed_ms(stage_start);
        stage_start = std::chrono::steady_clock::now();
        const ComplexVector reference = model.solve_fine_reference(load);
        record.reference_ms = elapsed_ms(stage_start);
        const ComplexVector error = reference - lod.fine_values;
        record.energy_error = discrete_energy_norm(model.operators(), error);
        record.relative_energy_error = relative_norm(
            record.energy_error, discrete_energy_norm(model.operators(), reference));
        record.l2_error = discrete_l2_norm(model.operators(), error);
        record.relative_l2_error = relative_norm(
            record.l2_error, discrete_l2_norm(model.operators(), reference));
        if (exact) {
            const HelmholtzError lod_exact = compute_helmholtz_error(
                model.problem().fine,
                lod.fine_values,
                model.config().wavenumber,
                exact,
                exact_gradient,
                model.config().quadrature,
                model.config().quadrature_context);
            const HelmholtzError fine_exact = compute_helmholtz_error(
                model.problem().fine,
                reference,
                model.config().wavenumber,
                exact,
                exact_gradient,
                model.config().quadrature,
                model.config().quadrature_context);
            record.exact_energy_error = lod_exact.energy;
            record.exact_l2_error = lod_exact.l2;
            record.fine_exact_energy_error = fine_exact.energy;
            record.fine_exact_l2_error = fine_exact.l2;
        }

        stage_start = std::chrono::steady_clock::now();
        const diagnostics::HelmholtzResidualContributions contributions =
            diagnostics::assemble_helmholtz_residual_contributions(
                model.problem(), model.operators(), lod.fine_values, load, source,
                model.config().quadrature, model.config().quadrature_context);
        const diagnostics::HelmholtzIndicatorSet indicators = diagnostics::build_helmholtz_indicators(
            model.problem(), contributions);
        record.estimate_ms = elapsed_ms(stage_start);
        record.eta_fine = indicators.fine;
        record.eta_mixed = indicators.mixed;
        record.eta_macro = indicators.macro;
        record.selected_estimator = indicators.global(config.estimator);
        record.selected_effectivity = record.selected_estimator
            / std::max(record.energy_error, 1e-30);
        if (exact) {
            record.exact_effectivity = record.selected_estimator
                / std::max(record.exact_energy_error, 1e-30);
        }
        record.residual_identity_error = contributions.algebraic_relative_difference;
        record.q_max = q_max_value(model.problem(), contributions.fine_element_parent);
        record.q_effective = q_effective_value(
            model.problem(),
            contributions.fine_element_parent,
            indicators.squared(config.estimator));

        const std::vector<double> &selected_squared = indicators.squared(config.estimator);
        result.final_indicator_squared = selected_squared;
        std::vector<int> marked;
        if (config.compute_dual_calibration) {
            stage_start = std::chrono::steady_clock::now();
            const std::vector<double> dual = diagnostics::build_local_dual_indicators(
                model.problem(), model.operators(), contributions, config.dual_patch_layers);
            record.dual_ms = elapsed_ms(stage_start);
            std::vector<double> dual_squared(dual.size());
            std::transform(dual.begin(), dual.end(), dual_squared.begin(),
                           [](double value) { return value * value; });
            record.dual_spearman = spearman_rank_correlation(selected_squared, dual_squared);
            const std::vector<int> selected_marked = mark_doerfler(selected_squared, config.theta);
            const std::vector<int> dual_marked = mark_doerfler(dual_squared, config.theta);
            record.dual_marking_overlap = doerfler_energy_overlap(
                selected_marked, dual_marked, dual_squared);
        }

        record.petrov_residual = lod.petrov_residual;
        record.corrector_residual = model.correctors().diagnostics.max_primal_residual;
        record.constraint_residual = model.correctors().diagnostics.max_constraint_residual;
        record.inf_sup = energy_inf_sup(model, config.stability_max_dofs);

        bool stop = false;
        if (record.selected_estimator <= config.tolerance && config.tolerance > 0.0) {
            result.stop_reason = "estimator tolerance reached";
            stop = true;
        } else if (record.coarse_nodes >= config.max_coarse_dofs) {
            result.stop_reason = "coarse dof limit reached";
            stop = true;
        } else if (record.q_max > config.q_limit) {
            result.stop_reason = "coarse/fine scale separation limit reached";
            stop = true;
        } else if (iteration + 1 >= config.max_iterations) {
            result.stop_reason = "maximum adaptive iterations reached";
            stop = true;
        }

        if (!stop) {
            stage_start = std::chrono::steady_clock::now();
            const int maximum_allowed_level = std::min(
                config.fine_level - 1,
                static_cast<int>(std::floor(
                    config.fine_level + 2.0 * std::log2(config.q_limit) + 1e-12)));
            std::vector<char> eligible(hierarchy.coarse_levels().size(), false);
            double eligible_energy = 0.0;
            const double total_energy = std::accumulate(
                selected_squared.begin(), selected_squared.end(), 0.0);
            for (int t = 0; t < static_cast<int>(eligible.size()); ++t)
                if ((eligible[t] =
                        hierarchy.coarse_levels()[t] < maximum_allowed_level)) {
                    eligible_energy += selected_squared[t];
                }
            if (!(total_energy > 0.0)) {
                result.stop_reason = "selected estimator is zero";
                stop = true;
            } else if (eligible_energy + 1e-14 * total_energy
                       < config.theta * total_energy) {
                result.stop_reason =
                    "eligible indicator energy cannot satisfy Doerfler marking";
                stop = true;
            } else {
                marked = mark_doerfler(selected_squared, config.theta, eligible);
            }
            if (!stop && marked.empty()) {
                result.stop_reason = "no eligible elements within scale separation limit";
                stop = true;
            }
            if (!stop) {
                AdaptiveMeshHierarchy candidate = hierarchy;
                candidate.refine(marked);
                const int candidate_max_level = *std::max_element(
                    candidate.coarse_levels().begin(), candidate.coarse_levels().end());
                if (candidate_max_level > maximum_allowed_level) {
                    result.stop_reason =
                        "NVB closure would exceed scale separation limit";
                    stop = true;
                } else {
                    record.marked_elements = static_cast<int>(marked.size());
                    const int old_elements =
                        static_cast<int>(hierarchy.coarse_mesh().elems.size());
                    hierarchy = std::move(candidate);
                    const int new_elements =
                        static_cast<int>(hierarchy.coarse_mesh().elems.size());
                    record.closure_added_elements = std::max(
                        0, new_elements - old_elements - record.marked_elements);
                }
            }
            record.mark_refine_ms = elapsed_ms(stage_start);
        }

        record.total_ms = elapsed_ms(iteration_start);
        result.history.push_back(record);
        if (stop) break;
    }

    result.final_coarse_mesh = hierarchy.coarse_mesh();
    result.final_coarse_levels = hierarchy.coarse_levels();
    result.final_coarse_element_ids = hierarchy.coarse_element_ids();
    if (result.stop_reason.empty()) result.stop_reason = "adaptive loop completed";
    return result;
}

} // namespace lod2d::helmholtz::adaptive

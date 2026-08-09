#include "helmholtz/adaptive/error_control.h"
#include "helmholtz/adaptive/estimator.h"
#include "helmholtz/adaptive/kernel_residual.h"
#include "helmholtz/benchmarks/paper_cases.h"
#include "helmholtz/model.h"

#include <Eigen/QR>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <vector>

using namespace lod2d;
using namespace lod2d::helmholtz;
using namespace lod2d::helmholtz::adaptive;
using namespace lod2d::helmholtz::benchmarks;
using namespace lod2d::helmholtz::experiments;

namespace {

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

void require_close(double first, double second, double tolerance, const char *message) {
    const double scale = std::max({1.0, std::abs(first), std::abs(second)});
    if (std::abs(first - second) > tolerance * scale)
        throw std::runtime_error(message);
}

struct AuditProblem {
    PaperCaseData data;
    AdaptiveMeshHierarchy hierarchy;
    HelmholtzOperators operators;
    ComplexVector load;

    AuditProblem(PaperCase id, int coarse_level, int audit_level)
        : data(make_paper_case(id, 4.0)),
          hierarchy(data.initial_mesh, coarse_level, audit_level),
          operators(assemble_helmholtz_operators(
              hierarchy.cert_audit_mesh(), data.wavenumber)),
          load(assemble_helmholtz_load(
              hierarchy.cert_audit_mesh(), data.source, {}, data.quadrature_context)) {}
};

struct LodAuditCandidate {
    ComplexVector values;
    double petrov_residual = 0.0;
};

LodAuditCandidate solve_lod_candidate(
    const AuditProblem &problem,
    int ell = 1) {
    HelmholtzProblemConfig config;
    config.H = *std::min_element(
        problem.hierarchy.coarse_levels().begin(),
        problem.hierarchy.coarse_levels().end());
    config.h = problem.hierarchy.fine_level();
    config.ell = ell;
    config.wavenumber = problem.data.wavenumber;
    config.initial_mesh = problem.data.initial_mesh;
    config.quadrature_context = problem.data.quadrature_context;
    HelmholtzLodModel model = HelmholtzLodModel::build_adaptive(
        config,
        problem.hierarchy.coarse_mesh(),
        problem.hierarchy.coarse_levels(),
        problem.hierarchy.fine_mesh(),
        problem.hierarchy.fine_element_levels());
    const HelmholtzLodSolution solution = model.solve_source(problem.data.source);
    LodAuditCandidate result;
    result.values = problem.hierarchy.fine_to_cert_audit().cast<Complex>()
        * solution.fine_values;
    result.petrov_residual = solution.petrov_residual;
    return result;
}

void verify_patch_constraints_and_solver_agreement() {
    static_assert(std::is_copy_constructible_v<AuditKernelResidualEvidence>);
    static_assert(!std::is_aggregate_v<AuditKernelResidualEvidence>);
    const AuditKernelResidualEvidence uninitialized_evidence;
    require(!uninitialized_evidence.verified(),
            "default eta_H evidence was promoted to verified");
    require(uninitialized_evidence.diagnostic_fingerprint().empty(),
            "default eta_H evidence has a fabricated fingerprint");
    AuditKernelResidualEstimate uninitialized_estimate;
    require(!uninitialized_evidence.matches_result(uninitialized_estimate),
            "default eta_H evidence matched an uninitialized result");

    AuditProblem problem(PaperCase::R1, 1, 3);
    const ComplexVector candidate = ComplexVector::Zero(problem.load.size());
    const AuditKernelResidualEstimate saddle = estimate_audit_kernel_residual(
        problem.hierarchy,
        problem.operators,
        problem.load,
        candidate,
        0.5,
        KernelRieszSolver::SaddlePoint);
    const AuditKernelResidualEstimate basis = estimate_audit_kernel_residual(
        problem.hierarchy,
        problem.operators,
        problem.load,
        candidate,
        0.5,
        KernelRieszSolver::KernelBasisReference);

    require(saddle.policy.patch_layers
                == saddle.policy.interpolation_support_layers + 1,
            "D_z layer count is not p_I+1");
    require(saddle.policy.enlargement_layers == 1,
            "D_z+ is not the frozen one-layer enlargement");
    require(!saddle.policy.hash.empty(), "patch policy hash is empty");
    require(saddle.patches.size() == problem.hierarchy.coarse_mesh().nodes.size(),
            "there is not one audit-kernel patch per geometric coarse node");
    require(saddle.global_residual.isApprox(problem.load, 2e-13),
            "zero candidate did not retain the single global audit residual");
    require(saddle.eta > 0.0, "nonzero audit residual produced a zero eta_H");
    require(saddle.marked_elements.size() > 0,
            "positive eta_H produced an empty Doerfler set");
    require(!saddle.evidence.verified()
                && saddle.evidence.level()
                    == KernelResidualEvidenceLevel::Diagnostic,
            "ordinary Eigen eta_H was promoted to verified");
    require(saddle.evidence.backend().find("Eigen/") == 0,
            "eta_H evidence does not identify the Eigen backend");
    require(!saddle.evidence.source().empty()
                && !saddle.evidence.failure_reason().empty(),
            "unverified eta_H evidence has incomplete provenance");
    require(!saddle.evidence.context_fingerprint().empty()
                && !saddle.evidence.diagnostic_fingerprint().empty()
                && !saddle.evidence.result_fingerprint().empty(),
            "eta_H evidence fingerprints are empty");
    require(saddle.evidence.matches_result(saddle),
            "eta_H evidence does not match its estimator result");
    AuditKernelResidualEstimate modified = saddle;
    modified.eta = std::nextafter(modified.eta,
                                  std::numeric_limits<double>::infinity());
    require(!modified.evidence.matches_result(modified),
            "eta_H evidence accepted a modified eta value");
    modified = saddle;
    modified.element_eta_squared.front() = std::nextafter(
        modified.element_eta_squared.front(),
        std::numeric_limits<double>::infinity());
    require(!modified.evidence.matches_result(modified),
            "eta_H evidence accepted a modified element allocation");
    modified = saddle;
    modified.local_results.front().eta = std::nextafter(
        modified.local_results.front().eta,
        std::numeric_limits<double>::infinity());
    require(!modified.evidence.matches_result(modified),
            "eta_H evidence accepted a modified local Riesz result");
    modified = saddle;
    ++modified.patches.front().coarse_node;
    require(!modified.evidence.matches_result(modified),
            "eta_H evidence accepted modified patch metadata");
    require(saddle.evidence.context_fingerprint()
                == basis.evidence.context_fingerprint(),
            "eta_H context fingerprint depends on the diagnostic solver");
    require(saddle.evidence.diagnostic_fingerprint()
                != basis.evidence.diagnostic_fingerprint(),
            "eta_H diagnostic fingerprint ignores the Riesz solver");

    const AuditKernelResidualEstimate repeated = estimate_audit_kernel_residual(
        problem.hierarchy,
        problem.operators,
        problem.load,
        candidate,
        0.5,
        KernelRieszSolver::SaddlePoint);
    require(repeated.evidence.diagnostic_fingerprint()
                == saddle.evidence.diagnostic_fingerprint(),
            "identical eta_H inputs produced a non-deterministic fingerprint");
    ComplexVector perturbed_candidate = candidate;
    perturbed_candidate(0) = Complex(1e-6, -2e-6);
    const AuditKernelResidualEstimate perturbed = estimate_audit_kernel_residual(
        problem.hierarchy,
        problem.operators,
        problem.load,
        perturbed_candidate,
        0.5,
        KernelRieszSolver::SaddlePoint);
    require(perturbed.evidence.context_fingerprint()
                == saddle.evidence.context_fingerprint(),
            "eta_H context fingerprint depends on the candidate solution");
    require(perturbed.evidence.diagnostic_fingerprint()
                != saddle.evidence.diagnostic_fingerprint(),
            "eta_H diagnostic fingerprint ignores the candidate solution");

    for (int index = 0; index < static_cast<int>(saddle.patches.size()); ++index) {
        const AuditKernelPatch &patch = saddle.patches[index];
        require(std::includes(
                    patch.coarse_elements.begin(), patch.coarse_elements.end(),
                    patch.coarse_hat_support.begin(), patch.coarse_hat_support.end()),
                "D_z lost support(lambda_z)");
        require(std::includes(
                    patch.coarse_elements.begin(), patch.coarse_elements.end(),
                    patch.interpolation_support.begin(), patch.interpolation_support.end()),
                "D_z lost interpolation-propagated support");
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(
            patch.constraints.transpose());
        qr.setThreshold(2e-12);
        require(qr.rank() == patch.constraints.rows(),
                "dependent local constraint rows survived reduction");

        const LocalKernelRieszResult &saddle_local = saddle.local_results[index];
        const LocalKernelRieszResult &basis_local = basis.local_results[index];
        require(saddle_local.constraint_relative_residual <= 1e-10,
                "saddle solution violates I_H xi_z = 0");
        require(saddle_local.riesz_relative_residual <= 1e-10,
                "saddle local Riesz equation residual is too large");
        require(saddle_local.energy_identity_relative_error <= 1e-10,
                "saddle local Riesz energy identity failed");
        require(basis_local.constraint_relative_residual <= 1e-10,
                "kernel-basis solution violates I_H xi_z = 0");
        require(basis_local.riesz_relative_residual <= 1e-10,
                "kernel-basis local Riesz equation residual is too large");
        require(basis_local.energy_identity_relative_error <= 1e-10,
                "kernel-basis local Riesz energy identity failed");
        ComplexVector saddle_ambient = ComplexVector::Zero(problem.load.size());
        ComplexVector basis_ambient = ComplexVector::Zero(problem.load.size());
        for (int local = 0; local < static_cast<int>(patch.audit_dofs.size()); ++local) {
            saddle_ambient(patch.audit_dofs[local]) = saddle_local.local_values(local);
            basis_ambient(patch.audit_dofs[local]) = basis_local.local_values(local);
        }
        const ComplexVector saddle_full_constraint =
            problem.hierarchy.cert_audit_quasi_interpolation().cast<Complex>()
            * saddle_ambient;
        const ComplexVector basis_full_constraint =
            problem.hierarchy.cert_audit_quasi_interpolation().cast<Complex>()
            * basis_ambient;
        require(saddle_full_constraint.norm() <= 1e-10 * std::max(1.0, saddle_ambient.norm()),
                "saddle zero extension violates the full global I_H constraint");
        require(basis_full_constraint.norm() <= 1e-10 * std::max(1.0, basis_ambient.norm()),
                "kernel-basis zero extension violates the full global I_H constraint");
        require_close(
            saddle_local.eta,
            basis_local.eta,
            2e-10,
            "saddle and kernel-basis eta_H,z disagree");
    }

    const double node_sum = std::accumulate(
        saddle.node_eta_squared.begin(), saddle.node_eta_squared.end(), 0.0);
    const double element_sum = std::accumulate(
        saddle.element_eta_squared.begin(), saddle.element_eta_squared.end(), 0.0);
    require_close(node_sum, element_sum, 1e-12,
                  "node-to-element eta allocation is not conservative");
    require(saddle.allocation_relative_error <= 1e-12,
            "reported eta allocation residual is too large");
}

void verify_explicit_global_kernel_bounds() {
    AuditProblem problem(PaperCase::R1, 1, 2);
    const AuditKernelResidualEstimate estimate = estimate_audit_kernel_residual(
        problem.hierarchy,
        problem.operators,
        problem.load,
        ComplexVector::Zero(problem.load.size()),
        0.5);
    const ExplicitGlobalKernelComparison comparison =
        compare_with_explicit_global_kernel(
            problem.hierarchy, problem.operators, estimate, 512);
    require(comparison.global_kernel_dimension > 0,
            "explicit test mesh unexpectedly has a trivial global kernel");
    require(comparison.local_kernel_columns >= comparison.global_kernel_dimension,
            "local kernels cannot span the explicit global kernel");
    require(comparison.span_relative_residual <= 2e-10,
            "local audit kernels do not span the explicit global kernel");
    require(comparison.lower_direction_holds,
            "eta_H/C_ov <= ||r||_W' failed on the explicit kernel");
    require(comparison.upper_direction_holds,
            "||r||_W' <= C_sd eta_H failed on the explicit kernel");
}

void verify_effectivity_distribution(PaperCase id) {
    AuditProblem problem(id, 0, 2);
    const LodAuditCandidate lod = solve_lod_candidate(problem);
    const ComplexVector &candidate = lod.values;
    const ComplexVector certification = solve_helmholtz_fem(
        problem.operators, problem.load);
    require(lod.petrov_residual <= 1e-10,
            "paper-case LOD candidate has a large Petrov residual");
    require(candidate.norm() > 1e-12,
            "paper-case effectivity used a zero LOD candidate");
    require(discrete_energy_norm(
                problem.operators, certification - candidate) > 1e-12,
            "paper-case effectivity did not use an independent LOD candidate");
    const AuditKernelResidualEstimate estimate = estimate_audit_kernel_residual(
        problem.hierarchy,
        problem.operators,
        problem.load,
        candidate,
        0.5);
    const LocalKernelEffectivityReport report = evaluate_local_kernel_effectivity(
        problem.hierarchy,
        problem.operators,
        estimate,
        certification,
        candidate);
    require(report.local.size() == problem.hierarchy.coarse_mesh().nodes.size(),
            "local effectivity count does not match coarse nodes");
    require(report.distribution.count > 0,
            "paper case produced no valid local effectivity ratios");
    require(report.distribution.minimum <= report.distribution.first_quartile
                && report.distribution.first_quartile <= report.distribution.median
                && report.distribution.median <= report.distribution.third_quartile
                && report.distribution.third_quartile <= report.distribution.percentile_90
                && report.distribution.percentile_90 <= report.distribution.maximum,
            "local effectivity quantiles are not ordered");
}

void verify_equal_mesh_degeneracy() {
    AuditProblem problem(PaperCase::R1, 0, 0);
    const LodAuditCandidate lod = solve_lod_candidate(problem);
    const ComplexVector &candidate = lod.values;
    const ComplexVector certification = solve_helmholtz_fem(
        problem.operators, problem.load);
    const AuditKernelResidualEstimate estimate = estimate_audit_kernel_residual(
        problem.hierarchy,
        problem.operators,
        problem.load,
        candidate,
        0.5);
    const HelmholtzError error = compute_discrete_helmholtz_error(
        problem.hierarchy.cert_audit_mesh(),
        problem.operators,
        certification,
        candidate);
    const ExplicitGlobalKernelComparison comparison =
        compare_with_explicit_global_kernel(
            problem.hierarchy, problem.operators, estimate, 128);
    require(comparison.global_kernel_dimension == 0,
            "T_H=T_h=T_audit did not collapse ker(I_H)");
    require(estimate.eta <= 2e-11,
            "equal-space audit-kernel residual is not at machine precision");
    require(lod.petrov_residual <= 1e-12,
            "equal-space LOD solve has a large Petrov residual");
    require(error.energy <= 2e-12,
            "independent equal-space FEM and LOD solutions do not agree");
    // Strong/broken residuals are diagnostics and intentionally do not enter
    // this degeneracy gate.
}

} // namespace

int main() {
    try {
        verify_patch_constraints_and_solver_agreement();
        verify_explicit_global_kernel_bounds();
        verify_effectivity_distribution(PaperCase::R1);
        verify_effectivity_distribution(PaperCase::R2a);
        verify_effectivity_distribution(PaperCase::R2b);
        verify_effectivity_distribution(PaperCase::S);
        verify_equal_mesh_degeneracy();
        std::cout << "Audit-kernel residual estimator passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_kernel_residual failed: "
                  << error.what() << '\n';
        return 1;
    }
}

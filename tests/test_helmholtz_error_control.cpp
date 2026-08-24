#include "helmholtz/adaptive/error_control.h"

#include "helmholtz/benchmarks/paper_cases.h"
#include "mesh/refine.h"

#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
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

struct FemResult {
    TriMesh mesh;
    HelmholtzOperators operators;
    ComplexVector solution;
};

FemResult solve_case(const PaperCaseData &data, int level) {
    FemResult result;
    result.mesh = refine_mesh_nvb(data.initial_mesh, level).mesh;
    result.operators = assemble_helmholtz_operators(result.mesh, data.wavenumber);
    const ComplexVector load = assemble_helmholtz_load(
        result.mesh, data.source, {}, data.quadrature_context);
    result.solution = solve_helmholtz_fem(result.operators, load);
    return result;
}

std::pair<std::vector<int>, std::vector<int>> split_elements(const TriMesh &mesh) {
    std::vector<int> first;
    std::vector<int> second;
    for (int element = 0; element < static_cast<int>(mesh.elems.size()); ++element) {
        (element % 2 == 0 ? first : second).push_back(element);
    }
    return {first, second};
}

void verify_global_and_local_errors() {
    const PaperCaseData data = make_paper_case(PaperCase::R1, 4.0);
    const FemResult fem = solve_case(data, 2);
    const ComplexVector zero = ComplexVector::Zero(fem.solution.size());
    const HelmholtzError global = compute_discrete_helmholtz_error(
        fem.mesh, fem.operators, fem.solution, zero);
    const auto [first_elements, second_elements] = split_elements(fem.mesh);
    const HelmholtzError first = compute_local_discrete_helmholtz_error(
        fem.mesh, fem.operators, fem.solution, zero, first_elements);
    const HelmholtzError second = compute_local_discrete_helmholtz_error(
        fem.mesh, fem.operators, fem.solution, zero, second_elements);
    require_close(global.energy * global.energy,
                  first.energy * first.energy + second.energy * second.energy,
                  2e-13, "local discrete energy errors do not sum to the global error");
    require_close(global.l2 * global.l2,
                  first.l2 * first.l2 + second.l2 * second.l2,
                  2e-13, "local discrete L2 errors do not sum to the global error");

    const ErrorReference exact = ErrorReference::exact(
        data.exact, data.exact_gradient, {}, data.quadrature_context);
    require(exact.role() == ErrorReferenceRole::Exact,
            "exact error reference reports the wrong role");
    const HelmholtzError exact_global = exact.evaluate(
        fem.mesh, fem.operators, fem.solution);
    const HelmholtzError exact_first = compute_local_exact_helmholtz_error(
        fem.mesh, fem.operators, fem.solution, first_elements,
        data.exact, data.exact_gradient, {}, data.quadrature_context);
    const HelmholtzError exact_second = compute_local_exact_helmholtz_error(
        fem.mesh, fem.operators, fem.solution, second_elements,
        data.exact, data.exact_gradient, {}, data.quadrature_context);
    require_close(exact_global.energy * exact_global.energy,
                  exact_first.energy * exact_first.energy
                      + exact_second.energy * exact_second.energy,
                  2e-12, "local exact energy errors do not sum to the global error");
    require_close(exact_global.l2 * exact_global.l2,
                  exact_first.l2 * exact_first.l2 + exact_second.l2 * exact_second.l2,
                  2e-12, "local exact L2 errors do not sum to the global error");

    const ErrorReference two_level = ErrorReference::two_level(fem.solution);
    require(two_level.role() == ErrorReferenceRole::TwoLevelSaturation,
            "two-level error reference reports the wrong role");
    const HelmholtzError two_level_error = two_level.evaluate(
        fem.mesh, fem.operators, zero);
    require_close(two_level_error.energy, global.energy, 1e-14,
                  "two-level role does not use its private finer values");

    const ErrorReference external = ErrorReference::external(
        [](const TriMesh &, const ComplexVector &) { return HelmholtzError{3.0, 4.0}; });
    require(external.role() == ErrorReferenceRole::ExternalEstimator,
            "external error reference reports the wrong role");
    const HelmholtzError external_error = external.evaluate(
        fem.mesh, fem.operators, fem.solution);
    require(external_error.energy == 3.0 && external_error.l2 == 4.0,
            "external error estimator callback was not used");

    const ErrorReference invalid_external = ErrorReference::external(
        [](const TriMesh &, const ComplexVector &) {
            return HelmholtzError{-1.0, 0.0};
        });
    bool rejected = false;
    try {
        (void)invalid_external.evaluate(fem.mesh, fem.operators, fem.solution);
    } catch (const std::runtime_error &) {
        rejected = true;
    }
    require(rejected, "negative external error estimate was accepted");
}

void verify_separate_factorization_and_timing() {
    const PaperCaseData data = make_paper_case(PaperCase::R1, 4.0);
    const TriMesh mesh = refine_mesh_nvb(data.initial_mesh, 2).mesh;

    CertificationAuditService certification(mesh, data.wavenumber);
    const CertificationAuditSolution first = certification.solve_source(data.source);
    const CertificationAuditSolution second = certification.solve_source(data.source);
    require((first.values - second.values).norm() < 1e-14,
            "certification cache changed a repeated solve");
    require(certification.timings().factorization_builds == 1,
            "certification factorization was not cached");
    require(certification.timings().solve_count == 2,
            "certification solve count is incorrect");
    require(!certification.timings().excluded_from_method_time,
            "certification cost was excluded from CALOD method time");
    require(certification.timings().included_method_time_ms()
                == certification.timings().total_ms(),
            "certification included time is inconsistent");

    EvaluationReferenceService evaluation(mesh, data.wavenumber);
    require(!evaluation.ready(), "evaluation reference is ready before a solve");
    evaluation.prepare_source(data.source);
    require(evaluation.ready(), "evaluation reference did not become ready");
    const EvaluationReferenceError report = evaluation.evaluate_candidate_on_reference(
        ComplexVector::Zero(mesh.nodes.size()));
    require_close(report.relative_energy, 1.0, 1e-13,
                  "zero candidate has the wrong relative reference energy error");
    require_close(report.relative_l2, 1.0, 1e-13,
                  "zero candidate has the wrong relative reference L2 error");
    evaluation.prepare_source(data.source);
    require(evaluation.timings().factorization_builds == 1,
            "evaluation-reference factorization was not cached");
    require(evaluation.timings().solve_count == 2,
            "evaluation-reference solve count is incorrect");
    require(evaluation.timings().excluded_from_method_time,
            "evaluation-reference cost was included in method time");
    require(evaluation.timings().included_method_time_ms() == 0.0,
            "excluded evaluation-reference time leaked into method time");
    require(certification.timings().factorization_builds
                == evaluation.timings().factorization_builds,
            "test setup did not exercise two independent one-build caches");
}

void verify_empirical_saturation_audit_interval_and_marking() {
    const PaperCaseData data = make_paper_case(PaperCase::R1, 4.0);
    const FemResult audit = solve_case(data, 1);
    std::vector<int> audit_parent_fine(audit.mesh.elems.size());
    std::iota(audit_parent_fine.begin(), audit_parent_fine.end(), 0);
    const double saturation_factor = 0.05;
    const EmpiricalSaturationAuditEstimate estimate =
        estimate_empirical_saturation_audit_error(
            audit.mesh, audit.solution, data.source, data.wavenumber,
            audit_parent_fine, static_cast<int>(audit.mesh.elems.size()),
            saturation_factor, 0.5, {}, {}, 1.0, {},
            data.quadrature_context);
    require(!estimate.assumption_verified,
            "empirical saturation assumption was promoted to verified");
    require(estimate.two_level_energy_error > 0.0,
            "two-level audit probe produced a zero difference unexpectedly");
    require_close(
        estimate.audit_error_lower,
        estimate.two_level_energy_error / (1.0 + saturation_factor),
        1e-14, "saturation lower interval is incorrect");
    require_close(
        estimate.audit_error_upper,
        estimate.two_level_energy_error / (1.0 - saturation_factor),
        1e-14, "saturation upper interval is incorrect");
    require(estimate.audit_error_lower <= estimate.two_level_energy_error
                && estimate.two_level_energy_error
                    <= estimate.audit_error_upper,
            "two-level difference is outside its saturation interval");
    require(estimate.allocation_relative_error < 1e-12,
            "audit element contributions do not conserve the probe error");
    require(!estimate.marked_fine_elements.empty(),
            "positive audit probe error produced an empty fine marking");
}

void verify_exact_reference_convergence(PaperCase id) {
    const PaperCaseData data = make_paper_case(id, 4.0);
    const FemResult candidate = solve_case(data, 1);
    const ErrorReference exact = ErrorReference::exact(
        data.exact, data.exact_gradient, {}, data.quadrature_context);
    const HelmholtzError exact_error = exact.evaluate(
        candidate.mesh, candidate.operators, candidate.solution);

    std::array<double, 2> energy_differences{};
    std::array<double, 2> l2_differences{};
    int slot = 0;
    // The revised R1 load is strongly localized; two adjacent very coarse
    // levels can show pre-asymptotic L2 oscillation. Compare against a
    // genuinely finer evaluation reference while retaining the same check.
    for (int reference_level : {2, 4}) {
        const TriMesh reference_mesh =
            refine_mesh_nvb(data.initial_mesh, reference_level).mesh;
        EvaluationReferenceService evaluation(reference_mesh, data.wavenumber);
        evaluation.prepare_source(data.source, {}, data.quadrature_context);
        const EvaluationReferenceError reference_error =
            evaluation.evaluate_nested_candidate(candidate.mesh, candidate.solution);
        energy_differences[slot] = std::abs(
            reference_error.absolute.energy - exact_error.energy);
        l2_differences[slot++] = std::abs(
            reference_error.absolute.l2 - exact_error.l2);
    }
    std::cout << "case=" << static_cast<int>(id)
              << " exact/reference energy differences=" << energy_differences[0]
              << ',' << energy_differences[1]
              << " L2 differences=" << l2_differences[0]
              << ',' << l2_differences[1] << '\n';
    require(energy_differences[1] < energy_differences[0],
            "exact and evaluation-reference energy errors did not converge");
    require(l2_differences[1] < l2_differences[0],
            "exact and evaluation-reference L2 errors did not converge");
}

} // namespace

int main() {
    try {
        verify_global_and_local_errors();
        verify_separate_factorization_and_timing();
        verify_empirical_saturation_audit_interval_and_marking();
        verify_exact_reference_convergence(PaperCase::R1);
        verify_exact_reference_convergence(PaperCase::S);
        std::cout << "Helmholtz error roles, local errors, and fine-space caches passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_error_control failed: " << error.what() << '\n';
        return 1;
    }
}

#include "helmholtz/adaptive/practical_driver.h"
#include "helmholtz/benchmarks/paper_cases.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace lod2d;
using namespace lod2d::helmholtz;
using namespace lod2d::helmholtz::adaptive;
using namespace lod2d::helmholtz::benchmarks;
using namespace lod2d::helmholtz::experiments;

namespace {

void require(const bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

PracticalDriverProblem r1_problem() {
    const PaperCaseData data = make_paper_case(PaperCase::R1, 4.0);
    PracticalDriverProblem problem;
    problem.initial_mesh = data.initial_mesh;
    problem.source = data.source;
    problem.quadrature_context = data.quadrature_context;
    return problem;
}

PracticalDriverConfig base_config() {
    PracticalDriverConfig config;
    config.initial_coarse_level = 2;
    config.reference_level = 4;
    config.ell0 = 1;
    config.ell_max = 3;
    config.wavenumber = 4.0;
    config.c_H = 100.0;
    config.theta_loc = 100.0;
    config.theta_H = 0.5;
    config.rho_star = 0.25;
    config.tolerance_reference = 1e6;
    config.limits.maximum_iterations = 20;
    config.limits.maximum_H_steps = 4;
    config.limits.maximum_unknowns = 100000;
    config.limits.maximum_coarse_elements = 10000;
    config.limits.maximum_ambient_elements = 100000;
    return config;
}

std::size_t count_action(
    const PracticalDriverResult &result,
    const PracticalDriverAction action) {
    return static_cast<std::size_t>(std::count_if(
        result.journal.begin(),
        result.journal.end(),
        [action](const PracticalIterationRecord &record) {
            return record.action == action;
        }));
}

const PracticalIterationRecord &first_action(
    const PracticalDriverResult &result,
    const PracticalDriverAction action) {
    const auto found = std::find_if(
        result.journal.begin(),
        result.journal.end(),
        [action](const PracticalIterationRecord &record) {
            return record.action == action;
        });
    if (found == result.journal.end()) {
        throw std::runtime_error("expected practical-driver action was not recorded");
    }
    return *found;
}

void verify_real_reference_chain_converges() {
    const PracticalDriverConfig config = base_config();
    PracticalAdaptiveDriver driver(r1_problem(), config);
    const std::uint64_t reference_version =
        driver.hierarchy().reference_mesh_version();
    const std::size_t reference_nodes =
        driver.hierarchy().reference_mesh().nodes.size();
    const PracticalDriverResult result = driver.run();

    require(result.state == PracticalDriverState::Converged,
            "real practical reference chain did not converge at a loose tolerance");
    require(count_action(result, PracticalDriverAction::AcceptCoarse) == 1,
            "real chain did not accept coarse admissibility");
    require(count_action(result, PracticalDriverAction::AcceptLocalization) == 1,
            "real chain did not accept the localization certificate");
    require(count_action(result, PracticalDriverAction::Complete) == 1,
            "real chain did not record practical completion");
    require(std::isfinite(result.eta_H) && result.eta_H >= 0.0,
            "real chain returned an invalid reference residual estimator");
    require(std::isfinite(result.theta_loc) && result.theta_loc >= 0.0,
            "real chain returned an invalid localization certificate");
    const PracticalIterationRecord &complete =
        first_action(result, PracticalDriverAction::Complete);
    require(complete.evaluation_candidate.size() ==
                static_cast<int>(driver.hierarchy().reference_mesh().nodes.size()) &&
                complete.evaluation_candidate.allFinite(),
            "evaluation-only candidate was not exported after the decision");
    require(result.final_element_eta_squared.size() ==
                driver.hierarchy().coarse_mesh().elems.size(),
            "final eta_H,T field was not exported on the final coarse mesh");
    require(driver.hierarchy().reference_mesh_version() == reference_version &&
                driver.hierarchy().reference_mesh().nodes.size() == reference_nodes,
            "a practical iteration changed the frozen reference space");
}

void verify_g4_localization_only_increases_global_ell() {
    PracticalDriverConfig config = base_config();
    config.initial_coarse_level = 4;
    config.reference_level = 6;
    config.theta_loc = 1e-30;
    config.ell0 = 1;
    config.ell_max = 2;
    PracticalAdaptiveDriver driver(r1_problem(), config);
    const PracticalDriverResult result = driver.run();

    if (result.state != PracticalDriverState::WorkLimitReached) {
        throw std::runtime_error(
            std::string("ell_max localization failure stopped in ") +
            practical_driver_state_name(result.state) + ": " + result.stop_reason);
    }
    require(count_action(result, PracticalDriverAction::IncreaseGlobalEll) == 1,
            "localization failure did not increase the single global ell");
    require(result.ell == 2,
            "global ell was not advanced exactly to the configured ell_max");
    require(count_action(result, PracticalDriverAction::FormCoarseMarking) == 0 &&
                count_action(result, PracticalDriverAction::RefineCoarse) == 0 &&
                count_action(
                    result,
                    PracticalDriverAction::RefineCoarseForAdmissibility) == 0,
            "localization failure leaked into an H/fine refinement action");

    const PracticalIterationRecord &increase =
        first_action(result, PracticalDriverAction::IncreaseGlobalEll);
    require(increase.state_before == PracticalDriverState::LocalizationCheck &&
                increase.state_after == PracticalDriverState::LocalizationCheck,
            "G4 ell increase left the localization state");
    const auto next = std::next(
        result.journal.begin(), static_cast<std::ptrdiff_t>(increase.sequence + 1));
    require(next != result.journal.end() &&
                next->state_before == PracticalDriverState::LocalizationCheck,
            "the action after a failed localization observation was not another localization check");
}

void verify_fixed_ell_skips_localization_certificate() {
    PracticalDriverConfig config = base_config();
    config.localization_policy = PracticalLocalizationPolicy::FixedGlobalEll;
    config.ell0 = 2;
    config.ell_max = 2;
    config.theta_loc = 1e-30;
    PracticalAdaptiveDriver driver(r1_problem(), config);
    const PracticalDriverResult result = driver.run();

    require(result.state == PracticalDriverState::Converged,
            "fixed-ell HLOD chain did not complete");
    require(result.ell == 2,
            "fixed-ell HLOD changed the calibrated oversampling level");
    require(count_action(result, PracticalDriverAction::IncreaseGlobalEll) == 0,
            "fixed-ell HLOD entered the PALOD localization branch");
    const PracticalIterationRecord &accepted =
        first_action(result, PracticalDriverAction::AcceptFixedEll);
    require(count_action(result, PracticalDriverAction::AcceptLocalization) == 0,
            "fixed-ell HLOD falsely recorded localization acceptance");
    require(accepted.theta_loc == 0.0
                && accepted.time_certificate_seconds <= 1e-12,
            "fixed-ell HLOD computed a localization certificate");
    require(driver.hierarchy().ambient_mesh().nodes.size()
                == driver.hierarchy().reference_mesh().nodes.size()
                && std::all_of(
                    result.journal.begin(), result.journal.end(),
                    [](const PracticalIterationRecord &record) {
                        return record.ambient_refined_elements == 0;
                    }),
            "fixed-ell HLOD paid for an unused ambient shadow refinement");
}

void verify_estimator_driven_H_refinement_preserves_epoch() {
    PracticalDriverConfig config = base_config();
    config.reference_level = 5;
    config.tolerance_reference = 1e-14;
    config.limits.maximum_iterations = 5;
    config.limits.maximum_H_steps = 2;
    PracticalAdaptiveDriver driver(r1_problem(), config);
    const std::uint64_t epoch = driver.hierarchy().reference_epoch();
    const std::uint64_t reference_version =
        driver.hierarchy().reference_mesh_version();
    const std::size_t reference_nodes =
        driver.hierarchy().reference_mesh().nodes.size();
    const PracticalDriverResult result = driver.run();

    require(count_action(result, PracticalDriverAction::FormCoarseMarking) == 1,
            "reference residual did not form a coarse marking");
    require(count_action(result, PracticalDriverAction::RefineCoarse) == 1,
            "reference residual marking was not applied to H");
    require(result.H_steps == 1,
            "estimator-driven test performed an unexpected number of H steps");
    require(driver.hierarchy().reference_epoch() == epoch &&
                driver.hierarchy().reference_mesh_version() == reference_version &&
                driver.hierarchy().reference_mesh().nodes.size() == reference_nodes,
            "H refinement mutated the frozen reference epoch");
    require(driver.hierarchy().ambient_ratio() <= config.rho_star * (1.0 + 1e-12),
            "H refinement did not restore the ambient/coarse ratio gate");
}

void verify_reference_capacity_stops_transactionally() {
    PracticalDriverConfig config = base_config();
    config.initial_coarse_level = 1;
    config.reference_level = 2;
    config.tolerance_reference = 1e-14;
    config.theta_loc = 100.0;
    config.limits.maximum_H_steps = 10;
    config.limits.maximum_iterations = 30;
    PracticalAdaptiveDriver driver(r1_problem(), config);
    const std::uint64_t epoch = driver.hierarchy().reference_epoch();
    const std::uint64_t reference_version =
        driver.hierarchy().reference_mesh_version();
    const std::size_t reference_nodes =
        driver.hierarchy().reference_mesh().nodes.size();
    const PracticalDriverResult result = driver.run();

    require(result.state == PracticalDriverState::ReferenceRefreshRequired,
            "reference capacity exhaustion did not return ReferenceRefreshRequired");
    require(count_action(
                result,
                PracticalDriverAction::StopReferenceRefreshRequired) == 1,
            "reference capacity exhaustion lacked its structured stop action");
    require(driver.hierarchy().reference_epoch() == epoch &&
                driver.hierarchy().reference_mesh_version() == reference_version &&
                driver.hierarchy().reference_mesh().nodes.size() == reference_nodes,
            "capacity exhaustion silently refreshed the frozen reference epoch");
    require(driver.hierarchy().reference_embedding_holds(),
            "failed H refinement damaged the reference embedding");
    require(result.journal.back().coarse_elements ==
                result.journal[result.journal.size() - 2].coarse_elements,
            "failed H refinement partially committed a coarse mesh");
}

void verify_one_trajectory_multiple_target_extraction() {
    std::vector<PracticalIterationRecord> journal(4);
    journal[0].reference_energy_error = 0.4;
    journal[1].reference_energy_error = 0.08;
    journal[2].reference_energy_error = 0.018;
    journal[3].reference_energy_error = 0.009;
    const std::vector<PracticalTargetHit> hits =
        extract_practical_target_hits(
            journal, {0.1, 0.05, 0.02, 0.01, 0.001});
    require(hits.size() == 5, "multi-target extraction lost a target");
    require(hits[0].journal_index == 1 && hits[1].journal_index == 2 &&
                hits[2].journal_index == 2 && hits[3].journal_index == 3 &&
                !hits[4].journal_index,
            "multi-target extraction did not select first hits from one trajectory");
}

void verify_posterior_convergence_diagnostic() {
    std::vector<PracticalIterationRecord> stable(3);
    stable[0].coarse_nodes = 10;
    stable[0].reference_energy_error = 0.4;
    stable[1].coarse_nodes = 20;
    stable[1].reference_energy_error = 0.2;
    stable[2].coarse_nodes = 40;
    stable[2].reference_energy_error = 0.1;
    const PracticalConvergenceDiagnostic stable_result =
        diagnose_practical_convergence_regime(stable);
    require(stable_result.regime
                == PracticalConvergenceRegime::ObservedStableDecay
                && stable_result.distinct_points == 3
                && stable_result.previous_log_slope
                && stable_result.last_log_slope,
            "stable posterior error decay was not recognized");

    stable[1].reference_energy_error = 0.5;
    require(diagnose_practical_convergence_regime(stable).regime
                == PracticalConvergenceRegime::PreAsymptotic,
            "nonmonotone posterior error was not marked pre-asymptotic");

    stable.pop_back();
    require(diagnose_practical_convergence_regime(stable).regime
                == PracticalConvergenceRegime::InsufficientData,
            "two posterior error points were treated as a convergence regime");
}

} // namespace

int main() {
    try {
        verify_real_reference_chain_converges();
        verify_g4_localization_only_increases_global_ell();
        verify_fixed_ell_skips_localization_certificate();
        verify_estimator_driven_H_refinement_preserves_epoch();
        verify_reference_capacity_stops_transactionally();
        verify_one_trajectory_multiple_target_extraction();
        verify_posterior_convergence_diagnostic();
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_practical_driver failed: "
                  << error.what() << '\n';
        return 1;
    }
    std::cout << "test_helmholtz_practical_driver passed\n";
    return 0;
}

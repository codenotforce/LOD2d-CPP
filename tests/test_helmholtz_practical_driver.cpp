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
    PracticalAdaptiveDriver driver(r1_problem(), base_config());
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
    config.c_H = 1e-6;
    config.limits.maximum_H_steps = 10;
    config.limits.maximum_iterations = 20;
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

} // namespace

int main() {
    try {
        verify_real_reference_chain_converges();
        verify_g4_localization_only_increases_global_ell();
        verify_estimator_driven_H_refinement_preserves_epoch();
        verify_reference_capacity_stops_transactionally();
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_practical_driver failed: "
                  << error.what() << '\n';
        return 1;
    }
    std::cout << "test_helmholtz_practical_driver passed\n";
    return 0;
}

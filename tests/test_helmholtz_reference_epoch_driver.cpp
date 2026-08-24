#include "helmholtz/adaptive/reference_epoch_driver.h"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace lod2d::helmholtz::adaptive;

namespace {

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

struct ScriptedBackend final : ReferenceEpochDriverBackend {
    std::vector<double> corrector_bounds;
    std::vector<double> error_bounds;
    std::vector<double> dual_gaps;
    std::vector<bool> structural;
    std::size_t corrector_index = 0;
    std::size_t solve_index = 0;
    std::size_t dual_index = 0;
    std::size_t cycle_index = 0;
    std::vector<std::string> calls;
    ReferenceEpochResourceSnapshot snapshot{100, 100};
    int reference_level_gap = std::numeric_limits<int>::max();
    bool stop_inside_refresh = false;
    bool reserve_unavailable_inside_refresh = false;

    void begin_epoch() override { calls.push_back("begin"); }

    ReferenceEpochCorrectorObservation corrector_check(int ell) override {
        calls.push_back("corrector:" + std::to_string(ell));
        if (corrector_index >= corrector_bounds.size())
            throw std::runtime_error("corrector script exhausted");
        const double value = corrector_bounds[corrector_index++];
        ReferenceEpochCorrectorObservation observation;
        observation.theta_loc = value;
        observation.delta_loc_hat = value;
        observation.corrector_parallel_threads = 3;
        observation.corrector_patch_assembly_work_seconds = 1.25;
        observation.corrector_patch_solve_work_seconds = 2.5;
        observation.corrector_patch_pack_work_seconds = 0.125;
        observation.corrector_maximum_patch_dofs = 41;
        observation.corrector_maximum_patch_constraints = 7;
        observation.corrector_maximum_patch_rhs = 3;
        observation.gram_factor_cache_hits = 11;
        observation.gram_factor_cache_misses = 13;
        return observation;
    }

    ReferenceEpochSolveObservation solve_and_estimate() override {
        calls.push_back("solve");
        if (solve_index >= error_bounds.size())
            throw std::runtime_error("solve script exhausted");
        const double value = error_bounds[solve_index++];
        return {0.5 * value, value, {0}};
    }

    void propose_coarse_refinement(const std::vector<int> &marked) override {
        require(!marked.empty(), "driver passed empty coarse marking");
        calls.push_back("propose");
    }

    ReferenceEpochCandidateObservation enrich_candidate() override {
        calls.push_back("enrich");
        ReferenceEpochCandidateObservation result;
        result.eta_eq_c = 1.0;
        result.marked_c = {0};
        return result;
    }

    bool proposal_contained_in_reference() const override {
        return cycle_index >= structural.size() || !structural[cycle_index];
    }

    int minimum_reference_level_gap() const override {
        return reference_level_gap;
    }

    ReferenceEpochDualObservation candidate_dual_check(double) override {
        calls.push_back("dual");
        if (dual_index >= dual_gaps.size())
            throw std::runtime_error("dual script exhausted");
        return {2.0 * dual_gaps[dual_index], dual_gaps[dual_index++]};
    }

    void commit_coarse_refinement() override {
        calls.push_back("commit");
        ++cycle_index;
    }

    void refresh_reference(int minimum_post_refresh_level_gap) override {
        calls.push_back(
            "refresh:" + std::to_string(minimum_post_refresh_level_gap));
        if (stop_inside_refresh) {
            throw ReferenceEpochWorkLimitExceeded(
                "maximum_candidate_unknowns reached");
        }
        if (reserve_unavailable_inside_refresh) {
            throw ReferenceEpochReserveUnavailable(
                "hybrid graded reserve unavailable");
        }
        if (minimum_post_refresh_level_gap > 0)
            reference_level_gap = minimum_post_refresh_level_gap;
    }

    ReferenceEpochResourceSnapshot resources() const override {
        return snapshot;
    }
};

std::size_t action_count(
    const ReferenceEpochDriverResult &result,
    ReferenceEpochDriverAction action) {
    return std::count_if(
        result.journal.begin(), result.journal.end(),
        [&](const ReferenceEpochDriverRecord &record) {
            return record.action == action;
        });
}

ReferenceEpochDriverConfig base_config() {
    ReferenceEpochDriverConfig config;
    config.ell0 = 2;
    config.ell_max = 4;
    config.tau_loc = 0.5;
    config.q_dual = 0.5;
    config.m_dual = 3;
    config.tau_ep = 0.5;
    config.tolerance_reference = 0.01;
    config.limits.maximum_H_steps = 8;
    config.limits.maximum_epochs = 4;
    config.limits.maximum_dual_checks = 4;
    return config;
}

void verify_corrector_loop_lazy_skip_and_forced_termination() {
    ScriptedBackend backend;
    backend.corrector_bounds = {1.0, 0.1, 0.1};
    backend.error_bounds = {0.2, 0.005};
    backend.dual_gaps = {0.0};
    backend.structural = {false, false};
    ReferenceEpochPracticalDriver driver(backend, base_config());
    const ReferenceEpochDriverResult result = driver.run();
    require(result.state == ReferenceEpochDriverState::Converged,
            "forced termination dual check did not converge");
    require(action_count(result, ReferenceEpochDriverAction::IncreaseGlobalEll) == 1,
            "corrector failure did not produce exactly one ell increase");
    const auto increase = std::find_if(
        result.journal.begin(), result.journal.end(),
        [](const ReferenceEpochDriverRecord &record) {
            return record.action
                == ReferenceEpochDriverAction::IncreaseGlobalEll;
        });
    require(increase != result.journal.end() && increase->ell == 2
                && increase->theta_loc == 1.0,
            "failed corrector observation was labeled with the next ell");
    require(action_count(result, ReferenceEpochDriverAction::SkipCandidateDual) == 1,
            "lazy dual path did not skip the first check");
    require(result.dual_checks == 1 && result.H_steps == 2,
            "forced termination check/solve counts are wrong");
    require(std::find(backend.calls.begin(), backend.calls.end(), "refresh")
                == backend.calls.end(),
            "termination without gap unexpectedly refreshed reference");
}

void verify_corrector_performance_diagnostics_are_copied() {
    ScriptedBackend backend;
    backend.corrector_bounds = {0.1};
    backend.error_bounds = {0.2};
    backend.structural = {false};
    ReferenceEpochDriverConfig config = base_config();
    config.limits.maximum_H_steps = 1;
    ReferenceEpochPracticalDriver driver(backend, config);
    const ReferenceEpochDriverResult result = driver.run();
    const auto corrector = std::find_if(
        result.journal.begin(), result.journal.end(),
        [](const ReferenceEpochDriverRecord &record) {
            return record.action == ReferenceEpochDriverAction::AcceptCorrector;
        });
    require(corrector != result.journal.end(),
            "accepted corrector record is missing");
    require(corrector->corrector_parallel_threads == 3
                && corrector->corrector_patch_assembly_work_seconds == 1.25
                && corrector->corrector_patch_solve_work_seconds == 2.5
                && corrector->corrector_patch_pack_work_seconds == 0.125
                && corrector->corrector_maximum_patch_dofs == 41
                && corrector->corrector_maximum_patch_constraints == 7
                && corrector->corrector_maximum_patch_rhs == 3,
            "corrector performance diagnostics were not copied to the journal");
    require(corrector->gram_factor_cache_hits == 11
                && corrector->gram_factor_cache_misses == 13,
            "Gram factor-cache diagnostics were not copied to the journal");
}

void verify_structural_refresh_precedes_commit() {
    ScriptedBackend backend;
    backend.corrector_bounds = {0.1, 0.1};
    backend.error_bounds = {0.2, 0.2};
    backend.dual_gaps = {0.0};
    backend.structural = {true};
    ReferenceEpochDriverConfig config = base_config();
    config.limits.maximum_H_steps = 2;
    ReferenceEpochPracticalDriver driver(backend, config);
    const ReferenceEpochDriverResult result = driver.run();
    require(result.state == ReferenceEpochDriverState::WorkLimitReached,
            "structural refresh horizon did not stop structurally");
    require(action_count(result, ReferenceEpochDriverAction::RefreshReference) == 1,
            "structural hierarchy trigger did not refresh reference");
    require(result.dual_checks == 0,
            "structural refresh performed a redundant candidate dual solve");
    const auto refresh = std::find(
        backend.calls.begin(), backend.calls.end(), "refresh:0");
    const auto commit = std::find(
        backend.calls.begin(), backend.calls.end(), "commit");
    require(refresh != backend.calls.end() && commit != backend.calls.end()
                && refresh < commit,
            "prospective coarse mesh committed before structural refresh");
    const auto committed = std::find_if(
        result.journal.begin(), result.journal.end(),
        [](const ReferenceEpochDriverRecord &record) {
            return record.action
                == ReferenceEpochDriverAction::CommitCoarseRefinement;
        });
    const auto next_begin = std::find_if(
        committed == result.journal.end() ? result.journal.end()
                                          : std::next(committed),
        result.journal.end(),
        [](const ReferenceEpochDriverRecord &record) {
            return record.action == ReferenceEpochDriverAction::BeginEpoch;
        });
    require(committed != result.journal.end() && committed->epoch == 0
                && next_begin != result.journal.end() && next_begin->epoch == 1,
            "refresh commit and next BeginEpoch have inconsistent ownership");
}

void verify_numerical_gap_overrides_termination() {
    ScriptedBackend backend;
    backend.corrector_bounds = {0.1, 0.1};
    backend.error_bounds = {0.005, 0.005};
    backend.dual_gaps = {0.004};
    backend.structural = {false};
    ReferenceEpochDriverConfig config = base_config();
    config.limits.maximum_H_steps = 2;
    ReferenceEpochPracticalDriver driver(backend, config);
    const ReferenceEpochDriverResult result = driver.run();
    require(result.state == ReferenceEpochDriverState::WorkLimitReached,
            "gap-triggered refresh incorrectly terminated the run");
    require(action_count(result, ReferenceEpochDriverAction::RefreshReference) == 1,
            "positive candidate gap did not refresh reference");
}

void verify_ell_is_inherited_across_refresh() {
    ScriptedBackend backend;
    backend.corrector_bounds = {1.0, 0.1, 0.1};
    backend.error_bounds = {0.2, 0.2};
    backend.dual_gaps = {0.0};
    backend.structural = {true};
    ReferenceEpochDriverConfig config = base_config();
    config.limits.maximum_H_steps = 2;
    ReferenceEpochPracticalDriver driver(backend, config);
    const ReferenceEpochDriverResult result = driver.run();
    require(result.state == ReferenceEpochDriverState::WorkLimitReached,
            "ell inheritance trajectory did not reach its fixed horizon");
    const std::vector<std::string> expected{
        "corrector:2", "corrector:3", "corrector:3"};
    std::vector<std::string> actual;
    for (const std::string &call : backend.calls) {
        if (call.rfind("corrector:", 0) == 0) actual.push_back(call);
    }
    require(actual == expected,
            "reference refresh did not inherit ell");
}

void verify_refresh_restores_target_gap() {
    ScriptedBackend backend;
    backend.corrector_bounds = {0.1, 0.1};
    backend.error_bounds = {0.2, 0.2};
    backend.structural = {true};
    backend.reference_level_gap = 2;
    ReferenceEpochDriverConfig config = base_config();
    config.reference_refresh_level_gap = 4;
    config.reference_refresh_target_gap = 6;
    config.limits.maximum_H_steps = 2;
    ReferenceEpochPracticalDriver driver(backend, config);
    const ReferenceEpochDriverResult result = driver.run();
    require(result.state == ReferenceEpochDriverState::WorkLimitReached,
            "target-gap trajectory did not reach its fixed horizon");
    require(std::find(
                backend.calls.begin(), backend.calls.end(), "refresh:6")
            != backend.calls.end(),
            "driver did not forward the post-refresh target gap");
    const auto begin = std::find_if(
        result.journal.rbegin(), result.journal.rend(),
        [](const ReferenceEpochDriverRecord &record) {
            return record.action == ReferenceEpochDriverAction::BeginEpoch;
        });
    require(begin != result.journal.rend()
                && begin->minimum_reference_level_gap == 6,
            "new epoch did not record the restored target gap");
}

void verify_refresh_requires_enough_future_solve_points() {
    ScriptedBackend backend;
    backend.corrector_bounds = {0.1, 0.1, 0.1};
    backend.error_bounds = {0.3, 0.2, 0.1};
    backend.structural = {false, false, true};
    ReferenceEpochDriverConfig config = base_config();
    config.minimum_solved_points_per_new_epoch = 4;
    config.limits.maximum_H_steps = 5;
    ReferenceEpochPracticalDriver driver(backend, config);
    const ReferenceEpochDriverResult result = driver.run();
    require(result.state == ReferenceEpochDriverState::WorkLimitReached
                && result.stop_reason
                    == "insufficient remaining H-step budget for a new reference epoch",
            "short trailing epoch was not stopped before refresh");
    require(action_count(result, ReferenceEpochDriverAction::RefreshReference)
                == 0
                && result.H_steps == 3,
            "budget guard opened or committed an unusably short epoch");
    const auto terminal = result.journal.empty()
        ? result.journal.end() : std::prev(result.journal.end());
    require(terminal != result.journal.end()
                && terminal->action
                    == ReferenceEpochDriverAction::StopWorkLimit
                && terminal->detail == result.stop_reason,
            "remaining-budget stop has no structured terminal journal row");
}

void verify_level_gap_guard_skips_redundant_dual() {
    ScriptedBackend backend;
    backend.corrector_bounds = {0.1, 0.1};
    backend.error_bounds = {0.2, 0.2};
    backend.structural = {false};
    backend.reference_level_gap = 4;
    ReferenceEpochDriverConfig config = base_config();
    config.reference_refresh_level_gap = 4;
    config.limits.maximum_H_steps = 2;
    ReferenceEpochPracticalDriver driver(backend, config);
    const ReferenceEpochDriverResult result = driver.run();
    require(result.state == ReferenceEpochDriverState::WorkLimitReached,
            "level-gap guard trajectory did not reach its fixed horizon");
    require(result.dual_checks == 0,
            "level-gap guard performed a redundant candidate dual solve");
    require(action_count(result, ReferenceEpochDriverAction::RefreshReference)
                == 1,
            "level-gap guard did not force a reference refresh at zero dual gap");
}

void verify_level_gap_guard_is_mandatory() {
    ScriptedBackend backend;
    backend.corrector_bounds = {0.1, 0.1};
    backend.error_bounds = {0.3, 0.3};
    backend.structural = {false};
    backend.reference_level_gap = 4;
    ReferenceEpochDriverConfig config = base_config();
    config.reference_refresh_level_gap = 4;
    // Retain a nonzero legacy value to prove that it cannot suppress the
    // mandatory reserve trigger of the revised paper.
    config.minimum_H_steps_per_epoch = 2;
    config.limits.maximum_H_steps = 2;
    ReferenceEpochPracticalDriver driver(backend, config);
    const ReferenceEpochDriverResult result = driver.run();
    require(result.state == ReferenceEpochDriverState::WorkLimitReached,
            "mandatory level-gap trajectory did not reach its fixed horizon");
    require(result.H_steps == 2 && result.dual_checks == 0,
            "legacy epoch minimum suppressed or dualized a reserve refresh");
    require(action_count(result, ReferenceEpochDriverAction::RefreshReference)
                == 1,
            "mandatory reserve trigger did not refresh immediately");
}

void verify_first_epoch_solve_initializes_decrease_trigger() {
    ScriptedBackend backend;
    // The fixed-work guard is reached in SolveEstimate, after the next
    // corrector check.  Supply that final accepted check so this fixture
    // exercises the lazy-dual transition rather than script exhaustion.
    backend.corrector_bounds = {0.1, 0.1, 0.1};
    backend.error_bounds = {0.2, 0.14, 0.13};
    backend.dual_gaps = {0.0};
    backend.structural = {false, false};
    ReferenceEpochDriverConfig config = base_config();
    config.q_dual = 0.75;
    config.m_dual = 100;
    config.limits.maximum_H_steps = 3;
    ReferenceEpochPracticalDriver driver(backend, config);
    const ReferenceEpochDriverResult result = driver.run();
    require(result.state == ReferenceEpochDriverState::WorkLimitReached,
            "lazy-dual initialization trajectory did not reach its horizon");
    require(result.dual_checks == 1,
            "first epoch solve did not initialize the decrease baseline");
    const bool saw_decrease = std::any_of(
        result.journal.begin(), result.journal.end(),
        [](const ReferenceEpochDriverRecord &record) {
            return record.action
                    == ReferenceEpochDriverAction::RequestCandidateDual
                && record.decrease_dual_trigger
                && !record.interval_dual_trigger;
        });
    require(saw_decrease,
            "decrease trigger was not evaluated against the first solved point");
}

void verify_solve_horizon_stops_before_postprocessing() {
    ScriptedBackend backend;
    backend.corrector_bounds = {0.1};
    backend.error_bounds = {0.2};
    backend.structural = {false};
    ReferenceEpochDriverConfig config = base_config();
    config.limits.maximum_H_steps = 1;
    ReferenceEpochPracticalDriver driver(backend, config);
    const ReferenceEpochDriverResult result = driver.run();
    require(result.state == ReferenceEpochDriverState::WorkLimitReached
                && result.stop_reason == "maximum_H_steps reached"
                && result.H_steps == 1,
            "one-point solve horizon was not reported structurally");
    require(action_count(
                result, ReferenceEpochDriverAction::SolveAndEstimate) == 1,
            "one-point solve horizon did not record exactly one solve");
    require(action_count(
                result,
                ReferenceEpochDriverAction::ProposeCoarseRefinement) == 0
                && action_count(
                    result, ReferenceEpochDriverAction::EnrichCandidate) == 0
                && action_count(
                    result,
                    ReferenceEpochDriverAction::CommitCoarseRefinement) == 0,
            "one-point solve horizon performed post-solve mesh work");
    require(std::find(backend.calls.begin(), backend.calls.end(), "propose")
                == backend.calls.end()
                && std::find(backend.calls.begin(), backend.calls.end(), "enrich")
                    == backend.calls.end()
                && std::find(backend.calls.begin(), backend.calls.end(), "commit")
                    == backend.calls.end(),
            "one-point solve horizon called a post-solve backend action");
    const auto solved_record = std::find_if(
        result.journal.begin(), result.journal.end(),
        [](const ReferenceEpochDriverRecord &record) {
            return record.action
                == ReferenceEpochDriverAction::SolveAndEstimate;
        });
    require(solved_record != result.journal.end()
                && solved_record->H_step == 0,
            "solve record no longer uses the committed mesh index");
}

void verify_resource_and_ell_limits_are_structured() {
    ScriptedBackend resource_backend;
    resource_backend.snapshot.reference_unknowns = 1000001;
    ReferenceEpochPracticalDriver resource_driver(
        resource_backend, base_config());
    const ReferenceEpochDriverResult resource = resource_driver.run();
    require(resource.state == ReferenceEpochDriverState::WorkLimitReached
                && resource.stop_reason == "maximum_reference_unknowns reached",
            "reference resource limit was not structured");

    ScriptedBackend ell_backend;
    ell_backend.corrector_bounds = {1.0};
    ReferenceEpochDriverConfig config = base_config();
    config.ell0 = 2;
    config.ell_max = 2;
    ReferenceEpochPracticalDriver ell_driver(ell_backend, config);
    const ReferenceEpochDriverResult ell = ell_driver.run();
    require(ell.state == ReferenceEpochDriverState::WorkLimitReached
                && ell.stop_reason == "corrector threshold failed at ell_max",
            "ell_max failure was not structured");
    require(!ell.journal.empty()
                && ell.journal.back().action
                    == ReferenceEpochDriverAction::StopWorkLimit
                && ell.journal.back().ell == 2
                && ell.journal.back().theta_loc == 1.0
                && ell.journal.back().delta_loc_hat == 1.0,
            "ell_max stop discarded or mislabeled the final corrector diagnostics");
}

void verify_internal_refresh_limit_is_structured() {
    ScriptedBackend backend;
    backend.corrector_bounds = {0.1};
    backend.error_bounds = {0.2};
    backend.structural = {true};
    backend.stop_inside_refresh = true;
    ReferenceEpochPracticalDriver driver(backend, base_config());
    const ReferenceEpochDriverResult result = driver.run();
    require(result.state == ReferenceEpochDriverState::WorkLimitReached
                && result.stop_reason
                    == "maximum_candidate_unknowns reached",
            "in-action refresh resource stop was reported as a failure");
    require(action_count(result, ReferenceEpochDriverAction::StopWorkLimit) == 1
                && action_count(result, ReferenceEpochDriverAction::Fail) == 0,
            "in-action refresh resource stop has the wrong journal action");
    require(std::find(backend.calls.begin(), backend.calls.end(), "commit")
                == backend.calls.end(),
            "driver committed a proposal after an internal refresh limit");
}

void verify_unavailable_reserve_is_structured() {
    ScriptedBackend backend;
    backend.corrector_bounds = {0.1};
    backend.error_bounds = {0.2};
    backend.structural = {true};
    backend.reserve_unavailable_inside_refresh = true;
    ReferenceEpochPracticalDriver driver(backend, base_config());
    const ReferenceEpochDriverResult result = driver.run();
    require(result.state == ReferenceEpochDriverState::WorkLimitReached
                && result.stop_reason
                    == "hybrid graded reserve unavailable",
            "unavailable reserve was reported as a numerical failure");
    require(action_count(result, ReferenceEpochDriverAction::StopWorkLimit) == 1
                && action_count(result, ReferenceEpochDriverAction::Fail) == 0,
            "unavailable reserve has the wrong structured journal action");
}

} // namespace

int main() {
    try {
        verify_corrector_loop_lazy_skip_and_forced_termination();
        verify_corrector_performance_diagnostics_are_copied();
        verify_structural_refresh_precedes_commit();
        verify_numerical_gap_overrides_termination();
        verify_ell_is_inherited_across_refresh();
        verify_refresh_restores_target_gap();
        verify_refresh_requires_enough_future_solve_points();
        verify_level_gap_guard_skips_redundant_dual();
        verify_level_gap_guard_is_mandatory();
        verify_first_epoch_solve_initializes_decrease_trigger();
        verify_solve_horizon_stops_before_postprocessing();
        verify_resource_and_ell_limits_are_structured();
        verify_internal_refresh_limit_is_structured();
        verify_unavailable_reserve_is_structured();
        std::cout << "Reference-epoch practical state machine passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_reference_epoch_driver failed: "
                  << error.what() << '\n';
        return 1;
    }
}

#include "helmholtz/adaptive/practical_driver.h"

#include "helmholtz/boundary.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace lod2d::helmholtz::adaptive {
namespace {

double elapsed_seconds(
    const std::chrono::steady_clock::time_point begin,
    const std::chrono::steady_clock::time_point end) {
    return std::chrono::duration<double>(end - begin).count();
}

double triangle_diameter(const TriMesh &mesh, const Triangle &triangle) {
    const Point2 &a = mesh.nodes[static_cast<std::size_t>(triangle[0])];
    const Point2 &b = mesh.nodes[static_cast<std::size_t>(triangle[1])];
    const Point2 &c = mesh.nodes[static_cast<std::size_t>(triangle[2])];
    const auto distance = [](const Point2 &x, const Point2 &y) {
        return std::hypot(x[0] - y[0], x[1] - y[1]);
    };
    return std::max({distance(a, b), distance(b, c), distance(c, a)});
}

double maximum_kappa_H(const TriMesh &mesh, const double wavenumber) {
    double value = 0.0;
    for (const Triangle &triangle : mesh.elems) {
        value = std::max(value, wavenumber * triangle_diameter(mesh, triangle));
    }
    return value;
}

std::vector<int> inadmissible_coarse_elements(
    const TriMesh &mesh,
    const double wavenumber,
    const double threshold) {
    std::vector<int> marked;
    const double guard = 64.0 * std::numeric_limits<double>::epsilon() *
                         std::max(1.0, std::abs(threshold));
    for (std::size_t element = 0; element < mesh.elems.size(); ++element) {
        const double value =
            wavenumber * triangle_diameter(mesh, mesh.elems[element]);
        if (value > threshold + guard) {
            marked.push_back(static_cast<int>(element));
        }
    }
    return marked;
}

std::vector<int> free_coarse_nodes(const TriMesh &mesh) {
    const std::vector<int> boundary = dirichlet_nodes(mesh);
    std::vector<char> is_dirichlet(mesh.nodes.size(), false);
    for (const int node : boundary) {
        is_dirichlet[static_cast<std::size_t>(node)] = true;
    }
    std::vector<int> free_nodes;
    free_nodes.reserve(mesh.nodes.size());
    for (std::size_t node = 0; node < mesh.nodes.size(); ++node) {
        if (!is_dirichlet[node]) {
            free_nodes.push_back(static_cast<int>(node));
        }
    }
    return free_nodes;
}

std::vector<double> pull_element_values(
    const std::vector<double> &reference_values,
    const std::vector<int> &ambient_parents) {
    std::vector<double> ambient_values(ambient_parents.size(), 0.0);
    for (std::size_t element = 0; element < ambient_parents.size(); ++element) {
        const int parent = ambient_parents[element];
        if (parent < 0 || static_cast<std::size_t>(parent) >= reference_values.size()) {
            throw std::invalid_argument("Ambient/reference element ancestry is invalid.");
        }
        ambient_values[element] = reference_values[static_cast<std::size_t>(parent)];
    }
    return ambient_values;
}

template <typename Enum>
[[noreturn]] void throw_unknown_enum(const char *kind, const Enum value) {
    std::ostringstream stream;
    stream << "Unknown " << kind << " value " << static_cast<int>(value) << ".";
    throw std::invalid_argument(stream.str());
}

} // namespace

const char *practical_driver_state_name(const PracticalDriverState state) {
    switch (state) {
    case PracticalDriverState::CoarseAdmissibility:
        return "CoarseAdmissibility";
    case PracticalDriverState::LocalizationCheck:
        return "LocalizationCheck";
    case PracticalDriverState::SolveAndEstimate:
        return "SolveAndEstimate";
    case PracticalDriverState::RefineCoarse:
        return "RefineCoarse";
    case PracticalDriverState::Converged:
        return "Converged";
    case PracticalDriverState::ReferenceRefreshRequired:
        return "ReferenceRefreshRequired";
    case PracticalDriverState::WorkLimitReached:
        return "WorkLimitReached";
    case PracticalDriverState::Failed:
        return "Failed";
    }
    throw_unknown_enum("practical driver state", state);
}

const char *practical_driver_action_name(const PracticalDriverAction action) {
    switch (action) {
    case PracticalDriverAction::InitializeReferenceEpoch:
        return "InitializeReferenceEpoch";
    case PracticalDriverAction::RefineCoarseForAdmissibility:
        return "RefineCoarseForAdmissibility";
    case PracticalDriverAction::AcceptCoarse:
        return "AcceptCoarse";
    case PracticalDriverAction::IncreaseGlobalEll:
        return "IncreaseGlobalEll";
    case PracticalDriverAction::AcceptLocalization:
        return "AcceptLocalization";
    case PracticalDriverAction::AcceptFixedEll:
        return "AcceptFixedEll";
    case PracticalDriverAction::SolveStandardLod:
        return "SolveStandardLod";
    case PracticalDriverAction::RefineUniformLod:
        return "RefineUniformLod";
    case PracticalDriverAction::SolveUniformFem:
        return "SolveUniformFem";
    case PracticalDriverAction::RefineUniformFem:
        return "RefineUniformFem";
    case PracticalDriverAction::FormCoarseMarking:
        return "FormCoarseMarking";
    case PracticalDriverAction::RefineCoarse:
        return "RefineCoarse";
    case PracticalDriverAction::Complete:
        return "Complete";
    case PracticalDriverAction::StopReferenceRefreshRequired:
        return "StopReferenceRefreshRequired";
    case PracticalDriverAction::StopWorkLimit:
        return "StopWorkLimit";
    case PracticalDriverAction::Fail:
        return "Fail";
    }
    throw_unknown_enum("practical driver action", action);
}

std::vector<PracticalTargetHit> extract_practical_target_hits(
    const std::vector<PracticalIterationRecord> &journal,
    const std::vector<double> &targets) {
    std::vector<PracticalTargetHit> result;
    result.reserve(targets.size());
    for (const double target : targets) {
        if (!std::isfinite(target) || !(target > 0.0)) {
            throw std::invalid_argument(
                "Practical empirical-error targets must be positive and finite.");
        }
        PracticalTargetHit hit;
        hit.target = target;
        for (std::size_t index = 0; index < journal.size(); ++index) {
            const auto error = journal[index].reference_energy_error;
            if (error && std::isfinite(*error) && *error <= target) {
                hit.journal_index = index;
                break;
            }
        }
        result.push_back(hit);
    }
    return result;
}

PracticalAdaptiveDriver::PracticalAdaptiveDriver(
    PracticalDriverProblem problem,
    PracticalDriverConfig config)
    : problem_(std::move(problem)), config_(std::move(config)), ell_(config_.ell0) {
    start_ = std::chrono::steady_clock::now();
    validate_config();
    hierarchy_ = std::make_unique<ReferenceEpochHierarchy>(
        problem_.initial_mesh,
        config_.initial_coarse_level,
        config_.reference_level);
    reference_load_ = assemble_helmholtz_load(
        hierarchy_->reference_mesh(),
        problem_.source,
        problem_.quadrature,
        problem_.quadrature_context);
}

void PracticalAdaptiveDriver::validate_config() const {
    if (problem_.initial_mesh.nodes.empty() || problem_.initial_mesh.elems.empty()) {
        throw std::invalid_argument("Practical driver initial mesh must be nonempty.");
    }
    if (!problem_.source) {
        throw std::invalid_argument("Practical driver source must be callable.");
    }
    if (config_.initial_coarse_level < 0 ||
        config_.reference_level <= config_.initial_coarse_level) {
        throw std::invalid_argument(
            "Practical driver requires 0 <= initial coarse level < reference level.");
    }
    if (config_.ell0 < 0 || config_.ell_max < config_.ell0) {
        throw std::invalid_argument("Practical driver requires 0 <= ell0 <= ell_max.");
    }
    if (config_.localization_policy
            == PracticalLocalizationPolicy::FixedGlobalEll
        && config_.ell0 != config_.ell_max) {
        throw std::invalid_argument(
            "Fixed-ell practical driver requires ell0 == ell_max.");
    }
    const auto positive_finite = [](const double value) {
        return std::isfinite(value) && value > 0.0;
    };
    if (!positive_finite(config_.wavenumber) ||
        !positive_finite(config_.boundary_beta) ||
        !positive_finite(config_.c_H) ||
        !positive_finite(config_.theta_loc) ||
        !positive_finite(config_.C0_usr) ||
        !positive_finite(config_.C1_usr) ||
        !positive_finite(config_.tolerance_reference)) {
        throw std::invalid_argument(
            "Practical driver physical parameters and tolerances must be positive and finite.");
    }
    if (!(config_.theta_H > 0.0 && config_.theta_H <= 1.0) ||
        !(config_.rho_star > 0.0 && config_.rho_star <= 1.0) ||
        !std::isfinite(config_.theta_H) || !std::isfinite(config_.rho_star)) {
        throw std::invalid_argument("Practical driver theta_H and rho_star must lie in (0,1].");
    }
    if (config_.limits.maximum_iterations == 0 ||
        config_.limits.maximum_unknowns == 0 ||
        config_.limits.maximum_coarse_elements == 0 ||
        config_.limits.maximum_ambient_elements == 0 ||
        !std::isfinite(config_.limits.maximum_wall_seconds) ||
        config_.limits.maximum_wall_seconds < 0.0) {
        throw std::invalid_argument("Practical driver work limits are invalid.");
    }
}

void PracticalAdaptiveDriver::invalidate_discrete_cache() {
    model_.reset();
    localization_.reset();
    estimator_.reset();
}

void PracticalAdaptiveDriver::append_record(PracticalIterationRecord record) {
    record.sequence = journal_.size();
    record.reference_epoch = hierarchy_->reference_epoch();
    record.ell = ell_;
    record.coarse_nodes = hierarchy_->coarse_mesh().nodes.size();
    record.reference_nodes = hierarchy_->reference_mesh().nodes.size();
    record.ambient_nodes = hierarchy_->ambient_mesh().nodes.size();
    record.coarse_elements = hierarchy_->coarse_mesh().elems.size();
    record.ambient_elements = hierarchy_->ambient_mesh().elems.size();
    record.kappa_H_max = maximum_kappa_H(
        hierarchy_->coarse_mesh(), config_.wavenumber);
    record.rho_ambient = hierarchy_->ambient_ratio();
    record.eta_H = eta_H_;
    record.theta_loc = theta_loc_;
    record.U_practical = U_practical_;
    record.time_total_cumulative_seconds =
        elapsed_seconds(start_, std::chrono::steady_clock::now());
    journal_.push_back(std::move(record));
}

bool PracticalAdaptiveDriver::work_limit_exceeded(std::string &reason) const {
    if (journal_.size() >= config_.limits.maximum_iterations) {
        reason = "maximum practical-driver iterations reached";
        return true;
    }
    if (hierarchy_->coarse_mesh().elems.size() >
        config_.limits.maximum_coarse_elements) {
        reason = "maximum coarse element count exceeded";
        return true;
    }
    if (hierarchy_->ambient_mesh().elems.size() >
        config_.limits.maximum_ambient_elements) {
        reason = "maximum ambient element count exceeded";
        return true;
    }
    if (hierarchy_->ambient_mesh().nodes.size() >
        config_.limits.maximum_unknowns) {
        reason = "maximum unknown count exceeded";
        return true;
    }
    if (config_.limits.maximum_wall_seconds > 0.0 &&
        elapsed_seconds(start_, std::chrono::steady_clock::now()) >=
            config_.limits.maximum_wall_seconds) {
        reason = "maximum wall time reached";
        return true;
    }
    return false;
}

PracticalDriverResult PracticalAdaptiveDriver::run() {
    if (!journal_.empty()) {
        throw std::logic_error("PracticalAdaptiveDriver::run may only be called once.");
    }
    PracticalIterationRecord initialization;
    initialization.state_before = PracticalDriverState::CoarseAdmissibility;
    initialization.state_after = PracticalDriverState::CoarseAdmissibility;
    initialization.action = PracticalDriverAction::InitializeReferenceEpoch;
    initialization.detail = "initialized fixed reference epoch and ambient shadow";
    append_record(std::move(initialization));

    std::string stop_reason;
    try {
        while (state_ != PracticalDriverState::Converged &&
               state_ != PracticalDriverState::ReferenceRefreshRequired &&
               state_ != PracticalDriverState::WorkLimitReached &&
               state_ != PracticalDriverState::Failed) {
            if (work_limit_exceeded(stop_reason)) {
                const PracticalDriverState before = state_;
                state_ = PracticalDriverState::WorkLimitReached;
                PracticalIterationRecord record;
                record.state_before = before;
                record.state_after = state_;
                record.action = PracticalDriverAction::StopWorkLimit;
                record.detail = stop_reason;
                append_record(std::move(record));
                break;
            }

            if (state_ == PracticalDriverState::CoarseAdmissibility) {
                const std::vector<int> marked = inadmissible_coarse_elements(
                    hierarchy_->coarse_mesh(), config_.wavenumber, config_.c_H);
                if (!marked.empty()) {
                    if (H_steps_ >= config_.limits.maximum_H_steps) {
                        stop_reason = "maximum coarse-refinement steps reached during admissibility";
                        state_ = PracticalDriverState::WorkLimitReached;
                        PracticalIterationRecord record;
                        record.state_before = PracticalDriverState::CoarseAdmissibility;
                        record.state_after = state_;
                        record.action = PracticalDriverAction::StopWorkLimit;
                        record.marked_H = marked.size();
                        record.detail = stop_reason;
                        append_record(std::move(record));
                        break;
                    }
                    const auto begin = std::chrono::steady_clock::now();
                    const ReferenceEpochRefinementResult refined =
                        hierarchy_->refine_coarse_preserving_reference(marked);
                    const auto end = std::chrono::steady_clock::now();
                    if (refined.status ==
                        ReferenceEpochRefinementStatus::ReferenceRefreshRequired) {
                        stop_reason = refined.detail;
                        state_ = PracticalDriverState::ReferenceRefreshRequired;
                        PracticalIterationRecord record;
                        record.state_before = PracticalDriverState::CoarseAdmissibility;
                        record.state_after = state_;
                        record.action = PracticalDriverAction::StopReferenceRefreshRequired;
                        record.marked_H = marked.size();
                        record.time_mesh_seconds = elapsed_seconds(begin, end);
                        record.detail = stop_reason;
                        append_record(std::move(record));
                        break;
                    }
                    if (!refined.changed()) {
                        throw std::logic_error(
                            "Nonempty admissibility marking produced no H refinement.");
                    }
                    ++H_steps_;
                    AmbientRatioEnforcementResult ambient;
                    if (config_.localization_policy
                        == PracticalLocalizationPolicy::AdaptiveGlobalEll) {
                        ambient = hierarchy_->enforce_ambient_ratio(
                            config_.rho_star);
                    }
                    invalidate_discrete_cache();
                    localization_warm_start_.resize(0);
                    PracticalIterationRecord record;
                    record.state_before = PracticalDriverState::CoarseAdmissibility;
                    record.state_after = PracticalDriverState::CoarseAdmissibility;
                    record.action = PracticalDriverAction::RefineCoarseForAdmissibility;
                    record.marked_H = marked.size();
                    record.ambient_refined_elements = ambient.refined_elements;
                    record.time_mesh_seconds =
                        elapsed_seconds(begin, std::chrono::steady_clock::now());
                    record.detail = refined.detail;
                    append_record(std::move(record));
                    continue;
                }

                const auto begin = std::chrono::steady_clock::now();
                AmbientRatioEnforcementResult ambient;
                if (config_.localization_policy
                    == PracticalLocalizationPolicy::AdaptiveGlobalEll) {
                    ambient = hierarchy_->enforce_ambient_ratio(
                        config_.rho_star);
                }
                state_ = PracticalDriverState::LocalizationCheck;
                PracticalIterationRecord record;
                record.state_before = PracticalDriverState::CoarseAdmissibility;
                record.state_after = state_;
                record.action = PracticalDriverAction::AcceptCoarse;
                record.ambient_refined_elements = ambient.refined_elements;
                record.time_mesh_seconds =
                    elapsed_seconds(begin, std::chrono::steady_clock::now());
                record.detail = "coarse admissibility and ambient-ratio gates accepted";
                append_record(std::move(record));
                continue;
            }

            if (state_ == PracticalDriverState::LocalizationCheck) {
                const auto corrector_begin = std::chrono::steady_clock::now();
                HelmholtzProblemConfig model_config;
                model_config.H = *std::min_element(
                    hierarchy_->coarse_levels().begin(),
                    hierarchy_->coarse_levels().end());
                model_config.h = config_.reference_level;
                model_config.ell = ell_;
                model_config.wavenumber = config_.wavenumber;
                model_config.boundary_beta = config_.boundary_beta;
                model_config.mode = config_.mode;
                model_config.initial_mesh = problem_.initial_mesh;
                model_config.patch_solver = config_.patch_solver;
                model_config.quadrature = problem_.quadrature;
                model_config.quadrature_context = problem_.quadrature_context;
                model_ = std::make_unique<HelmholtzLodModel>(
                    HelmholtzLodModel::build_adaptive(
                        model_config,
                        hierarchy_->coarse_mesh(),
                        hierarchy_->coarse_levels(),
                        hierarchy_->reference_mesh(),
                        hierarchy_->reference_element_levels()));
                const auto corrector_end = std::chrono::steady_clock::now();

                const auto certificate_begin = std::chrono::steady_clock::now();
                if (config_.localization_policy
                    == PracticalLocalizationPolicy::AdaptiveGlobalEll) {
                    const HelmholtzOperators ambient_operators =
                        assemble_helmholtz_operators(
                            hierarchy_->ambient_mesh(),
                            config_.wavenumber,
                            pull_element_values(
                                model_->operators().diffusion,
                                hierarchy_->ambient_parent_reference_elements()),
                            pull_element_values(
                                model_->operators().refractive_index,
                                hierarchy_->ambient_parent_reference_elements()),
                            config_.boundary_beta);
                    localization_ =
                        std::make_unique<ReferenceLocalizationCertificate>(
                            compute_reference_localization_certificate(
                                *hierarchy_,
                                model_->operators(),
                                ambient_operators,
                                model_->corrected_test_basis(),
                                free_coarse_nodes(hierarchy_->coarse_mesh()),
                                config_.riesz_solver,
                                [&] {
                                    LocalizationEigenConfig eigen =
                                        config_.localization_eigen;
                                    if (localization_warm_start_.size() > 0) {
                                        eigen.warm_start =
                                            localization_warm_start_;
                                    }
                                    return eigen;
                                }()));
                    theta_loc_ = localization_->theta_loc;
                } else {
                    localization_.reset();
                    theta_loc_ = 0.0;
                }
                const auto certificate_end = std::chrono::steady_clock::now();

                PracticalIterationRecord record;
                record.state_before = PracticalDriverState::LocalizationCheck;
                record.time_corrector_seconds =
                    elapsed_seconds(corrector_begin, corrector_end);
                record.time_certificate_seconds =
                    config_.localization_policy
                            == PracticalLocalizationPolicy::AdaptiveGlobalEll
                        ? elapsed_seconds(certificate_begin, certificate_end)
                        : 0.0;
                record.rebuilt_correctors = model_->correctors().primal.size();
                std::string operation_limit;
                if (work_limit_exceeded(operation_limit)) {
                    stop_reason = operation_limit;
                    state_ = PracticalDriverState::WorkLimitReached;
                    record.state_after = state_;
                    record.action = PracticalDriverAction::StopWorkLimit;
                    record.detail = stop_reason;
                    append_record(std::move(record));
                    break;
                }
                if (config_.localization_policy
                        == PracticalLocalizationPolicy::AdaptiveGlobalEll
                    && theta_loc_ > config_.theta_loc) {
                    if (ell_ >= config_.ell_max) {
                        stop_reason = "localization threshold failed at ell_max";
                        state_ = PracticalDriverState::WorkLimitReached;
                        record.state_after = state_;
                        record.action = PracticalDriverAction::StopWorkLimit;
                        record.detail = stop_reason;
                        append_record(std::move(record));
                        break;
                    }
                    localization_warm_start_ =
                        localization_->spectrum.dominant_vector;
                    ++ell_;
                    invalidate_discrete_cache();
                    record.state_after = PracticalDriverState::LocalizationCheck;
                    record.action = PracticalDriverAction::IncreaseGlobalEll;
                    record.detail = "localization threshold failed; increased global ell";
                    append_record(std::move(record));
                    continue;
                }
                if (localization_) {
                    localization_warm_start_ =
                        localization_->spectrum.dominant_vector;
                }
                state_ = PracticalDriverState::SolveAndEstimate;
                record.state_after = state_;
                if (config_.localization_policy
                    == PracticalLocalizationPolicy::AdaptiveGlobalEll) {
                    record.action = PracticalDriverAction::AcceptLocalization;
                    record.detail =
                        "reference localization certificate accepted";
                } else {
                    record.action = PracticalDriverAction::AcceptFixedEll;
                    record.detail =
                        "fixed global ell accepted without localization certification";
                }
                append_record(std::move(record));
                continue;
            }

            if (state_ == PracticalDriverState::SolveAndEstimate) {
                if (!model_ ||
                    (config_.localization_policy
                         == PracticalLocalizationPolicy::AdaptiveGlobalEll
                     && !localization_)) {
                    throw std::logic_error(
                        "SolveAndEstimate requires the accepted current LOD model and localization certificate.");
                }
                const auto solve_begin = std::chrono::steady_clock::now();
                const HelmholtzLodSolution solution = model_->solve_load(reference_load_);
                const ComplexVector &candidate = solution.fine_values;
                const auto solve_end = std::chrono::steady_clock::now();
                const auto estimator_begin = std::chrono::steady_clock::now();
                estimator_ = std::make_unique<ReferenceResidualRiesz>(
                    compute_reference_residual_riesz(
                        *hierarchy_,
                        model_->operators(),
                        reference_load_,
                        candidate,
                        config_.theta_H,
                        config_.riesz_solver));
                const auto estimator_end = std::chrono::steady_clock::now();
                eta_H_ = estimator_->eta;
                U_practical_ =
                    (config_.C0_usr + config_.C1_usr * theta_loc_) * eta_H_;
                if (!std::isfinite(U_practical_)) {
                    throw std::runtime_error("Practical upper indicator is nonfinite.");
                }

                PracticalIterationRecord record;
                record.state_before = PracticalDriverState::SolveAndEstimate;
                record.time_solve_seconds = elapsed_seconds(solve_begin, solve_end);
                record.time_estimator_seconds =
                    elapsed_seconds(estimator_begin, estimator_end);
                record.evaluation_candidate = candidate;
                std::string operation_limit;
                if (work_limit_exceeded(operation_limit)) {
                    stop_reason = operation_limit;
                    state_ = PracticalDriverState::WorkLimitReached;
                    record.state_after = state_;
                    record.action = PracticalDriverAction::StopWorkLimit;
                    record.detail = stop_reason;
                    append_record(std::move(record));
                    break;
                }
                if (U_practical_ <= config_.tolerance_reference) {
                    state_ = PracticalDriverState::Converged;
                    record.state_after = state_;
                    record.action = PracticalDriverAction::Complete;
                    record.detail = "practical reference tolerance reached";
                    append_record(std::move(record));
                    break;
                }
                pending_marking_ = estimator_->marked_elements;
                if (pending_marking_.empty()) {
                    throw std::runtime_error(
                        "Reference residual estimator produced empty coarse marking above tolerance.");
                }
                state_ = PracticalDriverState::RefineCoarse;
                record.state_after = state_;
                record.action = PracticalDriverAction::FormCoarseMarking;
                record.marked_H = pending_marking_.size();
                record.detail = "formed reference-residual coarse marking";
                append_record(std::move(record));
                continue;
            }

            if (state_ == PracticalDriverState::RefineCoarse) {
                if (H_steps_ >= config_.limits.maximum_H_steps) {
                    stop_reason = "maximum coarse-refinement steps reached";
                    state_ = PracticalDriverState::WorkLimitReached;
                    PracticalIterationRecord record;
                    record.state_before = PracticalDriverState::RefineCoarse;
                    record.state_after = state_;
                    record.action = PracticalDriverAction::StopWorkLimit;
                    record.marked_H = pending_marking_.size();
                    record.detail = stop_reason;
                    append_record(std::move(record));
                    break;
                }
                const std::size_t marked_count = pending_marking_.size();
                const auto begin = std::chrono::steady_clock::now();
                const ReferenceEpochRefinementResult refined =
                    hierarchy_->refine_coarse_preserving_reference(pending_marking_);
                if (refined.status ==
                    ReferenceEpochRefinementStatus::ReferenceRefreshRequired) {
                    stop_reason = refined.detail;
                    state_ = PracticalDriverState::ReferenceRefreshRequired;
                    PracticalIterationRecord record;
                    record.state_before = PracticalDriverState::RefineCoarse;
                    record.state_after = state_;
                    record.action = PracticalDriverAction::StopReferenceRefreshRequired;
                    record.marked_H = marked_count;
                    record.time_mesh_seconds =
                        elapsed_seconds(begin, std::chrono::steady_clock::now());
                    record.detail = stop_reason;
                    append_record(std::move(record));
                    break;
                }
                if (!refined.changed()) {
                    throw std::logic_error(
                        "Nonempty residual marking produced no H refinement.");
                }
                ++H_steps_;
                AmbientRatioEnforcementResult ambient;
                if (config_.localization_policy
                    == PracticalLocalizationPolicy::AdaptiveGlobalEll) {
                    ambient = hierarchy_->enforce_ambient_ratio(
                        config_.rho_star);
                }
                pending_marking_.clear();
                invalidate_discrete_cache();
                localization_warm_start_.resize(0);
                state_ = PracticalDriverState::LocalizationCheck;
                PracticalIterationRecord record;
                record.state_before = PracticalDriverState::RefineCoarse;
                record.state_after = state_;
                record.action = PracticalDriverAction::RefineCoarse;
                record.marked_H = marked_count;
                record.ambient_refined_elements = ambient.refined_elements;
                record.time_mesh_seconds =
                    elapsed_seconds(begin, std::chrono::steady_clock::now());
                record.detail = refined.detail;
                append_record(std::move(record));
                continue;
            }

            throw std::logic_error("Practical driver reached an unhandled state.");
        }
    } catch (const std::exception &error) {
        const PracticalDriverState before = state_;
        state_ = PracticalDriverState::Failed;
        stop_reason = error.what();
        PracticalIterationRecord record;
        record.state_before = before;
        record.state_after = state_;
        record.action = PracticalDriverAction::Fail;
        record.detail = stop_reason;
        append_record(std::move(record));
    }

    if (stop_reason.empty()) {
        stop_reason = state_ == PracticalDriverState::Converged
                          ? "practical reference tolerance reached"
                          : practical_driver_state_name(state_);
    }
    PracticalDriverResult result;
    result.state = state_;
    result.stop_reason = stop_reason;
    result.ell = ell_;
    result.H_steps = H_steps_;
    result.eta_H = eta_H_;
    result.theta_loc = theta_loc_;
    result.U_practical = U_practical_;
    result.final_marked_H = pending_marking_;
    if (estimator_ && estimator_->element_eta_squared.size() ==
                          hierarchy_->coarse_mesh().elems.size()) {
        result.final_element_eta_squared = estimator_->element_eta_squared;
    }
    result.journal = journal_;
    return result;
}

} // namespace lod2d::helmholtz::adaptive

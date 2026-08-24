#include "helmholtz/adaptive/practical_driver.h"
#include "helmholtz/adaptive/estimator.h"
#include "helmholtz/benchmarks/paper_cases.h"
#include "helmholtz/boundary.h"
#include "helmholtz/experiments/paper_config.h"
#include "helmholtz/experiments/reference_solution_cache.h"
#include "helmholtz/experiments/reference_epoch_runner.h"
#include "helmholtz/operators.h"
#include "io/vtk_writer.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <optional>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

using namespace lod2d;
using namespace lod2d::helmholtz;
using namespace lod2d::helmholtz::adaptive;
using namespace lod2d::helmholtz::benchmarks;
using namespace lod2d::helmholtz::experiments;

namespace {

struct Arguments {
    std::filesystem::path config;
    std::filesystem::path output_directory;
    std::filesystem::path manuscript_baseline;
    std::filesystem::path reference_cache_directory;
    bool check = false;
    bool validate_only = false;
};

struct PaperExecution {
    PracticalDriverResult result;
    TriMesh final_mesh;
};

struct PostprocessErrors {
    std::optional<double> relative_reference_energy;
    std::optional<double> relative_reference_L2;
    std::optional<double> exact_energy;
    std::optional<double> exact_L2;
    std::optional<double> relative_exact_energy;
    std::optional<double> relative_exact_L2;
};

struct EvaluationReference {
    TriMesh mesh;
    HelmholtzOperators operators;
    ComplexVector solution;
    std::optional<HelmholtzError> exact_norm;
    bool cache_hit = false;
    std::string cache_key;
};

const char *patch_solver_label(const HelmholtzPatchSolverKind kind) {
    switch (kind) {
    case HelmholtzPatchSolverKind::DirectSaddle:
        return "direct_saddle";
    case HelmholtzPatchSolverKind::DirectSchur:
        return "direct_schur";
    case HelmholtzPatchSolverKind::ShiftedGmres:
        return "shifted_gmres";
    }
    return "unknown";
}

std::string read_text(const std::filesystem::path &path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("cannot open input file: " + path.string());
    return std::string(
        std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>());
}

Arguments parse_arguments(const int argc, char **argv) {
    Arguments result;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        const auto value = [&](const std::string &prefix) {
            return argument.substr(prefix.size());
        };
        if (argument.starts_with("--config=")) {
            result.config = value("--config=");
        } else if (argument.starts_with("--output-dir=")) {
            result.output_directory = value("--output-dir=");
        } else if (argument.starts_with("--manuscript-baseline=")) {
            result.manuscript_baseline = value("--manuscript-baseline=");
        } else if (argument.starts_with("--reference-cache-dir=")) {
            result.reference_cache_directory = value("--reference-cache-dir=");
        } else if (argument == "--check") {
            result.check = true;
        } else if (argument == "--validate-only") {
            result.validate_only = true;
        } else {
            throw std::invalid_argument("unknown argument: " + argument);
        }
    }
    if (result.config.empty() || result.manuscript_baseline.empty()
        || (!result.validate_only && result.output_directory.empty())) {
        throw std::invalid_argument(
            "required: --config=FILE --output-dir=DIR --manuscript-baseline=FILE "
            "[--reference-cache-dir=DIR] [--validate-only]");
    }
    if (!result.validate_only && result.reference_cache_directory.empty())
        result.reference_cache_directory =
            result.output_directory / "_reference_cache";
    return result;
}

std::string first_token(const std::string &text) {
    std::istringstream stream(text);
    std::string token;
    stream >> token;
    return token;
}

std::string json_string(const std::string &value) {
    std::ostringstream out;
    out << '"';
    for (const unsigned char character : value) {
        switch (character) {
        case '"': out << "\\\""; break;
        case '\\': out << "\\\\"; break;
        case '\n': out << "\\n"; break;
        case '\r': out << "\\r"; break;
        case '\t': out << "\\t"; break;
        default:
            if (character < 0x20) {
                throw std::invalid_argument("JSON string contains a control character");
            }
            out << static_cast<char>(character);
        }
    }
    out << '"';
    return out.str();
}

std::string csv_string(const std::string &value) {
    if (value.find_first_of(",\"\r\n") == std::string::npos) return value;
    std::string escaped = "\"";
    for (const char character : value) {
        if (character == '"') escaped += '"';
        escaped += character;
    }
    escaped += '"';
    return escaped;
}

template <class Value>
std::string numeric(const Value value) {
    std::ostringstream out;
    out << std::setprecision(17) << value;
    return out.str();
}

std::string optional_numeric(const std::optional<double> &value) {
    return value ? numeric(*value) : std::string{};
}

double elapsed_seconds(
    const std::chrono::steady_clock::time_point begin,
    const std::chrono::steady_clock::time_point end) {
    return std::chrono::duration<double>(end - begin).count();
}

double peak_memory_mb() {
#if defined(__unix__) || defined(__APPLE__)
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
#if defined(__APPLE__)
        return static_cast<double>(usage.ru_maxrss) / (1024.0 * 1024.0);
#else
        return static_cast<double>(usage.ru_maxrss) / 1024.0;
#endif
    }
#endif
    return 0.0;
}

void stream_evaluation_candidate(
    PracticalIterationRecord &record,
    const ReferenceEpochHierarchy &hierarchy,
    const PracticalEvaluationSink &sink) {
    if (sink && record.evaluation_candidate.size() > 0) {
        sink(record.sequence, hierarchy, record.evaluation_candidate);
        record.evaluation_candidate.resize(0);
    }
}

PaperExecution run_uniform_fem_trajectory(
    const PracticalPaperConfig &config,
    const PaperCaseData &data,
    const PracticalEvaluationSink &evaluation_sink,
    const double &evaluation_seconds_excluded) {
    ReferenceEpochHierarchy hierarchy(
        data.initial_mesh, config.initial_coarse_level, config.reference_level,
        config.reference_epoch);
    PaperExecution execution;
    std::size_t refinements = 0;
    const auto start = std::chrono::steady_clock::now();
    auto cumulative_seconds = [&] {
        return std::max(
            0.0, elapsed_seconds(start, std::chrono::steady_clock::now())
                - evaluation_seconds_excluded);
    };
    auto fill_mesh_fields = [&](PracticalIterationRecord &record) {
        record.H_step = refinements;
        record.reference_epoch = hierarchy.reference_epoch();
        record.coarse_nodes = hierarchy.coarse_mesh().nodes.size();
        record.coarse_dofs = record.coarse_nodes
            - dirichlet_nodes(hierarchy.coarse_mesh()).size();
        record.reference_nodes = hierarchy.reference_mesh().nodes.size();
        record.reference_dofs = record.reference_nodes
            - dirichlet_nodes(hierarchy.reference_mesh()).size();
        record.coarse_elements = hierarchy.coarse_mesh().elems.size();
        record.ell = 0;
        record.kappa_H_max = config.wavenumber
            * max_element_diameter(hierarchy.coarse_mesh());
        record.time_total_cumulative_seconds = cumulative_seconds();
    };

    PracticalIterationRecord initialization;
    initialization.sequence = execution.result.journal.size();
    initialization.state_before = PracticalDriverState::CoarseAdmissibility;
    initialization.state_after = PracticalDriverState::SolveAndEstimate;
    initialization.action = PracticalDriverAction::InitializeReferenceEpoch;
    initialization.detail =
        "initialized fixed evaluation-reference epoch for uniform FEM";
    fill_mesh_fields(initialization);
    execution.result.journal.push_back(std::move(initialization));

    while (true) {
        const double wall_seconds = cumulative_seconds();
        const bool wall_limited = config.work_limits.maximum_wall_seconds > 0.0
            && wall_seconds >= config.work_limits.maximum_wall_seconds;
        if (execution.result.journal.size()
                >= config.work_limits.maximum_iterations
            || hierarchy.coarse_mesh().nodes.size()
                > config.work_limits.maximum_unknowns
            || hierarchy.coarse_mesh().elems.size()
                > config.work_limits.maximum_coarse_elements
            || wall_limited) {
            execution.result.state = PracticalDriverState::WorkLimitReached;
            execution.result.stop_reason =
                "uniform FEM trajectory reached a configured work limit";
            break;
        }

        const auto solve_begin = std::chrono::steady_clock::now();
        const HelmholtzOperators operators = assemble_helmholtz_operators(
            hierarchy.coarse_mesh(), config.wavenumber, {}, {},
            config.boundary_beta);
        const ComplexVector load = assemble_helmholtz_load(
            hierarchy.coarse_mesh(), data.source, config.quadrature,
            data.quadrature_context);
        const ComplexVector solution = solve_helmholtz_fem(operators, load);
        const ComplexVector candidate =
            hierarchy.coarse_to_reference().cast<Complex>() * solution;
        const auto solve_end = std::chrono::steady_clock::now();

        PracticalIterationRecord solved;
        solved.sequence = execution.result.journal.size();
        solved.state_before = PracticalDriverState::SolveAndEstimate;
        solved.state_after = refinements >= config.work_limits.maximum_H_steps
            ? (config.trajectory_policy
                       == PracticalTrajectoryPolicy::FixedWorkHorizon
                   ? PracticalDriverState::TrajectoryComplete
                   : PracticalDriverState::WorkLimitReached)
            : PracticalDriverState::RefineCoarse;
        solved.action = PracticalDriverAction::SolveUniformFem;
        solved.evaluation_candidate = candidate;
        solved.time_solve_seconds = elapsed_seconds(solve_begin, solve_end);
        solved.detail = "solved conforming P1 FEM on the current uniform mesh";
        fill_mesh_fields(solved);
        stream_evaluation_candidate(solved, hierarchy, evaluation_sink);
        execution.result.journal.push_back(std::move(solved));

        if (refinements >= config.work_limits.maximum_H_steps) {
            execution.result.state = config.trajectory_policy
                    == PracticalTrajectoryPolicy::FixedWorkHorizon
                ? PracticalDriverState::TrajectoryComplete
                : PracticalDriverState::WorkLimitReached;
            execution.result.stop_reason = config.trajectory_policy
                    == PracticalTrajectoryPolicy::FixedWorkHorizon
                ? "fixed H-step trajectory complete"
                : "completed configured uniform FEM refinement trajectory";
            break;
        }

        std::vector<int> marked(hierarchy.coarse_mesh().elems.size());
        std::iota(marked.begin(), marked.end(), 0);
        const auto mesh_begin = std::chrono::steady_clock::now();
        const ReferenceEpochRefinementResult refined =
            hierarchy.refine_coarse_preserving_reference(marked);
        const auto mesh_end = std::chrono::steady_clock::now();
        if (refined.status
            == ReferenceEpochRefinementStatus::ReferenceRefreshRequired) {
            execution.result.state =
                PracticalDriverState::ReferenceRefreshRequired;
            execution.result.stop_reason = refined.detail;
            break;
        }
        if (!refined.changed()) {
            throw std::runtime_error(
                "uniform FEM marking produced no coarse refinement");
        }
        ++refinements;
        PracticalIterationRecord refined_record;
        refined_record.sequence = execution.result.journal.size();
        refined_record.state_before = PracticalDriverState::RefineCoarse;
        refined_record.state_after = PracticalDriverState::SolveAndEstimate;
        refined_record.action = PracticalDriverAction::RefineUniformFem;
        refined_record.marked_H = marked.size();
        refined_record.time_mesh_seconds = elapsed_seconds(mesh_begin, mesh_end);
        refined_record.detail = "uniformly refined every conforming P1 element";
        fill_mesh_fields(refined_record);
        execution.result.journal.push_back(std::move(refined_record));
    }
    execution.result.H_steps = refinements;
    execution.final_mesh = hierarchy.coarse_mesh();
    return execution;
}

PaperExecution run_adaptive_fem_trajectory(
    const PracticalPaperConfig &config,
    const PaperCaseData &data,
    const PracticalEvaluationSink &evaluation_sink,
    const std::function<void(
        std::size_t, const TriMesh &, const ComplexVector &)> &
        detached_exact_evaluation_sink,
    const double &evaluation_seconds_excluded) {
    const bool exact_only = data.exact && data.exact_gradient;
    std::unique_ptr<ReferenceEpochHierarchy> hierarchy;
    std::optional<TriMesh> detached_mesh;
    std::vector<int> detached_levels;
    if (exact_only) {
        RefineOutput initial =
            refine_mesh_nvb(data.initial_mesh, config.initial_coarse_level);
        detached_mesh = std::move(initial.mesh);
        detached_levels.assign(
            detached_mesh->elems.size(), config.initial_coarse_level);
    } else {
        hierarchy = std::make_unique<ReferenceEpochHierarchy>(
            data.initial_mesh, config.initial_coarse_level,
            config.reference_level, config.reference_epoch);
    }
    PaperExecution execution;
    std::size_t refinements = 0;
    const auto start = std::chrono::steady_clock::now();
    auto cumulative_seconds = [&] {
        return std::max(
            0.0, elapsed_seconds(start, std::chrono::steady_clock::now())
                - evaluation_seconds_excluded);
    };
    auto current_mesh = [&]() -> const TriMesh & {
        return detached_mesh ? *detached_mesh : hierarchy->coarse_mesh();
    };
    auto fill_mesh_fields = [&](PracticalIterationRecord &record) {
        record.H_step = refinements;
        record.reference_epoch = detached_mesh
            ? config.reference_epoch : hierarchy->reference_epoch();
        record.coarse_nodes = current_mesh().nodes.size();
        record.coarse_dofs = record.coarse_nodes
            - dirichlet_nodes(current_mesh()).size();
        record.reference_nodes = detached_mesh
            ? 0 : hierarchy->reference_mesh().nodes.size();
        record.coarse_elements = current_mesh().elems.size();
        record.ell = 0;
        record.kappa_H_max = config.wavenumber
            * max_element_diameter(current_mesh());
        record.time_total_cumulative_seconds = cumulative_seconds();
    };

    PracticalIterationRecord initialization;
    initialization.sequence = execution.result.journal.size();
    initialization.state_before = PracticalDriverState::CoarseAdmissibility;
    initialization.state_after = PracticalDriverState::SolveAndEstimate;
    initialization.action = PracticalDriverAction::InitializeReferenceEpoch;
    initialization.detail = exact_only
        ? "initialized manufactured-exact-only conforming P1 AFEM trajectory"
        : "initialized fixed evaluation-reference epoch for conforming P1 AFEM";
    fill_mesh_fields(initialization);
    execution.result.journal.push_back(std::move(initialization));

    while (true) {
        const double wall_seconds = cumulative_seconds();
        const bool wall_limited = config.work_limits.maximum_wall_seconds > 0.0
            && wall_seconds >= config.work_limits.maximum_wall_seconds;
        const bool iteration_limited = execution.result.journal.size()
            >= config.work_limits.maximum_iterations;
        const bool unknown_limited = current_mesh().nodes.size()
            > config.work_limits.maximum_unknowns;
        const bool coarse_element_limited = current_mesh().elems.size()
            > config.work_limits.maximum_coarse_elements;
        if (iteration_limited || unknown_limited || coarse_element_limited
            || wall_limited) {
            execution.result.state = PracticalDriverState::WorkLimitReached;
            if (iteration_limited) {
                execution.result.stop_reason =
                    "adaptive FEM trajectory reached maximum_iterations";
            } else if (unknown_limited) {
                execution.result.stop_reason =
                    "adaptive FEM trajectory reached maximum_unknowns";
            } else if (coarse_element_limited) {
                execution.result.stop_reason =
                    "adaptive FEM trajectory reached maximum_coarse_elements";
            } else {
                execution.result.stop_reason =
                    "adaptive FEM trajectory reached maximum_wall_seconds";
            }
            execution.result.final_element_eta_squared.clear();
            break;
        }

        const auto solve_begin = std::chrono::steady_clock::now();
        const HelmholtzOperators operators = assemble_helmholtz_operators(
            current_mesh(), config.wavenumber, {}, {},
            config.boundary_beta);
        const ComplexVector load = assemble_helmholtz_load(
            current_mesh(), data.source, config.quadrature,
            data.quadrature_context);
        const ComplexVector solution = solve_helmholtz_fem(operators, load);
        const ComplexVector candidate = detached_mesh
            ? ComplexVector{}
            : hierarchy->coarse_to_reference().cast<Complex>() * solution;
        const auto solve_end = std::chrono::steady_clock::now();

        const auto estimator_begin = std::chrono::steady_clock::now();
        const diagnostics::HelmholtzP1ResidualEstimate estimate =
            diagnostics::estimate_conforming_p1_residual(
                current_mesh(), operators, solution, load,
                data.source, config.quadrature, data.quadrature_context);
        std::vector<int> marked;
        if (estimate.eta > 0.0) {
            marked = mark_doerfler(estimate.element_squared, config.theta_H);
        }
        const auto estimator_end = std::chrono::steady_clock::now();

        execution.result.eta_H = estimate.eta;
        execution.result.final_marked_H = marked;
        execution.result.final_element_eta_squared = estimate.element_squared;
        const bool maximum_detached_level_reached = detached_mesh
            && !detached_levels.empty()
            && *std::max_element(
                detached_levels.begin(), detached_levels.end())
                >= config.reference_level;
        const bool completed_trajectory = maximum_detached_level_reached
            || refinements >= config.work_limits.maximum_H_steps
            || marked.empty();
        PracticalIterationRecord solved;
        solved.sequence = execution.result.journal.size();
        solved.state_before = PracticalDriverState::SolveAndEstimate;
        solved.state_after = completed_trajectory
            ? (config.trajectory_policy
                       == PracticalTrajectoryPolicy::FixedWorkHorizon
                   ? PracticalDriverState::TrajectoryComplete
                   : PracticalDriverState::WorkLimitReached)
            : PracticalDriverState::RefineCoarse;
        solved.action = PracticalDriverAction::SolveAdaptiveFem;
        solved.evaluation_candidate = candidate;
        solved.eta_H = estimate.eta;
        solved.marked_H = marked.size();
        solved.time_solve_seconds = elapsed_seconds(solve_begin, solve_end);
        solved.time_estimator_seconds =
            elapsed_seconds(estimator_begin, estimator_end);
        solved.detail =
            "solved conforming P1 FEM and formed body/jump/impedance residual marking"
            "; algebraic residual difference="
            + numeric(estimate.algebraic_relative_difference);
        fill_mesh_fields(solved);
        if (detached_mesh) {
            if (!detached_exact_evaluation_sink) {
                throw std::logic_error(
                    "detached AFEM requires a manufactured-exact evaluation sink");
            }
            detached_exact_evaluation_sink(
                solved.sequence, current_mesh(), solution);
        } else {
            stream_evaluation_candidate(solved, *hierarchy, evaluation_sink);
        }
        execution.result.journal.push_back(std::move(solved));

        if (completed_trajectory) {
            execution.result.state = config.trajectory_policy
                    == PracticalTrajectoryPolicy::FixedWorkHorizon
                ? PracticalDriverState::TrajectoryComplete
                : PracticalDriverState::WorkLimitReached;
            execution.result.stop_reason = config.trajectory_policy
                    == PracticalTrajectoryPolicy::FixedWorkHorizon
                ? (maximum_detached_level_reached
                    ? "maximum adaptive FEM level reached"
                    : "fixed H-step trajectory complete")
                : (marked.empty()
                    ? "adaptive FEM residual indicator vanished"
                    : "completed configured adaptive FEM refinement trajectory");
            break;
        }

        if (!detached_mesh && config.minimum_reference_level_gap > 0
            && config.reference_level
                    - *std::max_element(
                        hierarchy->coarse_levels().begin(),
                        hierarchy->coarse_levels().end())
                <= config.minimum_reference_level_gap) {
            execution.result.state = PracticalDriverState::TrajectoryComplete;
            execution.result.stop_reason =
                "minimum reference/coarse level gap reached";
            execution.result.journal.back().state_after =
                PracticalDriverState::TrajectoryComplete;
            execution.result.journal.back().detail +=
                "; minimum reference/coarse level gap reached";
            break;
        }

        const auto mesh_begin = std::chrono::steady_clock::now();
        ReferenceEpochRefinementResult refined;
        if (detached_mesh) {
            const TriMesh parent_mesh = *detached_mesh;
            const std::vector<int> parent_levels = detached_levels;
            RefineOutput detached =
                bisect_newest_vertex(parent_mesh, marked);
            detached_levels = refinement_child_levels(
                parent_mesh, parent_levels, detached);
            detached_mesh = std::move(detached.mesh);
            refined.status = ReferenceEpochRefinementStatus::Refined;
            refined.detail =
                "locally refined detached exact-error AFEM mesh";
        } else {
            refined = hierarchy->refine_coarse_preserving_reference(marked);
        }
        const auto mesh_end = std::chrono::steady_clock::now();
        if (refined.status
            == ReferenceEpochRefinementStatus::ReferenceRefreshRequired) {
            execution.result.state =
                PracticalDriverState::ReferenceRefreshRequired;
            execution.result.stop_reason = refined.detail;
            break;
        }
        if (!refined.changed()) {
            throw std::runtime_error(
                "adaptive FEM Doerfler marking produced no coarse refinement");
        }
        ++refinements;
        execution.result.final_element_eta_squared.clear();
        PracticalIterationRecord refined_record;
        refined_record.sequence = execution.result.journal.size();
        refined_record.state_before = PracticalDriverState::RefineCoarse;
        refined_record.state_after = PracticalDriverState::SolveAndEstimate;
        refined_record.action = PracticalDriverAction::RefineAdaptiveFem;
        refined_record.marked_H = marked.size();
        refined_record.time_mesh_seconds = elapsed_seconds(mesh_begin, mesh_end);
        refined_record.detail =
            "locally refined the conforming P1 mesh by residual Doerfler marking";
        fill_mesh_fields(refined_record);
        execution.result.journal.push_back(std::move(refined_record));
    }
    execution.result.H_steps = refinements;
    execution.final_mesh = current_mesh();
    return execution;
}

PaperExecution run_standard_lod_trajectory(
    const PracticalPaperConfig &config,
    const PaperCaseData &data,
    const PracticalEvaluationSink &evaluation_sink,
    const double &evaluation_seconds_excluded) {
    const int prior_ell = config.ell0;
    const auto start = std::chrono::steady_clock::now();
    const auto hierarchy_begin = start;
    ReferenceEpochHierarchy hierarchy(
        data.initial_mesh, config.initial_coarse_level, config.reference_level,
        config.reference_epoch);
    const auto hierarchy_end = std::chrono::steady_clock::now();
    const auto initial_load_begin = hierarchy_end;
    ComplexVector reference_load = assemble_helmholtz_load(
        hierarchy.reference_mesh(), data.source, config.quadrature,
        data.quadrature_context);
    const auto initial_load_end = std::chrono::steady_clock::now();
    const double fixed_fine_to_coarse_ratio = hierarchy.ambient_ratio();
    PaperExecution execution;
    std::size_t refinements = 0;
    auto cumulative_seconds = [&] {
        return std::max(
            0.0, elapsed_seconds(start, std::chrono::steady_clock::now())
                - evaluation_seconds_excluded);
    };
    auto fill_mesh_fields = [&](PracticalIterationRecord &record) {
        record.H_step = refinements;
        record.reference_epoch = hierarchy.reference_epoch();
        record.coarse_nodes = hierarchy.coarse_mesh().nodes.size();
        record.coarse_dofs = record.coarse_nodes
            - dirichlet_nodes(hierarchy.coarse_mesh()).size();
        record.reference_nodes = hierarchy.reference_mesh().nodes.size();
        record.reference_dofs = record.reference_nodes
            - dirichlet_nodes(hierarchy.reference_mesh()).size();
        record.coarse_elements = hierarchy.coarse_mesh().elems.size();
        record.ell = prior_ell;
        record.kappa_H_max = config.wavenumber
            * max_element_diameter(hierarchy.coarse_mesh());
        record.time_total_cumulative_seconds = cumulative_seconds();
    };

    PracticalIterationRecord initialization;
    initialization.sequence = execution.result.journal.size();
    initialization.state_before = PracticalDriverState::CoarseAdmissibility;
    initialization.state_after = PracticalDriverState::SolveAndEstimate;
    initialization.action = PracticalDriverAction::InitializeReferenceEpoch;
    initialization.detail =
        "initialized standard LOD with a fixed fine-to-coarse mesh ratio";
    initialization.time_mesh_seconds = elapsed_seconds(
        hierarchy_begin, hierarchy_end);
    initialization.time_load_assembly_seconds = elapsed_seconds(
        initial_load_begin, initial_load_end);
    fill_mesh_fields(initialization);
    execution.result.journal.push_back(std::move(initialization));

    while (true) {
        const double wall_seconds = cumulative_seconds();
        const bool wall_limited = config.work_limits.maximum_wall_seconds > 0.0
            && wall_seconds >= config.work_limits.maximum_wall_seconds;
        if (execution.result.journal.size()
                >= config.work_limits.maximum_iterations
            || hierarchy.reference_mesh().nodes.size()
                > config.work_limits.maximum_unknowns
            || hierarchy.coarse_mesh().elems.size()
                > config.work_limits.maximum_coarse_elements
            || wall_limited) {
            execution.result.state = PracticalDriverState::WorkLimitReached;
            execution.result.stop_reason =
                "standard LOD trajectory reached a configured work limit";
            break;
        }

        HelmholtzProblemConfig model_config;
        model_config.H = *std::min_element(
            hierarchy.coarse_levels().begin(), hierarchy.coarse_levels().end());
        model_config.h = config.reference_level;
        model_config.ell = prior_ell;
        model_config.wavenumber = config.wavenumber;
        model_config.boundary_beta = config.boundary_beta;
        model_config.mode = config.petrov_mode;
        model_config.initial_mesh = data.initial_mesh;
        const std::size_t reference_dofs =
            hierarchy.reference_mesh().nodes.size()
            - dirichlet_nodes(hierarchy.reference_mesh()).size();
        const bool auto_direct_schur =
            config.slod_direct_schur_min_reference_dofs > 0
            && reference_dofs
                >= config.slod_direct_schur_min_reference_dofs;
        model_config.patch_solver.kind = auto_direct_schur
            ? HelmholtzPatchSolverKind::DirectSchur
            : config.patch_solver_kind;
        model_config.patch_solver.symbolic_cache_slots =
            config.patch_symbolic_cache_slots;
        model_config.patch_solver.reuse_identical_factorization =
            config.patch_reuse_identical_factorization;
        model_config.patch_solver.maximum_parallel_solves =
            config.maximum_patch_threads;
        model_config.patch_solver.gmres.relative_tolerance =
            config.tolerances.linear_relative_residual;
        model_config.quadrature = config.quadrature;
        model_config.quadrature_context = data.quadrature_context;
        if (std::getenv("LOD2D_PROGRESS") != nullptr) {
            std::cerr
                << "LOD2D_SLOD_PROGRESS stage=model_begin"
                << " H_step=" << refinements
                << " coarse_nodes=" << hierarchy.coarse_mesh().nodes.size()
                << " reference_nodes=" << hierarchy.reference_mesh().nodes.size()
                << " reference_dofs=" << reference_dofs
                << " solver="
                << patch_solver_label(model_config.patch_solver.kind)
                << " auto_direct_schur=" << (auto_direct_schur ? 1 : 0)
                << std::endl;
        }
        const auto build_begin = std::chrono::steady_clock::now();
        HelmholtzLodModel model = HelmholtzLodModel::build_adaptive(
            model_config, hierarchy);
        const auto build_end = std::chrono::steady_clock::now();
        if (std::getenv("LOD2D_PROFILE_MODEL_STAGES") != nullptr) {
            const HelmholtzBuildTimings &timings = model.build_timings();
            const HelmholtzCorrectorDiagnostics &diagnostics =
                model.correctors().diagnostics;
            std::cerr
                << "LOD2D_MODEL_STAGES method=SLOD"
                << " refinement=" << refinements
                << " coarse_nodes=" << hierarchy.coarse_mesh().nodes.size()
                << " reference_nodes=" << hierarchy.reference_mesh().nodes.size()
                << " patch_count=" << diagnostics.patch_count
                << " parallel_threads=" << diagnostics.parallel_threads
                << " solver="
                << patch_solver_label(model_config.patch_solver.kind)
                << " auto_direct_schur=" << (auto_direct_schur ? 1 : 0)
                << " symbolic_analyses=" << diagnostics.symbolic_analyses
                << " symbolic_reuses=" << diagnostics.symbolic_reuses
                << " factorization_reuses="
                << diagnostics.factorization_reuses
                << " max_patch_dofs=" << diagnostics.maximum_patch_dofs
                << " max_patch_constraints="
                << diagnostics.maximum_patch_constraints
                << " patch_assembly_work_ms="
                << 1000.0 * diagnostics.patch_assembly_work_seconds
                << " patch_solve_work_ms="
                << 1000.0 * diagnostics.patch_solve_work_seconds
                << " patch_pack_work_ms="
                << 1000.0 * diagnostics.patch_pack_work_seconds
                << " mesh_interpolation_ms="
                << timings.mesh_and_interpolation_ms
                << " operators_ms=" << timings.operators_ms
                << " correctors_ms=" << timings.correctors_ms
                << " corrected_basis_ms=" << timings.corrected_basis_ms
                << " coarse_operator_ms=" << timings.coarse_operator_ms
                << " coarse_factorization_ms="
                << timings.coarse_factorization_ms
                << " total_ms=" << timings.total_ms
                << '\n';
        }

        const auto solve_begin = std::chrono::steady_clock::now();
        const ComplexVector candidate = model.solve_load(reference_load).fine_values;
        const auto solve_end = std::chrono::steady_clock::now();

        PracticalIterationRecord solved;
        solved.sequence = execution.result.journal.size();
        solved.state_before = PracticalDriverState::SolveAndEstimate;
        solved.state_after = refinements >= config.work_limits.maximum_H_steps
            ? (config.trajectory_policy
                       == PracticalTrajectoryPolicy::FixedWorkHorizon
                   ? PracticalDriverState::TrajectoryComplete
                   : PracticalDriverState::WorkLimitReached)
            : PracticalDriverState::RefineCoarse;
        solved.action = PracticalDriverAction::SolveStandardLod;
        solved.evaluation_candidate = candidate;
        solved.rebuilt_correctors = model.correctors().primal.size();
        const HelmholtzBuildTimings &timings = model.build_timings();
        const HelmholtzCorrectorDiagnostics &diagnostics =
            model.correctors().diagnostics;
        solved.patch_solver_kind_used = model_config.patch_solver.kind;
        solved.slod_auto_direct_schur = auto_direct_schur;
        solved.corrector_parallel_threads = diagnostics.parallel_threads;
        solved.corrector_symbolic_analyses = diagnostics.symbolic_analyses;
        solved.corrector_symbolic_reuses = diagnostics.symbolic_reuses;
        solved.corrector_factorization_reuses =
            diagnostics.factorization_reuses;
        solved.corrector_maximum_patch_dofs =
            diagnostics.maximum_patch_dofs;
        solved.corrector_maximum_patch_constraints =
            diagnostics.maximum_patch_constraints;
        solved.corrector_patch_assembly_work_seconds =
            diagnostics.patch_assembly_work_seconds;
        solved.corrector_patch_solve_work_seconds =
            diagnostics.patch_solve_work_seconds;
        solved.corrector_patch_pack_work_seconds =
            diagnostics.patch_pack_work_seconds;
        solved.time_model_mesh_interpolation_seconds =
            timings.mesh_and_interpolation_ms / 1000.0;
        solved.time_operator_assembly_seconds =
            timings.operators_ms / 1000.0;
        solved.time_corrector_seconds =
            timings.correctors_ms / 1000.0;
        solved.time_basis_assembly_seconds =
            timings.corrected_basis_ms / 1000.0;
        solved.time_coarse_operator_seconds =
            timings.coarse_operator_ms / 1000.0;
        solved.time_coarse_factorization_seconds =
            timings.coarse_factorization_ms / 1000.0;
        solved.time_model_total_seconds =
            elapsed_seconds(build_begin, build_end);
        solved.time_solve_seconds = elapsed_seconds(solve_begin, solve_end);
        solved.detail =
            "rebuilt and solved standard LOD on the current uniform coarse mesh"
            "; full model build seconds="
            + numeric(elapsed_seconds(build_begin, build_end));
        fill_mesh_fields(solved);
        stream_evaluation_candidate(solved, hierarchy, evaluation_sink);
        execution.result.journal.push_back(std::move(solved));

        if (refinements >= config.work_limits.maximum_H_steps) {
            execution.result.state = config.trajectory_policy
                    == PracticalTrajectoryPolicy::FixedWorkHorizon
                ? PracticalDriverState::TrajectoryComplete
                : PracticalDriverState::WorkLimitReached;
            execution.result.stop_reason = config.trajectory_policy
                    == PracticalTrajectoryPolicy::FixedWorkHorizon
                ? "fixed H-step trajectory complete"
                : "completed configured standard LOD refinement trajectory";
            break;
        }

        std::vector<int> marked(hierarchy.coarse_mesh().elems.size());
        std::iota(marked.begin(), marked.end(), 0);
        const int coarse_level_before = *std::min_element(
            hierarchy.coarse_levels().begin(), hierarchy.coarse_levels().end());
        const int fine_level_before = *std::min_element(
            hierarchy.reference_element_levels().begin(),
            hierarchy.reference_element_levels().end());
        const auto mesh_begin = std::chrono::steady_clock::now();
        const ReferenceEpochRefinementResult refined =
            hierarchy.refine_coarse_preserving_reference(marked);
        if (refined.status
            == ReferenceEpochRefinementStatus::ReferenceRefreshRequired) {
            execution.result.state =
                PracticalDriverState::ReferenceRefreshRequired;
            execution.result.stop_reason = refined.detail;
            break;
        }
        if (!refined.changed()) {
            throw std::runtime_error(
                "standard LOD uniform marking produced no coarse refinement");
        }
        const AmbientRatioEnforcementResult fine_refined =
            hierarchy.enforce_ambient_ratio(fixed_fine_to_coarse_ratio);
        if (!fine_refined.changed) {
            throw std::runtime_error(
                "standard LOD H refinement did not trigger the matching fine refinement");
        }
        hierarchy.refresh_reference_from_ambient();
        const int coarse_level_after = *std::min_element(
            hierarchy.coarse_levels().begin(), hierarchy.coarse_levels().end());
        const int fine_level_after = *std::min_element(
            hierarchy.reference_element_levels().begin(),
            hierarchy.reference_element_levels().end());
        if (coarse_level_after != coarse_level_before + 1
            || fine_level_after != fine_level_before + 1
            || fine_level_after - coarse_level_after
                != fine_level_before - coarse_level_before) {
            throw std::runtime_error(
                "standard LOD failed to preserve its fixed h/H level difference");
        }
        const auto load_begin = std::chrono::steady_clock::now();
        reference_load = assemble_helmholtz_load(
            hierarchy.reference_mesh(), data.source, config.quadrature,
            data.quadrature_context);
        const auto load_end = std::chrono::steady_clock::now();
        const auto mesh_end = std::chrono::steady_clock::now();
        ++refinements;
        PracticalIterationRecord refined_record;
        refined_record.sequence = execution.result.journal.size();
        refined_record.state_before = PracticalDriverState::RefineCoarse;
        refined_record.state_after = PracticalDriverState::SolveAndEstimate;
        refined_record.action = PracticalDriverAction::RefineUniformLod;
        refined_record.marked_H = marked.size();
        refined_record.ambient_refined_elements =
            fine_refined.refined_elements;
        refined_record.time_load_assembly_seconds =
            elapsed_seconds(load_begin, load_end);
        refined_record.time_mesh_seconds = std::max(
            0.0, elapsed_seconds(mesh_begin, mesh_end)
                - refined_record.time_load_assembly_seconds);
        refined_record.detail =
            "uniformly refined H and h by one level while preserving fixed h/H";
        fill_mesh_fields(refined_record);
        execution.result.journal.push_back(std::move(refined_record));
    }
    execution.result.ell = prior_ell;
    execution.result.H_steps = refinements;
    execution.final_mesh = hierarchy.coarse_mesh();
    return execution;
}

std::pair<double, double> relative_reference_errors(
    const HelmholtzOperators &operators,
    const ComplexVector &reference,
    const ComplexVector &candidate) {
    if (reference.size() != candidate.size() || !reference.allFinite() ||
        !candidate.allFinite()) {
        throw std::invalid_argument("evaluation reference and PALOD candidate mismatch");
    }
    const ComplexVector difference = reference - candidate;
    const Eigen::SparseMatrix<double> energy = operators.stiffness +
        (operators.wavenumber * operators.wavenumber) * operators.mass;
    const auto square = [](const ComplexVector &values,
                           const Eigen::SparseMatrix<double> &matrix) {
        return std::max(0.0, std::real(values.dot(matrix.cast<Complex>() * values)));
    };
    const double energy_denominator = square(reference, energy);
    const double l2_denominator = square(reference, operators.mass);
    if (!(energy_denominator > 0.0) || !(l2_denominator > 0.0)) {
        throw std::runtime_error("evaluation reference has a zero norm");
    }
    return {
        std::sqrt(square(difference, energy) / energy_denominator),
        std::sqrt(square(difference, operators.mass) / l2_denominator),
    };
}

PostprocessErrors postprocess_errors(
    const PracticalPaperConfig &config,
    const PaperCaseData &data,
    const TriMesh &reference_mesh,
    const HelmholtzOperators &reference_operators,
    const ComplexVector &reference_solution,
    const std::optional<HelmholtzError> &exact_norm,
    const ComplexVector &candidate,
    double &reference_error_seconds,
    double &exact_error_seconds) {
    const auto reference_begin = std::chrono::steady_clock::now();
    const auto reference = relative_reference_errors(
        reference_operators, reference_solution, candidate);
    reference_error_seconds += elapsed_seconds(
        reference_begin, std::chrono::steady_clock::now());
    PostprocessErrors result;
    result.relative_reference_energy = reference.first;
    result.relative_reference_L2 = reference.second;
    if (!exact_norm) return result;

    const auto exact_begin = std::chrono::steady_clock::now();
    const HelmholtzError exact_error = compute_helmholtz_error(
        reference_mesh, candidate, config.wavenumber,
        data.exact, data.exact_gradient,
        config.quadrature, data.quadrature_context);
    if (!(exact_norm->energy > 0.0) || !(exact_norm->l2 > 0.0)) {
        throw std::runtime_error(
            "manufactured exact solution has a zero norm");
    }
    result.exact_energy = exact_error.energy;
    result.exact_L2 = exact_error.l2;
    result.relative_exact_energy = exact_error.energy / exact_norm->energy;
    result.relative_exact_L2 = exact_error.l2 / exact_norm->l2;
    exact_error_seconds += elapsed_seconds(
        exact_begin, std::chrono::steady_clock::now());
    return result;
}

void assign_postprocess_errors(
    PracticalIterationRecord &record,
    const PostprocessErrors &errors) {
    record.reference_energy_error = errors.relative_reference_energy;
    record.reference_L2_error = errors.relative_reference_L2;
    record.exact_energy_error = errors.exact_energy;
    record.exact_L2_error = errors.exact_L2;
    record.relative_exact_energy_error = errors.relative_exact_energy;
    record.relative_exact_L2_error = errors.relative_exact_L2;
}

PostprocessErrors postprocess_exact_errors(
    const PracticalPaperConfig &config,
    const PaperCaseData &data,
    const TriMesh &mesh,
    const ComplexVector &candidate,
    double &exact_error_seconds) {
    if (!data.exact || !data.exact_gradient) {
        throw std::logic_error(
            "exact-only postprocessing requires a manufactured solution");
    }
    const auto begin = std::chrono::steady_clock::now();
    const HelmholtzError candidate_error = compute_helmholtz_error(
        mesh, candidate, config.wavenumber, data.exact, data.exact_gradient,
        config.quadrature, data.quadrature_context);
    const HelmholtzError exact_norm = compute_helmholtz_error(
        mesh, ComplexVector::Zero(candidate.size()), config.wavenumber,
        data.exact, data.exact_gradient, config.quadrature,
        data.quadrature_context);
    if (!(exact_norm.energy > 0.0) || !(exact_norm.l2 > 0.0)) {
        throw std::runtime_error(
            "manufactured exact solution has a zero norm");
    }
    PostprocessErrors result;
    result.exact_energy = candidate_error.energy;
    result.exact_L2 = candidate_error.l2;
    result.relative_exact_energy = candidate_error.energy / exact_norm.energy;
    result.relative_exact_L2 = candidate_error.l2 / exact_norm.l2;
    exact_error_seconds += elapsed_seconds(
        begin, std::chrono::steady_clock::now());
    return result;
}

PaperRunStatus paper_status(const PracticalDriverState state,
                            const std::string &reason) {
    switch (state) {
    case PracticalDriverState::Converged:
    case PracticalDriverState::TrajectoryComplete:
        return PaperRunStatus::Success;
    case PracticalDriverState::WorkLimitReached:
        if (reason.find("wall") != std::string::npos)
            return PaperRunStatus::CensoredTimeLimit;
        if (reason.find("iteration") != std::string::npos)
            return PaperRunStatus::CensoredIterationLimit;
        return PaperRunStatus::CensoredWorkLimit;
    case PracticalDriverState::ReferenceRefreshRequired:
        return PaperRunStatus::Interrupted;
    case PracticalDriverState::Failed:
        return PaperRunStatus::LinearAlgebraFailure;
    default:
        return PaperRunStatus::Interrupted;
    }
}

void write_iterations(
    const std::filesystem::path &path,
    const PracticalPaperConfig &config,
    const std::string &run_id,
    const PracticalDriverResult &result,
    const double peak_mb) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "schema_version,case,method,kappa,run_id,reference_epoch,H_step,iteration,action,stop_reason,"
           "N_H,DoF_H,N_ref,DoF_ref,N_amb,ell,kappa_H_max,mu_H,rho_amb,eta_H,Theta_loc,"
           "localization_eigen_iterations,localization_eigen_relative_residual,"
           "localization_sparse_generalized,localization_used_warm_start,"
           "localization_patch_threads,localization_ambient_patch_count,"
           "localization_ambient_patch_factorizations,"
           "localization_ambient_rhs_solves,"
           "localization_ambient_max_active_columns,"
           "patch_solver_used,slod_auto_direct_schur,"
           "corrector_parallel_threads,corrector_symbolic_analyses,"
           "corrector_symbolic_reuses,corrector_factorization_reuses,"
           "corrector_maximum_patch_dofs,corrector_maximum_patch_constraints,"
           "corrector_patch_assembly_work,corrector_patch_solve_work,"
           "corrector_patch_pack_work,U_prac,"
           "reference_energy_error,reference_L2_error,exact_energy_error,exact_L2_error,"
           "relative_exact_energy_error,relative_exact_L2_error,reference_error_ratio,"
           "reference_log_improvement,reference_dof_log_slope,marked_H,rebuilt_correctors,"
           "ambient_refined_elements,time_mesh,time_load_assembly,"
           "time_model_mesh_interpolation,time_operator_assembly,time_corrector,"
           "time_basis_assembly,time_coarse_operator,time_coarse_factorization,"
           "time_model_total,time_localization_ambient_operator_assembly,"
           "time_localization_retraction,time_localization_defect_rhs,"
           "time_localization_ambient_riesz,"
           "time_localization_ambient_patch_solve,"
           "time_localization_ambient_gram_reduction,"
           "time_localization_coarse_energy,time_localization_spectrum,"
           "time_certificate,time_solve,estimator_patch_threads,"
           "estimator_patch_factorizations,time_estimator_prepare,"
           "time_estimator_patch_solve,time_estimator_reduction,"
           "time_estimator,time_total_cumulative,peak_memory_mb\n";
    const bool ell_method = config.method_id == PracticalPaperMethod::Palod
        || config.method_id == PracticalPaperMethod::HlodFixed
        || config.method_id == PracticalPaperMethod::Slod;
    const bool practical_bound_method =
        config.method_id == PracticalPaperMethod::Palod
        || config.method_id == PracticalPaperMethod::HlodFixed;
    const bool eta_method = practical_bound_method
        || config.method_id == PracticalPaperMethod::Afem;
    const bool localization_method =
        config.method_id == PracticalPaperMethod::Palod;
    const bool slod_method = config.method_id == PracticalPaperMethod::Slod;
    std::optional<std::pair<std::size_t, double>> previous_reference_point;
    for (const PracticalIterationRecord &record : result.journal) {
        const bool lod_model_profile =
            record.action == PracticalDriverAction::AcceptLocalization
            || record.action == PracticalDriverAction::AcceptFixedEll
            || record.action == PracticalDriverAction::SolveStandardLod;
        std::optional<double> error_ratio;
        std::optional<double> log_improvement;
        std::optional<double> dof_log_slope;
        if (record.reference_energy_error && previous_reference_point
            && record.coarse_nodes > previous_reference_point->first
            && *record.reference_energy_error > 0.0
            && previous_reference_point->second > 0.0) {
            error_ratio = *record.reference_energy_error
                / previous_reference_point->second;
            if (std::isfinite(*error_ratio) && *error_ratio > 0.0) {
                log_improvement = -std::log(*error_ratio);
                dof_log_slope = *log_improvement
                    / std::log(static_cast<double>(record.coarse_nodes)
                               / static_cast<double>(
                                   previous_reference_point->first));
            }
        }
        out << config.schema_version << ',' << to_string(config.case_id) << ','
            << to_string(config.method_id) << ',' << numeric(config.wavenumber) << ','
            << run_id << ',' << record.reference_epoch << ',' << record.H_step
            << ',' << record.sequence << ','
            << practical_driver_action_name(record.action) << ','
            << csv_string(record.detail) << ',' << record.coarse_nodes << ','
            << record.coarse_dofs << ',' << record.reference_nodes << ','
            << (record.reference_dofs > 0
                    ? std::to_string(record.reference_dofs) : "") << ','
            << (localization_method ? std::to_string(record.ambient_nodes) : "")
            << ',' << (ell_method ? std::to_string(record.ell) : "") << ','
            << numeric(record.kappa_H_max) << ",,"
            << (localization_method ? numeric(record.rho_ambient) : "") << ','
            << (eta_method ? numeric(record.eta_H) : "") << ','
            << (localization_method ? numeric(record.theta_loc) : "") << ','
            << (localization_method
                    ? std::to_string(record.localization_eigen_iterations)
                    : "") << ','
            << (localization_method
                    ? numeric(record.localization_eigen_relative_residual)
                    : "") << ','
            << (localization_method
                    ? std::to_string(
                        record.localization_sparse_generalized ? 1 : 0)
                    : "") << ','
            << (localization_method
                    ? std::to_string(
                        record.localization_used_warm_start ? 1 : 0)
                    : "") << ','
            << (localization_method
                    ? std::to_string(record.localization_patch_threads)
                    : "") << ','
            << (localization_method
                    ? std::to_string(record.localization_ambient_patch_count)
                    : "") << ','
            << (localization_method
                    ? std::to_string(
                        record.localization_ambient_patch_factorizations)
                    : "") << ','
            << (localization_method
                    ? std::to_string(record.localization_ambient_rhs_solves)
                    : "") << ','
            << (localization_method
                    ? std::to_string(
                        record.localization_ambient_max_active_columns)
                    : "") << ','
            << (lod_model_profile
                    ? patch_solver_label(record.patch_solver_kind_used) : "")
            << ','
            << (slod_method
                    && record.action == PracticalDriverAction::SolveStandardLod
                    ? std::to_string(record.slod_auto_direct_schur ? 1 : 0)
                    : "") << ','
            << (lod_model_profile
                    ? std::to_string(record.corrector_parallel_threads) : "")
            << ','
            << (lod_model_profile
                    ? std::to_string(record.corrector_symbolic_analyses) : "")
            << ','
            << (lod_model_profile
                    ? std::to_string(record.corrector_symbolic_reuses) : "")
            << ','
            << (lod_model_profile
                    ? std::to_string(record.corrector_factorization_reuses) : "")
            << ','
            << (lod_model_profile
                    ? std::to_string(record.corrector_maximum_patch_dofs) : "")
            << ','
            << (lod_model_profile
                    ? std::to_string(
                        record.corrector_maximum_patch_constraints) : "")
            << ','
            << (lod_model_profile
                    ? numeric(record.corrector_patch_assembly_work_seconds) : "")
            << ','
            << (lod_model_profile
                    ? numeric(record.corrector_patch_solve_work_seconds) : "")
            << ','
            << (lod_model_profile
                    ? numeric(record.corrector_patch_pack_work_seconds) : "")
            << ','
            << (practical_bound_method ? numeric(record.U_practical) : "") << ','
            << optional_numeric(record.reference_energy_error) << ','
            << optional_numeric(record.reference_L2_error) << ','
            << optional_numeric(record.exact_energy_error) << ','
            << optional_numeric(record.exact_L2_error) << ','
            << optional_numeric(record.relative_exact_energy_error) << ','
            << optional_numeric(record.relative_exact_L2_error) << ','
            << optional_numeric(error_ratio) << ','
            << optional_numeric(log_improvement) << ','
            << optional_numeric(dof_log_slope) << ','
            << record.marked_H << ',' << record.rebuilt_correctors << ','
            << record.ambient_refined_elements << ','
            << numeric(record.time_mesh_seconds) << ','
            << numeric(record.time_load_assembly_seconds) << ','
            << numeric(record.time_model_mesh_interpolation_seconds) << ','
            << numeric(record.time_operator_assembly_seconds) << ','
            << numeric(record.time_corrector_seconds) << ','
            << numeric(record.time_basis_assembly_seconds) << ','
            << numeric(record.time_coarse_operator_seconds) << ','
            << numeric(record.time_coarse_factorization_seconds) << ','
            << numeric(record.time_model_total_seconds) << ','
            << numeric(
                record.time_localization_ambient_operator_assembly_seconds)
            << ','
            << numeric(record.time_localization_retraction_seconds) << ','
            << numeric(record.time_localization_defect_rhs_seconds) << ','
            << numeric(record.time_localization_ambient_riesz_seconds) << ','
            << numeric(
                record.time_localization_ambient_patch_solve_seconds) << ','
            << numeric(
                record.time_localization_ambient_gram_reduction_seconds) << ','
            << numeric(record.time_localization_coarse_energy_seconds) << ','
            << numeric(record.time_localization_spectrum_seconds) << ','
            << numeric(record.time_certificate_seconds) << ','
            << numeric(record.time_solve_seconds) << ','
            << record.estimator_patch_threads << ','
            << record.estimator_patch_factorizations << ','
            << numeric(record.time_estimator_prepare_seconds) << ','
            << numeric(record.time_estimator_patch_solve_seconds) << ','
            << numeric(record.time_estimator_reduction_seconds) << ','
            << numeric(record.time_estimator_seconds) << ','
            << numeric(record.time_total_cumulative_seconds) << ','
            << numeric(peak_mb) << '\n';
        if (record.reference_energy_error && record.coarse_nodes > 0
            && std::isfinite(*record.reference_energy_error)
            && *record.reference_energy_error > 0.0) {
            previous_reference_point = {
                record.coarse_nodes, *record.reference_energy_error};
        }
    }
}

void write_summary(
    const std::filesystem::path &path,
    const PracticalPaperConfig &config,
    const PracticalDriverResult &result) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "target,status,first_iteration,reference_energy_error,U_prac,N_H,DoF_H,ell,"
           "rebuilt_correctors_cumulative,ambient_refined_elements_cumulative,"
           "time_total_cumulative\n";
    const std::vector<double> targets(
        config.relative_energy_targets.begin(),
        config.relative_energy_targets.end());
    for (const PracticalTargetHit &target_hit :
         extract_practical_target_hits(result.journal, targets)) {
        out << numeric(target_hit.target) << ',';
        if (!target_hit.journal_index) {
            out << "not_reached,,,,,,,,,\n";
        } else {
            const auto hit = result.journal.begin() +
                static_cast<std::ptrdiff_t>(*target_hit.journal_index);
            std::size_t rebuilt = 0;
            std::size_t ambient_refined = 0;
            for (auto record = result.journal.begin(); record != std::next(hit); ++record) {
                rebuilt += record->rebuilt_correctors;
                ambient_refined += record->ambient_refined_elements;
            }
            const bool adaptive_lod_method =
                config.method_id == PracticalPaperMethod::Palod
                || config.method_id == PracticalPaperMethod::HlodFixed;
            const bool ell_method = adaptive_lod_method
                || config.method_id == PracticalPaperMethod::Slod;
            out << "reached," << hit->sequence << ','
                << numeric(*hit->reference_energy_error) << ','
                << (adaptive_lod_method ? numeric(hit->U_practical) : "") << ','
                << hit->coarse_nodes << ',' << hit->coarse_dofs << ','
                << (ell_method ? std::to_string(hit->ell) : "") << ','
                << rebuilt << ',' << ambient_refined << ','
                << numeric(hit->time_total_cumulative_seconds)
                << '\n';
        }
    }
}

void write_ell_history(
    const std::filesystem::path &path,
    const PracticalPaperConfig &config,
    const PracticalDriverResult &result) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "H_step,iteration,action,ell,Theta_loc,rebuilt_correctors,"
           "time_model_total,time_model_mesh_interpolation,time_operator_assembly,"
           "time_corrector,time_basis_assembly,time_coarse_operator,"
           "time_coarse_factorization,time_localization_ambient_operator_assembly,"
           "time_localization_retraction,time_localization_defect_rhs,"
           "time_localization_ambient_riesz,"
           "time_localization_ambient_patch_solve,"
           "time_localization_ambient_gram_reduction,"
           "time_localization_coarse_energy,time_localization_spectrum,"
           "time_certificate\n";
    for (const PracticalIterationRecord &record : result.journal) {
        if (record.state_before != PracticalDriverState::LocalizationCheck
            && record.action != PracticalDriverAction::SolveStandardLod) {
            continue;
        }
        out << record.H_step << ',' << record.sequence << ','
            << practical_driver_action_name(record.action) << ','
            << record.ell << ','
            << (config.method_id == PracticalPaperMethod::Palod
                    ? numeric(record.theta_loc) : "") << ','
            << record.rebuilt_correctors << ','
            << numeric(record.time_model_total_seconds) << ','
            << numeric(record.time_model_mesh_interpolation_seconds) << ','
            << numeric(record.time_operator_assembly_seconds) << ','
            << numeric(record.time_corrector_seconds) << ','
            << numeric(record.time_basis_assembly_seconds) << ','
            << numeric(record.time_coarse_operator_seconds) << ','
            << numeric(record.time_coarse_factorization_seconds) << ','
            << numeric(
                record.time_localization_ambient_operator_assembly_seconds)
            << ','
            << numeric(record.time_localization_retraction_seconds) << ','
            << numeric(record.time_localization_defect_rhs_seconds) << ','
            << numeric(record.time_localization_ambient_riesz_seconds) << ','
            << numeric(
                record.time_localization_ambient_patch_solve_seconds) << ','
            << numeric(
                record.time_localization_ambient_gram_reduction_seconds) << ','
            << numeric(record.time_localization_coarse_energy_seconds) << ','
            << numeric(record.time_localization_spectrum_seconds) << ','
            << numeric(record.time_certificate_seconds) << '\n';
    }
}

void write_final_mesh(
    const std::filesystem::path &path,
    const TriMesh &mesh,
    const PracticalDriverResult &result) {
    std::vector<double> eta(mesh.elems.size(), 0.0);
    std::vector<int> eta_available(mesh.elems.size(), 0);
    if (result.final_element_eta_squared.size() == mesh.elems.size()) {
        for (std::size_t element = 0; element < eta.size(); ++element) {
            eta[element] = std::sqrt(std::max(
                0.0, result.final_element_eta_squared[element]));
            eta_available[element] = 1;
        }
    }
    const std::array<io::VtkDoubleFieldView, 1> doubles{{{"eta_H_T", eta}}};
    const std::array<io::VtkIntFieldView, 1> integers{{{"eta_H_T_available", eta_available}}};
    io::VtuDataView fields;
    fields.cell_double = doubles;
    fields.cell_int = integers;
    io::write_vtu(path, mesh, fields);
}

void write_run_json(
    const std::filesystem::path &path,
    const PracticalPaperConfig &config,
    const std::string &run_id,
    const PracticalDriverResult &result,
    const double method_seconds,
    const double reference_seconds,
    const double exact_seconds,
    const bool reference_cache_hit,
    const std::string &reference_cache_key,
    const bool exact_only_evaluation,
    const double peak_mb) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    const PracticalConvergenceDiagnostic convergence =
        diagnose_practical_convergence_regime(
            result.journal, config.plateau_diagnostic);
    const auto optional_json_number = [](const std::optional<double> value) {
        return value ? numeric(*value) : std::string("null");
    };
    const auto accumulated = [&](double PracticalIterationRecord::*member) {
        double total = 0.0;
        for (const PracticalIterationRecord &record : result.journal)
            total += record.*member;
        return total;
    };
    out << "{\n  \"schema_version\":" << config.schema_version << ",\n"
        << "  \"run_id\":" << json_string(run_id) << ",\n"
        << "  \"config_hash\":" << json_string(canonical_config_hash(config)) << ",\n"
        << "  \"case\":" << json_string(std::string(to_string(config.case_id))) << ",\n"
        << "  \"method\":" << json_string(std::string(to_string(config.method_id))) << ",\n"
        << "  \"status\":" << json_string(std::string(to_string(
               paper_status(result.state, result.stop_reason)))) << ",\n"
        << "  \"claim\":\"implementation-study\",\n"
        << "  \"stop_reason\":" << json_string(result.stop_reason) << ",\n"
        << "  \"driver_state\":"
        << json_string(practical_driver_state_name(result.state)) << ",\n"
        << "  \"error_evaluation_mode\":"
        << json_string(exact_only_evaluation
                ? "manufactured-exact-only" : "fixed-reference")
        << ",\n"
        << "  \"config\":" << canonical_json(config) << ",\n"
        << "  \"files\":{\"iterations\":\"iterations.csv\","
           "\"summary\":\"summary.csv\",\"ell_history\":\"ell_history.csv\","
           "\"final_mesh\":\"final_mesh.vtu\"},\n"
        << "  \"convergence_diagnostic\":{\"status\":"
        << json_string(practical_convergence_regime_name(convergence.regime))
        << ",\"distinct_points\":" << convergence.distinct_points
        << ",\"criterion\":\"last two positive error-vs-DOF log slopes differ by at most factor 2\""
        << ",\"previous_log_slope\":"
        << optional_json_number(convergence.previous_log_slope)
        << ",\"last_log_slope\":"
        << optional_json_number(convergence.last_log_slope)
        << ",\"last_error_ratio\":"
        << optional_json_number(convergence.last_error_ratio)
        << ",\"last_log_improvement\":"
        << optional_json_number(convergence.last_log_improvement)
        << ",\"plateau_minimum_geometric_mean_ratio\":"
        << numeric(config.plateau_diagnostic.minimum_geometric_mean_ratio)
        << ",\"plateau_maximum_relative_oscillation\":"
        << numeric(config.plateau_diagnostic.maximum_relative_oscillation)
        << ",\"plateau_window_steps\":"
        << config.plateau_diagnostic.window_steps
        << ",\"window_geometric_mean_ratio\":"
        << optional_json_number(convergence.window_geometric_mean_ratio)
        << ",\"window_relative_oscillation\":"
        << optional_json_number(convergence.window_relative_oscillation)
        << ",\"plateau_observed\":"
        << (convergence.plateau_observed ? "true" : "false") << "},\n"
        << "  \"timing\":{\"method_seconds\":" << numeric(method_seconds)
        << ",\"evaluation_reference_seconds\":" << numeric(reference_seconds)
        << ",\"evaluation_exact_seconds\":" << numeric(exact_seconds)
        << ",\"stage_totals_seconds\":{\"mesh\":"
        << numeric(accumulated(&PracticalIterationRecord::time_mesh_seconds))
        << ",\"load_assembly\":" << numeric(accumulated(
            &PracticalIterationRecord::time_load_assembly_seconds))
        << ",\"model_total\":" << numeric(accumulated(
            &PracticalIterationRecord::time_model_total_seconds))
        << ",\"model_mesh_interpolation\":" << numeric(accumulated(
            &PracticalIterationRecord::time_model_mesh_interpolation_seconds))
        << ",\"operator_assembly\":" << numeric(accumulated(
            &PracticalIterationRecord::time_operator_assembly_seconds))
        << ",\"corrector\":" << numeric(accumulated(
            &PracticalIterationRecord::time_corrector_seconds))
        << ",\"basis_assembly\":" << numeric(accumulated(
            &PracticalIterationRecord::time_basis_assembly_seconds))
        << ",\"coarse_operator\":" << numeric(accumulated(
            &PracticalIterationRecord::time_coarse_operator_seconds))
        << ",\"coarse_factorization\":" << numeric(accumulated(
            &PracticalIterationRecord::time_coarse_factorization_seconds))
        << ",\"localization_ambient_operator_assembly\":"
        << numeric(accumulated(
            &PracticalIterationRecord::
                time_localization_ambient_operator_assembly_seconds))
        << ",\"localization_retraction\":" << numeric(accumulated(
            &PracticalIterationRecord::time_localization_retraction_seconds))
        << ",\"localization_defect_rhs\":" << numeric(accumulated(
            &PracticalIterationRecord::time_localization_defect_rhs_seconds))
        << ",\"localization_ambient_riesz\":" << numeric(accumulated(
            &PracticalIterationRecord::
                time_localization_ambient_riesz_seconds))
        << ",\"localization_ambient_patch_solve\":"
        << numeric(accumulated(
            &PracticalIterationRecord::
                time_localization_ambient_patch_solve_seconds))
        << ",\"localization_ambient_gram_reduction\":"
        << numeric(accumulated(
            &PracticalIterationRecord::
                time_localization_ambient_gram_reduction_seconds))
        << ",\"localization_coarse_energy\":" << numeric(accumulated(
            &PracticalIterationRecord::
                time_localization_coarse_energy_seconds))
        << ",\"localization_spectrum\":" << numeric(accumulated(
            &PracticalIterationRecord::time_localization_spectrum_seconds))
        << ",\"certificate_total\":" << numeric(accumulated(
            &PracticalIterationRecord::time_certificate_seconds))
        << ",\"solve\":" << numeric(accumulated(
            &PracticalIterationRecord::time_solve_seconds))
        << ",\"estimator\":" << numeric(accumulated(
            &PracticalIterationRecord::time_estimator_seconds)) << "}"
        << ",\"evaluation_reference_excluded_from_method_time\":true"
        << ",\"evaluation_exact_excluded_from_method_time\":true"
        << ",\"reference_solution_cache_hit\":"
        << (reference_cache_hit ? "true" : "false")
        << ",\"reference_solution_cache_key\":"
        << json_string(reference_cache_key) << "},\n"
        << "  \"hardware\":{\"hardware_threads\":"
        << std::thread::hardware_concurrency()
        << ",\"compiler\":" << json_string(__VERSION__)
        << ",\"peak_memory_mb\":" << numeric(peak_mb) << "}\n}\n";
}

} // namespace

int main(const int argc, char **argv) {
    try {
        const Arguments arguments = parse_arguments(argc, argv);
        const std::string config_text = read_text(arguments.config);
        if (is_reference_epoch_paper_config(config_text)) {
            return run_reference_epoch_paper(
                config_text, arguments.output_directory,
                arguments.manuscript_baseline, arguments.check,
                arguments.validate_only);
        }
        const PracticalPaperConfig config =
            parse_practical_paper_config(config_text);
        if (config.method_id != PracticalPaperMethod::Palod
            && config.method_id != PracticalPaperMethod::HlodFixed
            && config.method_id != PracticalPaperMethod::Slod
            && config.method_id != PracticalPaperMethod::Ufem
            && config.method_id != PracticalPaperMethod::Afem) {
            throw std::invalid_argument(
                "paper runner received an unsupported practical method");
        }
        const std::string frozen_hash = first_token(read_text(arguments.manuscript_baseline));
        if (frozen_hash != config.manuscript_sha256) {
            throw std::invalid_argument(
                "config manuscript_sha256 does not match the frozen baseline artifact");
        }
        const PaperCaseData data = make_paper_case(
            config.case_id, config.wavenumber,
            config.singular_oscillatory_fraction,
            config.singular_cutoff_outer_radius,
            config.singular_quintic_cutoff,
            config.smooth_wave_amplitude);
        if (arguments.validate_only) {
            std::cout << "config=valid\n"
                      << "run_id=" << make_run_id(config) << '\n';
            return 0;
        }

        const std::string reference_identity =
            "reference-sparse-lu-v1/" + config.git_commit + "/"
            + config.build_hash;
        std::map<std::uint64_t, std::unique_ptr<EvaluationReference>>
            evaluation_references;
        double reference_construction_seconds = 0.0;
        double exact_norm_seconds = 0.0;
        bool reference_cache_hit = true;
        std::string reference_cache_key = "not-applicable";
        const auto ensure_evaluation_reference =
            [&](const ReferenceEpochHierarchy &hierarchy)
                -> const EvaluationReference & {
            const std::uint64_t epoch = hierarchy.reference_epoch();
            const auto found = evaluation_references.find(epoch);
            if (found != evaluation_references.end()) return *found->second;

            const auto reference_begin = std::chrono::steady_clock::now();
            auto evaluation = std::make_unique<EvaluationReference>();
            evaluation->mesh = hierarchy.reference_mesh();
            evaluation->operators = assemble_helmholtz_operators(
                evaluation->mesh, config.wavenumber, {}, {},
                config.boundary_beta);
            const ComplexVector reference_load = assemble_helmholtz_load(
                evaluation->mesh, data.source, config.quadrature,
                data.quadrature_context);
            evaluation->cache_key = reference_solution_cache_key(
                evaluation->mesh, evaluation->operators, reference_load,
                reference_identity);
            ReferenceSolutionCacheLookup reference_cache =
                load_reference_solution_cache(
                    arguments.reference_cache_directory,
                    evaluation->cache_key, reference_load.size());
            evaluation->cache_hit = reference_cache.solution.has_value();
            if (evaluation->cache_hit) {
                evaluation->solution = std::move(*reference_cache.solution);
            } else {
                evaluation->solution = solve_helmholtz_fem(
                    evaluation->operators, reference_load);
                store_reference_solution_cache(
                    arguments.reference_cache_directory,
                    evaluation->cache_key, evaluation->solution);
            }
            reference_construction_seconds += elapsed_seconds(
                reference_begin, std::chrono::steady_clock::now());
            if (data.exact && data.exact_gradient) {
                const auto exact_norm_begin = std::chrono::steady_clock::now();
                evaluation->exact_norm = compute_helmholtz_error(
                    evaluation->mesh,
                    ComplexVector::Zero(evaluation->solution.size()),
                    config.wavenumber, data.exact, data.exact_gradient,
                    config.quadrature, data.quadrature_context);
                exact_norm_seconds += elapsed_seconds(
                    exact_norm_begin, std::chrono::steady_clock::now());
            }
            reference_cache_hit = reference_cache_hit && evaluation->cache_hit;
            reference_cache_key = evaluation->cache_key;
            const auto inserted = evaluation_references.emplace(
                epoch, std::move(evaluation));
            return *inserted.first->second;
        };

        const bool exact_only_evaluation =
            (config.method_id == PracticalPaperMethod::Afem
             || config.manufactured_exact_only_errors)
            && data.exact && data.exact_gradient;
        if (config.manufactured_exact_only_errors
            && !exact_only_evaluation) {
            throw std::invalid_argument(
                "manufactured_exact_only_errors requires a manufactured exact solution");
        }
        if (exact_only_evaluation) reference_cache_hit = false;
        if (!exact_only_evaluation) {
            ReferenceEpochHierarchy initial_evaluation_hierarchy(
                data.initial_mesh, config.initial_coarse_level,
                config.reference_level, config.reference_epoch);
            (void)ensure_evaluation_reference(initial_evaluation_hierarchy);
        }

        std::vector<std::optional<PostprocessErrors>> streamed_errors;
        double streamed_evaluation_seconds = 0.0;
        double reference_error_seconds = 0.0;
        double exact_error_seconds = 0.0;
        const PracticalEvaluationSink evaluation_sink =
            [&](const std::size_t sequence,
                const ReferenceEpochHierarchy &hierarchy,
                const ComplexVector &candidate) {
                const auto begin = std::chrono::steady_clock::now();
                if (streamed_errors.size() <= sequence)
                    streamed_errors.resize(sequence + 1);
                if (exact_only_evaluation) {
                    streamed_errors[sequence] = postprocess_exact_errors(
                        config, data, hierarchy.reference_mesh(), candidate,
                        exact_error_seconds);
                } else {
                    const EvaluationReference &evaluation =
                        ensure_evaluation_reference(hierarchy);
                    streamed_errors[sequence] = postprocess_errors(
                        config, data, evaluation.mesh, evaluation.operators,
                        evaluation.solution, evaluation.exact_norm,
                        candidate, reference_error_seconds,
                        exact_error_seconds);
                }
                streamed_evaluation_seconds += elapsed_seconds(
                    begin, std::chrono::steady_clock::now());
            };
        const auto detached_exact_evaluation_sink =
            [&](const std::size_t sequence,
                const TriMesh &mesh,
                const ComplexVector &candidate) {
                const auto begin = std::chrono::steady_clock::now();
                if (streamed_errors.size() <= sequence)
                    streamed_errors.resize(sequence + 1);
                streamed_errors[sequence] = postprocess_exact_errors(
                    config, data, mesh, candidate, exact_error_seconds);
                streamed_evaluation_seconds += elapsed_seconds(
                    begin, std::chrono::steady_clock::now());
            };

        const auto method_begin = std::chrono::steady_clock::now();
        PaperExecution execution;
        if (config.method_id == PracticalPaperMethod::Ufem) {
            execution = run_uniform_fem_trajectory(
                config, data, evaluation_sink, streamed_evaluation_seconds);
        } else if (config.method_id == PracticalPaperMethod::Afem) {
            execution = run_adaptive_fem_trajectory(
                config, data, evaluation_sink,
                detached_exact_evaluation_sink,
                streamed_evaluation_seconds);
        } else if (config.method_id == PracticalPaperMethod::Slod) {
            execution = run_standard_lod_trajectory(
                config, data, evaluation_sink, streamed_evaluation_seconds);
        } else {
            PracticalDriverProblem problem;
            problem.initial_mesh = data.initial_mesh;
            problem.source = data.source;
            problem.quadrature = config.quadrature;
            problem.quadrature_context = data.quadrature_context;
            PracticalAdaptiveDriver driver(
                std::move(problem), make_practical_driver_config(config),
                evaluation_sink);
            execution.result = driver.run();
            execution.final_mesh = driver.hierarchy().coarse_mesh();
        }
        PracticalDriverResult &result = execution.result;
        const auto method_end = std::chrono::steady_clock::now();
        for (PracticalIterationRecord &record : result.journal) {
            if (record.sequence < streamed_errors.size()
                && streamed_errors[record.sequence]) {
                assign_postprocess_errors(
                    record, *streamed_errors[record.sequence]);
            } else if (record.evaluation_candidate.size() > 0) {
                if (exact_only_evaluation) {
                    throw std::runtime_error(
                        "manufactured-exact-only evaluation was not streamed on its matching mesh");
                }
                const EvaluationReference &evaluation =
                    *evaluation_references.at(record.reference_epoch);
                const PostprocessErrors errors = postprocess_errors(
                    config, data, evaluation.mesh, evaluation.operators,
                    evaluation.solution, evaluation.exact_norm,
                    record.evaluation_candidate,
                    reference_error_seconds, exact_error_seconds);
                assign_postprocess_errors(record, errors);
                record.evaluation_candidate.resize(0);
            }
        }
        const bool fixed_empirical_trajectory =
            config.method_id == PracticalPaperMethod::Ufem
            || config.method_id == PracticalPaperMethod::Slod
            || config.method_id == PracticalPaperMethod::Afem;
        const bool fixed_work_horizon =
            config.trajectory_policy
            == PracticalTrajectoryPolicy::FixedWorkHorizon;
        if (fixed_work_horizon
            && result.state == PracticalDriverState::WorkLimitReached
            && (result.H_steps >= config.work_limits.maximum_H_steps
                || result.stop_reason
                    == "adaptive FEM residual indicator vanished")) {
            result.state = PracticalDriverState::TrajectoryComplete;
            result.stop_reason = "fixed H-step trajectory complete";
        }
        if (fixed_empirical_trajectory && !exact_only_evaluation
            && result.state == PracticalDriverState::WorkLimitReached) {
            const double smallest_target = config.relative_energy_targets.back();
            const bool reached = std::any_of(
                result.journal.begin(), result.journal.end(),
                [smallest_target](const PracticalIterationRecord &record) {
                    return record.reference_energy_error
                        && *record.reference_energy_error <= smallest_target;
                });
            if (reached) {
                result.state = PracticalDriverState::Converged;
                result.stop_reason =
                    "uniform FEM fixed trajectory reached the smallest empirical target";
            }
        }

        const std::string run_id = make_run_id(config);
        const std::filesystem::path run_directory =
            arguments.output_directory / run_id;
        std::filesystem::create_directories(run_directory);
        const double peak_mb = peak_memory_mb();
        write_iterations(run_directory / "iterations.csv", config, run_id, result, peak_mb);
        write_summary(run_directory / "summary.csv", config, result);
        write_ell_history(run_directory / "ell_history.csv", config, result);
        write_final_mesh(
            run_directory / "final_mesh.vtu", execution.final_mesh, result);
        write_run_json(
            run_directory / "run.json", config, run_id, result,
            std::max(
                0.0, elapsed_seconds(method_begin, method_end)
                    - streamed_evaluation_seconds),
            reference_construction_seconds + reference_error_seconds,
            exact_norm_seconds + exact_error_seconds,
            reference_cache_hit,
            reference_cache_key,
            exact_only_evaluation,
            peak_mb);

        if (arguments.check) {
            const bool acceptable_state =
                result.state == PracticalDriverState::Converged
                || result.state == PracticalDriverState::TrajectoryComplete
                || ((fixed_empirical_trajectory || fixed_work_horizon)
                    && result.state == PracticalDriverState::WorkLimitReached);
            if (!acceptable_state ||
                std::none_of(result.journal.begin(), result.journal.end(),
                             [](const PracticalIterationRecord &record) {
                                 return record.action == PracticalDriverAction::Complete
                                     || record.action
                                         == PracticalDriverAction::CompleteTrajectory
                                     || record.action
                                         == PracticalDriverAction::StopWorkLimit
                                     || record.action
                                         == PracticalDriverAction::SolveUniformFem
                                     || record.action
                                         == PracticalDriverAction::SolveStandardLod
                                     || record.action
                                         == PracticalDriverAction::SolveAdaptiveFem;
                             })) {
                throw std::runtime_error(
                    "paper smoke did not complete a real method chain");
            }
            const bool has_reference_error = std::any_of(
                result.journal.begin(), result.journal.end(),
                [](const PracticalIterationRecord &record) {
                    return record.reference_energy_error
                        && record.reference_L2_error;
                });
            if (has_reference_error == exact_only_evaluation) {
                throw std::runtime_error(
                    exact_only_evaluation
                        ? "manufactured-exact-only run unexpectedly used fixed-reference errors"
                        : "WP5 smoke did not compute post-run reference errors");
            }
            const bool expects_exact = data.exact && data.exact_gradient;
            const bool has_exact = std::any_of(
                result.journal.begin(), result.journal.end(),
                [](const PracticalIterationRecord &record) {
                    return record.exact_energy_error
                        && record.exact_L2_error
                        && record.relative_exact_energy_error
                        && record.relative_exact_L2_error;
                });
            if (expects_exact != has_exact) {
                throw std::runtime_error(
                    expects_exact
                        ? "manufactured paper case did not compute exact errors"
                        : "non-manufactured paper case unexpectedly reported exact errors");
            }
            const auto has_action = [&](const PracticalDriverAction action) {
                return std::any_of(
                    result.journal.begin(), result.journal.end(),
                    [action](const PracticalIterationRecord &record) {
                        return record.action == action;
                    });
            };
            if (config.method_id == PracticalPaperMethod::Palod) {
                const bool invalid_completion = fixed_work_horizon
                    ? (!has_action(PracticalDriverAction::CompleteTrajectory)
                       || has_action(PracticalDriverAction::Complete)
                       || result.H_steps
                           != config.work_limits.maximum_H_steps
                       || result.stop_reason
                           != "fixed H-step trajectory complete")
                    : !has_action(PracticalDriverAction::Complete);
                const bool missing_stage_profile = std::none_of(
                    result.journal.begin(), result.journal.end(),
                    [](const PracticalIterationRecord &record) {
                        return record.state_before
                                == PracticalDriverState::LocalizationCheck
                            && record.time_model_total_seconds > 0.0
                            && record.time_model_mesh_interpolation_seconds > 0.0
                            && record.time_operator_assembly_seconds > 0.0
                            && record.time_corrector_seconds > 0.0
                            && record.time_basis_assembly_seconds > 0.0
                            && record.time_coarse_operator_seconds > 0.0
                            && record.time_coarse_factorization_seconds > 0.0
                            && record.time_certificate_seconds > 0.0
                            && record.time_localization_retraction_seconds > 0.0
                            && record.time_localization_defect_rhs_seconds > 0.0
                            && record.time_localization_ambient_riesz_seconds > 0.0
                            && record.time_localization_ambient_patch_solve_seconds
                                > 0.0
                            && record.time_localization_coarse_energy_seconds > 0.0
                            && record.time_localization_spectrum_seconds > 0.0
                            && record.localization_ambient_patch_count > 0
                            && record.localization_ambient_rhs_solves > 0;
                    });
                if (!has_action(PracticalDriverAction::AcceptLocalization)
                    || invalid_completion || missing_stage_profile) {
                    throw std::runtime_error(
                        "PALOD smoke did not execute its trajectory with a complete stage profile");
                }
            }
            if (config.method_id == PracticalPaperMethod::HlodFixed) {
                const bool paid_certificate_cost = std::any_of(
                    result.journal.begin(), result.journal.end(),
                    [](const PracticalIterationRecord &record) {
                        return record.time_certificate_seconds != 0.0
                            || record.ambient_refined_elements != 0;
                    });
                if (!has_action(PracticalDriverAction::AcceptFixedEll)
                    || has_action(PracticalDriverAction::AcceptLocalization)
                    || has_action(PracticalDriverAction::IncreaseGlobalEll)
                    || paid_certificate_cost) {
                    throw std::runtime_error(
                        "HLOD-fixed smoke entered or paid for PALOD localization");
                }
            }
            if (config.method_id == PracticalPaperMethod::Ufem
                && (!has_action(PracticalDriverAction::SolveUniformFem)
                    || !has_action(PracticalDriverAction::RefineUniformFem)
                    || has_action(PracticalDriverAction::AcceptLocalization)
                    || has_action(PracticalDriverAction::AcceptFixedEll)
                    || has_action(PracticalDriverAction::FormCoarseMarking))) {
                throw std::runtime_error(
                    "UFEM smoke did not remain on the uniform conforming FEM path");
            }
            if (config.method_id == PracticalPaperMethod::Slod) {
                const int expected_ell = config.ell0;
                const bool wrong_ell = std::any_of(
                    result.journal.begin(), result.journal.end(),
                    [expected_ell](const PracticalIterationRecord &record) {
                        return record.action
                                   == PracticalDriverAction::SolveStandardLod
                            && record.ell != expected_ell;
                    });
                const bool missed_synchronized_fine_refinement = std::any_of(
                    result.journal.begin(), result.journal.end(),
                    [](const PracticalIterationRecord &record) {
                        return record.action
                                   == PracticalDriverAction::RefineUniformLod
                            && (record.ambient_refined_elements == 0
                                || record.reference_epoch == 0);
                    });
                const bool missed_solver_switch =
                    config.slod_direct_schur_min_reference_dofs > 0
                    && std::any_of(
                        result.journal.begin(), result.journal.end(),
                        [&](const PracticalIterationRecord &record) {
                            if (record.action
                                != PracticalDriverAction::SolveStandardLod)
                                return false;
                            const bool expected_switch = record.reference_dofs
                                >= config.slod_direct_schur_min_reference_dofs;
                            const HelmholtzPatchSolverKind expected_solver =
                                expected_switch
                                    ? HelmholtzPatchSolverKind::DirectSchur
                                    : HelmholtzPatchSolverKind::DirectSaddle;
                            return record.slod_auto_direct_schur
                                    != expected_switch
                                || record.patch_solver_kind_used
                                    != expected_solver
                                || record.corrector_parallel_threads <= 0
                                || record.corrector_maximum_patch_dofs <= 0
                                || !(record.time_corrector_seconds > 0.0);
                        });
                if (!has_action(PracticalDriverAction::SolveStandardLod)
                    || !has_action(PracticalDriverAction::RefineUniformLod)
                    || has_action(PracticalDriverAction::AcceptLocalization)
                    || has_action(PracticalDriverAction::AcceptFixedEll)
                    || has_action(PracticalDriverAction::FormCoarseMarking)
                    || wrong_ell || missed_synchronized_fine_refinement
                    || missed_solver_switch) {
                    throw std::runtime_error(
                        "SLOD smoke did not preserve its empirical ell and synchronized uniform H/h refinement");
                }
            }
            if (config.method_id == PracticalPaperMethod::Afem) {
                const bool invalid_estimate = std::any_of(
                    result.journal.begin(), result.journal.end(),
                    [](const PracticalIterationRecord &record) {
                        return record.action
                                   == PracticalDriverAction::SolveAdaptiveFem
                            && (!(record.eta_H > 0.0) || record.marked_H == 0);
                    });
                const bool paid_lod_cost = std::any_of(
                    result.journal.begin(), result.journal.end(),
                    [](const PracticalIterationRecord &record) {
                        return record.time_corrector_seconds != 0.0
                            || record.time_certificate_seconds != 0.0
                            || record.ambient_refined_elements != 0;
                    });
                if (!has_action(PracticalDriverAction::SolveAdaptiveFem)
                    || !has_action(PracticalDriverAction::RefineAdaptiveFem)
                    || has_action(PracticalDriverAction::SolveUniformFem)
                    || has_action(PracticalDriverAction::SolveStandardLod)
                    || has_action(PracticalDriverAction::AcceptLocalization)
                    || has_action(PracticalDriverAction::AcceptFixedEll)
                    || has_action(PracticalDriverAction::FormCoarseMarking)
                    || invalid_estimate || paid_lod_cost) {
                    throw std::runtime_error(
                        "AFEM smoke did not remain on the conforming residual-adaptive path");
                }
            }
            for (const char *file : {
                     "iterations.csv", "summary.csv", "ell_history.csv",
                     "final_mesh.vtu", "run.json"}) {
                if (!std::filesystem::is_regular_file(run_directory / file)) {
                    throw std::runtime_error(std::string("missing WP5 artifact: ") + file);
                }
            }
        }

        std::cout << "run_id=" << run_id << '\n'
                  << "state=" << practical_driver_state_name(result.state) << '\n'
                  << "convergence_regime="
                  << practical_convergence_regime_name(
                         diagnose_practical_convergence_regime(
                             result.journal,
                             config.plateau_diagnostic).regime)
                  << '\n'
                  << "reference_cache="
                  << (reference_cache_hit ? "hit" : "miss") << '\n'
                  << "claim=implementation-study\n"
                  << "output=" << run_directory.string() << '\n';
        return result.state == PracticalDriverState::Failed ? 2 : 0;
    } catch (const std::exception &error) {
        std::cerr << "bench_helmholtz_adaptive_paper failed: "
                  << error.what() << '\n';
        return 1;
    }
}

#include "helmholtz/adaptive/practical_driver.h"
#include "helmholtz/adaptive/estimator.h"
#include "helmholtz/benchmarks/paper_cases.h"
#include "helmholtz/experiments/paper_config.h"
#include "helmholtz/operators.h"
#include "io/vtk_writer.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
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
    bool check = false;
};

struct PaperExecution {
    PracticalDriverResult result;
    TriMesh final_mesh;
};

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
        } else if (argument == "--check") {
            result.check = true;
        } else {
            throw std::invalid_argument("unknown argument: " + argument);
        }
    }
    if (result.config.empty() || result.output_directory.empty() ||
        result.manuscript_baseline.empty()) {
        throw std::invalid_argument(
            "required: --config=FILE --output-dir=DIR --manuscript-baseline=FILE");
    }
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

PaperExecution run_uniform_fem_trajectory(
    const PracticalPaperConfig &config,
    const PaperCaseData &data) {
    ReferenceEpochHierarchy hierarchy(
        data.initial_mesh, config.initial_coarse_level, config.reference_level);
    PaperExecution execution;
    const auto start = std::chrono::steady_clock::now();
    auto cumulative_seconds = [&] {
        return elapsed_seconds(start, std::chrono::steady_clock::now());
    };
    auto fill_mesh_fields = [&](PracticalIterationRecord &record) {
        record.reference_epoch = hierarchy.reference_epoch();
        record.coarse_nodes = hierarchy.coarse_mesh().nodes.size();
        record.reference_nodes = hierarchy.reference_mesh().nodes.size();
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

    std::size_t refinements = 0;
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
            ? PracticalDriverState::WorkLimitReached
            : PracticalDriverState::RefineCoarse;
        solved.action = PracticalDriverAction::SolveUniformFem;
        solved.evaluation_candidate = candidate;
        solved.time_solve_seconds = elapsed_seconds(solve_begin, solve_end);
        solved.detail = "solved conforming P1 FEM on the current uniform mesh";
        fill_mesh_fields(solved);
        execution.result.journal.push_back(std::move(solved));

        if (refinements >= config.work_limits.maximum_H_steps) {
            execution.result.state = PracticalDriverState::WorkLimitReached;
            execution.result.stop_reason =
                "completed configured uniform FEM refinement trajectory";
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
    const PaperCaseData &data) {
    ReferenceEpochHierarchy hierarchy(
        data.initial_mesh, config.initial_coarse_level, config.reference_level);
    PaperExecution execution;
    const auto start = std::chrono::steady_clock::now();
    auto cumulative_seconds = [&] {
        return elapsed_seconds(start, std::chrono::steady_clock::now());
    };
    auto fill_mesh_fields = [&](PracticalIterationRecord &record) {
        record.reference_epoch = hierarchy.reference_epoch();
        record.coarse_nodes = hierarchy.coarse_mesh().nodes.size();
        record.reference_nodes = hierarchy.reference_mesh().nodes.size();
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
        "initialized fixed evaluation-reference epoch for conforming P1 AFEM";
    fill_mesh_fields(initialization);
    execution.result.journal.push_back(std::move(initialization));

    std::size_t refinements = 0;
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
                "adaptive FEM trajectory reached a configured work limit";
            execution.result.final_element_eta_squared.clear();
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

        const auto estimator_begin = std::chrono::steady_clock::now();
        const diagnostics::HelmholtzP1ResidualEstimate estimate =
            diagnostics::estimate_conforming_p1_residual(
                hierarchy.coarse_mesh(), operators, solution, load,
                data.source, config.quadrature, data.quadrature_context);
        std::vector<int> marked;
        if (estimate.eta > 0.0) {
            marked = mark_doerfler(estimate.element_squared, config.theta_H);
        }
        const auto estimator_end = std::chrono::steady_clock::now();

        execution.result.eta_H = estimate.eta;
        execution.result.final_marked_H = marked;
        execution.result.final_element_eta_squared = estimate.element_squared;
        const bool completed_trajectory =
            refinements >= config.work_limits.maximum_H_steps || marked.empty();
        PracticalIterationRecord solved;
        solved.sequence = execution.result.journal.size();
        solved.state_before = PracticalDriverState::SolveAndEstimate;
        solved.state_after = completed_trajectory
            ? PracticalDriverState::WorkLimitReached
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
        execution.result.journal.push_back(std::move(solved));

        if (completed_trajectory) {
            execution.result.state = PracticalDriverState::WorkLimitReached;
            execution.result.stop_reason = marked.empty()
                ? "adaptive FEM residual indicator vanished"
                : "completed configured adaptive FEM refinement trajectory";
            break;
        }

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
    execution.final_mesh = hierarchy.coarse_mesh();
    return execution;
}

PaperExecution run_standard_lod_trajectory(
    const PracticalPaperConfig &config,
    const PaperCaseData &data) {
    const int prior_ell = standard_lod_prior_ell(config.wavenumber);
    ReferenceEpochHierarchy hierarchy(
        data.initial_mesh, config.initial_coarse_level, config.reference_level);
    const ComplexVector reference_load = assemble_helmholtz_load(
        hierarchy.reference_mesh(), data.source, config.quadrature,
        data.quadrature_context);
    PaperExecution execution;
    const auto start = std::chrono::steady_clock::now();
    auto cumulative_seconds = [&] {
        return elapsed_seconds(start, std::chrono::steady_clock::now());
    };
    auto fill_mesh_fields = [&](PracticalIterationRecord &record) {
        record.reference_epoch = hierarchy.reference_epoch();
        record.coarse_nodes = hierarchy.coarse_mesh().nodes.size();
        record.reference_nodes = hierarchy.reference_mesh().nodes.size();
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
        "initialized fixed evaluation-reference epoch for uniform standard LOD";
    fill_mesh_fields(initialization);
    execution.result.journal.push_back(std::move(initialization));

    std::size_t refinements = 0;
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
        model_config.patch_solver.kind = config.patch_solver_kind;
        model_config.patch_solver.gmres.relative_tolerance =
            config.tolerances.linear_relative_residual;
        model_config.quadrature = config.quadrature;
        model_config.quadrature_context = data.quadrature_context;
        const auto build_begin = std::chrono::steady_clock::now();
        HelmholtzLodModel model = HelmholtzLodModel::build_adaptive(
            model_config,
            hierarchy.coarse_mesh(), hierarchy.coarse_levels(),
            hierarchy.reference_mesh(), hierarchy.reference_element_levels());
        const auto build_end = std::chrono::steady_clock::now();

        const auto solve_begin = std::chrono::steady_clock::now();
        const ComplexVector candidate = model.solve_load(reference_load).fine_values;
        const auto solve_end = std::chrono::steady_clock::now();

        PracticalIterationRecord solved;
        solved.sequence = execution.result.journal.size();
        solved.state_before = PracticalDriverState::SolveAndEstimate;
        solved.state_after = refinements >= config.work_limits.maximum_H_steps
            ? PracticalDriverState::WorkLimitReached
            : PracticalDriverState::RefineCoarse;
        solved.action = PracticalDriverAction::SolveStandardLod;
        solved.evaluation_candidate = candidate;
        solved.rebuilt_correctors = model.correctors().primal.size();
        solved.time_corrector_seconds =
            model.build_timings().correctors_ms / 1000.0;
        solved.time_solve_seconds = elapsed_seconds(solve_begin, solve_end);
        solved.detail =
            "rebuilt and solved standard LOD on the current uniform coarse mesh"
            "; full model build seconds="
            + numeric(elapsed_seconds(build_begin, build_end));
        fill_mesh_fields(solved);
        execution.result.journal.push_back(std::move(solved));

        if (refinements >= config.work_limits.maximum_H_steps) {
            execution.result.state = PracticalDriverState::WorkLimitReached;
            execution.result.stop_reason =
                "completed configured standard LOD refinement trajectory";
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
                "standard LOD uniform marking produced no coarse refinement");
        }
        ++refinements;
        PracticalIterationRecord refined_record;
        refined_record.sequence = execution.result.journal.size();
        refined_record.state_before = PracticalDriverState::RefineCoarse;
        refined_record.state_after = PracticalDriverState::SolveAndEstimate;
        refined_record.action = PracticalDriverAction::RefineUniformLod;
        refined_record.marked_H = marked.size();
        refined_record.time_mesh_seconds = elapsed_seconds(mesh_begin, mesh_end);
        refined_record.detail =
            "uniformly refined every standard LOD coarse element";
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

PaperRunStatus paper_status(const PracticalDriverState state,
                            const std::string &reason) {
    switch (state) {
    case PracticalDriverState::Converged:
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
    out << "schema_version,case,method,kappa,run_id,reference_epoch,iteration,action,stop_reason,"
           "N_H,N_ref,N_amb,ell,kappa_H_max,mu_H,rho_amb,eta_H,Theta_loc,U_prac,"
           "reference_energy_error,reference_L2_error,marked_H,rebuilt_correctors,"
           "ambient_refined_elements,time_mesh,time_corrector,time_certificate,time_solve,"
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
    for (const PracticalIterationRecord &record : result.journal) {
        out << practical_paper_schema_version << ',' << to_string(config.case_id) << ','
            << to_string(config.method_id) << ',' << numeric(config.wavenumber) << ','
            << run_id << ',' << record.reference_epoch << ',' << record.sequence << ','
            << practical_driver_action_name(record.action) << ','
            << csv_string(record.detail) << ',' << record.coarse_nodes << ','
            << record.reference_nodes << ','
            << (localization_method ? std::to_string(record.ambient_nodes) : "")
            << ',' << (ell_method ? std::to_string(record.ell) : "") << ','
            << numeric(record.kappa_H_max) << ",,"
            << (localization_method ? numeric(record.rho_ambient) : "") << ','
            << (eta_method ? numeric(record.eta_H) : "") << ','
            << (localization_method ? numeric(record.theta_loc) : "") << ','
            << (practical_bound_method ? numeric(record.U_practical) : "") << ','
            << optional_numeric(record.reference_energy_error) << ','
            << optional_numeric(record.reference_L2_error) << ','
            << record.marked_H << ',' << record.rebuilt_correctors << ','
            << record.ambient_refined_elements << ','
            << numeric(record.time_mesh_seconds) << ','
            << numeric(record.time_corrector_seconds) << ','
            << numeric(record.time_certificate_seconds) << ','
            << numeric(record.time_solve_seconds) << ','
            << numeric(record.time_estimator_seconds) << ','
            << numeric(record.time_total_cumulative_seconds) << ','
            << numeric(peak_mb) << '\n';
    }
}

void write_summary(
    const std::filesystem::path &path,
    const PracticalPaperConfig &config,
    const PracticalDriverResult &result) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "target,status,first_iteration,reference_energy_error,U_prac,N_H,ell,"
           "rebuilt_correctors_cumulative,ambient_refined_elements_cumulative,"
           "time_total_cumulative\n";
    const std::vector<double> targets(
        config.relative_energy_targets.begin(),
        config.relative_energy_targets.end());
    for (const PracticalTargetHit &target_hit :
         extract_practical_target_hits(result.journal, targets)) {
        out << numeric(target_hit.target) << ',';
        if (!target_hit.journal_index) {
            out << "not_reached,,,,,,,,\n";
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
                << hit->coarse_nodes << ','
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
    out << "iteration,action,ell,Theta_loc,rebuilt_correctors,time_corrector,time_certificate\n";
    for (const PracticalIterationRecord &record : result.journal) {
        if (record.state_before != PracticalDriverState::LocalizationCheck
            && record.action != PracticalDriverAction::SolveStandardLod) {
            continue;
        }
        out << record.sequence << ',' << practical_driver_action_name(record.action) << ','
            << record.ell << ','
            << (config.method_id == PracticalPaperMethod::Palod
                    ? numeric(record.theta_loc) : "") << ','
            << record.rebuilt_correctors << ','
            << numeric(record.time_corrector_seconds) << ','
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
    const double peak_mb) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "{\n  \"schema_version\":2,\n"
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
        << "  \"config\":" << canonical_json(config) << ",\n"
        << "  \"files\":{\"iterations\":\"iterations.csv\","
           "\"summary\":\"summary.csv\",\"ell_history\":\"ell_history.csv\","
           "\"final_mesh\":\"final_mesh.vtu\"},\n"
        << "  \"timing\":{\"method_seconds\":" << numeric(method_seconds)
        << ",\"evaluation_reference_seconds\":" << numeric(reference_seconds)
        << ",\"evaluation_reference_excluded_from_method_time\":true},\n"
        << "  \"hardware\":{\"hardware_threads\":"
        << std::thread::hardware_concurrency()
        << ",\"compiler\":" << json_string(__VERSION__)
        << ",\"peak_memory_mb\":" << numeric(peak_mb) << "}\n}\n";
}

} // namespace

int main(const int argc, char **argv) {
    try {
        const Arguments arguments = parse_arguments(argc, argv);
        const PracticalPaperConfig config =
            parse_practical_paper_config(read_text(arguments.config));
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
        const PaperCaseData data = make_paper_case(config.case_id, config.wavenumber);

        const auto reference_begin = std::chrono::steady_clock::now();
        ReferenceEpochHierarchy reference_hierarchy(
            data.initial_mesh, config.initial_coarse_level, config.reference_level);
        const HelmholtzOperators reference_operators = assemble_helmholtz_operators(
            reference_hierarchy.reference_mesh(), config.wavenumber, {}, {},
            config.boundary_beta);
        const ComplexVector reference_load = assemble_helmholtz_load(
            reference_hierarchy.reference_mesh(), data.source,
            config.quadrature, data.quadrature_context);
        const ComplexVector reference_solution =
            solve_helmholtz_fem(reference_operators, reference_load);
        const auto reference_end = std::chrono::steady_clock::now();

        const auto method_begin = std::chrono::steady_clock::now();
        PaperExecution execution;
        if (config.method_id == PracticalPaperMethod::Ufem) {
            execution = run_uniform_fem_trajectory(config, data);
        } else if (config.method_id == PracticalPaperMethod::Afem) {
            execution = run_adaptive_fem_trajectory(config, data);
        } else if (config.method_id == PracticalPaperMethod::Slod) {
            execution = run_standard_lod_trajectory(config, data);
        } else {
            PracticalDriverProblem problem;
            problem.initial_mesh = data.initial_mesh;
            problem.source = data.source;
            problem.quadrature = config.quadrature;
            problem.quadrature_context = data.quadrature_context;
            PracticalAdaptiveDriver driver(
                std::move(problem), make_practical_driver_config(config));
            execution.result = driver.run();
            execution.final_mesh = driver.hierarchy().coarse_mesh();
        }
        PracticalDriverResult &result = execution.result;
        const auto method_end = std::chrono::steady_clock::now();
        const auto evaluation_begin = std::chrono::steady_clock::now();
        for (PracticalIterationRecord &record : result.journal) {
            if (record.evaluation_candidate.size() == 0) continue;
            const auto errors = relative_reference_errors(
                reference_operators, reference_solution,
                record.evaluation_candidate);
            record.reference_energy_error = errors.first;
            record.reference_L2_error = errors.second;
            record.evaluation_candidate.resize(0);
        }
        const auto evaluation_end = std::chrono::steady_clock::now();
        const bool fixed_empirical_trajectory =
            config.method_id == PracticalPaperMethod::Ufem
            || config.method_id == PracticalPaperMethod::Slod
            || config.method_id == PracticalPaperMethod::Afem;
        if (fixed_empirical_trajectory
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
            elapsed_seconds(method_begin, method_end),
            elapsed_seconds(reference_begin, reference_end) +
                elapsed_seconds(evaluation_begin, evaluation_end),
            peak_mb);

        if (arguments.check) {
            const bool acceptable_state =
                result.state == PracticalDriverState::Converged
                || (fixed_empirical_trajectory
                    && result.state == PracticalDriverState::WorkLimitReached);
            if (!acceptable_state ||
                std::none_of(result.journal.begin(), result.journal.end(),
                             [](const PracticalIterationRecord &record) {
                                 return record.action == PracticalDriverAction::Complete
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
            if (std::none_of(result.journal.begin(), result.journal.end(),
                             [](const PracticalIterationRecord &record) {
                                 return record.reference_energy_error &&
                                        record.reference_L2_error;
                             })) {
                throw std::runtime_error(
                    "WP5 smoke did not compute post-run reference errors");
            }
            const auto has_action = [&](const PracticalDriverAction action) {
                return std::any_of(
                    result.journal.begin(), result.journal.end(),
                    [action](const PracticalIterationRecord &record) {
                        return record.action == action;
                    });
            };
            if (config.method_id == PracticalPaperMethod::Palod
                && (!has_action(PracticalDriverAction::AcceptLocalization)
                    || !has_action(PracticalDriverAction::Complete))) {
                throw std::runtime_error(
                    "PALOD smoke did not execute localization and completion");
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
                const int expected_ell =
                    standard_lod_prior_ell(config.wavenumber);
                const bool wrong_ell = std::any_of(
                    result.journal.begin(), result.journal.end(),
                    [expected_ell](const PracticalIterationRecord &record) {
                        return record.action
                                   == PracticalDriverAction::SolveStandardLod
                            && record.ell != expected_ell;
                    });
                const bool paid_adaptive_cost = std::any_of(
                    result.journal.begin(), result.journal.end(),
                    [](const PracticalIterationRecord &record) {
                        return record.time_certificate_seconds != 0.0
                            || record.time_estimator_seconds != 0.0
                            || record.ambient_refined_elements != 0;
                    });
                if (!has_action(PracticalDriverAction::SolveStandardLod)
                    || !has_action(PracticalDriverAction::RefineUniformLod)
                    || has_action(PracticalDriverAction::AcceptLocalization)
                    || has_action(PracticalDriverAction::AcceptFixedEll)
                    || has_action(PracticalDriverAction::FormCoarseMarking)
                    || wrong_ell || paid_adaptive_cost) {
                    throw std::runtime_error(
                        "SLOD smoke did not remain on the frozen-prior uniform LOD path");
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
                  << "claim=implementation-study\n"
                  << "output=" << run_directory.string() << '\n';
        return result.state == PracticalDriverState::Failed ? 2 : 0;
    } catch (const std::exception &error) {
        std::cerr << "bench_helmholtz_adaptive_paper failed: "
                  << error.what() << '\n';
        return 1;
    }
}

#include "helmholtz/adaptive/practical_driver.h"
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
    for (const PracticalIterationRecord &record : result.journal) {
        out << practical_paper_schema_version << ',' << to_string(config.case_id) << ','
            << to_string(config.method_id) << ',' << numeric(config.wavenumber) << ','
            << run_id << ',' << record.reference_epoch << ',' << record.sequence << ','
            << practical_driver_action_name(record.action) << ','
            << csv_string(record.detail) << ',' << record.coarse_nodes << ','
            << record.reference_nodes << ',' << record.ambient_nodes << ','
            << record.ell << ',' << numeric(record.kappa_H_max) << ",,"
            << numeric(record.rho_ambient) << ',' << numeric(record.eta_H) << ','
            << numeric(record.theta_loc) << ',' << numeric(record.U_practical) << ','
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
            out << "reached," << hit->sequence << ','
                << numeric(*hit->reference_energy_error) << ','
                << numeric(hit->U_practical) << ',' << hit->coarse_nodes << ','
                << hit->ell << ',' << rebuilt << ',' << ambient_refined << ','
                << numeric(hit->time_total_cumulative_seconds)
                << '\n';
        }
    }
}

void write_ell_history(
    const std::filesystem::path &path,
    const PracticalDriverResult &result) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "iteration,action,ell,Theta_loc,rebuilt_correctors,time_corrector,time_certificate\n";
    for (const PracticalIterationRecord &record : result.journal) {
        if (record.state_before != PracticalDriverState::LocalizationCheck) continue;
        out << record.sequence << ',' << practical_driver_action_name(record.action) << ','
            << record.ell << ',' << numeric(record.theta_loc) << ','
            << record.rebuilt_correctors << ','
            << numeric(record.time_corrector_seconds) << ','
            << numeric(record.time_certificate_seconds) << '\n';
    }
}

void write_final_mesh(
    const std::filesystem::path &path,
    const PracticalAdaptiveDriver &driver,
    const PracticalDriverResult &result) {
    const TriMesh &mesh = driver.hierarchy().coarse_mesh();
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
        if (config.method_id != PracticalPaperMethod::Palod) {
            throw std::invalid_argument(
                "WP5 runner currently executes PALOD only; baseline methods remain pending");
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

        PracticalDriverProblem problem;
        problem.initial_mesh = data.initial_mesh;
        problem.source = data.source;
        problem.quadrature = config.quadrature;
        problem.quadrature_context = data.quadrature_context;
        const auto method_begin = std::chrono::steady_clock::now();
        PracticalAdaptiveDriver driver(
            std::move(problem), make_practical_driver_config(config));
        PracticalDriverResult result = driver.run();
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

        const std::string run_id = make_run_id(config);
        const std::filesystem::path run_directory =
            arguments.output_directory / run_id;
        std::filesystem::create_directories(run_directory);
        const double peak_mb = peak_memory_mb();
        write_iterations(run_directory / "iterations.csv", config, run_id, result, peak_mb);
        write_summary(run_directory / "summary.csv", config, result);
        write_ell_history(run_directory / "ell_history.csv", result);
        write_final_mesh(run_directory / "final_mesh.vtu", driver, result);
        write_run_json(
            run_directory / "run.json", config, run_id, result,
            elapsed_seconds(method_begin, method_end),
            elapsed_seconds(reference_begin, reference_end) +
                elapsed_seconds(evaluation_begin, evaluation_end),
            peak_mb);

        if (arguments.check) {
            if (result.state != PracticalDriverState::Converged ||
                std::none_of(result.journal.begin(), result.journal.end(),
                             [](const PracticalIterationRecord &record) {
                                 return record.action == PracticalDriverAction::Complete;
                             })) {
                throw std::runtime_error("WP5 smoke did not complete the real PALOD chain");
            }
            if (std::none_of(result.journal.begin(), result.journal.end(),
                             [](const PracticalIterationRecord &record) {
                                 return record.reference_energy_error &&
                                        record.reference_L2_error;
                             })) {
                throw std::runtime_error(
                    "WP5 smoke did not compute post-run reference errors");
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

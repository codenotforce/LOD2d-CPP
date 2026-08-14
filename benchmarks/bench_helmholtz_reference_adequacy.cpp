#include "helmholtz/adaptive/hierarchy.h"
#include "helmholtz/benchmarks/paper_cases.h"
#include "helmholtz/experiments/paper_config.h"
#include "helmholtz/experiments/reference_solution_cache.h"
#include "helmholtz/operators.h"
#include "mesh/refine.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using namespace lod2d;
using namespace lod2d::helmholtz;
using namespace lod2d::helmholtz::adaptive;
using namespace lod2d::helmholtz::benchmarks;
using namespace lod2d::helmholtz::experiments;

namespace {

struct Arguments {
    std::filesystem::path config;
    std::filesystem::path source_iterations;
    std::filesystem::path output_directory;
    std::filesystem::path reference_cache_directory;
    std::filesystem::path manuscript_baseline;
    bool check = false;
};

struct SourceTrajectoryPoint {
    int schema_version = 0;
    PaperCase case_id = PaperCase::R1;
    PracticalPaperMethod method_id = PracticalPaperMethod::Palod;
    double wavenumber = 0.0;
    std::string run_id;
    std::string content_fingerprint;
    double terminal_reference_error = 0.0;
};

std::string read_text(const std::filesystem::path &path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("cannot read " + path.string());
    return {std::istreambuf_iterator<char>(input),
            std::istreambuf_iterator<char>()};
}

std::string first_token(const std::string &text) {
    std::istringstream input(text);
    std::string token;
    input >> token;
    return token;
}

std::string json_string(const std::string &value) {
    std::string result = "\"";
    for (const char character : value) {
        if (character == '\\' || character == '"') result += '\\';
        result += character;
    }
    result += '"';
    return result;
}

template <class Value>
std::string numeric(const Value value) {
    std::ostringstream out;
    out << std::setprecision(17) << value;
    return out.str();
}

double elapsed_seconds(
    const std::chrono::steady_clock::time_point begin,
    const std::chrono::steady_clock::time_point end) {
    return std::chrono::duration<double>(end - begin).count();
}

Arguments parse_arguments(const int argc, char **argv) {
    Arguments result;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        const auto value = [&](const std::string &prefix) {
            return std::filesystem::path(argument.substr(prefix.size()));
        };
        if (argument.rfind("--config=", 0) == 0)
            result.config = value("--config=");
        else if (argument.rfind("--source-iterations=", 0) == 0)
            result.source_iterations = value("--source-iterations=");
        else if (argument.rfind("--output-dir=", 0) == 0)
            result.output_directory = value("--output-dir=");
        else if (argument.rfind("--reference-cache-dir=", 0) == 0)
            result.reference_cache_directory = value("--reference-cache-dir=");
        else if (argument.rfind("--manuscript-baseline=", 0) == 0)
            result.manuscript_baseline = value("--manuscript-baseline=");
        else if (argument == "--check")
            result.check = true;
        else
            throw std::invalid_argument("unknown reference-audit argument: " + argument);
    }
    if (result.config.empty() || result.source_iterations.empty()
        || result.output_directory.empty() || result.manuscript_baseline.empty()) {
        throw std::invalid_argument(
            "reference audit requires config, source iterations, output dir, and manuscript baseline");
    }
    return result;
}

std::vector<std::string> parse_csv_line(const std::string &line) {
    std::vector<std::string> fields;
    std::string field;
    bool quoted = false;
    for (std::size_t index = 0; index < line.size(); ++index) {
        const char character = line[index];
        if (character == '"') {
            if (quoted && index + 1 < line.size() && line[index + 1] == '"') {
                field += '"';
                ++index;
            } else {
                quoted = !quoted;
            }
        } else if (character == ',' && !quoted) {
            fields.push_back(field);
            field.clear();
        } else if (character != '\r') {
            field += character;
        }
    }
    if (quoted) throw std::runtime_error("unterminated quoted CSV field");
    fields.push_back(field);
    return fields;
}

SourceTrajectoryPoint read_source_trajectory(
    const std::filesystem::path &path) {
    const std::string content = read_text(path);
    std::uint64_t hash = 14695981039346656037ULL;
    for (const unsigned char byte : content) {
        hash ^= byte;
        hash *= 1099511628211ULL;
    }
    std::ifstream input(path);
    std::string line;
    if (!std::getline(input, line))
        throw std::runtime_error("source iterations CSV is empty");
    const std::vector<std::string> header = parse_csv_line(line);
    std::unordered_map<std::string, std::size_t> column;
    for (std::size_t index = 0; index < header.size(); ++index)
        column.emplace(header[index], index);
    for (const char *required : {
             "schema_version", "case", "method", "kappa", "run_id",
             "reference_energy_error"}) {
        if (!column.contains(required))
            throw std::runtime_error(
                "source iterations CSV is missing column " + std::string(required));
    }

    SourceTrajectoryPoint result;
    {
        std::ostringstream encoded;
        encoded << "fnv1a64:" << std::hex << std::setfill('0')
                << std::setw(16) << hash;
        result.content_fingerprint = encoded.str();
    }
    bool initialized = false;
    bool found_error = false;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> fields = parse_csv_line(line);
        if (fields.size() != header.size())
            throw std::runtime_error("source iterations CSV row width changed");
        if (!initialized) {
            result.schema_version = std::stoi(fields[column.at("schema_version")]);
            result.case_id = parse_paper_case(fields[column.at("case")]);
            result.method_id = parse_practical_paper_method(fields[column.at("method")]);
            result.wavenumber = std::stod(fields[column.at("kappa")]);
            result.run_id = fields[column.at("run_id")];
            initialized = true;
        } else if (std::stoi(fields[column.at("schema_version")])
                       != result.schema_version
                   || parse_paper_case(fields[column.at("case")])
                       != result.case_id
                   || parse_practical_paper_method(fields[column.at("method")])
                       != result.method_id
                   || std::stod(fields[column.at("kappa")])
                       != result.wavenumber
                   || fields[column.at("run_id")] != result.run_id) {
            throw std::runtime_error(
                "source iterations CSV changes run identity between rows");
        }
        const std::string &error = fields[column.at("reference_energy_error")];
        if (!error.empty()) {
            result.terminal_reference_error = std::stod(error);
            found_error = true;
        }
    }
    if (!initialized || !found_error
        || !std::isfinite(result.terminal_reference_error)
        || !(result.terminal_reference_error > 0.0)) {
        throw std::runtime_error(
            "source iterations CSV has no positive terminal reference error");
    }
    return result;
}

double relative_energy_difference(
    const HelmholtzOperators &fine_operators,
    const ComplexVector &fine_solution,
    const ComplexVector &prolongated_reference) {
    if (fine_solution.size() != prolongated_reference.size())
        throw std::invalid_argument("reference audit prolongation size mismatch");
    const Eigen::SparseMatrix<double> energy = fine_operators.stiffness
        + fine_operators.wavenumber * fine_operators.wavenumber
            * fine_operators.mass;
    const auto square = [&](const ComplexVector &value) {
        return std::max(
            0.0, std::real(value.dot(energy.cast<Complex>() * value)));
    };
    const double denominator = square(fine_solution);
    if (!(denominator > 0.0))
        throw std::runtime_error("fine reference solution has zero energy norm");
    return std::sqrt(square(fine_solution - prolongated_reference) / denominator);
}

struct CachedReference {
    ComplexVector solution;
    std::string key;
    bool cache_hit = false;
};

CachedReference load_or_solve_reference(
    const std::filesystem::path &cache_directory,
    const TriMesh &mesh,
    const HelmholtzOperators &operators,
    const ComplexVector &load,
    const std::string &identity) {
    CachedReference result;
    result.key = reference_solution_cache_key(mesh, operators, load, identity);
    ReferenceSolutionCacheLookup lookup = load_reference_solution_cache(
        cache_directory, result.key, load.size());
    result.cache_hit = lookup.solution.has_value();
    if (lookup.solution) {
        result.solution = std::move(*lookup.solution);
    } else {
        result.solution = solve_helmholtz_fem(operators, load);
        store_reference_solution_cache(
            cache_directory, result.key, result.solution);
    }
    return result;
}

} // namespace

int main(const int argc, char **argv) {
    try {
        const Arguments arguments = parse_arguments(argc, argv);
        const PracticalPaperConfig config =
            parse_practical_paper_config(read_text(arguments.config));
        if (!config.reference_adequacy.enabled)
            throw std::invalid_argument("reference adequacy audit is disabled in config");
        if (first_token(read_text(arguments.manuscript_baseline))
            != config.manuscript_sha256) {
            throw std::invalid_argument("reference audit manuscript hash mismatch");
        }
        const SourceTrajectoryPoint source =
            read_source_trajectory(arguments.source_iterations);
        if (source.case_id != config.case_id
            || source.method_id != config.method_id
            || source.wavenumber != config.wavenumber) {
            throw std::invalid_argument(
                "reference audit source trajectory does not match config");
        }

        const PaperCaseData data = make_paper_case(
            config.case_id, config.wavenumber,
            config.singular_oscillatory_fraction,
            config.singular_cutoff_outer_radius,
            config.singular_quintic_cutoff,
            config.smooth_wave_amplitude);
        ReferenceEpochHierarchy hierarchy(
            data.initial_mesh, config.initial_coarse_level,
            config.reference_level, config.reference_epoch);
        const TriMesh &base_mesh = hierarchy.reference_mesh();
        const HelmholtzOperators base_operators = assemble_helmholtz_operators(
            base_mesh, config.wavenumber, {}, {}, config.boundary_beta);
        const ComplexVector base_load = assemble_helmholtz_load(
            base_mesh, data.source, config.quadrature,
            data.quadrature_context);
        const std::string identity = "reference-adequacy-sparse-lu-v1/"
            + config.git_commit + "/" + config.build_hash;

        const auto begin = std::chrono::steady_clock::now();
        const CachedReference base = load_or_solve_reference(
            arguments.reference_cache_directory, base_mesh, base_operators,
            base_load, identity + "/base");
        const RefineOutput refined = refine_mesh_nvb(
            base_mesh, config.reference_adequacy.refinement_levels);
        ReferenceEpochHierarchy expected_fine_hierarchy(
            data.initial_mesh, config.initial_coarse_level,
            config.reference_level + config.reference_adequacy.refinement_levels,
            config.reference_epoch + 1);
        const TriMesh &expected_fine = expected_fine_hierarchy.reference_mesh();
        if (refined.mesh.nodes.size() != expected_fine.nodes.size()
            || refined.mesh.elems != expected_fine.elems) {
            throw std::runtime_error(
                "reference audit refinement does not match the next uniform NVB level");
        }
        for (std::size_t node = 0; node < refined.mesh.nodes.size(); ++node) {
            if ((refined.mesh.nodes[node] - expected_fine.nodes[node]).norm()
                > 1e-14) {
                throw std::runtime_error(
                    "reference audit node coordinates drifted from the next epoch");
            }
        }
        const HelmholtzOperators fine_operators = assemble_helmholtz_operators(
            refined.mesh, config.wavenumber, {}, {}, config.boundary_beta);
        const ComplexVector fine_load = assemble_helmholtz_load(
            refined.mesh, data.source, config.quadrature,
            data.quadrature_context);
        const CachedReference fine = load_or_solve_reference(
            arguments.reference_cache_directory, refined.mesh, fine_operators,
            fine_load, identity + "/fine");
        const ComplexVector prolongated =
            refined.P_node.cast<Complex>() * base.solution;
        const double relative_difference = relative_energy_difference(
            fine_operators, fine.solution, prolongated);
        const double terminal_fraction = relative_difference
            / source.terminal_reference_error;
        const bool refresh_recommended = terminal_fraction
            > config.reference_adequacy.maximum_terminal_error_fraction;
        const double seconds = elapsed_seconds(
            begin, std::chrono::steady_clock::now());

        std::filesystem::create_directories(arguments.output_directory);
        const std::filesystem::path output =
            arguments.output_directory / "reference_adequacy.json";
        std::ofstream out(output);
        if (!out) throw std::runtime_error("cannot write " + output.string());
        out << "{\n  \"schema_version\":1,\n"
            << "  \"status\":"
            << json_string(refresh_recommended
                    ? "reference_refresh_recommended" : "reference_adequate")
            << ",\n  \"source_run_id\":" << json_string(source.run_id)
            << ",\n  \"source_schema_version\":" << source.schema_version
            << ",\n  \"source_iterations_fingerprint\":"
            << json_string(source.content_fingerprint)
            << ",\n  \"config_hash\":"
            << json_string(canonical_config_hash(config))
            << ",\n  \"case\":" << json_string(std::string(to_string(config.case_id)))
            << ",\n  \"method\":" << json_string(std::string(to_string(config.method_id)))
            << ",\n  \"wavenumber\":" << numeric(config.wavenumber)
            << ",\n  \"base_reference_level\":" << config.reference_level
            << ",\n  \"audit_reference_level\":"
            << config.reference_level + config.reference_adequacy.refinement_levels
            << ",\n  \"base_reference_nodes\":" << base_mesh.nodes.size()
            << ",\n  \"audit_reference_nodes\":" << refined.mesh.nodes.size()
            << ",\n  \"terminal_reference_error\":"
            << numeric(source.terminal_reference_error)
            << ",\n  \"relative_reference_difference\":"
            << numeric(relative_difference)
            << ",\n  \"terminal_error_fraction\":" << numeric(terminal_fraction)
            << ",\n  \"maximum_terminal_error_fraction\":"
            << numeric(config.reference_adequacy.maximum_terminal_error_fraction)
            << ",\n  \"reference_refresh_recommended\":"
            << (refresh_recommended ? "true" : "false")
            << ",\n  \"base_cache_hit\":" << (base.cache_hit ? "true" : "false")
            << ",\n  \"audit_cache_hit\":" << (fine.cache_hit ? "true" : "false")
            << ",\n  \"base_cache_key\":" << json_string(base.key)
            << ",\n  \"audit_cache_key\":" << json_string(fine.key)
            << ",\n  \"evaluation_seconds\":" << numeric(seconds) << "\n}\n";
        out.close();

        if (arguments.check
            && (!(relative_difference >= 0.0)
                || !std::isfinite(relative_difference)
                || !std::filesystem::is_regular_file(output))) {
            throw std::runtime_error("reference adequacy smoke failed");
        }
        std::cout << "source_run_id=" << source.run_id << '\n'
                  << "relative_reference_difference=" << relative_difference << '\n'
                  << "terminal_error_fraction=" << terminal_fraction << '\n'
                  << "reference_refresh_recommended="
                  << (refresh_recommended ? "true" : "false") << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "bench_helmholtz_reference_adequacy failed: "
                  << error.what() << '\n';
        return 1;
    }
}

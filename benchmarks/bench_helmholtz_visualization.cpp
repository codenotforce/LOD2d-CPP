#include "helmholtz/manufactured.h"
#include "helmholtz/model.h"
#include "helmholtz/operators.h"
#include "io/vtk_writer.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <span>
#include <stdexcept>
#include <string>

using namespace lod2d;
using namespace lod2d::helmholtz;
using namespace lod2d::io;

namespace {

struct Options {
    double wavenumber = 4.0;
    int coarse_level = 4;
    int fine_level = 8;
    int oversampling = 3;
    std::filesystem::path output_dir =
        "results/visualization/helmholtz_manufactured_k4_H4_h8_ell3";
};

int parse_int(const std::string &text, const char *name) {
    std::size_t consumed = 0;
    const int value = std::stoi(text, &consumed);
    if (consumed != text.size())
        throw std::invalid_argument(std::string("invalid ") + name + ": " + text);
    return value;
}

double parse_double(const std::string &text, const char *name) {
    std::size_t consumed = 0;
    const double value = std::stod(text, &consumed);
    if (consumed != text.size() || !std::isfinite(value))
        throw std::invalid_argument(std::string("invalid ") + name + ": " + text);
    return value;
}

Options parse_options(int argc, char **argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        auto value_after = [&](const std::string &prefix) {
            return argument.substr(prefix.size());
        };
        if (argument.rfind("--k=", 0) == 0) {
            options.wavenumber = parse_double(value_after("--k="), "wavenumber");
        } else if (argument.rfind("--H=", 0) == 0) {
            options.coarse_level = parse_int(value_after("--H="), "coarse level");
        } else if (argument.rfind("--h=", 0) == 0) {
            options.fine_level = parse_int(value_after("--h="), "fine level");
        } else if (argument.rfind("--ell=", 0) == 0) {
            options.oversampling = parse_int(value_after("--ell="), "oversampling");
        } else if (argument.rfind("--output-dir=", 0) == 0) {
            options.output_dir = value_after("--output-dir=");
        } else if (argument == "--help") {
            std::cout
                << "Usage: bench_helmholtz_visualization [--k=4] [--H=4] "
                   "[--h=8] [--ell=3] [--output-dir=PATH]\n"
                   "This explicit benchmark performs the extra fine reference solve "
                   "and writes VTU/JSON files; normal solver benchmarks do neither.\n";
            std::exit(0);
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (!(options.wavenumber > 0.0))
        throw std::invalid_argument("k must be positive");
    if (options.coarse_level < 0 || options.fine_level < options.coarse_level)
        throw std::invalid_argument("levels must satisfy 0 <= H <= h");
    if (options.oversampling < 0)
        throw std::invalid_argument("ell must be nonnegative");
    return options;
}

std::string json_escape(const std::string &text) {
    std::string result;
    result.reserve(text.size());
    for (char ch : text) {
        switch (ch) {
        case '\\': result += "\\\\"; break;
        case '"': result += "\\\""; break;
        case '\n': result += "\\n"; break;
        case '\r': result += "\\r"; break;
        case '\t': result += "\\t"; break;
        default: result += ch; break;
        }
    }
    return result;
}

void write_manifest(
    const std::filesystem::path &path,
    const Options &options,
    const HelmholtzLodModel &model,
    const HelmholtzLodSolution &lod,
    const HelmholtzError &reference_error,
    const HelmholtzError &lod_error) {
    std::ofstream output(path);
    if (!output)
        throw std::runtime_error("cannot open manifest: " + path.string());
    const auto &problem = model.problem();
    const auto &diagnostics = model.correctors().diagnostics;
    output << std::setprecision(17)
           << "{\n"
           << "  \"schema\": \"lod2d.helmholtz.visualization.v1\",\n"
           << "  \"problem\": \"manufactured_polynomial_plane_wave\",\n"
           << "  \"norm\": \"absolute weighted energy norm: "
              "sqrt(|e|_H1_seminorm^2 + k^2 ||e||_L2^2)\",\n"
           << "  \"parameters\": {\"k\": " << options.wavenumber
           << ", \"H_level\": " << options.coarse_level
           << ", \"h_level\": " << options.fine_level
           << ", \"ell\": " << options.oversampling << "},\n"
           << "  \"mesh\": {\"H_max\": " << max_element_diameter(problem.coarse)
           << ", \"h_max\": " << max_element_diameter(problem.fine)
           << ", \"coarse_nodes\": " << problem.coarse.nodes.size()
           << ", \"coarse_elements\": " << problem.coarse.elems.size()
           << ", \"fine_nodes\": " << problem.fine.nodes.size()
           << ", \"fine_elements\": " << problem.fine.elems.size() << "},\n"
           << "  \"errors\": {\"reference_energy_abs\": " << reference_error.energy
           << ", \"reference_l2_abs\": " << reference_error.l2
           << ", \"lod_energy_abs\": " << lod_error.energy
           << ", \"lod_l2_abs\": " << lod_error.l2 << "},\n"
           << "  \"residuals\": {\"petrov\": " << lod.petrov_residual
           << ", \"corrector\": " << diagnostics.max_primal_residual
           << ", \"constraint\": " << diagnostics.max_constraint_residual << "},\n"
           << "  \"fields\": {\"complex_encoding\": \"paired_real_imag\", "
              "\"point\": [\"u_exact\", \"u_reference\", \"u_lod\", "
              "\"u_coarse\", \"u_fine_scale\", \"error_lod\"], "
              "\"cell\": [\"diffusion\", \"refractive_index\"]},\n"
           << "  \"files\": {\"coarse_mesh\": \""
           << json_escape("coarse_mesh.vtu")
           << "\", \"fine_solution\": \"" << json_escape("fine_solution.vtu")
           << "\"}\n"
           << "}\n";
    if (!output)
        throw std::runtime_error("failed while writing manifest: " + path.string());
}

std::span<const Complex> view(const ComplexVector &values) {
    return {values.data(), static_cast<std::size_t>(values.size())};
}

} // namespace

int main(int argc, char **argv) {
    try {
        const Options options = parse_options(argc, argv);
        const HelmholtzManufacturedSolution manufactured =
            make_polynomial_plane_wave_solution(options.wavenumber);

        HelmholtzProblemConfig config;
        config.H = options.coarse_level;
        config.h = options.fine_level;
        config.ell = options.oversampling;
        config.wavenumber = options.wavenumber;
        config.mode = HelmholtzPetrovMode::TwoSided;
        HelmholtzLodModel model = HelmholtzLodModel::build(config);

        const ComplexVector load =
            assemble_helmholtz_load(model.problem().fine, manufactured.source);
        const HelmholtzLodSolution lod = model.solve_load(load);
        const ComplexVector reference = model.solve_fine_reference(load);
        const ComplexVector coarse =
            model.problem().coarse_to_fine.cast<Complex>() * lod.coarse_coefficients;
        const ComplexVector fine_scale = lod.fine_values - coarse;
        const ComplexVector lod_error_values = lod.fine_values - reference;

        ComplexVector exact(static_cast<Eigen::Index>(model.problem().fine.nodes.size()));
        for (Eigen::Index node = 0; node < exact.size(); ++node)
            exact(node) = manufactured.value(
                model.problem().fine.nodes[static_cast<std::size_t>(node)]);

        const HelmholtzError reference_error = compute_helmholtz_error(
            model.problem().fine,
            reference,
            options.wavenumber,
            manufactured.value,
            manufactured.gradient);
        const HelmholtzError lod_error = compute_helmholtz_error(
            model.problem().fine,
            lod.fine_values,
            options.wavenumber,
            manufactured.value,
            manufactured.gradient);

        std::filesystem::create_directories(options.output_dir);
        write_vtu(options.output_dir / "coarse_mesh.vtu", model.problem().coarse);

        const std::array point_complex{
            VtkComplexFieldView{"u_exact", view(exact)},
            VtkComplexFieldView{"u_reference", view(reference)},
            VtkComplexFieldView{"u_lod", view(lod.fine_values)},
            VtkComplexFieldView{"u_coarse", view(coarse)},
            VtkComplexFieldView{"u_fine_scale", view(fine_scale)},
            VtkComplexFieldView{"error_lod", view(lod_error_values)}};
        const std::array cell_double{
            VtkDoubleFieldView{"diffusion", model.operators().diffusion},
            VtkDoubleFieldView{
                "refractive_index", model.operators().refractive_index}};
        VtuDataView fields;
        fields.point_complex = point_complex;
        fields.cell_double = cell_double;
        write_vtu(
            options.output_dir / "fine_solution.vtu",
            model.problem().fine,
            fields);
        write_manifest(
            options.output_dir / "run.json",
            options,
            model,
            lod,
            reference_error,
            lod_error);

        std::cout << "wrote " << options.output_dir.string() << '\n'
                  << "fine nodes=" << model.problem().fine.nodes.size()
                  << " elements=" << model.problem().fine.elems.size()
                  << " lod_energy_abs=" << lod_error.energy << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "bench_helmholtz_visualization failed: "
                  << error.what() << '\n';
        return 1;
    }
}

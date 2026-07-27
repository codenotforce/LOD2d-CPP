#include "helmholtz/manufactured.h"
#include "helmholtz/model.h"
#include "helmholtz/operators.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace lod2d;
using namespace lod2d::helmholtz;

namespace {

namespace fs = std::filesystem;

struct Options {
    double wavenumber = 16.0;
    int fine_level = 15;
    std::vector<int> coarse_levels{7, 8, 9, 10, 11, 12};
    int ell = 4;
    int threads = 0;
    std::string solver = "saddle";
    int symbolic_cache_slots = 1;
    bool reuse_identical_factorization = false;
    HelmholtzPetrovMode mode = HelmholtzPetrovMode::TwoSided;
    bool solve_fine_reference = false;
    bool export_fields = false;
    bool csv = false;
    bool csv_header = false;
    bool check = false;
    fs::path export_dir;
    fs::path summary_out;
};

struct Row {
    double k = 0.0;
    int H_level = 0;
    int h_level = 0;
    int ell = 0;
    std::string solver;
    std::string mode;
    int coarse_nodes = 0;
    int coarse_elements = 0;
    int fine_nodes = 0;
    int fine_elements = 0;
    double H_max = 0.0;
    double H_min = 0.0;
    double h_max = 0.0;
    double h_min = 0.0;
    double kH = 0.0;
    double kh = 0.0;
    double h_over_H = 0.0;
    HelmholtzError p1_exact;
    HelmholtzError lod_exact;
    HelmholtzError fine_exact{
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()};
    double p1_energy_rate = std::numeric_limits<double>::quiet_NaN();
    double lod_energy_rate = std::numeric_limits<double>::quiet_NaN();
    double p1_l2_rate = std::numeric_limits<double>::quiet_NaN();
    double lod_l2_rate = std::numeric_limits<double>::quiet_NaN();
    double p1_residual = 0.0;
    double petrov_residual = 0.0;
    double corrector_residual = 0.0;
    double constraint_residual = 0.0;
    double nesting_coordinate_residual = 0.0;
    int corrector_threads = 1;
    int symbolic_analyses = 0;
    int symbolic_reuses = 0;
    int factorization_reuses = 0;
    double schur_residual = 0.0;
    double schur_rcond = 1.0;
    double mesh_ms = 0.0;
    double operators_ms = 0.0;
    double correctors_ms = 0.0;
    double basis_ms = 0.0;
    double build_ms = 0.0;
    double load_ms = 0.0;
    double lod_solve_ms = 0.0;
    double p1_solve_ms = 0.0;
    double fine_solve_ms = 0.0;
    double error_ms = 0.0;
    double export_ms = 0.0;
    double total_ms = 0.0;
    std::string coarse_nodes_file;
    std::string coarse_elements_file;
    std::string fine_nodes_file;
    std::string fine_elements_file;
    std::string fields_file;
};

double elapsed_ms(const std::chrono::steady_clock::time_point &start) {
    return std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - start).count();
}

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

std::vector<int> parse_levels(const std::string &text) {
    std::vector<int> result;
    std::stringstream stream(text);
    std::string token;
    while (std::getline(stream, token, ',')) {
        if (token.empty()) throw std::invalid_argument("empty H level");
        result.push_back(parse_int(token, "H level"));
    }
    if (result.empty()) throw std::invalid_argument("H levels must not be empty");
    return result;
}

Options parse_options(int argc, char **argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        auto value_after = [&](const std::string &prefix) {
            return argument.substr(prefix.size());
        };
        if (argument.rfind("--k=", 0) == 0)
            options.wavenumber = parse_double(value_after("--k="), "k");
        else if (argument.rfind("--h=", 0) == 0)
            options.fine_level = parse_int(value_after("--h="), "h level");
        else if (argument.rfind("--H-levels=", 0) == 0)
            options.coarse_levels = parse_levels(value_after("--H-levels="));
        else if (argument.rfind("--ell=", 0) == 0)
            options.ell = parse_int(value_after("--ell="), "ell");
        else if (argument.rfind("--threads=", 0) == 0)
            options.threads = parse_int(value_after("--threads="), "threads");
        else if (argument.rfind("--solver=", 0) == 0)
            options.solver = value_after("--solver=");
        else if (argument.rfind("--symbolic-cache-slots=", 0) == 0)
            options.symbolic_cache_slots =
                parse_int(value_after("--symbolic-cache-slots="), "symbolic cache slots");
        else if (argument.rfind("--factorization-reuse=", 0) == 0) {
            const std::string value = value_after("--factorization-reuse=");
            if (value == "identical") options.reuse_identical_factorization = true;
            else if (value == "none") options.reuse_identical_factorization = false;
            else throw std::invalid_argument(
                "factorization reuse must be none or identical");
        } else if (argument.rfind("--mode=", 0) == 0) {
            const std::string value = value_after("--mode=");
            if (value == "two-sided") options.mode = HelmholtzPetrovMode::TwoSided;
            else if (value == "test-only")
                options.mode = HelmholtzPetrovMode::CorrectedTestOnly;
            else throw std::invalid_argument("mode must be two-sided or test-only");
        } else if (argument.rfind("--export-dir=", 0) == 0)
            options.export_dir = value_after("--export-dir=");
        else if (argument.rfind("--summary-out=", 0) == 0)
            options.summary_out = value_after("--summary-out=");
        else if (argument == "--export-fields")
            options.export_fields = true;
        else if (argument == "--fine-reference")
            options.solve_fine_reference = true;
        else if (argument == "--format=csv")
            options.csv = true;
        else if (argument == "--csv-header")
            options.csv_header = true;
        else if (argument == "--check")
            options.check = true;
        else if (argument == "--help") {
            std::cout
                << "Usage: bench_helmholtz_H_convergence [options]\n"
                << "  --k=16 --h=15 --H-levels=7,8,9,10,11,12 --ell=4\n"
                << "  --solver=saddle|schur --threads=N\n"
                << "  --symbolic-cache-slots=N --factorization-reuse=none|identical\n"
                << "  --mode=two-sided|test-only --fine-reference\n"
                << "  --export-dir=PATH [--export-fields] --summary-out=PATH\n"
                << "  --format=csv | --csv-header | --check\n\n"
                << "H and h are global NVB sweep counts. P1 and LOD errors are both\n"
                << "absolute errors integrated on the same fixed fine mesh. Mesh CSVs\n"
                << "are always exported when --export-dir is set; nodal fields require\n"
                << "--export-fields.\n";
            std::exit(0);
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (!(options.wavenumber > 0.0)) throw std::invalid_argument("k must be positive");
    if (options.fine_level < 0) throw std::invalid_argument("h must be nonnegative");
    if (options.ell < 0) throw std::invalid_argument("ell must be nonnegative");
    if (options.threads < 0) throw std::invalid_argument("threads must be nonnegative");
    if (options.symbolic_cache_slots <= 0)
        throw std::invalid_argument("symbolic cache slots must be positive");
    if (options.solver != "saddle" && options.solver != "schur")
        throw std::invalid_argument("solver must be saddle or schur");
    std::sort(options.coarse_levels.begin(), options.coarse_levels.end());
    options.coarse_levels.erase(
        std::unique(options.coarse_levels.begin(), options.coarse_levels.end()),
        options.coarse_levels.end());
    for (int level : options.coarse_levels) {
        if (level < 0 || level >= options.fine_level)
            throw std::invalid_argument("each H level must satisfy 0 <= H < h");
    }
    return options;
}

double triangle_diameter(const TriMesh &mesh, const Triangle &triangle) {
    return std::max({
        (mesh.nodes[triangle[0]] - mesh.nodes[triangle[1]]).norm(),
        (mesh.nodes[triangle[1]] - mesh.nodes[triangle[2]]).norm(),
        (mesh.nodes[triangle[2]] - mesh.nodes[triangle[0]]).norm()});
}

double min_element_diameter(const TriMesh &mesh) {
    double result = std::numeric_limits<double>::infinity();
    for (const Triangle &triangle : mesh.elems)
        result = std::min(result, triangle_diameter(mesh, triangle));
    return result;
}

double triangle_area(const TriMesh &mesh, const Triangle &triangle) {
    const Point2 a = mesh.nodes[triangle[1]] - mesh.nodes[triangle[0]];
    const Point2 b = mesh.nodes[triangle[2]] - mesh.nodes[triangle[0]];
    return 0.5 * std::abs(a.x() * b.y() - a.y() * b.x());
}

double convergence_rate(
    double previous_error,
    double error,
    double previous_meshwidth,
    double meshwidth) {
    if (!(previous_error > 0.0 && error > 0.0
          && previous_meshwidth > meshwidth))
        return std::numeric_limits<double>::quiet_NaN();
    return std::log(previous_error / error)
         / std::log(previous_meshwidth / meshwidth);
}

double relative_residual(
    const ComplexSparseMatrix &matrix,
    const ComplexVector &solution,
    const ComplexVector &right_hand_side) {
    return (matrix * solution - right_hand_side).norm()
         / std::max(right_hand_side.norm(), 1e-30);
}

void verify_manufactured_solution(
    const HelmholtzManufacturedSolution &manufactured,
    double wavenumber) {
    for (double t : {0.0, 0.17, 0.53, 1.0}) {
        const std::array<std::pair<Point2, Eigen::Vector2d>, 4> points{{
            {Point2(0.0, t), Eigen::Vector2d(-1.0, 0.0)},
            {Point2(1.0, t), Eigen::Vector2d(1.0, 0.0)},
            {Point2(t, 0.0), Eigen::Vector2d(0.0, -1.0)},
            {Point2(t, 1.0), Eigen::Vector2d(0.0, 1.0)}}};
        for (const auto &[point, normal] : points) {
            const Complex robin =
                normal.cast<Complex>().dot(manufactured.gradient(point))
                - Complex(0.0, wavenumber) * manufactured.value(point);
            if (std::abs(robin) > 1e-11)
                throw std::runtime_error(
                    "manufactured solution violates homogeneous Robin data");
        }
    }
}

std::string level_tag(int level) {
    std::ostringstream stream;
    stream << "H_L" << std::setw(2) << std::setfill('0') << level;
    return stream.str();
}

fs::path temporary_path(const fs::path &path) {
    return fs::path(path.string() + ".tmp");
}

void publish_file(const fs::path &temporary, const fs::path &path) {
    if (fs::exists(path)) fs::remove(path);
    fs::rename(temporary, path);
}

void write_nodes(const fs::path &path, const TriMesh &mesh) {
    const fs::path temporary = temporary_path(path);
    std::ofstream output(temporary);
    if (!output) throw std::runtime_error("cannot write " + path.string());
    output << "node_id,x,y\n" << std::setprecision(17);
    for (int node = 0; node < static_cast<int>(mesh.nodes.size()); ++node)
        output << node << ',' << mesh.nodes[node].x() << ',' << mesh.nodes[node].y() << '\n';
    output.close();
    publish_file(temporary, path);
}

void write_elements(const fs::path &path, const TriMesh &mesh) {
    const fs::path temporary = temporary_path(path);
    std::ofstream output(temporary);
    if (!output) throw std::runtime_error("cannot write " + path.string());
    output << "element_id,n0,n1,n2,cx,cy,area,diameter\n" << std::setprecision(17);
    for (int element = 0; element < static_cast<int>(mesh.elems.size()); ++element) {
        const Triangle &triangle = mesh.elems[element];
        const Point2 centroid =
            (mesh.nodes[triangle[0]] + mesh.nodes[triangle[1]]
             + mesh.nodes[triangle[2]]) / 3.0;
        output << element << ',' << triangle[0] << ',' << triangle[1] << ','
               << triangle[2] << ',' << centroid.x() << ',' << centroid.y() << ','
               << triangle_area(mesh, triangle) << ','
               << triangle_diameter(mesh, triangle) << '\n';
    }
    output.close();
    publish_file(temporary, path);
}

void write_coarse_values(
    const fs::path &path,
    const TriMesh &mesh,
    const ComplexVector &p1,
    const ComplexVector &lod_coefficients) {
    const fs::path temporary = temporary_path(path);
    std::ofstream output(temporary);
    if (!output) throw std::runtime_error("cannot write " + path.string());
    output << "node_id,x,y,p1_real,p1_imag,lod_coefficient_real,"
              "lod_coefficient_imag\n"
           << std::setprecision(17);
    for (int node = 0; node < static_cast<int>(mesh.nodes.size()); ++node) {
        output << node << ',' << mesh.nodes[node].x() << ',' << mesh.nodes[node].y()
               << ',' << p1(node).real() << ',' << p1(node).imag()
               << ',' << lod_coefficients(node).real()
               << ',' << lod_coefficients(node).imag() << '\n';
    }
    output.close();
    publish_file(temporary, path);
}

void write_fields(
    const fs::path &path,
    const HelmholtzProblemData &problem,
    const HelmholtzLodSolution &lod,
    const ComplexVector &p1_on_fine,
    const HelmholtzManufacturedSolution &manufactured,
    const ComplexVector *fine_reference) {
    const ComplexVector coarse_scale =
        problem.coarse_to_fine.cast<Complex>() * lod.coarse_coefficients;
    const fs::path temporary = temporary_path(path);
    std::ofstream output(temporary);
    if (!output) throw std::runtime_error("cannot write " + path.string());
    output << "node_id,x,y,lod_real,lod_imag,lod_abs,coarse_scale_real,"
              "coarse_scale_imag,correction_real,correction_imag,p1_real,p1_imag,"
              "exact_real,exact_imag,lod_point_error_abs,p1_point_error_abs,"
              "fine_reference_real,fine_reference_imag\n"
           << std::setprecision(17);
    for (int node = 0; node < static_cast<int>(problem.fine.nodes.size()); ++node) {
        const Complex exact = manufactured.value(problem.fine.nodes[node]);
        const Complex correction = lod.fine_values(node) - coarse_scale(node);
        const Complex reference = fine_reference
            ? (*fine_reference)(node)
            : Complex(
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN());
        output << node << ',' << problem.fine.nodes[node].x() << ','
               << problem.fine.nodes[node].y() << ','
               << lod.fine_values(node).real() << ',' << lod.fine_values(node).imag()
               << ',' << std::abs(lod.fine_values(node)) << ','
               << coarse_scale(node).real() << ',' << coarse_scale(node).imag() << ','
               << correction.real() << ',' << correction.imag() << ','
               << p1_on_fine(node).real() << ',' << p1_on_fine(node).imag() << ','
               << exact.real() << ',' << exact.imag() << ','
               << std::abs(lod.fine_values(node) - exact) << ','
               << std::abs(p1_on_fine(node) - exact) << ','
               << reference.real() << ',' << reference.imag() << '\n';
    }
    output.close();
    publish_file(temporary, path);
}

double nesting_coordinate_residual(const HelmholtzProblemData &problem) {
    Eigen::VectorXd coarse_x(problem.coarse.nodes.size());
    Eigen::VectorXd coarse_y(problem.coarse.nodes.size());
    Eigen::VectorXd fine_x(problem.fine.nodes.size());
    Eigen::VectorXd fine_y(problem.fine.nodes.size());
    for (int i = 0; i < coarse_x.size(); ++i) {
        coarse_x(i) = problem.coarse.nodes[i].x();
        coarse_y(i) = problem.coarse.nodes[i].y();
    }
    for (int i = 0; i < fine_x.size(); ++i) {
        fine_x(i) = problem.fine.nodes[i].x();
        fine_y(i) = problem.fine.nodes[i].y();
    }
    return std::max(
        (problem.coarse_to_fine * coarse_x - fine_x).cwiseAbs().maxCoeff(),
        (problem.coarse_to_fine * coarse_y - fine_y).cwiseAbs().maxCoeff());
}

const char *csv_header() {
    return "k,H_level,h_level,ell,solver,mode,coarse_nodes,coarse_elements,"
           "fine_nodes,fine_elements,H_max,H_min,h_max,h_min,kH,kh,h_over_H,"
           "p1_energy_abs,lod_energy_abs,p1_l2_abs,lod_l2_abs,"
           "fine_energy_abs,fine_l2_abs,p1_energy_rate,lod_energy_rate,"
           "p1_l2_rate,lod_l2_rate,p1_residual,petrov_residual,"
           "corrector_residual,constraint_residual,nesting_coordinate_residual,"
           "corrector_threads,symbolic_analyses,symbolic_reuses,"
           "factorization_reuses,schur_residual,schur_rcond,"
           "mesh_ms,operators_ms,correctors_ms,basis_ms,build_ms,load_ms,"
           "lod_solve_ms,p1_solve_ms,fine_solve_ms,error_ms,export_ms,total_ms,"
           "coarse_nodes_file,coarse_elements_file,fine_nodes_file,"
           "fine_elements_file,fields_file";
}

void print_csv_row(std::ostream &output, const Row &row) {
    output << std::setprecision(17)
           << row.k << ',' << row.H_level << ',' << row.h_level << ',' << row.ell << ','
           << row.solver << ',' << row.mode << ',' << row.coarse_nodes << ','
           << row.coarse_elements << ',' << row.fine_nodes << ',' << row.fine_elements << ','
           << row.H_max << ',' << row.H_min << ',' << row.h_max << ',' << row.h_min << ','
           << row.kH << ',' << row.kh << ',' << row.h_over_H << ','
           << row.p1_exact.energy << ',' << row.lod_exact.energy << ','
           << row.p1_exact.l2 << ',' << row.lod_exact.l2 << ','
           << row.fine_exact.energy << ',' << row.fine_exact.l2 << ','
           << row.p1_energy_rate << ',' << row.lod_energy_rate << ','
           << row.p1_l2_rate << ',' << row.lod_l2_rate << ','
           << row.p1_residual << ',' << row.petrov_residual << ','
           << row.corrector_residual << ',' << row.constraint_residual << ','
           << row.nesting_coordinate_residual << ',' << row.corrector_threads << ','
           << row.symbolic_analyses << ',' << row.symbolic_reuses << ','
           << row.factorization_reuses << ',' << row.schur_residual << ','
           << row.schur_rcond << ',' << row.mesh_ms << ',' << row.operators_ms << ','
           << row.correctors_ms << ',' << row.basis_ms << ',' << row.build_ms << ','
           << row.load_ms << ',' << row.lod_solve_ms << ',' << row.p1_solve_ms << ','
           << row.fine_solve_ms << ',' << row.error_ms << ',' << row.export_ms << ','
           << row.total_ms << ',' << row.coarse_nodes_file << ','
           << row.coarse_elements_file << ',' << row.fine_nodes_file << ','
           << row.fine_elements_file << ',' << row.fields_file << '\n';
}

void print_human(const std::vector<Row> &rows) {
    std::cout << std::setprecision(8) << std::scientific
              << "Absolute manufactured-solution errors in "
                 "||v||_{1,k}=(||grad v||^2+k^2||v||^2)^(1/2)\n"
              << "Hlev H h kH h/H Nc Nf P1_E rate LOD_E rate "
                 "P1_L2 LOD_L2 build_ms\n";
    for (const Row &row : rows) {
        std::cout << row.H_level << ' ' << row.H_max << ' ' << row.h_max << ' '
                  << row.kH << ' ' << row.h_over_H << ' ' << row.coarse_nodes << ' '
                  << row.fine_nodes << ' ' << row.p1_exact.energy << ' '
                  << row.p1_energy_rate << ' ' << row.lod_exact.energy << ' '
                  << row.lod_energy_rate << ' ' << row.p1_exact.l2 << ' '
                  << row.lod_exact.l2 << ' ' << row.build_ms << '\n';
    }
}

void export_snapshot(
    const Options &options,
    Row &row,
    const HelmholtzLodModel &model,
    const ComplexVector &p1_values,
    const ComplexVector &p1_on_fine,
    const HelmholtzLodSolution &lod,
    const HelmholtzManufacturedSolution &manufactured,
    const ComplexVector *fine_reference) {
    if (options.export_dir.empty()) return;
    fs::create_directories(options.export_dir);
    const std::string tag = level_tag(row.H_level);
    row.coarse_nodes_file = tag + "_coarse_nodes.csv";
    row.coarse_elements_file = tag + "_coarse_elements.csv";
    row.fine_nodes_file = "fine_L" + std::to_string(row.h_level) + "_nodes.csv";
    row.fine_elements_file = "fine_L" + std::to_string(row.h_level) + "_elements.csv";
    row.fields_file = options.export_fields ? tag + "_fields.csv" : "";

    write_coarse_values(
        options.export_dir / row.coarse_nodes_file,
        model.problem().coarse, p1_values, lod.coarse_coefficients);
    write_elements(
        options.export_dir / row.coarse_elements_file,
        model.problem().coarse);
    if (!fs::exists(options.export_dir / row.fine_nodes_file))
        write_nodes(options.export_dir / row.fine_nodes_file, model.problem().fine);
    if (!fs::exists(options.export_dir / row.fine_elements_file))
        write_elements(options.export_dir / row.fine_elements_file, model.problem().fine);
    if (options.export_fields) {
        write_fields(
            options.export_dir / row.fields_file,
            model.problem(), lod, p1_on_fine, manufactured, fine_reference);
    }
}

Row run_level(
    const Options &options,
    int coarse_level,
    const HelmholtzManufacturedSolution &manufactured) {
    const auto total_start = std::chrono::steady_clock::now();
    Row row;
    row.k = options.wavenumber;
    row.H_level = coarse_level;
    row.h_level = options.fine_level;
    row.ell = options.ell;
    row.solver = options.solver;
    row.mode = options.mode == HelmholtzPetrovMode::TwoSided
        ? "two-sided" : "test-only";

    HelmholtzProblemConfig config;
    config.H = coarse_level;
    config.h = options.fine_level;
    config.ell = options.ell;
    config.wavenumber = options.wavenumber;
    config.mode = options.mode;
    config.patch_solver.kind = options.solver == "schur"
        ? HelmholtzPatchSolverKind::DirectSchur
        : HelmholtzPatchSolverKind::DirectSaddle;
    config.patch_solver.symbolic_cache_slots = options.symbolic_cache_slots;
    config.patch_solver.reuse_identical_factorization =
        options.reuse_identical_factorization;
    config.patch_solver.fallback_to_direct = false;

    HelmholtzLodModel model = HelmholtzLodModel::build(config);
    const HelmholtzBuildTimings &build = model.build_timings();
    row.mesh_ms = build.mesh_and_interpolation_ms;
    row.operators_ms = build.operators_ms;
    row.correctors_ms = build.correctors_ms;
    row.basis_ms = build.basis_and_factorization_ms;
    row.build_ms = build.total_ms;
    row.coarse_nodes = static_cast<int>(model.problem().coarse.nodes.size());
    row.coarse_elements = static_cast<int>(model.problem().coarse.elems.size());
    row.fine_nodes = static_cast<int>(model.problem().fine.nodes.size());
    row.fine_elements = static_cast<int>(model.problem().fine.elems.size());
    row.H_max = max_element_diameter(model.problem().coarse);
    row.H_min = min_element_diameter(model.problem().coarse);
    row.h_max = max_element_diameter(model.problem().fine);
    row.h_min = min_element_diameter(model.problem().fine);
    row.kH = row.k * row.H_max;
    row.kh = row.k * row.h_max;
    row.h_over_H = row.h_max / row.H_max;
    row.nesting_coordinate_residual = nesting_coordinate_residual(model.problem());

    auto stage_start = std::chrono::steady_clock::now();
    const ComplexVector fine_load =
        assemble_helmholtz_load(model.problem().fine, manufactured.source);
    row.load_ms = elapsed_ms(stage_start);

    stage_start = std::chrono::steady_clock::now();
    const HelmholtzLodSolution lod = model.solve_load(fine_load);
    row.lod_solve_ms = elapsed_ms(stage_start);
    row.petrov_residual = lod.petrov_residual;

    stage_start = std::chrono::steady_clock::now();
    const HelmholtzOperators p1_operators =
        assemble_helmholtz_operators(model.problem().coarse, row.k);
    const ComplexVector p1_load =
        assemble_helmholtz_load(model.problem().coarse, manufactured.source);
    const ComplexVector p1_values = solve_helmholtz_fem(p1_operators, p1_load);
    const ComplexVector p1_on_fine =
        model.problem().coarse_to_fine.cast<Complex>() * p1_values;
    row.p1_solve_ms = elapsed_ms(stage_start);
    row.p1_residual = relative_residual(
        p1_operators.system, p1_values, p1_load);

    ComplexVector fine_reference;
    const ComplexVector *fine_reference_pointer = nullptr;
    if (options.solve_fine_reference) {
        stage_start = std::chrono::steady_clock::now();
        fine_reference = model.solve_fine_reference(fine_load);
        row.fine_solve_ms = elapsed_ms(stage_start);
        fine_reference_pointer = &fine_reference;
    }

    stage_start = std::chrono::steady_clock::now();
    row.p1_exact = compute_helmholtz_error(
        model.problem().fine, p1_on_fine, row.k,
        manufactured.value, manufactured.gradient);
    row.lod_exact = compute_helmholtz_error(
        model.problem().fine, lod.fine_values, row.k,
        manufactured.value, manufactured.gradient);
    if (fine_reference_pointer) {
        row.fine_exact = compute_helmholtz_error(
            model.problem().fine, fine_reference, row.k,
            manufactured.value, manufactured.gradient);
    }
    row.error_ms = elapsed_ms(stage_start);

    const HelmholtzCorrectorDiagnostics &diagnostics =
        model.correctors().diagnostics;
    row.corrector_residual = diagnostics.max_primal_residual;
    row.constraint_residual = diagnostics.max_constraint_residual;
    row.corrector_threads = diagnostics.parallel_threads;
    row.symbolic_analyses = diagnostics.symbolic_analyses;
    row.symbolic_reuses = diagnostics.symbolic_reuses;
    row.factorization_reuses = diagnostics.factorization_reuses;
    row.schur_residual = diagnostics.max_schur_residual;
    row.schur_rcond = diagnostics.min_schur_reciprocal_condition;

    stage_start = std::chrono::steady_clock::now();
    export_snapshot(
        options, row, model, p1_values, p1_on_fine, lod,
        manufactured, fine_reference_pointer);
    row.export_ms = elapsed_ms(stage_start);
    row.total_ms = elapsed_ms(total_start);
    return row;
}

void check_row(const Options &options, const Row &row) {
    const std::array<double, 7> positive{{
        row.H_max, row.H_min, row.h_max, row.h_min,
        row.p1_exact.energy, row.lod_exact.energy, row.lod_exact.l2}};
    for (double value : positive) {
        if (!(std::isfinite(value) && value > 0.0))
            throw std::runtime_error("non-finite or non-positive convergence metric");
    }
    if (!(row.p1_residual < 1e-8
          && row.petrov_residual < 1e-8
          && row.corrector_residual < 1e-8
          && row.constraint_residual < 1e-8))
        throw std::runtime_error("algebraic residual check failed");
    if (!(row.nesting_coordinate_residual < 1e-13))
        throw std::runtime_error("V_H is not represented exactly in V_h");
    if (options.solver == "schur"
        && (!(row.schur_residual < 1e-8) || !(row.schur_rcond > 1e-14)))
        throw std::runtime_error("DirectSchur diagnostics check failed");
}

} // namespace

int main(int argc, char **argv) {
    try {
        const Options options = parse_options(argc, argv);
        if (options.csv_header) {
            std::cout << csv_header() << '\n';
            return 0;
        }
#ifdef _OPENMP
        if (options.threads > 0) omp_set_num_threads(options.threads);
#endif
        const HelmholtzManufacturedSolution manufactured =
            make_polynomial_plane_wave_solution(options.wavenumber);
        verify_manufactured_solution(manufactured, options.wavenumber);

        std::vector<Row> rows;
        rows.reserve(options.coarse_levels.size());
        for (int level : options.coarse_levels) {
            Row row = run_level(options, level, manufactured);
            if (!rows.empty()) {
                row.p1_energy_rate = convergence_rate(
                    rows.back().p1_exact.energy, row.p1_exact.energy,
                    rows.back().H_max, row.H_max);
                row.lod_energy_rate = convergence_rate(
                    rows.back().lod_exact.energy, row.lod_exact.energy,
                    rows.back().H_max, row.H_max);
                row.p1_l2_rate = convergence_rate(
                    rows.back().p1_exact.l2, row.p1_exact.l2,
                    rows.back().H_max, row.H_max);
                row.lod_l2_rate = convergence_rate(
                    rows.back().lod_exact.l2, row.lod_exact.l2,
                    rows.back().H_max, row.H_max);
            }
            if (options.check) check_row(options, row);
            rows.push_back(std::move(row));
        }
        if (options.check && rows.size() > 1
            && (!(rows.back().lod_exact.energy < rows.front().lod_exact.energy)
                || !(rows.back().p1_exact.energy < rows.front().p1_exact.energy)))
            throw std::runtime_error("absolute energy errors did not decrease");

        if (!options.summary_out.empty()) {
            if (!options.summary_out.parent_path().empty())
                fs::create_directories(options.summary_out.parent_path());
            std::ofstream output(options.summary_out);
            if (!output)
                throw std::runtime_error(
                    "cannot write " + options.summary_out.string());
            output << csv_header() << '\n';
            for (const Row &row : rows) print_csv_row(output, row);
        }
        if (options.csv) {
            std::cout << csv_header() << '\n';
            for (const Row &row : rows) print_csv_row(std::cout, row);
        } else {
            print_human(rows);
        }
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "bench_helmholtz_H_convergence failed: "
                  << error.what() << '\n';
        return 1;
    }
}

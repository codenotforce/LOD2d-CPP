#include "helmholtz/manufactured.h"
#include "helmholtz/model.h"
#include "helmholtz/operators.h"
#include "mesh/refine.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace lod2d;
using namespace lod2d::helmholtz;

namespace {

struct Options {
    double wavenumber = 4.0;
    int fine_level = 10;
    std::vector<int> coarse_levels{3, 4, 5, 6};
    std::vector<int> oversampling_levels{2, 3};
    std::vector<int> fem_levels{6, 8, 10};
    HelmholtzPetrovMode petrov_mode = HelmholtzPetrovMode::TwoSided;
    bool csv = false;
    bool check = false;
};

struct DiscreteError {
    double energy = 0.0;
    double l2 = 0.0;
};

struct FemRow {
    int level = 0;
    int nodes = 0;
    int elements = 0;
    double meshwidth = 0.0;
    HelmholtzError exact;
    double energy_rate = std::numeric_limits<double>::quiet_NaN();
    double l2_rate = std::numeric_limits<double>::quiet_NaN();
    double relative_residual = 0.0;
};

struct LodRow {
    int coarse_level = 0;
    int fine_level = 0;
    int ell = 0;
    int coarse_nodes = 0;
    int coarse_elements = 0;
    int fine_nodes = 0;
    int fine_elements = 0;
    double H = 0.0;
    double h = 0.0;
    HelmholtzError fem_exact;
    HelmholtzError lod_exact;
    DiscreteError lod_fem;
    double exact_energy_rate = std::numeric_limits<double>::quiet_NaN();
    double exact_l2_rate = std::numeric_limits<double>::quiet_NaN();
    double petrov_residual = 0.0;
    double corrector_residual = 0.0;
    double constraint_residual = 0.0;
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

std::vector<int> parse_int_list(const std::string &text, const char *name) {
    std::vector<int> values;
    std::stringstream stream(text);
    std::string token;
    while (std::getline(stream, token, ',')) {
        if (token.empty()) throw std::invalid_argument(std::string("empty entry in ") + name);
        values.push_back(parse_int(token, name));
    }
    if (values.empty()) throw std::invalid_argument(std::string(name) + " must not be empty");
    return values;
}

Options parse_options(int argc, char **argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        auto value_after = [&](const std::string &prefix) {
            return argument.substr(prefix.size());
        };
        if (argument.rfind("--k=", 0) == 0) {
            options.wavenumber = parse_double(value_after("--k="), "wavenumber");
        } else if (argument.rfind("--h=", 0) == 0) {
            options.fine_level = parse_int(value_after("--h="), "fine level");
        } else if (argument.rfind("--H-levels=", 0) == 0) {
            options.coarse_levels = parse_int_list(value_after("--H-levels="), "H levels");
        } else if (argument.rfind("--ell-levels=", 0) == 0) {
            options.oversampling_levels =
                parse_int_list(value_after("--ell-levels="), "ell levels");
        } else if (argument.rfind("--fem-levels=", 0) == 0) {
            options.fem_levels = parse_int_list(value_after("--fem-levels="), "FEM levels");
        } else if (argument.rfind("--mode=", 0) == 0) {
            const std::string mode = value_after("--mode=");
            if (mode == "two-sided") options.petrov_mode = HelmholtzPetrovMode::TwoSided;
            else if (mode == "test-only")
                options.petrov_mode = HelmholtzPetrovMode::CorrectedTestOnly;
            else throw std::invalid_argument("mode must be two-sided or test-only");
        } else if (argument == "--format=csv") {
            options.csv = true;
        } else if (argument == "--check") {
            options.check = true;
        } else if (argument == "--help") {
            std::cout
                << "Usage: bench_helmholtz_manufactured [--k=4] [--h=10] "
                   "[--H-levels=3,4,5,6] [--ell-levels=2,3] "
                   "[--fem-levels=6,8,10] [--mode=two-sided|test-only] "
                   "[--format=csv] [--check]\n\n"
                   "All level arguments are counts of global NVB sweeps. The output also "
                   "reports the measured maximum element diameters H and h.\n";
            std::exit(0);
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (!(options.wavenumber > 0.0)) throw std::invalid_argument("k must be positive");
    if (options.fine_level < 0) throw std::invalid_argument("h must be nonnegative");
    for (int level : options.coarse_levels) {
        if (level < 0 || level > options.fine_level)
            throw std::invalid_argument("each H level must satisfy 0 <= H <= h");
    }
    for (int ell : options.oversampling_levels) {
        if (ell < 0) throw std::invalid_argument("ell levels must be nonnegative");
    }
    for (int level : options.fem_levels) {
        if (level < 0) throw std::invalid_argument("FEM levels must be nonnegative");
    }
    std::sort(options.coarse_levels.begin(), options.coarse_levels.end());
    std::sort(options.oversampling_levels.begin(), options.oversampling_levels.end());
    std::sort(options.fem_levels.begin(), options.fem_levels.end());
    return options;
}

double quadratic_form(
    const Eigen::SparseMatrix<double> &matrix,
    const ComplexVector &vector) {
    return std::max(0.0, std::real(vector.dot(matrix.cast<Complex>() * vector)));
}

DiscreteError discrete_error(
    const HelmholtzOperators &operators,
    const ComplexVector &approximation,
    const ComplexVector &reference) {
    const ComplexVector difference = approximation - reference;
    DiscreteError result;
    result.l2 = std::sqrt(quadratic_form(operators.mass, difference));
    result.energy = std::sqrt(
        quadratic_form(operators.stiffness, difference)
        + operators.wavenumber * operators.wavenumber * result.l2 * result.l2);
    return result;
}

double convergence_rate(double previous_error, double error, double previous_h, double h) {
    if (!(previous_error > 0.0 && error > 0.0 && previous_h > h))
        return std::numeric_limits<double>::quiet_NaN();
    return std::log(previous_error / error) / std::log(previous_h / h);
}

double relative_linear_residual(
    const ComplexSparseMatrix &matrix,
    const ComplexVector &solution,
    const ComplexVector &right_hand_side) {
    return (matrix * solution - right_hand_side).norm()
         / std::max(right_hand_side.norm(), 1e-30);
}

void verify_homogeneous_robin_data(
    const HelmholtzManufacturedSolution &manufactured,
    double wavenumber) {
    for (double t : {0.0, 0.17, 0.53, 1.0}) {
        const std::array<std::pair<Point2, Eigen::Vector2d>, 4> points{{
            {Point2(0.0, t), Eigen::Vector2d(-1.0, 0.0)},
            {Point2(1.0, t), Eigen::Vector2d(1.0, 0.0)},
            {Point2(t, 0.0), Eigen::Vector2d(0.0, -1.0)},
            {Point2(t, 1.0), Eigen::Vector2d(0.0, 1.0)}}};
        for (const auto &[point, normal] : points) {
            const Complex robin = normal.cast<Complex>().dot(manufactured.gradient(point))
                                - Complex(0.0, wavenumber) * manufactured.value(point);
            if (std::abs(robin) > 1e-11)
                throw std::runtime_error("manufactured solution violates homogeneous Robin data");
        }
    }
}

std::vector<FemRow> run_fem_study(
    const Options &options,
    const HelmholtzManufacturedSolution &manufactured) {
    std::vector<FemRow> rows;
    const TriMesh initial = make_helmholtz_unit_square_mesh();
    for (int level : options.fem_levels) {
        const TriMesh mesh = refine_mesh_nvb(initial, level).mesh;
        const HelmholtzOperators operators =
            assemble_helmholtz_operators(mesh, options.wavenumber);
        const ComplexVector load = assemble_helmholtz_load(mesh, manufactured.source);
        const ComplexVector solution = solve_helmholtz_fem(operators, load);
        FemRow row;
        row.level = level;
        row.nodes = static_cast<int>(mesh.nodes.size());
        row.elements = static_cast<int>(mesh.elems.size());
        row.meshwidth = max_element_diameter(mesh);
        row.exact = compute_helmholtz_error(
            mesh, solution, options.wavenumber,
            manufactured.value, manufactured.gradient);
        row.relative_residual = relative_linear_residual(operators.system, solution, load);
        if (!rows.empty()) {
            row.energy_rate = convergence_rate(
                rows.back().exact.energy, row.exact.energy,
                rows.back().meshwidth, row.meshwidth);
            row.l2_rate = convergence_rate(
                rows.back().exact.l2, row.exact.l2,
                rows.back().meshwidth, row.meshwidth);
        }
        rows.push_back(row);
    }
    return rows;
}

std::vector<LodRow> run_lod_study(
    const Options &options,
    const HelmholtzManufacturedSolution &manufactured) {
    std::vector<LodRow> rows;
    for (int ell : options.oversampling_levels) {
        LodRow previous;
        bool have_previous = false;
        for (int coarse_level : options.coarse_levels) {
            HelmholtzProblemConfig config;
            config.H = coarse_level;
            config.h = options.fine_level;
            config.ell = ell;
            config.wavenumber = options.wavenumber;
            config.mode = options.petrov_mode;
            HelmholtzLodModel model = HelmholtzLodModel::build(config);

            const ComplexVector load =
                assemble_helmholtz_load(model.problem().fine, manufactured.source);
            const HelmholtzLodSolution lod = model.solve_load(load);
            const ComplexVector fem = model.solve_fine_reference(load);

            LodRow row;
            row.coarse_level = coarse_level;
            row.fine_level = options.fine_level;
            row.ell = ell;
            row.coarse_nodes = static_cast<int>(model.problem().coarse.nodes.size());
            row.coarse_elements = static_cast<int>(model.problem().coarse.elems.size());
            row.fine_nodes = static_cast<int>(model.problem().fine.nodes.size());
            row.fine_elements = static_cast<int>(model.problem().fine.elems.size());
            row.H = max_element_diameter(model.problem().coarse);
            row.h = max_element_diameter(model.problem().fine);
            row.fem_exact = compute_helmholtz_error(
                model.problem().fine, fem, options.wavenumber,
                manufactured.value, manufactured.gradient);
            row.lod_exact = compute_helmholtz_error(
                model.problem().fine, lod.fine_values, options.wavenumber,
                manufactured.value, manufactured.gradient);
            row.lod_fem = discrete_error(model.operators(), lod.fine_values, fem);
            row.petrov_residual = lod.petrov_residual;
            row.corrector_residual = model.correctors().diagnostics.max_primal_residual;
            row.constraint_residual = model.correctors().diagnostics.max_constraint_residual;
            if (have_previous) {
                row.exact_energy_rate = convergence_rate(
                    previous.lod_exact.energy, row.lod_exact.energy,
                    previous.H, row.H);
                row.exact_l2_rate = convergence_rate(
                    previous.lod_exact.l2, row.lod_exact.l2,
                    previous.H, row.H);
            }
            previous = row;
            have_previous = true;
            rows.push_back(row);
        }
    }
    return rows;
}

void print_human(const std::vector<FemRow> &fem, const std::vector<LodRow> &lod) {
    std::cout << std::setprecision(6) << std::scientific;
    std::cout << "FEM manufactured-solution convergence (global NVB)\n"
              << "level nodes elems h energy_error energy_rate l2_error l2_rate residual\n";
    for (const FemRow &row : fem) {
        std::cout << row.level << ' ' << row.nodes << ' ' << row.elements << ' '
                  << row.meshwidth << ' ' << row.exact.energy << ' ' << row.energy_rate << ' '
                  << row.exact.l2 << ' ' << row.l2_rate << ' ' << row.relative_residual << '\n';
    }
    std::cout << "\nLOD manufactured-solution convergence (fixed fine grid)\n"
              << "Hlev hlev ell H h Nc Nf FEM_E LOD_E rate_E LOD_L2 rate_L2 "
                 "LOD-FEM_E LOD-FEM_L2 petrov corrector constraint\n";
    for (const LodRow &row : lod) {
        std::cout << row.coarse_level << ' ' << row.fine_level << ' ' << row.ell << ' '
                  << row.H << ' ' << row.h << ' ' << row.coarse_nodes << ' ' << row.fine_nodes << ' '
                  << row.fem_exact.energy << ' ' << row.lod_exact.energy << ' '
                  << row.exact_energy_rate << ' ' << row.lod_exact.l2 << ' '
                  << row.exact_l2_rate << ' ' << row.lod_fem.energy << ' '
                  << row.lod_fem.l2 << ' ' << row.petrov_residual << ' '
                  << row.corrector_residual << ' ' << row.constraint_residual << '\n';
    }
}

void print_csv(
    double wavenumber,
    const std::vector<FemRow> &fem,
    const std::vector<LodRow> &lod) {
    std::cout << "study,k,H_level,h_level,ell,H,h,coarse_nodes,coarse_elements,fine_nodes,"
                 "fine_elements,fem_exact_energy,fem_exact_l2,lod_exact_energy,lod_exact_l2,"
                 "energy_rate,l2_rate,lod_fem_energy,lod_fem_l2,petrov_residual,"
                 "corrector_residual,constraint_residual,linear_residual\n";
    std::cout << std::setprecision(17);
    for (const FemRow &row : fem) {
        std::cout << "fem," << wavenumber << ',' << -1 << ',' << row.level << ',' << -1 << ','
                  << -1 << ',' << row.meshwidth << ',' << -1 << ',' << -1 << ','
                  << row.nodes << ',' << row.elements << ',' << row.exact.energy << ','
                  << row.exact.l2 << ",nan,nan," << row.energy_rate << ',' << row.l2_rate
                  << ",nan,nan,nan,nan,nan," << row.relative_residual << '\n';
    }
    for (const LodRow &row : lod) {
        std::cout << "lod," << wavenumber << ',' << row.coarse_level << ',' << row.fine_level << ','
                  << row.ell << ',' << row.H << ',' << row.h << ',' << row.coarse_nodes << ','
                  << row.coarse_elements << ',' << row.fine_nodes << ',' << row.fine_elements << ','
                  << row.fem_exact.energy << ',' << row.fem_exact.l2 << ','
                  << row.lod_exact.energy << ',' << row.lod_exact.l2 << ','
                  << row.exact_energy_rate << ',' << row.exact_l2_rate << ','
                  << row.lod_fem.energy << ',' << row.lod_fem.l2 << ','
                  << row.petrov_residual << ',' << row.corrector_residual << ','
                  << row.constraint_residual << ",nan\n";
    }
}

void check_results(const std::vector<FemRow> &fem, const std::vector<LodRow> &lod) {
    if (fem.size() >= 2) {
        if (!(fem.back().exact.energy < fem.front().exact.energy
              && fem.back().exact.l2 < fem.front().exact.l2))
            throw std::runtime_error("fine FEM error did not decrease under global NVB");
        const double energy_rate = convergence_rate(
            fem.front().exact.energy, fem.back().exact.energy,
            fem.front().meshwidth, fem.back().meshwidth);
        const double l2_rate = convergence_rate(
            fem.front().exact.l2, fem.back().exact.l2,
            fem.front().meshwidth, fem.back().meshwidth);
        if (!(energy_rate > 0.7 && l2_rate > 1.5))
            throw std::runtime_error("fine FEM manufactured-solution rates are too low");
    }
    for (const FemRow &row : fem) {
        if (!(row.relative_residual < 1e-9))
            throw std::runtime_error("fine FEM algebraic residual is too large");
    }
    for (const LodRow &row : lod) {
        if (!(row.petrov_residual < 1e-8
              && row.corrector_residual < 1e-8
              && row.constraint_residual < 1e-8))
            throw std::runtime_error("LOD algebraic residual is too large");
    }
    for (std::size_t first = 0; first < lod.size();) {
        std::size_t last = first;
        while (last + 1 < lod.size() && lod[last + 1].ell == lod[first].ell) ++last;
        if (last > first
            && !(lod[last].lod_exact.energy < lod[first].lod_exact.energy
                 && lod[last].lod_exact.l2 < lod[first].lod_exact.l2
                 && lod[last].lod_fem.energy < lod[first].lod_fem.energy
                 && lod[last].lod_fem.l2 < lod[first].lod_fem.l2))
            throw std::runtime_error("LOD error did not decrease as H was refined");
        first = last + 1;
    }
}

} // namespace

int main(int argc, char **argv) {
    try {
        const Options options = parse_options(argc, argv);
        const HelmholtzManufacturedSolution manufactured =
            make_polynomial_plane_wave_solution(options.wavenumber);
        verify_homogeneous_robin_data(manufactured, options.wavenumber);
        const std::vector<FemRow> fem = run_fem_study(options, manufactured);
        const std::vector<LodRow> lod = run_lod_study(options, manufactured);
        if (options.csv) print_csv(options.wavenumber, fem, lod);
        else print_human(fem, lod);
        if (options.check) check_results(fem, lod);
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "bench_helmholtz_manufactured failed: " << error.what() << '\n';
        return 1;
    }
}

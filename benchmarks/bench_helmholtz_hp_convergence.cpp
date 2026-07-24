#include "helmholtz/hp_model.h"
#include "helmholtz/manufactured.h"
#include "mesh/refine.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace lod2d;
using namespace lod2d::helmholtz;

namespace {

struct Options {
    std::string study = "fem";
    std::string solver = "saddle";
    double wavenumber = 4.0;
    std::vector<int> degrees{1, 2, 3};
    std::vector<int> fine_levels{4, 6, 8};
    std::vector<int> coarse_levels{2, 4, 6};
    std::vector<int> ell_levels{1, 2, 3};
    std::vector<int> ell_by_degree{2, 2, 2};
    int coarse_level = 4;
    int fine_level = 8;
    int gap = 4;
    int threads = 1;
    int progress_interval = 0;
    bool stream = false;
    bool check = false;
};

struct ErrorPair {
    double energy = 0.0;
    double l2 = 0.0;
};

struct FineReference {
    ComplexVector solution;
    ErrorPair exact;
};

struct Row {
    std::string study;
    std::string solver;
    int degree = 0;
    int H_level = -1;
    int h_level = -1;
    int ell = -1;
    int fine_dofs = 0;
    double H = std::numeric_limits<double>::quiet_NaN();
    double h = std::numeric_limits<double>::quiet_NaN();
    ErrorPair exact;
    ErrorPair fine_exact;
    ErrorPair lod_fine;
    double delta = std::numeric_limits<double>::quiet_NaN();
    double energy_rate = std::numeric_limits<double>::quiet_NaN();
    double l2_rate = std::numeric_limits<double>::quiet_NaN();
    double petrov_residual = 0.0;
    double corrector_residual = 0.0;
    double constraint_residual = 0.0;
    double schur_residual = 0.0;
    double schur_rcond = 1.0;
    int direct_fallbacks = 0;
    int corrector_threads = 1;
};

void print_header();
void print_row(const Row &row);

int parse_int(const std::string &text, const char *name) {
    std::size_t consumed = 0;
    const int result = std::stoi(text, &consumed);
    if (consumed != text.size())
        throw std::invalid_argument(std::string("invalid ") + name);
    return result;
}

double parse_double(const std::string &text, const char *name) {
    std::size_t consumed = 0;
    const double result = std::stod(text, &consumed);
    if (consumed != text.size() || !std::isfinite(result))
        throw std::invalid_argument(std::string("invalid ") + name);
    return result;
}

std::vector<int> parse_list(const std::string &text, const char *name) {
    std::vector<int> result;
    std::stringstream stream(text);
    std::string token;
    while (std::getline(stream, token, ','))
        result.push_back(parse_int(token, name));
    if (result.empty())
        throw std::invalid_argument(std::string(name) + " must not be empty");
    return result;
}

Options parse_options(int argc, char **argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        auto after = [&](const std::string &prefix) {
            return argument.substr(prefix.size());
        };
        if (argument.rfind("--study=", 0) == 0)
            options.study = after("--study=");
        else if (argument.rfind("--solver=", 0) == 0)
            options.solver = after("--solver=");
        else if (argument.rfind("--k=", 0) == 0)
            options.wavenumber = parse_double(after("--k="), "k");
        else if (argument.rfind("--p=", 0) == 0)
            options.degrees = parse_list(after("--p="), "p");
        else if (argument.rfind("--h-levels=", 0) == 0)
            options.fine_levels =
                parse_list(after("--h-levels="), "h levels");
        else if (argument.rfind("--H-levels=", 0) == 0)
            options.coarse_levels =
                parse_list(after("--H-levels="), "H levels");
        else if (argument.rfind("--ell-levels=", 0) == 0)
            options.ell_levels =
                parse_list(after("--ell-levels="), "ell levels");
        else if (argument.rfind("--ell-by-p=", 0) == 0)
            options.ell_by_degree =
                parse_list(after("--ell-by-p="), "ell by p");
        else if (argument.rfind("--H=", 0) == 0)
            options.coarse_level = parse_int(after("--H="), "H");
        else if (argument.rfind("--h=", 0) == 0)
            options.fine_level = parse_int(after("--h="), "h");
        else if (argument.rfind("--gap=", 0) == 0)
            options.gap = parse_int(after("--gap="), "gap");
        else if (argument.rfind("--threads=", 0) == 0)
            options.threads = parse_int(after("--threads="), "threads");
        else if (argument.rfind("--progress=", 0) == 0)
            options.progress_interval =
                parse_int(after("--progress="), "progress");
        else if (argument == "--stream")
            options.stream = true;
        else if (argument == "--check")
            options.check = true;
        else if (argument == "--help") {
            std::cout
                << "Usage: bench_helmholtz_hp_convergence "
                   "--study=fem|ell|H|coupled "
                   "[--solver=saddle|schur] [--k=4] [--p=1,2,3] "
                   "[--h-levels=4,6,8] [--H-levels=2,4,6] "
                   "[--ell-levels=1,2,3] [--ell-by-p=2,2,2] "
                   "[--H=4] [--h=8] [--gap=4] [--threads=1] "
                   "[--progress=0] [--stream] [--check]\n";
            std::exit(0);
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.study != "fem" && options.study != "ell"
        && options.study != "H" && options.study != "coupled")
        throw std::invalid_argument("study must be fem, ell, H, or coupled");
    if (options.solver != "saddle" && options.solver != "schur")
        throw std::invalid_argument("solver must be saddle or schur");
    if (options.ell_by_degree.size() != options.degrees.size())
        throw std::invalid_argument(
            "ell-by-p must have one entry for each requested p");
    for (int degree : options.degrees)
        if (degree < 1 || degree > 3)
            throw std::invalid_argument("p must be 1, 2, or 3");
    if (options.threads < 1)
        throw std::invalid_argument("threads must be positive");
    if (options.progress_interval < 0)
        throw std::invalid_argument("progress must be nonnegative");
    return options;
}

double quadratic_form(
    const Eigen::SparseMatrix<double> &matrix,
    const ComplexVector &vector) {
    return std::max(
        0.0, std::real(vector.dot(matrix.cast<Complex>() * vector)));
}

ErrorPair discrete_error(
    const HelmholtzHpOperators &operators,
    const ComplexVector &approximation,
    const ComplexVector &reference) {
    const ComplexVector difference = approximation - reference;
    ErrorPair result;
    result.l2 = std::sqrt(quadratic_form(operators.mass, difference));
    result.energy = std::sqrt(
        quadratic_form(operators.stiffness, difference)
        + operators.wavenumber * operators.wavenumber
            * result.l2 * result.l2);
    return result;
}

FineReference build_fine_reference(
    const Options &options,
    const HelmholtzManufacturedSolution &manufactured,
    int degree,
    int h_level) {
    TriMesh mesh = refine_mesh_nvb(
        make_helmholtz_unit_square_mesh(), h_level).mesh;
    HpTriSpace space(mesh, degree);
    const HelmholtzHpOperators operators =
        assemble_helmholtz_hp_operators(space, options.wavenumber);
    const ComplexVector load =
        assemble_helmholtz_hp_load(space, manufactured.source);
    FineReference reference;
    reference.solution = solve_helmholtz_hp_fem(operators, load);
    const HelmholtzError exact = compute_helmholtz_hp_error(
        space, reference.solution, options.wavenumber,
        manufactured.value, manufactured.gradient);
    reference.exact = {exact.energy, exact.l2};
    return reference;
}
double rate(double previous, double current, double previous_h, double h) {
    if (!(previous > 0.0 && current > 0.0 && previous_h > h))
        return std::numeric_limits<double>::quiet_NaN();
    return std::log(previous / current) / std::log(previous_h / h);
}

int ell_for(const Options &options, int degree) {
    const auto it =
        std::find(options.degrees.begin(), options.degrees.end(), degree);
    return options.ell_by_degree[
        static_cast<std::size_t>(it - options.degrees.begin())];
}

Row run_lod_case(
    const Options &options,
    const HelmholtzManufacturedSolution &manufactured,
    int degree,
    int H_level,
    int h_level,
    int ell,
    const std::string &study,
    const FineReference *fine_reference = nullptr,
    const ComplexVector *previous_solution = nullptr,
    ComplexVector *solution_output = nullptr) {
    HelmholtzHpProblemConfig config;
    config.H = H_level;
    config.h = h_level;
    config.ell = ell;
    config.degree = degree;
    config.wavenumber = options.wavenumber;
    config.patch_solver.kind = options.solver == "schur"
        ? HelmholtzPatchSolverKind::DirectSchur
        : HelmholtzPatchSolverKind::DirectSaddle;
    config.patch_solver.fallback_to_direct = false;
    config.corrector_threads = options.threads;
    config.progress_interval = options.progress_interval;
    if (options.progress_interval > 0) {
        config.corrector_progress =
            [degree, H_level, h_level, ell, solver = options.solver]
            (int done, int total) {
                std::cerr << "[hp] p=" << degree
                          << " H=" << H_level
                          << " h=" << h_level
                          << " ell=" << ell
                          << " solver=" << solver
                          << " correctors=" << done << '/' << total
                          << '\n' << std::flush;
            };
    }
    std::cerr << "[hp] start p=" << degree
              << " H=" << H_level << " h=" << h_level
              << " ell=" << ell << " solver=" << options.solver
              << " threads=" << options.threads
              << '\n' << std::flush;
    HelmholtzHpLodModel model = HelmholtzHpLodModel::build(config);
    const ComplexVector load =
        assemble_helmholtz_hp_load(model.fine_space(), manufactured.source);
    ComplexVector owned_fine;
    if (!fine_reference)
        owned_fine = model.solve_fine_reference(load);
    const ComplexVector &fine = fine_reference
        ? fine_reference->solution : owned_fine;
    if (fine.size() != model.fine_space().dof_count())
        throw std::runtime_error("fine reference size mismatch");
    const HelmholtzLodSolution lod = model.solve_load(load);
    if (solution_output)
        *solution_output = lod.fine_values;

    Row row;
    row.study = study;
    row.solver = options.solver;
    row.degree = degree;
    row.H_level = H_level;
    row.h_level = h_level;
    row.ell = ell;
    row.fine_dofs = model.fine_space().dof_count();
    row.H = max_element_diameter(model.problem().coarse);
    row.h = max_element_diameter(model.problem().fine);
    const HelmholtzError exact = compute_helmholtz_hp_error(
        model.fine_space(), lod.fine_values, options.wavenumber,
        manufactured.value, manufactured.gradient);
    row.exact = {exact.energy, exact.l2};
    if (fine_reference) {
        row.fine_exact = fine_reference->exact;
    } else {
        const HelmholtzError fine_exact = compute_helmholtz_hp_error(
            model.fine_space(), fine, options.wavenumber,
            manufactured.value, manufactured.gradient);
        row.fine_exact = {fine_exact.energy, fine_exact.l2};
    }
    row.lod_fine = discrete_error(model.operators(), lod.fine_values, fine);
    if (previous_solution)
        row.delta = discrete_error(
            model.operators(), lod.fine_values, *previous_solution).energy;
    row.petrov_residual = lod.petrov_residual;
    row.corrector_residual = std::max(
        model.corrector_diagnostics().max_primal_residual,
        model.corrector_diagnostics().max_adjoint_residual);
    row.constraint_residual =
        model.corrector_diagnostics().max_constraint_residual;
    row.schur_residual =
        model.corrector_diagnostics().max_schur_residual;
    row.schur_rcond =
        model.corrector_diagnostics().min_schur_rcond;
    row.direct_fallbacks =
        model.corrector_diagnostics().direct_fallback_count;
    row.corrector_threads =
        model.corrector_diagnostics().parallel_threads;
    std::cerr << "[hp] complete p=" << degree
              << " H=" << H_level << " h=" << h_level
              << " ell=" << ell << '\n' << std::flush;
    return row;
}

std::vector<Row> run_fem(
    const Options &options,
    const HelmholtzManufacturedSolution &manufactured) {
    std::vector<Row> rows;
    for (int degree : options.degrees) {
        Row previous;
        bool have_previous = false;
        for (int level : options.fine_levels) {
            TriMesh mesh = refine_mesh_nvb(
                make_helmholtz_unit_square_mesh(), level).mesh;
            HpTriSpace space(mesh, degree);
            const auto operators = assemble_helmholtz_hp_operators(
                space, options.wavenumber);
            const ComplexVector load =
                assemble_helmholtz_hp_load(space, manufactured.source);
            const ComplexVector solution =
                solve_helmholtz_hp_fem(operators, load);
            const HelmholtzError error = compute_helmholtz_hp_error(
                space, solution, options.wavenumber,
                manufactured.value, manufactured.gradient);
            Row row;
            row.study = "fem";
            row.solver = "none";
            row.degree = degree;
            row.h_level = level;
            row.fine_dofs = space.dof_count();
            row.h = max_element_diameter(mesh);
            row.exact = {error.energy, error.l2};
            row.fine_exact = row.exact;
            const ComplexVector residual = operators.system * solution - load;
            row.petrov_residual =
                residual.norm() / std::max(1.0, load.norm());
            if (have_previous) {
                row.energy_rate = rate(
                    previous.exact.energy, row.exact.energy,
                    previous.h, row.h);
                row.l2_rate = rate(
                    previous.exact.l2, row.exact.l2,
                    previous.h, row.h);
            }
            rows.push_back(row);
            if (options.stream) print_row(row);
            previous = row;
            have_previous = true;
        }
    }
    return rows;
}

std::vector<Row> run_ell(
    const Options &options,
    const HelmholtzManufacturedSolution &manufactured) {
    std::vector<Row> rows;
    for (int degree : options.degrees) {
        const FineReference fine_reference = build_fine_reference(
            options, manufactured, degree, options.fine_level);
        ComplexVector previous_solution;
        bool have_previous = false;
        for (int ell : options.ell_levels) {
            ComplexVector current_solution;
            Row row = run_lod_case(
                options, manufactured, degree, options.coarse_level,
                options.fine_level, ell, "ell", &fine_reference,
                have_previous ? &previous_solution : nullptr,
                &current_solution);
            previous_solution = std::move(current_solution);
            have_previous = true;
            rows.push_back(row);
            if (options.stream) print_row(row);
        }
    }
    return rows;
}

std::vector<Row> run_H_like(
    const Options &options,
    const HelmholtzManufacturedSolution &manufactured,
    bool coupled) {
    std::vector<Row> rows;
    for (int degree : options.degrees) {
        FineReference fixed_reference;
        if (!coupled)
            fixed_reference = build_fine_reference(
                options, manufactured, degree, options.fine_level);
        Row previous;
        bool have_previous = false;
        for (int H_level : options.coarse_levels) {
            const int h_level =
                coupled ? H_level + options.gap : options.fine_level;
            FineReference coupled_reference;
            const FineReference *fine_reference = &fixed_reference;
            if (coupled) {
                coupled_reference = build_fine_reference(
                    options, manufactured, degree, h_level);
                fine_reference = &coupled_reference;
            }
            Row row = run_lod_case(
                options, manufactured, degree, H_level, h_level,
                ell_for(options, degree), coupled ? "coupled" : "H",
                fine_reference);
            if (have_previous) {
                row.energy_rate = rate(
                    previous.exact.energy, row.exact.energy,
                    previous.H, row.H);
                row.l2_rate = rate(
                    previous.exact.l2, row.exact.l2,
                    previous.H, row.H);
            }
            rows.push_back(row);
            if (options.stream) print_row(row);
            previous = row;
            have_previous = true;
        }
    }
    return rows;
}

void print_header() {
    std::cout
        << "study,p,H_level,h_level,ell,fine_dofs,H,h,"
           "exact_energy,exact_l2,fine_energy,fine_l2,"
           "lod_fine_energy,lod_fine_l2,delta,rate_energy,rate_l2,"
           "petrov_residual,corrector_residual,constraint_residual,"
           "corrector_threads,solver,schur_residual,schur_rcond,"
           "direct_fallbacks\n"
        << std::flush;
}

void print_row(const Row &row) {
    std::cout << std::setprecision(12)
              << row.study << ',' << row.degree << ','
              << row.H_level << ',' << row.h_level << ',' << row.ell << ','
              << row.fine_dofs << ',' << row.H << ',' << row.h << ','
              << row.exact.energy << ',' << row.exact.l2 << ','
              << row.fine_exact.energy << ',' << row.fine_exact.l2 << ','
              << row.lod_fine.energy << ',' << row.lod_fine.l2 << ','
              << row.delta << ',' << row.energy_rate << ',' << row.l2_rate
              << ',' << row.petrov_residual << ','
              << row.corrector_residual << ','
              << row.constraint_residual << ','
              << row.corrector_threads << ',' << row.solver << ','
              << row.schur_residual << ',' << row.schur_rcond << ','
              << row.direct_fallbacks << '\n'
              << std::flush;
}

void print_rows(const std::vector<Row> &rows) {
    print_header();
    for (const Row &row : rows) print_row(row);
}
void check_rows(const Options &options, const std::vector<Row> &rows) {
    for (const Row &row : rows) {
        if (!(std::isfinite(row.exact.energy)
              && std::isfinite(row.exact.l2)))
            throw std::runtime_error("nonfinite hp convergence error");
        if (row.petrov_residual > 1e-8
            || row.corrector_residual > 1e-8
            || row.constraint_residual > 1e-8
            || row.schur_residual > 1e-8
            || row.direct_fallbacks != 0)
            throw std::runtime_error("hp convergence residual check failed");
    }
    if (options.study == "fem") {
        for (int degree : options.degrees) {
            std::vector<const Row *> selected;
            for (const Row &row : rows)
                if (row.degree == degree) selected.push_back(&row);
            if (selected.size() >= 2
                && !(selected.back()->exact.energy
                     < selected.front()->exact.energy))
                throw std::runtime_error("fine-hp error did not decrease");
        }
    }
}

} // namespace

int main(int argc, char **argv) {
    try {
        const Options options = parse_options(argc, argv);
        const auto manufactured =
            make_polynomial_plane_wave_solution(options.wavenumber);
        std::vector<Row> rows;
        if (options.stream) print_header();
        if (options.study == "fem")
            rows = run_fem(options, manufactured);
        else if (options.study == "ell")
            rows = run_ell(options, manufactured);
        else
            rows = run_H_like(
                options, manufactured, options.study == "coupled");
        if (options.check) check_rows(options, rows);
        if (!options.stream) print_rows(rows);
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "bench_helmholtz_hp_convergence failed: "
                  << error.what() << '\n';
        return 1;
    }
}

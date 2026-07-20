#include "helmholtz/manufactured.h"
#include "helmholtz/operators.h"
#include "helmholtz/two_level_schwarz.h"
#include "solver/right_gmres.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace lod2d::helmholtz;
namespace solver = lod2d::solver;

namespace {

using Clock = std::chrono::steady_clock;

struct Options {
    int H = 5;
    int h = 10;
    int ell = 3;
    int threads = 8;
    int restart = 100;
    int max_iterations = 2000;
    double k = 4.0;
    double tolerance = 1e-10;
    double impedance_beta = 1.0;
    double local_shift_alpha = 0.2;
    double local_tolerance = 1e-10;
    int local_restart = 30;
    int local_max_iterations = 200;
    int local_pre_smooth = 2;
    int local_post_smooth = 2;
    int local_coarse_max_dofs = 256;
    double local_jacobi_weight = 0.6;
    std::string local_solver = "direct";
    std::string local_inverse = "lu";
    std::string solver = "all";
    std::string boundary = "dirichlet";
    std::string extension = "weighted";
    std::string factorization_reuse = "none";
    std::string source = "gaussian";
};

int parse_int(const std::string &text, const char *name) {
    std::size_t consumed = 0;
    const int value = std::stoi(text, &consumed);
    if (consumed != text.size())
        throw std::invalid_argument(std::string("invalid ") + name);
    return value;
}

double parse_double(const std::string &text, const char *name) {
    std::size_t consumed = 0;
    const double value = std::stod(text, &consumed);
    if (consumed != text.size() || !std::isfinite(value))
        throw std::invalid_argument(std::string("invalid ") + name);
    return value;
}

Options parse_options(int argc, char **argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        const auto after = [&](const std::string &prefix) {
            return argument.substr(prefix.size());
        };
        if (argument.rfind("--H=", 0) == 0)
            options.H = parse_int(after("--H="), "H");
        else if (argument.rfind("--h=", 0) == 0)
            options.h = parse_int(after("--h="), "h");
        else if (argument.rfind("--ell=", 0) == 0)
            options.ell = parse_int(after("--ell="), "ell");
        else if (argument.rfind("--k=", 0) == 0)
            options.k = parse_double(after("--k="), "k");
        else if (argument.rfind("--threads=", 0) == 0)
            options.threads = parse_int(after("--threads="), "threads");
        else if (argument.rfind("--restart=", 0) == 0)
            options.restart = parse_int(after("--restart="), "restart");
        else if (argument.rfind("--max-iters=", 0) == 0)
            options.max_iterations =
                parse_int(after("--max-iters="), "max iterations");
        else if (argument.rfind("--tol=", 0) == 0)
            options.tolerance = parse_double(after("--tol="), "tolerance");
        else if (argument.rfind("--solver=", 0) == 0)
            options.solver = after("--solver=");
        else if (argument.rfind("--source=", 0) == 0)
            options.source = after("--source=");
        else if (argument.rfind("--boundary=", 0) == 0)
            options.boundary = after("--boundary=");
        else if (argument.rfind("--extension=", 0) == 0)
            options.extension = after("--extension=");
        else if (argument.rfind("--factorization-reuse=", 0) == 0)
            options.factorization_reuse =
                after("--factorization-reuse=");
        else if (argument.rfind("--impedance-beta=", 0) == 0)
            options.impedance_beta =
                parse_double(after("--impedance-beta="), "impedance beta");
        else if (argument.rfind("--local-solver=", 0) == 0)
            options.local_solver = after("--local-solver=");
        else if (argument.rfind("--local-alpha=", 0) == 0)
            options.local_shift_alpha =
                parse_double(after("--local-alpha="), "local shift alpha");
        else if (argument.rfind("--local-tol=", 0) == 0)
            options.local_tolerance =
                parse_double(after("--local-tol="), "local tolerance");
        else if (argument.rfind("--local-restart=", 0) == 0)
            options.local_restart =
                parse_int(after("--local-restart="), "local restart");
        else if (argument.rfind("--local-max-iters=", 0) == 0)
            options.local_max_iterations =
                parse_int(after("--local-max-iters="), "local max iterations");
        else if (argument.rfind("--local-inverse=", 0) == 0)
            options.local_inverse = after("--local-inverse=");
        else if (argument.rfind("--local-pre-smooth=", 0) == 0)
            options.local_pre_smooth =
                parse_int(after("--local-pre-smooth="), "local pre smooth");
        else if (argument.rfind("--local-post-smooth=", 0) == 0)
            options.local_post_smooth =
                parse_int(after("--local-post-smooth="), "local post smooth");
        else if (argument.rfind("--local-coarse-max=", 0) == 0)
            options.local_coarse_max_dofs =
                parse_int(after("--local-coarse-max="), "local coarse max");
        else if (argument.rfind("--local-omega=", 0) == 0)
            options.local_jacobi_weight =
                parse_double(after("--local-omega="), "local omega");
        else if (argument == "--help") {
            std::cout
                << "Usage: bench_helmholtz_two_level_schwarz [--H=5] "
                   "[--h=10] [--ell=3] [--k=4] [--threads=8] "
                   "[--restart=100] [--max-iters=2000] [--tol=1e-10] "
                   "[--solver=all|identity|local|additive|hybrid] "
                   "[--source=gaussian|manufactured] "
                   "[--boundary=dirichlet|impedance] "
                   "[--extension=weighted|restricted] [--impedance-beta=1] "
                   "[--factorization-reuse=none|identical] "
                   "[--local-solver=direct|shifted-gmres] [--local-alpha=0.2] "
                   "[--local-tol=1e-10] [--local-restart=30] "
                   "[--local-max-iters=200] [--local-inverse=lu|vcycle] "
                   "[--local-pre-smooth=2] [--local-post-smooth=2] "
                   "[--local-coarse-max=256] [--local-omega=0.6]\n";
            std::exit(0);
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.H < 0 || options.h < options.H || options.ell < 0)
        throw std::invalid_argument("require 0 <= H <= h and ell >= 0");
    if (!(options.k > 0.0) || !(options.tolerance > 0.0)
        || !(options.impedance_beta >= 0.0)
        || !(options.local_shift_alpha >= 0.0)
        || !(options.local_tolerance > 0.0))
        throw std::invalid_argument(
            "k and tolerances must be positive; shift parameters nonnegative");
    if (options.threads <= 0 || options.restart <= 0
        || options.max_iterations <= 0 || options.local_restart <= 0
        || options.local_max_iterations <= 0
        || options.local_pre_smooth < 0 || options.local_post_smooth < 0
        || options.local_coarse_max_dofs <= 0
        || !(options.local_jacobi_weight > 0.0))
        throw std::invalid_argument(
            "thread/coarse limits must be positive; smoothing nonnegative");
    if (options.solver != "all" && options.solver != "identity"
        && options.solver != "local" && options.solver != "additive"
        && options.solver != "hybrid")
        throw std::invalid_argument(
            "solver must be all, identity, local, additive, or hybrid");
    if (options.source != "gaussian" && options.source != "manufactured")
        throw std::invalid_argument(
            "source must be gaussian or manufactured");
    if (options.boundary != "dirichlet" && options.boundary != "impedance")
        throw std::invalid_argument(
            "boundary must be dirichlet or impedance");
    if (options.extension != "weighted" && options.extension != "restricted")
        throw std::invalid_argument(
            "extension must be weighted or restricted");
    if (options.factorization_reuse != "none"
        && options.factorization_reuse != "identical")
        throw std::invalid_argument(
            "factorization reuse must be none or identical");
    if (options.local_solver != "direct"
        && options.local_solver != "shifted-gmres")
        throw std::invalid_argument(
            "local solver must be direct or shifted-gmres");
    if (options.local_inverse != "lu"
        && options.local_inverse != "vcycle")
        throw std::invalid_argument(
            "local inverse must be lu or vcycle");
    if (options.factorization_reuse == "identical"
        && options.local_solver != "direct")
        throw std::invalid_argument(
            "identical factorization reuse requires --local-solver=direct");
    if (options.local_inverse == "vcycle"
        && options.local_solver != "shifted-gmres")
        throw std::invalid_argument(
            "local vcycle requires --local-solver=shifted-gmres");
    if (options.local_inverse == "vcycle"
        && options.boundary != "dirichlet")
        throw std::invalid_argument(
            "local vcycle currently requires --boundary=dirichlet");
    return options;
}

double milliseconds(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

double energy_norm(
    const HelmholtzOperators &operators,
    const ComplexVector &vector) {
    const ComplexVector weighted =
        operators.stiffness.cast<Complex>() * vector
        + operators.wavenumber * operators.wavenumber
            * (operators.mass.cast<Complex>() * vector);
    return std::sqrt(std::max(0.0, std::real(vector.dot(weighted))));
}

struct RunResult {
    std::string name;
    bool converged = false;
    double solve_ms = 0.0;
    int iterations = 0;
    int restarts = 0;
    int preconditioner_applications = 0;
    double true_residual = 0.0;
    double energy_error = 0.0;
    double fine_error = 0.0;
    double petrov_residual = 0.0;
    double exact_energy_error = std::numeric_limits<double>::quiet_NaN();
    double exact_l2_error = std::numeric_limits<double>::quiet_NaN();
    double exact_energy_relative = std::numeric_limits<double>::quiet_NaN();
    double exact_l2_relative = std::numeric_limits<double>::quiet_NaN();
    std::string message;
};

template <class Preconditioner>
RunResult run_gmres(
    const std::string &name,
    const HelmholtzLodModel &model,
    const ComplexVector &right_hand_side,
    const ComplexVector &reference,
    const Options &options,
    const HelmholtzManufacturedSolution *manufactured,
    Preconditioner preconditioner) {
    const ComplexSparseMatrix &matrix = model.operators().system;
    solver::RightGmresConfig config;
    config.restart = options.restart;
    config.max_iterations = options.max_iterations;
    config.relative_tolerance = options.tolerance;

    int preconditioner_applications = 0;
    const auto apply_operator = [&](const solver::ComplexVector &vector) {
        return solver::ComplexVector(matrix * vector);
    };
    const auto apply_preconditioner = [&](const solver::ComplexVector &vector) {
        ++preconditioner_applications;
        return solver::ComplexVector(preconditioner(vector));
    };

    const auto start = Clock::now();
    const bool flexible =
        options.local_solver == "shifted-gmres";
    const solver::RightGmresResult result = flexible
        ? solver::solve_right_preconditioned_fgmres(
            matrix.rows(), apply_operator, apply_preconditioner,
            right_hand_side, config)
        : solver::solve_right_preconditioned_gmres(
            matrix.rows(), apply_operator, apply_preconditioner,
            right_hand_side, config);
    const auto end = Clock::now();

    const double minimum_scale = std::numeric_limits<double>::min();
    RunResult report;
    report.name = name;
    report.converged = result.converged;
    report.solve_ms = milliseconds(start, end);
    report.iterations = result.iterations;
    report.restarts = result.restarts;
    report.preconditioner_applications = preconditioner_applications;
    report.true_residual =
        (right_hand_side - matrix * result.solution).norm()
        / std::max(minimum_scale, right_hand_side.norm());
    report.energy_error =
        energy_norm(model.operators(), result.solution - reference)
        / std::max(minimum_scale, energy_norm(model.operators(), reference));
    report.fine_error = (result.solution - reference).norm()
        / std::max(minimum_scale, reference.norm());
    const ComplexVector coarse_rhs =
        model.test_basis().adjoint() * right_hand_side;
    const ComplexVector petrov = model.test_basis().adjoint()
        * (matrix * result.solution - right_hand_side);
    report.petrov_residual =
        petrov.norm() / std::max(1.0, coarse_rhs.norm());
    if (manufactured) {
        const HelmholtzError exact = compute_helmholtz_error(
            model.problem().fine,
            result.solution,
            model.operators().wavenumber,
            manufactured->value,
            manufactured->gradient);
        report.exact_energy_error = exact.energy;
        report.exact_l2_error = exact.l2;
        const HelmholtzError exact_norm = compute_helmholtz_error(
            model.problem().fine,
            ComplexVector::Zero(model.operators().system.rows()),
            model.operators().wavenumber,
            manufactured->value,
            manufactured->gradient);
        report.exact_energy_relative =
            exact.energy / std::max(minimum_scale, exact_norm.energy);
        report.exact_l2_relative =
            exact.l2 / std::max(minimum_scale, exact_norm.l2);
    }
    report.message = result.message;
    return report;
}

void print_result(const RunResult &result) {
    std::cout << std::left << std::setw(9) << result.name
              << " status=" << std::setw(10)
              << (result.converged ? "converged" : "failed")
              << " solve_ms=" << std::setw(11) << result.solve_ms
              << " iterations=" << std::setw(5) << result.iterations
              << " restarts=" << std::setw(4) << result.restarts
              << " prec_apps=" << std::setw(5)
              << result.preconditioner_applications
              << " residual=" << std::setw(12) << result.true_residual
              << " energy_rel=" << std::setw(12) << result.energy_error
              << " fine_rel=" << std::setw(12) << result.fine_error
              << " petrov=" << std::setw(12) << result.petrov_residual
              << " exact_E=" << std::setw(12) << result.exact_energy_error
              << " exact_L2=" << std::setw(12) << result.exact_l2_error
              << " exact_E_rel=" << std::setw(12)
              << result.exact_energy_relative
              << " exact_L2_rel=" << result.exact_l2_relative;
    if (!result.converged) std::cout << " message=\"" << result.message << '"';
    std::cout << '\n';
}

} // namespace

int main(int argc, char **argv) {
    try {
        const auto benchmark_start = Clock::now();
        const Options options = parse_options(argc, argv);
#ifdef _OPENMP
        omp_set_num_threads(options.threads);
#else
        (void)options.threads;
#endif

        HelmholtzProblemConfig config;
        config.H = options.H;
        config.h = options.h;
        config.ell = options.ell;
        config.wavenumber = options.k;
        config.patch_solver.kind = HelmholtzPatchSolverKind::DirectSaddle;

        const auto build_start = Clock::now();
        HelmholtzLodModel model = HelmholtzLodModel::build(config);
        const auto build_end = Clock::now();

        HelmholtzManufacturedSolution manufactured;
        const HelmholtzManufacturedSolution *manufactured_ptr = nullptr;
        ComplexFunction source;
        if (options.source == "manufactured") {
            manufactured =
                make_polynomial_plane_wave_solution(options.k);
            manufactured_ptr = &manufactured;
            source = manufactured.source;
        } else {
            source = [](const lod2d::Point2 &point) {
                const double dx = point.x() - 0.35;
                const double dy = point.y() - 0.55;
                return Complex(std::exp(-40.0 * (dx * dx + dy * dy)), 0.0);
            };
        }
        const ComplexVector right_hand_side =
            assemble_helmholtz_load(model.problem().fine, source);

        const auto direct_start = Clock::now();
        const ComplexVector reference =
            model.solve_fine_reference(right_hand_side);
        const auto direct_end = Clock::now();

        const auto lod_start = Clock::now();
        const HelmholtzLodSolution lod_baseline =
            model.solve_load(right_hand_side);
        const auto lod_end = Clock::now();
        const double minimum_scale = std::numeric_limits<double>::min();
        const double lod_energy_relative = energy_norm(
            model.operators(), lod_baseline.fine_values - reference)
            / std::max(
                minimum_scale,
                energy_norm(model.operators(), reference));
        const double lod_fine_relative =
            (lod_baseline.fine_values - reference).norm()
            / std::max(minimum_scale, reference.norm());
        HelmholtzError exact_norm;
        HelmholtzError fem_exact;
        HelmholtzError lod_exact;
        if (manufactured_ptr) {
            exact_norm = compute_helmholtz_error(
                model.problem().fine,
                ComplexVector::Zero(model.operators().system.rows()),
                options.k, manufactured.value, manufactured.gradient);
            fem_exact = compute_helmholtz_error(
                model.problem().fine, reference, options.k,
                manufactured.value, manufactured.gradient);
            lod_exact = compute_helmholtz_error(
                model.problem().fine, lod_baseline.fine_values, options.k,
                manufactured.value, manufactured.gradient);
        } else {
            exact_norm.energy = exact_norm.l2 =
                std::numeric_limits<double>::quiet_NaN();
            fem_exact.energy = fem_exact.l2 =
                std::numeric_limits<double>::quiet_NaN();
            lod_exact.energy = lod_exact.l2 =
                std::numeric_limits<double>::quiet_NaN();
        }

        HelmholtzTwoLevelSchwarzConfig schwarz_config;
        schwarz_config.artificial_boundary =
            options.boundary == "impedance"
                ? HelmholtzSchwarzArtificialBoundary::Impedance
                : HelmholtzSchwarzArtificialBoundary::HomogeneousDirichlet;
        schwarz_config.artificial_impedance_beta = options.impedance_beta;
        schwarz_config.extension = options.extension == "restricted"
            ? HelmholtzSchwarzExtension::RestrictedCore
            : HelmholtzSchwarzExtension::WeightedOverlap;
        schwarz_config.factorization_reuse =
            options.factorization_reuse == "identical"
                ? HelmholtzSchwarzFactorizationReuse::IdenticalMatrix
                : HelmholtzSchwarzFactorizationReuse::None;
        schwarz_config.local_solver.kind =
            options.local_solver == "shifted-gmres"
                ? HelmholtzSchwarzLocalSolverKind::ShiftedGmres
                : HelmholtzSchwarzLocalSolverKind::SparseLu;
        schwarz_config.local_solver.shifted_inverse =
            options.local_inverse == "vcycle"
                ? HelmholtzSchwarzShiftedInverseKind::GeometricVcycle
                : HelmholtzSchwarzShiftedInverseKind::SparseLu;
        schwarz_config.local_solver.shift_alpha = options.local_shift_alpha;
        schwarz_config.local_solver.vcycle_pre_smooth =
            options.local_pre_smooth;
        schwarz_config.local_solver.vcycle_post_smooth =
            options.local_post_smooth;
        schwarz_config.local_solver.vcycle_coarse_max_dofs =
            options.local_coarse_max_dofs;
        schwarz_config.local_solver.vcycle_jacobi_weight =
            options.local_jacobi_weight;
        schwarz_config.local_solver.gmres.relative_tolerance =
            options.local_tolerance;
        schwarz_config.local_solver.gmres.restart = options.local_restart;
        schwarz_config.local_solver.gmres.max_iterations =
            options.local_max_iterations;
        const auto setup_start = Clock::now();
        HelmholtzTwoLevelSchwarzPreconditioner preconditioner(
            model, schwarz_config);
        const auto setup_end = Clock::now();

        std::cout << "Helmholtz two-level Schwarz benchmark: H=" << options.H
                  << " h=" << options.h << " ell=" << options.ell
                  << " k=" << options.k << " threads=" << options.threads
                  << " source=" << options.source
                  << " boundary=" << options.boundary
                  << " extension=" << options.extension
                  << " factorization_reuse="
                  << options.factorization_reuse
                  << " impedance_beta=" << options.impedance_beta
                  << " local_solver=" << options.local_solver
                  << " local_inverse=" << options.local_inverse
                  << " local_alpha=" << options.local_shift_alpha
                  << "\ncoarse_dofs=" << model.coarse_operator().rows()
                  << " fine_dofs=" << model.operators().system.rows()
                  << " coarse_meshwidth="
                  << max_element_diameter(model.problem().coarse)
                  << " fine_meshwidth="
                  << max_element_diameter(model.problem().fine)
                  << " kH=" << options.k
                    * max_element_diameter(model.problem().coarse)
                  << " build_ms=" << milliseconds(build_start, build_end)
                  << " fine_sparse_lu_ms="
                  << milliseconds(direct_start, direct_end)
                  << " lod_solve_ms=" << milliseconds(lod_start, lod_end)
                  << " schwarz_setup_ms="
                  << milliseconds(setup_start, setup_end) << '\n'
                  << "LOD baseline petrov="
                  << lod_baseline.petrov_residual
                  << " energy_rel=" << lod_energy_relative
                  << " fine_rel=" << lod_fine_relative
                  << " fem_exact_E=" << fem_exact.energy
                  << " fem_exact_L2=" << fem_exact.l2
                  << " fem_exact_E_rel="
                  << fem_exact.energy / exact_norm.energy
                  << " fem_exact_L2_rel="
                  << fem_exact.l2 / exact_norm.l2
                  << " lod_exact_E=" << lod_exact.energy
                  << " lod_exact_L2=" << lod_exact.l2
                  << " lod_exact_E_rel="
                  << lod_exact.energy / exact_norm.energy
                  << " lod_exact_L2_rel="
                  << lod_exact.l2 / exact_norm.l2
                  << " corrector_residual="
                  << model.correctors().diagnostics.max_primal_residual
                  << " constraint_residual="
                  << model.correctors().diagnostics.max_constraint_residual
                  << '\n';
        const auto &diagnostics = preconditioner.diagnostics();
        std::cout << "Schwarz subdomains=" << diagnostics.subdomains
                  << " local_dofs=" << diagnostics.min_local_dofs
                  << ".." << diagnostics.max_local_dofs
                  << " owned_dofs=" << diagnostics.min_owned_dofs
                  << ".." << diagnostics.max_owned_dofs
                  << " uncovered=" << diagnostics.uncovered_dofs
                  << " partition_error="
                  << diagnostics.partition_unity_error
                  << " artificial_edges="
                  << diagnostics.artificial_boundary_edges
                  << " physical_edges="
                  << diagnostics.physical_boundary_edges
                  << " solver_groups="
                  << diagnostics.local_solver_groups
                  << " reused_factorizations="
                  << diagnostics.reused_factorizations
                  << " max_reuse_group="
                  << diagnostics.max_reuse_group << '\n';

        const auto should_run = [&](const char *name) {
            return options.solver == "all" || options.solver == name;
        };
        if (should_run("identity"))
            print_result(run_gmres(
                "identity", model, right_hand_side, reference, options,
                manufactured_ptr,
                [](const ComplexVector &vector) { return vector; }));
        if (should_run("local"))
            print_result(run_gmres(
                "local", model, right_hand_side, reference, options,
                manufactured_ptr,
                [&](const ComplexVector &vector) {
                    return preconditioner.apply_local(vector);
                }));
        if (should_run("additive"))
            print_result(run_gmres(
                "additive", model, right_hand_side, reference, options,
                manufactured_ptr,
                [&](const ComplexVector &vector) {
                    return preconditioner.apply(
                        vector, HelmholtzTwoLevelSchwarzMode::Additive);
                }));
        if (should_run("hybrid"))
            print_result(run_gmres(
                "hybrid", model, right_hand_side, reference, options,
                manufactured_ptr,
                [&](const ComplexVector &vector) {
                    return preconditioner.apply(
                        vector, HelmholtzTwoLevelSchwarzMode::Hybrid);
                }));
        const HelmholtzSchwarzLocalSolverDiagnostics local_diagnostics =
            preconditioner.local_solver_diagnostics();
        const double average_local_iterations =
            local_diagnostics.solve_calls > 0
                ? static_cast<double>(local_diagnostics.total_iterations)
                    / static_cast<double>(local_diagnostics.solve_calls)
                : 0.0;
        std::cout << "Local solver calls=" << local_diagnostics.solve_calls
                  << " avg_iterations=" << average_local_iterations
                  << " max_iterations=" << local_diagnostics.max_iterations
                  << " restarts=" << local_diagnostics.total_restarts
                  << " max_true_residual="
                  << local_diagnostics.max_relative_residual
                  << " vcycle_levels="
                  << local_diagnostics.max_vcycle_levels
                  << " vcycle_coarse_dofs="
                  << local_diagnostics.min_vcycle_coarse_dofs
                  << " vcycle_finest_dofs="
                  << local_diagnostics.max_vcycle_finest_dofs << '\n';
        std::cout << "Benchmark total_ms="
                  << milliseconds(benchmark_start, Clock::now()) << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "bench_helmholtz_two_level_schwarz failed: "
                  << error.what() << '\n';
        return 1;
    }
}

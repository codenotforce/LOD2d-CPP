#include "helmholtz/model.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace lod2d::helmholtz;
#ifdef _OPENMP
#include <omp.h>
#endif


namespace {

struct Options {
    int H = 1;
    int h = 5;
    int ell = 1;
    int threads = 1;
    int restart = 30;
    int max_iterations = 200;
    double k = 4.0;
    double alpha = 0.2;
    double tolerance = 1e-10;
    double jacobi_weight = 0.6;
    int pre_smooth = 2;
    int post_smooth = 2;
    int coarse_max_dofs = 200;
    std::string inverse = "lu";
    std::string solver = "all";
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
        else if (argument.rfind("--solver=", 0) == 0)
            options.solver = after("--solver=");
        else if (argument.rfind("--alpha=", 0) == 0)
            options.alpha = parse_double(after("--alpha="), "alpha");
        else if (argument.rfind("--inverse=", 0) == 0)
            options.inverse = after("--inverse=");
        else if (argument.rfind("--pre-smooth=", 0) == 0)
            options.pre_smooth = parse_int(after("--pre-smooth="), "pre smooth");
        else if (argument.rfind("--post-smooth=", 0) == 0)
            options.post_smooth = parse_int(after("--post-smooth="), "post smooth");
        else if (argument.rfind("--coarse-max=", 0) == 0)
            options.coarse_max_dofs =
                parse_int(after("--coarse-max="), "coarse maximum DOFs");
        else if (argument.rfind("--omega=", 0) == 0)
            options.jacobi_weight = parse_double(after("--omega="), "omega");
        else if (argument.rfind("--tol=", 0) == 0)
            options.tolerance = parse_double(after("--tol="), "tolerance");
        else if (argument.rfind("--restart=", 0) == 0)
            options.restart = parse_int(after("--restart="), "restart");
        else if (argument.rfind("--max-iters=", 0) == 0)
            options.max_iterations = parse_int(after("--max-iters="), "max iterations");
        else if (argument == "--help") {
            std::cout
                << "Usage: bench_helmholtz_patch_solver [--H=1] [--h=5] "
                   "[--ell=1] [--k=4] [--threads=1] "
                   "[--solver=saddle|schur|gmres|all] [--alpha=0.2] "
                   "[--inverse=none|lu|vcycle] [--pre-smooth=2] [--post-smooth=2] "
                   "[--coarse-max=200] [--omega=0.6] [--tol=1e-10] "
                   "[--restart=30] [--max-iters=200]\n";
            std::exit(0);
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.H < 0 || options.h < options.H || options.ell < 0)
        throw std::invalid_argument("require 0 <= H <= h and ell >= 0");
    if (!(options.k > 0.0) || !(options.alpha >= 0.0)
        || !(options.tolerance >= 0.0))
        throw std::invalid_argument("k must be positive; alpha and tol nonnegative");
    if (options.threads <= 0 || options.restart <= 0 || options.max_iterations <= 0)
        throw std::invalid_argument("thread and iteration counts must be positive");
    if (options.pre_smooth < 0 || options.post_smooth < 0
        || options.coarse_max_dofs <= 0 || !(options.jacobi_weight > 0.0))
        throw std::invalid_argument(
            "smoothing counts must be nonnegative; coarse-max and omega positive");
    if (options.inverse != "none" && options.inverse != "lu" && options.inverse != "vcycle")
        throw std::invalid_argument("inverse must be none, lu, or vcycle");
    if (options.solver != "saddle" && options.solver != "schur"
        && options.solver != "gmres" && options.solver != "all")
        throw std::invalid_argument("solver must be saddle, schur, gmres, or all");
    return options;
}

struct RunResult {
    std::string name;
    double corrector_ms = 0.0;
    double total_ms = 0.0;
    ComplexMatrix basis;
    HelmholtzCorrectorDiagnostics diagnostics;
};

RunResult run(
    const Options &options,
    HelmholtzPatchSolverKind kind,
    const std::string &name) {
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
    config.patch_solver.kind = kind;
    config.patch_solver.shifted.alpha = options.alpha;
    config.patch_solver.gmres.relative_tolerance = options.tolerance;
    config.patch_solver.gmres.restart = options.restart;
    config.patch_solver.shifted.inverse = options.inverse == "vcycle"
        ? HelmholtzShiftedInverseKind::GeometricVcycle
        : options.inverse == "none"
            ? HelmholtzShiftedInverseKind::Identity
            : HelmholtzShiftedInverseKind::SparseLu;
    config.patch_solver.shifted.pre_smooth = options.pre_smooth;
    config.patch_solver.shifted.post_smooth = options.post_smooth;
    config.patch_solver.shifted.coarse_max_dofs = options.coarse_max_dofs;
    config.patch_solver.shifted.jacobi_weight = options.jacobi_weight;
    config.patch_solver.gmres.max_iterations = options.max_iterations;

    HelmholtzLodModel model = HelmholtzLodModel::build(config);
    RunResult result;
    result.name = name;
    result.corrector_ms = model.build_timings().correctors_ms;
    result.total_ms = model.build_timings().total_ms;
    result.basis = ComplexMatrix(model.corrected_trial_basis());
    result.diagnostics = model.correctors().diagnostics;
    return result;
}

void print_result(const RunResult &result, const ComplexMatrix *reference) {
    const auto &d = result.diagnostics;
    double difference = 0.0;
    if (reference)
        difference = (result.basis - *reference).norm()
                   / std::max(1.0, reference->norm());
    std::cout << std::left << std::setw(8) << result.name
              << " corrector_ms=" << std::setw(11) << result.corrector_ms
              << " total_ms=" << std::setw(11) << result.total_ms
              << " basis_rel=" << std::setw(12) << difference
              << " primal=" << std::setw(12) << d.max_primal_residual
              << " constraint=" << std::setw(12) << d.max_constraint_residual;
    if (d.gmres_right_hand_sides > 0) {
        const double average = static_cast<double>(d.gmres_iterations)
                             / d.gmres_right_hand_sides;
        std::cout << " gmres_rhs=" << d.gmres_right_hand_sides
                  << " avg_it=" << average
                  << " max_it=" << d.gmres_max_iterations
                  << " true_res=" << d.max_gmres_relative_residual;
    }
    if (d.max_vcycle_levels > 0) {
        std::cout << " vcycle_rho=" << d.max_vcycle_relative_residual
                  << " levels=" << d.max_vcycle_levels
                  << " coarse=" << d.max_vcycle_coarse_dofs
                  << " fine=" << d.max_vcycle_finest_dofs;
    }
    std::cout << " fallback=" << d.direct_fallbacks << '\n';
}

} // namespace

int main(int argc, char **argv) {
    try {
        const Options options = parse_options(argc, argv);
        std::cout << "Helmholtz patch solver benchmark: H=" << options.H
                  << " h=" << options.h << " ell=" << options.ell
                  << " k=" << options.k << " threads=" << options.threads
                  << " alpha=" << options.alpha
                  << " inverse=" << options.inverse << '\n';

        if (options.solver == "all") {
            const RunResult saddle =
                run(options, HelmholtzPatchSolverKind::DirectSaddle, "saddle");
            print_result(saddle, &saddle.basis);
            const RunResult schur =
                run(options, HelmholtzPatchSolverKind::DirectSchur, "schur");
            print_result(schur, &saddle.basis);
            const RunResult gmres =
                run(options, HelmholtzPatchSolverKind::ShiftedGmres, "gmres");
            print_result(gmres, &saddle.basis);
        } else {
            HelmholtzPatchSolverKind kind = HelmholtzPatchSolverKind::DirectSaddle;
            if (options.solver == "schur") kind = HelmholtzPatchSolverKind::DirectSchur;
            if (options.solver == "gmres") kind = HelmholtzPatchSolverKind::ShiftedGmres;
            print_result(run(options, kind, options.solver), nullptr);
        }
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "bench_helmholtz_patch_solver failed: "
                  << error.what() << '\n';
        return 1;
    }
}

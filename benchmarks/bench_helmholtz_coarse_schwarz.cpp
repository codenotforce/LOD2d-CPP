#include "helmholtz/coarse_schwarz.h"
#include "helmholtz/model.h"
#include "helmholtz/operators.h"
#include "solver/right_gmres.h"

#include <Eigen/SparseLU>
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

namespace {
namespace solver = lod2d::solver;

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
        else if (argument == "--help") {
            std::cout
                << "Usage: bench_helmholtz_coarse_schwarz [--H=5] [--h=10] "
                   "[--ell=3] [--k=4] [--threads=8] [--restart=100] "
                   "[--max-iters=2000] [--tol=1e-10]\n";
            std::exit(0);
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.H < 0 || options.h < options.H || options.ell < 0)
        throw std::invalid_argument("require 0 <= H <= h and ell >= 0");
    if (!(options.k > 0.0) || !(options.tolerance > 0.0))
        throw std::invalid_argument("k and tolerance must be positive");
    if (options.threads <= 0 || options.restart <= 0
        || options.max_iterations <= 0)
        throw std::invalid_argument("thread and iteration counts must be positive");
    return options;
}

double milliseconds(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

struct IterativeResult {
    std::string name;
    double solve_ms = 0.0;
    int iterations = 0;
    int restarts = 0;
    int preconditioner_applications = 0;
    double true_residual = 0.0;
    double coefficient_error = 0.0;
    double fine_error = 0.0;
};

template <class Preconditioner>
IterativeResult run_gmres(
    const std::string &name,
    const ComplexSparseMatrix &matrix,
    const ComplexVector &right_hand_side,
    const ComplexVector &reference,
    const ComplexSparseMatrix &trial_basis,
    const Options &options,
    Preconditioner preconditioner) {
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
    const solver::RightGmresResult result =
        solver::solve_right_preconditioned_gmres(
            matrix.rows(), apply_operator, apply_preconditioner,
            right_hand_side, config);
    const auto end = Clock::now();
    if (!result.converged)
        throw std::runtime_error(name + " GMRES failed: " + result.message);

    const double minimum_scale = std::numeric_limits<double>::min();
    const double rhs_scale =
        std::max(minimum_scale, right_hand_side.norm());
    const double reference_scale =
        std::max(minimum_scale, reference.norm());
    const ComplexVector reference_fine = trial_basis * reference;
    const ComplexVector computed_fine = trial_basis * result.solution;

    IterativeResult report;
    report.name = name;
    report.solve_ms = milliseconds(start, end);
    report.iterations = result.iterations;
    report.restarts = result.restarts;
    report.preconditioner_applications = preconditioner_applications;
    report.true_residual =
        (right_hand_side - matrix * result.solution).norm() / rhs_scale;
    report.coefficient_error =
        (result.solution - reference).norm() / reference_scale;
    report.fine_error =
        (computed_fine - reference_fine).norm()
        / std::max(minimum_scale, reference_fine.norm());
    return report;
}

void print_result(const IterativeResult &result) {
    std::cout << std::left << std::setw(9) << result.name
              << " solve_ms=" << std::setw(11) << result.solve_ms
              << " iterations=" << std::setw(5) << result.iterations
              << " restarts=" << std::setw(4) << result.restarts
              << " prec_apps=" << std::setw(5)
              << result.preconditioner_applications
              << " residual=" << std::setw(12) << result.true_residual
              << " coeff_rel=" << std::setw(12) << result.coefficient_error
              << " fine_rel=" << result.fine_error << '\n';
}

} // namespace

int main(int argc, char **argv) {
    try {
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

        const ComplexFunction source = [](const lod2d::Point2 &point) {
            const double x = point.x();
            const double y = point.y();
            return Complex(std::exp(-40.0 * (
                (x - 0.35) * (x - 0.35) + (y - 0.55) * (y - 0.55))), 0.0);
        };
        const ComplexVector fine_load =
            assemble_helmholtz_load(model.problem().fine, source);
        const ComplexVector right_hand_side =
            model.test_basis().adjoint() * fine_load;
        const ComplexSparseMatrix &matrix = model.coarse_operator();

        const auto direct_start = Clock::now();
        Eigen::SparseLU<ComplexSparseMatrix> direct;
        direct.analyzePattern(matrix);
        direct.factorize(matrix);
        if (direct.info() != Eigen::Success)
            throw std::runtime_error("coarse LOD SparseLU factorization failed");
        const ComplexVector reference = direct.solve(right_hand_side);
        if (direct.info() != Eigen::Success || !reference.allFinite())
            throw std::runtime_error("coarse LOD SparseLU solve failed");
        const auto direct_end = Clock::now();

        const auto ras_start = Clock::now();
        HelmholtzCoarseRasPreconditioner ras(
            matrix, model.problem().coarse, model.problem().patches);
        const auto ras_end = Clock::now();

        ComplexVector inverse_diagonal(matrix.rows());
        for (int index = 0; index < matrix.rows(); ++index) {
            const Complex diagonal = matrix.coeff(index, index);
            if (std::abs(diagonal) <= 1e-14)
                throw std::runtime_error("coarse LOD matrix has a zero diagonal");
            inverse_diagonal(index) = Complex(1.0, 0.0) / diagonal;
        }

        std::cout << "Helmholtz coarse LOD Schwarz benchmark: H=" << options.H
                  << " h=" << options.h << " ell=" << options.ell
                  << " k=" << options.k << " threads=" << options.threads
                  << "\ncoarse_dofs=" << matrix.rows()
                  << " fine_dofs=" << model.problem().fine.nodes.size()
                  << " build_ms=" << milliseconds(build_start, build_end)
                  << " direct_ms=" << milliseconds(direct_start, direct_end)
                  << " ras_setup_ms=" << milliseconds(ras_start, ras_end)
                  << '\n';
        const auto &diagnostics = ras.diagnostics();
        std::cout << "RAS subdomains=" << diagnostics.subdomains
                  << " local_dofs=" << diagnostics.min_local_dofs
                  << ".." << diagnostics.max_local_dofs
                  << " partition_error="
                  << diagnostics.partition_unity_error << '\n';

        print_result(run_gmres(
            "identity", matrix, right_hand_side, reference,
            model.trial_basis(), options,
            [](const ComplexVector &vector) { return vector; }));
        print_result(run_gmres(
            "jacobi", matrix, right_hand_side, reference,
            model.trial_basis(), options,
            [&](const ComplexVector &vector) {
                return ComplexVector(inverse_diagonal.array() * vector.array());
            }));
        print_result(run_gmres(
            "ras", matrix, right_hand_side, reference,
            model.trial_basis(), options,
            [&](const ComplexVector &vector) { return ras.apply(vector); }));
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "bench_helmholtz_coarse_schwarz failed: "
                  << error.what() << '\n';
        return 1;
    }
}

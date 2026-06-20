/// Benchmark repeated LOD solves for the same A/mesh/H/h/ell and different RHS.
#include "lod/lod_model.h"
#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace lod2d;
namespace chr = std::chrono;

namespace {
struct Options {
    CorrectorSolver solver = CorrectorSolver::Cholmod;
    int rhs_count = 8;
    int threads = -1;
};

const char *solver_name(CorrectorSolver s) {
    if (s == CorrectorSolver::Cholmod) return "cholmod";
    if (s == CorrectorSolver::CholmodCached) return "cholmod_cached";
    return "eigen";
}

Options parse_options(int argc, char **argv) {
    Options opt;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.rfind("--solver=", 0) == 0) {
            std::string value = arg.substr(std::string("--solver=").size());
            if (value == "eigen") opt.solver = CorrectorSolver::EigenLLT;
            else if (value == "cholmod") opt.solver = CorrectorSolver::Cholmod;
            else if (value == "cholmod_cached") opt.solver = CorrectorSolver::CholmodCached;
            else if (value == "auto") opt.solver = CorrectorSolver::Cholmod;
            else throw std::invalid_argument("unknown solver: " + value);
        } else if (arg.rfind("--rhs=", 0) == 0) {
            opt.rhs_count = std::stoi(arg.substr(std::string("--rhs=").size()));
        } else if (arg.rfind("--threads=", 0) == 0) {
            std::string value = arg.substr(std::string("--threads=").size());
            if (value == "auto") opt.threads = -1;
            else if (value == "env") opt.threads = 0;
            else opt.threads = std::stoi(value);
        } else {
            throw std::invalid_argument("usage: bench_reuse_rhs [--solver=eigen|cholmod|cholmod_cached|auto] [--rhs=N] [--threads=auto|env|N]");
        }
    }
    if (opt.rhs_count <= 0) throw std::invalid_argument("--rhs must be positive");
    return opt;
}

#ifdef _OPENMP
void apply_thread_option(const Options &opt) {
    if (opt.threads > 0) {
        omp_set_num_threads(opt.threads);
        return;
    }
    if (opt.threads == 0) return;
    if ((opt.solver == CorrectorSolver::Cholmod || opt.solver == CorrectorSolver::CholmodCached) && omp_get_max_threads() > 8)
        omp_set_num_threads(8);
}
#else
void apply_thread_option(const Options &) {}
#endif

std::vector<double> read_coefficients(const std::string &path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("failed to open coefficient file: " + path);
    int n = 0;
    in >> n;
    std::vector<double> Ah(n);
    for (int i = 0; i < n; ++i) in >> Ah[i];
    if (!in) throw std::runtime_error("failed to read coefficient file: " + path);
    return Ah;
}

Eigen::VectorXd make_rhs(const TriMesh &coarse, int m) {
    Eigen::VectorXd f(coarse.nodes.size());
    const double a = static_cast<double>(m + 1);
    const double b = static_cast<double>(m + 2);
    for (size_t i = 0; i < coarse.nodes.size(); ++i) {
        const double x = coarse.nodes[i][0];
        const double y = coarse.nodes[i][1];
        f(static_cast<int>(i)) = 1.0 + 0.2 * std::sin(a * M_PI * x) * std::cos(b * M_PI * y);
    }
    return f;
}

} // namespace

int main(int argc, char **argv) {
    Options opt;
    try {
        opt = parse_options(argc, argv);
    } catch (const std::exception &e) {
        std::cerr << e.what() << "\n";
        return 2;
    }
    apply_thread_option(opt);

    LodProblemConfig config;
    config.H = 4;
    config.h = 10;
    config.ell = 2;
    config.d = 2;
    config.solver = opt.solver;
    config.initial_mesh = make_unit_square_mesh();

    std::cout << "=== Reusable LOD RHS benchmark H=" << config.H << " h=" << config.h
              << " solver=" << solver_name(opt.solver)
              << " rhs=" << opt.rhs_count << " ===\n";

    try {
        std::vector<double> Ah = read_coefficients("benchmarks/data_H4h10.txt");

        auto t0 = chr::high_resolution_clock::now();
        LodModel model = LodModel::build(config, Ah);
        auto t_setup = chr::high_resolution_clock::now();
        Ah.clear();
        Ah.shrink_to_fit();

        const LodProblemData &problem = model.problem();
        double checksum = 0.0;
        auto t_rhs0 = chr::high_resolution_clock::now();
        for (int r = 0; r < opt.rhs_count; ++r) {
            Eigen::VectorXd f = make_rhs(problem.coarse, r);
            LodReuseSolution sol = model.solve_from_coarse_values(f);
            checksum += sol.uHms.squaredNorm();
        }
        auto t_rhs1 = chr::high_resolution_clock::now();

        const double setup_ms = chr::duration<double, std::milli>(t_setup - t0).count();
        const double rhs_ms = chr::duration<double, std::milli>(t_rhs1 - t_rhs0).count();

        std::cout << "Coarse:" << problem.NH << "v " << problem.NTH << "t  Fine:"
                  << problem.Nh << "v " << problem.NTh << "t\n";
        std::cout << "Reusable setup: " << setup_ms << " ms";
#ifdef _OPENMP
        std::cout << " (" << omp_get_max_threads() << " threads)";
#endif
        std::cout << "\n";
        std::cout << "Repeated RHS total: " << rhs_ms << " ms\n";
        std::cout << "Repeated RHS avg: " << rhs_ms / opt.rhs_count << " ms\n";
        std::cout << "Checksum: " << checksum << "\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << "\n";
        return 1;
    }
    return 0;
}

#include "helmholtz/adaptive/driver.h"
#include "helmholtz/manufactured.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace lod2d;
using namespace lod2d::helmholtz;
using namespace lod2d::helmholtz::adaptive;

namespace {

struct Options {
    double wavenumber = 8.0;
    int initial_H = 2;
    int fine_h = 7;
    int ell = 4;
    int iterations = 6;
    int max_dofs = 20000;
    int threads = 0;
    double theta = 0.5;
    double tolerance = 0.0;
    double q_limit = 0.5;
    int dual_patch_layers = 1;
    diagnostics::ResidualEstimatorKind estimator =
        diagnostics::ResidualEstimatorKind::Mixed;
    bool dual = true;
    bool csv = false;
    std::string source = "gaussian";
    std::string mesh_output;
};

int parse_int(const std::string &text, const char *name) {
    size_t consumed = 0;
    const int value = std::stoi(text, &consumed);
    if (consumed != text.size()) throw std::invalid_argument(std::string("invalid ") + name);
    return value;
}

double parse_double(const std::string &text, const char *name) {
    size_t consumed = 0;
    const double value = std::stod(text, &consumed);
    if (consumed != text.size() || !std::isfinite(value))
        throw std::invalid_argument(std::string("invalid ") + name);
    return value;
}

Options parse_options(int argc, char **argv) {
    Options result;
    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        auto value = [&](const std::string &prefix) { return argument.substr(prefix.size()); };
        if (argument.rfind("--k=", 0) == 0) result.wavenumber = parse_double(value("--k="), "k");
        else if (argument.rfind("--H=", 0) == 0) result.initial_H = parse_int(value("--H="), "H");
        else if (argument.rfind("--h=", 0) == 0) result.fine_h = parse_int(value("--h="), "h");
        else if (argument.rfind("--ell=", 0) == 0) result.ell = parse_int(value("--ell="), "ell");
        else if (argument.rfind("--iterations=", 0) == 0)
            result.iterations = parse_int(value("--iterations="), "iterations");
        else if (argument.rfind("--max-dofs=", 0) == 0)
            result.max_dofs = parse_int(value("--max-dofs="), "max dofs");
        else if (argument.rfind("--threads=", 0) == 0)
            result.threads = parse_int(value("--threads="), "threads");
        else if (argument.rfind("--theta=", 0) == 0)
            result.theta = parse_double(value("--theta="), "theta");
        else if (argument.rfind("--tol=", 0) == 0)
            result.tolerance = parse_double(value("--tol="), "tolerance");
        else if (argument.rfind("--q-limit=", 0) == 0)
            result.q_limit = parse_double(value("--q-limit="), "q limit");
        else if (argument.rfind("--dual-patch=", 0) == 0)
            result.dual_patch_layers = parse_int(value("--dual-patch="), "dual patch");
        else if (argument.rfind("--estimator=", 0) == 0) {
            const std::string name = value("--estimator=");
            if (name == "fine") result.estimator = diagnostics::ResidualEstimatorKind::Fine;
            else if (name == "mixed") result.estimator = diagnostics::ResidualEstimatorKind::Mixed;
            else if (name == "macro") result.estimator = diagnostics::ResidualEstimatorKind::Macro;
            else throw std::invalid_argument("estimator must be fine, mixed, or macro");
        } else if (argument == "--no-dual") result.dual = false;
        else if (argument.rfind("--source=", 0) == 0) {
            result.source = value("--source=");
            if (result.source != "gaussian" && result.source != "manufactured")
                throw std::invalid_argument("source must be gaussian or manufactured");
        }
        else if (argument == "--format=csv") result.csv = true;
        else if (argument.rfind("--mesh-out=", 0) == 0) result.mesh_output = value("--mesh-out=");
        else if (argument == "--help") {
            std::cout
                << "Usage: bench_helmholtz_adaptive [--k=8] [--H=2] [--h=7] [--ell=4] "
                   "[--iterations=6] [--theta=0.5] [--estimator=fine|mixed|macro] "
                   "[--q-limit=0.5] [--tol=0] [--max-dofs=20000] [--threads=8] "
                   "[--source=gaussian|manufactured] "
                   "[--dual-patch=1|--no-dual] [--format=csv] [--mesh-out=path]\n";
            std::exit(0);
        } else throw std::invalid_argument("unknown option: " + argument);
    }
    return result;
}

void print_csv_header() {
    std::cout
        << "iteration,coarse_nodes,coarse_elements,fine_nodes,fine_elements,min_level,max_level,"
           "marked,closure_added,H_max,q_max,q_effective,energy_error,energy_error_rel,"
           "l2_error,l2_error_rel,exact_energy_error,exact_l2_error,"
           "fine_exact_energy_error,fine_exact_l2_error,eta_fine,eta_mixed,eta_macro,"
           "eta_selected,effectivity,exact_effectivity,"
           "residual_identity,dual_spearman,dual_overlap,inf_sup,petrov_residual,corrector_residual,"
           "constraint_residual,build_ms,solve_ms,reference_ms,estimate_ms,dual_ms,"
           "mark_refine_ms,total_ms\n";
}

void print_csv(const AdaptiveIterationRecord &r) {
    std::cout << std::setprecision(17)
              << r.iteration << ',' << r.coarse_nodes << ',' << r.coarse_elements << ','
              << r.fine_nodes << ',' << r.fine_elements << ','
              << r.min_coarse_level << ',' << r.max_coarse_level << ','
              << r.marked_elements << ',' << r.closure_added_elements << ','
              << r.H_max << ',' << r.q_max << ',' << r.q_effective << ','
              << r.energy_error << ',' << r.relative_energy_error << ','
              << r.l2_error << ',' << r.relative_l2_error << ','
              << r.exact_energy_error << ',' << r.exact_l2_error << ','
              << r.fine_exact_energy_error << ',' << r.fine_exact_l2_error << ','
              << r.eta_fine << ',' << r.eta_mixed << ',' << r.eta_macro << ','
              << r.selected_estimator << ',' << r.selected_effectivity << ','
              << r.exact_effectivity << ','
              << r.residual_identity_error << ',' << r.dual_spearman << ','
              << r.dual_marking_overlap << ',' << r.inf_sup << ',' << r.petrov_residual << ','
              << r.corrector_residual << ',' << r.constraint_residual << ','
              << r.build_ms << ',' << r.solve_ms << ',' << r.reference_ms << ','
              << r.estimate_ms << ',' << r.dual_ms << ',' << r.mark_refine_ms << ','
              << r.total_ms << '\n';
}

void print_human(const AdaptiveIterationRecord &r) {
    std::cout << std::setprecision(6)
              << "iter=" << r.iteration << " coarse=" << r.coarse_elements
              << " levels=" << r.min_coarse_level << ':' << r.max_coarse_level
              << " marked/closure=" << r.marked_elements << '/' << r.closure_added_elements
              << " q=" << r.q_max << '\n'
              << "  error(E/L2)=" << r.relative_energy_error << '/' << r.relative_l2_error
              << " eta(fine/mix/macro)=" << r.eta_fine << '/' << r.eta_mixed
              << '/' << r.eta_macro << " eff=" << r.selected_effectivity << '\n'
              << "  exact LOD(E/L2)=" << r.exact_energy_error << '/' << r.exact_l2_error
              << " fine(E/L2)=" << r.fine_exact_energy_error << '/' << r.fine_exact_l2_error
              << " exact-eff=" << r.exact_effectivity << '\n'
              << "  residual identity=" << r.residual_identity_error
              << " dual rho/overlap=" << r.dual_spearman << '/' << r.dual_marking_overlap
              << " inf-sup=" << r.inf_sup << " Petrov=" << r.petrov_residual << '\n'
              << "  time build/estimate/dual/total(ms)=" << r.build_ms << '/'
              << r.estimate_ms << '/' << r.dual_ms << '/' << r.total_ms << '\n';
}

void write_mesh(const std::string &path, const AdaptiveHelmholtzResult &result) {
    if (path.empty()) return;
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot open adaptive mesh output: " + path);
    output << "element_id,level,n0,n1,n2,x0,y0,x1,y1,x2,y2,indicator_squared\n";
    for (int t = 0; t < static_cast<int>(result.final_coarse_mesh.elems.size()); ++t) {
        const Triangle &tri = result.final_coarse_mesh.elems[t];
        output << result.final_coarse_element_ids[t] << ',' << result.final_coarse_levels[t] << ','
               << tri[0] << ',' << tri[1] << ',' << tri[2];
        for (int node : tri)
            output << ',' << result.final_coarse_mesh.nodes[node].x()
                   << ',' << result.final_coarse_mesh.nodes[node].y();
        output << ',' << result.final_indicator_squared[t] << '\n';
    }
}

} // namespace

int main(int argc, char **argv) {
    try {
        const Options options = parse_options(argc, argv);
#ifdef _OPENMP
        if (options.threads > 0) omp_set_num_threads(options.threads);
#endif
        AdaptiveHelmholtzConfig config;
        config.problem.wavenumber = options.wavenumber;
        config.problem.ell = options.ell;
        config.problem.mode = HelmholtzPetrovMode::TwoSided;
        config.initial_coarse_level = options.initial_H;
        config.fine_level = options.fine_h;
        config.max_iterations = options.iterations;
        config.max_coarse_dofs = options.max_dofs;
        config.theta = options.theta;
        config.tolerance = options.tolerance;
        config.q_limit = options.q_limit;
        config.dual_patch_layers = options.dual_patch_layers;
        config.estimator = options.estimator;
        config.compute_dual_calibration = options.dual;

        ComplexFunction source;
        ComplexFunction exact;
        ComplexGradientFunction exact_gradient;
        if (options.source == "manufactured") {
            const HelmholtzManufacturedSolution manufactured =
                make_polynomial_plane_wave_solution(options.wavenumber);
            source = manufactured.source;
            exact = manufactured.value;
            exact_gradient = manufactured.gradient;
        } else {
            source = [](const Point2 &point) {
                const double dx = point.x() - 0.35;
                const double dy = point.y() - 0.55;
                return Complex(std::exp(-80.0 * (dx * dx + dy * dy)), 0.0);
            };
        }
        const AdaptiveHelmholtzResult result =
            run_adaptive_helmholtz(config, source, exact, exact_gradient);
        if (options.csv) print_csv_header();
        else std::cout << "=== Adaptive Helmholtz LOD: estimator="
                       << residual_estimator_name(options.estimator) << " ===\n";
        for (const AdaptiveIterationRecord &record : result.history) {
            if (options.csv) print_csv(record);
            else print_human(record);
        }
        if (!options.csv) std::cout << "stop: " << result.stop_reason << '\n';
        write_mesh(options.mesh_output, result);
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "bench_helmholtz_adaptive failed: " << error.what() << '\n';
        return 1;
    }
}

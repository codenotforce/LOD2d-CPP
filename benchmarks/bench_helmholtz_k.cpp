#include "helmholtz/manufactured.h"
#include "helmholtz/model.h"
#include "helmholtz/operators.h"
#include "mesh/refine.h"

#include <Eigen/Cholesky>
#include <Eigen/SVD>
#include <algorithm>
#include <array>
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

using namespace lod2d;
using namespace lod2d::helmholtz;

namespace {

struct Options {
    double wavenumber = 4.0;
    int coarse_level = -1;
    int fine_level = -1;
    int fine_gap = 8;
    int ell = -1;
    double coarse_kh_target = 1.0;
    int stability_max_dofs = 512;
    int threads = 0;
    int symbolic_cache_slots = 1;
    std::string source = "gaussian";
    std::string solver = "saddle";
    std::string factorization_reuse = "none";
    HelmholtzPetrovMode mode = HelmholtzPetrovMode::TwoSided;
    bool csv = false;
    bool csv_header = false;
    bool check = false;
};

struct ErrorMetrics {
    double energy_absolute = 0.0;
    double energy_relative = 0.0;
    double l2_absolute = 0.0;
    double l2_relative = 0.0;
};

struct Result {
    double k = 0.0;
    int H_level = 0;
    int h_level = 0;
    int ell = 0;
    std::string mode;
    int coarse_nodes = 0;
    int coarse_elements = 0;
    int fine_nodes = 0;
    int fine_elements = 0;
    double H = 0.0;
    double h = 0.0;
    double kH = 0.0;
    double kh = 0.0;
    std::string source;
    std::string solver;
    double source_l2 = 0.0;
    ErrorMetrics fem;
    ErrorMetrics lod;
    double exact_energy_norm = std::numeric_limits<double>::quiet_NaN();
    double exact_l2_norm = std::numeric_limits<double>::quiet_NaN();
    double fem_exact_energy_rel = std::numeric_limits<double>::quiet_NaN();
    double fem_exact_l2_rel = std::numeric_limits<double>::quiet_NaN();
    double lod_exact_energy_rel = std::numeric_limits<double>::quiet_NaN();
    double lod_exact_l2_rel = std::numeric_limits<double>::quiet_NaN();
    double fine_exact_energy_rel = std::numeric_limits<double>::quiet_NaN();
    double fine_exact_l2_rel = std::numeric_limits<double>::quiet_NaN();
    double normalized_lod_energy = 0.0;
    double inf_sup = std::numeric_limits<double>::quiet_NaN();
    double petrov_residual = 0.0;
    double corrector_residual = 0.0;
    double constraint_residual = 0.0;
    int corrector_threads = 1;
    int symbolic_analyses = 0;
    int symbolic_reuses = 0;
    int factorization_reuses = 0;
    double schur_residual = 0.0;
    double schur_rcond = 1.0;
    int direct_fallbacks = 0;
    HelmholtzBuildTimings build;
    double load_ms = 0.0;
    double lod_solve_ms = 0.0;
    double coarse_fem_ms = 0.0;
    double reference_ms = 0.0;
    double total_ms = 0.0;
};

double elapsed_ms(const std::chrono::steady_clock::time_point &start) {
    return std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - start).count();
}

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
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        auto value_after = [&](const std::string &prefix) {
            return argument.substr(prefix.size());
        };
        if (argument.rfind("--k=", 0) == 0) {
            options.wavenumber = parse_double(value_after("--k="), "wavenumber");
        } else if (argument.rfind("--H=", 0) == 0) {
            const std::string value = value_after("--H=");
            options.coarse_level = value == "auto" ? -1 : parse_int(value, "H level");
        } else if (argument.rfind("--h=", 0) == 0) {
            const std::string value = value_after("--h=");
            options.fine_level = value == "auto" ? -1 : parse_int(value, "h level");
        } else if (argument.rfind("--fine-gap=", 0) == 0) {
            options.fine_gap = parse_int(value_after("--fine-gap="), "fine gap");
        } else if (argument.rfind("--ell=", 0) == 0) {
            const std::string value = value_after("--ell=");
            options.ell = value == "auto" ? -1 : parse_int(value, "ell");
        } else if (argument.rfind("--kH-target=", 0) == 0) {
            options.coarse_kh_target = parse_double(value_after("--kH-target="), "kH target");
        } else if (argument.rfind("--stability-max-dofs=", 0) == 0) {
            options.stability_max_dofs = parse_int(
                value_after("--stability-max-dofs="), "stability max dofs");
        } else if (argument.rfind("--threads=", 0) == 0) {
            options.threads = parse_int(value_after("--threads="), "threads");
        } else if (argument.rfind("--symbolic-cache-slots=", 0) == 0) {
            options.symbolic_cache_slots = parse_int(
                value_after("--symbolic-cache-slots="), "symbolic cache slots");
        } else if (argument.rfind("--source=", 0) == 0) {
            options.source = value_after("--source=");
        } else if (argument.rfind("--solver=", 0) == 0) {
            options.solver = value_after("--solver=");
        } else if (argument.rfind("--factorization-reuse=", 0) == 0) {
            options.factorization_reuse = value_after("--factorization-reuse=");
        } else if (argument.rfind("--mode=", 0) == 0) {
            const std::string value = value_after("--mode=");
            if (value == "two-sided") options.mode = HelmholtzPetrovMode::TwoSided;
            else if (value == "test-only") options.mode = HelmholtzPetrovMode::CorrectedTestOnly;
            else throw std::invalid_argument("mode must be two-sided or test-only");
        } else if (argument == "--format=csv") {
            options.csv = true;
        } else if (argument == "--csv-header") {
            options.csv_header = true;
        } else if (argument == "--check") {
            options.check = true;
        } else if (argument == "--help") {
            std::cout
                << "Usage: bench_helmholtz_k [--k=4] [--H=auto] [--h=auto] "
                   "[--fine-gap=8] [--ell=auto] [--kH-target=1] "
                   "[--source=gaussian|manufactured] "
                   "[--solver=saddle|schur] [--threads=0] "
                   "[--symbolic-cache-slots=1] "
                   "[--factorization-reuse=none|identical] "
                   "[--mode=two-sided|test-only] [--stability-max-dofs=512] "
                   "[--format=csv|--csv-header] [--check]\n";
            std::exit(0);
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.wavenumber <= 0.0) throw std::invalid_argument("k must be positive");
    if (options.coarse_kh_target <= 0.0) throw std::invalid_argument("kH target must be positive");
    if (options.fine_gap < 0) throw std::invalid_argument("fine gap must be nonnegative");
    if (options.stability_max_dofs < 0)
        throw std::invalid_argument("stability max dofs must be nonnegative");
    if (options.threads < 0)
        throw std::invalid_argument("threads must be nonnegative");
    if (options.symbolic_cache_slots <= 0)
        throw std::invalid_argument("symbolic cache slots must be positive");
    if (options.source != "gaussian" && options.source != "manufactured")
        throw std::invalid_argument("source must be gaussian or manufactured");
    if (options.solver != "saddle" && options.solver != "schur")
        throw std::invalid_argument("solver must be saddle or schur");
    if (options.factorization_reuse != "none"
        && options.factorization_reuse != "identical")
        throw std::invalid_argument(
            "factorization reuse must be none or identical");
    return options;
}

int automatic_coarse_level(double k, double target) {
    TriMesh mesh = make_helmholtz_unit_square_mesh();
    for (int level = 0; level <= 24; ++level) {
        if (k * max_element_diameter(mesh) <= target) return level;
        mesh = refine_nvb(mesh).mesh;
    }
    throw std::runtime_error("automatic H selection exceeded 24 NVB levels");
}

double quadratic_form(
    const Eigen::SparseMatrix<double> &matrix,
    const ComplexVector &vector) {
    const ComplexVector product = matrix.cast<Complex>() * vector;
    return std::max(0.0, std::real(vector.dot(product)));
}

ErrorMetrics discrete_error(
    const HelmholtzOperators &operators,
    const ComplexVector &approximation,
    const ComplexVector &reference) {
    const ComplexVector difference = reference - approximation;
    const double k2 = operators.wavenumber * operators.wavenumber;
    const double reference_l2_squared = quadratic_form(operators.mass, reference);
    const double error_l2_squared = quadratic_form(operators.mass, difference);
    const double reference_energy_squared = quadratic_form(operators.stiffness, reference)
                                          + k2 * reference_l2_squared;
    const double error_energy_squared = quadratic_form(operators.stiffness, difference)
                                      + k2 * error_l2_squared;
    ErrorMetrics result;
    result.energy_absolute = std::sqrt(error_energy_squared);
    result.energy_relative = result.energy_absolute
                           / std::max(std::sqrt(reference_energy_squared), 1e-30);
    result.l2_absolute = std::sqrt(error_l2_squared);
    result.l2_relative = result.l2_absolute
                       / std::max(std::sqrt(reference_l2_squared), 1e-30);
    return result;
}

ErrorMetrics exact_error(
    const TriMesh &mesh,
    const ComplexVector &approximation,
    double wavenumber,
    const HelmholtzManufacturedSolution &manufactured,
    double exact_energy_norm,
    double exact_l2_norm) {
    const HelmholtzError error = compute_helmholtz_error(
        mesh, approximation, wavenumber,
        manufactured.value, manufactured.gradient);
    ErrorMetrics result;
    result.energy_absolute = error.energy;
    result.energy_relative = error.energy / std::max(exact_energy_norm, 1e-30);
    result.l2_absolute = error.l2;
    result.l2_relative = error.l2 / std::max(exact_l2_norm, 1e-30);
    return result;
}

double source_l2_norm(
    const TriMesh &mesh,
    const Eigen::SparseMatrix<double> &mass,
    const ComplexFunction &source) {
    ComplexVector values(mesh.nodes.size());
    for (int i = 0; i < values.size(); ++i) values(i) = source(mesh.nodes[i]);
    return std::sqrt(quadratic_form(mass, values));
}

double energy_inf_sup(const HelmholtzLodModel &model, int max_dofs) {
    const int dofs = model.coarse_operator().rows();
    if (dofs == 0 || dofs > max_dofs) return std::numeric_limits<double>::quiet_NaN();

    Eigen::SparseMatrix<double> energy = model.operators().stiffness;
    energy += model.config().wavenumber * model.config().wavenumber
            * model.operators().mass;
    const ComplexSparseMatrix complex_energy = energy.cast<Complex>();
    const ComplexMatrix trial_gram(
        model.trial_basis().adjoint() * complex_energy * model.trial_basis());
    const ComplexMatrix test_gram(
        model.test_basis().adjoint() * complex_energy * model.test_basis());
    Eigen::LLT<ComplexMatrix> trial_llt(trial_gram);
    Eigen::LLT<ComplexMatrix> test_llt(test_gram);
    if (trial_llt.info() != Eigen::Success || test_llt.info() != Eigen::Success)
        return std::numeric_limits<double>::quiet_NaN();

    const ComplexMatrix coarse(model.coarse_operator());
    const ComplexMatrix left_scaled = test_llt.matrixL().solve(coarse);
    const ComplexMatrix scaled = trial_llt.matrixL()
        .solve(left_scaled.adjoint()).adjoint();
    Eigen::JacobiSVD<ComplexMatrix> svd(scaled, Eigen::ComputeThinU | Eigen::ComputeThinV);
    return svd.singularValues().minCoeff();
}

std::string mode_name(HelmholtzPetrovMode mode) {
    return mode == HelmholtzPetrovMode::TwoSided ? "two-sided" : "test-only";
}

const char *csv_header() {
    return "k,H_level,h_level,ell,mode,source,solver,coarse_nodes,coarse_elements,fine_nodes,fine_elements,"
           "H,h,kH,kh,source_l2,fem_energy_rel,fem_l2_rel,lod_energy_rel,lod_l2_rel,"
           "exact_energy_norm,exact_l2_norm,fem_exact_energy_rel,fem_exact_l2_rel,"
           "lod_exact_energy_rel,lod_exact_l2_rel,fine_exact_energy_rel,fine_exact_l2_rel,"
           "lod_energy_over_Hf,inf_sup,petrov_residual,corrector_residual,constraint_residual,"
           "corrector_threads,symbolic_analyses,symbolic_reuses,factorization_reuses,"
           "schur_residual,schur_rcond,direct_fallbacks,"
           "mesh_interp_ms,operators_ms,correctors_ms,basis_factor_ms,build_total_ms,"
           "load_ms,lod_solve_ms,coarse_fem_ms,reference_ms,total_ms";
}

void print_csv(const Result &r) {
    std::cout << std::setprecision(17)
              << r.k << ',' << r.H_level << ',' << r.h_level << ',' << r.ell << ',' << r.mode << ','
              << r.source << ',' << r.solver << ','
              << r.coarse_nodes << ',' << r.coarse_elements << ','
              << r.fine_nodes << ',' << r.fine_elements << ','
              << r.H << ',' << r.h << ',' << r.kH << ',' << r.kh << ',' << r.source_l2 << ','
              << r.fem.energy_relative << ',' << r.fem.l2_relative << ','
              << r.lod.energy_relative << ',' << r.lod.l2_relative << ','
              << r.exact_energy_norm << ',' << r.exact_l2_norm << ','
              << r.fem_exact_energy_rel << ',' << r.fem_exact_l2_rel << ','
              << r.lod_exact_energy_rel << ',' << r.lod_exact_l2_rel << ','
              << r.fine_exact_energy_rel << ',' << r.fine_exact_l2_rel << ','
              << r.normalized_lod_energy << ',' << r.inf_sup << ',' << r.petrov_residual << ','
              << r.corrector_residual << ',' << r.constraint_residual << ','
              << r.corrector_threads << ',' << r.symbolic_analyses << ','
              << r.symbolic_reuses << ',' << r.factorization_reuses << ','
              << r.schur_residual << ','
              << r.schur_rcond << ',' << r.direct_fallbacks << ','
              << r.build.mesh_and_interpolation_ms << ',' << r.build.operators_ms << ','
              << r.build.correctors_ms << ',' << r.build.basis_and_factorization_ms << ','
              << r.build.total_ms << ',' << r.load_ms << ',' << r.lod_solve_ms << ','
              << r.coarse_fem_ms << ',' << r.reference_ms << ',' << r.total_ms << '\n';
}

void print_human(const Result &r) {
    std::cout << std::setprecision(8)
              << "=== Helmholtz wave-number benchmark ===\n"
              << "k=" << r.k << " H-level=" << r.H_level << " h-level=" << r.h_level
              << " ell=" << r.ell << " mode=" << r.mode << " source=" << r.source
              << " solver=" << r.solver << '\n'
              << "coarse: " << r.coarse_nodes << " vertices, " << r.coarse_elements
              << " triangles; fine: " << r.fine_nodes << " vertices, " << r.fine_elements
              << " triangles\n"
              << "H=" << r.H << " h=" << r.h << " kH=" << r.kH << " kh=" << r.kh << '\n'
              << "\nRelative errors against fine FEM\n"
              << "standard FEM: energy=" << r.fem.energy_relative << " L2=" << r.fem.l2_relative << '\n'
              << "Petrov LOD : energy=" << r.lod.energy_relative << " L2=" << r.lod.l2_relative << '\n'
              << "\nRelative errors against manufactured exact solution\n"
              << "standard FEM: energy=" << r.fem_exact_energy_rel
              << " L2=" << r.fem_exact_l2_rel << '\n'
              << "Petrov LOD : energy=" << r.lod_exact_energy_rel
              << " L2=" << r.lod_exact_l2_rel << '\n'
              << "fine FEM   : energy=" << r.fine_exact_energy_rel
              << " L2=" << r.fine_exact_l2_rel << '\n'
              << "LOD energy/(H||f||)=" << r.normalized_lod_energy
              << " inf-sup=" << r.inf_sup << '\n'
              << "Petrov residual=" << r.petrov_residual
              << " corrector=" << r.corrector_residual
              << " constraint=" << r.constraint_residual << '\n'
              << "corrector threads=" << r.corrector_threads << " symbolic analyze/reuse="
              << r.symbolic_analyses << '/' << r.symbolic_reuses
              << " numeric reuse=" << r.factorization_reuses
              << " Schur residual/rcond=" << r.schur_residual << '/'
              << r.schur_rcond << " fallbacks=" << r.direct_fallbacks << '\n'
              << "\nTimings (ms)\n"
              << "mesh+I_H=" << r.build.mesh_and_interpolation_ms
              << " operators=" << r.build.operators_ms
              << " correctors=" << r.build.correctors_ms
              << " basis+factor=" << r.build.basis_and_factorization_ms << '\n'
              << "LOD build=" << r.build.total_ms << " LOD solve=" << r.lod_solve_ms
              << " coarse FEM=" << r.coarse_fem_ms << " reference=" << r.reference_ms
              << " total=" << r.total_ms << '\n';
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

        const auto total_start = std::chrono::steady_clock::now();
        Result result;
        result.k = options.wavenumber;
        result.H_level = options.coarse_level >= 0
            ? options.coarse_level
            : automatic_coarse_level(options.wavenumber, options.coarse_kh_target);
        result.h_level = options.fine_level >= 0
            ? options.fine_level
            : result.H_level + options.fine_gap;
        if (result.h_level < result.H_level)
            throw std::invalid_argument("h level must be at least H level");
        result.ell = options.ell >= 0
            ? options.ell
            : std::max(1, static_cast<int>(std::ceil(std::log2(options.wavenumber))));
        result.mode = mode_name(options.mode);
        result.source = options.source;
        result.solver = options.solver;

        HelmholtzProblemConfig config;
        config.H = result.H_level;
        config.h = result.h_level;
        config.ell = result.ell;
        config.wavenumber = options.wavenumber;
        config.mode = options.mode;
        config.patch_solver.kind = options.solver == "schur"
            ? HelmholtzPatchSolverKind::DirectSchur
            : HelmholtzPatchSolverKind::DirectSaddle;
        config.patch_solver.symbolic_cache_slots = options.symbolic_cache_slots;
        config.patch_solver.reuse_identical_factorization =
            options.factorization_reuse == "identical";
        config.patch_solver.fallback_to_direct = false;
        HelmholtzLodModel model = HelmholtzLodModel::build(config);
        result.build = model.build_timings();
        result.coarse_nodes = model.problem().coarse.nodes.size();
        result.coarse_elements = model.problem().coarse.elems.size();
        result.fine_nodes = model.problem().fine.nodes.size();
        result.fine_elements = model.problem().fine.elems.size();
        result.H = max_element_diameter(model.problem().coarse);
        result.h = max_element_diameter(model.problem().fine);
        result.kH = result.k * result.H;
        result.kh = result.k * result.h;

        const ComplexFunction gaussian_source = [](const Point2 &point) {
            const double dx = point.x() - 0.35;
            const double dy = point.y() - 0.55;
            return Complex(std::exp(-40.0 * (dx * dx + dy * dy)), 0.0);
        };
        HelmholtzManufacturedSolution manufactured;
        const HelmholtzManufacturedSolution *manufactured_pointer = nullptr;
        ComplexFunction source = gaussian_source;
        if (options.source == "manufactured") {
            manufactured = make_polynomial_plane_wave_solution(options.wavenumber);
            manufactured_pointer = &manufactured;
            source = manufactured.source;
        }

        auto stage_start = std::chrono::steady_clock::now();
        const ComplexVector fine_load = assemble_helmholtz_load(model.problem().fine, source);
        result.load_ms = elapsed_ms(stage_start);
        result.source_l2 = source_l2_norm(
            model.problem().fine, model.operators().mass, source);

        stage_start = std::chrono::steady_clock::now();
        const HelmholtzLodSolution lod_solution = model.solve_load(fine_load);
        result.lod_solve_ms = elapsed_ms(stage_start);
        result.petrov_residual = lod_solution.petrov_residual;

        stage_start = std::chrono::steady_clock::now();
        const HelmholtzOperators coarse_operators = assemble_helmholtz_operators(
            model.problem().coarse, options.wavenumber);
        const ComplexVector coarse_load = assemble_helmholtz_load(model.problem().coarse, source);
        const ComplexVector coarse_values = solve_helmholtz_fem(coarse_operators, coarse_load);
        const ComplexVector coarse_on_fine = model.problem().coarse_to_fine.cast<Complex>()
                                           * coarse_values;
        result.coarse_fem_ms = elapsed_ms(stage_start);

        stage_start = std::chrono::steady_clock::now();
        const ComplexVector reference = model.solve_fine_reference(fine_load);
        result.reference_ms = elapsed_ms(stage_start);

        result.fem = discrete_error(model.operators(), coarse_on_fine, reference);
        result.lod = discrete_error(model.operators(), lod_solution.fine_values, reference);
        if (manufactured_pointer) {
            const ComplexVector zero = ComplexVector::Zero(reference.size());
            const HelmholtzError exact_norm = compute_helmholtz_error(
                model.problem().fine, zero, options.wavenumber,
                manufactured.value, manufactured.gradient);
            result.exact_energy_norm = exact_norm.energy;
            result.exact_l2_norm = exact_norm.l2;
            const ErrorMetrics fem_exact = exact_error(
                model.problem().fine, coarse_on_fine, options.wavenumber,
                manufactured, result.exact_energy_norm, result.exact_l2_norm);
            const ErrorMetrics lod_exact = exact_error(
                model.problem().fine, lod_solution.fine_values, options.wavenumber,
                manufactured, result.exact_energy_norm, result.exact_l2_norm);
            const ErrorMetrics fine_exact = exact_error(
                model.problem().fine, reference, options.wavenumber,
                manufactured, result.exact_energy_norm, result.exact_l2_norm);
            result.fem_exact_energy_rel = fem_exact.energy_relative;
            result.fem_exact_l2_rel = fem_exact.l2_relative;
            result.lod_exact_energy_rel = lod_exact.energy_relative;
            result.lod_exact_l2_rel = lod_exact.l2_relative;
            result.fine_exact_energy_rel = fine_exact.energy_relative;
            result.fine_exact_l2_rel = fine_exact.l2_relative;
        }
        result.normalized_lod_energy = result.lod.energy_absolute
            / std::max(result.H * result.source_l2, 1e-30);
        result.inf_sup = energy_inf_sup(model, options.stability_max_dofs);
        result.corrector_residual = model.correctors().diagnostics.max_primal_residual;
        result.constraint_residual = model.correctors().diagnostics.max_constraint_residual;
        result.corrector_threads = model.correctors().diagnostics.parallel_threads;
        result.symbolic_analyses = model.correctors().diagnostics.symbolic_analyses;
        result.symbolic_reuses = model.correctors().diagnostics.symbolic_reuses;
        result.factorization_reuses =
            model.correctors().diagnostics.factorization_reuses;
        result.schur_residual = model.correctors().diagnostics.max_schur_residual;
        result.schur_rcond =
            model.correctors().diagnostics.min_schur_reciprocal_condition;
        result.direct_fallbacks = model.correctors().diagnostics.direct_fallbacks;
        result.total_ms = elapsed_ms(total_start);

        if (options.check) {
            if (!(result.petrov_residual < 1e-8
                  && result.corrector_residual < 1e-8
                  && result.constraint_residual < 1e-8))
                throw std::runtime_error("Helmholtz pollution residual check failed");
            if (options.solver == "schur"
                && (!(result.schur_residual < 1e-8)
                    || !(result.schur_rcond > 1e-14)
                    || result.direct_fallbacks != 0))
                throw std::runtime_error("Helmholtz pollution Schur check failed");
            if (manufactured_pointer) {
                const std::array<double, 8> exact_values{{
                    result.exact_energy_norm, result.exact_l2_norm,
                    result.fem_exact_energy_rel, result.fem_exact_l2_rel,
                    result.lod_exact_energy_rel, result.lod_exact_l2_rel,
                    result.fine_exact_energy_rel, result.fine_exact_l2_rel}};
                for (double value : exact_values)
                    if (!(std::isfinite(value) && value > 0.0))
                        throw std::runtime_error(
                            "manufactured pollution metric is not finite and positive");
            }
        }

        if (options.csv) print_csv(result);
        else print_human(result);
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "bench_helmholtz_k failed: " << error.what() << '\n';
        return 1;
    }
}

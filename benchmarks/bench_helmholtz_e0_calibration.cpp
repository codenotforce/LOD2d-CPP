#include "helmholtz/adaptive/error_control.h"
#include "helmholtz/adaptive/kernel_residual.h"
#include "helmholtz/adaptive/practical_driver.h"
#include "helmholtz/adaptive/reference_retraction.h"
#include "helmholtz/benchmarks/paper_cases.h"
#include "helmholtz/boundary.h"
#include "helmholtz/model.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace lod2d;
using namespace lod2d::helmholtz;
using namespace lod2d::helmholtz::adaptive;
using namespace lod2d::helmholtz::benchmarks;
using namespace lod2d::helmholtz::experiments;

namespace {

constexpr double kWavenumber = 16.0;
constexpr int kCoarseLevel = 6;
constexpr int kReferenceLevel = 8;
constexpr double kRhoStar = 0.25;
constexpr double kDoerflerTheta = 0.5;
constexpr double kSafetyFactor = 1.25;

struct Arguments {
    std::filesystem::path output_directory;
    bool check = false;
};

struct CalibrationRow {
    int ell = 0;
    double theta_loc = 0.0;
    double direct_delta = 0.0;
    double upper_certificate = 0.0;
    double reference_energy_error = 0.0;
    double reference_relative_energy_error = 0.0;
    double eta_H = 0.0;
    double required_multiplier = 0.0;
    double maximum_constraint_residual = 0.0;
    double maximum_riesz_residual = 0.0;
    double maximum_energy_identity_error = 0.0;
    double allocation_relative_error = 0.0;
    std::size_t effectivity_count = 0;
    std::size_t effectivity_excluded = 0;
    double effectivity_minimum = 0.0;
    double effectivity_median = 0.0;
    double effectivity_percentile_90 = 0.0;
    double effectivity_maximum = 0.0;
    SmallMatrixLocalizationValidation direct_validation;
};

struct FrozenParameters {
    int ell = 2;
    double theta_loc = 0.0;
    double C0_usr = 0.0;
    double C1_usr = 0.0;
    double rho_star = kRhoStar;
};

struct DriverSmokeEvidence {
    std::size_t journal_size = 0;
    int final_ell = 0;
    double diagnostic_c_H = 0.0;
    double diagnostic_tolerance = 0.0;
    bool increased_global_ell = false;
    bool accepted_localization = false;
    bool completed = false;
};

Arguments parse_arguments(const int argc, char **argv) {
    Arguments result;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument.starts_with("--output-dir=")) {
            result.output_directory = argument.substr(
                std::string("--output-dir=").size());
        } else if (argument == "--check") {
            result.check = true;
        } else {
            throw std::invalid_argument("unknown argument: " + argument);
        }
    }
    if (result.output_directory.empty()) {
        throw std::invalid_argument("required: --output-dir=DIR [--check]");
    }
    return result;
}

std::string numeric(const double value) {
    std::ostringstream out;
    out << std::setprecision(17) << value;
    return out.str();
}

std::vector<int> free_coarse_nodes(const TriMesh &mesh) {
    std::vector<char> is_dirichlet(mesh.nodes.size(), false);
    for (const int node : dirichlet_nodes(mesh)) is_dirichlet[node] = true;
    std::vector<int> result;
    for (int node = 0; node < static_cast<int>(mesh.nodes.size()); ++node) {
        if (!is_dirichlet[node]) result.push_back(node);
    }
    return result;
}

HelmholtzLodModel build_reference_model(
    const PaperCaseData &data,
    const ReferenceEpochHierarchy &hierarchy,
    const int ell) {
    HelmholtzProblemConfig config;
    config.H = *std::min_element(
        hierarchy.coarse_levels().begin(), hierarchy.coarse_levels().end());
    config.h = hierarchy.reference_level();
    config.ell = ell;
    config.wavenumber = data.wavenumber;
    config.initial_mesh = data.initial_mesh;
    config.quadrature_context = data.quadrature_context;
    return HelmholtzLodModel::build_adaptive(
        config,
        hierarchy.coarse_mesh(),
        hierarchy.coarse_levels(),
        hierarchy.reference_mesh(),
        hierarchy.reference_element_levels());
}

double quantile(const std::vector<double> &sorted, const double probability) {
    if (sorted.empty()) return 0.0;
    const double position = probability * static_cast<double>(sorted.size() - 1);
    const auto lower = static_cast<std::size_t>(std::floor(position));
    const auto upper = static_cast<std::size_t>(std::ceil(position));
    const double fraction = position - static_cast<double>(lower);
    return (1.0 - fraction) * sorted[lower] + fraction * sorted[upper];
}

void populate_effectivity(
    CalibrationRow &row,
    const ReferenceEpochHierarchy &hierarchy,
    const HelmholtzOperators &operators,
    const ReferenceResidualRiesz &estimate,
    const ComplexVector &reference,
    const ComplexVector &candidate) {
    std::vector<double> ratios;
    const double global_scale = std::max({
        1.0,
        compute_discrete_helmholtz_error(
            hierarchy.reference_mesh(), operators, reference,
            ComplexVector::Zero(reference.size())).energy,
        compute_discrete_helmholtz_error(
            hierarchy.reference_mesh(), operators, candidate,
            ComplexVector::Zero(candidate.size())).energy});
    const double zero_tolerance = 64.0
        * std::numeric_limits<double>::epsilon() * global_scale;
    for (std::size_t index = 0; index < estimate.patches.size(); ++index) {
        const double local_error = compute_local_discrete_helmholtz_error(
            hierarchy.reference_mesh(), operators, reference, candidate,
            estimate.patches[index].enlarged_discrete_elements).energy;
        if (local_error <= zero_tolerance) {
            ++row.effectivity_excluded;
            continue;
        }
        const double ratio = estimate.local_results[index].eta / local_error;
        if (std::isfinite(ratio)) ratios.push_back(ratio);
    }
    std::sort(ratios.begin(), ratios.end());
    row.effectivity_count = ratios.size();
    if (ratios.empty()) return;
    row.effectivity_minimum = ratios.front();
    row.effectivity_median = quantile(ratios, 0.5);
    row.effectivity_percentile_90 = quantile(ratios, 0.9);
    row.effectivity_maximum = ratios.back();
}

FrozenParameters derive_frozen_parameters(
    const std::vector<CalibrationRow> &rows) {
    if (rows.size() != 4) {
        throw std::invalid_argument("E0 calibration requires ell=1,2,3,4");
    }
    FrozenParameters frozen;
    // Keep ell=2 accepted with a small observation margin while rejecting
    // the visibly under-localized ell=1 result.
    frozen.theta_loc = 1.10 * rows[1].theta_loc;
    double asymptotic_multiplier = std::numeric_limits<double>::infinity();
    for (const CalibrationRow &row : rows) {
        asymptotic_multiplier = std::min(
            asymptotic_multiplier, row.required_multiplier);
    }
    // C0 covers the observed well-localized asymptote.  C1 covers only the
    // excess multiplier correlated with localization defect, avoiding the
    // ill-conditioned split obtained by forcing C0 to half the baseline.
    frozen.C0_usr = kSafetyFactor * asymptotic_multiplier;
    for (const CalibrationRow &row : rows) {
        frozen.C1_usr = std::max(
            frozen.C1_usr,
            kSafetyFactor
                * (row.required_multiplier - asymptotic_multiplier)
                / row.theta_loc);
    }
    frozen.C1_usr = std::max(0.0, frozen.C1_usr);
    return frozen;
}

void check_acceptance(
    const std::vector<CalibrationRow> &rows,
    const FrozenParameters &frozen,
    const double rho_ambient) {
    if (!(rho_ambient <= kRhoStar + 1e-12)) {
        throw std::runtime_error("E0 ambient ratio gate failed");
    }
    if (!(rows[0].theta_loc > rows[1].theta_loc
          && rows[1].theta_loc >= rows[2].theta_loc * (1.0 - 2e-8)
          && rows[2].theta_loc >= rows[3].theta_loc * (1.0 - 2e-8))) {
        throw std::runtime_error("E0 Theta_loc is not decreasing with ell");
    }
    if (!(rows[0].theta_loc > frozen.theta_loc
          && rows[1].theta_loc < frozen.theta_loc)) {
        throw std::runtime_error("E0 theta_loc does not separate ell=1 from ell=2");
    }
    for (const CalibrationRow &row : rows) {
        if (!row.direct_validation.one_sided_control_holds
            || row.direct_validation.decomposition_relative_residual > 2e-9
            || row.maximum_constraint_residual > 2e-9
            || row.maximum_riesz_residual > 2e-9
            || row.maximum_energy_identity_error > 2e-9
            || row.allocation_relative_error > 2e-12) {
            throw std::runtime_error("E0 certificate or Riesz identity gate failed");
        }
        if (row.effectivity_count == 0
            || !(row.effectivity_median > 0.0)
            || row.effectivity_percentile_90
                > 20.0 * row.effectivity_median) {
            throw std::runtime_error(
                "E0 smooth-case local effectivity has anomalous concentration");
        }
        const double practical_multiplier =
            frozen.C0_usr + frozen.C1_usr * row.theta_loc;
        if (practical_multiplier + 1e-14
            < kSafetyFactor * row.required_multiplier) {
            throw std::runtime_error("E0 frozen constants do not cover observations");
        }
    }
}

DriverSmokeEvidence run_driver_smoke(
    const PaperCaseData &data,
    const ReferenceEpochHierarchy &calibration_hierarchy,
    const std::vector<CalibrationRow> &rows,
    const FrozenParameters &frozen) {
    PracticalDriverProblem problem;
    problem.initial_mesh = data.initial_mesh;
    problem.source = data.source;
    problem.quadrature_context = data.quadrature_context;
    PracticalDriverConfig config;
    config.initial_coarse_level = kCoarseLevel;
    config.reference_level = kReferenceLevel;
    config.ell0 = 1;
    config.ell_max = 4;
    config.wavenumber = kWavenumber;
    // E0 keeps H fixed.  This diagnostic-only value admits exactly the
    // chosen calibration mesh and is not a frozen production parameter.
    config.c_H = 1.01 * kWavenumber
        * max_element_diameter(calibration_hierarchy.coarse_mesh());
    config.theta_loc = frozen.theta_loc;
    config.C0_usr = frozen.C0_usr;
    config.C1_usr = frozen.C1_usr;
    config.theta_H = kDoerflerTheta;
    config.rho_star = frozen.rho_star;
    config.tolerance_reference = 1.01
        * (frozen.C0_usr + frozen.C1_usr * rows[1].theta_loc)
        * rows[1].eta_H;
    config.localization_eigen.maximum_iterations = 500;
    config.localization_eigen.relative_tolerance = 2e-12;
    config.localization_eigen.dense_cross_check_max_dimension = 64;
    config.limits.maximum_iterations = 20;
    config.limits.maximum_H_steps = 0;
    config.limits.maximum_unknowns = 100000;
    config.limits.maximum_coarse_elements = 10000;
    config.limits.maximum_ambient_elements = 100000;
    config.limits.maximum_wall_seconds = 300.0;

    PracticalAdaptiveDriver driver(std::move(problem), config);
    const PracticalDriverResult result = driver.run();
    DriverSmokeEvidence evidence;
    evidence.journal_size = result.journal.size();
    evidence.final_ell = result.ell;
    evidence.diagnostic_c_H = config.c_H;
    evidence.diagnostic_tolerance = config.tolerance_reference;
    for (const PracticalIterationRecord &record : result.journal) {
        evidence.increased_global_ell = evidence.increased_global_ell
            || record.action == PracticalDriverAction::IncreaseGlobalEll;
        evidence.accepted_localization = evidence.accepted_localization
            || record.action == PracticalDriverAction::AcceptLocalization;
        evidence.completed = evidence.completed
            || record.action == PracticalDriverAction::Complete;
        if (record.action == PracticalDriverAction::RefineCoarse
            || record.action
                == PracticalDriverAction::StopReferenceRefreshRequired) {
            throw std::runtime_error(
                "E0 fixed-grid driver smoke changed H or requested refresh");
        }
    }
    if (result.state != PracticalDriverState::Converged
        || result.ell != frozen.ell
        || !evidence.increased_global_ell
        || !evidence.accepted_localization
        || !evidence.completed) {
        throw std::runtime_error(
            "E0 driver smoke did not execute ell=1 -> ell=2 -> Complete");
    }
    return evidence;
}

void write_csv(
    const std::filesystem::path &path,
    const std::vector<CalibrationRow> &rows) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "schema_version,case,kappa,coarse_level,reference_level,rho_star,ell,"
           "Theta_loc,direct_delta,upper_certificate,reference_energy_error,"
           "reference_relative_energy_error,eta_H,required_multiplier,"
           "max_constraint_residual,max_riesz_residual,max_energy_identity_error,"
           "allocation_relative_error,effectivity_count,effectivity_excluded,"
           "effectivity_min,effectivity_median,effectivity_p90,effectivity_max,"
           "ambient_kernel_dimension,local_kernel_columns,C_J,C_ret,c_W,C_sd,"
           "decomposition_relative_residual,one_sided_control_holds\n";
    for (const CalibrationRow &row : rows) {
        out << "1,R1," << numeric(kWavenumber) << ',' << kCoarseLevel << ','
            << kReferenceLevel << ',' << numeric(kRhoStar) << ',' << row.ell << ','
            << numeric(row.theta_loc) << ',' << numeric(row.direct_delta) << ','
            << numeric(row.upper_certificate) << ','
            << numeric(row.reference_energy_error) << ','
            << numeric(row.reference_relative_energy_error) << ','
            << numeric(row.eta_H) << ',' << numeric(row.required_multiplier) << ','
            << numeric(row.maximum_constraint_residual) << ','
            << numeric(row.maximum_riesz_residual) << ','
            << numeric(row.maximum_energy_identity_error) << ','
            << numeric(row.allocation_relative_error) << ','
            << row.effectivity_count << ',' << row.effectivity_excluded << ','
            << numeric(row.effectivity_minimum) << ','
            << numeric(row.effectivity_median) << ','
            << numeric(row.effectivity_percentile_90) << ','
            << numeric(row.effectivity_maximum) << ','
            << row.direct_validation.ambient_kernel_dimension << ','
            << row.direct_validation.local_kernel_columns << ','
            << numeric(row.direct_validation.projector_stability_constant) << ','
            << numeric(row.direct_validation.retraction_stability_constant) << ','
            << numeric(row.direct_validation.ambient_kernel_coercivity) << ','
            << numeric(row.direct_validation.stable_decomposition_constant) << ','
            << numeric(row.direct_validation.decomposition_relative_residual) << ','
            << (row.direct_validation.one_sided_control_holds ? "true" : "false")
            << '\n';
    }
}

void write_frozen_json(
    const std::filesystem::path &path,
    const FrozenParameters &frozen,
    const double rho_ambient,
    const DriverSmokeEvidence &smoke) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "{\n"
        << "  \"schema_version\": 1,\n"
        << "  \"scope\": \"E0-R1-k16-calibration\",\n"
        << "  \"claim\": \"implementation-study\",\n"
        << "  \"status\": \"passed\",\n"
        << "  \"calibration\": {\"case\": \"R1\", \"kappa\": 16, "
           "\"coarse_level\": 6, \"reference_level\": 8, "
           "\"ell_sweep\": [1, 2, 3, 4], \"safety_factor\": 1.25},\n"
        << "  \"frozen_parameters\": {\n"
        << "    \"ell\": " << frozen.ell << ",\n"
        << "    \"theta_loc\": " << numeric(frozen.theta_loc) << ",\n"
        << "    \"C0_usr\": " << numeric(frozen.C0_usr) << ",\n"
        << "    \"C1_usr\": " << numeric(frozen.C1_usr) << ",\n"
        << "    \"rho_star\": " << numeric(frozen.rho_star) << "\n"
        << "  },\n"
        << "  \"observed_rho_ambient\": " << numeric(rho_ambient) << ",\n"
        << "  \"driver_smoke\": {\n"
        << "    \"status\": \"passed\",\n"
        << "    \"journal_size\": " << smoke.journal_size << ",\n"
        << "    \"final_ell\": " << smoke.final_ell << ",\n"
        << "    \"diagnostic_c_H\": " << numeric(smoke.diagnostic_c_H) << ",\n"
        << "    \"diagnostic_tolerance\": "
        << numeric(smoke.diagnostic_tolerance) << ",\n"
        << "    \"increase_global_ell\": "
        << (smoke.increased_global_ell ? "true" : "false") << ",\n"
        << "    \"accept_localization\": "
        << (smoke.accepted_localization ? "true" : "false") << ",\n"
        << "    \"complete\": " << (smoke.completed ? "true" : "false")
        << "\n  },\n"
        << "  \"raw_data\": \"e0-calibration.csv\",\n"
        << "  \"notes\": [\n"
        << "    \"Constants are conservative observed-effectivity parameters, not theorem constants.\",\n"
        << "    \"Reference errors are post-processing diagnostics and do not enter marking or stopping.\"\n"
        << "  ]\n"
        << "}\n";
}

} // namespace

int main(const int argc, char **argv) {
    try {
        const Arguments arguments = parse_arguments(argc, argv);
        std::filesystem::create_directories(arguments.output_directory);
        const PaperCaseData data = make_paper_case(PaperCase::R1, kWavenumber);
        ReferenceEpochHierarchy hierarchy(
            data.initial_mesh, kCoarseLevel, kReferenceLevel);
        hierarchy.enforce_ambient_ratio(kRhoStar);
        const double rho_ambient = hierarchy.ambient_ratio();
        const HelmholtzOperators ambient_operators = assemble_helmholtz_operators(
            hierarchy.ambient_mesh(), kWavenumber);
        const std::vector<int> basis_nodes = free_coarse_nodes(
            hierarchy.coarse_mesh());
        if (basis_nodes.empty()) {
            throw std::runtime_error("E0 coarse space has no free basis nodes");
        }

        const HelmholtzOperators reference_operators = assemble_helmholtz_operators(
            hierarchy.reference_mesh(), kWavenumber);
        const ComplexVector reference_load = assemble_helmholtz_load(
            hierarchy.reference_mesh(), data.source, {}, data.quadrature_context);
        const ComplexVector reference_solution = solve_helmholtz_fem(
            reference_operators, reference_load);
        const HelmholtzError reference_norm = compute_discrete_helmholtz_error(
            hierarchy.reference_mesh(), reference_operators, reference_solution,
            ComplexVector::Zero(reference_solution.size()));
        if (!(reference_norm.energy > 0.0)) {
            throw std::runtime_error("E0 reference solution has zero energy norm");
        }

        const HelmholtzLodModel ideal = build_reference_model(
            data, hierarchy,
            static_cast<int>(hierarchy.coarse_mesh().elems.size()) + 2);
        std::vector<CalibrationRow> rows;
        ComplexVector warm_start;
        for (const int ell : std::array{1, 2, 3, 4}) {
            HelmholtzLodModel model = build_reference_model(data, hierarchy, ell);
            LocalizationEigenConfig eigen_config;
            eigen_config.maximum_iterations = 500;
            eigen_config.relative_tolerance = 2e-12;
            eigen_config.dense_cross_check_max_dimension = 64;
            eigen_config.warm_start = warm_start;
            const ReferenceLocalizationCertificate certificate =
                compute_reference_localization_certificate(
                    hierarchy, model.operators(), ambient_operators,
                    model.corrected_test_basis(), basis_nodes,
                    KernelRieszSolver::SaddlePoint, eigen_config);
            if (!certificate.spectrum.converged) {
                throw std::runtime_error("E0 Theta_loc eigensolver did not converge");
            }
            warm_start = certificate.spectrum.dominant_vector;

            CalibrationRow row;
            row.ell = ell;
            row.theta_loc = certificate.theta_loc;
            if (rows.empty()) {
                row.direct_validation =
                    validate_reference_localization_certificate_small_matrix(
                        hierarchy, model.operators(), ambient_operators,
                        model.corrected_test_basis(), ideal.corrected_test_basis(),
                        basis_nodes, certificate, 2048);
            } else {
                row.direct_validation = rows.front().direct_validation;
                row.direct_validation.direct_delta =
                    compute_reference_localization_direct_delta(
                        hierarchy, model.operators(),
                        model.corrected_test_basis(), ideal.corrected_test_basis(),
                        basis_nodes, certificate);
                row.direct_validation.upper_certificate =
                    row.direct_validation.stable_decomposition_constant
                    / row.direct_validation.ambient_kernel_coercivity
                    * certificate.theta_loc;
                const double tolerance = 2e-9 * std::max({
                    1.0, row.direct_validation.direct_delta,
                    row.direct_validation.upper_certificate});
                row.direct_validation.one_sided_control_holds =
                    row.direct_validation.decomposition_relative_residual <= 2e-10
                    && row.direct_validation.direct_delta
                        <= row.direct_validation.upper_certificate + tolerance;
            }
            row.direct_delta = row.direct_validation.direct_delta;
            row.upper_certificate = row.direct_validation.upper_certificate;

            const HelmholtzLodSolution solution = model.solve_load(reference_load);
            const HelmholtzError error = compute_discrete_helmholtz_error(
                hierarchy.reference_mesh(), reference_operators,
                reference_solution, solution.fine_values);
            const ReferenceResidualRiesz estimate =
                compute_reference_residual_riesz(
                    hierarchy, reference_operators, reference_load,
                    solution.fine_values, kDoerflerTheta,
                    KernelRieszSolver::SaddlePoint);
            row.reference_energy_error = error.energy;
            row.reference_relative_energy_error = error.energy / reference_norm.energy;
            row.eta_H = estimate.eta;
            if (!(row.eta_H > 0.0)) {
                throw std::runtime_error("E0 produced a zero eta_H");
            }
            row.required_multiplier = error.energy / row.eta_H;
            row.allocation_relative_error = estimate.allocation_relative_error;
            for (const LocalKernelRieszResult &local : estimate.local_results) {
                row.maximum_constraint_residual = std::max(
                    row.maximum_constraint_residual,
                    local.constraint_relative_residual);
                row.maximum_riesz_residual = std::max(
                    row.maximum_riesz_residual,
                    local.riesz_relative_residual);
                row.maximum_energy_identity_error = std::max(
                    row.maximum_energy_identity_error,
                    local.energy_identity_relative_error);
            }
            populate_effectivity(
                row, hierarchy, reference_operators, estimate,
                reference_solution, solution.fine_values);
            rows.push_back(row);
            std::cout << "ell=" << ell
                      << " theta_loc=" << numeric(row.theta_loc)
                      << " direct_delta=" << numeric(row.direct_delta)
                      << " upper=" << numeric(row.upper_certificate)
                      << " rel_error=" << numeric(row.reference_relative_energy_error)
                      << " eta_H=" << numeric(row.eta_H)
                      << " required_multiplier="
                      << numeric(row.required_multiplier) << '\n';
        }

        const FrozenParameters frozen = derive_frozen_parameters(rows);
        check_acceptance(rows, frozen, rho_ambient);
        const DriverSmokeEvidence smoke = run_driver_smoke(
            data, hierarchy, rows, frozen);
        write_csv(arguments.output_directory / "e0-calibration.csv", rows);
        write_frozen_json(
            arguments.output_directory / "frozen-parameters.json",
            frozen, rho_ambient, smoke);
        if (arguments.check) {
            for (const char *file : {
                     "e0-calibration.csv", "frozen-parameters.json"}) {
                if (!std::filesystem::is_regular_file(
                        arguments.output_directory / file)) {
                    throw std::runtime_error(
                        std::string("missing E0 artifact: ") + file);
                }
            }
        }
        std::cout << "status=passed\n"
                  << "theta_loc=" << numeric(frozen.theta_loc) << '\n'
                  << "C0_usr=" << numeric(frozen.C0_usr) << '\n'
                  << "C1_usr=" << numeric(frozen.C1_usr) << '\n'
                  << "rho_star=" << numeric(frozen.rho_star) << '\n'
                  << "output=" << arguments.output_directory.string() << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "bench_helmholtz_e0_calibration failed: "
                  << error.what() << '\n';
        return 1;
    }
}

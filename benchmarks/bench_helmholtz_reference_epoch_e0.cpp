#include "helmholtz/adaptive/candidate_dual.h"
#include "helmholtz/adaptive/candidate_flux.h"
#include "helmholtz/adaptive/certificates.h"
#include "helmholtz/adaptive/error_control.h"
#include "helmholtz/adaptive/kernel_residual.h"
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
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace lod2d;
using namespace lod2d::helmholtz;
using namespace lod2d::helmholtz::adaptive;
using namespace lod2d::helmholtz::benchmarks;
using namespace lod2d::helmholtz::experiments;

namespace {

constexpr double kKappa = 16.0;
constexpr int kCoarseLevel = 3;
constexpr int kReferenceLevel = 5;
constexpr double kThetaH = 0.5;
constexpr double kSafety = 1.25;

struct Arguments {
    std::filesystem::path output;
    bool check = false;
};

struct CorrectorRow {
    int ell = 0;
    double theta = 0.0;
    double direct = 0.0;
    double lower = 0.0;
    double upper = 0.0;
};

void require(bool condition, const std::string &message) {
    if (!condition) throw std::runtime_error(message);
}

std::string number(double value) {
    std::ostringstream out;
    out << std::setprecision(17) << value;
    return out.str();
}

Arguments parse_arguments(int argc, char **argv) {
    Arguments result;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument.starts_with("--output-dir="))
            result.output = argument.substr(13);
        else if (argument == "--check") result.check = true;
        else throw std::invalid_argument("unknown argument: " + argument);
    }
    if (result.output.empty())
        throw std::invalid_argument("required: --output-dir=DIR [--check]");
    return result;
}

std::vector<int> free_nodes(const TriMesh &mesh) {
    std::vector<char> fixed(mesh.nodes.size(), false);
    for (int node : dirichlet_nodes(mesh)) fixed[node] = true;
    std::vector<int> result;
    for (int node = 0; node < static_cast<int>(mesh.nodes.size()); ++node)
        if (!fixed[node]) result.push_back(node);
    return result;
}

HelmholtzLodModel build_model(
    const PaperCaseData &data,
    const ReferenceEpochHierarchy &hierarchy,
    int ell) {
    HelmholtzProblemConfig config;
    config.ell = ell;
    config.wavenumber = data.wavenumber;
    config.initial_mesh = data.initial_mesh;
    config.quadrature_context = data.quadrature_context;
    return HelmholtzLodModel::build_adaptive(config, hierarchy);
}

} // namespace

int main(int argc, char **argv) {
    try {
        const Arguments arguments = parse_arguments(argc, argv);
        std::filesystem::create_directories(arguments.output);
        const PaperCaseData data = make_paper_case(PaperCase::R1, kKappa);
        ReferenceEpochHierarchy hierarchy(
            data.initial_mesh, kCoarseLevel, kReferenceLevel);
        const std::vector<int> basis = free_nodes(hierarchy.coarse_mesh());
        const HelmholtzLodModel ideal = build_model(data, hierarchy, 8);

        std::vector<CorrectorRow> corrector_rows;
        ComplexVector warm;
        for (int ell : std::array{1, 2, 3, 4}) {
            const HelmholtzLodModel localized = build_model(data, hierarchy, ell);
            LocalizationEigenConfig eigen;
            eigen.maximum_iterations = 500;
            eigen.relative_tolerance = 2e-11;
            eigen.dense_cross_check_max_dimension = 256;
            eigen.warm_start = warm;
            const ReferenceCorrectorCertificate certificate =
                build_reference_corrector_certificate(
                    hierarchy, localized.operators(),
                    localized.corrected_test_basis(), basis,
                    KernelRieszSolver::SaddlePoint, eigen);
            require(certificate.spectrum.converged,
                    "E0 localization eigensolve did not converge");
            const ReferenceCorrectorDirectValidation direct =
                validate_reference_corrector_certificate_small_matrix(
                    hierarchy, localized.operators(),
                    localized.corrected_test_basis(),
                    ideal.corrected_test_basis(), basis, certificate,
                    8.0, 8.0, 8.0, 0.125, 512);
            require(direct.bracket_holds,
                    "E0 direct corrector defect escaped its bracket");
            corrector_rows.push_back(
                {ell, certificate.theta_loc, direct.direct_delta,
                 direct.theorem_lower, direct.theorem_upper});
            warm = certificate.spectrum.dominant_vector;
        }
        for (const CorrectorRow &row : corrector_rows) {
            require(std::isfinite(row.theta) && row.theta >= 0.0,
                    "E0 Theta_loc is invalid");
            std::cout << "ell=" << row.ell
                      << " Theta_loc=" << number(row.theta)
                      << " direct_delta=" << number(row.direct) << '\n';
        }
        {
            std::ofstream out(arguments.output / "01-corrector-localization.csv");
            out << "schema_version,case,kappa,coarse_level,reference_level,ell,Theta_loc,direct_delta,theorem_lower,theorem_upper,bracket_holds\n";
            for (const CorrectorRow &row : corrector_rows)
                out << "5,R1,16," << kCoarseLevel << ',' << kReferenceLevel
                    << ',' << row.ell << ',' << number(row.theta) << ','
                    << number(row.direct) << ',' << number(row.lower) << ','
                    << number(row.upper) << ",true\n";
        }

        const HelmholtzLodModel calibrated = build_model(data, hierarchy, 2);
        const ComplexVector reference_load = assemble_helmholtz_load(
            hierarchy.reference_mesh(), data.source, {}, data.quadrature_context);
        const ComplexVector reference_solution = solve_helmholtz_fem(
            calibrated.operators(), reference_load);
        const HelmholtzLodSolution lod = calibrated.solve_load(reference_load);
        const ReferenceResidualRiesz residual = compute_reference_residual_riesz(
            hierarchy, calibrated.operators(), reference_load,
            lod.fine_values, kThetaH, KernelRieszSolver::SaddlePoint);
        double maximum_constraint = 0.0;
        double maximum_riesz = 0.0;
        double maximum_energy_identity = 0.0;
        for (const LocalKernelRieszResult &local : residual.local_results) {
            maximum_constraint = std::max(
                maximum_constraint, local.constraint_relative_residual);
            maximum_riesz = std::max(
                maximum_riesz, local.riesz_relative_residual);
            maximum_energy_identity = std::max(
                maximum_energy_identity, local.energy_identity_relative_error);
        }
        require(residual.local_square_sum_relative_error < 2e-12
                    && residual.allocation_relative_error < 2e-12
                    && maximum_constraint < 2e-9
                    && maximum_riesz < 2e-9
                    && maximum_energy_identity < 2e-9,
                "E0 reference residual Riesz identity failed");
        {
            std::ofstream out(arguments.output / "02-reference-riesz.csv");
            out << "schema_version,eta_H,local_square_sum_relative_error,nodal_to_element_allocation_relative_error,max_constraint_residual,max_riesz_residual,max_energy_identity_error,skipped_zero_kernel_patches,status\n"
                << "5," << number(residual.eta) << ','
                << number(residual.local_square_sum_relative_error) << ','
                << number(residual.allocation_relative_error) << ','
                << number(maximum_constraint) << ',' << number(maximum_riesz)
                << ',' << number(maximum_energy_identity) << ','
                << residual.skipped_zero_kernel_patches << ",passed\n";
        }

        ReferenceEpochHierarchy candidate_hierarchy(data.initial_mesh, 2, 3);
        candidate_hierarchy.begin_reference_epoch();
        (void)candidate_hierarchy.enrich_candidate({0, 1});
        const HelmholtzOperators candidate_operators =
            assemble_helmholtz_operators(
                candidate_hierarchy.candidate_mesh(), kKappa);
        const ComplexVector candidate_load = assemble_helmholtz_load(
            candidate_hierarchy.candidate_mesh(), data.source, {},
            data.quadrature_context);
        const ComplexVector candidate_solution = solve_helmholtz_fem(
            candidate_operators, candidate_load);
        CandidateFluxConfig flux_config;
        flux_config.doerfler_theta = kThetaH;
        flux_config.quadrature_context = data.quadrature_context;
        flux_config.compute_discrete_residual_audit = true;
        const CandidateFluxRT2Result flux = reconstruct_candidate_flux_rt2(
            candidate_hierarchy.candidate_mesh(), candidate_operators,
            data.source, candidate_solution, flux_config);
        require(flux.maximum_patch_compatibility_error < 2e-9
                    && flux.maximum_element_divergence_residual < 2e-9
                    && flux.maximum_normal_continuity_residual < 2e-9
                    && flux.maximum_boundary_flux_residual < 2e-9
                    && flux.eta_eq + 2e-9 >= flux.discrete_residual_dual_norm,
                "E0 RT2/P2 identity failed");
        {
            std::ofstream out(arguments.output / "03-candidate-rt2-p2.csv");
            out << "schema_version,eta_eq,discrete_residual_dual_norm,max_patch_compatibility,max_divergence_residual,max_normal_jump,max_boundary_flux_residual,status\n"
                << "5," << number(flux.eta_eq) << ','
                << number(flux.discrete_residual_dual_norm) << ','
                << number(flux.maximum_patch_compatibility_error) << ','
                << number(flux.maximum_element_divergence_residual) << ','
                << number(flux.maximum_normal_continuity_residual) << ','
                << number(flux.maximum_boundary_flux_residual) << ",passed\n";
        }

        const HelmholtzOperators small_reference_operators =
            assemble_helmholtz_operators(
                candidate_hierarchy.reference_mesh(), kKappa);
        const ComplexVector small_reference_load = assemble_helmholtz_load(
            candidate_hierarchy.reference_mesh(), data.source, {},
            data.quadrature_context);
        const ComplexVector small_reference_solution = solve_helmholtz_fem(
            small_reference_operators, small_reference_load);
        const ComplexVector reference_on_candidate =
            candidate_hierarchy.reference_to_candidate().cast<Complex>()
            * small_reference_solution;
        CandidateDualGapConfig gap_config;
        gap_config.continuity_constant = 8.0;
        gap_config.overlap_constant = 8.0;
        gap_config.reference_upper_bound = 0.0;
        gap_config.epoch_switch_ratio = 0.5;
        gap_config.evidence_mode = CandidateGapEvidenceMode::Certified;
        const CandidateDualGapResult gap = build_candidate_dual_gap(
            candidate_hierarchy, candidate_operators, candidate_load,
            reference_on_candidate, gap_config);
        const CandidateDualGapValidation dual_validation =
            validate_candidate_dual_gap_small_mesh(
                candidate_operators, candidate_solution,
                reference_on_candidate, reference_on_candidate, gap);
        require(dual_validation.L_c_lower_bound_holds
                    && dual_validation.L_gap_lower_bound_holds,
                "E0 candidate dual lower bound failed");
        {
            std::ofstream out(arguments.output / "04-candidate-dual-gap.csv");
            out << "schema_version,eta_dual_c,L_c,candidate_error,L_gap_c,reference_candidate_gap,L_c_lower_bound_holds,L_gap_lower_bound_holds,status\n"
                << "5," << number(gap.eta_dual_c) << ',' << number(gap.L_c)
                << ',' << number(dual_validation.candidate_error) << ','
                << number(gap.L_gap_c) << ','
                << number(dual_validation.reference_candidate_gap)
                << ",true,true,passed\n";
        }

        const HelmholtzError reference_norm = compute_discrete_helmholtz_error(
            hierarchy.reference_mesh(), calibrated.operators(),
            reference_solution, ComplexVector::Zero(reference_solution.size()));
        const HelmholtzError lod_error = compute_discrete_helmholtz_error(
            hierarchy.reference_mesh(), calibrated.operators(),
            reference_solution, lod.fine_values);
        const double C_rel_usr = kSafety * lod_error.energy / residual.eta;
        double theta_loc_usr = kSafety * corrector_rows[1].theta;
        if (corrector_rows[0].theta > corrector_rows[1].theta)
            theta_loc_usr = 0.5
                * (corrector_rows[0].theta + corrector_rows[1].theta);
        require(C_rel_usr > 0.0 && std::isfinite(C_rel_usr)
                    && theta_loc_usr >= corrector_rows[1].theta,
                "E0 frozen parameters are invalid");
        {
            std::ofstream out(
                arguments.output / "05-frozen-reference-epoch-parameters.json");
            out << "{\n  \"schema_version\": 5,\n"
                << "  \"scope\": \"E0-localized-smooth-R1-k16\",\n"
                << "  \"claim\": \"implementation-study\",\n"
                << "  \"status\": \"passed\",\n"
                << "  \"theta_loc_usr\": " << number(theta_loc_usr) << ",\n"
                << "  \"C_rel_usr\": " << number(C_rel_usr) << ",\n"
                << "  \"safety_factor\": " << number(kSafety) << ",\n"
                << "  \"observed_relative_reference_energy\": "
                << number(lod_error.energy / reference_norm.energy) << ",\n"
                << "  \"notes\": [\"Frozen once for later experiments; not a rigorous theorem constant.\"]\n}\n";
        }

        if (arguments.check) {
            for (const char *file : {"01-corrector-localization.csv",
                     "02-reference-riesz.csv", "03-candidate-rt2-p2.csv",
                     "04-candidate-dual-gap.csv",
                     "05-frozen-reference-epoch-parameters.json"})
                require(std::filesystem::is_regular_file(arguments.output / file),
                        std::string("missing E0 artifact: ") + file);
        }
        std::cout << "status=passed\n"
                  << "theta_loc_usr=" << number(theta_loc_usr) << '\n'
                  << "C_rel_usr=" << number(C_rel_usr) << '\n'
                  << "output=" << arguments.output.string() << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "bench_helmholtz_reference_epoch_e0 failed: "
                  << error.what() << '\n';
        return 1;
    }
}

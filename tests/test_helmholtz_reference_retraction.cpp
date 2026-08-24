#include "helmholtz/adaptive/reference_retraction.h"
#include "helmholtz/benchmarks/paper_cases.h"
#include "helmholtz/boundary.h"
#include "helmholtz/model.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace lod2d;
using namespace lod2d::helmholtz;
using namespace lod2d::helmholtz::adaptive;
using namespace lod2d::helmholtz::benchmarks;
using namespace lod2d::helmholtz::experiments;

namespace {

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

std::vector<int> free_coarse_nodes(const TriMesh &mesh) {
    std::vector<char> is_dirichlet(mesh.nodes.size(), false);
    for (int node : dirichlet_nodes(mesh)) is_dirichlet[node] = true;
    std::vector<int> result;
    for (int node = 0; node < static_cast<int>(mesh.nodes.size()); ++node) {
        if (!is_dirichlet[node]) result.push_back(node);
    }
    return result;
}

HelmholtzLodModel build_reference_model(
    const PaperCaseData &data,
    const ReferenceEpochHierarchy &hierarchy,
    int ell) {
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

void verify_cached_hierarchy_model_build() {
    const PaperCaseData data = make_paper_case(PaperCase::R2a, 4.0);
    ReferenceEpochHierarchy hierarchy(data.initial_mesh, 1, 3);
    HelmholtzProblemConfig config;
    config.H = 1;
    config.h = 3;
    config.ell = 1;
    config.wavenumber = data.wavenumber;
    config.initial_mesh = data.initial_mesh;
    config.quadrature_context = data.quadrature_context;

    HelmholtzLodModel rebuilt = HelmholtzLodModel::build_adaptive(
        config,
        hierarchy.coarse_mesh(), hierarchy.coarse_levels(),
        hierarchy.reference_mesh(), hierarchy.reference_element_levels());
    HelmholtzLodModel reused = HelmholtzLodModel::build_adaptive(
        config, hierarchy);

    require(rebuilt.problem().coarse_to_fine.isApprox(
                reused.problem().coarse_to_fine, 0.0),
            "cached hierarchy changed the nodal embedding");
    require(rebuilt.problem().fine_element_prolongation.isApprox(
                reused.problem().fine_element_prolongation, 0.0),
            "cached hierarchy changed the element embedding");
    require(rebuilt.problem().fine_dg_prolongation.isApprox(
                reused.problem().fine_dg_prolongation, 0.0),
            "cached hierarchy changed the DG embedding");
    require(rebuilt.problem().quasi_interpolation.isApprox(
                reused.problem().quasi_interpolation, 2e-14),
            "cached hierarchy changed quasi-interpolation");
    require(rebuilt.coarse_operator().isApprox(
                reused.coarse_operator(), 2e-13),
            "cached hierarchy changed the SLOD coarse operator");
}

void verify_retraction_identities() {
    const PaperCaseData data = make_paper_case(PaperCase::R2a, 4.0);
    ReferenceEpochHierarchy hierarchy(data.initial_mesh, 1, 3);
    const AmbientRatioEnforcementResult ambient_update =
        hierarchy.enforce_ambient_ratio(0.2);
    require(ambient_update.changed,
            "retraction test did not create a distinct ambient mesh");

    const ReferenceRetraction retraction =
        build_reference_retraction(hierarchy);
    const ReferenceRetractionDiagnostics &diagnostics =
        retraction.diagnostics();
    require(diagnostics.projector_identity_relative_error <= 2e-13,
            "J_ref is not the identity on the embedded reference space");
    require(diagnostics.reference_kernel_identity_relative_error <= 2e-13,
            "R_ref is not the identity on the embedded reference kernel");
    require(diagnostics.kernel_constraint_relative_error <= 2e-13,
            "the range of R_ref is not contained in W_ref");
    require(diagnostics.maximum_column_nonzeros > 0,
            "reference retraction has empty support");

    ComplexVector reference_values(hierarchy.reference_mesh().nodes.size());
    for (int node = 0; node < reference_values.size(); ++node) {
        reference_values(node) = Complex(
            std::sin(0.37 * (node + 1)),
            std::cos(0.23 * (node + 1)));
    }
    for (int node : dirichlet_nodes(hierarchy.reference_mesh()))
        reference_values(node) = Complex(0.0, 0.0);
    const ComplexVector reference_kernel = reference_values
        - hierarchy.coarse_to_reference().cast<Complex>()
            * (hierarchy.reference_quasi_interpolation().cast<Complex>()
               * reference_values);
    const ComplexVector embedded =
        hierarchy.reference_to_ambient().cast<Complex>() * reference_kernel;
    const ComplexVector recovered = retraction.apply(embedded);
    require((recovered - reference_kernel).norm()
                <= 2e-12 * std::max(1.0, reference_kernel.norm()),
            "R_ref did not reproduce a reference-kernel vector");

    ComplexVector arbitrary_ambient(hierarchy.ambient_mesh().nodes.size());
    for (int node = 0; node < arbitrary_ambient.size(); ++node) {
        arbitrary_ambient(node) = Complex(
            std::cos(0.19 * (node + 1)),
            std::sin(0.31 * (node + 1)));
    }
    for (int node : dirichlet_nodes(hierarchy.ambient_mesh()))
        arbitrary_ambient(node) = Complex(0.0, 0.0);
    const ComplexVector retracted = retraction.apply(arbitrary_ambient);
    const ComplexVector constraint =
        hierarchy.reference_quasi_interpolation().cast<Complex>() * retracted;
    require(constraint.norm()
                <= 2e-12 * std::max(1.0, retracted.norm()),
            "an arbitrary ambient input was not retracted into W_ref");
    for (int node : dirichlet_nodes(hierarchy.reference_mesh())) {
        require(std::abs(retracted(node)) <= 2e-13,
                "reference retraction did not preserve homogeneous trace");
    }
}

void verify_clustered_sparse_generalized_spectrum() {
    constexpr int dimension = 64;
    Eigen::SparseMatrix<double> denominator(dimension, dimension);
    std::vector<Eigen::Triplet<double>> denominator_entries;
    denominator_entries.reserve(dimension);
    ComplexMatrix numerator = ComplexMatrix::Zero(dimension, dimension);
    for (int index = 0; index < dimension; ++index) {
        const double energy_weight = 1.0 + 0.01 * index;
        denominator_entries.emplace_back(index, index, energy_weight);
        const double generalized_eigenvalue = index == 0
            ? 1.0
            : (index == 1 ? 1.0 - 1e-8 : 0.1);
        numerator(index, index) = energy_weight * generalized_eigenvalue;
    }
    denominator.setFromTriplets(
        denominator_entries.begin(), denominator_entries.end());

    LocalizationEigenConfig config;
    config.maximum_iterations = 80;
    config.relative_tolerance = 1e-12;
    config.dense_cross_check_max_dimension = dimension;
    config.dense_fallback_max_dimension = 0;
    config.sparse_generalized_min_dimension = 0;
    const LocalizationSpectrum spectrum = compute_localization_spectrum(
        numerator, denominator, config);
    require(spectrum.used_sparse_generalized_solver
                && spectrum.converged
                && spectrum.relative_residual <= config.relative_tolerance
                && spectrum.dense_cross_checked
                && spectrum.dense_relative_difference <= 2e-12
                && std::abs(spectrum.lambda_max - 1.0) <= 2e-12
                && spectrum.iterations < config.maximum_iterations,
            "preconditioned Ritz iteration failed on a clustered dominant spectrum");
}

void verify_localization_certificate() {
    const PaperCaseData data = make_paper_case(PaperCase::R2a, 4.0);
    ReferenceEpochHierarchy hierarchy(data.initial_mesh, 4, 6);
    const AmbientRatioEnforcementResult ambient_update =
        hierarchy.enforce_ambient_ratio(0.2);
    require(ambient_update.changed,
            "certificate test did not create a distinct ambient mesh");
    const HelmholtzOperators ambient_operators = assemble_helmholtz_operators(
        hierarchy.ambient_mesh(), data.wavenumber);
    const std::vector<int> basis_nodes = free_coarse_nodes(
        hierarchy.coarse_mesh());
    require(!basis_nodes.empty(),
            "certificate test has no free coarse basis functions");

    std::vector<double> theta;
    std::vector<ReferenceLocalizationCertificate> certificates;
    ComplexVector warm_start;
    for (int ell : {1, 2, 3}) {
        const HelmholtzLodModel model = build_reference_model(
            data, hierarchy, ell);
        LocalizationEigenConfig eigen_config;
        eigen_config.maximum_iterations = 400;
        eigen_config.relative_tolerance = 2e-12;
        eigen_config.dense_cross_check_max_dimension = 64;
        eigen_config.warm_start = warm_start;
        ReferenceLocalizationCertificate certificate =
            compute_reference_localization_certificate(
                hierarchy,
                model.operators(),
                ambient_operators,
                model.corrected_test_basis(),
                basis_nodes,
                KernelRieszSolver::SaddlePoint,
                eigen_config);
        require(certificate.spectrum.converged,
                "Theta_loc eigensolver did not converge");
        require(certificate.status
                    == LocalizationCertificateStatus::ImplementationStudy,
                "ordinary floating-point WP3 result claimed rigorous certification");
        require(certificate.spectrum.relative_residual <= 3e-11,
                "Theta_loc eigen residual is too large");
        require(certificate.spectrum.dense_cross_checked
                    && certificate.spectrum.dense_relative_difference <= 3e-10,
                "iterative Theta_loc disagrees with the dense cross-check");
        require(certificate.coarse_energy_operator.rows()
                    == certificate.coarse_energy.rows()
                    && ComplexMatrix(
                        certificate.coarse_energy_operator.cast<Complex>())
                        .isApprox(certificate.coarse_energy, 2e-14),
                "sparse and retained dense coarse energy matrices disagree");
        require(certificate.spectrum.used_warm_start == (ell > 1),
                "localization eigensolver warm-start accounting is wrong");
        require(certificate.ambient_riesz.space
                    == KernelRieszSpace::AmbientDefect,
                "localization certificate did not use ambient defect Riesz");
        require(certificate.ambient_riesz.local_square_sum_relative_error
                    <= 3e-11,
                "certificate ambient local square sum is inconsistent");
        require(!certificate.ambient_riesz.local_details_stored
                    && certificate.ambient_riesz.patches.empty()
                    && certificate.ambient_riesz.local_results.empty()
                    && certificate.ambient_riesz.patch_factorizations
                        <= certificate.ambient_riesz.patch_count
                    && certificate.ambient_riesz.right_hand_side_solves
                        <= certificate.ambient_riesz.patch_count
                            * basis_nodes.size()
                    && certificate.ambient_riesz.maximum_active_columns
                        <= basis_nodes.size(),
                "production ambient defect Riesz retained details or exceeded dense work");
        require(certificate.ambient_riesz.patch_solve_seconds > 0.0
                    && certificate.ambient_riesz.gram_reduction_seconds >= 0.0,
                "ambient defect Riesz stage timings were not recorded");
        const double timed_certificate_components =
            certificate.timings.retraction_seconds
            + certificate.timings.defect_rhs_seconds
            + certificate.timings.ambient_riesz_seconds
            + certificate.timings.coarse_energy_seconds
            + certificate.timings.spectrum_seconds;
        require(certificate.timings.retraction_seconds > 0.0
                    && certificate.timings.defect_rhs_seconds > 0.0
                    && certificate.timings.ambient_riesz_seconds > 0.0
                    && certificate.timings.coarse_energy_seconds > 0.0
                    && certificate.timings.spectrum_seconds > 0.0
                    && certificate.timings.total_seconds > 0.0
                    && timed_certificate_components
                        <= certificate.timings.total_seconds * (1.0 + 1e-9),
                "localization certificate stage timings are incomplete or inconsistent");
        require(certificate.retraction.diagnostics()
                    .kernel_constraint_relative_error <= 2e-13,
                "certificate used a retraction outside W_ref");
        theta.push_back(certificate.theta_loc);
        std::cout << "ell=" << ell
                  << " theta_loc=" << certificate.theta_loc
                  << " defect_norm=" << certificate.defect_rhs.norm()
                  << '\n';
        warm_start = certificate.spectrum.dominant_vector;
        certificates.push_back(std::move(certificate));
    }
    require(theta.front() > 1e-12,
            "localized corrector produced a trivial initial certificate");
    for (int index = 1; index < static_cast<int>(theta.size()); ++index) {
        require(theta[index]
                    <= theta[index - 1] * (1.0 + 2e-8) + 2e-11,
                "Theta_loc did not decrease when ell increased");
    }

    const HelmholtzLodModel fallback_model = build_reference_model(
        data, hierarchy, 1);
    LocalizationEigenConfig fallback_config;
    fallback_config.maximum_iterations = 1;
    fallback_config.relative_tolerance = 1e-10;
    fallback_config.dense_cross_check_max_dimension = 0;
    fallback_config.dense_fallback_max_dimension = 512;
    const ReferenceLocalizationCertificate fallback_certificate =
        compute_reference_localization_certificate(
            hierarchy,
            fallback_model.operators(),
            ambient_operators,
            fallback_model.corrected_test_basis(),
            basis_nodes,
            KernelRieszSolver::SaddlePoint,
            fallback_config);
    require(fallback_certificate.spectrum.used_dense_fallback
                && fallback_certificate.spectrum.converged
                && fallback_certificate.spectrum.relative_residual <= 1e-10
                && std::abs(
                    fallback_certificate.theta_loc - theta.front())
                    <= 3e-10 * std::max(1.0, theta.front()),
            "dense localization fallback did not recover the same spectrum");

    LocalizationEigenConfig sparse_config;
    sparse_config.maximum_iterations = 2000;
    sparse_config.relative_tolerance = 1e-10;
    sparse_config.dense_cross_check_max_dimension = 512;
    sparse_config.sparse_generalized_min_dimension = 0;
    sparse_config.warm_start = certificates.front().spectrum.dominant_vector;
    const ReferenceLocalizationCertificate sparse_certificate =
        compute_reference_localization_certificate(
            hierarchy,
            fallback_model.operators(),
            ambient_operators,
            fallback_model.corrected_test_basis(),
            basis_nodes,
            KernelRieszSolver::SaddlePoint,
            sparse_config);
    require(sparse_certificate.spectrum.used_sparse_generalized_solver
                && sparse_certificate.spectrum.used_warm_start
                && sparse_certificate.spectrum.converged
                && sparse_certificate.spectrum.dense_cross_checked
                && sparse_certificate.spectrum.dense_relative_difference <= 3e-9
                && std::abs(sparse_certificate.theta_loc - theta.front())
                    <= 3e-9 * std::max(1.0, theta.front()),
            "sparse generalized localization spectrum disagrees with dense validation");

    fallback_config.dense_fallback_max_dimension = 0;
    bool disabled_fallback_rejected = false;
    try {
        (void)compute_reference_localization_certificate(
            hierarchy,
            fallback_model.operators(),
            ambient_operators,
            fallback_model.corrected_test_basis(),
            basis_nodes,
            KernelRieszSolver::SaddlePoint,
            fallback_config);
    } catch (const std::runtime_error &) {
        disabled_fallback_rejected = true;
    }
    require(disabled_fallback_rejected,
            "disabled dense fallback accepted a nonconverged eigen iteration");
    require(LocalizationEigenConfig{}.maximum_iterations == 1000,
            "production localization iteration budget regressed");
    require(LocalizationEigenConfig{}.dense_fallback_max_dimension == 1024,
            "production dense localization fallback range regressed");
    require(LocalizationEigenConfig{}.sparse_generalized_min_dimension == 1025,
            "production sparse generalized threshold regressed");

    const HelmholtzLodModel localized = build_reference_model(
        data, hierarchy, 1);
    const HelmholtzLodModel ideal = build_reference_model(
        data, hierarchy,
        static_cast<int>(hierarchy.coarse_mesh().elems.size()) + 2);
    HelmholtzOperators stale_ambient_operators = ambient_operators;
    stale_ambient_operators.mass.coeffRef(0, 0) += 1e-5;
    bool stale_operators_rejected = false;
    try {
        (void)compute_reference_localization_certificate(
            hierarchy,
            localized.operators(),
            stale_ambient_operators,
            localized.corrected_test_basis(),
            basis_nodes);
    } catch (const std::invalid_argument &) {
        stale_operators_rejected = true;
    }
    require(stale_operators_rejected,
            "localization certificate accepted stale ambient operators");
    const SmallMatrixLocalizationValidation validation =
        validate_reference_localization_certificate_small_matrix(
            hierarchy,
            localized.operators(),
            ambient_operators,
            localized.corrected_test_basis(),
            ideal.corrected_test_basis(),
            basis_nodes,
            certificates.front(),
            2048);
    require(validation.ambient_kernel_dimension > 0,
            "small-matrix ambient kernel is unexpectedly trivial");
    require(validation.local_kernel_columns
                >= validation.ambient_kernel_dimension,
            "local ambient kernels cannot span the global kernel");
    require(validation.decomposition_relative_residual <= 2e-10,
            "ambient local kernels do not span the global kernel");
    require(validation.ambient_kernel_coercivity > 0.0,
            "small-matrix ambient kernel is not coercive");
    require(std::isfinite(validation.projector_stability_constant)
                && validation.projector_stability_constant > 0.0
                && validation.projector_stability_constant < 20.0,
            "local J_ref energy-stability diagnostic is invalid");
    require(std::isfinite(validation.retraction_stability_constant)
                && validation.retraction_stability_constant > 0.0
                && validation.retraction_stability_constant < 50.0,
            "R_ref energy-stability diagnostic is invalid");
    require(validation.direct_delta > 0.0,
            "ideal/localized reference correctors unexpectedly agree");
    require(validation.one_sided_control_holds,
            "retracted ambient certificate failed the theorem upper direction");
    std::cout << "direct_delta=" << validation.direct_delta
              << " upper_certificate=" << validation.upper_certificate
              << " C_J=" << validation.projector_stability_constant
              << " C_ret=" << validation.retraction_stability_constant
              << " c_W=" << validation.ambient_kernel_coercivity
              << " C_sd=" << validation.stable_decomposition_constant
              << '\n';

    ComplexSparseMatrix stale_basis = localized.corrected_test_basis();
    stale_basis *= Complex(1.0 + 1e-5, 0.0);
    bool stale_rejected = false;
    try {
        (void)validate_reference_localization_certificate_small_matrix(
            hierarchy,
            localized.operators(),
            ambient_operators,
            stale_basis,
            ideal.corrected_test_basis(),
            basis_nodes,
            certificates.front(),
            2048);
    } catch (const std::invalid_argument &) {
        stale_rejected = true;
    }
    require(stale_rejected,
            "small-matrix validation accepted a stale localization certificate");

    const ReferenceLocalizationCertificate ideal_certificate =
        compute_reference_localization_certificate(
            hierarchy,
            ideal.operators(),
            ambient_operators,
            ideal.corrected_test_basis(),
            basis_nodes);
    require(ideal_certificate.theta_loc
                <= 2e-8 * std::max(1.0, certificates.front().theta_loc),
            "retracted defect did not vanish for the ideal reference corrector");
}

} // namespace

int main() {
    try {
        verify_cached_hierarchy_model_build();
        verify_retraction_identities();
        verify_clustered_sparse_generalized_spectrum();
        verify_localization_certificate();
        std::cout << "Ambient-to-reference retraction and certificate passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_reference_retraction failed: "
                  << error.what() << '\n';
        return 1;
    }
}

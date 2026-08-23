#include "helmholtz/adaptive/certificates.h"
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
    std::vector<char> constrained(mesh.nodes.size(), false);
    for (int node : dirichlet_nodes(mesh)) constrained[node] = true;
    std::vector<int> result;
    for (int node = 0; node < static_cast<int>(mesh.nodes.size()); ++node)
        if (!constrained[node]) result.push_back(node);
    return result;
}

HelmholtzLodModel build_model(
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
    return HelmholtzLodModel::build_adaptive(config, hierarchy);
}

} // namespace

int main() {
    try {
        const PaperCaseData data = make_paper_case(PaperCase::R1, 16.0);
        ReferenceEpochHierarchy hierarchy(data.initial_mesh, 3, 5);
        const std::vector<int> basis_nodes = free_coarse_nodes(
            hierarchy.coarse_mesh());
        const HelmholtzLodModel ideal = build_model(data, hierarchy, 8);

        std::vector<double> theta;
        std::vector<double> direct;
        ComplexVector warm_start;
        for (int ell : {1, 2, 3}) {
            const HelmholtzLodModel localized = build_model(
                data, hierarchy, ell);
            LocalizationEigenConfig eigen;
            eigen.maximum_iterations = 300;
            eigen.relative_tolerance = 2e-11;
            eigen.dense_cross_check_max_dimension = 256;
            eigen.warm_start = warm_start;
            const ReferenceCorrectorCertificate certificate =
                build_reference_corrector_certificate(
                    hierarchy, localized.operators(),
                    localized.corrected_test_basis(), basis_nodes,
                    KernelRieszSolver::SaddlePoint, eigen);
            require(certificate.status
                        == ReferenceCorrectorCertificateStatus::ImplementationStudy,
                    "floating-point reference certificate claimed rigorous status");
            require(certificate.matrix_free_action_used,
                    "reference certificate did not use the matrix-free action");
            require(certificate.spectrum.converged
                        && certificate.spectrum.dense_cross_checked,
                    "reference certificate eigensolve was not cross-checked");
            require(certificate.defect_rhs.rows()
                        == static_cast<int>(hierarchy.reference_mesh().nodes.size()),
                    "reference defect left W_h dimensions");
            const ReferenceCorrectorDirectValidation validation =
                validate_reference_corrector_certificate_small_matrix(
                    hierarchy, localized.operators(),
                    localized.corrected_test_basis(),
                    ideal.corrected_test_basis(), basis_nodes, certificate,
                    8.0, 8.0, 8.0, 0.125, 512);
            require(validation.bracket_holds,
                    "direct corrector defect is outside the theorem bracket");
            require(certificate.theta_loc_diagnostic_lower
                        <= certificate.theta_loc
                    && certificate.theta_loc
                        <= certificate.theta_loc_diagnostic_upper,
                    "Theta_loc is outside its algebraic diagnostic interval");
            if (ell == 1) {
                LocalizationEigenConfig block_eigen = eigen;
                block_eigen.relative_tolerance = 1e-8;
                block_eigen.dense_cross_check_max_dimension = 0;
                block_eigen.dense_fallback_max_dimension = 0;
                block_eigen.warm_start.resize(0);
                const ReferenceCorrectorCertificate block_certificate =
                    build_reference_corrector_certificate(
                        hierarchy, localized.operators(),
                        localized.corrected_test_basis(), basis_nodes,
                        KernelRieszSolver::SaddlePoint, block_eigen);
                require(block_certificate.spectrum.converged
                            && !block_certificate.spectrum.dense_cross_checked,
                        "matrix-free block eigensolve did not converge");
                require(std::abs(block_certificate.theta_loc
                                     - certificate.theta_loc)
                            <= 1e-7 * std::max(1.0, certificate.theta_loc),
                        "matrix-free block eigensolve disagrees with dense spectrum");
            }
            if (ell == 1) {
                ComplexVector probe(basis_nodes.size());
                for (int index = 0; index < probe.size(); ++index)
                    probe(index) = Complex(0.25 * (index + 1), -0.125 * index);
                ReferenceDefectGramOperator serial(
                    hierarchy, localized.operators(), certificate.defect_rhs,
                    KernelRieszSolver::SaddlePoint, 1);
                ReferenceDefectGramOperator parallel(
                    hierarchy, localized.operators(), certificate.defect_rhs,
                    KernelRieszSolver::SaddlePoint, 4);
                ReferenceDefectGramFactorCache factor_cache(1024);
                ReferenceDefectGramOperator cached_first(
                    hierarchy, localized.operators(), certificate.defect_rhs,
                    KernelRieszSolver::SaddlePoint, 1, &factor_cache);
                ReferenceDefectGramOperator cached_second(
                    hierarchy, localized.operators(), certificate.defect_rhs,
                    KernelRieszSolver::SaddlePoint, 1, &factor_cache);
                const ComplexVector serial_first = serial.apply(probe);
                const ComplexVector serial_second = serial.apply(probe);
                const ComplexVector parallel_value = parallel.apply(probe);
                ComplexMatrix probe_block(probe.size(), 2);
                probe_block.col(0) = probe;
                probe_block.col(1) = Complex(0.25, -0.5) * probe;
                const ComplexMatrix parallel_block = parallel.apply(probe_block);
                const ReferenceDefectRiesz assembled =
                    compute_reference_defect_riesz(
                        hierarchy, localized.operators(), certificate.defect_rhs,
                        KernelRieszSolver::SaddlePoint,
                        AmbientDefectDetail::SummaryOnly, 4);
                const ComplexVector assembled_value = assembled.gram * probe;
                require((serial_first - serial_second).norm() <= 1e-14
                            * std::max(1.0, serial_first.norm()),
                        "prepared Gram action is not repeatable");
                require((serial_first - parallel_value).norm() <= 1e-13
                            * std::max(1.0, serial_first.norm()),
                        "parallel prepared Gram action changed the mathematics");
                require((parallel_block.col(0) - serial_first).norm() <= 1e-13
                            * std::max(1.0, serial_first.norm())
                            && (parallel_block.col(1)
                                - Complex(0.25, -0.5) * serial_first).norm()
                                <= 1e-13 * std::max(1.0, serial_first.norm()),
                        "block prepared Gram action changed the mathematics");
                require((serial_first - assembled_value).norm() <= 1e-12
                            * std::max(1.0, assembled_value.norm()),
                        "prepared Gram action disagrees with assembled Gram");
                require(cached_first.diagnostics().factor_cache_misses
                                == cached_first.diagnostics().patch_count
                            && cached_second.diagnostics().factor_cache_hits
                                == cached_second.diagnostics().patch_count
                            && cached_second.diagnostics().patch_factorizations
                                == 0
                            && (cached_second.apply(probe) - serial_first).norm()
                                <= 1e-13 * std::max(1.0, serial_first.norm()),
                        "cross-check Gram factor LRU missed or changed the action");
                require(serial.diagnostics().action_calls == 2
                            && serial.diagnostics().patch_factorizations
                                == serial.diagnostics().patch_count,
                        "prepared Gram action repeated or missed a patch factorization");
            }
            theta.push_back(certificate.theta_loc);
            direct.push_back(validation.direct_delta);
            warm_start = certificate.spectrum.dominant_vector;
            std::cout << "ell=" << ell
                      << " theta_loc=" << certificate.theta_loc
                      << " direct_delta=" << validation.direct_delta
                      << " lower=" << validation.theorem_lower
                      << " upper=" << validation.theorem_upper << '\n';
        }
        require(theta[1] <= theta[0] * (1.0 + 2e-8)
                    && theta[2] <= theta[1] * (1.0 + 2e-8),
                "reference Theta_loc did not decrease with ell");
        require(direct[1] <= direct[0] * (1.0 + 2e-8)
                    && direct[2] <= direct[1] * (1.0 + 2e-8),
                "direct reference corrector defect did not decrease with ell");
        std::cout << "Reference-kernel corrector certificate passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_reference_corrector_certificate failed: "
                  << error.what() << '\n';
        return 1;
    }
}

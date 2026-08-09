#include "helmholtz/adaptive/certificates.h"

#include "helmholtz/boundary.h"
#include "helmholtz/patch_system.h"
#include "lod/patches.h"

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace lod2d::helmholtz::adaptive {
namespace {

using Enclosure = solver::ComplexMatrixEnclosure;

constexpr double kInfinity = std::numeric_limits<double>::infinity();

struct RieszBlockResult {
    Enclosure solution;
    Enclosure gram;
    RieszCertificateDiagnostics diagnostics;
};

struct FineAssemblyResult {
    Enclosure primal;
    Enclosure adjoint;
    std::vector<Enclosure> element_primal;
    std::vector<Enclosure> element_adjoint;
    RieszCertificateDiagnostics primal_diagnostics;
    RieszCertificateDiagnostics adjoint_diagnostics;
};

bool direction_satisfies(CertificateBoundDirection supplied,
                         CertificateBoundDirection required) {
    return supplied == CertificateBoundDirection::Exact || supplied == required;
}

const std::vector<std::pair<std::string, CertificateBoundDirection>> &
required_constants() {
    static const std::vector<std::pair<std::string, CertificateBoundDirection>> values = {
        {"C_app", CertificateBoundDirection::Upper},
        {"C_st", CertificateBoundDirection::Upper},
        {"C_sd", CertificateBoundDirection::Upper},
        {"C_ov", CertificateBoundDirection::Upper},
        {"C_a", CertificateBoundDirection::Upper},
        {"C_Fort", CertificateBoundDirection::Upper},
        {"c_W", CertificateBoundDirection::Lower},
        {"C_Pi", CertificateBoundDirection::Upper},
        {"C_loc", CertificateBoundDirection::Upper},
        {"C_ol(ell)", CertificateBoundDirection::Upper},
        {"C_ol(ell+s)", CertificateBoundDirection::Upper},
        {"beta", CertificateBoundDirection::Upper},
        {"s", CertificateBoundDirection::Upper}};
    return values;
}

double triangle_diameter(const TriMesh &mesh, const Triangle &triangle) {
    double diameter = 0.0;
    for (int first = 0; first < 3; ++first) {
        for (int second = first + 1; second < 3; ++second) {
            diameter = std::max(
                diameter,
                (mesh.nodes[triangle[first]] - mesh.nodes[triangle[second]]).norm());
        }
    }
    return diameter;
}

double maximum_coarse_diameter(const TriMesh &mesh) {
    double diameter = 0.0;
    for (const Triangle &triangle : mesh.elems)
        diameter = std::max(diameter, triangle_diameter(mesh, triangle));
    return diameter;
}

double exact_patch_overlap(const TriMesh &mesh, int layers) {
    if (layers < 0) throw std::invalid_argument("patch overlap layers must be nonnegative");
    const Eigen::SparseMatrix<double> patches = build_patches(mesh, layers);
    std::vector<int> multiplicity(mesh.elems.size(), 0);
    for (int target = 0; target < patches.outerSize(); ++target) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(patches, target); it; ++it) {
            if (it.value() != 0.0) ++multiplicity[it.row()];
        }
    }
    return static_cast<double>(
        *std::max_element(multiplicity.begin(), multiplicity.end()));
}

double rounding_gamma(Eigen::Index operation_count) {
    const double unit_roundoff = 0.5 * std::numeric_limits<double>::epsilon();
    const double count = static_cast<double>(std::max<Eigen::Index>(1, operation_count));
    const double product = count * unit_roundoff;
    if (product >= 0.5) return kInfinity;
    return product / (1.0 - product);
}

double outward_nonnegative(double value) {
    if (!(value >= 0.0) || !std::isfinite(value)) return value;
    return std::nextafter(value, kInfinity);
}

Enclosure point_enclosure(const ComplexMatrix &matrix) {
    Enclosure result;
    result.midpoint = matrix;
    result.radius = Eigen::MatrixXd::Zero(matrix.rows(), matrix.cols());
    result.entries_verified = true;
    return result;
}

Enclosure uniform_enclosure(const ComplexMatrix &matrix, double radius,
                            bool verified) {
    if (!std::isfinite(radius) || radius < 0.0)
        throw std::invalid_argument("uniform matrix enclosure radius is invalid");
    Enclosure result;
    result.midpoint = matrix;
    result.radius = Eigen::MatrixXd::Constant(matrix.rows(), matrix.cols(), radius);
    result.entries_verified = verified;
    return result;
}

void validate_enclosure_dimensions(const Enclosure &matrix) {
    if (matrix.radius.rows() != matrix.midpoint.rows() ||
        matrix.radius.cols() != matrix.midpoint.cols())
        throw std::invalid_argument("matrix enclosure dimensions are inconsistent");
}

Enclosure enclosure_adjoint(const Enclosure &matrix) {
    validate_enclosure_dimensions(matrix);
    Enclosure result;
    result.midpoint = matrix.midpoint.adjoint();
    result.radius = matrix.radius.transpose();
    result.entries_verified = matrix.entries_verified;
    return result;
}

Enclosure enclosure_conjugate(const Enclosure &matrix) {
    validate_enclosure_dimensions(matrix);
    Enclosure result;
    result.midpoint = matrix.midpoint.conjugate();
    result.radius = matrix.radius;
    result.entries_verified = matrix.entries_verified;
    return result;
}

Enclosure enclosure_add(const Enclosure &left, const Enclosure &right,
                        double right_sign = 1.0) {
    validate_enclosure_dimensions(left);
    validate_enclosure_dimensions(right);
    if (left.midpoint.rows() != right.midpoint.rows() ||
        left.midpoint.cols() != right.midpoint.cols())
        throw std::invalid_argument("matrix enclosure addition dimensions differ");
    Enclosure result;
    result.midpoint = left.midpoint + right_sign * right.midpoint;
    result.radius.resize(left.midpoint.rows(), left.midpoint.cols());
    const double gamma = rounding_gamma(2);
    for (Eigen::Index i = 0; i < result.midpoint.rows(); ++i) {
        for (Eigen::Index j = 0; j < result.midpoint.cols(); ++j) {
            const double magnitude =
                std::abs(left.midpoint(i, j).real()) +
                std::abs(left.midpoint(i, j).imag()) +
                std::abs(right.midpoint(i, j).real()) +
                std::abs(right.midpoint(i, j).imag());
            result.radius(i, j) = outward_nonnegative(
                left.radius(i, j) + right.radius(i, j) + gamma * magnitude);
        }
    }
    result.entries_verified = left.entries_verified && right.entries_verified;
    return result;
}

Enclosure enclosure_multiply(const Enclosure &left, const Enclosure &right) {
    validate_enclosure_dimensions(left);
    validate_enclosure_dimensions(right);
    if (left.midpoint.cols() != right.midpoint.rows())
        throw std::invalid_argument("matrix enclosure multiplication dimensions differ");
    Enclosure result;
    result.midpoint = left.midpoint * right.midpoint;
    result.radius = Eigen::MatrixXd::Zero(
        left.midpoint.rows(), right.midpoint.cols());
    const double gamma = rounding_gamma(4 * left.midpoint.cols() + 2);
    for (Eigen::Index i = 0; i < left.midpoint.rows(); ++i) {
        for (Eigen::Index j = 0; j < right.midpoint.cols(); ++j) {
            double propagated = 0.0;
            double product_magnitude = 0.0;
            for (Eigen::Index k = 0; k < left.midpoint.cols(); ++k) {
                const double left_l1 =
                    std::abs(left.midpoint(i, k).real()) +
                    std::abs(left.midpoint(i, k).imag());
                const double right_l1 =
                    std::abs(right.midpoint(k, j).real()) +
                    std::abs(right.midpoint(k, j).imag());
                const double left_radius = left.radius(i, k);
                const double right_radius = right.radius(k, j);
                product_magnitude += left_l1 * right_l1;
                propagated += left_l1 * 2.0 * right_radius
                    + 2.0 * left_radius * right_l1
                    + 4.0 * left_radius * right_radius;
            }
            result.radius(i, j) = outward_nonnegative(
                propagated + gamma * product_magnitude);
        }
    }
    result.entries_verified = left.entries_verified && right.entries_verified;
    return result;
}

Enclosure restrict_rows(const Enclosure &matrix, const std::vector<int> &rows) {
    validate_enclosure_dimensions(matrix);
    Enclosure result;
    result.midpoint.resize(rows.size(), matrix.midpoint.cols());
    result.radius.resize(rows.size(), matrix.midpoint.cols());
    for (int local = 0; local < static_cast<int>(rows.size()); ++local) {
        if (rows[local] < 0 || rows[local] >= matrix.midpoint.rows())
            throw std::out_of_range("matrix enclosure row restriction is invalid");
        result.midpoint.row(local) = matrix.midpoint.row(rows[local]);
        result.radius.row(local) = matrix.radius.row(rows[local]);
    }
    result.entries_verified = matrix.entries_verified;
    return result;
}

Enclosure select_columns(const Enclosure &matrix, const std::vector<int> &columns) {
    validate_enclosure_dimensions(matrix);
    Enclosure result;
    result.midpoint.resize(matrix.midpoint.rows(), columns.size());
    result.radius.resize(matrix.midpoint.rows(), columns.size());
    for (int local = 0; local < static_cast<int>(columns.size()); ++local) {
        if (columns[local] < 0 || columns[local] >= matrix.midpoint.cols())
            throw std::out_of_range("matrix enclosure column selection is invalid");
        result.midpoint.col(local) = matrix.midpoint.col(columns[local]);
        result.radius.col(local) = matrix.radius.col(columns[local]);
    }
    result.entries_verified = matrix.entries_verified;
    return result;
}

Enclosure restrict_square(const Enclosure &matrix, const std::vector<int> &dofs) {
    return select_columns(restrict_rows(matrix, dofs), dofs);
}

Enclosure hermitianize(const Enclosure &matrix) {
    if (matrix.midpoint.rows() != matrix.midpoint.cols())
        throw std::invalid_argument("Hermitian enclosure must be square");
    Enclosure result;
    result.midpoint = 0.5 * (matrix.midpoint + matrix.midpoint.adjoint());
    result.radius.resize(matrix.midpoint.rows(), matrix.midpoint.cols());
    const double gamma = rounding_gamma(4);
    for (Eigen::Index i = 0; i < matrix.midpoint.rows(); ++i) {
        for (Eigen::Index j = 0; j < matrix.midpoint.cols(); ++j) {
            const Complex mismatch = matrix.midpoint(i, j)
                - std::conj(matrix.midpoint(j, i));
            const double recenter = 0.5 * std::max(
                std::abs(mismatch.real()), std::abs(mismatch.imag()));
            const double magnitude =
                std::abs(matrix.midpoint(i, j).real()) +
                std::abs(matrix.midpoint(i, j).imag()) +
                std::abs(matrix.midpoint(j, i).real()) +
                std::abs(matrix.midpoint(j, i).imag());
            result.radius(i, j) = outward_nonnegative(
                std::max(matrix.radius(i, j), matrix.radius(j, i))
                + recenter + gamma * magnitude);
        }
    }
    result.entries_verified = matrix.entries_verified;
    return result;
}

double vector_enclosure_norm_upper(const Enclosure &matrix, Eigen::Index column) {
    long double squared = 0.0L;
    for (Eigen::Index row = 0; row < matrix.midpoint.rows(); ++row) {
        const double magnitude = std::abs(matrix.midpoint(row, column))
            + std::sqrt(2.0) * matrix.radius(row, column);
        squared += static_cast<long double>(magnitude)
            * static_cast<long double>(magnitude);
    }
    return outward_nonnegative(std::sqrt(static_cast<double>(squared)));
}

RieszBlockResult solve_constrained_riesz(
    const Enclosure &energy,
    const Eigen::MatrixXd &constraints,
    const Enclosure &rhs,
    int precision_bits) {
    if (energy.midpoint.rows() != energy.midpoint.cols() ||
        energy.midpoint.rows() != rhs.midpoint.rows() ||
        constraints.cols() != energy.midpoint.rows())
        throw std::invalid_argument("constrained Riesz dimensions are inconsistent");
    const int unknowns = static_cast<int>(energy.midpoint.rows());
    const int constraint_count = static_cast<int>(constraints.rows());
    const int columns = static_cast<int>(rhs.midpoint.cols());
    RieszBlockResult result;
    result.diagnostics.solve_count = 1;

    Enclosure saddle;
    saddle.midpoint = ComplexMatrix::Zero(
        unknowns + constraint_count, unknowns + constraint_count);
    saddle.radius = Eigen::MatrixXd::Zero(
        unknowns + constraint_count, unknowns + constraint_count);
    saddle.midpoint.topLeftCorner(unknowns, unknowns) = energy.midpoint;
    saddle.radius.topLeftCorner(unknowns, unknowns) = energy.radius;
    if (constraint_count > 0) {
        saddle.midpoint.topRightCorner(unknowns, constraint_count) =
            constraints.transpose().cast<Complex>();
        saddle.midpoint.bottomLeftCorner(constraint_count, unknowns) =
            constraints.cast<Complex>();
    }
    saddle.entries_verified = energy.entries_verified;

    Enclosure saddle_rhs;
    saddle_rhs.midpoint = ComplexMatrix::Zero(
        unknowns + constraint_count, columns);
    saddle_rhs.radius = Eigen::MatrixXd::Zero(
        unknowns + constraint_count, columns);
    saddle_rhs.midpoint.topRows(unknowns) = rhs.midpoint;
    saddle_rhs.radius.topRows(unknowns) = rhs.radius;
    saddle_rhs.entries_verified = rhs.entries_verified;

    Eigen::FullPivLU<ComplexMatrix> factorization(saddle.midpoint);
    if (!factorization.isInvertible())
        throw std::runtime_error("constrained Riesz midpoint saddle is singular");
    const ComplexMatrix approximate_solution = factorization.solve(saddle_rhs.midpoint);
    if (!approximate_solution.allFinite())
        throw std::runtime_error("constrained Riesz solve returned non-finite values");

    const Enclosure approximate = point_enclosure(approximate_solution);
    const Enclosure residual = enclosure_add(
        saddle_rhs, enclosure_multiply(saddle, approximate), -1.0);
    const solver::VerifiedScalarResult singular =
        solver::verified_minimum_singular_value(saddle, precision_bits);
    const double denominator = singular.metadata.verified
        ? singular.enclosure.lower
        : singular.approximation;

    Enclosure enclosed_solution = approximate;
    enclosed_solution.entries_verified = singular.metadata.verified
        && saddle_rhs.entries_verified;
    for (int column = 0; column < columns; ++column) {
        const double residual_bound = vector_enclosure_norm_upper(residual, column);
        const double error_bound = denominator > 0.0
            ? outward_nonnegative(residual_bound / denominator)
            : kInfinity;
        result.diagnostics.max_solution_error_bound = std::max(
            result.diagnostics.max_solution_error_bound, error_bound);
        for (int row = 0; row < unknowns + constraint_count; ++row)
            enclosed_solution.radius(row, column) = error_bound;
    }
    if (enclosed_solution.entries_verified)
        result.diagnostics.verified_solve_count = 1;

    result.solution.midpoint = enclosed_solution.midpoint.topRows(unknowns);
    result.solution.radius = enclosed_solution.radius.topRows(unknowns);
    result.solution.entries_verified = enclosed_solution.entries_verified;

    const Enclosure energy_times_solution =
        enclosure_multiply(energy, result.solution);
    result.gram = hermitianize(enclosure_multiply(
        enclosure_adjoint(result.solution), energy_times_solution));

    const ComplexMatrix local_solution = approximate_solution.topRows(unknowns);
    const ComplexMatrix multipliers = approximate_solution.bottomRows(constraint_count);
    const ComplexMatrix constraint_residual = constraints.cast<Complex>() * local_solution;
    ComplexMatrix stationarity = energy.midpoint * local_solution - rhs.midpoint;
    if (constraint_count > 0)
        stationarity += constraints.transpose().cast<Complex>() * multipliers;
    result.diagnostics.max_constraint_relative_residual =
        constraint_residual.norm() /
        std::max(1.0, constraints.norm() * local_solution.norm());
    result.diagnostics.max_stationarity_relative_residual =
        stationarity.norm() / std::max(1.0, rhs.midpoint.norm());
    const ComplexMatrix energy_gram =
        local_solution.adjoint() * energy.midpoint * local_solution;
    const ComplexMatrix action = local_solution.adjoint() * rhs.midpoint;
    result.diagnostics.max_energy_identity_relative_error =
        (energy_gram - action).norm() /
        std::max({1.0, energy_gram.norm(), action.norm()});
    return result;
}

void accumulate_riesz_diagnostics(const RieszCertificateDiagnostics &source,
                                  RieszCertificateDiagnostics &target) {
    target.solve_count += source.solve_count;
    target.verified_solve_count += source.verified_solve_count;
    target.max_constraint_relative_residual = std::max(
        target.max_constraint_relative_residual,
        source.max_constraint_relative_residual);
    target.max_stationarity_relative_residual = std::max(
        target.max_stationarity_relative_residual,
        source.max_stationarity_relative_residual);
    target.max_energy_identity_relative_error = std::max(
        target.max_energy_identity_relative_error,
        source.max_energy_identity_relative_error);
    target.max_solution_error_bound = std::max(
        target.max_solution_error_bound,
        source.max_solution_error_bound);
}

CertificateMatrixDiagnostics matrix_diagnostics(Enclosure enclosure) {
    CertificateMatrixDiagnostics result;
    result.enclosure = hermitianize(enclosure);
    const ComplexMatrix &midpoint = result.enclosure.midpoint;
    result.hermitian_relative_error =
        (midpoint - midpoint.adjoint()).norm() / std::max(1.0, midpoint.norm());
    Eigen::SelfAdjointEigenSolver<ComplexMatrix> eigensolver(midpoint);
    result.minimum_midpoint_eigenvalue = eigensolver.info() == Eigen::Success
        ? eigensolver.eigenvalues().minCoeff()
        : std::numeric_limits<double>::quiet_NaN();
    return result;
}

std::vector<int> free_nodes(const TriMesh &mesh) {
    std::vector<char> constrained(mesh.nodes.size(), false);
    for (int node : dirichlet_nodes(mesh)) constrained[node] = true;
    std::vector<int> result;
    for (int node = 0; node < static_cast<int>(mesh.nodes.size()); ++node) {
        if (!constrained[node]) result.push_back(node);
    }
    return result;
}

Enclosure zero_square_enclosure(int size) {
    return point_enclosure(ComplexMatrix::Zero(size, size));
}

Enclosure expand_local_gram(
    const Enclosure &local,
    const Triangle &triangle,
    const std::vector<int> &free_index,
    int free_count) {
    Enclosure result = zero_square_enclosure(free_count);
    result.entries_verified = local.entries_verified;
    for (int local_row = 0; local_row < 3; ++local_row) {
        const int row = free_index[triangle[local_row]];
        if (row < 0) continue;
        for (int local_column = 0; local_column < 3; ++local_column) {
            const int column = free_index[triangle[local_column]];
            if (column < 0) continue;
            result.midpoint(row, column) = local.midpoint(local_row, local_column);
            result.radius(row, column) = local.radius(local_row, local_column);
        }
    }
    return result;
}

GeneralizedSpectrumCertificate generalized_spectrum(
    const Enclosure &numerator,
    const Enclosure &denominator,
    const CorrectorCertificateConfig &config) {
    GeneralizedSpectrumCertificate result;
    result.verified_lambda = solver::verified_generalized_largest_eigenvalue(
        numerator, denominator, config.precision_bits);
    result.lambda_enclosure = result.verified_lambda.enclosure;
    result.lambda_max_approximation = result.verified_lambda.approximation;
    result.theta_approximation = std::sqrt(std::max(
        0.0, result.lambda_max_approximation));
    if (result.verified_lambda.metadata.verified) {
        result.theta_enclosure.lower = std::nextafter(
            std::sqrt(std::max(0.0, result.lambda_enclosure.lower)), 0.0);
        result.theta_enclosure.upper = std::nextafter(
            std::sqrt(std::max(0.0, result.lambda_enclosure.upper)), kInfinity);
    } else {
        result.theta_enclosure = {
            result.theta_approximation, result.theta_approximation};
    }

    const ComplexMatrix hermitian_denominator =
        0.5 * (denominator.midpoint + denominator.midpoint.adjoint());
    Eigen::LLT<ComplexMatrix> factor(hermitian_denominator);
    if (factor.info() != Eigen::Success)
        throw std::runtime_error("coarse energy midpoint is not positive definite");
    const ComplexMatrix inverse_lower = factor.matrixL().solve(
        ComplexMatrix::Identity(denominator.midpoint.rows(),
                                denominator.midpoint.cols()));
    const ComplexMatrix transformed = inverse_lower
        * (0.5 * (numerator.midpoint + numerator.midpoint.adjoint()))
        * inverse_lower.adjoint();
    Eigen::SelfAdjointEigenSolver<ComplexMatrix> eigensolver(transformed);
    if (eigensolver.info() != Eigen::Success)
        throw std::runtime_error("coarse generalized eigensolve failed");
    result.eigenvalues = eigensolver.eigenvalues();
    const double top = result.eigenvalues.maxCoeff();
    const double threshold = std::max(
        config.cluster_absolute_gap,
        config.cluster_relative_gap * std::max(1.0, std::abs(top)));
    int first = result.eigenvalues.size() - 1;
    while (first > 0 && top - result.eigenvalues(first - 1) <= threshold)
        --first;
    result.dominant_cluster_size = result.eigenvalues.size() - first;
    result.dominant_cluster_basis = inverse_lower.adjoint()
        * eigensolver.eigenvectors().rightCols(result.dominant_cluster_size);
    for (Eigen::Index column = 0;
         column < result.dominant_cluster_basis.cols(); ++column) {
        ComplexVector vector = result.dominant_cluster_basis.col(column);
        const double norm_squared = std::real(
            vector.dot(denominator.midpoint * vector));
        if (!(norm_squared > 0.0))
            throw std::runtime_error("dominant generalized eigenvector has zero energy");
        result.dominant_cluster_basis.col(column) /= std::sqrt(norm_squared);
    }
    const ComplexVector dominant = result.dominant_cluster_basis.rightCols(1);
    const double dominant_lambda = result.eigenvalues.maxCoeff();
    result.dominant_residual =
        (numerator.midpoint * dominant
         - dominant_lambda * denominator.midpoint * dominant).norm()
        / std::max({1.0,
                    (numerator.midpoint * dominant).norm(),
                    std::abs(dominant_lambda)
                        * (denominator.midpoint * dominant).norm()});
    return result;
}

GeneralizedSpectrumCertificate worse_spectrum(
    const GeneralizedSpectrumCertificate &primal,
    const GeneralizedSpectrumCertificate &adjoint) {
    const bool use_primal = primal.lambda_max_approximation
        >= adjoint.lambda_max_approximation;
    GeneralizedSpectrumCertificate result = use_primal ? primal : adjoint;
    result.lambda_max_approximation = std::max(
        primal.lambda_max_approximation, adjoint.lambda_max_approximation);
    result.theta_approximation = std::sqrt(std::max(
        0.0, result.lambda_max_approximation));
    const bool verified = primal.verified_lambda.metadata.verified
        && adjoint.verified_lambda.metadata.verified;
    if (verified) {
        result.lambda_enclosure.lower = std::max(
            primal.lambda_enclosure.lower, adjoint.lambda_enclosure.lower);
        result.lambda_enclosure.upper = std::max(
            primal.lambda_enclosure.upper, adjoint.lambda_enclosure.upper);
        result.theta_enclosure.lower = std::nextafter(
            std::sqrt(std::max(0.0, result.lambda_enclosure.lower)), 0.0);
        result.theta_enclosure.upper = std::nextafter(
            std::sqrt(std::max(0.0, result.lambda_enclosure.upper)), kInfinity);
    } else {
        result.lambda_enclosure = {
            result.lambda_max_approximation, result.lambda_max_approximation};
        result.theta_enclosure = {
            result.theta_approximation, result.theta_approximation};
    }
    result.verified_lambda.approximation = result.lambda_max_approximation;
    result.verified_lambda.enclosure = result.lambda_enclosure;
    result.verified_lambda.metadata.verified = verified;
    result.verified_lambda.metadata.backend_available =
        primal.verified_lambda.metadata.backend_available
        && adjoint.verified_lambda.metadata.backend_available;
    result.verified_lambda.metadata.backend =
        primal.verified_lambda.metadata.backend;
    result.verified_lambda.metadata.precision_bits =
        primal.verified_lambda.metadata.precision_bits;
    result.verified_lambda.metadata.rounding_mode =
        primal.verified_lambda.metadata.rounding_mode;
    if (!verified) {
        result.verified_lambda.metadata.failure_reason =
            "both primal and adjoint spectral enclosures are required";
    }
    return result;
}

double conjugation_relative_error(const Enclosure &primal,
                                  const Enclosure &adjoint) {
    return (adjoint.midpoint - primal.midpoint.conjugate()).norm()
        / std::max({1.0, primal.midpoint.norm(), adjoint.midpoint.norm()});
}

Enclosure assemble_total_gram(
    const Enclosure &energy,
    const Enclosure &rhs,
    const std::vector<AuditKernelPatch> &patches,
    int precision_bits,
    RieszCertificateDiagnostics &diagnostics) {
    const int columns = static_cast<int>(rhs.midpoint.cols());
    Enclosure gram = zero_square_enclosure(columns);
    for (const AuditKernelPatch &patch : patches) {
        if (patch.audit_dofs.empty()) continue;
        const Enclosure local_energy = restrict_square(energy, patch.audit_dofs);
        const Enclosure local_rhs = restrict_rows(rhs, patch.audit_dofs);
        const RieszBlockResult solved = solve_constrained_riesz(
            local_energy, patch.constraints, local_rhs, precision_bits);
        gram = enclosure_add(gram, solved.gram);
        accumulate_riesz_diagnostics(solved.diagnostics, diagnostics);
    }
    return hermitianize(gram);
}

FineAssemblyResult assemble_fine_grams(
    const AdaptiveMeshHierarchy &hierarchy,
    const HelmholtzLodModel &model,
    const HelmholtzOperators &audit_operators,
    bool algebraic_verified,
    const CertificateAssemblyEvidence &evidence,
    const std::vector<int> &free_coarse,
    int precision_bits) {
    const TriMesh &coarse = hierarchy.coarse_mesh();
    const int free_count = static_cast<int>(free_coarse.size());
    std::vector<int> free_index(coarse.nodes.size(), -1);
    for (int index = 0; index < free_count; ++index)
        free_index[free_coarse[index]] = index;

    const std::vector<TriMesh> no_meshes;
    const std::vector<Eigen::SparseMatrix<double>> no_prolongations;
    HelmholtzPatchAssembler assembler(
        coarse,
        hierarchy.cert_audit_mesh(),
        hierarchy.coarse_elements_to_cert_audit(),
        hierarchy.coarse_dg_to_cert_audit(),
        hierarchy.cert_audit_quasi_interpolation(),
        model.problem().patches,
        no_meshes,
        no_prolongations,
        no_prolongations,
        audit_operators);

    const Enclosure fine_to_audit = point_enclosure(
        ComplexMatrix(hierarchy.fine_to_cert_audit().cast<Complex>()));
    FineAssemblyResult result;
    result.primal = zero_square_enclosure(free_count);
    result.adjoint = zero_square_enclosure(free_count);
    result.element_primal.reserve(coarse.elems.size());
    result.element_adjoint.reserve(coarse.elems.size());

    for (int target = 0; target < static_cast<int>(coarse.elems.size()); ++target) {
        const HelmholtzPatchSystem system = assembler.assemble(target);
        ComplexMatrix local_energy = ComplexMatrix(system.stiffness);
        local_energy += audit_operators.wavenumber * audit_operators.wavenumber
            * ComplexMatrix(system.mass);
        const Enclosure energy = uniform_enclosure(
            local_energy, evidence.energy_entry_radius, algebraic_verified);
        const Enclosure local_operator = uniform_enclosure(
            ComplexMatrix(system.helmholtz), evidence.system_entry_radius,
            algebraic_verified);
        const Enclosure source = uniform_enclosure(
            system.rhs, evidence.local_source_entry_radius,
            algebraic_verified);

        ComplexMatrix corrector_fine = ComplexMatrix::Zero(
            model.problem().fine.nodes.size(), 3);
        for (const HelmholtzCorrectorEntry &entry :
             model.correctors().primal[target]) {
            corrector_fine(entry.row, entry.local_coarse_vertex) += entry.value;
        }
        const Enclosure corrector_audit = enclosure_multiply(
            fine_to_audit, point_enclosure(corrector_fine));
        const Enclosure corrector_local = restrict_rows(
            corrector_audit, system.local_vertices);

        const Enclosure primal_rhs = enclosure_add(
            source, enclosure_multiply(local_operator, corrector_local), -1.0);
        const Enclosure adjoint_rhs = enclosure_add(
            enclosure_conjugate(source),
            enclosure_multiply(
                enclosure_adjoint(local_operator),
                enclosure_conjugate(corrector_local)),
            -1.0);
        const RieszBlockResult primal = solve_constrained_riesz(
            energy, system.constraints, primal_rhs, precision_bits);
        const RieszBlockResult adjoint = solve_constrained_riesz(
            energy, system.constraints, adjoint_rhs, precision_bits);
        accumulate_riesz_diagnostics(
            primal.diagnostics, result.primal_diagnostics);
        accumulate_riesz_diagnostics(
            adjoint.diagnostics, result.adjoint_diagnostics);

        const Enclosure expanded_primal = expand_local_gram(
            primal.gram, coarse.elems[target], free_index, free_count);
        const Enclosure expanded_adjoint = expand_local_gram(
            adjoint.gram, coarse.elems[target], free_index, free_count);
        result.primal = enclosure_add(result.primal, expanded_primal);
        result.adjoint = enclosure_add(result.adjoint, expanded_adjoint);
        result.element_primal.push_back(expanded_primal);
        result.element_adjoint.push_back(expanded_adjoint);
    }
    result.primal = hermitianize(result.primal);
    result.adjoint = hermitianize(result.adjoint);
    return result;
}

double constant_value(const CertificateConstantRegistry &registry,
                      const std::string &name) {
    const CertificateConstant *constant = registry.find(name);
    return constant ? constant->value
                    : std::numeric_limits<double>::quiet_NaN();
}

void add_reason(std::vector<std::string> &reasons, std::string reason) {
    if (std::find(reasons.begin(), reasons.end(), reason) == reasons.end())
        reasons.push_back(std::move(reason));
}

std::string join_reasons(const std::vector<std::string> &reasons) {
    std::ostringstream message;
    for (std::size_t index = 0; index < reasons.size(); ++index) {
        if (index > 0) message << "; ";
        message << reasons[index];
    }
    return message.str();
}

std::vector<double> cluster_contributions(
    const GeneralizedSpectrumCertificate &spectrum,
    const std::vector<Enclosure> &element_grams,
    double &allocation_relative_error) {
    std::vector<double> result(element_grams.size(), 0.0);
    if (spectrum.dominant_cluster_basis.cols() == 0) {
        allocation_relative_error = 0.0;
        return result;
    }
    const double inverse_cluster =
        1.0 / static_cast<double>(spectrum.dominant_cluster_basis.cols());
    double sum = 0.0;
    for (std::size_t element = 0; element < element_grams.size(); ++element) {
        const ComplexMatrix projected = spectrum.dominant_cluster_basis.adjoint()
            * element_grams[element].midpoint
            * spectrum.dominant_cluster_basis;
        result[element] = std::max(0.0, projected.trace().real() * inverse_cluster);
        sum += result[element];
    }
    const double target = std::max(0.0, spectrum.lambda_max_approximation);
    if (sum > 0.0) {
        const double scale = target / sum;
        for (double &value : result) value *= scale;
    }
    double allocated = 0.0;
    for (double value : result) allocated += value;
    allocation_relative_error = std::abs(allocated - target)
        / std::max(1.0, target);
    return result;
}

} // namespace

void CertificateConstantRegistry::set(CertificateConstant constant) {
    if (constant.name.empty())
        throw std::invalid_argument("certificate constant name must not be empty");
    if (!std::isfinite(constant.value))
        throw std::invalid_argument("certificate constant value must be finite");
    entries_[constant.name] = std::move(constant);
}

const CertificateConstant *CertificateConstantRegistry::find(
    const std::string &name) const {
    const auto found = entries_.find(name);
    return found == entries_.end() ? nullptr : &found->second;
}

bool CertificateConstantRegistry::has_verified(
    const std::string &name,
    CertificateBoundDirection required_direction) const {
    const CertificateConstant *constant = find(name);
    return constant && constant->verified
        && direction_satisfies(constant->direction, required_direction)
        && !constant->source.empty()
        && !constant->mesh_class.empty();
}

std::vector<std::string>
CertificateConstantRegistry::missing_or_unverified_required() const {
    std::vector<std::string> result;
    for (const auto &[name, direction] : required_constants()) {
        if (!has_verified(name, direction)) result.push_back(name);
    }
    return result;
}

bool CertificateAssemblyEvidence::valid() const {
    return verified
        && std::isfinite(energy_entry_radius) && energy_entry_radius >= 0.0
        && std::isfinite(system_entry_radius) && system_entry_radius >= 0.0
        && std::isfinite(local_source_entry_radius)
        && local_source_entry_radius >= 0.0
        && !source.empty() && !hash.empty();
}

void derive_certificate_constants(
    CertificateConstantRegistry &registry,
    const AdaptiveMeshHierarchy &hierarchy,
    double wavenumber,
    int ell,
    int localization_shift,
    const KernelPatchPolicy &patch_policy) {
    if (!(wavenumber > 0.0) || ell < 0 || localization_shift < 0)
        throw std::invalid_argument("constant derivation parameters are invalid");
    const std::string mesh_class = "current conforming NVB hierarchy";
    const std::string hash = patch_policy.hash;
    registry.set({
        "C_ol(ell)", exact_patch_overlap(hierarchy.coarse_mesh(), ell),
        CertificateBoundDirection::Exact,
        "exact patch-incidence enumeration",
        "maximum number of ell-layer patches containing one coarse element",
        mesh_class, hash, true});
    registry.set({
        "C_ol(ell+s)",
        exact_patch_overlap(hierarchy.coarse_mesh(), ell + localization_shift),
        CertificateBoundDirection::Exact,
        "exact patch-incidence enumeration",
        "maximum number of (ell+s)-layer patches containing one coarse element",
        mesh_class, hash, true});

    if (!registry.find("s")) {
        registry.set({
            "s", static_cast<double>(localization_shift),
            CertificateBoundDirection::Exact,
            "caller-supplied localization shift",
            "the theorem-dependent shift has not been independently verified",
            mesh_class, hash, false});
    }

    const CertificateConstant *c_app = registry.find("C_app");
    const CertificateConstant *c_st = registry.find("C_st");
    const CertificateConstant *c_a = registry.find("C_a");
    const double kappa_h = wavenumber
        * maximum_coarse_diameter(hierarchy.coarse_mesh());
    if (c_app) {
        const double mu = c_app->value * kappa_h;
        if (mu < 1.0) {
            const double c_w = (1.0 - mu * mu) / (1.0 + mu * mu);
            registry.set({
                "c_W", c_w, CertificateBoundDirection::Lower,
                "paper equation (kernel-coercivity-new)",
                "(1-mu^2)/(1+mu^2), mu=C_app*kappa*max(H_T)",
                mesh_class, hash,
                c_app->verified
                    && direction_satisfies(
                        c_app->direction, CertificateBoundDirection::Upper)});
        }
    }
    if (c_app && c_st) {
        const double mu = c_app->value * kappa_h;
        registry.set({
            "C_Pi", 1.0 + c_st->value + mu,
            CertificateBoundDirection::Upper,
            "paper equation (Pi-kappa-new)",
            "1+C_st+C_app*kappa*max(H_T)",
            mesh_class, hash,
            c_app->verified && c_st->verified
                && direction_satisfies(
                    c_app->direction, CertificateBoundDirection::Upper)
                && direction_satisfies(
                    c_st->direction, CertificateBoundDirection::Upper)});
    }
    const CertificateConstant *c_w = registry.find("c_W");
    if (c_a && c_w && c_w->value > 0.0) {
        registry.set({
            "C_Fort", 1.0 + c_a->value / c_w->value,
            CertificateBoundDirection::Upper,
            "paper equation (CF-new)",
            "1+C_a/c_W",
            mesh_class, hash,
            c_a->verified && c_w->verified
                && direction_satisfies(
                    c_a->direction, CertificateBoundDirection::Upper)
                && direction_satisfies(
                    c_w->direction, CertificateBoundDirection::Lower)});
    }
}

CorrectorCertificateResult build_corrector_certificates(
    const AdaptiveMeshHierarchy &hierarchy,
    const HelmholtzLodModel &model,
    const HelmholtzOperators &audit_operators,
    CertificateConstantRegistry constants,
    double eta_H,
    bool eta_H_verified,
    const CertificateAssemblyEvidence &assembly_evidence,
    const CorrectorCertificateConfig &config) {
    if (!(eta_H >= 0.0)
        || !std::isfinite(eta_H))
        throw std::invalid_argument("coarse kernel residual estimator is invalid");
    if (config.precision_bits < 64 || config.cluster_relative_gap < 0.0
        || config.cluster_absolute_gap < 0.0
        || config.conjugation_tolerance < 0.0 || config.q0 < 0.0)
        throw std::invalid_argument("corrector certificate configuration is invalid");
    if (model.problem().coarse.nodes.size() != hierarchy.coarse_mesh().nodes.size()
        || model.problem().fine.nodes.size() != hierarchy.fine_mesh().nodes.size()
        || audit_operators.system.rows()
            != static_cast<int>(hierarchy.cert_audit_mesh().nodes.size()))
        throw std::invalid_argument("model, hierarchy, and audit operator dimensions differ");

    CorrectorCertificateResult result;
    result.assembly_evidence = assembly_evidence;
    result.eta_H = eta_H;
    result.eta_H_verified = eta_H_verified;
    result.patch_policy = audit_kernel_patch_policy(hierarchy);
    int localization_shift = 0;
    if (const CertificateConstant *shift = constants.find("s"))
        localization_shift = std::max(0, static_cast<int>(std::ceil(shift->value)));
    derive_certificate_constants(
        constants, hierarchy, model.config().wavenumber,
        model.config().ell, localization_shift, result.patch_policy);
    result.constants = constants;

    const bool algebraic_verified = assembly_evidence.valid();
    if (!algebraic_verified) {
        add_reason(result.conditional_reasons,
                   "finite-element assembly enclosure is not verified");
    }
    if (!solver::verified_spectrum_backend_available()) {
        add_reason(result.conditional_reasons,
                   "MPFR/MPFI verified spectrum backend is unavailable");
    }

    const std::vector<int> free_coarse = free_nodes(hierarchy.coarse_mesh());
    if (free_coarse.empty())
        throw std::runtime_error("coarse certificate space has no free nodes");
    const Enclosure coarse_to_audit_all = point_enclosure(
        ComplexMatrix(hierarchy.coarse_to_cert_audit().cast<Complex>()));
    const Enclosure coarse_to_audit = select_columns(
        coarse_to_audit_all, free_coarse);
    const Enclosure fine_to_audit = point_enclosure(
        ComplexMatrix(hierarchy.fine_to_cert_audit().cast<Complex>()));

    const ComplexSparseMatrix full_corrector_sparse =
        build_helmholtz_corrector_matrix(
            model.problem().coarse,
            static_cast<int>(model.problem().fine.nodes.size()),
            model.correctors().primal);
    const Enclosure full_corrector = point_enclosure(
        ComplexMatrix(full_corrector_sparse));
    const Enclosure selected_corrector = select_columns(
        full_corrector, free_coarse);
    const Enclosure audit_corrector = enclosure_multiply(
        fine_to_audit, selected_corrector);
    const Enclosure primal_remainder = enclosure_add(
        coarse_to_audit, audit_corrector, -1.0);
    const Enclosure adjoint_remainder = enclosure_add(
        coarse_to_audit, enclosure_conjugate(audit_corrector), -1.0);

    ComplexMatrix energy_midpoint = ComplexMatrix(audit_operators.stiffness);
    energy_midpoint += audit_operators.wavenumber * audit_operators.wavenumber
        * ComplexMatrix(audit_operators.mass);
    const Enclosure audit_energy = uniform_enclosure(
        energy_midpoint, assembly_evidence.energy_entry_radius,
        algebraic_verified);
    const Enclosure audit_system = uniform_enclosure(
        ComplexMatrix(audit_operators.system),
        assembly_evidence.system_entry_radius,
        algebraic_verified);

    const Enclosure coarse_energy = hermitianize(enclosure_multiply(
        enclosure_adjoint(coarse_to_audit),
        enclosure_multiply(audit_energy, coarse_to_audit)));
    result.matrices.coarse_energy = matrix_diagnostics(coarse_energy);

    const Enclosure total_primal_rhs = enclosure_multiply(
        audit_system, primal_remainder);
    const Enclosure total_adjoint_rhs = enclosure_multiply(
        enclosure_adjoint(audit_system), adjoint_remainder);
    const std::vector<AuditKernelPatch> audit_patches =
        build_audit_kernel_patches(hierarchy, result.patch_policy);
    const Enclosure total_primal = assemble_total_gram(
        audit_energy, total_primal_rhs, audit_patches,
        config.precision_bits, result.total_primal_riesz);
    const Enclosure total_adjoint = assemble_total_gram(
        audit_energy, total_adjoint_rhs, audit_patches,
        config.precision_bits, result.total_adjoint_riesz);
    result.matrices.total_primal = matrix_diagnostics(total_primal);
    result.matrices.total_adjoint = matrix_diagnostics(total_adjoint);

    FineAssemblyResult fine = assemble_fine_grams(
        hierarchy, model, audit_operators, algebraic_verified,
        assembly_evidence, free_coarse, config.precision_bits);
    result.matrices.fine_primal = matrix_diagnostics(fine.primal);
    result.matrices.fine_adjoint = matrix_diagnostics(fine.adjoint);
    result.matrices.fine_element_primal = std::move(fine.element_primal);
    result.matrices.fine_element_adjoint = std::move(fine.element_adjoint);
    result.fine_primal_riesz = fine.primal_diagnostics;
    result.fine_adjoint_riesz = fine.adjoint_diagnostics;

    result.total_conjugation_relative_error = conjugation_relative_error(
        total_primal, total_adjoint);
    result.fine_conjugation_relative_error = conjugation_relative_error(
        fine.primal, fine.adjoint);
    result.conjugation_passed =
        result.total_conjugation_relative_error <= config.conjugation_tolerance
        && result.fine_conjugation_relative_error <= config.conjugation_tolerance;
    result.used_independent_worse_side = !result.conjugation_passed;

    result.total_primal_spectrum = generalized_spectrum(
        total_primal, coarse_energy, config);
    result.total_adjoint_spectrum = generalized_spectrum(
        total_adjoint, coarse_energy, config);
    result.fine_primal_spectrum = generalized_spectrum(
        fine.primal, coarse_energy, config);
    result.fine_adjoint_spectrum = generalized_spectrum(
        fine.adjoint, coarse_energy, config);
    result.total_spectrum = worse_spectrum(
        result.total_primal_spectrum, result.total_adjoint_spectrum);
    result.fine_spectrum = worse_spectrum(
        result.fine_primal_spectrum, result.fine_adjoint_spectrum);

    const std::vector<int> free_audit = free_nodes(hierarchy.cert_audit_mesh());
    const Enclosure free_audit_system = restrict_square(audit_system, free_audit);
    const Enclosure free_audit_energy = restrict_square(audit_energy, free_audit);
    result.audit_infsup = solver::verified_energy_scaled_minimum_singular_value(
        free_audit_system, free_audit_energy, config.precision_bits);
    result.gamma_audit_approximation = result.audit_infsup.approximation;
    result.gamma_audit_lower = result.audit_infsup.metadata.verified
        ? result.audit_infsup.enclosure.lower
        : result.audit_infsup.approximation;

    result.theta_total_lower = result.total_spectrum.theta_enclosure.lower;
    result.theta_total_upper = result.total_spectrum.theta_enclosure.upper;
    result.theta_h_upper = result.fine_spectrum.theta_enclosure.upper;

    const double c_app = constant_value(constants, "C_app");
    const double c_sd = constant_value(constants, "C_sd");
    const double c_ov = constant_value(constants, "C_ov");
    const double c_a = constant_value(constants, "C_a");
    const double c_fort = constant_value(constants, "C_Fort");
    const double c_w = constant_value(constants, "c_W");
    const double c_pi = constant_value(constants, "C_Pi");
    const double c_loc = constant_value(constants, "C_loc");
    const double c_ol_ell = constant_value(constants, "C_ol(ell)");
    const double c_ol_shift = constant_value(constants, "C_ol(ell+s)");
    const double beta = constant_value(constants, "beta");
    result.mu = c_app * model.config().wavenumber
        * maximum_coarse_diameter(hierarchy.coarse_mesh());

    const bool formula_inputs_valid =
        std::isfinite(c_sd) && c_sd > 0.0
        && std::isfinite(c_ov) && c_ov > 0.0
        && std::isfinite(c_a) && c_a > 0.0
        && std::isfinite(c_fort) && c_fort > 0.0
        && std::isfinite(c_w) && c_w > 0.0
        && std::isfinite(c_pi) && c_pi > 0.0
        && std::isfinite(c_ol_ell) && c_ol_ell > 0.0
        && std::isfinite(result.gamma_audit_lower)
        && result.gamma_audit_lower > 0.0;
    if (formula_inputs_valid) {
        result.delta_total_lower =
            result.theta_total_lower / (c_a * c_ov);
        result.delta_total_upper =
            (c_sd / c_w) * result.theta_total_upper;
        result.delta_h_upper =
            std::sqrt(c_ol_ell) / c_w * result.theta_h_upper;
        result.delta_ell_lower = std::max(
            0.0, result.delta_total_lower - result.delta_h_upper);
        result.delta_ell_upper =
            result.delta_total_upper + result.delta_h_upper;
        if (std::isfinite(c_loc) && c_loc > 0.0
            && std::isfinite(c_ol_shift) && c_ol_shift > 0.0
            && std::isfinite(beta) && beta > 0.0 && beta < 1.0) {
            const double decay = c_loc * std::sqrt(c_ol_shift)
                * std::pow(beta, model.config().ell);
            result.delta_ell_upper = std::min(result.delta_ell_upper, decay);
            result.localization_decay_bound_used = true;
        }
        const double perturbation_factor =
            c_pi * c_fort / result.gamma_audit_lower;
        result.q_total = perturbation_factor * result.delta_total_upper;
        result.q_h = perturbation_factor * result.delta_h_upper;
        result.q_ell = perturbation_factor * result.delta_ell_upper;
        result.q_acceptable = result.q_total <= config.q0;
        const double epsilon = c_pi * result.delta_total_upper;
        result.stability_margin = result.gamma_audit_lower / (2.0 * c_fort)
            - c_a * epsilon * (2.0 + epsilon);
        result.lod_error_lower = eta_H / (c_a * c_ov);
        result.lod_error_upper = c_sd * (
            1.0 / c_w + perturbation_factor * result.delta_total_upper)
            * eta_H;
    } else {
        result.delta_total_lower = result.delta_total_upper =
            result.delta_h_upper = result.delta_ell_lower =
            result.delta_ell_upper = result.q_total = result.q_h =
            result.q_ell = result.stability_margin = result.lod_error_lower =
            result.lod_error_upper =
                std::numeric_limits<double>::quiet_NaN();
        add_reason(result.conditional_reasons,
                   "certificate formulas have missing or nonpositive inputs");
    }

    const bool use_primal_fine = result.fine_primal_spectrum.lambda_max_approximation
        >= result.fine_adjoint_spectrum.lambda_max_approximation;
    result.fine_marker_side = use_primal_fine ? "primal" : "adjoint";
    result.eta_h_element_squared = cluster_contributions(
        use_primal_fine ? result.fine_primal_spectrum
                        : result.fine_adjoint_spectrum,
        use_primal_fine ? result.matrices.fine_element_primal
                        : result.matrices.fine_element_adjoint,
        result.eta_h_allocation_relative_error);

    for (const std::string &name : constants.missing_or_unverified_required())
        add_reason(result.conditional_reasons,
                   "missing or unverified constant: " + name);
    if (!eta_H_verified)
        add_reason(result.conditional_reasons,
                   "coarse kernel residual eta_H is not verified");
    if (!result.total_spectrum.verified_lambda.metadata.verified)
        add_reason(result.conditional_reasons,
                   "total-corrector spectral enclosure is not verified");
    if (!result.fine_spectrum.verified_lambda.metadata.verified)
        add_reason(result.conditional_reasons,
                   "fine-corrector spectral enclosure is not verified");
    if (!result.audit_infsup.metadata.verified)
        add_reason(result.conditional_reasons,
                   "audit inf-sup lower enclosure is not verified");
    if (!(result.mu < 1.0))
        add_reason(result.conditional_reasons,
                   "coarse admissibility mu < 1 is not satisfied");
    result.stability_verified = result.conditional_reasons.empty()
        && std::isfinite(result.stability_margin)
        && result.stability_margin >= 0.0;
    if (std::isfinite(result.stability_margin) && result.stability_margin < 0.0)
        add_reason(result.conditional_reasons,
                   "verified stability margin is negative");

    result.status = result.conditional_reasons.empty()
        ? CorrectorCertificateStatus::Certified
        : CorrectorCertificateStatus::Conditional;
    result.verification_metadata.backend_available =
        solver::verified_spectrum_backend_available();
    result.verification_metadata.verified =
        result.status == CorrectorCertificateStatus::Certified;
    result.verification_metadata.backend = result.audit_infsup.metadata.backend;
    result.verification_metadata.precision_bits = config.precision_bits;
    result.verification_metadata.rounding_mode =
        result.audit_infsup.metadata.rounding_mode;
    result.verification_metadata.failure_reason =
        join_reasons(result.conditional_reasons);
    return result;
}

const char *certificate_bound_direction_name(
    CertificateBoundDirection direction) {
    switch (direction) {
    case CertificateBoundDirection::Lower:
        return "lower";
    case CertificateBoundDirection::Upper:
        return "upper";
    case CertificateBoundDirection::Exact:
        return "exact";
    }
    return "unknown";
}

const char *corrector_certificate_status_name(
    CorrectorCertificateStatus status) {
    switch (status) {
    case CorrectorCertificateStatus::Certified:
        return "certified";
    case CorrectorCertificateStatus::Conditional:
        return "conditional";
    case CorrectorCertificateStatus::Invalid:
        return "invalid";
    }
    return "unknown";
}

} // namespace lod2d::helmholtz::adaptive

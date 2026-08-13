#include "helmholtz/adaptive/certificates.h"

#include "helmholtz/boundary.h"
#include "helmholtz/patch_system.h"
#include "lod/patches.h"

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>

#include <algorithm>
#include <bit>
#include <cmath>
#include <complex>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <type_traits>
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

class FingerprintBuilder {
  public:
    void add_u64(std::uint64_t value) {
        for (int shift = 0; shift < 64; shift += 8) {
            state_ ^= static_cast<unsigned char>((value >> shift) & 0xffU);
            state_ *= 1099511628211ULL;
        }
    }

    void add_i64(std::int64_t value) {
        add_u64(std::bit_cast<std::uint64_t>(value));
    }

    void add_bool(bool value) { add_u64(value ? 1U : 0U); }

    void add_double(double value) {
        add_u64(std::bit_cast<std::uint64_t>(value));
    }

    void add_complex(Complex value) {
        add_double(value.real());
        add_double(value.imag());
    }

    void add_string(std::string_view value) {
        add_u64(static_cast<std::uint64_t>(value.size()));
        for (unsigned char byte : value) {
            state_ ^= byte;
            state_ *= 1099511628211ULL;
        }
    }

    std::string finish(std::string_view kind) const {
        std::ostringstream encoded;
        encoded << kind << ":fnv1a64:" << std::hex << std::setw(16)
                << std::setfill('0') << state_;
        return encoded.str();
    }

  private:
    std::uint64_t state_ = 14695981039346656037ULL;
};

void add_mesh(FingerprintBuilder &builder, const TriMesh &mesh) {
    builder.add_u64(mesh.nodes.size());
    for (const Point2 &node : mesh.nodes) {
        builder.add_double(node.x());
        builder.add_double(node.y());
    }
    builder.add_u64(mesh.elems.size());
    for (const Triangle &triangle : mesh.elems) {
        for (int node : triangle) builder.add_i64(node);
    }
    builder.add_u64(mesh.dirichlet.size());
    for (int node : mesh.dirichlet) builder.add_i64(node);
    builder.add_u64(mesh.boundary_edges.size());
    for (const BoundaryEdge &edge : mesh.boundary_edges) {
        builder.add_i64(edge.nodes[0]);
        builder.add_i64(edge.nodes[1]);
        builder.add_u64(static_cast<std::uint64_t>(edge.tag));
    }
}

template <class Scalar>
void add_scalar(FingerprintBuilder &builder, const Scalar &value) {
    if constexpr (std::is_same_v<Scalar, Complex>) {
        builder.add_complex(value);
    } else {
        builder.add_double(static_cast<double>(value));
    }
}

template <class Scalar, int Options, class StorageIndex>
void add_sparse_matrix(
    FingerprintBuilder &builder,
    const Eigen::SparseMatrix<Scalar, Options, StorageIndex> &matrix) {
    builder.add_i64(matrix.rows());
    builder.add_i64(matrix.cols());
    builder.add_i64(matrix.nonZeros());
    builder.add_i64(matrix.outerSize());
    for (int outer = 0; outer < matrix.outerSize(); ++outer) {
        builder.add_i64(outer);
        for (typename Eigen::SparseMatrix<Scalar, Options, StorageIndex>::InnerIterator
                 entry(matrix, outer);
             entry; ++entry) {
            builder.add_i64(entry.row());
            builder.add_i64(entry.col());
            add_scalar(builder, entry.value());
        }
    }
}

void add_double_vector(FingerprintBuilder &builder,
                       const std::vector<double> &values) {
    builder.add_u64(values.size());
    for (double value : values) builder.add_double(value);
}

void add_int_vector(FingerprintBuilder &builder,
                    const std::vector<int> &values) {
    builder.add_u64(values.size());
    for (int value : values) builder.add_i64(value);
}

void add_u64_vector(FingerprintBuilder &builder,
                    const std::vector<std::uint64_t> &values) {
    builder.add_u64(values.size());
    for (std::uint64_t value : values) builder.add_u64(value);
}

void add_operators(FingerprintBuilder &builder,
                   const HelmholtzOperators &operators) {
    builder.add_double(operators.wavenumber);
    builder.add_double(operators.boundary_beta);
    add_double_vector(builder, operators.diffusion);
    add_double_vector(builder, operators.refractive_index);
    add_int_vector(builder, operators.dirichlet_nodes);
    builder.add_u64(operators.element_blocks.size());
    for (const Eigen::Matrix3cd &block : operators.element_blocks) {
        for (Eigen::Index column = 0; column < block.cols(); ++column) {
            for (Eigen::Index row = 0; row < block.rows(); ++row)
                builder.add_complex(block(row, column));
        }
    }
    add_sparse_matrix(builder, operators.stiffness);
    add_sparse_matrix(builder, operators.mass);
    add_sparse_matrix(builder, operators.boundary_mass);
    add_sparse_matrix(builder, operators.system);
}

bool constant_matches_context(
    const CertificateConstant &constant,
    const CertificateContextFingerprint &context) {
    return context.complete()
        && constant.mesh_fingerprint == context.mesh
        && constant.pde_fingerprint == context.pde
        && constant.patch_policy_hash == context.patch_policy
        && constant.operator_fingerprint == context.operators;
}

void bind_constant_to_context(
    CertificateConstant &constant,
    const CertificateContextFingerprint &context) {
    constant.mesh_fingerprint = context.mesh;
    constant.pde_fingerprint = context.pde;
    constant.patch_policy_hash = context.patch_policy;
    constant.operator_fingerprint = context.operators;
}

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
        {"C_ol(ell)", CertificateBoundDirection::Upper}};
    return values;
}

bool required_constant_name(const std::string &name) {
    return std::any_of(
        required_constants().begin(), required_constants().end(),
        [&](const auto &required) { return required.first == name; });
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

// A point enclosure is sound only when the represented floating-point value is
// itself the mathematical input. This is used for exact zero accumulators and
// for a freely chosen residual-correction centre; it must never be used to
// attest a numerically assembled operator or a computed corrector.
Enclosure stored_value_point_enclosure(const ComplexMatrix &matrix) {
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
    const Enclosure &constraints,
    const Enclosure &rhs,
    int precision_bits) {
    validate_enclosure_dimensions(constraints);
    if (energy.midpoint.rows() != energy.midpoint.cols() ||
        energy.midpoint.rows() != rhs.midpoint.rows() ||
        constraints.midpoint.cols() != energy.midpoint.rows())
        throw std::invalid_argument("constrained Riesz dimensions are inconsistent");
    const int unknowns = static_cast<int>(energy.midpoint.rows());
    const int constraint_count = static_cast<int>(constraints.midpoint.rows());
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
            constraints.midpoint.adjoint();
        saddle.radius.topRightCorner(unknowns, constraint_count) =
            constraints.radius.transpose();
        saddle.midpoint.bottomLeftCorner(constraint_count, unknowns) =
            constraints.midpoint;
        saddle.radius.bottomLeftCorner(constraint_count, unknowns) =
            constraints.radius;
    }
    saddle.entries_verified = energy.entries_verified
        && constraints.entries_verified;

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

    const Enclosure approximate = stored_value_point_enclosure(approximate_solution);
    const Enclosure residual = enclosure_add(
        saddle_rhs, enclosure_multiply(saddle, approximate), -1.0);
    const solver::VerifiedScalarResult singular =
        solver::verified_minimum_singular_value(saddle, precision_bits);
    const double denominator = singular.metadata.verified
        ? singular.enclosure.lower
        : singular.approximation;

    Enclosure enclosed_solution = approximate;
    bool solution_verified = singular.metadata.verified
        && saddle.entries_verified
        && saddle_rhs.entries_verified
        && residual.entries_verified
        && std::isfinite(denominator) && denominator > 0.0;
    for (int column = 0; column < columns; ++column) {
        const double residual_bound = vector_enclosure_norm_upper(residual, column);
        const double error_bound = denominator > 0.0
            ? outward_nonnegative(residual_bound / denominator)
            : kInfinity;
        solution_verified = solution_verified && std::isfinite(error_bound);
        result.diagnostics.max_solution_error_bound = std::max(
            result.diagnostics.max_solution_error_bound, error_bound);
        for (int row = 0; row < unknowns + constraint_count; ++row)
            enclosed_solution.radius(row, column) = error_bound;
    }
    enclosed_solution.entries_verified = solution_verified;
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
    const ComplexMatrix constraint_residual =
        constraints.midpoint * local_solution;
    ComplexMatrix stationarity = energy.midpoint * local_solution - rhs.midpoint;
    if (constraint_count > 0)
        stationarity += constraints.midpoint.adjoint() * multipliers;
    result.diagnostics.max_constraint_relative_residual =
        constraint_residual.norm() /
        std::max(1.0, constraints.midpoint.norm() * local_solution.norm());
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
    return stored_value_point_enclosure(ComplexMatrix::Zero(size, size));
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
    const CertificateAssemblyEvidence &evidence,
    bool algebraic_verified,
    int precision_bits,
    RieszCertificateDiagnostics &diagnostics) {
    const int columns = static_cast<int>(rhs.midpoint.cols());
    Enclosure gram = zero_square_enclosure(columns);
    for (const AuditKernelPatch &patch : patches) {
        if (patch.discrete_dofs.empty()) continue;
        const Enclosure local_energy = restrict_square(energy, patch.discrete_dofs);
        const Enclosure local_rhs = restrict_rows(rhs, patch.discrete_dofs);
        const Enclosure constraints = uniform_enclosure(
            patch.constraints.cast<Complex>(),
            evidence.constraint_entry_radius,
            algebraic_verified);
        const RieszBlockResult solved = solve_constrained_riesz(
            local_energy, constraints, local_rhs, precision_bits);
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

    const Enclosure fine_to_audit = uniform_enclosure(
        ComplexMatrix(hierarchy.fine_to_cert_audit().cast<Complex>()),
        evidence.prolongation_entry_radius,
        algebraic_verified);
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
        const Enclosure enclosed_corrector = uniform_enclosure(
            corrector_fine,
            evidence.corrector_entry_radius,
            algebraic_verified);
        const Enclosure corrector_audit = enclosure_multiply(
            fine_to_audit, enclosed_corrector);
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
        const Enclosure constraints = uniform_enclosure(
            system.constraints.cast<Complex>(),
            evidence.constraint_entry_radius,
            algebraic_verified);
        const RieszBlockResult primal = solve_constrained_riesz(
            energy, constraints, primal_rhs, precision_bits);
        const RieszBlockResult adjoint = solve_constrained_riesz(
            energy, constraints, adjoint_rhs, precision_bits);
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

bool identical_mesh(const TriMesh &left, const TriMesh &right) {
    if (left.nodes.size() != right.nodes.size()
        || left.elems != right.elems
        || left.dirichlet != right.dirichlet
        || left.boundary_edges.size() != right.boundary_edges.size()) {
        return false;
    }
    for (std::size_t node = 0; node < left.nodes.size(); ++node) {
        if (left.nodes[node].x() != right.nodes[node].x()
            || left.nodes[node].y() != right.nodes[node].y()) {
            return false;
        }
    }
    for (std::size_t edge = 0; edge < left.boundary_edges.size(); ++edge) {
        if (left.boundary_edges[edge].nodes
                != right.boundary_edges[edge].nodes
            || left.boundary_edges[edge].tag
                != right.boundary_edges[edge].tag) {
            return false;
        }
    }
    return true;
}

template <class Scalar>
bool identical_sparse(const Eigen::SparseMatrix<Scalar> &left,
                      const Eigen::SparseMatrix<Scalar> &right) {
    return left.rows() == right.rows() && left.cols() == right.cols()
        && (left - right).norm() == 0.0;
}

bool identical_operators(const HelmholtzOperators &left,
                         const HelmholtzOperators &right) {
    if (left.wavenumber != right.wavenumber
        || left.boundary_beta != right.boundary_beta
        || left.diffusion != right.diffusion
        || left.refractive_index != right.refractive_index
        || left.dirichlet_nodes != right.dirichlet_nodes
        || left.element_blocks.size() != right.element_blocks.size()
        || !identical_sparse(left.stiffness, right.stiffness)
        || !identical_sparse(left.mass, right.mass)
        || !identical_sparse(left.boundary_mass, right.boundary_mass)
        || !identical_sparse(left.system, right.system)) {
        return false;
    }
    for (std::size_t element = 0; element < left.element_blocks.size(); ++element) {
        if ((left.element_blocks[element] - right.element_blocks[element]).norm()
            != 0.0) {
            return false;
        }
    }
    return true;
}

void validate_certificate_inputs(
    const AdaptiveMeshHierarchy &hierarchy,
    const HelmholtzLodModel &model,
    const HelmholtzOperators &audit_operators) {
    const HelmholtzProblemData &problem = model.problem();
    if (!identical_mesh(problem.coarse, hierarchy.coarse_mesh())
        || !identical_mesh(problem.fine, hierarchy.fine_mesh())
        || problem.coarse_element_levels != hierarchy.coarse_levels()
        || problem.fine_element_levels != hierarchy.fine_element_levels()
        || !identical_sparse(problem.coarse_to_fine,
                             hierarchy.coarse_to_fine())
        || !identical_sparse(problem.fine_element_prolongation,
                             hierarchy.coarse_elements_to_fine())
        || !identical_sparse(problem.fine_dg_prolongation,
                             hierarchy.coarse_dg_to_fine())
        || !identical_sparse(problem.quasi_interpolation,
                             hierarchy.fine_quasi_interpolation())
        || !identical_sparse(
            problem.patches,
            build_patches(hierarchy.coarse_mesh(), model.config().ell))) {
        throw std::invalid_argument(
            "LOD model is not built from the supplied adaptive hierarchy");
    }

    const HelmholtzOperators rebuilt_model = assemble_helmholtz_operators(
        hierarchy.fine_mesh(), model.config().wavenumber,
        model.operators().diffusion, model.operators().refractive_index,
        model.config().boundary_beta);
    const HelmholtzOperators rebuilt_audit = assemble_helmholtz_operators(
        hierarchy.cert_audit_mesh(), audit_operators.wavenumber,
        audit_operators.diffusion, audit_operators.refractive_index,
        audit_operators.boundary_beta);
    if (!identical_operators(model.operators(), rebuilt_model)
        || !identical_operators(audit_operators, rebuilt_audit)
        || model.config().wavenumber != audit_operators.wavenumber
        || model.config().boundary_beta != audit_operators.boundary_beta) {
        throw std::invalid_argument(
            "LOD and audit operators do not represent the supplied PDE data");
    }

    const Eigen::VectorXd fine_diffusion = Eigen::Map<const Eigen::VectorXd>(
        model.operators().diffusion.data(), model.operators().diffusion.size());
    const Eigen::VectorXd fine_refractive = Eigen::Map<const Eigen::VectorXd>(
        model.operators().refractive_index.data(),
        model.operators().refractive_index.size());
    const Eigen::VectorXd audit_diffusion = Eigen::Map<const Eigen::VectorXd>(
        audit_operators.diffusion.data(), audit_operators.diffusion.size());
    const Eigen::VectorXd audit_refractive = Eigen::Map<const Eigen::VectorXd>(
        audit_operators.refractive_index.data(),
        audit_operators.refractive_index.size());
    if ((hierarchy.fine_elements_to_cert_audit() * fine_diffusion
             - audit_diffusion).norm() != 0.0
        || (hierarchy.fine_elements_to_cert_audit() * fine_refractive
                - audit_refractive).norm() != 0.0) {
        throw std::invalid_argument(
            "LOD and audit coefficient fields are inconsistent across the hierarchy");
    }
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

bool CertificateContextFingerprint::complete() const {
    return !mesh.empty() && !pde.empty() && !patch_policy.empty()
        && !operators.empty();
}

CertificateContextFingerprint certificate_context_fingerprint(
    const AdaptiveMeshHierarchy &hierarchy,
    const HelmholtzLodModel &model,
    const HelmholtzOperators &audit_operators,
    const KernelPatchPolicy &patch_policy) {
    CertificateContextFingerprint result;

    FingerprintBuilder mesh;
    mesh.add_string("certificate-mesh-context-v1");
    add_mesh(mesh, hierarchy.initial_mesh());
    add_mesh(mesh, hierarchy.coarse_mesh());
    add_mesh(mesh, hierarchy.fine_mesh());
    add_mesh(mesh, hierarchy.cert_audit_mesh());
    add_int_vector(mesh, hierarchy.coarse_levels());
    add_int_vector(mesh, hierarchy.fine_element_levels());
    add_int_vector(mesh, hierarchy.cert_audit_element_levels());
    add_u64_vector(mesh, hierarchy.coarse_element_ids());
    add_u64_vector(mesh, hierarchy.coarse_parent_ids());
    mesh.add_i64(hierarchy.fine_level());
    mesh.add_u64(hierarchy.coarse_mesh_version());
    mesh.add_u64(hierarchy.fine_mesh_version());
    mesh.add_u64(hierarchy.cert_audit_mesh_version());
    mesh.add_u64(hierarchy.interpolation_version());
    mesh.add_u64(hierarchy.boundary_version());
    mesh.add_u64(hierarchy.corrector_space_version());
    add_sparse_matrix(mesh, hierarchy.coarse_to_fine());
    add_sparse_matrix(mesh, hierarchy.fine_to_cert_audit());
    add_sparse_matrix(mesh, hierarchy.coarse_to_cert_audit());
    add_sparse_matrix(mesh, hierarchy.fine_quasi_interpolation());
    add_sparse_matrix(mesh, hierarchy.cert_audit_quasi_interpolation());
    result.mesh = mesh.finish("mesh");

    const HelmholtzProblemConfig &config = model.config();
    FingerprintBuilder pde;
    pde.add_string("certificate-pde-context-v1");
    pde.add_double(config.wavenumber);
    pde.add_double(config.boundary_beta);
    pde.add_u64(static_cast<std::uint64_t>(config.mode));
    add_double_vector(pde, config.diffusion);
    add_double_vector(pde, config.refractive_index);
    pde.add_i64(config.quadrature.base_triangle_order);
    pde.add_i64(config.quadrature.gaussian_triangle_order);
    pde.add_i64(config.quadrature.singular_triangle_order);
    pde.add_i64(config.quadrature.max_recursive_subdivisions);
    pde.add_u64(static_cast<std::uint64_t>(
        config.quadrature_context.integrand_class));
    pde.add_double(config.quadrature_context.feature_point.x());
    pde.add_double(config.quadrature_context.feature_point.y());
    pde.add_double(config.quadrature_context.feature_scale);
    pde.add_u64(static_cast<std::uint64_t>(config.patch_solver.kind));
    pde.add_i64(config.patch_solver.symbolic_cache_slots);
    pde.add_bool(config.patch_solver.reuse_identical_factorization);
    pde.add_i64(config.patch_solver.gmres.restart);
    pde.add_i64(config.patch_solver.gmres.max_iterations);
    pde.add_double(config.patch_solver.gmres.relative_tolerance);
    pde.add_double(config.patch_solver.gmres.absolute_tolerance);
    pde.add_bool(config.patch_solver.gmres.reorthogonalize);
    pde.add_u64(static_cast<std::uint64_t>(config.patch_solver.shifted.rule));
    pde.add_double(config.patch_solver.shifted.alpha);
    pde.add_double(config.patch_solver.shifted.absolute_epsilon);
    pde.add_u64(static_cast<std::uint64_t>(
        config.patch_solver.shifted.inverse));
    pde.add_i64(config.patch_solver.shifted.pre_smooth);
    pde.add_i64(config.patch_solver.shifted.post_smooth);
    pde.add_i64(config.patch_solver.shifted.coarse_max_dofs);
    pde.add_double(config.patch_solver.shifted.jacobi_weight);
    pde.add_bool(config.patch_solver.fallback_to_direct);
    result.pde = pde.finish("pde");

    result.patch_policy = patch_policy.hash;

    FingerprintBuilder operators;
    operators.add_string("certificate-operator-context-v1");
    add_operators(operators, model.operators());
    add_operators(operators, audit_operators);
    add_sparse_matrix(operators, model.problem().coarse_to_fine);
    add_sparse_matrix(operators, model.problem().fine_element_prolongation);
    add_sparse_matrix(operators, model.problem().fine_dg_prolongation);
    add_sparse_matrix(operators, model.problem().quasi_interpolation);
    add_sparse_matrix(operators, model.problem().patches);
    operators.add_u64(model.correctors().primal.size());
    for (const HelmholtzElementCorrector &element : model.correctors().primal) {
        operators.add_u64(element.size());
        for (const HelmholtzCorrectorEntry &entry : element) {
            operators.add_i64(entry.row);
            operators.add_i64(entry.local_coarse_vertex);
            operators.add_complex(entry.value);
        }
    }
    result.operators = operators.finish("operators");
    return result;
}

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

bool CertificateConstantRegistry::has_verified(
    const std::string &name,
    CertificateBoundDirection required_direction,
    const CertificateContextFingerprint &context) const {
    const CertificateConstant *constant = find(name);
    return constant && has_verified(name, required_direction)
        && constant_matches_context(*constant, context);
}

std::vector<std::string>
CertificateConstantRegistry::missing_or_unverified_required() const {
    std::vector<std::string> result;
    for (const auto &[name, direction] : required_constants()) {
        if (!has_verified(name, direction)) result.push_back(name);
    }
    return result;
}

std::vector<std::string>
CertificateConstantRegistry::missing_or_unverified_required(
    const CertificateContextFingerprint &context) const {
    std::vector<std::string> result;
    for (const auto &[name, direction] : required_constants()) {
        if (!has_verified(name, direction, context)) result.push_back(name);
    }
    return result;
}

bool CertificateAssemblyEvidence::valid() const {
    const auto positive_finite = [](double value) {
        return std::isfinite(value) && value > 0.0;
    };
    return verified
        && positive_finite(energy_entry_radius)
        && positive_finite(system_entry_radius)
        && positive_finite(local_source_entry_radius)
        && positive_finite(prolongation_entry_radius)
        && positive_finite(corrector_entry_radius)
        && positive_finite(constraint_entry_radius)
        && !source.empty() && !hash.empty()
        && !mesh_fingerprint.empty() && !pde_fingerprint.empty()
        && !patch_policy_hash.empty() && !operator_fingerprint.empty();
}

bool CertificateAssemblyEvidence::valid_for(
    const CertificateContextFingerprint &context) const {
    return valid() && context.complete()
        && mesh_fingerprint == context.mesh
        && pde_fingerprint == context.pde
        && patch_policy_hash == context.patch_policy
        && operator_fingerprint == context.operators;
}

void derive_certificate_constants(
    CertificateConstantRegistry &registry,
    const AdaptiveMeshHierarchy &hierarchy,
    double wavenumber,
    int ell,
    int localization_shift,
    const KernelPatchPolicy &patch_policy,
    const CertificateContextFingerprint &context) {
    if (!(wavenumber > 0.0) || ell < 0 || localization_shift < 0)
        throw std::invalid_argument("constant derivation parameters are invalid");
    const std::string mesh_class = "current conforming NVB hierarchy";
    const std::string hash = patch_policy.hash;
    const auto set_current = [&](CertificateConstant constant) {
        bind_constant_to_context(constant, context);
        registry.set(std::move(constant));
    };
    set_current({
        "C_ol(ell)", exact_patch_overlap(hierarchy.coarse_mesh(), ell),
        CertificateBoundDirection::Exact,
        "exact patch-incidence enumeration",
        "maximum number of ell-layer patches containing one coarse element",
        mesh_class, hash, true});
    set_current({
        "C_ol(ell+s)",
        exact_patch_overlap(hierarchy.coarse_mesh(), ell + localization_shift),
        CertificateBoundDirection::Exact,
        "exact patch-incidence enumeration",
        "maximum number of (ell+s)-layer patches containing one coarse element",
        mesh_class, hash, true});

    if (!registry.find("s")) {
        set_current({
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
            set_current({
                "c_W", c_w, CertificateBoundDirection::Lower,
                "paper equation (kernel-coercivity-new)",
                "(1-mu^2)/(1+mu^2), mu=C_app*kappa*max(H_T)",
                mesh_class, hash, false});
        }
    }
    if (c_app && c_st) {
        const double mu = c_app->value * kappa_h;
        set_current({
            "C_Pi", 1.0 + c_st->value + mu,
            CertificateBoundDirection::Upper,
            "paper equation (Pi-kappa-new)",
            "1+C_st+C_app*kappa*max(H_T)",
            mesh_class, hash, false});
    }
    const CertificateConstant *c_w = registry.find("c_W");
    if (c_a && c_w && c_w->value > 0.0) {
        set_current({
            "C_Fort", 1.0 + c_a->value / c_w->value,
            CertificateBoundDirection::Upper,
            "paper equation (CF-new)",
            "1+C_a/c_W",
            mesh_class, hash, false});
    }
}

CorrectorCertificateResult build_corrector_certificates(
    const AdaptiveMeshHierarchy &hierarchy,
    const HelmholtzLodModel &model,
    const HelmholtzOperators &audit_operators,
    CertificateConstantRegistry constants,
    double eta_H,
    const AuditKernelResidualEvidence &eta_H_evidence,
    const CertificateAssemblyEvidence &assembly_evidence,
    const CorrectorCertificateConfig &config) {
    if (!(eta_H >= 0.0)
        || !std::isfinite(eta_H))
        throw std::invalid_argument("coarse kernel residual estimator is invalid");
    if (config.precision_bits < 64 || config.cluster_relative_gap < 0.0
        || config.cluster_absolute_gap < 0.0
        || config.conjugation_tolerance < 0.0 || config.q0 < 0.0)
        throw std::invalid_argument("corrector certificate configuration is invalid");
    validate_certificate_inputs(hierarchy, model, audit_operators);

    CorrectorCertificateResult result;
    result.assembly_evidence = assembly_evidence;
    result.eta_H = eta_H;
    result.eta_H_evidence = eta_H_evidence;
    result.patch_policy = audit_kernel_patch_policy(hierarchy);
    result.context_fingerprint = certificate_context_fingerprint(
        hierarchy, model, audit_operators, result.patch_policy);
    const std::string expected_eta_H_context =
        audit_kernel_residual_context_fingerprint(
            hierarchy, audit_operators, result.patch_policy);
    const bool eta_H_evidence_matches_context =
        !eta_H_evidence.context_fingerprint().empty()
        && !eta_H_evidence.diagnostic_fingerprint().empty()
        && eta_H_evidence.context_fingerprint() == expected_eta_H_context
        && eta_H_evidence.matches_eta(eta_H);

    CertificateConstantRegistry context_scoped_constants;
    for (const auto &[name, supplied] : constants.entries()) {
        CertificateConstant constant = supplied;
        if (constant.verified
            && !constant_matches_context(
                constant, result.context_fingerprint)) {
            constant.verified = false;
            if (required_constant_name(name)) {
                add_reason(result.conditional_reasons,
                           "constant context fingerprint mismatch: " + name);
            }
        }
        context_scoped_constants.set(std::move(constant));
    }
    constants = std::move(context_scoped_constants);

    int localization_shift = 0;
    if (constants.has_verified(
            "s", CertificateBoundDirection::Upper,
            result.context_fingerprint)) {
        const CertificateConstant *shift = constants.find("s");
        localization_shift = std::max(0, static_cast<int>(std::ceil(shift->value)));
    }
    derive_certificate_constants(
        constants, hierarchy, model.config().wavenumber,
        model.config().ell, localization_shift, result.patch_policy,
        result.context_fingerprint);
    result.constants = constants;

    const bool assembly_evidence_matches_context = assembly_evidence.valid_for(
        result.context_fingerprint);
    // Input radii alone do not make the subsequent ordinary-double matrix
    // additions and products rigorous.  The current code records useful
    // diagnostic radii, but must not mark them verified until every radius
    // accumulation is performed by a directed interval backend.
    result.matrix_enclosure_arithmetic_verified = false;
    const bool algebraic_verified = false;
    if (!assembly_evidence.valid()) {
        add_reason(result.conditional_reasons,
                   "finite-element assembly enclosure is incomplete or unverified");
    } else if (!assembly_evidence_matches_context) {
        add_reason(result.conditional_reasons,
                   "finite-element assembly evidence context fingerprint mismatch");
    } else {
        add_reason(
            result.conditional_reasons,
            "matrix enclosure propagation lacks directed interval arithmetic");
    }
    if (!solver::verified_spectrum_backend_available()) {
        add_reason(result.conditional_reasons,
                   "MPFR/MPFI verified spectrum backend is unavailable");
    }

    const std::vector<int> free_coarse = free_nodes(hierarchy.coarse_mesh());
    if (free_coarse.empty())
        throw std::runtime_error("coarse certificate space has no free nodes");
    const Enclosure coarse_to_audit_all = uniform_enclosure(
        ComplexMatrix(hierarchy.coarse_to_cert_audit().cast<Complex>()),
        assembly_evidence.prolongation_entry_radius,
        algebraic_verified);
    const Enclosure coarse_to_audit = select_columns(
        coarse_to_audit_all, free_coarse);
    const Enclosure fine_to_audit = uniform_enclosure(
        ComplexMatrix(hierarchy.fine_to_cert_audit().cast<Complex>()),
        assembly_evidence.prolongation_entry_radius,
        algebraic_verified);

    const ComplexSparseMatrix full_corrector_sparse =
        build_helmholtz_corrector_matrix(
            model.problem().coarse,
            static_cast<int>(model.problem().fine.nodes.size()),
            model.correctors().primal);
    const Enclosure full_corrector = uniform_enclosure(
        ComplexMatrix(full_corrector_sparse),
        assembly_evidence.corrector_entry_radius,
        algebraic_verified);
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
        assembly_evidence, algebraic_verified,
        config.precision_bits, result.total_primal_riesz);
    const Enclosure total_adjoint = assemble_total_gram(
        audit_energy, total_adjoint_rhs, audit_patches,
        assembly_evidence, algebraic_verified,
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
    // The current model retains only a primal corrector. The adjoint objects
    // above are conjugation-derived diagnostics, not an independently solved
    // fallback, so a failed invariance gate must remain fail-closed.
    result.used_independent_worse_side = false;
    if (!result.conjugation_passed) {
        add_reason(
            result.conditional_reasons,
            "conjugation invariance failed and no independently enclosed adjoint corrector is available");
    }

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
        const bool localization_decay_verified =
            constants.has_verified(
                "C_loc", CertificateBoundDirection::Upper,
                result.context_fingerprint)
            && constants.has_verified(
                "C_ol(ell+s)", CertificateBoundDirection::Upper,
                result.context_fingerprint)
            && constants.has_verified(
                "beta", CertificateBoundDirection::Upper,
                result.context_fingerprint)
            && constants.has_verified(
                "s", CertificateBoundDirection::Upper,
                result.context_fingerprint);
        if (localization_decay_verified
            && std::isfinite(c_loc) && c_loc > 0.0
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

    for (const std::string &name :
         constants.missing_or_unverified_required(result.context_fingerprint))
        add_reason(result.conditional_reasons,
                   "missing or unverified constant: " + name);
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
    if (std::isfinite(result.stability_margin) && result.stability_margin < 0.0)
        add_reason(result.conditional_reasons,
                   "verified stability margin is negative");
    // The current implementation rigorously encloses matrix and spectral
    // stages, but combines theorem constants, delta/q values, stability, and
    // LOD error endpoints with ordinary double arithmetic.  Until that
    // scalar chain uses directed interval rounding, certification must stop
    // here even in an MPFR-enabled matrix build.
    result.scalar_formula_enclosures_verified = false;
    add_reason(
        result.conditional_reasons,
        "scalar certificate formulas lack directed-rounding enclosures");

    result.corrector_conditional_reasons = result.conditional_reasons;
    result.corrector_status = result.corrector_conditional_reasons.empty()
        ? CorrectorCertificateStatus::Certified
        : CorrectorCertificateStatus::Conditional;
    result.stability_verified =
        result.corrector_status == CorrectorCertificateStatus::Certified
        && std::isfinite(result.stability_margin)
        && result.stability_margin >= 0.0;

    if (!eta_H_evidence.verified()) {
        std::string reason = "coarse kernel residual eta_H is not verified";
        if (!eta_H_evidence.failure_reason().empty())
            reason += ": " + eta_H_evidence.failure_reason();
        add_reason(result.conditional_reasons, std::move(reason));
    } else if (!eta_H_evidence_matches_context) {
        add_reason(
            result.conditional_reasons,
            "coarse kernel residual eta_H evidence context fingerprint mismatch");
    }

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

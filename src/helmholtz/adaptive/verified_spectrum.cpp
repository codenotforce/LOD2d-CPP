#include "helmholtz/adaptive/verified_spectrum.h"

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <Eigen/SVD>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifdef LOD_ENABLE_VERIFIED_CERTIFICATES
#include <mpfi.h>
#include <mpfr.h>
#endif

namespace lod2d::solver {
namespace {

using helmholtz::ComplexMatrix;

constexpr double kInfinity = std::numeric_limits<double>::infinity();

bool valid_enclosure(const ComplexMatrixEnclosure &matrix, bool require_square,
                     std::string &reason) {
    if (matrix.midpoint.rows() == 0 || matrix.midpoint.cols() == 0) {
        reason = "matrix enclosure must be nonempty";
        return false;
    }
    if (require_square && matrix.midpoint.rows() != matrix.midpoint.cols()) {
        reason = "matrix enclosure must be square";
        return false;
    }
    if (matrix.radius.rows() != matrix.midpoint.rows() ||
        matrix.radius.cols() != matrix.midpoint.cols()) {
        reason = "matrix radius has incompatible dimensions";
        return false;
    }
    for (Eigen::Index i = 0; i < matrix.midpoint.rows(); ++i) {
        for (Eigen::Index j = 0; j < matrix.midpoint.cols(); ++j) {
            const auto value = matrix.midpoint(i, j);
            const double radius = matrix.radius(i, j);
            if (!std::isfinite(value.real()) || !std::isfinite(value.imag()) ||
                !std::isfinite(radius) || radius < 0.0) {
                reason = "matrix enclosure contains a non-finite or negative radius";
                return false;
            }
        }
    }
    return true;
}

bool encloses_hermitian_matrix(const ComplexMatrixEnclosure &matrix,
                               std::string &reason) {
    for (Eigen::Index i = 0; i < matrix.midpoint.rows(); ++i) {
        if (std::abs(matrix.midpoint(i, i).imag()) > matrix.radius(i, i)) {
            reason = "diagonal enclosure excludes a real Hermitian entry";
            return false;
        }
        for (Eigen::Index j = i + 1; j < matrix.midpoint.cols(); ++j) {
            const auto mismatch = matrix.midpoint(i, j) -
                                  std::conj(matrix.midpoint(j, i));
            const double allowance = matrix.radius(i, j) + matrix.radius(j, i);
            if (std::abs(mismatch.real()) > allowance ||
                std::abs(mismatch.imag()) > allowance) {
                reason = "entry enclosures exclude Hermitian conjugacy";
                return false;
            }
        }
    }
    return true;
}

#ifdef LOD_ENABLE_VERIFIED_CERTIFICATES

Eigen::MatrixXd complex_real_block(const ComplexMatrix &matrix) {
    Eigen::MatrixXd block(2 * matrix.rows(), 2 * matrix.cols());
    block.topLeftCorner(matrix.rows(), matrix.cols()) = matrix.real();
    block.topRightCorner(matrix.rows(), matrix.cols()) = -matrix.imag();
    block.bottomLeftCorner(matrix.rows(), matrix.cols()) = matrix.imag();
    block.bottomRightCorner(matrix.rows(), matrix.cols()) = matrix.real();
    return block;
}

Eigen::MatrixXd complex_real_block_radius(const Eigen::MatrixXd &radius) {
    Eigen::MatrixXd block(2 * radius.rows(), 2 * radius.cols());
    block.topLeftCorner(radius.rows(), radius.cols()) = radius;
    block.topRightCorner(radius.rows(), radius.cols()) = radius;
    block.bottomLeftCorner(radius.rows(), radius.cols()) = radius;
    block.bottomRightCorner(radius.rows(), radius.cols()) = radius;
    return block;
}

#endif

double generalized_largest_eigenvalue_approximation(
    const ComplexMatrix &numerator, const ComplexMatrix &denominator) {
    const ComplexMatrix hermitian_denominator =
        0.5 * (denominator + denominator.adjoint());
    Eigen::LLT<ComplexMatrix> factor(hermitian_denominator);
    if (factor.info() != Eigen::Success) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const ComplexMatrix identity = ComplexMatrix::Identity(
        denominator.rows(), denominator.cols());
    const ComplexMatrix inverse_lower = factor.matrixL().solve(identity);
    const ComplexMatrix transformed = inverse_lower *
        (0.5 * (numerator + numerator.adjoint())) * inverse_lower.adjoint();
    Eigen::SelfAdjointEigenSolver<ComplexMatrix> eigensolver(transformed);
    if (eigensolver.info() != Eigen::Success) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return std::max(0.0, eigensolver.eigenvalues().maxCoeff());
}

double energy_scaled_minimum_singular_value_approximation(
    const ComplexMatrix &operator_matrix, const ComplexMatrix &energy) {
    const ComplexMatrix hermitian_energy = 0.5 * (energy + energy.adjoint());
    Eigen::LLT<ComplexMatrix> factor(hermitian_energy);
    if (factor.info() != Eigen::Success) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const ComplexMatrix identity = ComplexMatrix::Identity(energy.rows(), energy.cols());
    const ComplexMatrix inverse_lower = factor.matrixL().solve(identity);
    const ComplexMatrix transformed =
        inverse_lower * operator_matrix * inverse_lower.adjoint();
    Eigen::JacobiSVD<ComplexMatrix> svd(transformed);
    if (svd.singularValues().size() == 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return svd.singularValues().minCoeff();
}

VerifiedScalarResult unavailable_result(double approximation,
                                        const std::string &reason,
                                        int precision_bits) {
    VerifiedScalarResult result;
    result.approximation = approximation;
    result.enclosure = {0.0, kInfinity};
    result.metadata.backend_available = false;
    result.metadata.verified = false;
    result.metadata.backend = "Eigen diagnostic only";
    result.metadata.precision_bits = precision_bits;
    result.metadata.rounding_mode = "nearest; no directed rounding";
    result.metadata.failure_reason = reason;
    return result;
}

#ifdef LOD_ENABLE_VERIFIED_CERTIFICATES

class MpfiValue {
  public:
    explicit MpfiValue(mpfr_prec_t precision) {
        mpfi_init2(value_, precision);
        mpfi_set_d(value_, 0.0);
    }

    MpfiValue(const MpfiValue &other) {
        mpfi_init2(value_, mpfi_get_prec(other.value_));
        mpfi_set(value_, other.value_);
    }

    MpfiValue(MpfiValue &&other) noexcept {
        mpfi_init2(value_, mpfi_get_prec(other.value_));
        mpfi_set(value_, other.value_);
    }

    MpfiValue &operator=(const MpfiValue &other) {
        if (this != &other) {
            mpfi_set(value_, other.value_);
        }
        return *this;
    }

    MpfiValue &operator=(MpfiValue &&other) noexcept {
        if (this != &other) {
            mpfi_set(value_, other.value_);
        }
        return *this;
    }

    ~MpfiValue() { mpfi_clear(value_); }

    mpfi_ptr get() { return value_; }
    mpfi_srcptr get() const { return value_; }

  private:
    mpfi_t value_;
};

class IntervalMatrix {
  public:
    IntervalMatrix(Eigen::Index rows, Eigen::Index cols, mpfr_prec_t precision)
        : rows_(rows), cols_(cols), precision_(precision) {
        entries_.reserve(static_cast<std::size_t>(rows * cols));
        for (Eigen::Index k = 0; k < rows * cols; ++k) {
            entries_.emplace_back(precision);
        }
    }

    Eigen::Index rows() const { return rows_; }
    Eigen::Index cols() const { return cols_; }
    mpfr_prec_t precision() const { return precision_; }

    MpfiValue &operator()(Eigen::Index row, Eigen::Index col) {
        return entries_[static_cast<std::size_t>(row * cols_ + col)];
    }

    const MpfiValue &operator()(Eigen::Index row, Eigen::Index col) const {
        return entries_[static_cast<std::size_t>(row * cols_ + col)];
    }

  private:
    Eigen::Index rows_;
    Eigen::Index cols_;
    mpfr_prec_t precision_;
    std::vector<MpfiValue> entries_;
};

mpfr_prec_t checked_precision(int precision_bits) {
    return static_cast<mpfr_prec_t>(std::max(64, precision_bits));
}

double interval_lower(const MpfiValue &value) {
    mpfr_t endpoint;
    mpfr_init2(endpoint, mpfi_get_prec(value.get()));
    mpfi_get_left(endpoint, value.get());
    const double result = mpfr_get_d(endpoint, MPFR_RNDD);
    mpfr_clear(endpoint);
    return result;
}

double interval_upper(const MpfiValue &value) {
    mpfr_t endpoint;
    mpfr_init2(endpoint, mpfi_get_prec(value.get()));
    mpfi_get_right(endpoint, value.get());
    const double result = mpfr_get_d(endpoint, MPFR_RNDU);
    mpfr_clear(endpoint);
    return result;
}

bool finite_interval(const MpfiValue &value) {
    return mpfi_bounded_p(value.get()) != 0 && mpfi_nan_p(value.get()) == 0;
}

IntervalMatrix interval_matrix(const Eigen::MatrixXd &midpoint,
                               const Eigen::MatrixXd &radius,
                               mpfr_prec_t precision) {
    IntervalMatrix result(midpoint.rows(), midpoint.cols(), precision);
    MpfiValue perturbation(precision);
    for (Eigen::Index i = 0; i < midpoint.rows(); ++i) {
        for (Eigen::Index j = 0; j < midpoint.cols(); ++j) {
            mpfi_set_d(result(i, j).get(), midpoint(i, j));
            if (radius(i, j) > 0.0) {
                mpfi_interv_d(perturbation.get(), -radius(i, j), radius(i, j));
                mpfi_add(result(i, j).get(), result(i, j).get(),
                         perturbation.get());
            }
        }
    }
    return result;
}

IntervalMatrix point_interval_matrix(const Eigen::MatrixXd &matrix,
                                     mpfr_prec_t precision) {
    return interval_matrix(matrix, Eigen::MatrixXd::Zero(matrix.rows(), matrix.cols()),
                           precision);
}

IntervalMatrix interval_identity(Eigen::Index size, mpfr_prec_t precision) {
    return point_interval_matrix(Eigen::MatrixXd::Identity(size, size), precision);
}

IntervalMatrix multiply(const IntervalMatrix &left,
                        const IntervalMatrix &right) {
    if (left.cols() != right.rows() || left.precision() != right.precision()) {
        throw std::invalid_argument("incompatible interval matrix multiplication");
    }
    IntervalMatrix result(left.rows(), right.cols(), left.precision());
    MpfiValue product(left.precision());
    for (Eigen::Index i = 0; i < left.rows(); ++i) {
        for (Eigen::Index j = 0; j < right.cols(); ++j) {
            mpfi_set_d(result(i, j).get(), 0.0);
            for (Eigen::Index k = 0; k < left.cols(); ++k) {
                mpfi_mul(product.get(), left(i, k).get(), right(k, j).get());
                mpfi_add(result(i, j).get(), result(i, j).get(), product.get());
            }
        }
    }
    return result;
}

IntervalMatrix subtract(const IntervalMatrix &left,
                        const IntervalMatrix &right) {
    if (left.rows() != right.rows() || left.cols() != right.cols() ||
        left.precision() != right.precision()) {
        throw std::invalid_argument("incompatible interval matrix subtraction");
    }
    IntervalMatrix result(left.rows(), left.cols(), left.precision());
    for (Eigen::Index i = 0; i < left.rows(); ++i) {
        for (Eigen::Index j = 0; j < left.cols(); ++j) {
            mpfi_sub(result(i, j).get(), left(i, j).get(), right(i, j).get());
        }
    }
    return result;
}

IntervalMatrix congruence(const Eigen::MatrixXd &transform,
                          const IntervalMatrix &matrix) {
    const auto left = point_interval_matrix(transform, matrix.precision());
    const auto right = point_interval_matrix(transform.transpose(), matrix.precision());
    return multiply(multiply(left, matrix), right);
}

IntervalMatrix quadratic_form(const Eigen::VectorXd &vector,
                              const IntervalMatrix &matrix) {
    Eigen::MatrixXd row(1, vector.size());
    row.row(0) = vector.transpose();
    Eigen::MatrixXd column(vector.size(), 1);
    column.col(0) = vector;
    return multiply(
        multiply(point_interval_matrix(row, matrix.precision()), matrix),
        point_interval_matrix(column, matrix.precision()));
}

double gershgorin_upper(const IntervalMatrix &matrix) {
    MpfiValue off_diagonal(matrix.precision());
    MpfiValue absolute_entry(matrix.precision());
    MpfiValue disc(matrix.precision());
    double upper = -kInfinity;
    for (Eigen::Index i = 0; i < matrix.rows(); ++i) {
        mpfi_set_d(off_diagonal.get(), 0.0);
        for (Eigen::Index j = 0; j < matrix.cols(); ++j) {
            if (i == j) {
                continue;
            }
            mpfi_abs(absolute_entry.get(), matrix(i, j).get());
            mpfi_add(off_diagonal.get(), off_diagonal.get(), absolute_entry.get());
        }
        mpfi_add(disc.get(), matrix(i, i).get(), off_diagonal.get());
        upper = std::max(upper, interval_upper(disc));
    }
    return upper;
}

double gershgorin_lower(const IntervalMatrix &matrix) {
    MpfiValue off_diagonal(matrix.precision());
    MpfiValue absolute_entry(matrix.precision());
    MpfiValue disc(matrix.precision());
    double lower = kInfinity;
    for (Eigen::Index i = 0; i < matrix.rows(); ++i) {
        mpfi_set_d(off_diagonal.get(), 0.0);
        for (Eigen::Index j = 0; j < matrix.cols(); ++j) {
            if (i == j) {
                continue;
            }
            mpfi_abs(absolute_entry.get(), matrix(i, j).get());
            mpfi_add(off_diagonal.get(), off_diagonal.get(), absolute_entry.get());
        }
        mpfi_sub(disc.get(), matrix(i, i).get(), off_diagonal.get());
        lower = std::min(lower, interval_lower(disc));
    }
    return lower;
}

double matrix_infinity_norm_upper(const IntervalMatrix &matrix) {
    MpfiValue row_sum(matrix.precision());
    MpfiValue absolute_entry(matrix.precision());
    double norm = 0.0;
    for (Eigen::Index i = 0; i < matrix.rows(); ++i) {
        mpfi_set_d(row_sum.get(), 0.0);
        for (Eigen::Index j = 0; j < matrix.cols(); ++j) {
            mpfi_abs(absolute_entry.get(), matrix(i, j).get());
            mpfi_add(row_sum.get(), row_sum.get(), absolute_entry.get());
        }
        norm = std::max(norm, interval_upper(row_sum));
    }
    return norm;
}

double matrix_one_norm_upper(const IntervalMatrix &matrix) {
    MpfiValue column_sum(matrix.precision());
    MpfiValue absolute_entry(matrix.precision());
    double norm = 0.0;
    for (Eigen::Index j = 0; j < matrix.cols(); ++j) {
        mpfi_set_d(column_sum.get(), 0.0);
        for (Eigen::Index i = 0; i < matrix.rows(); ++i) {
            mpfi_abs(absolute_entry.get(), matrix(i, j).get());
            mpfi_add(column_sum.get(), column_sum.get(), absolute_entry.get());
        }
        norm = std::max(norm, interval_upper(column_sum));
    }
    return norm;
}

double directed_sqrt_product_upper(double left, double right,
                                   mpfr_prec_t precision) {
    mpfr_t a;
    mpfr_t b;
    mpfr_t product;
    mpfr_init2(a, precision);
    mpfr_init2(b, precision);
    mpfr_init2(product, precision);
    mpfr_set_d(a, left, MPFR_RNDU);
    mpfr_set_d(b, right, MPFR_RNDU);
    mpfr_mul(product, a, b, MPFR_RNDU);
    mpfr_sqrt(product, product, MPFR_RNDU);
    const double result = mpfr_get_d(product, MPFR_RNDU);
    mpfr_clear(product);
    mpfr_clear(b);
    mpfr_clear(a);
    return result;
}

double spectral_norm_upper(const IntervalMatrix &matrix) {
    return directed_sqrt_product_upper(matrix_one_norm_upper(matrix),
                                       matrix_infinity_norm_upper(matrix),
                                       matrix.precision());
}

double directed_division_upper(double numerator, double denominator,
                               mpfr_prec_t precision) {
    mpfr_t a;
    mpfr_t b;
    mpfr_init2(a, precision);
    mpfr_init2(b, precision);
    mpfr_set_d(a, numerator, MPFR_RNDU);
    mpfr_set_d(b, denominator, MPFR_RNDD);
    mpfr_div(a, a, b, MPFR_RNDU);
    const double result = mpfr_get_d(a, MPFR_RNDU);
    mpfr_clear(b);
    mpfr_clear(a);
    return result;
}

double directed_division_lower(double numerator, double denominator,
                               mpfr_prec_t precision) {
    mpfr_t a;
    mpfr_t b;
    mpfr_init2(a, precision);
    mpfr_init2(b, precision);
    mpfr_set_d(a, numerator, MPFR_RNDD);
    mpfr_set_d(b, denominator, MPFR_RNDU);
    mpfr_div(a, a, b, MPFR_RNDD);
    const double result = mpfr_get_d(a, MPFR_RNDD);
    mpfr_clear(b);
    mpfr_clear(a);
    return result;
}

double inverse_residual_lower_bound(double residual_norm, double inverse_norm,
                                    double energy_upper,
                                    mpfr_prec_t precision) {
    mpfr_t residual;
    mpfr_t numerator;
    mpfr_t denominator;
    mpfr_init2(residual, precision);
    mpfr_init2(numerator, precision);
    mpfr_init2(denominator, precision);
    mpfr_set_d(residual, residual_norm, MPFR_RNDU);
    mpfr_set_ui(numerator, 1, MPFR_RNDD);
    mpfr_sub(numerator, numerator, residual, MPFR_RNDD);
    mpfr_set_d(denominator, inverse_norm, MPFR_RNDU);
    mpfr_mul_d(denominator, denominator, energy_upper, MPFR_RNDU);
    mpfr_div(numerator, numerator, denominator, MPFR_RNDD);
    const double result = mpfr_get_d(numerator, MPFR_RNDD);
    mpfr_clear(denominator);
    mpfr_clear(numerator);
    mpfr_clear(residual);
    return result;
}

bool all_finite(const IntervalMatrix &matrix) {
    for (Eigen::Index i = 0; i < matrix.rows(); ++i) {
        for (Eigen::Index j = 0; j < matrix.cols(); ++j) {
            if (!finite_interval(matrix(i, j))) {
                return false;
            }
        }
    }
    return true;
}

VerificationMetadata verified_metadata(int precision_bits, bool inputs_verified,
                                       const std::string &failure_reason) {
    VerificationMetadata metadata;
    metadata.backend_available = true;
    metadata.verified = inputs_verified && failure_reason.empty();
    metadata.backend = "MPFI interval arithmetic with MPFR endpoints";
    metadata.precision_bits = std::max(64, precision_bits);
    metadata.rounding_mode = "MPFI outward; MPFR RNDD/RNDU";
    if (!failure_reason.empty()) {
        metadata.failure_reason = failure_reason;
    } else if (!inputs_verified) {
        metadata.failure_reason =
            "input entry radii are diagnostic and do not verify upstream errors";
    }
    return metadata;
}

#endif

} // namespace

bool verified_spectrum_backend_available() {
#ifdef LOD_ENABLE_VERIFIED_CERTIFICATES
    return true;
#else
    return false;
#endif
}

VerifiedScalarResult verified_generalized_largest_eigenvalue(
    const ComplexMatrixEnclosure &hermitian_numerator,
    const ComplexMatrixEnclosure &hermitian_positive_denominator,
    int precision_bits) {
    std::string reason;
    if (!valid_enclosure(hermitian_numerator, true, reason) ||
        !valid_enclosure(hermitian_positive_denominator, true, reason) ||
        hermitian_numerator.midpoint.rows() !=
            hermitian_positive_denominator.midpoint.rows()) {
        if (reason.empty()) {
            reason = "generalized eigenvalue matrices have incompatible dimensions";
        }
        return unavailable_result(std::numeric_limits<double>::quiet_NaN(), reason,
                                  precision_bits);
    }
    if (!encloses_hermitian_matrix(hermitian_numerator, reason) ||
        !encloses_hermitian_matrix(hermitian_positive_denominator, reason)) {
        return unavailable_result(std::numeric_limits<double>::quiet_NaN(), reason,
                                  precision_bits);
    }

    const double approximation = generalized_largest_eigenvalue_approximation(
        hermitian_numerator.midpoint,
        hermitian_positive_denominator.midpoint);

#ifndef LOD_ENABLE_VERIFIED_CERTIFICATES
    return unavailable_result(
        approximation,
        "LOD_ENABLE_VERIFIED_CERTIFICATES is OFF; Eigen values are diagnostics only",
        precision_bits);
#else
    VerifiedScalarResult result;
    result.approximation = approximation;
    result.enclosure = {0.0, kInfinity};
    if (!hermitian_numerator.entries_verified
        || !hermitian_positive_denominator.entries_verified) {
        result.metadata = verified_metadata(precision_bits, false, "");
        return result;
    }
    const auto precision = checked_precision(precision_bits);
    try {
        const Eigen::MatrixXd numerator_midpoint =
            complex_real_block(hermitian_numerator.midpoint);
        const Eigen::MatrixXd numerator_radius =
            complex_real_block_radius(hermitian_numerator.radius);
        const Eigen::MatrixXd denominator_midpoint =
            complex_real_block(hermitian_positive_denominator.midpoint);
        const Eigen::MatrixXd denominator_radius = complex_real_block_radius(
            hermitian_positive_denominator.radius);

        Eigen::LLT<Eigen::MatrixXd> factor(
            0.5 * (denominator_midpoint + denominator_midpoint.transpose()));
        if (factor.info() != Eigen::Success) {
            result.metadata = verified_metadata(
                precision_bits, false,
                "midpoint energy Cholesky factorization failed");
            return result;
        }
        const Eigen::MatrixXd inverse_lower = factor.matrixL().solve(
            Eigen::MatrixXd::Identity(denominator_midpoint.rows(),
                                      denominator_midpoint.cols()));

        const auto numerator_interval = interval_matrix(
            numerator_midpoint, numerator_radius, precision);
        const auto denominator_interval = interval_matrix(
            denominator_midpoint, denominator_radius, precision);
        const auto transformed_numerator =
            congruence(inverse_lower, numerator_interval);
        const auto transformed_denominator =
            congruence(inverse_lower, denominator_interval);
        if (!all_finite(transformed_numerator) ||
            !all_finite(transformed_denominator)) {
            result.metadata = verified_metadata(
                precision_bits, false, "non-finite interval during congruence");
            return result;
        }

        const double denominator_lower =
            gershgorin_lower(transformed_denominator);
        if (!(denominator_lower > 0.0)) {
            result.metadata = verified_metadata(
                precision_bits, false,
                "positive definiteness was not verified by the energy enclosure");
            return result;
        }
        const double numerator_upper =
            std::max(0.0, gershgorin_upper(transformed_numerator));
        double upper = directed_division_upper(numerator_upper,
                                               denominator_lower, precision);
        if (std::isfinite(approximation)) {
            upper = std::max(upper,
                             std::nextafter(approximation, kInfinity));
        }
        double lower = 0.0;
        const Eigen::MatrixXd transformed_numerator_midpoint =
            inverse_lower * numerator_midpoint * inverse_lower.transpose();
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(
            0.5 * (transformed_numerator_midpoint +
                   transformed_numerator_midpoint.transpose()));
        if (eigensolver.info() == Eigen::Success) {
            const Eigen::VectorXd trial =
                eigensolver.eigenvectors().rightCols(1);
            const auto numerator_rayleigh =
                quadratic_form(trial, transformed_numerator);
            const auto denominator_rayleigh =
                quadratic_form(trial, transformed_denominator);
            const double rayleigh_numerator_lower =
                interval_lower(numerator_rayleigh(0, 0));
            const double rayleigh_denominator_upper =
                interval_upper(denominator_rayleigh(0, 0));
            if (rayleigh_numerator_lower > 0.0 &&
                rayleigh_denominator_upper > 0.0) {
                lower = directed_division_lower(
                    rayleigh_numerator_lower,
                    rayleigh_denominator_upper, precision);
                lower = std::max(0.0, lower);
                if (std::isfinite(approximation)) {
                    lower = std::min(
                        lower, std::nextafter(approximation, 0.0));
                }
            }
        }
        result.enclosure = {lower, upper};
        result.verification_residual_bound =
            std::max(0.0, 1.0 - denominator_lower);
        const bool inputs_verified = hermitian_numerator.entries_verified &&
                                     hermitian_positive_denominator.entries_verified;
        result.metadata = verified_metadata(precision_bits, inputs_verified, "");
        return result;
    } catch (const std::exception &error) {
        result.metadata = verified_metadata(precision_bits, false, error.what());
        return result;
    }
#endif
}

VerifiedScalarResult verified_energy_scaled_minimum_singular_value(
    const ComplexMatrixEnclosure &operator_matrix,
    const ComplexMatrixEnclosure &hermitian_positive_energy,
    int precision_bits) {
    std::string reason;
    if (!valid_enclosure(operator_matrix, true, reason) ||
        !valid_enclosure(hermitian_positive_energy, true, reason) ||
        operator_matrix.midpoint.rows() !=
            hermitian_positive_energy.midpoint.rows()) {
        if (reason.empty()) {
            reason = "operator and energy matrices have incompatible dimensions";
        }
        return unavailable_result(std::numeric_limits<double>::quiet_NaN(), reason,
                                  precision_bits);
    }
    if (!encloses_hermitian_matrix(hermitian_positive_energy, reason)) {
        return unavailable_result(std::numeric_limits<double>::quiet_NaN(), reason,
                                  precision_bits);
    }

    const double approximation = energy_scaled_minimum_singular_value_approximation(
        operator_matrix.midpoint, hermitian_positive_energy.midpoint);

#ifndef LOD_ENABLE_VERIFIED_CERTIFICATES
    return unavailable_result(
        approximation,
        "LOD_ENABLE_VERIFIED_CERTIFICATES is OFF; Eigen values are diagnostics only",
        precision_bits);
#else
    VerifiedScalarResult result;
    result.approximation = approximation;
    result.enclosure = {0.0, kInfinity};
    if (!operator_matrix.entries_verified
        || !hermitian_positive_energy.entries_verified) {
        result.metadata = verified_metadata(precision_bits, false, "");
        return result;
    }
    const auto precision = checked_precision(precision_bits);
    try {
        const Eigen::MatrixXd operator_midpoint =
            complex_real_block(operator_matrix.midpoint);
        const Eigen::MatrixXd operator_radius =
            complex_real_block_radius(operator_matrix.radius);
        const Eigen::MatrixXd energy_midpoint =
            complex_real_block(hermitian_positive_energy.midpoint);
        const Eigen::MatrixXd energy_radius =
            complex_real_block_radius(hermitian_positive_energy.radius);

        Eigen::LLT<Eigen::MatrixXd> factor(
            0.5 * (energy_midpoint + energy_midpoint.transpose()));
        if (factor.info() != Eigen::Success) {
            result.metadata = verified_metadata(
                precision_bits, false,
                "midpoint energy Cholesky factorization failed");
            return result;
        }
        const Eigen::MatrixXd inverse_lower = factor.matrixL().solve(
            Eigen::MatrixXd::Identity(energy_midpoint.rows(),
                                      energy_midpoint.cols()));

        const auto operator_interval = interval_matrix(
            operator_midpoint, operator_radius, precision);
        const auto energy_interval = interval_matrix(
            energy_midpoint, energy_radius, precision);
        const auto transformed_operator =
            congruence(inverse_lower, operator_interval);
        const auto transformed_energy =
            congruence(inverse_lower, energy_interval);
        if (!all_finite(transformed_operator) ||
            !all_finite(transformed_energy)) {
            result.metadata = verified_metadata(
                precision_bits, false, "non-finite interval during congruence");
            return result;
        }

        const double energy_lower = gershgorin_lower(transformed_energy);
        const double energy_upper = gershgorin_upper(transformed_energy);
        if (!(energy_lower > 0.0) || !(energy_upper > 0.0)) {
            result.metadata = verified_metadata(
                precision_bits, false,
                "positive definiteness was not verified by the energy enclosure");
            return result;
        }

        const Eigen::MatrixXd transformed_midpoint =
            inverse_lower * operator_midpoint * inverse_lower.transpose();
        Eigen::FullPivLU<Eigen::MatrixXd> inverse_solver(transformed_midpoint);
        if (!inverse_solver.isInvertible()) {
            result.metadata = verified_metadata(
                precision_bits, false, "midpoint operator is singular");
            return result;
        }
        const Eigen::MatrixXd approximate_inverse = inverse_solver.inverse();
        const auto inverse_interval =
            point_interval_matrix(approximate_inverse, precision);
        const auto residual_interval = subtract(
            interval_identity(transformed_operator.rows(), precision),
            multiply(inverse_interval, transformed_operator));
        const double residual_norm = spectral_norm_upper(residual_interval);
        result.verification_residual_bound = residual_norm;
        if (!(residual_norm < 1.0)) {
            result.metadata = verified_metadata(
                precision_bits, false,
                "inverse residual norm is not strictly below one");
            return result;
        }
        const double inverse_norm = spectral_norm_upper(inverse_interval);
        if (!(inverse_norm > 0.0) || !std::isfinite(inverse_norm) ||
            !std::isfinite(energy_upper)) {
            result.metadata = verified_metadata(
                precision_bits, false, "invalid inverse or energy norm enclosure");
            return result;
        }
        double lower = inverse_residual_lower_bound(
            residual_norm, inverse_norm, energy_upper, precision);
        lower = std::max(0.0, lower);
        if (std::isfinite(approximation)) {
            lower = std::min(lower,
                             std::nextafter(approximation, 0.0));
        }
        result.enclosure = {lower, kInfinity};
        const bool inputs_verified = operator_matrix.entries_verified &&
                                     hermitian_positive_energy.entries_verified;
        result.metadata = verified_metadata(precision_bits, inputs_verified, "");
        return result;
    } catch (const std::exception &error) {
        result.metadata = verified_metadata(precision_bits, false, error.what());
        return result;
    }
#endif
}

VerifiedScalarResult verified_minimum_singular_value(
    const ComplexMatrixEnclosure &operator_matrix, int precision_bits) {
    ComplexMatrixEnclosure identity;
    identity.midpoint = ComplexMatrix::Identity(operator_matrix.midpoint.rows(),
                                                operator_matrix.midpoint.rows());
    identity.radius = Eigen::MatrixXd::Zero(operator_matrix.midpoint.rows(),
                                            operator_matrix.midpoint.rows());
    identity.entries_verified = true;
    return verified_energy_scaled_minimum_singular_value(
        operator_matrix, identity, precision_bits);
}

} // namespace lod2d::solver

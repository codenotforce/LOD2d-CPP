#pragma once

#include "helmholtz/types.h"

#include <Eigen/Dense>

#include <string>

namespace lod2d::solver {

struct ScalarInterval {
    double lower = 0.0;
    double upper = 0.0;

    bool contains(double value) const {
        return value >= lower && value <= upper;
    }
};

struct VerificationMetadata {
    bool backend_available = false;
    bool verified = false;
    std::string backend;
    int precision_bits = 0;
    std::string rounding_mode;
    std::string failure_reason;
};

// One nonnegative radius encloses both the real and imaginary component of
// each midpoint entry. entries_verified must be true only when these radii
// rigorously include assembly and upstream solve errors.
struct ComplexMatrixEnclosure {
    helmholtz::ComplexMatrix midpoint;
    Eigen::MatrixXd radius;
    bool entries_verified = false;
};

struct VerifiedScalarResult {
    double approximation = 0.0;
    ScalarInterval enclosure;
    double verification_residual_bound = 0.0;
    VerificationMetadata metadata;
};

bool verified_spectrum_backend_available();

VerifiedScalarResult verified_generalized_largest_eigenvalue(
    const ComplexMatrixEnclosure &hermitian_numerator,
    const ComplexMatrixEnclosure &hermitian_positive_denominator,
    int precision_bits = 128);

VerifiedScalarResult verified_energy_scaled_minimum_singular_value(
    const ComplexMatrixEnclosure &operator_matrix,
    const ComplexMatrixEnclosure &hermitian_positive_energy,
    int precision_bits = 128);

VerifiedScalarResult verified_minimum_singular_value(
    const ComplexMatrixEnclosure &operator_matrix,
    int precision_bits = 128);

} // namespace lod2d::solver

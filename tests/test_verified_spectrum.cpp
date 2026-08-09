#include "helmholtz/adaptive/verified_spectrum.h"

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace lod2d;
using namespace lod2d::helmholtz;
using namespace lod2d::solver;

namespace {

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

void require_close(double first, double second, double tolerance,
                   const char *message) {
    const double scale = std::max({1.0, std::abs(first), std::abs(second)});
    if (std::abs(first - second) > tolerance * scale)
        throw std::runtime_error(message);
}

ComplexMatrixEnclosure enclosure(const ComplexMatrix &matrix,
                                 double radius,
                                 bool verified = true) {
    ComplexMatrixEnclosure result;
    result.midpoint = matrix;
    result.radius = Eigen::MatrixXd::Constant(
        matrix.rows(), matrix.cols(), radius);
    result.entries_verified = verified;
    return result;
}

double generalized_reference(const ComplexMatrix &numerator,
                             const ComplexMatrix &denominator) {
    Eigen::LLT<ComplexMatrix> factor(denominator);
    const ComplexMatrix inverse_lower = factor.matrixL().solve(
        ComplexMatrix::Identity(denominator.rows(), denominator.cols()));
    const ComplexMatrix transformed =
        inverse_lower * numerator * inverse_lower.adjoint();
    Eigen::SelfAdjointEigenSolver<ComplexMatrix> solver(transformed);
    return solver.eigenvalues().maxCoeff();
}

double gamma_reference(const ComplexMatrix &matrix,
                       const ComplexMatrix &energy) {
    Eigen::LLT<ComplexMatrix> factor(energy);
    const ComplexMatrix inverse_lower = factor.matrixL().solve(
        ComplexMatrix::Identity(energy.rows(), energy.cols()));
    Eigen::JacobiSVD<ComplexMatrix> svd(
        inverse_lower * matrix * inverse_lower.adjoint());
    return svd.singularValues().minCoeff();
}

void verify_directional_enclosures() {
    ComplexMatrix energy(2, 2);
    energy << Complex(2.0, 0.0), Complex(0.25, 0.1),
              Complex(0.25, -0.1), Complex(1.4, 0.0);
    ComplexMatrix gram(2, 2);
    gram << Complex(3.0, 0.0), Complex(0.2, -0.15),
            Complex(0.2, 0.15), Complex(0.8, 0.0);
    ComplexMatrix system(2, 2);
    system << Complex(1.7, -0.2), Complex(-0.3, 0.4),
              Complex(0.25, -0.1), Complex(1.1, 0.3);

    const double lambda = generalized_reference(gram, energy);
    const double gamma = gamma_reference(system, energy);
    const auto lambda_128 = verified_generalized_largest_eigenvalue(
        enclosure(gram, 2e-16), enclosure(energy, 2e-16), 128);
    const auto gamma_128 = verified_energy_scaled_minimum_singular_value(
        enclosure(system, 2e-16), enclosure(energy, 2e-16), 128);

    require_close(lambda_128.approximation, lambda, 2e-13,
                  "generalized eigenvalue diagnostic disagrees with dense reference");
    require_close(gamma_128.approximation, gamma, 2e-13,
                  "energy-scaled singular diagnostic disagrees with dense reference");

    if (!verified_spectrum_backend_available()) {
        require(!lambda_128.metadata.verified,
                "Eigen-only generalized eigenvalue was marked verified");
        require(!gamma_128.metadata.verified,
                "Eigen-only singular value was marked verified");
        require(std::isinf(lambda_128.enclosure.upper),
                "Eigen-only path emitted a finite certified eigenvalue bound");
        require(gamma_128.enclosure.lower == 0.0,
                "Eigen-only path emitted a positive certified singular lower bound");
        return;
    }

    require(lambda_128.metadata.verified,
            "MPFR/MPFI generalized eigenvalue enclosure is not verified");
    require(gamma_128.metadata.verified,
            "MPFR/MPFI singular lower enclosure is not verified");
    require(lambda_128.enclosure.lower <= lambda
                && lambda <= lambda_128.enclosure.upper,
            "verified generalized eigenvalue enclosure has the wrong direction");
    require(gamma_128.enclosure.lower <= gamma,
            "verified singular lower enclosure exceeds the reference value");
    require(lambda_128.enclosure.lower > 0.0,
            "verified Rayleigh lower enclosure is unexpectedly trivial");
    require(gamma_128.enclosure.lower > 0.0,
            "verified inverse-residual singular lower enclosure is trivial");

    const auto lambda_256 = verified_generalized_largest_eigenvalue(
        enclosure(gram, 2e-16), enclosure(energy, 2e-16), 256);
    const auto gamma_256 = verified_energy_scaled_minimum_singular_value(
        enclosure(system, 2e-16), enclosure(energy, 2e-16), 256);
    require(lambda_256.metadata.verified && gamma_256.metadata.verified,
            "higher precision verification unexpectedly failed");
    std::cout << "lambda[128]=[" << lambda_128.enclosure.lower << ','
              << lambda_128.enclosure.upper << "] lambda[256]=["
              << lambda_256.enclosure.lower << ','
              << lambda_256.enclosure.upper << "] gamma_lower[128/256]="
              << gamma_128.enclosure.lower << '/'
              << gamma_256.enclosure.lower << '\n';
    require(lambda_256.enclosure.lower >= lambda_128.enclosure.lower,
            "higher precision generalized lower bound is not nested");
    require(lambda_256.enclosure.upper <= lambda_128.enclosure.upper,
            "higher precision generalized upper bound is not nested");
    require(gamma_256.enclosure.lower >= gamma_128.enclosure.lower,
            "higher precision singular lower bound is not nested");

    auto unverified_gram = enclosure(gram, 2e-16, false);
    const auto conditional = verified_generalized_largest_eigenvalue(
        unverified_gram, enclosure(energy, 2e-16), 128);
    require(!conditional.metadata.verified,
            "unverified matrix entries were promoted to a verified spectrum");
}

} // namespace

int main() {
    try {
        verify_directional_enclosures();
        std::cout << "Verified spectrum direction and backend gates passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_verified_spectrum failed: " << error.what() << '\n';
        return 1;
    }
}

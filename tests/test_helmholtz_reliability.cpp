#include "helmholtz/adaptive/reliability.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace lod2d::helmholtz::adaptive;

namespace {

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

} // namespace

int main() {
    try {
        const double rho = estimate_localization_decay(
            {1, 2, 3, 4}, {0.4, 0.2, 0.1, 0.05});
        require(std::abs(rho - 0.5) < 1e-12, "localization decay fit is incorrect");

        ReliabilityConstants constants;
        constants.c1 = 1.0;
        constants.c4 = 0.1;
        ReliabilitySample sample{4.0, 2, 0.125, 0.5};
        const ReliabilityEvaluation valid = evaluate_reliability_envelope(
            constants, sample, rho);
        require(valid.applicable && valid.predicted >= 1.0,
                "valid reliability envelope was rejected");

        constants.c6 = 1.0;
        const ReliabilityEvaluation invalid = evaluate_reliability_envelope(
            constants, sample, rho);
        require(!invalid.applicable, "invalid reliability denominator was accepted");

        const std::vector<ReliabilitySample> training{
            {4.0, 1, 0.125, 1.4},
            {4.0, 2, 0.125, 1.2},
            {8.0, 2, 0.0625, 1.6},
            {8.0, 3, 0.0625, 1.3},
            {16.0, 4, 0.03125, 1.8}};
        const ReliabilityFitResult fit = fit_reliability_envelope(training, rho);
        require(fit.constants.c7 == 1.0, "reliability fit did not enforce C7 normalization");
        require(fit.max_observed_over_predicted <= 1.0 + 1e-9,
                "reliability fit does not envelope its training data");

        require(std::abs(spearman_rank_correlation(
                    {1.0, 2.0, 3.0, 4.0}, {2.0, 4.0, 6.0, 8.0}) - 1.0) < 1e-12,
                "Spearman correlation is incorrect");
        const double overlap = doerfler_energy_overlap({0, 2}, {0, 1}, {4.0, 3.0, 2.0});
        require(std::abs(overlap - 4.0 / 7.0) < 1e-12,
                "Doerfler energy overlap is incorrect");

        std::cout << "Adaptive Helmholtz reliability calibration tools passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_reliability failed: " << error.what() << '\n';
        return 1;
    }
}

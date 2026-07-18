#pragma once

#include <vector>

namespace lod2d::helmholtz::adaptive {

struct ReliabilityConstants {
    double c1 = 0.0;
    double c2 = 0.0;
    double c3 = 0.0;
    double c4 = 0.0;
    double c5 = 0.0;
    double c6 = 0.0;
    double c7 = 1.0;
};

struct ReliabilitySample {
    double wavenumber = 0.0;
    int ell = 0;
    double q_max = 0.0;
    double observed = 0.0;
};

struct ReliabilityEvaluation {
    bool applicable = false;
    double predicted = 0.0;
    double c_ell_k = 0.0;
    double denominator_margin = 0.0;
    double stability_margin = 0.0;
};

struct ReliabilityFitResult {
    ReliabilityConstants constants;
    double objective = 0.0;
    double max_observed_over_predicted = 0.0;
    double median_predicted_over_observed = 0.0;
};

double estimate_localization_decay(
    const std::vector<int> &ell_values,
    const std::vector<double> &corrector_differences);

ReliabilityEvaluation evaluate_reliability_envelope(
    const ReliabilityConstants &constants,
    const ReliabilitySample &sample,
    double rho,
    double epsilon = 0.01,
    double stability_delta = 1e-6);

ReliabilityFitResult fit_reliability_envelope(
    const std::vector<ReliabilitySample> &training_samples,
    double rho,
    double epsilon = 0.01,
    double stability_delta = 1e-6);

double spearman_rank_correlation(
    const std::vector<double> &left,
    const std::vector<double> &right);

double doerfler_energy_overlap(
    const std::vector<int> &marked,
    const std::vector<int> &reference_marked,
    const std::vector<double> &reference_indicator_squared);

} // namespace lod2d::helmholtz::adaptive

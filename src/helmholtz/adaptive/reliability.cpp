#include "helmholtz/adaptive/reliability.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <unordered_set>

namespace lod2d::helmholtz::adaptive {
namespace {

std::vector<double> average_ranks(const std::vector<double> &values) {
    std::vector<int> order(values.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int left, int right) {
        if (values[left] != values[right]) return values[left] < values[right];
        return left < right;
    });
    std::vector<double> ranks(values.size());
    for (int begin = 0; begin < static_cast<int>(order.size());) {
        int end = begin + 1;
        while (end < static_cast<int>(order.size())
               && values[order[end]] == values[order[begin]]) ++end;
        const double rank = 0.5 * (begin + end - 1);
        for (int i = begin; i < end; ++i) ranks[order[i]] = rank;
        begin = end;
    }
    return ranks;
}

double median(std::vector<double> values) {
    if (values.empty()) return 0.0;
    std::sort(values.begin(), values.end());
    const size_t middle = values.size() / 2;
    if (values.size() % 2 == 1) return values[middle];
    return 0.5 * (values[middle - 1] + values[middle]);
}

struct CandidateScore {
    bool valid = false;
    ReliabilityConstants constants;
    double objective = std::numeric_limits<double>::infinity();
    double max_ratio = std::numeric_limits<double>::infinity();
    double median_slack = std::numeric_limits<double>::infinity();
};

CandidateScore score_candidate(
    ReliabilityConstants constants,
    const std::vector<ReliabilitySample> &samples,
    double rho,
    double epsilon,
    double stability_delta) {
    constants.c1 = 0.0;
    double required_c1 = 0.0;
    for (const ReliabilitySample &sample : samples) {
        ReliabilityEvaluation evaluation = evaluate_reliability_envelope(
            constants, sample, rho, epsilon, stability_delta);
        if (!evaluation.applicable) return {};
        const double stability_factor = (1.0 + evaluation.c_ell_k * evaluation.c_ell_k)
                                      / (1.0 - evaluation.c_ell_k * evaluation.c_ell_k);
        const double rho_term = std::pow(rho, sample.ell);
        const double q_term = std::pow(sample.q_max, 0.5 - epsilon);
        required_c1 = std::max(
            required_c1,
            sample.observed / stability_factor
                - constants.c2 * rho_term - constants.c3 * q_term);
    }
    constants.c1 = std::max(0.0, required_c1) * (1.0 + 1e-12);

    std::vector<double> slacks;
    slacks.reserve(samples.size());
    double log_sum = 0.0;
    double max_ratio = 0.0;
    for (const ReliabilitySample &sample : samples) {
        const ReliabilityEvaluation evaluation = evaluate_reliability_envelope(
            constants, sample, rho, epsilon, stability_delta);
        if (!evaluation.applicable || !(evaluation.predicted > 0.0)) return {};
        const double ratio = sample.observed / evaluation.predicted;
        if (ratio > 1.0 + 1e-9) return {};
        max_ratio = std::max(max_ratio, ratio);
        const double slack = evaluation.predicted / sample.observed;
        slacks.push_back(slack);
        log_sum += std::log(std::max(slack, 1.0));
    }

    CandidateScore result;
    result.valid = true;
    result.constants = constants;
    result.max_ratio = max_ratio;
    result.median_slack = median(slacks);
    result.objective = log_sum / samples.size()
                     + 0.05 * std::log(std::max(1.0, result.median_slack));
    return result;
}

} // namespace

double estimate_localization_decay(
    const std::vector<int> &ell_values,
    const std::vector<double> &corrector_differences) {
    if (ell_values.size() != corrector_differences.size() || ell_values.size() < 2)
        throw std::invalid_argument("decay fit requires matching data with at least two points");
    const double mean_ell = std::accumulate(ell_values.begin(), ell_values.end(), 0.0)
                          / ell_values.size();
    double mean_log = 0.0;
    for (double value : corrector_differences) {
        if (!(value > 0.0) || !std::isfinite(value))
            throw std::invalid_argument("decay differences must be finite and positive");
        mean_log += std::log(value);
    }
    mean_log /= corrector_differences.size();
    double numerator = 0.0;
    double denominator = 0.0;
    for (size_t i = 0; i < ell_values.size(); ++i) {
        const double centered = ell_values[i] - mean_ell;
        numerator += centered * (std::log(corrector_differences[i]) - mean_log);
        denominator += centered * centered;
    }
    if (denominator == 0.0) throw std::invalid_argument("decay fit ell values are constant");
    const double rho = std::exp(numerator / denominator);
    if (!(rho > 0.0 && rho < 1.0))
        throw std::runtime_error("corrector differences do not exhibit geometric decay");
    return rho;
}

ReliabilityEvaluation evaluate_reliability_envelope(
    const ReliabilityConstants &c,
    const ReliabilitySample &sample,
    double rho,
    double epsilon,
    double stability_delta) {
    ReliabilityEvaluation result;
    if (!(sample.wavenumber > 0.0) || sample.ell < 0
        || !(sample.q_max >= 0.0) || !(rho > 0.0 && rho < 1.0)
        || !(epsilon > 0.0 && epsilon < 0.5) || !(c.c7 > 0.0)) return result;
    const double rho_term = std::pow(rho, sample.ell);
    const double q_term = std::pow(sample.q_max, 0.5 - epsilon);
    const double k3 = sample.wavenumber * sample.wavenumber * sample.wavenumber;
    const double perturbation = c.c5 * rho_term + c.c6 * q_term;
    result.denominator_margin = c.c7 - perturbation * k3;
    if (!(result.denominator_margin > 0.0)) return result;
    result.c_ell_k = (c.c4 + perturbation * k3) / result.denominator_margin;
    result.stability_margin = 1.0 - result.c_ell_k;
    if (!(result.c_ell_k >= 0.0 && result.stability_margin > stability_delta)) return result;
    const double prefactor = c.c1 + c.c2 * rho_term + c.c3 * q_term;
    if (!(prefactor >= 0.0)) return result;
    result.predicted = prefactor
        * (1.0 + result.c_ell_k * result.c_ell_k)
        / (1.0 - result.c_ell_k * result.c_ell_k);
    result.applicable = std::isfinite(result.predicted);
    return result;
}

ReliabilityFitResult fit_reliability_envelope(
    const std::vector<ReliabilitySample> &samples,
    double rho,
    double epsilon,
    double stability_delta) {
    if (samples.empty()) throw std::invalid_argument("reliability fit needs training samples");
    double max_observed = 0.0;
    double max_k3 = 0.0;
    for (const ReliabilitySample &sample : samples) {
        if (!(sample.observed > 0.0) || !std::isfinite(sample.observed))
            throw std::invalid_argument("reliability observations must be finite and positive");
        max_observed = std::max(max_observed, sample.observed);
        max_k3 = std::max(max_k3, std::pow(sample.wavenumber, 3));
    }

    const std::vector<double> amplitude{
        0.0, 0.01 * max_observed, 0.05 * max_observed, 0.1 * max_observed,
        0.25 * max_observed, 0.5 * max_observed, max_observed, 2.0 * max_observed};
    const std::vector<double> c4_values{0.0, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8};
    std::vector<double> perturbation_values;
    for (double fraction : {0.0, 1e-4, 1e-3, 0.005, 0.01, 0.025, 0.05, 0.1})
        perturbation_values.push_back(fraction / std::max(max_k3, 1.0));

    ReliabilityConstants current;
    CandidateScore best = score_candidate(
        current, samples, rho, epsilon, stability_delta);
    if (!best.valid) throw std::runtime_error("failed to initialize reliability envelope fit");
    current = best.constants;

    for (int pass = 0; pass < 4; ++pass) {
        const std::array<std::vector<double>, 5> candidates{
            amplitude, amplitude, c4_values, perturbation_values, perturbation_values};
        for (int parameter = 0; parameter < 5; ++parameter) {
            CandidateScore local_best = best;
            for (double value : candidates[parameter]) {
                ReliabilityConstants trial = current;
                if (parameter == 0) trial.c2 = value;
                else if (parameter == 1) trial.c3 = value;
                else if (parameter == 2) trial.c4 = value;
                else if (parameter == 3) trial.c5 = value;
                else trial.c6 = value;
                CandidateScore score = score_candidate(
                    trial, samples, rho, epsilon, stability_delta);
                if (score.valid && score.objective < local_best.objective) local_best = score;
            }
            best = local_best;
            current = best.constants;
        }
    }

    ReliabilityFitResult result;
    result.constants = best.constants;
    result.objective = best.objective;
    result.max_observed_over_predicted = best.max_ratio;
    result.median_predicted_over_observed = best.median_slack;
    return result;
}

double spearman_rank_correlation(
    const std::vector<double> &left,
    const std::vector<double> &right) {
    if (left.size() != right.size() || left.size() < 2)
        throw std::invalid_argument("Spearman correlation requires matching vectors of size at least two");
    const std::vector<double> left_rank = average_ranks(left);
    const std::vector<double> right_rank = average_ranks(right);
    const double left_mean = std::accumulate(left_rank.begin(), left_rank.end(), 0.0) / left_rank.size();
    const double right_mean = std::accumulate(right_rank.begin(), right_rank.end(), 0.0) / right_rank.size();
    double covariance = 0.0;
    double left_variance = 0.0;
    double right_variance = 0.0;
    for (size_t i = 0; i < left_rank.size(); ++i) {
        const double a = left_rank[i] - left_mean;
        const double b = right_rank[i] - right_mean;
        covariance += a * b;
        left_variance += a * a;
        right_variance += b * b;
    }
    if (left_variance == 0.0 || right_variance == 0.0) return 0.0;
    return covariance / std::sqrt(left_variance * right_variance);
}

double doerfler_energy_overlap(
    const std::vector<int> &marked,
    const std::vector<int> &reference_marked,
    const std::vector<double> &reference_indicator_squared) {
    std::unordered_set<int> selected(marked.begin(), marked.end());
    double reference_energy = 0.0;
    double overlap_energy = 0.0;
    for (int element : reference_marked) {
        if (element < 0 || element >= static_cast<int>(reference_indicator_squared.size()))
            throw std::out_of_range("reference marked element is out of range");
        reference_energy += reference_indicator_squared[element];
        if (selected.count(element) != 0)
            overlap_energy += reference_indicator_squared[element];
    }
    return reference_energy > 0.0 ? overlap_energy / reference_energy : 1.0;
}

} // namespace lod2d::helmholtz::adaptive

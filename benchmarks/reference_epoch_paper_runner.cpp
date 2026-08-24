#include "helmholtz/experiments/reference_epoch_runner.h"

#include "helmholtz/adaptive/candidate_dual.h"
#include "helmholtz/adaptive/candidate_flux.h"
#include "helmholtz/adaptive/certificates.h"
#include "helmholtz/adaptive/error_control.h"
#include "helmholtz/adaptive/estimator.h"
#include "helmholtz/adaptive/kernel_residual.h"
#include "helmholtz/adaptive/singularity_hybrid.h"
#include "helmholtz/benchmarks/paper_cases.h"
#include "helmholtz/boundary.h"
#include "helmholtz/experiments/paper_config.h"
#include "helmholtz/model.h"
#include "io/vtk_writer.h"
#include "lod/patches.h"
#include "mesh/refine.h"

#include <Eigen/SparseLU>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <thread>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

namespace lod2d::helmholtz::experiments {
namespace {

using namespace adaptive;
using benchmarks::PaperCaseData;
using Clock = std::chrono::steady_clock;

struct ReferenceEpochMeshSnapshot {
    std::size_t epoch = 0;
    std::size_t H_step = 0;
    std::uint64_t reference_mesh_version = 0;
    std::string stage;
    ReferenceEpochDriverAction anchor_action =
        ReferenceEpochDriverAction::BeginEpoch;
    std::size_t anchor_occurrence = 0;
    TriMesh coarse;
    TriMesh candidate;
    std::optional<TriMesh> reference;
};

struct HybridReserveDiagnostic {
    std::size_t epoch = 0;
    std::size_t H_step = 0;
    std::size_t refresh_index = 0;
    int trial_index = -1;
    std::string row_type;
    std::string status;
    std::string reject_reason;
    int requested_target_gap = 0;
    int collar = -1;
    int maximum_graph_distance = 0;
    int ell = 0;
    int ell_s = -1;
    double physical_radius = std::numeric_limits<double>::quiet_NaN();
    std::size_t omega_s_elements = 0;
    std::size_t omega_f_elements = 0;
    std::size_t full_target_elements = 0;
    bool target_satisfied = false;
    std::size_t matching_spill = 0;
    int profile_margin_before = std::numeric_limits<int>::max();
    int profile_margin_after = std::numeric_limits<int>::max();
    int far_gap_before = std::numeric_limits<int>::max();
    int far_gap_after = std::numeric_limits<int>::max();
    std::size_t candidate_elements = 0;
    std::size_t candidate_nodes = 0;
    std::size_t deepened_elements = 0;
    double time_matching = 0.0;
    double time_probe = 0.0;
    double time_deepen = 0.0;
    double time_total = 0.0;
};

double seconds_since(const Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

std::string read_text(const std::filesystem::path &path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("cannot open " + path.string());
    return std::string(
        std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>());
}

std::string first_token(const std::string &text) {
    std::istringstream input(text);
    std::string token;
    input >> token;
    return token;
}

std::string number(double value) {
    if (!std::isfinite(value)) return "NA";
    std::ostringstream out;
    out << std::setprecision(17) << value;
    return out.str();
}

std::string csv(std::string value) {
    if (value.find_first_of(",\"\r\n") == std::string::npos) return value;
    std::string result = "\"";
    for (char character : value) {
        if (character == '\"') result.push_back('\"');
        result.push_back(character);
    }
    result.push_back('\"');
    return result;
}

std::string json_string(std::string_view value) {
    std::string result = "\"";
    for (char character : value) {
        switch (character) {
        case '\"': result += "\\\""; break;
        case '\\': result += "\\\\"; break;
        case '\n': result += "\\n"; break;
        case '\r': result += "\\r"; break;
        case '\t': result += "\\t"; break;
        default: result.push_back(character); break;
        }
    }
    result.push_back('\"');
    return result;
}

double peak_memory_mb() {
#if defined(__unix__) || defined(__APPLE__)
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
#if defined(__APPLE__)
        return static_cast<double>(usage.ru_maxrss) / (1024.0 * 1024.0);
#else
        return static_cast<double>(usage.ru_maxrss) / 1024.0;
#endif
    }
#endif
    return 0.0;
}

std::vector<int> free_nodes(const TriMesh &mesh) {
    std::vector<char> fixed(mesh.nodes.size(), false);
    for (int node : dirichlet_nodes(mesh)) fixed[node] = true;
    std::vector<int> result;
    for (int node = 0; node < static_cast<int>(mesh.nodes.size()); ++node)
        if (!fixed[node]) result.push_back(node);
    return result;
}

std::size_t unknowns(const TriMesh &mesh) {
    return mesh.nodes.size() - dirichlet_nodes(mesh).size();
}

void verify_discrete_helmholtz_factorization(
    const HelmholtzOperators &operators) {
    const int global_size = operators.system.rows();
    if (operators.system.cols() != global_size)
        throw std::invalid_argument(
            "Helmholtz well-posedness audit received a nonsquare operator");
    std::vector<char> fixed(global_size, false);
    for (const int node : operators.dirichlet_nodes) {
        if (node < 0 || node >= global_size)
            throw std::invalid_argument(
                "Helmholtz well-posedness audit has an invalid Dirichlet node");
        fixed[node] = true;
    }
    std::vector<int> global_to_free(global_size, -1);
    int free_count = 0;
    for (int node = 0; node < global_size; ++node) {
        if (!fixed[node]) global_to_free[node] = free_count++;
    }
    if (free_count == 0)
        throw std::runtime_error(
            "Helmholtz well-posedness audit has no free degrees of freedom");
    std::vector<ComplexTriplet> triplets;
    triplets.reserve(operators.system.nonZeros());
    for (int global_column = 0; global_column < global_size; ++global_column) {
        const int local_column = global_to_free[global_column];
        if (local_column < 0) continue;
        for (ComplexSparseMatrix::InnerIterator it(
                 operators.system, global_column); it; ++it) {
            const int local_row = global_to_free[it.row()];
            if (local_row >= 0)
                triplets.emplace_back(local_row, local_column, it.value());
        }
    }
    ComplexSparseMatrix reduced(free_count, free_count);
    reduced.setFromTriplets(triplets.begin(), triplets.end());
    reduced.makeCompressed();
    Eigen::SparseLU<ComplexSparseMatrix> solver;
    solver.analyzePattern(reduced);
    solver.factorize(reduced);
    if (solver.info() != Eigen::Success)
        throw std::runtime_error(
            "candidate Helmholtz Galerkin factorization failed");
}

int prospective_reference_level_gap(
    const ReferenceEpochHierarchy &hierarchy,
    const std::vector<char> *included_coarse_elements = nullptr) {
    if (hierarchy.has_proposed_coarse_refinement()) {
        if (included_coarse_elements)
            return hierarchy.minimum_proposed_reference_level_gap(
                *included_coarse_elements);
        return hierarchy.minimum_proposed_reference_level_gap();
    }
    const std::vector<int> *coarse_levels = &hierarchy.coarse_levels();
    const std::vector<int> *parents =
        &hierarchy.reference_parent_coarse_elements();
    const std::vector<int> &reference_levels =
        hierarchy.reference_element_levels();
    if (reference_levels.size() != parents->size()) {
        throw std::logic_error("reference level-gap metadata is incomplete");
    }
    int gap = std::numeric_limits<int>::max();
    for (std::size_t child = 0; child < parents->size(); ++child) {
        const int parent = (*parents)[child];
        if (parent < 0
            || static_cast<std::size_t>(parent) >= coarse_levels->size()) {
            throw std::logic_error("reference level-gap parent is invalid");
        }
        if (included_coarse_elements
            && !(*included_coarse_elements)[parent]) continue;
        const int local_gap =
            reference_levels[child] - (*coarse_levels)[parent];
        gap = std::min(gap, local_gap);
    }
    return gap;
}

struct ClosureAwareHybridMarking {
    std::vector<int> marked_elements;
    std::optional<ReferenceEpochCoarseRefinementPreview> accepted_preview;
    double regular_mass = 0.0;
    double admissible_mass = 0.0;
    double marked_mass = 0.0;
    double time_preview = 0.0;
    double time_preview_nvb = 0.0;
    double time_preview_reference_embedding = 0.0;
    std::size_t preview_attempts = 0;
    int conformity_collar_layers = -1;
    bool closure_safe = false;
    bool full_regular_doerfler = false;
};

ClosureAwareHybridMarking mark_hybrid_regular_region_closure_aware(
    const ReferenceEpochHierarchy &hierarchy,
    const std::vector<double> &element_eta_squared,
    const SingularRegionClassification &regions,
    const double theta,
    const ReferenceEpochRefinementGuard &resource_guard) {
    if (element_eta_squared.size() != hierarchy.coarse_mesh().elems.size()
        || element_eta_squared.size() != regions.in_regular.size()) {
        throw std::invalid_argument(
            "closure-aware hybrid indicators do not match the coarse mesh");
    }
    if (!(theta > 0.0 && theta <= 1.0)) {
        throw std::invalid_argument(
            "closure-aware hybrid Doerfler theta must lie in (0,1]");
    }

    ClosureAwareHybridMarking result;
    for (int element = 0;
         element < static_cast<int>(element_eta_squared.size()); ++element) {
        const double value = element_eta_squared[element];
        if (!(value >= 0.0) || !std::isfinite(value)) {
            throw std::invalid_argument(
                "closure-aware hybrid indicators must be finite and nonnegative");
        }
        if (regions.in_regular[element]) result.regular_mass += value;
    }
    if (!(result.regular_mass > 0.0)) {
        result.closure_safe = true;
        result.full_regular_doerfler = true;
        result.conformity_collar_layers = 0;
        return result;
    }

    const double required_mass = theta * result.regular_mass;
    const HybridGradedReserveProfile distances =
        make_hybrid_graded_reserve_profile(
            hierarchy.coarse_mesh(), regions, 1, 0);
    std::vector<int> fallback = mark_hybrid_regular_region(
        element_eta_squared, regions, theta);
    double fallback_mass = 0.0;
    for (const int element : fallback)
        fallback_mass += element_eta_squared[element];

    for (int collar = 0; collar <= distances.maximum_graph_distance;
         ++collar) {
        const HybridGradedReserveProfile profile =
            make_hybrid_graded_reserve_profile(
                hierarchy.coarse_mesh(), regions, 1, collar);
        std::vector<char> eligible(element_eta_squared.size(), false);
        double admissible_mass = 0.0;
        for (int element = 0;
             element < static_cast<int>(element_eta_squared.size()); ++element) {
            eligible[element] = regions.in_regular[element]
                && profile.target_level_gaps[element] > 0;
            if (eligible[element])
                admissible_mass += element_eta_squared[element];
        }
        if (admissible_mass + 1e-14 * result.regular_mass < required_mass)
            break;

        const double active_theta = std::min(
            1.0, required_mass / admissible_mass);
        std::vector<double> admissible_indicators = element_eta_squared;
        for (int element = 0;
             element < static_cast<int>(admissible_indicators.size()); ++element) {
            if (!eligible[element]) admissible_indicators[element] = 0.0;
        }
        std::vector<int> marked = mark_doerfler(
            admissible_indicators, active_theta, eligible);
        double marked_mass = 0.0;
        for (const int element : marked)
            marked_mass += element_eta_squared[element];
        const bool full_doerfler =
            marked_mass + 1e-14 * result.regular_mass >= required_mass;
        if (!full_doerfler) continue;
        const Clock::time_point preview_start = Clock::now();
        ReferenceEpochCoarseRefinementPreview preview =
            hierarchy.preview_coarse_refinement(marked, resource_guard);
        result.time_preview += seconds_since(preview_start);
        result.time_preview_nvb += preview.time_nvb_refine;
        result.time_preview_reference_embedding +=
            preview.time_reference_embedding_update;
        ++result.preview_attempts;
        if (!preview.reference_contained())
            continue;

        result.marked_elements = std::move(marked);
        result.accepted_preview = std::move(preview);
        result.admissible_mass = admissible_mass;
        result.marked_mass = marked_mass;
        result.conformity_collar_layers = collar;
        result.closure_safe = true;
        result.full_regular_doerfler = true;
        return result;
    }

    // The paper's full regular-region bulk condition takes precedence.  If
    // no closure-safe subset carries enough mass, retain the original set and
    // let the structural branch refresh the reference as prescribed.
    result.marked_elements = std::move(fallback);
    result.admissible_mass = std::numeric_limits<double>::quiet_NaN();
    result.marked_mass = fallback_mass;
    result.full_regular_doerfler =
        fallback_mass + 1e-14 * result.regular_mass >= required_mass;
    return result;
}

std::size_t embedding_children(
    const Eigen::SparseMatrix<double> &embedding,
    const int parent_element) {
    std::size_t count = 0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(
             embedding, parent_element); it; ++it) {
        if (std::abs(it.value()) > 1e-14) ++count;
    }
    return count;
}

enum class HybridMatchingTarget { Reference, Candidate };

struct ProposedHybridMatchingResult {
    SingularRegionClassification regions;
    std::size_t refinement_rounds = 0;
    std::size_t added_elements = 0;
    double time_classification = 0.0;
    double time_nvb = 0.0;
    double time_embedding_update = 0.0;
};

ProposedHybridMatchingResult restore_proposed_hybrid_matching(
    ReferenceEpochHierarchy &hierarchy,
    const int ell,
    const double physical_radius,
    const HybridMatchingTarget target,
    const std::size_t maximum_rounds = 64,
    const ReferenceEpochRefinementGuard &resource_guard = {}) {
    if (!hierarchy.has_proposed_coarse_refinement() || maximum_rounds == 0) {
        throw std::invalid_argument(
            "prospective hybrid matching requires a pending proposal and a positive limit");
    }
    ProposedHybridMatchingResult result;
    for (;;) {
        const bool contained = target == HybridMatchingTarget::Reference
            ? hierarchy.reference_contains_proposed_coarse()
            : hierarchy.candidate_contains_proposed_coarse();
        if (!contained) {
            throw std::runtime_error(
                "prospective hybrid matching target does not contain the proposal");
        }
        const Clock::time_point classification_start = Clock::now();
        result.regions = classify_singular_regions_with_physical_radius(
            hierarchy.proposed_coarse_mesh(), {Point2(0.0, 0.0)}, ell,
            physical_radius);
        result.time_classification += seconds_since(classification_start);
        const Eigen::SparseMatrix<double> &embedding =
            target == HybridMatchingTarget::Reference
            ? hierarchy.proposed_coarse_elements_to_reference()
            : hierarchy.proposed_coarse_elements_to_candidate();
        std::vector<int> marked;
        for (const int element : result.regions.omega_f_elements) {
            if (embedding_children(embedding, element) > 1)
                marked.push_back(element);
        }
        if (marked.empty()) return result;
        if (result.refinement_rounds >= maximum_rounds) {
            throw ReferenceEpochWorkLimitExceeded(
                "prospective hybrid matching refinement limit reached");
        }
        const std::size_t before =
            hierarchy.proposed_coarse_mesh().elems.size();
        const ReferenceEpochRefinementResult refined =
            hierarchy.refine_proposed_coarse(marked, resource_guard);
        if (!refined.changed()) {
            throw std::runtime_error(
                "prospective hybrid matching did not refine the proposal");
        }
        result.added_elements +=
            hierarchy.proposed_coarse_mesh().elems.size() - before;
        result.time_nvb += refined.time_nvb_refine;
        result.time_embedding_update +=
            refined.time_parent_map_update
            + refined.time_embedding_composition
            + refined.time_proposed_embedding_update;
        ++result.refinement_rounds;
    }
}

HybridMatchingResult restore_committed_hybrid_matching_incrementally(
    ReferenceEpochHierarchy &hierarchy,
    const int ell,
    const double physical_radius,
    const std::size_t maximum_rounds = 64,
    const ReferenceEpochRefinementGuard &resource_guard = {}) {
    if (hierarchy.has_proposed_coarse_refinement()) {
        throw std::invalid_argument(
            "committed hybrid matching cannot start with a pending proposal");
    }
    const std::uint64_t reference_version = hierarchy.reference_mesh_version();
    const std::size_t reference_elements = hierarchy.reference_mesh().elems.size();
    HybridMatchingResult result;
    result.regions = classify_singular_regions_with_physical_radius(
        hierarchy.coarse_mesh(), {Point2(0.0, 0.0)}, ell,
        physical_radius);
    std::vector<int> marked;
    for (const int element : result.regions.omega_f_elements) {
        if (embedding_children(
                hierarchy.coarse_elements_to_reference(), element) > 1) {
            marked.push_back(element);
        }
    }
    if (!marked.empty()) {
        const std::size_t before = hierarchy.coarse_mesh().elems.size();
        const ReferenceEpochRefinementResult first =
            hierarchy.propose_coarse_refinement(marked, resource_guard);
        if (!first.changed() || !hierarchy.reference_contains_proposed_coarse()) {
            throw std::runtime_error(
                "hybrid matching proposal is not contained in the fixed reference");
        }
        ProposedHybridMatchingResult closure =
            restore_proposed_hybrid_matching(
                hierarchy, ell, physical_radius,
                HybridMatchingTarget::Reference, maximum_rounds,
                resource_guard);
        const ReferenceEpochRefinementResult committed =
            hierarchy.commit_coarse_refinement();
        if (!committed.changed()) {
            throw std::runtime_error("hybrid matching coarse commit failed");
        }
        result.refinement_rounds = 1 + closure.refinement_rounds;
        result.refined_coarse_elements =
            hierarchy.coarse_mesh().elems.size() - before;
    }
    result.regions = classify_singular_regions_with_physical_radius(
        hierarchy.coarse_mesh(), {Point2(0.0, 0.0)}, ell,
        physical_radius);
    result.matching_holds = hybrid_reference_matching_holds(
        hierarchy, result.regions);
    result.reference_unchanged =
        hierarchy.reference_mesh_version() == reference_version
        && hierarchy.reference_mesh().elems.size() == reference_elements;
    if (!result.matching_holds || !result.reference_unchanged) {
        throw std::runtime_error(
            "incremental hybrid matching invariant was not restored");
    }
    return result;
}

struct HybridCorrectorPatchAudit {
    std::size_t active_patches = 0;
    std::size_t maximum_fine_elements = 0;
    std::size_t p95_fine_elements = 0;
    int maximum_target = -1;
};

HybridCorrectorPatchAudit audit_hybrid_corrector_patch_costs(
    const ReferenceEpochHierarchy &hierarchy, const int ell,
    const SingularRegionClassification &regions) {
    const std::size_t coarse_elements = hierarchy.coarse_mesh().elems.size();
    const Eigen::SparseMatrix<double> &embedding =
        hierarchy.coarse_elements_to_reference();
    if (embedding.cols() != static_cast<int>(coarse_elements)
        || regions.in_omega_s.size() != coarse_elements) {
        throw std::logic_error(
            "hybrid corrector preflight metadata has inconsistent dimensions");
    }

    std::vector<std::size_t> reference_children(coarse_elements, 0);
    for (int coarse = 0; coarse < embedding.outerSize(); ++coarse) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(embedding, coarse);
             it; ++it) {
            if (it.value() != 0.0) ++reference_children[coarse];
        }
    }

    const Eigen::SparseMatrix<double> patches =
        build_patches(hierarchy.coarse_mesh(), ell);
    std::vector<std::size_t> costs;
    costs.reserve(coarse_elements);
    HybridCorrectorPatchAudit audit;
    for (int target = 0; target < patches.outerSize(); ++target) {
        if (regions.in_omega_s[target]) continue;
        std::size_t cost = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(patches, target);
             it; ++it) {
            if (it.value() != 0.0) cost += reference_children[it.row()];
        }
        costs.push_back(cost);
        if (cost > audit.maximum_fine_elements) {
            audit.maximum_fine_elements = cost;
            audit.maximum_target = target;
        }
    }
    audit.active_patches = costs.size();
    if (!costs.empty()) {
        std::sort(costs.begin(), costs.end());
        const std::size_t p95_index =
            (95 * costs.size() + 99) / 100 - 1;
        audit.p95_fine_elements = costs[p95_index];
    }
    return audit;
}

class NumericalReferenceEpochBackend final : public ReferenceEpochDriverBackend {
public:
    NumericalReferenceEpochBackend(
        ReferenceEpochPaperConfig config, PaperCaseData data)
        : config_(std::move(config)), data_(std::move(data)),
          hierarchy_(data_.initial_mesh, config_.initial_coarse_level,
                     config_.initial_reference_level) {}

    void begin_epoch() override {
        if (epoch_started_) ++epoch_index_;
        epoch_started_ = true;
        hierarchy_.begin_reference_epoch();
        reference_solution_.reset();
        reference_load_.reset();
        reference_discrete_norm_.reset();
        reference_exact_relative_energy_.reset();
        model_.reset();
        solution_.reset();
        regions_.reset();
        proposed_regions_.reset();
        accepted_coarse_preview_.reset();
        epoch_start_snapshot_pending_ = true;
    }

    ReferenceEpochCorrectorObservation corrector_check(int ell) override {
        current_ell_ = ell;
        const Clock::time_point total_start = Clock::now();
        std::vector<int> skipped;
        if (config_.singularity_hybrid) {
            std::cerr
                << "[hybrid-matching] begin ell=" << ell
                << " R=" << config_.hybrid_minimum_physical_radius
                << " coarse_elements=" << hierarchy_.coarse_mesh().elems.size()
                << " reference_elements="
                << hierarchy_.reference_mesh().elems.size() << std::endl;
            const Clock::time_point matching_start = Clock::now();
            HybridMatchingResult matching =
                restore_committed_hybrid_matching_incrementally(
                    hierarchy_, ell,
                    config_.hybrid_minimum_physical_radius, 64,
                    proposed_refinement_resource_guard());
            regions_ = matching.regions;
            skipped = matching.regions.omega_s_elements;
            std::cerr
                << "[hybrid-matching] end seconds="
                << seconds_since(matching_start)
                << " ell_S=" << matching.regions.l_s
                << " covered_radius=" << matching.regions.covered_physical_radius
                << " refinement_rounds=" << matching.refinement_rounds
                << " coarse_elements=" << hierarchy_.coarse_mesh().elems.size()
                << " omega_s=" << matching.regions.omega_s_elements.size()
                << " omega_f=" << matching.regions.omega_f_elements.size()
                << std::endl;
            const HybridCorrectorPatchAudit patch_audit =
                audit_hybrid_corrector_patch_costs(
                    hierarchy_, ell, matching.regions);
            std::cerr
                << "[hybrid-preflight] ell=" << ell
                << " ell_S=" << matching.regions.l_s
                << " active_patches=" << patch_audit.active_patches
                << " max_patch_fine_elements="
                << patch_audit.maximum_fine_elements
                << " p95_patch_fine_elements="
                << patch_audit.p95_fine_elements
                << " max_target=" << patch_audit.maximum_target
                << " guard="
                << config_.hybrid_maximum_corrector_patch_fine_elements
                << std::endl;
            if (patch_audit.maximum_fine_elements
                > config_.hybrid_maximum_corrector_patch_fine_elements) {
                std::ostringstream message;
                message
                    << "hybrid corrector preflight rejected pathological patch: "
                    << "target=" << patch_audit.maximum_target
                    << ", fine_elements="
                    << patch_audit.maximum_fine_elements
                    << ", limit="
                    << config_.hybrid_maximum_corrector_patch_fine_elements
                    << ", ell=" << ell
                    << ", ell_S=" << matching.regions.l_s
                    << ", R=" << config_.hybrid_minimum_physical_radius;
                throw ReferenceEpochWorkLimitExceeded(message.str());
            }
        }
        HelmholtzProblemConfig model_config;
        model_config.ell = ell;
        model_config.wavenumber = config_.wavenumber;
        model_config.initial_mesh = data_.initial_mesh;
        model_config.quadrature = config_.quadrature;
        model_config.quadrature_context = data_.quadrature_context;
        model_config.progress = config_.singularity_hybrid;
        if (config_.singularity_hybrid) {
            model_ = std::make_unique<HelmholtzLodModel>(
                HelmholtzLodModel::build_adaptive_hybrid(
                    model_config, hierarchy_, skipped, &corrector_cache_));
        } else {
            model_ = std::make_unique<HelmholtzLodModel>(
                HelmholtzLodModel::build_adaptive(
                    model_config, hierarchy_, &corrector_cache_));
        }
        double time_reference_stability = 0.0;
        if (!reference_solution_) {
            const Clock::time_point stability_start = Clock::now();
            reference_load_ = assemble_helmholtz_load(
                hierarchy_.reference_mesh(), data_.source,
                config_.quadrature, data_.quadrature_context);
            ensure_reference_solution(*reference_load_);
            time_reference_stability = seconds_since(stability_start);
            std::cerr
                << "[reference-stability] epoch=" << epoch_index_
                << " reference_unknowns="
                << unknowns(hierarchy_.reference_mesh())
                << " seconds=" << time_reference_stability
                << " status=solved" << std::endl;
        }
        const Clock::time_point theta_start = Clock::now();
        LocalizationEigenConfig eigen;
        eigen.maximum_iterations =
            config_.localization_eigen_maximum_iterations;
        eigen.relative_tolerance =
            config_.localization_eigen_relative_tolerance;
        // Small clustered spectra are both cheaper and more reliable through
        // explicit block Gram assembly. The matrix-free three-vector Ritz
        // path remains mandatory above this modest coarse dimension.
        eigen.dense_cross_check_max_dimension = 128;
        // Production runs must never silently replace a failed matrix-free
        // iteration by O(N_H) full Gram actions.  E0 exercises the explicit
        // dense path separately on its fixed small hierarchy.
        eigen.dense_fallback_max_dimension = 0;
        eigen.warm_start = localization_warm_start_;
        eigen.warm_start_block = localization_warm_start_block_;
        const ReferenceCorrectorCertificate certificate =
            build_reference_corrector_certificate(
                hierarchy_, model_->operators(), model_->corrected_test_basis(),
                free_nodes(hierarchy_.coarse_mesh()),
                KernelRieszSolver::SaddlePoint, eigen, &gram_factor_cache_);
        localization_warm_start_ = certificate.spectrum.dominant_vector;
        localization_warm_start_block_ =
            certificate.spectrum.dominant_subspace;
        ReferenceEpochCorrectorObservation result;
        result.theta_loc = certificate.theta_loc;
        // Compare the frozen implementation-study threshold against the
        // residual-enclosed upper endpoint, not merely the central Ritz value.
        result.delta_loc_hat = certificate.theta_loc_diagnostic_upper;
        const HelmholtzCorrectorDiagnostics &corrector_diagnostics =
            model_->correctors().diagnostics;
        result.active_correctors = static_cast<std::size_t>(
            corrector_diagnostics.patch_count);
        result.rebuilt_correctors = static_cast<std::size_t>(
            corrector_diagnostics.patch_cache_misses);
        result.reused_correctors = static_cast<std::size_t>(
            corrector_diagnostics.patch_cache_hits);
        result.corrector_cache_oversized_misses = static_cast<std::size_t>(
            corrector_diagnostics.patch_cache_oversized_misses);
        result.corrector_cache_budget_rejections = static_cast<std::size_t>(
            corrector_diagnostics.patch_cache_budget_rejections);
        const HelmholtzCorrectorPatchCache::Statistics cache_statistics =
            corrector_cache_.statistics();
        result.corrector_cache_entries = cache_statistics.entries;
        result.corrector_cache_current_bytes = cache_statistics.current_bytes;
        result.corrector_cache_peak_bytes = cache_statistics.peak_bytes;
        if (result.active_correctors
            != result.rebuilt_correctors + result.reused_correctors) {
            throw std::runtime_error(
                "corrector cache accounting does not recover active patches");
        }
        result.skipped_correctors = static_cast<std::size_t>(
            corrector_diagnostics.skipped_patch_count);
        result.skipped_corrector_work_units =
            corrector_diagnostics.skipped_patch_work_units;
        result.corrector_parallel_threads =
            corrector_diagnostics.parallel_threads;
        result.corrector_patch_assembly_work_seconds =
            corrector_diagnostics.patch_assembly_work_seconds;
        result.corrector_patch_solve_work_seconds =
            corrector_diagnostics.patch_solve_work_seconds;
        result.corrector_patch_pack_work_seconds =
            corrector_diagnostics.patch_pack_work_seconds;
        result.corrector_maximum_patch_dofs =
            corrector_diagnostics.maximum_patch_dofs;
        result.corrector_maximum_patch_constraints =
            corrector_diagnostics.maximum_patch_constraints;
        result.corrector_maximum_patch_rhs =
            corrector_diagnostics.maximum_patch_rhs;
        if (regions_) {
            result.hybrid_l_s = regions_->l_s;
            result.hybrid_minimum_physical_radius =
                regions_->minimum_physical_radius;
            result.hybrid_covered_physical_radius =
                regions_->covered_physical_radius;
            result.hybrid_omega_s_elements =
                regions_->omega_s_elements.size();
        result.hybrid_omega_f_elements =
                regions_->omega_f_elements.size();
        }
        result.time_reference_stability = time_reference_stability;
        const HelmholtzBuildTimings &build = model_->build_timings();
        result.time_lod_build_total = build.total_ms / 1000.0;
        result.time_lod_mesh_and_interpolation =
            build.mesh_and_interpolation_ms / 1000.0;
        result.time_lod_operators = build.operators_ms / 1000.0;
        result.time_corrector = build.correctors_ms / 1000.0;
        result.time_lod_corrected_basis = build.corrected_basis_ms / 1000.0;
        result.time_lod_coarse_operator = build.coarse_operator_ms / 1000.0;
        result.time_lod_coarse_factorization =
            build.coarse_factorization_ms / 1000.0;
        result.time_theta = seconds_since(theta_start);
        result.time_gram_prepare_structure =
            certificate.gram_operator.prepare_structure_seconds;
        result.time_gram_prepare_factorization =
            certificate.gram_operator.prepare_factorization_seconds;
        result.time_gram_action_rhs =
            certificate.gram_operator.action_rhs_seconds;
        result.time_gram_action_patch_solve =
            certificate.gram_operator.action_patch_solve_seconds;
        result.time_gram_action_scatter =
            certificate.gram_operator.action_scatter_seconds;
        result.gram_action_calls = certificate.gram_operator.action_calls;
        result.gram_patch_factorizations =
            certificate.gram_operator.patch_factorizations;
        result.gram_factor_cache_hits =
            certificate.gram_operator.factor_cache_hits;
        result.gram_factor_cache_misses =
            certificate.gram_operator.factor_cache_misses;
        result.gram_structure_parallel_threads =
            certificate.gram_operator.structure_parallel_threads;
        result.gram_parallel_threads =
            certificate.gram_operator.parallel_threads;
        result.localization_iterations = certificate.spectrum.iterations;
        result.localization_relative_residual =
            certificate.spectrum.relative_residual;
        result.localization_used_warm_start =
            certificate.spectrum.used_warm_start;
        if (epoch_start_snapshot_pending_
            && result.delta_loc_hat <= config_.theta_loc_usr) {
            const Clock::time_point artifact_start = Clock::now();
            ReferenceEpochMeshSnapshot snapshot;
            snapshot.epoch = epoch_index_;
            snapshot.H_step = committed_H_steps_;
            snapshot.reference_mesh_version =
                hierarchy_.reference_mesh_version();
            snapshot.stage = "epoch_start";
            snapshot.anchor_action = ReferenceEpochDriverAction::BeginEpoch;
            snapshot.anchor_occurrence = begin_occurrence_++;
            snapshot.coarse = hierarchy_.coarse_mesh();
            snapshot.candidate = hierarchy_.candidate_mesh();
            snapshot.reference = hierarchy_.reference_mesh();
            snapshots_.push_back(std::move(snapshot));
            artifact_capture_seconds_ += seconds_since(artifact_start);
            epoch_start_snapshot_pending_ = false;
        }
        result.time_corrector_check_total = seconds_since(total_start);
        return result;
    }

    ReferenceEpochSolveObservation solve_and_estimate() override {
        if (!model_) throw std::logic_error("solve requested before corrector check");
        if (!reference_load_)
            throw std::logic_error(
                "reference load was not prepared by the epoch stability check");
        const ComplexVector &load = *reference_load_;
        const Clock::time_point solve_start = Clock::now();
        solution_ = model_->solve_load(load);
        ReferenceEpochSolveObservation result;
        result.time_lod_solve = seconds_since(solve_start);
        const Clock::time_point riesz_start = Clock::now();
        const ReferenceResidualRiesz estimate = compute_reference_residual_riesz(
            hierarchy_, model_->operators(), load, solution_->fine_values,
            config_.theta_H, KernelRieszSolver::SaddlePoint);
        result.time_reference_riesz = seconds_since(riesz_start);
        result.eta_H = estimate.eta;
        result.U_practical = config_.C_rel_usr * estimate.eta;
        if (config_.singularity_hybrid) {
            const Clock::time_point marking_start = Clock::now();
            ClosureAwareHybridMarking marking =
                mark_hybrid_regular_region_closure_aware(
                    hierarchy_, estimate.element_eta_squared, *regions_,
                    config_.theta_H,
                    proposed_refinement_resource_guard());
            result.marked_H = marking.marked_elements;
            result.hybrid_regular_indicator_mass = marking.regular_mass;
            result.hybrid_admissible_indicator_mass = marking.admissible_mass;
            result.hybrid_marked_H_indicator_mass = marking.marked_mass;
            result.hybrid_coarse_conformity_collar =
                marking.conformity_collar_layers;
            result.hybrid_coarse_marking_closure_safe = marking.closure_safe;
            result.hybrid_full_regular_doerfler =
                marking.full_regular_doerfler;
            result.hybrid_coarse_preview_attempts = marking.preview_attempts;
            result.hybrid_coarse_preview_cached =
                marking.accepted_preview.has_value();
            result.time_hybrid_coarse_preview = marking.time_preview;
            result.time_hybrid_coarse_preview_nvb = marking.time_preview_nvb;
            result.time_hybrid_coarse_preview_reference_embedding =
                marking.time_preview_reference_embedding;
            accepted_coarse_preview_ = std::move(marking.accepted_preview);
            result.time_hybrid_coarse_marking = seconds_since(marking_start);
            std::cerr
                << "[hybrid-coarse-marking] epoch=" << epoch_index_
                << " H_step=" << committed_H_steps_
                << " regular_mass=" << marking.regular_mass
                << " admissible_mass=" << marking.admissible_mass
                << " marked_mass=" << marking.marked_mass
                << " collar=" << marking.conformity_collar_layers
                << " closure_safe="
                << (marking.closure_safe ? "true" : "false")
                << " full_doerfler="
                << (marking.full_regular_doerfler ? "true" : "false")
                << " preview_attempts=" << marking.preview_attempts
                << " preview_cached="
                << (result.hybrid_coarse_preview_cached ? "true" : "false")
                << " time_marking=" << result.time_hybrid_coarse_marking
                << " time_preview=" << result.time_hybrid_coarse_preview
                << std::endl;
        } else {
            accepted_coarse_preview_.reset();
            result.marked_H = estimate.marked_elements;
        }

        const Clock::time_point reference_validation_start = Clock::now();
        ensure_reference_solution(load);
        if (!reference_discrete_norm_) {
            reference_discrete_norm_ = compute_discrete_helmholtz_error(
                hierarchy_.reference_mesh(), model_->operators(),
                *reference_solution_,
                ComplexVector::Zero(reference_solution_->size()));
        }
        const HelmholtzError reference_error = compute_discrete_helmholtz_error(
            hierarchy_.reference_mesh(), model_->operators(),
            *reference_solution_, solution_->fine_values);
        result.relative_reference_energy = reference_error.energy
            / std::max(
                reference_discrete_norm_->energy,
                std::numeric_limits<double>::min());
        result.time_reference_validation =
            seconds_since(reference_validation_start);
        const Clock::time_point exact_validation_start = Clock::now();
        if (data_.exact && data_.exact_gradient) {
            const HelmholtzError exact_error = compute_helmholtz_error(
                hierarchy_.reference_mesh(), solution_->fine_values,
                config_.wavenumber, data_.exact, data_.exact_gradient,
                config_.quadrature, data_.quadrature_context);
            if (!reference_exact_norm_) {
                reference_exact_norm_ = compute_helmholtz_error(
                    hierarchy_.reference_mesh(),
                    ComplexVector::Zero(solution_->fine_values.size()),
                    config_.wavenumber, data_.exact, data_.exact_gradient,
                    config_.quadrature, data_.quadrature_context);
            }
            result.relative_exact_energy =
                exact_error.energy / reference_exact_norm_->energy;
            result.relative_exact_L2 =
                exact_error.l2 / reference_exact_norm_->l2;
            if (!reference_exact_relative_energy_) {
                const HelmholtzError exact_reference_error =
                    compute_helmholtz_error(
                        hierarchy_.reference_mesh(), *reference_solution_,
                        config_.wavenumber, data_.exact, data_.exact_gradient,
                        config_.quadrature, data_.quadrature_context);
                reference_exact_relative_energy_ =
                    exact_reference_error.energy / reference_exact_norm_->energy;
            }
            result.relative_exact_reference_energy =
                *reference_exact_relative_energy_;
        }
        result.time_exact_validation = seconds_since(exact_validation_start);
        {
            const Clock::time_point artifact_start = Clock::now();
            latest_solved_coarse_ = hierarchy_.coarse_mesh();
            artifact_capture_seconds_ += seconds_since(artifact_start);
        }
        return result;
    }

    void propose_coarse_refinement(const std::vector<int> &marked_H) override {
        proposed_regions_.reset();
        const ReferenceEpochRefinementGuard proposal_guard =
            proposed_refinement_resource_guard();
        ReferenceEpochRefinementResult proposed;
        if (marked_H.empty()) {
            accepted_coarse_preview_.reset();
            proposed = hierarchy_.propose_identity_coarse();
        } else if (accepted_coarse_preview_) {
            if (accepted_coarse_preview_->marked_elements != marked_H) {
                throw std::logic_error(
                    "accepted hybrid coarse preview does not match driver marking");
            }
            proposed = hierarchy_.propose_coarse_refinement(
                std::move(*accepted_coarse_preview_), proposal_guard);
            accepted_coarse_preview_.reset();
        } else {
            proposed = hierarchy_.propose_coarse_refinement(
                marked_H, proposal_guard);
        }
        if (!proposed.changed())
            throw std::runtime_error("coarse proposal transaction did not open");
        if (!config_.singularity_hybrid) return;

        const Clock::time_point matching_start = Clock::now();
        ProposedHybridMatchingResult matching;
        if (hierarchy_.reference_contains_proposed_coarse()) {
            matching = restore_proposed_hybrid_matching(
                hierarchy_, current_ell_,
                config_.hybrid_minimum_physical_radius,
                HybridMatchingTarget::Reference, 64, proposal_guard);
        } else {
            matching.regions =
                classify_singular_regions_with_physical_radius(
                    hierarchy_.proposed_coarse_mesh(),
                    {Point2(0.0, 0.0)}, current_ell_,
                    config_.hybrid_minimum_physical_radius);
        }
        proposed_regions_ = matching.regions;
        std::cerr
            << "[hybrid-prospective-matching] seconds="
            << seconds_since(matching_start)
            << " ell_S=" << matching.regions.l_s
            << " rounds=" << matching.refinement_rounds
            << " added_elements=" << matching.added_elements
            << " time_classification=" << matching.time_classification
            << " time_nvb=" << matching.time_nvb
            << " time_embedding_update=" << matching.time_embedding_update
            << " contained_in_reference="
            << (hierarchy_.reference_contains_proposed_coarse()
                    ? "true" : "false")
            << std::endl;
    }

    ReferenceEpochCandidateObservation enrich_candidate() override {
        if (!solution_) throw std::logic_error("candidate enrichment has no LOD value");
        const Clock::time_point total_start = Clock::now();
        const ReferenceEpochRefinementGuard candidate_guard =
            candidate_refinement_resource_guard();
        const Clock::time_point operator_start = Clock::now();
        HelmholtzOperators operators = assemble_helmholtz_operators(
            hierarchy_.candidate_mesh(), config_.wavenumber);
        const double time_operator = seconds_since(operator_start);
        const Clock::time_point prolongation_start = Clock::now();
        ComplexVector on_candidate = hierarchy_.reference_to_candidate().cast<Complex>()
            * solution_->fine_values;
        const double time_prolongation = seconds_since(prolongation_start);
        CandidateFluxConfig flux_config;
        flux_config.doerfler_theta = config_.theta_c;
        flux_config.quadrature = config_.quadrature;
        flux_config.quadrature_context = data_.quadrature_context;
        const Clock::time_point flux_start = Clock::now();
        const CandidateFluxRT2Result flux = reconstruct_candidate_flux_rt2(
            hierarchy_.candidate_mesh(), operators, data_.source,
            on_candidate, flux_config);
        const double time_flux = seconds_since(flux_start);
        ReferenceEpochCandidateObservation result;
        result.eta_eq_c = flux.eta_eq;
        std::vector<int> marked_candidate = flux.marked_elements;
        if (config_.singularity_hybrid) {
            if (!proposed_regions_)
                throw std::logic_error(
                    "hybrid candidate marking has no prospective regions");
            // Algorithm 2 constructs and marks the indicator on the current
            // candidate mesh, before hierarchy closure.  When the prospective
            // coarse mesh is already nested, classify through its cached
            // parent map.  A non-nested proposal is a structural-refresh
            // case, for which the paper retains the current physical regions.
            const bool classify_on_proposed =
                hierarchy_.candidate_contains_proposed_coarse();
            const std::vector<int> parents = classify_on_proposed
                ? hierarchy_.candidate_parent_proposed_coarse_elements()
                : hierarchy_.candidate_parent_coarse_elements();
            const SingularRegionClassification &marking_regions =
                classify_on_proposed ? *proposed_regions_ : *regions_;
            std::vector<char> in_omega_f(parents.size(), false);
            for (int element = 0;
                 element < static_cast<int>(parents.size()); ++element) {
                const int parent = parents[element];
                if (parent < 0
                    || parent >= static_cast<int>(
                        marking_regions.in_omega_f.size())) {
                    throw std::logic_error(
                        "candidate/prospective hybrid parent is invalid");
                }
                in_omega_f[element] =
                    marking_regions.in_omega_f[parent];
                if (in_omega_f[element]) ++result.candidate_cells_f;
                else ++result.candidate_cells_r;
            }
            const SplitRegionalDoerflerMarking regional =
                mark_split_regional_doerfler(
                    flux.element_eta_squared, in_omega_f, config_.theta_c);
            marked_candidate = regional.marked_elements;
            result.eta_eq_c_f = std::sqrt(regional.omega_f_mass);
            result.eta_eq_c_r = std::sqrt(regional.regular_mass);
            result.indicator_mass_c_f = regional.omega_f_mass;
            result.indicator_mass_c_r = regional.regular_mass;
            result.marked_mass_c_f = regional.marked_omega_f_mass;
            result.marked_mass_c_r = regional.marked_regular_mass;
            result.marked_c_f =
                regional.marked_omega_f_elements.size();
            result.marked_c_r =
                regional.marked_regular_elements.size();
            const double total_mass = regional.omega_f_mass
                + regional.regular_mass;
            const double flux_mass = flux.eta_eq * flux.eta_eq;
            if (std::abs(total_mass - flux_mass)
                > 1e-9 * std::max(1.0, flux_mass)) {
                throw std::runtime_error(
                    "hybrid regional candidate masses do not recover eta_eq_c");
            }
        }
        result.marked_c = marked_candidate;
        const Clock::time_point enrich_start = Clock::now();
        ReferenceEpochRefinementResult enrich_result;
        if (!marked_candidate.empty())
            enrich_result = hierarchy_.enrich_candidate(
                marked_candidate, candidate_guard);
        result.time_candidate_enrich = seconds_since(enrich_start);
        const Clock::time_point close_start = Clock::now();
        const ReferenceEpochRefinementResult close_result =
            hierarchy_.close_candidate_over_proposed_coarse(candidate_guard);
        const double time_close = seconds_since(close_start);
        result.time_candidate_close = time_close;
        result.time_candidate_operator_assembly = time_operator;
        result.time_candidate_prolongation = time_prolongation;
        result.time_candidate_flux_reconstruction = time_flux;
        result.time_candidate_flux_prepare = flux.time_prepare;
        result.time_candidate_flux_patch_solve = flux.time_patch_solve;
        result.time_candidate_flux_merge = flux.time_deterministic_merge;
        result.time_candidate_flux_audit = flux.time_estimator_and_audit;
        result.candidate_flux_parallel_threads = flux.parallel_threads;
        result.time_candidate_nvb_refine = close_result.time_nvb_refine
            + enrich_result.time_nvb_refine;
        result.time_candidate_embedding_composition =
            close_result.time_embedding_composition
            + enrich_result.time_embedding_composition;
        result.time_candidate_parent_map_update =
            close_result.time_parent_map_update
            + enrich_result.time_parent_map_update;
        result.time_candidate_quasi_interpolation =
            close_result.time_candidate_quasi_interpolation
            + enrich_result.time_candidate_quasi_interpolation;
        result.time_candidate_embedding_validation =
            close_result.time_embedding_validation
            + enrich_result.time_embedding_validation;
        result.time_candidate_flux = seconds_since(total_start);
        return result;
    }

    bool proposal_contained_in_reference() const override {
        return hierarchy_.reference_contains_proposed_coarse();
    }

    int minimum_reference_level_gap() const override {
        if (config_.singularity_hybrid && regions_) {
            if (hierarchy_.has_proposed_coarse_refinement()) {
                if (!proposed_regions_)
                    throw std::logic_error(
                        "hybrid level reserve has no prospective regions");
                const HybridGradedReserveProfile reserve =
                    make_hybrid_graded_reserve_profile(
                        hierarchy_.proposed_coarse_mesh(),
                        *proposed_regions_,
                        config_.reference_refresh_target_gap,
                        hybrid_conformity_collar_layers_);
                if (config_.reference_refresh_target_gap > 0
                    && reserve.full_target_elements == 0) {
                    return 0;
                }
                return prospective_reference_level_gap(
                    hierarchy_, &reserve.at_full_target);
            }
            const HybridGradedReserveProfile reserve =
                make_hybrid_graded_reserve_profile(
                    hierarchy_.coarse_mesh(), *regions_,
                    config_.reference_refresh_target_gap,
                    hybrid_conformity_collar_layers_);
            if (config_.reference_refresh_target_gap > 0
                && reserve.full_target_elements == 0) {
                return 0;
            }
            return prospective_reference_level_gap(
                hierarchy_, &reserve.at_full_target);
        }
        return prospective_reference_level_gap(hierarchy_);
    }

    ReferenceEpochDualObservation candidate_dual_check(
        double reference_upper_bound) override {
        if (!solution_) throw std::logic_error("candidate dual has no LOD value");
        const Clock::time_point total_start = Clock::now();
        const Clock::time_point operator_start = Clock::now();
        const HelmholtzOperators operators = assemble_helmholtz_operators(
            hierarchy_.candidate_mesh(), config_.wavenumber);
        const double time_operator = seconds_since(operator_start);
        const Clock::time_point wellposedness_start = Clock::now();
        verify_discrete_helmholtz_factorization(operators);
        const double time_wellposedness = seconds_since(wellposedness_start);
        const Clock::time_point load_start = Clock::now();
        const ComplexVector load = assemble_helmholtz_load(
            hierarchy_.candidate_mesh(), data_.source, config_.quadrature,
            data_.quadrature_context);
        const double time_load = seconds_since(load_start);
        const Clock::time_point prolongation_start = Clock::now();
        const ComplexVector on_candidate =
            hierarchy_.reference_to_candidate().cast<Complex>()
            * solution_->fine_values;
        const double time_prolongation = seconds_since(prolongation_start);
        CandidateDualGapConfig gap_config;
        gap_config.continuity_constant = config_.continuity_constant;
        gap_config.overlap_constant = config_.overlap_constant;
        gap_config.reference_upper_bound = reference_upper_bound;
        gap_config.epoch_switch_ratio = config_.tau_ep;
        gap_config.evidence_mode = CandidateGapEvidenceMode::Practical;
        const Clock::time_point solve_start = Clock::now();
        const CandidateDualGapResult gap = build_candidate_dual_gap(
            hierarchy_, operators, load, on_candidate, gap_config);
        ReferenceEpochDualObservation result;
        result.eta_dual_c = gap.eta_dual_c;
        result.L_gap_c = gap.L_gap_c;
        result.time_candidate_dual_solve = seconds_since(solve_start);
        result.time_candidate_dual_operator_assembly = time_operator;
        result.time_candidate_dual_load_assembly = time_load;
        result.time_candidate_dual_prolongation = time_prolongation;
        result.time_candidate_dual_prepare = gap.riesz.prepare_seconds;
        result.time_candidate_dual_patch_solve = gap.riesz.patch_solve_seconds;
        result.time_candidate_dual_reduction = gap.riesz.reduction_seconds;
        result.candidate_dual_patch_factorizations =
            gap.riesz.patch_factorizations;
        result.candidate_dual_parallel_threads = gap.riesz.parallel_threads;
        result.time_candidate_dual = seconds_since(total_start);
        result.time_candidate_wellposedness = time_wellposedness;
        return result;
    }

    void commit_coarse_refinement() override {
        ComplexVector prolonged_warm_start;
        ComplexMatrix prolonged_warm_start_block;
        const std::vector<int> old_free = free_nodes(hierarchy_.coarse_mesh());
        const TriMesh proposed = hierarchy_.proposed_coarse_mesh();
        const std::vector<int> proposed_free = free_nodes(proposed);
        if (localization_warm_start_.size()
                == static_cast<int>(old_free.size())
            && localization_warm_start_.allFinite()) {
            ComplexVector old_nodal = ComplexVector::Zero(
                hierarchy_.coarse_mesh().nodes.size());
            for (int index = 0; index < static_cast<int>(old_free.size()); ++index)
                old_nodal(old_free[index]) = localization_warm_start_(index);
            const ComplexVector proposed_nodal =
                hierarchy_.coarse_nodes_to_proposed_coarse().cast<Complex>()
                * old_nodal;
            prolonged_warm_start.resize(proposed_free.size());
            for (int index = 0;
                 index < static_cast<int>(proposed_free.size()); ++index) {
                prolonged_warm_start(index) = proposed_nodal(proposed_free[index]);
            }
        }
        if (localization_warm_start_block_.rows()
                == static_cast<int>(old_free.size())
            && localization_warm_start_block_.cols() > 0
            && localization_warm_start_block_.allFinite()) {
            const int columns = localization_warm_start_block_.cols();
            ComplexMatrix old_nodal = ComplexMatrix::Zero(
                hierarchy_.coarse_mesh().nodes.size(), columns);
            for (int index = 0; index < static_cast<int>(old_free.size()); ++index)
                old_nodal.row(old_free[index]) =
                    localization_warm_start_block_.row(index);
            const ComplexMatrix proposed_nodal =
                hierarchy_.coarse_nodes_to_proposed_coarse().cast<Complex>()
                * old_nodal;
            prolonged_warm_start_block.resize(proposed_free.size(), columns);
            for (int index = 0;
                 index < static_cast<int>(proposed_free.size()); ++index) {
                prolonged_warm_start_block.row(index) =
                    proposed_nodal.row(proposed_free[index]);
            }
        }
        const ReferenceEpochRefinementResult committed =
            hierarchy_.commit_coarse_refinement();
        if (!committed.changed())
            throw std::runtime_error("coarse proposal was not committable");
        localization_warm_start_ = std::move(prolonged_warm_start);
        localization_warm_start_block_ =
            std::move(prolonged_warm_start_block);
        model_.reset();
        solution_.reset();
        regions_.reset();
        proposed_regions_.reset();
        ++committed_H_steps_;
        const Clock::time_point artifact_start = Clock::now();
        ReferenceEpochMeshSnapshot snapshot;
        snapshot.epoch = epoch_index_;
        snapshot.H_step = committed_H_steps_;
        snapshot.reference_mesh_version =
            hierarchy_.reference_mesh_version();
        snapshot.stage = "committed";
        snapshot.anchor_action =
            ReferenceEpochDriverAction::CommitCoarseRefinement;
        snapshot.anchor_occurrence = commit_occurrence_++;
        snapshot.coarse = hierarchy_.coarse_mesh();
        snapshot.candidate = hierarchy_.candidate_mesh();
        snapshots_.push_back(std::move(snapshot));
        artifact_capture_seconds_ += seconds_since(artifact_start);
    }

    void refresh_reference(const int minimum_post_refresh_level_gap) override {
        const Clock::time_point refresh_start = Clock::now();
        std::optional<HybridReserveDiagnostic> achieved_hybrid_closure;
        enforce_current_mesh_work_limits();
        const ReferenceEpochRefinementGuard proposal_guard =
            proposed_refinement_resource_guard();
        const ReferenceEpochRefinementGuard candidate_guard =
            candidate_refinement_resource_guard();
        if (config_.singularity_hybrid) {
            // First close exact matching against the fixed candidate.  Freeze
            // that proposal and its Omega_F afterwards: alternating proposal
            // matching with uniform regular deepening moves the interface and
            // can create an unbounded chase.
            enforce_current_mesh_work_limits();
            const Clock::time_point closure_start = Clock::now();
            const ProposedHybridMatchingResult matching =
                restore_proposed_hybrid_matching(
                    hierarchy_, current_ell_,
                    config_.hybrid_minimum_physical_radius,
                    HybridMatchingTarget::Candidate, 64,
                    proposal_guard);
            const double matching_seconds = seconds_since(closure_start);
            const std::size_t diagnostic_refresh_index = refresh_occurrence_;
            proposed_regions_ = matching.regions;
            const auto matching_spill_count =
                [&](const ReferenceEpochHierarchy &state) {
                    std::size_t spill = 0;
                    for (const int element :
                         proposed_regions_->omega_f_elements) {
                        if (embedding_children(
                                state.proposed_coarse_elements_to_candidate(),
                                element) != 1) {
                            ++spill;
                        }
                    }
                    return spill;
                };

            HybridGradedReserveProfile reserve;
            ReferenceEpochRefinementResult deepened;
            int margin_before = std::numeric_limits<int>::max();
            int full_gap_before = std::numeric_limits<int>::max();
            int margin_after = std::numeric_limits<int>::max();
            int full_gap_after = std::numeric_limits<int>::max();
            int chosen_collar = -1;
            std::string trial_work_limit;
            const HybridGradedReserveProfile zero_collar_profile =
                make_hybrid_graded_reserve_profile(
                    hierarchy_.proposed_coarse_mesh(), *proposed_regions_,
                    minimum_post_refresh_level_gap, 0);
            const int last_collar = minimum_post_refresh_level_gap > 0
                ? zero_collar_profile.maximum_graph_distance
                    - minimum_post_refresh_level_gap
                : 0;
            if (last_collar < 0) {
                throw ReferenceEpochReserveUnavailable(
                    "hybrid graded reserve unavailable: no far-interior full-target cells");
            }
            for (int collar = 0; collar <= last_collar; ++collar) {
                HybridReserveDiagnostic diagnostic;
                diagnostic.epoch = epoch_index_;
                diagnostic.H_step = committed_H_steps_;
                diagnostic.refresh_index = diagnostic_refresh_index;
                diagnostic.trial_index = collar;
                diagnostic.row_type = "trial";
                diagnostic.requested_target_gap =
                    minimum_post_refresh_level_gap;
                diagnostic.collar = collar;
                diagnostic.ell = current_ell_;
                diagnostic.ell_s = proposed_regions_->l_s;
                diagnostic.physical_radius =
                    config_.hybrid_minimum_physical_radius;
                diagnostic.omega_s_elements =
                    proposed_regions_->omega_s_elements.size();
                diagnostic.omega_f_elements =
                    proposed_regions_->omega_f_elements.size();
                diagnostic.time_matching = matching_seconds;
                HybridGradedReserveProfile trial_profile =
                    make_hybrid_graded_reserve_profile(
                        hierarchy_.proposed_coarse_mesh(), *proposed_regions_,
                        minimum_post_refresh_level_gap, collar);
                diagnostic.maximum_graph_distance =
                    trial_profile.maximum_graph_distance;
                diagnostic.full_target_elements =
                    trial_profile.full_target_elements;
                if (minimum_post_refresh_level_gap > 0
                    && trial_profile.full_target_elements == 0) {
                    diagnostic.status = "rejected";
                    diagnostic.reject_reason = "empty_full_target";
                    hybrid_reserve_diagnostics_.push_back(
                        std::move(diagnostic));
                    continue;
                }
                const int trial_margin_before =
                    hierarchy_.minimum_proposed_candidate_level_gap_margin(
                        trial_profile.target_level_gaps);
                const int trial_full_gap_before =
                    hierarchy_.minimum_proposed_candidate_level_gap(
                        trial_profile.at_full_target);
                diagnostic.profile_margin_before = trial_margin_before;
                diagnostic.far_gap_before = trial_full_gap_before;
                if (trial_margin_before >= 0) {
                    reserve = std::move(trial_profile);
                    margin_before = trial_margin_before;
                    full_gap_before = trial_full_gap_before;
                    margin_after = trial_margin_before;
                    full_gap_after = trial_full_gap_before;
                    chosen_collar = collar;
                    diagnostic.status = "accepted_existing";
                    diagnostic.target_satisfied = true;
                    diagnostic.profile_margin_after = trial_margin_before;
                    diagnostic.far_gap_after = trial_full_gap_before;
                    diagnostic.candidate_elements =
                        hierarchy_.candidate_mesh().elems.size();
                    diagnostic.candidate_nodes =
                        hierarchy_.candidate_mesh().nodes.size();
                    diagnostic.time_total = seconds_since(closure_start);
                    hybrid_reserve_diagnostics_.push_back(
                        std::move(diagnostic));
                    break;
                }

                ReferenceEpochCandidateDeepeningProbe probe;
                const Clock::time_point probe_start = Clock::now();
                try {
                    probe = hierarchy_
                        .probe_candidate_deepening_over_proposed_coarse(
                            trial_profile.target_level_gaps,
                            proposed_regions_->in_omega_f,
                            trial_profile.at_full_target,
                            candidate_guard);
                } catch (const ReferenceEpochWorkLimitExceeded &error) {
                    diagnostic.time_probe = seconds_since(probe_start);
                    diagnostic.time_total = seconds_since(closure_start);
                    diagnostic.status = "rejected";
                    diagnostic.reject_reason = "work_limit:" + std::string(error.what());
                    hybrid_reserve_diagnostics_.push_back(
                        std::move(diagnostic));
                    trial_work_limit = error.what();
                    std::cerr
                        << "[hybrid-refresh-collar-trial] collar=" << collar
                        << " status=work_limit detail=" << error.what()
                        << std::endl;
                    continue;
                }
                diagnostic.time_probe = seconds_since(probe_start);
                diagnostic.target_satisfied = probe.target_satisfied;
                diagnostic.matching_spill = probe.protected_parent_spill;
                diagnostic.profile_margin_after = probe.minimum_gap_margin;
                diagnostic.far_gap_after = probe.minimum_full_target_gap;
                diagnostic.candidate_elements = probe.final_element_count;
                diagnostic.candidate_nodes = probe.final_node_count;
                std::cerr
                    << "[hybrid-refresh-collar-trial] collar=" << collar
                    << " matching_spill=" << probe.protected_parent_spill
                    << " profile_deficit="
                    << std::max(0, -probe.minimum_gap_margin)
                    << " far_regular_gap="
                    << probe.minimum_full_target_gap
                    << " full_target_elements="
                    << trial_profile.full_target_elements
                    << " candidate_elements="
                    << probe.final_element_count
                    << " candidate_nodes=" << probe.final_node_count
                    << std::endl;
                if (!probe.target_satisfied
                    || probe.protected_parent_spill != 0
                    || probe.minimum_gap_margin < 0
                    || probe.minimum_full_target_gap
                        < minimum_post_refresh_level_gap) {
                    diagnostic.status = "rejected";
                    diagnostic.reject_reason = !probe.target_satisfied
                        ? "target_unsatisfied"
                        : probe.protected_parent_spill != 0
                        ? "matching_spill"
                        : probe.minimum_gap_margin < 0
                        ? "profile_deficit"
                        : "far_gap_deficit";
                    diagnostic.time_total = seconds_since(closure_start);
                    hybrid_reserve_diagnostics_.push_back(
                        std::move(diagnostic));
                    continue;
                }

                reserve = std::move(trial_profile);
                margin_before = trial_margin_before;
                full_gap_before = trial_full_gap_before;
                chosen_collar = collar;
                diagnostic.status = "accepted_probe";
                diagnostic.time_total = seconds_since(closure_start);
                hybrid_reserve_diagnostics_.push_back(std::move(diagnostic));
                break;
            }
            if (chosen_collar < 0) {
                if (!trial_work_limit.empty())
                    throw ReferenceEpochWorkLimitExceeded(trial_work_limit);
                throw ReferenceEpochReserveUnavailable(
                    "hybrid graded reserve has no spill-free conformity collar");
            }
            hybrid_conformity_collar_layers_ = chosen_collar;
            const Clock::time_point deepen_start = Clock::now();
            if (margin_before < 0) {
                deepened = hierarchy_.deepen_candidate_over_proposed_coarse(
                    reserve.target_level_gaps, candidate_guard);
            }
            const double deepen_seconds = seconds_since(deepen_start);
            margin_after =
                hierarchy_.minimum_proposed_candidate_level_gap_margin(
                    reserve.target_level_gaps);
            full_gap_after =
                hierarchy_.minimum_proposed_candidate_level_gap(
                    reserve.at_full_target);
            if (matching_spill_count(hierarchy_) != 0) {
                throw std::runtime_error(
                    "hybrid graded reserve entered the frozen matching region");
            }
            if (minimum_post_refresh_level_gap > 0 && margin_after < 0) {
                throw std::runtime_error(
                    "hybrid graded reserve closure left a positive profile deficit");
            }
            HybridReserveDiagnostic closure;
            closure.epoch = epoch_index_;
            closure.H_step = committed_H_steps_;
            closure.refresh_index = diagnostic_refresh_index;
            closure.row_type = "closure";
            closure.status = "achieved";
            closure.requested_target_gap = minimum_post_refresh_level_gap;
            closure.collar = chosen_collar;
            closure.maximum_graph_distance = reserve.maximum_graph_distance;
            closure.ell = current_ell_;
            closure.ell_s = proposed_regions_->l_s;
            closure.physical_radius = config_.hybrid_minimum_physical_radius;
            closure.omega_s_elements = proposed_regions_->omega_s_elements.size();
            closure.omega_f_elements = proposed_regions_->omega_f_elements.size();
            closure.full_target_elements = reserve.full_target_elements;
            closure.target_satisfied = margin_after >= 0
                && full_gap_after >= minimum_post_refresh_level_gap;
            closure.matching_spill = matching_spill_count(hierarchy_);
            closure.profile_margin_before = margin_before;
            closure.profile_margin_after = margin_after;
            closure.far_gap_before = full_gap_before;
            closure.far_gap_after = full_gap_after;
            closure.candidate_elements = hierarchy_.candidate_mesh().elems.size();
            closure.candidate_nodes = hierarchy_.candidate_mesh().nodes.size();
            closure.deepened_elements = deepened.current_element_count
                - deepened.previous_element_count;
            closure.time_matching = matching_seconds;
            closure.time_deepen = deepen_seconds;
            closure.time_total = seconds_since(closure_start);
            achieved_hybrid_closure = std::move(closure);
            std::cerr
                << "[hybrid-refresh-closure] seconds="
                << seconds_since(closure_start)
                << " matching_rounds=" << matching.refinement_rounds
                << " ell_S=" << matching.regions.l_s
                << " conformity_collar=" << chosen_collar
                << " maximum_graph_distance="
                << reserve.maximum_graph_distance
                << " profile_deficit_before=" << std::max(0, -margin_before)
                << " profile_deficit_after=" << std::max(0, -margin_after)
                << " far_regular_gap_before=" << full_gap_before
                << " far_regular_gap_after=" << full_gap_after
                << " full_target_elements=" << reserve.full_target_elements
                << " deepened_elements="
                << deepened.current_element_count
                    - deepened.previous_element_count
                << " candidate_elements="
                << hierarchy_.candidate_mesh().elems.size()
                << " proposed_elements="
                << hierarchy_.proposed_coarse_mesh().elems.size()
                << std::endl;
        } else {
            hierarchy_.deepen_candidate_over_proposed_coarse(
                minimum_post_refresh_level_gap, candidate_guard);
        }
        enforce_reference_promotion_work_limit();
        const Clock::time_point artifact_start = Clock::now();
        ReferenceEpochMeshSnapshot pending_snapshot;
        pending_snapshot.epoch = epoch_index_;
        pending_snapshot.H_step = committed_H_steps_;
        pending_snapshot.reference_mesh_version =
            hierarchy_.reference_mesh_version();
        pending_snapshot.stage = "pre_switch";
        pending_snapshot.anchor_action =
            ReferenceEpochDriverAction::RefreshReference;
        pending_snapshot.anchor_occurrence = refresh_occurrence_;
        pending_snapshot.coarse = hierarchy_.proposed_coarse_mesh();
        pending_snapshot.candidate = hierarchy_.candidate_mesh();
        artifact_capture_seconds_ += seconds_since(artifact_start);
        hierarchy_.refresh_reference_from_candidate();
        if (config_.singularity_hybrid) {
            if (!proposed_regions_)
                throw std::logic_error(
                    "hybrid refresh lost its prospective regions");
            for (const int element : proposed_regions_->omega_f_elements) {
                if (embedding_children(
                        hierarchy_.proposed_coarse_elements_to_reference(),
                        element) != 1) {
                    throw std::runtime_error(
                        "hybrid refresh did not preserve exact matching");
                }
            }
            const HybridGradedReserveProfile reserve =
                make_hybrid_graded_reserve_profile(
                    hierarchy_.proposed_coarse_mesh(), *proposed_regions_,
                    minimum_post_refresh_level_gap,
                    hybrid_conformity_collar_layers_);
            const int refreshed_margin =
                hierarchy_.minimum_proposed_reference_level_gap_margin(
                    reserve.target_level_gaps);
            const int refreshed_full_gap =
                hierarchy_.minimum_proposed_reference_level_gap(
                    reserve.at_full_target);
            if (minimum_post_refresh_level_gap > 0
                && refreshed_margin != std::numeric_limits<int>::max()
                && refreshed_margin < 0) {
                throw std::runtime_error(
                    "hybrid refresh did not preserve the graded reserve target");
            }
            if (minimum_post_refresh_level_gap > 0
                && refreshed_full_gap != std::numeric_limits<int>::max()
                && refreshed_full_gap < minimum_post_refresh_level_gap) {
                throw std::runtime_error(
                    "hybrid refresh did not preserve the far-regular reserve target");
            }
            if (!achieved_hybrid_closure) {
                throw std::logic_error(
                    "hybrid refresh completed without a closure diagnostic");
            }
            achieved_hybrid_closure->time_total =
                seconds_since(refresh_start);
            hybrid_reserve_diagnostics_.push_back(
                std::move(*achieved_hybrid_closure));
        }
        // Publish the pre-switch artifact only after promotion and every
        // post-refresh audit succeed.  Otherwise the driver has no matching
        // RefreshReference journal row and an orphan snapshot would mask the
        // original failure while writing mesh_manifest.csv.
        const Clock::time_point artifact_publish_start = Clock::now();
        snapshots_.push_back(std::move(pending_snapshot));
        ++refresh_occurrence_;
        artifact_capture_seconds_ += seconds_since(artifact_publish_start);
        gram_factor_cache_.clear();
        reference_solution_.reset();
        reference_load_.reset();
        reference_discrete_norm_.reset();
        reference_exact_relative_energy_.reset();
        corrector_cache_.clear();
        // The reference operator changed, but the coarse coefficient space
        // did not. Keep the dominant coarse vector as an initial guess; the
        // subsequent prospective coarse commit prolongates it if needed.
    }

    ReferenceEpochResourceSnapshot resources() const override {
        ReferenceEpochResourceSnapshot result;
        result.coarse_unknowns = unknowns(hierarchy_.coarse_mesh());
        result.reference_unknowns = unknowns(hierarchy_.reference_mesh());
        result.candidate_unknowns = unknowns(hierarchy_.candidate_mesh());
        result.kappa_H_max = config_.wavenumber
            * max_element_diameter(hierarchy_.coarse_mesh());
        result.artifact_capture_seconds = artifact_capture_seconds_;
        return result;
    }

    const ReferenceEpochHierarchy &hierarchy() const { return hierarchy_; }
    const std::vector<ReferenceEpochMeshSnapshot> &snapshots() const {
        return snapshots_;
    }
    const std::optional<TriMesh> &latest_solved_coarse() const {
        return latest_solved_coarse_;
    }
    const std::vector<HybridReserveDiagnostic> &hybrid_reserve_diagnostics() const {
        return hybrid_reserve_diagnostics_;
    }

private:
    [[noreturn]] static void throw_reference_work_limit() {
        throw ReferenceEpochWorkLimitExceeded(
            "maximum_reference_unknowns reached");
    }

    [[noreturn]] static void throw_candidate_work_limit() {
        throw ReferenceEpochWorkLimitExceeded(
            "maximum_candidate_unknowns reached");
    }

    void enforce_current_mesh_work_limits() const {
        if (unknowns(hierarchy_.reference_mesh())
            > config_.work_limits.maximum_reference_unknowns) {
            throw_reference_work_limit();
        }
        if (unknowns(hierarchy_.candidate_mesh())
            > config_.work_limits.maximum_candidate_unknowns) {
            throw_candidate_work_limit();
        }
    }

    ReferenceEpochRefinementGuard
    proposed_refinement_resource_guard() const {
        return [this](
                   const ReferenceEpochRefinementGuardPoint,
                   const TriMesh &proposed_mesh) {
            enforce_current_mesh_work_limits();
            const std::size_t proposed_unknowns = unknowns(proposed_mesh);
            // Any future refreshed reference and candidate must contain the
            // proposal.  Stop before spending on embeddings when this is
            // already impossible within either configured budget.
            if (proposed_unknowns
                > config_.work_limits.maximum_reference_unknowns) {
                throw_reference_work_limit();
            }
            if (proposed_unknowns
                > config_.work_limits.maximum_candidate_unknowns) {
                throw_candidate_work_limit();
            }
        };
    }

    ReferenceEpochRefinementGuard
    candidate_refinement_resource_guard() const {
        return [this](
                   const ReferenceEpochRefinementGuardPoint point,
                   const TriMesh &candidate_mesh) {
            if (unknowns(hierarchy_.reference_mesh())
                > config_.work_limits.maximum_reference_unknowns) {
                throw_reference_work_limit();
            }
            const std::size_t candidate_unknowns = unknowns(candidate_mesh);
            const bool exhausted =
                point == ReferenceEpochRefinementGuardPoint::BeforeNvb
                ? candidate_unknowns
                    >= config_.work_limits.maximum_candidate_unknowns
                : candidate_unknowns
                    > config_.work_limits.maximum_candidate_unknowns;
            if (exhausted) throw_candidate_work_limit();
        };
    }

    void enforce_reference_promotion_work_limit() const {
        enforce_current_mesh_work_limits();
        if (unknowns(hierarchy_.candidate_mesh())
            > config_.work_limits.maximum_reference_unknowns) {
            throw_reference_work_limit();
        }
    }

    void ensure_reference_solution(const ComplexVector &load) {
        if (!reference_solution_)
            reference_solution_ = solve_helmholtz_fem(model_->operators(), load);
    }

    ReferenceEpochPaperConfig config_;
    PaperCaseData data_;
    ReferenceEpochHierarchy hierarchy_;
    // E1 patches grow well beyond the library default 4096-DoF admission
    // limit after the first refresh. Admit the observed patches, but bound
    // retained assembled systems and solutions by an explicit two-GiB LRU
    // budget; a count-only 64-entry pilot exhausted a 12-GiB WSL instance.
    HelmholtzCorrectorPatchCache corrector_cache_{
        256, 32768, std::size_t{2} * 1024 * 1024 * 1024};
    // Numeric constraints change on most patches after a coarse commit.
    // Keep only a small tail: a 512-entry pilot retained several GiB for
    // very few hits, while 32 bounds the optional cross-step memory cost.
    ReferenceDefectGramFactorCache gram_factor_cache_{32};
    ComplexVector localization_warm_start_;
    ComplexMatrix localization_warm_start_block_;
    std::unique_ptr<HelmholtzLodModel> model_;
    std::optional<HelmholtzLodSolution> solution_;
    std::optional<ComplexVector> reference_solution_;
    std::optional<ComplexVector> reference_load_;
    std::optional<HelmholtzError> reference_discrete_norm_;
    std::optional<double> reference_exact_relative_energy_;
    std::optional<HelmholtzError> reference_exact_norm_;
    std::optional<SingularRegionClassification> regions_;
    std::optional<SingularRegionClassification> proposed_regions_;
    std::optional<ReferenceEpochCoarseRefinementPreview>
        accepted_coarse_preview_;
    std::optional<TriMesh> latest_solved_coarse_;
    bool epoch_started_ = false;
    bool epoch_start_snapshot_pending_ = false;
    int current_ell_ = 0;
    int hybrid_conformity_collar_layers_ = 1;
    std::size_t epoch_index_ = 0;
    std::size_t committed_H_steps_ = 0;
    std::size_t begin_occurrence_ = 0;
    std::size_t commit_occurrence_ = 0;
    std::size_t refresh_occurrence_ = 0;
    std::vector<ReferenceEpochMeshSnapshot> snapshots_;
    std::vector<HybridReserveDiagnostic> hybrid_reserve_diagnostics_;
    double artifact_capture_seconds_ = 0.0;
};

void write_iterations(
    const std::filesystem::path &path,
    const ReferenceEpochPaperConfig &config,
    const std::string &run_id,
    const ReferenceEpochDriverResult &result) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "schema_version,manuscript_sha256,code_commit,case,method,kappa,run_id,repeat_index,"
           "epoch,event_sequence,H_step,state,action,stop_reason,decision_detail,"
           "N_H,N_reference,N_candidate,ell,kappa_H_max,"
           "eta_H,Theta_loc,delta_loc_hat,U_practical,eta_eq_c,eta_eq_c_F,eta_eq_c_R,"
           "indicator_mass_c_F,indicator_mass_c_R,marked_mass_c_F,marked_mass_c_R,"
           "eta_dual_c,L_gap_c,dual_check_performed,"
           "minimum_reference_level_gap,structural_dual_trigger,"
           "level_gap_dual_trigger,tolerance_dual_trigger,"
           "decrease_dual_trigger,interval_dual_trigger,"
           "relative_reference_energy,relative_exact_energy,relative_exact_L2,"
           "relative_exact_reference_energy,marked_H,"
           "hybrid_regular_indicator_mass,hybrid_admissible_indicator_mass,"
           "hybrid_marked_H_indicator_mass,hybrid_coarse_conformity_collar,"
           "hybrid_coarse_marking_closure_safe,hybrid_full_regular_doerfler,"
           "hybrid_coarse_preview_attempts,hybrid_coarse_preview_cached,"
           "marked_c,"
           "marked_c_F,marked_c_R,candidate_cells_F,candidate_cells_R,"
           "active_correctors,rebuilt_correctors,reused_correctors,"
           "corrector_cache_oversized_misses,corrector_cache_budget_rejections,"
           "corrector_cache_entries,corrector_cache_current_bytes,"
           "corrector_cache_peak_bytes,"
           "skipped_correctors,skipped_corrector_work_units,"
           "corrector_parallel_threads,corrector_patch_assembly_work,"
           "corrector_patch_solve_work,corrector_patch_pack_work,"
           "corrector_maximum_patch_dofs,corrector_maximum_patch_constraints,"
           "corrector_maximum_patch_rhs,"
           "hybrid_ell_S,hybrid_minimum_physical_radius,hybrid_covered_physical_radius,"
           "hybrid_omega_S_elements,hybrid_omega_F_elements,"
           "time_corrector_check_total,time_lod_build_total,"
           "time_lod_mesh_and_interpolation,time_lod_operators,time_corrector,"
           "time_lod_corrected_basis,time_lod_coarse_operator,"
           "time_lod_coarse_factorization,time_reference_stability,time_theta,"
           "time_gram_prepare_structure,time_gram_prepare_factorization,"
           "time_gram_action_rhs,time_gram_action_patch_solve,time_gram_action_scatter,"
           "gram_action_calls,gram_patch_factorizations,gram_factor_cache_hits,"
           "gram_factor_cache_misses,gram_structure_parallel_threads,"
           "gram_parallel_threads,"
           "localization_iterations,localization_relative_residual,"
           "localization_used_warm_start,time_lod_solve,"
           "time_reference_riesz,time_hybrid_coarse_marking,"
           "time_hybrid_coarse_preview,time_hybrid_coarse_preview_nvb,"
           "time_hybrid_coarse_preview_reference_embedding,"
           "time_reference_validation,"
           "time_exact_validation,time_candidate_enrich_total,time_candidate_close,"
           "time_candidate_operator_assembly,time_candidate_prolongation,"
           "time_candidate_flux_reconstruction,time_candidate_flux_prepare,"
           "time_candidate_flux_patch_solve,time_candidate_flux_merge,"
           "time_candidate_flux_audit,candidate_flux_parallel_threads,"
           "time_candidate_enrich,time_candidate_nvb_refine,"
           "time_candidate_embedding_composition,time_candidate_parent_map_update,"
           "time_candidate_quasi_interpolation,time_candidate_embedding_validation,"
           "time_candidate_dual,time_candidate_wellposedness,"
           "time_candidate_dual_operator_assembly,"
           "time_candidate_dual_load_assembly,time_candidate_dual_prolongation,"
           "time_candidate_dual_solve,time_candidate_dual_prepare,"
           "time_candidate_dual_patch_solve,time_candidate_dual_reduction,"
           "candidate_dual_patch_factorizations,candidate_dual_parallel_threads,"
           "time_mesh,time_validation_cumulative,"
           "time_artifact_capture_cumulative,time_method_cumulative,"
           "time_total_cumulative,peak_memory_mb\n";
    for (const ReferenceEpochDriverRecord &row : result.journal) {
        out << config.schema_version << ',' << config.manuscript_sha256 << ','
            << config.git_commit << ',' << to_string(config.case_id) << ','
            << config.method << ',' << number(config.wavenumber) << ',' << run_id
            << ',' << config.repeat_index << ',' << row.epoch << ',' << row.sequence
            << ',' << row.H_step
            << ',' << reference_epoch_driver_state_name(row.state_after) << ','
            << reference_epoch_driver_action_name(row.action) << ','
            << csv(row.action == ReferenceEpochDriverAction::Fail
                       || row.action == ReferenceEpochDriverAction::StopWorkLimit
                       ? row.detail : std::string{})
            << ',' << csv(row.detail)
            << ',' << row.coarse_unknowns << ',' << row.reference_unknowns << ','
            << row.candidate_unknowns << ',' << row.ell << ','
            << number(row.kappa_H_max) << ',' << number(row.eta_H) << ','
            << number(row.theta_loc) << ',' << number(row.delta_loc_hat) << ','
            << number(row.U_practical) << ','
            << number(row.eta_eq_c) << ',' << number(row.eta_eq_c_f) << ','
            << number(row.eta_eq_c_r) << ','
            << number(row.indicator_mass_c_f) << ','
            << number(row.indicator_mass_c_r) << ','
            << number(row.marked_mass_c_f) << ','
            << number(row.marked_mass_c_r) << ','
            << number(row.eta_dual_c) << ','
            << number(row.L_gap_c) << ','
            << (row.action == ReferenceEpochDriverAction::ComputeCandidateDual
                    ? "true" : "false")
            << ','
            << (row.minimum_reference_level_gap
                        == std::numeric_limits<int>::max()
                    ? "NA"
                    : std::to_string(row.minimum_reference_level_gap))
            << ',' << (row.structural_dual_trigger ? "true" : "false")
            << ',' << (row.level_gap_dual_trigger ? "true" : "false")
            << ',' << (row.tolerance_dual_trigger ? "true" : "false")
            << ',' << (row.decrease_dual_trigger ? "true" : "false")
            << ',' << (row.interval_dual_trigger ? "true" : "false")
            << ',' << number(row.relative_reference_energy) << ','
            << number(row.relative_exact_energy) << ','
            << number(row.relative_exact_L2) << ','
            << number(row.relative_exact_reference_energy) << ','
            << row.marked_H << ','
            << number(row.hybrid_regular_indicator_mass) << ','
            << number(row.hybrid_admissible_indicator_mass) << ','
            << number(row.hybrid_marked_H_indicator_mass) << ','
            << (row.hybrid_coarse_conformity_collar < 0
                    ? "NA"
                    : std::to_string(row.hybrid_coarse_conformity_collar))
            << ','
            << (row.hybrid_coarse_marking_closure_safe ? "true" : "false")
            << ',' << (row.hybrid_full_regular_doerfler ? "true" : "false")
            << ',' << row.hybrid_coarse_preview_attempts
            << ',' << (row.hybrid_coarse_preview_cached ? "true" : "false")
            << ','
            << row.marked_c << ',' << row.marked_c_f << ','
            << row.marked_c_r << ',' << row.candidate_cells_f << ','
            << row.candidate_cells_r << ',' << row.active_correctors << ','
            << row.rebuilt_correctors << ',' << row.reused_correctors << ','
            << row.corrector_cache_oversized_misses << ','
            << row.corrector_cache_budget_rejections << ','
            << row.corrector_cache_entries << ','
            << row.corrector_cache_current_bytes << ','
            << row.corrector_cache_peak_bytes << ','
            << row.skipped_correctors << ','
            << row.skipped_corrector_work_units << ','
            << row.corrector_parallel_threads << ','
            << number(row.corrector_patch_assembly_work_seconds) << ','
            << number(row.corrector_patch_solve_work_seconds) << ','
            << number(row.corrector_patch_pack_work_seconds) << ','
            << row.corrector_maximum_patch_dofs << ','
            << row.corrector_maximum_patch_constraints << ','
            << row.corrector_maximum_patch_rhs << ','
            << (row.hybrid_l_s < 0 ? "NA" : std::to_string(row.hybrid_l_s)) << ','
            << number(row.hybrid_minimum_physical_radius) << ','
            << number(row.hybrid_covered_physical_radius) << ','
            << row.hybrid_omega_s_elements << ','
            << row.hybrid_omega_f_elements << ','
            << number(row.time_corrector_check_total) << ','
            << number(row.time_lod_build_total) << ','
            << number(row.time_lod_mesh_and_interpolation) << ','
            << number(row.time_lod_operators) << ','
            << number(row.time_corrector) << ','
            << number(row.time_lod_corrected_basis) << ','
            << number(row.time_lod_coarse_operator) << ','
            << number(row.time_lod_coarse_factorization) << ','
            << number(row.time_reference_stability) << ','
            << number(row.time_theta) << ','
            << number(row.time_gram_prepare_structure) << ','
            << number(row.time_gram_prepare_factorization) << ','
            << number(row.time_gram_action_rhs) << ','
            << number(row.time_gram_action_patch_solve) << ','
            << number(row.time_gram_action_scatter) << ','
            << row.gram_action_calls << ',' << row.gram_patch_factorizations << ','
            << row.gram_factor_cache_hits << ','
            << row.gram_factor_cache_misses << ','
            << row.gram_structure_parallel_threads << ','
            << row.gram_parallel_threads << ',' << row.localization_iterations << ','
            << number(row.localization_relative_residual) << ','
            << (row.localization_used_warm_start ? "true" : "false") << ','
            << number(row.time_lod_solve) << ','
            << number(row.time_reference_riesz) << ','
            << number(row.time_hybrid_coarse_marking) << ','
            << number(row.time_hybrid_coarse_preview) << ','
            << number(row.time_hybrid_coarse_preview_nvb) << ','
            << number(row.time_hybrid_coarse_preview_reference_embedding) << ','
            << number(row.time_reference_validation) << ','
            << number(row.time_exact_validation) << ','
            << number(row.time_candidate_flux) << ','
            << number(row.time_candidate_close) << ','
            << number(row.time_candidate_operator_assembly) << ','
            << number(row.time_candidate_prolongation) << ','
            << number(row.time_candidate_flux_reconstruction) << ','
            << number(row.time_candidate_flux_prepare) << ','
            << number(row.time_candidate_flux_patch_solve) << ','
            << number(row.time_candidate_flux_merge) << ','
            << number(row.time_candidate_flux_audit) << ','
            << row.candidate_flux_parallel_threads << ','
            << number(row.time_candidate_enrich) << ','
            << number(row.time_candidate_nvb_refine) << ','
            << number(row.time_candidate_embedding_composition) << ','
            << number(row.time_candidate_parent_map_update) << ','
            << number(row.time_candidate_quasi_interpolation) << ','
            << number(row.time_candidate_embedding_validation) << ','
            << number(row.time_candidate_dual) << ','
            << number(row.time_candidate_wellposedness) << ','
            << number(row.time_candidate_dual_operator_assembly) << ','
            << number(row.time_candidate_dual_load_assembly) << ','
            << number(row.time_candidate_dual_prolongation) << ','
            << number(row.time_candidate_dual_solve) << ','
            << number(row.time_candidate_dual_prepare) << ','
            << number(row.time_candidate_dual_patch_solve) << ','
            << number(row.time_candidate_dual_reduction) << ','
            << row.candidate_dual_patch_factorizations << ','
            << row.candidate_dual_parallel_threads << ','
            << number(row.time_mesh) << ','
            << number(row.time_validation_cumulative) << ','
            << number(row.time_artifact_capture_cumulative) << ','
            << number(row.time_method_cumulative)
            << ',' << number(row.time_total_cumulative) << ','
            << number(peak_memory_mb()) << '\n';
    }
}

void write_auxiliary_outputs(
    const std::filesystem::path &directory,
    const ReferenceEpochPaperConfig &config,
    const std::string &run_id,
    const ReferenceEpochDriverResult &result,
    const ReferenceEpochHierarchy &hierarchy,
    const std::vector<ReferenceEpochMeshSnapshot> &snapshots,
    const std::optional<TriMesh> &latest_solved_coarse,
    const std::vector<HybridReserveDiagnostic> &hybrid_reserve_diagnostics) {
    {
        std::ofstream out(directory / "summary.csv");
        out << "schema_version,run_id,row_type,target,achieved,epoch,iteration,N_H,relative_energy,time_method_cumulative,time_total_cumulative,status,stop_reason\n";
        for (double target : {0.1, 0.05, 0.02, 0.01}) {
            const auto hit = std::find_if(
                result.journal.begin(), result.journal.end(),
                [&](const ReferenceEpochDriverRecord &row) {
                    const double error = std::isfinite(row.relative_exact_energy)
                        ? row.relative_exact_energy
                        : row.relative_reference_energy;
                    return std::isfinite(error) && error <= target;
                });
            out << config.schema_version << ',' << run_id << ",target,"
                << number(target) << ',' << (hit != result.journal.end() ? "true" : "false");
            if (hit != result.journal.end()) {
                const double error = std::isfinite(hit->relative_exact_energy)
                    ? hit->relative_exact_energy
                    : hit->relative_reference_energy;
                out << ',' << hit->epoch << ',' << hit->sequence << ','
                    << hit->coarse_unknowns << ',' << number(error) << ','
                    << number(hit->time_method_cumulative) << ','
                    << number(hit->time_total_cumulative);
            } else {
                out << ",NA,NA,NA,NA,NA,NA";
            }
            out << ',' << reference_epoch_driver_state_name(result.state)
                << ',' << csv(result.stop_reason) << '\n';
        }
        out << config.schema_version << ',' << run_id
            << ",run,NA,NA,NA,NA,NA,NA,NA,NA,"
            << reference_epoch_driver_state_name(result.state) << ','
            << csv(result.stop_reason) << '\n';
    }
    {
        std::ofstream out(directory / "epoch_history.csv");
        out << "schema_version,epoch,iteration,action,N_H,N_reference,N_candidate,L_gap_c\n";
        for (const auto &row : result.journal) {
            if (row.action == ReferenceEpochDriverAction::BeginEpoch
                || row.action == ReferenceEpochDriverAction::RefreshReference) {
                out << config.schema_version << ',' << row.epoch << ',' << row.sequence
                    << ',' << reference_epoch_driver_action_name(row.action) << ','
                    << row.coarse_unknowns << ',' << row.reference_unknowns << ','
                    << row.candidate_unknowns << ',' << number(row.L_gap_c) << '\n';
            }
        }
    }
    {
        std::ofstream out(directory / "corrector_work.csv");
        out << "schema_version,epoch,event_sequence,H_step,action,ell,Theta_loc,"
               "delta_loc_hat,decision_detail,hybrid_ell_S,"
               "hybrid_minimum_physical_radius,hybrid_covered_physical_radius,"
               "hybrid_omega_S_elements,hybrid_omega_F_elements,"
               "active,rebuilt,reused,corrector_cache_oversized_misses,"
               "corrector_cache_budget_rejections,"
               "corrector_cache_entries,corrector_cache_current_bytes,"
               "corrector_cache_peak_bytes,"
               "skipped,skipped_work_units,"
               "corrector_parallel_threads,corrector_patch_assembly_work,"
               "corrector_patch_solve_work,corrector_patch_pack_work,"
               "corrector_maximum_patch_dofs,corrector_maximum_patch_constraints,"
               "corrector_maximum_patch_rhs,"
               "time_corrector_check_total,time_lod_build_total,"
               "time_lod_mesh_and_interpolation,time_lod_operators,"
               "time_corrector,time_lod_corrected_basis,"
               "time_lod_coarse_operator,time_lod_coarse_factorization,"
               "time_reference_stability,time_theta,"
               "time_gram_prepare_structure,"
               "time_gram_prepare_factorization,time_gram_action_rhs,"
               "time_gram_action_patch_solve,time_gram_action_scatter,"
               "gram_action_calls,gram_patch_factorizations,"
               "gram_factor_cache_hits,gram_factor_cache_misses,"
               "gram_structure_parallel_threads,gram_parallel_threads,"
               "localization_iterations,localization_relative_residual,"
               "localization_used_warm_start\n";
        for (const auto &row : result.journal) {
            if (row.action == ReferenceEpochDriverAction::AcceptCorrector
                || row.action == ReferenceEpochDriverAction::IncreaseGlobalEll) {
                out << config.schema_version << ',' << row.epoch << ',' << row.sequence
                    << ',' << row.H_step << ','
                    << reference_epoch_driver_action_name(row.action) << ','
                    << row.ell << ',' << number(row.theta_loc) << ','
                    << number(row.delta_loc_hat) << ',' << csv(row.detail) << ','
                    << (row.hybrid_l_s < 0 ? "NA" : std::to_string(row.hybrid_l_s))
                    << ',' << number(row.hybrid_minimum_physical_radius)
                    << ',' << number(row.hybrid_covered_physical_radius)
                    << ',' << row.hybrid_omega_s_elements
                    << ',' << row.hybrid_omega_f_elements
                    << ',' << row.active_correctors
                    << ',' << row.rebuilt_correctors
                    << ',' << row.reused_correctors << ','
                    << row.corrector_cache_oversized_misses << ','
                    << row.corrector_cache_budget_rejections << ','
                    << row.corrector_cache_entries << ','
                    << row.corrector_cache_current_bytes << ','
                    << row.corrector_cache_peak_bytes << ','
                    << row.skipped_correctors << ',' << row.skipped_corrector_work_units
                    << ',' << row.corrector_parallel_threads
                    << ',' << number(row.corrector_patch_assembly_work_seconds)
                    << ',' << number(row.corrector_patch_solve_work_seconds)
                    << ',' << number(row.corrector_patch_pack_work_seconds)
                    << ',' << row.corrector_maximum_patch_dofs
                    << ',' << row.corrector_maximum_patch_constraints
                    << ',' << row.corrector_maximum_patch_rhs
                    << ',' << number(row.time_corrector_check_total)
                    << ',' << number(row.time_lod_build_total)
                    << ',' << number(row.time_lod_mesh_and_interpolation)
                    << ',' << number(row.time_lod_operators)
                    << ',' << number(row.time_corrector)
                    << ',' << number(row.time_lod_corrected_basis)
                    << ',' << number(row.time_lod_coarse_operator)
                    << ',' << number(row.time_lod_coarse_factorization) << ','
                    << number(row.time_reference_stability) << ','
                    << number(row.time_theta) << ','
                    << number(row.time_gram_prepare_structure) << ','
                    << number(row.time_gram_prepare_factorization) << ','
                    << number(row.time_gram_action_rhs) << ','
                    << number(row.time_gram_action_patch_solve) << ','
                    << number(row.time_gram_action_scatter) << ','
                    << row.gram_action_calls << ','
                    << row.gram_patch_factorizations << ','
                    << row.gram_factor_cache_hits << ','
                    << row.gram_factor_cache_misses << ','
                    << row.gram_structure_parallel_threads << ','
                    << row.gram_parallel_threads << ','
                    << row.localization_iterations << ','
                    << number(row.localization_relative_residual) << ','
                    << (row.localization_used_warm_start ? "true" : "false")
                    << '\n';
            }
        }
    }
    {
        std::ofstream out(directory / "hybrid_reserve.csv");
        out << "schema_version,epoch,H_step,refresh_index,trial_index,row_type,status,"
               "reject_reason,requested_target_gap,collar,maximum_graph_distance,"
               "ell,ell_S,R_star,omega_S_elements,omega_F_elements,"
               "full_target_elements,target_satisfied,matching_spill,"
               "profile_margin_before,profile_margin_after,far_gap_before,"
               "far_gap_after,candidate_elements,candidate_nodes,deepened_elements,"
               "time_matching,time_probe,time_deepen,time_total\n";
        const auto integer = [](const int value) {
            return value == std::numeric_limits<int>::max()
                ? std::string("NA") : std::to_string(value);
        };
        for (const HybridReserveDiagnostic &row : hybrid_reserve_diagnostics) {
            out << config.schema_version << ',' << row.epoch << ',' << row.H_step
                << ',' << row.refresh_index << ','
                << (row.trial_index < 0 ? "NA" : std::to_string(row.trial_index))
                << ',' << row.row_type << ',' << row.status << ','
                << csv(row.reject_reason) << ',' << row.requested_target_gap << ','
                << (row.collar < 0 ? "NA" : std::to_string(row.collar)) << ','
                << row.maximum_graph_distance << ',' << row.ell << ','
                << (row.ell_s < 0 ? "NA" : std::to_string(row.ell_s)) << ','
                << number(row.physical_radius) << ',' << row.omega_s_elements
                << ',' << row.omega_f_elements << ',' << row.full_target_elements
                << ',' << (row.target_satisfied ? "true" : "false") << ','
                << row.matching_spill << ',' << integer(row.profile_margin_before)
                << ',' << integer(row.profile_margin_after) << ','
                << integer(row.far_gap_before) << ',' << integer(row.far_gap_after)
                << ',' << row.candidate_elements << ',' << row.candidate_nodes
                << ',' << row.deepened_elements << ',' << number(row.time_matching)
                << ',' << number(row.time_probe) << ',' << number(row.time_deepen)
                << ',' << number(row.time_total) << '\n';
        }
    }
    const std::string coarse_name = "mesh_final_coarse.vtu";
    const std::string reference_name = "mesh_final_reference.vtu";
    const std::string candidate_name = "mesh_final_candidate.vtu";
    io::write_vtu(directory / coarse_name, hierarchy.coarse_mesh());
    io::write_vtu(directory / reference_name, hierarchy.reference_mesh());
    io::write_vtu(directory / candidate_name, hierarchy.candidate_mesh());
    {
        std::ofstream out(directory / "mesh_manifest.csv");
        out << "case,epoch,H_step,iteration,stage,mesh_role,filename,"
               "N_cells,N_dofs,reference_mesh_version\n";
        const auto anchor_sequence = [&](const ReferenceEpochMeshSnapshot &snapshot) {
            std::size_t occurrence = 0;
            for (const ReferenceEpochDriverRecord &row : result.journal) {
                if (row.action != snapshot.anchor_action) continue;
                if (occurrence++ == snapshot.anchor_occurrence) return row.sequence;
            }
            throw std::runtime_error("mesh snapshot has no matching journal action");
        };
        const auto experiment_name = [&] {
            return config.case_id == PaperCase::R1
                ? std::string("E1") : (config.case_id == PaperCase::S
                    ? std::string("E2")
                    : std::string(to_string(config.case_id)));
        };
        const auto snapshot_prefix = [&](const std::size_t epoch) {
            std::ostringstream prefix;
            prefix << "mesh_" << experiment_name() << "_e"
                   << std::setfill('0') << std::setw(3) << epoch;
            return prefix.str();
        };
        // Store a fixed reference mesh once per epoch.  Per-H-step manifest
        // rows below point to this same file and version, so plotting can show
        // the reference in every column without multiplying memory or I/O.
        for (const ReferenceEpochMeshSnapshot &snapshot : snapshots) {
            if (!snapshot.reference) continue;
            const std::string name = snapshot_prefix(snapshot.epoch)
                + "_reference.vtu";
            io::write_vtu(directory / name, *snapshot.reference);
        }
        for (const ReferenceEpochMeshSnapshot &snapshot : snapshots) {
            const bool selected = snapshot.epoch == 0
                || snapshot.stage == "epoch_start"
                || snapshot.stage == "pre_switch";
            if (!selected) continue;
            const std::size_t iteration = anchor_sequence(snapshot);
            const std::string prefix = snapshot_prefix(snapshot.epoch);
            const auto reference = std::find_if(
                snapshots.begin(), snapshots.end(),
                [&](const ReferenceEpochMeshSnapshot &entry) {
                    return entry.epoch == snapshot.epoch
                        && entry.reference.has_value();
                });
            if (reference == snapshots.end())
                throw std::runtime_error(
                    "E1/E2 mesh snapshot has no fixed epoch reference");
            if (reference->reference_mesh_version
                != snapshot.reference_mesh_version) {
                // A refresh commit is journalled as the final action of the
                // old epoch, although the backend has already promoted the
                // candidate to the next reference.  The following
                // BeginEpoch snapshot owns that same H-step and new version.
                if (snapshot.stage == "committed") continue;
                throw std::runtime_error(
                    "reference mesh version changed inside one epoch");
            }
            const std::string reference_snapshot = prefix + "_reference.vtu";
            out << to_string(config.case_id) << ',' << snapshot.epoch << ','
                << snapshot.H_step << ',' << iteration << ',' << snapshot.stage
                << ",reference," << reference_snapshot << ','
                << reference->reference->elems.size() << ','
                << unknowns(*reference->reference) << ','
                << snapshot.reference_mesh_version << '\n';
            std::ostringstream iteration_text;
            iteration_text << std::setfill('0') << std::setw(3) << iteration;
            const std::string coarse_snapshot = prefix + "_i"
                + iteration_text.str() + "_coarse.vtu";
            const std::string candidate_snapshot = prefix + "_i"
                + iteration_text.str() + "_candidate.vtu";
            io::write_vtu(directory / coarse_snapshot, snapshot.coarse);
            io::write_vtu(directory / candidate_snapshot, snapshot.candidate);
            out << to_string(config.case_id) << ',' << snapshot.epoch << ','
                << snapshot.H_step << ',' << iteration << ',' << snapshot.stage
                << ",coarse,"
                << coarse_snapshot << ',' << snapshot.coarse.elems.size() << ','
                << unknowns(snapshot.coarse) << ','
                << snapshot.reference_mesh_version << '\n';
            out << to_string(config.case_id) << ',' << snapshot.epoch << ','
                << snapshot.H_step << ',' << iteration << ',' << snapshot.stage
                << ",candidate,"
                << candidate_snapshot << ',' << snapshot.candidate.elems.size()
                << ',' << unknowns(snapshot.candidate) << ','
                << snapshot.reference_mesh_version << '\n';
        }
        const std::size_t iteration = result.journal.empty()
            ? 0 : result.journal.back().sequence;
        const std::size_t epoch = result.journal.empty()
            ? 0 : result.journal.back().epoch;
        const std::size_t H_step = result.journal.empty()
            ? 0 : result.journal.back().H_step;
        const auto row = [&](const char *role, const std::string &name,
                             const TriMesh &mesh) {
            out << to_string(config.case_id) << ',' << epoch << ',' << H_step
                << ',' << iteration << ",final," << role << ',' << name << ','
                << mesh.elems.size() << ',' << unknowns(mesh) << ','
                << hierarchy.reference_mesh_version() << '\n';
        };
        row("coarse", coarse_name, hierarchy.coarse_mesh());
        row("reference", reference_name, hierarchy.reference_mesh());
        row("candidate", candidate_name, hierarchy.candidate_mesh());
        if (config.singularity_hybrid && latest_solved_coarse) {
            const auto solved = std::find_if(
                result.journal.rbegin(), result.journal.rend(),
                [](const ReferenceEpochDriverRecord &entry) {
                    return entry.action
                        == ReferenceEpochDriverAction::SolveAndEstimate;
                });
            if (solved != result.journal.rend()) {
                const SingularRegionClassification regions =
                    classify_singular_regions_with_physical_radius(
                        *latest_solved_coarse, {Point2(0.0, 0.0)}, solved->ell,
                        config.hybrid_minimum_physical_radius);
                std::vector<int> region_code(
                    latest_solved_coarse->elems.size(), 0);
                for (int element : regions.omega_f_elements)
                    region_code[element] = 1;
                for (int element : regions.omega_s_elements)
                    region_code[element] = 2;
                const std::array<io::VtkIntFieldView, 1> fields{{
                    {"hybrid_region", std::span<const int>(region_code)}}};
                io::VtuDataView data;
                data.cell_int = fields;
                const std::string name = "mesh_E2_final_hybrid_regions.vtu";
                io::write_vtu(
                    directory / name, *latest_solved_coarse, data);
                out << to_string(config.case_id) << ',' << solved->epoch << ','
                    << solved->H_step << ',' << solved->sequence
                    << ",final_solved,hybrid_regions,"
                    << name << ',' << latest_solved_coarse->elems.size() << ','
                    << unknowns(*latest_solved_coarse) << ','
                    << hierarchy.reference_mesh_version() << '\n';
            }
        }
    }
    {
        std::ofstream out(directory / "run.json");
        out << "{\n  \"schema_version\": " << config.schema_version
            << ",\n  \"run_id\": " << json_string(run_id)
            << ",\n  \"status\": "
            << json_string(reference_epoch_driver_state_name(result.state))
            << ",\n  \"stop_reason\": " << json_string(result.stop_reason)
            << ",\n  \"claim\": \"implementation-study\",\n"
            << "  \"algorithm_variant\": {\"name\":"
            << json_string(
                   config.singularity_hybrid
                   ? "adaptive-conformity-collar-v1"
                   : "safeguarded-reference-epoch-v1")
            << ",\"manuscript_conformance\":"
            << json_string(
                   config.singularity_hybrid
                   ? "implementation-erratum"
                   : "implementation-study-variant")
            << ",\"coarse_marking\":"
            << json_string(
                   config.singularity_hybrid
                   ? "full-regular-Doerfler-with-closure-safe-selection-when-feasible"
                   : "full-Doerfler")
            << ",\"reserve_selection\":"
            << json_string(
                   config.singularity_hybrid
                   ? "smallest-spill-free-feasible-collar"
                   : "uniform")
            << ",\"reserve_trigger_scope\":"
            << json_string(
                   config.singularity_hybrid
                   ? "far-full-target-regular-cells"
                   : "not-applicable")
            << ",\"candidate_hierarchy_closure_order\":"
            << json_string("post-estimator")
            << ",\"coarse_admissibility\":"
            << json_string(
                   "not-enforced;pre-asymptotic-points-retained")
            << "},\n  \"algorithm_2_reserve_policy\": "
            << json_string(
                       config.singularity_hybrid
                       ? "graded-interface-reserve:min(g_tar,max(0,graph_distance_to_Omega_F-adaptive_conformity_collar))"
                       : "not-applicable")
            << ",\n  \"timing_semantics\": {"
               "\"time_total_cumulative\":\"raw driver wall time\","
               "\"time_method_cumulative\":\"raw wall time minus completed reference/exact validation and mesh artifact capture\","
               "\"time_hybrid_coarse_marking\":\"full closure-aware selection including all previews\","
               "\"time_hybrid_coarse_preview\":\"all marked NVB/reference-embedding previews; accepted preview is reused by proposal\","
               "\"wall_limit_clock\":\"raw driver wall including validation and artifact capture\","
               "\"validation_changes_estimators_or_marking\":false},\n"
            << "  \"manuscript_sha256\": "
            << json_string(config.manuscript_sha256)
            << ",\n  \"code_commit\": " << json_string(config.git_commit)
            << ",\n  \"build_hash\": " << json_string(config.build_hash)
            << ",\n  \"case_definition\": ";
        if (config.case_id == PaperCase::R1) {
            out << "{\"revision\":\"localized-smooth-revised-v1\","
                   "\"boundary_partition\":\"D-top-bottom,N-left,R-right\","
                   "\"localization_center\":[0.75,0.5],"
                   "\"localization_alpha\":80,"
                   "\"oscillatory_phase\":\"exp(i*kappa*x)\"}";
        } else if (config.case_id == PaperCase::S) {
            out << "{\"revision\":\"l-shape-additive-revised-v1\","
                   "\"boundary_partition\":\"D-reentrant-rays,N-empty,R-remainder\","
                   "\"singular_cutoff\":\"C-infinity-flat-step\","
                   "\"singular_cutoff_inner_radius\":0.25,"
                   "\"singular_cutoff_outer_radius\":"
                << number(config.singular_cutoff_outer_radius)
                << ",\"oscillatory_bump\":\"tensor-C-infinity\","
                   "\"oscillatory_support\":[-0.75,-0.25,0.25,0.75],"
                   "\"oscillatory_amplitude\":"
                << number(config.smooth_wave_amplitude)
                << ",\"singular_oscillatory_fraction\":"
                << number(config.singular_oscillatory_fraction)
                << ",\"oscillatory_phase\":\"exp(i*kappa*x)\","
                   "\"hybrid_physical_ball_geometric_tolerance\":1e-12}";
        } else {
            out << "{\"revision\":\"unchanged-paper-case\"}";
        }
        out
            << ",\n  \"hardware\": {\"compiler\": "
            << json_string(__VERSION__) << ", \"hardware_threads\": "
            << std::thread::hardware_concurrency()
            << ", \"peak_memory_mb\": " << number(peak_memory_mb()) << "}"
            << ",\n  \"config\": " << canonical_json(config) << "\n}\n";
    }
}

} // namespace

int run_reference_epoch_paper(
    const std::string_view config_json,
    const std::filesystem::path &output_directory,
    const std::filesystem::path &manuscript_baseline,
    bool check,
    bool validate_only) {
    const ReferenceEpochPaperConfig config =
        parse_reference_epoch_paper_config(config_json);
    const std::string baseline = first_token(read_text(manuscript_baseline));
    const std::string expected = config.manuscript_sha256.substr(7);
    if (baseline != expected)
        throw std::runtime_error(
            "manuscript baseline does not match schema-v6 config");
    if (validate_only) return 0;
    const std::string run_id = make_run_id(config);
    const std::filesystem::path run_directory = output_directory / run_id;
    std::filesystem::create_directories(run_directory);
    PaperCaseData data = benchmarks::make_paper_case(
        config.case_id, config.wavenumber,
        config.singular_oscillatory_fraction,
        config.singular_cutoff_outer_radius,
        config.singular_quintic_cutoff,
        config.smooth_wave_amplitude);
    NumericalReferenceEpochBackend backend(config, std::move(data));
    ReferenceEpochPracticalDriver driver(
        backend, make_reference_epoch_driver_config(config));
    const ReferenceEpochDriverResult result = driver.run();
    write_iterations(
        run_directory / "iterations.csv", config, run_id, result);
    write_auxiliary_outputs(
        run_directory, config, run_id, result, backend.hierarchy(),
        backend.snapshots(), backend.latest_solved_coarse(),
        backend.hybrid_reserve_diagnostics());
    if (check) {
        for (const char *file : {"iterations.csv", "run.json", "summary.csv",
                 "epoch_history.csv", "mesh_manifest.csv", "corrector_work.csv",
                 "hybrid_reserve.csv"}) {
            if (!std::filesystem::is_regular_file(run_directory / file))
                throw std::runtime_error(std::string("missing WP7 artifact: ") + file);
        }
        if (result.state == ReferenceEpochDriverState::Failed)
            throw std::runtime_error("reference-epoch smoke failed: " + result.stop_reason);
    }
    std::cout << "run_id=" << run_id << '\n'
              << "state=" << reference_epoch_driver_state_name(result.state) << '\n'
              << "stop_reason=" << result.stop_reason << '\n'
              << "output=" << run_directory.string() << '\n';
    return result.state == ReferenceEpochDriverState::Failed ? 1 : 0;
}

} // namespace lod2d::helmholtz::experiments

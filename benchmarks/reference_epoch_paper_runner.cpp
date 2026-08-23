#include "helmholtz/experiments/reference_epoch_runner.h"

#include "helmholtz/adaptive/candidate_dual.h"
#include "helmholtz/adaptive/candidate_flux.h"
#include "helmholtz/adaptive/certificates.h"
#include "helmholtz/adaptive/error_control.h"
#include "helmholtz/adaptive/kernel_residual.h"
#include "helmholtz/adaptive/singularity_hybrid.h"
#include "helmholtz/benchmarks/paper_cases.h"
#include "helmholtz/boundary.h"
#include "helmholtz/experiments/paper_config.h"
#include "helmholtz/model.h"
#include "io/vtk_writer.h"

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
    std::string stage;
    ReferenceEpochDriverAction anchor_action =
        ReferenceEpochDriverAction::BeginEpoch;
    std::size_t anchor_occurrence = 0;
    TriMesh coarse;
    TriMesh candidate;
    std::optional<TriMesh> reference;
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
        if (!included_coarse_elements || local_gap > 0)
            gap = std::min(gap, local_gap);
    }
    return gap;
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
        model_.reset();
        solution_.reset();
        regions_.reset();
        ReferenceEpochMeshSnapshot snapshot;
        snapshot.epoch = epoch_index_;
        snapshot.H_step = committed_H_steps_;
        snapshot.stage = "epoch_start";
        snapshot.anchor_action = ReferenceEpochDriverAction::BeginEpoch;
        snapshot.anchor_occurrence = begin_occurrence_++;
        snapshot.coarse = hierarchy_.coarse_mesh();
        snapshot.candidate = hierarchy_.candidate_mesh();
        snapshot.reference = hierarchy_.reference_mesh();
        snapshots_.push_back(std::move(snapshot));
    }

    ReferenceEpochCorrectorObservation corrector_check(int ell) override {
        const Clock::time_point total_start = Clock::now();
        std::vector<int> skipped;
        if (config_.singularity_hybrid) {
            HybridMatchingResult matching = restore_hybrid_reference_matching(
                hierarchy_, {Point2(0.0, 0.0)}, ell);
            regions_ = matching.regions;
            skipped = matching.regions.omega_s_elements;
        }
        HelmholtzProblemConfig model_config;
        model_config.ell = ell;
        model_config.wavenumber = config_.wavenumber;
        model_config.initial_mesh = data_.initial_mesh;
        model_config.quadrature = config_.quadrature;
        model_config.quadrature_context = data_.quadrature_context;
        if (config_.singularity_hybrid) {
            model_ = std::make_unique<HelmholtzLodModel>(
                HelmholtzLodModel::build_adaptive_hybrid(
                    model_config, hierarchy_, skipped, &corrector_cache_));
        } else {
            model_ = std::make_unique<HelmholtzLodModel>(
                HelmholtzLodModel::build_adaptive(
                    model_config, hierarchy_, &corrector_cache_));
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
        const ReferenceCorrectorCertificate certificate =
            build_reference_corrector_certificate(
                hierarchy_, model_->operators(), model_->corrected_test_basis(),
                free_nodes(hierarchy_.coarse_mesh()),
                KernelRieszSolver::SaddlePoint, eigen, &gram_factor_cache_);
        localization_warm_start_ = certificate.spectrum.dominant_vector;
        ReferenceEpochCorrectorObservation result;
        result.theta_loc = certificate.theta_loc;
        // Compare the frozen implementation-study threshold against the
        // residual-enclosed upper endpoint, not merely the central Ritz value.
        result.delta_loc_hat = certificate.theta_loc_diagnostic_upper;
        result.rebuilt_correctors = static_cast<std::size_t>(
            model_->correctors().diagnostics.patch_count);
        result.skipped_correctors = static_cast<std::size_t>(
            model_->correctors().diagnostics.skipped_patch_count);
        result.skipped_corrector_work_units =
            model_->correctors().diagnostics.skipped_patch_work_units;
        result.time_corrector = model_->build_timings().correctors_ms / 1000.0;
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
        result.gram_parallel_threads =
            certificate.gram_operator.parallel_threads;
        result.localization_iterations = certificate.spectrum.iterations;
        result.localization_relative_residual =
            certificate.spectrum.relative_residual;
        result.localization_used_warm_start =
            certificate.spectrum.used_warm_start;
        (void)total_start;
        return result;
    }

    ReferenceEpochSolveObservation solve_and_estimate() override {
        if (!model_) throw std::logic_error("solve requested before corrector check");
        const ComplexVector load = assemble_helmholtz_load(
            hierarchy_.reference_mesh(), data_.source, config_.quadrature,
            data_.quadrature_context);
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
        result.marked_H = config_.singularity_hybrid
            ? mark_hybrid_regular_region(
                  estimate.element_eta_squared, *regions_, config_.theta_H)
            : estimate.marked_elements;

        ensure_reference_solution(load);
        const HelmholtzError reference_norm = compute_discrete_helmholtz_error(
            hierarchy_.reference_mesh(), model_->operators(),
            *reference_solution_, ComplexVector::Zero(reference_solution_->size()));
        const HelmholtzError reference_error = compute_discrete_helmholtz_error(
            hierarchy_.reference_mesh(), model_->operators(),
            *reference_solution_, solution_->fine_values);
        result.relative_reference_energy = reference_error.energy
            / std::max(reference_norm.energy, std::numeric_limits<double>::min());
        if (data_.exact && data_.exact_gradient) {
            const HelmholtzError exact_error = compute_helmholtz_error(
                hierarchy_.reference_mesh(), solution_->fine_values,
                config_.wavenumber, data_.exact, data_.exact_gradient,
                config_.quadrature, data_.quadrature_context);
            const HelmholtzError exact_norm = compute_helmholtz_error(
                hierarchy_.reference_mesh(),
                ComplexVector::Zero(solution_->fine_values.size()),
                config_.wavenumber, data_.exact, data_.exact_gradient,
                config_.quadrature, data_.quadrature_context);
            result.relative_exact_energy = exact_error.energy / exact_norm.energy;
            result.relative_exact_L2 = exact_error.l2 / exact_norm.l2;
        }
        latest_solved_coarse_ = hierarchy_.coarse_mesh();
        return result;
    }

    void propose_coarse_refinement(const std::vector<int> &marked_H) override {
        const ReferenceEpochRefinementResult proposed =
            hierarchy_.propose_coarse_refinement(marked_H);
        if (!proposed.changed())
            throw std::runtime_error("coarse proposal did not refine the mesh");
    }

    ReferenceEpochCandidateObservation enrich_candidate() override {
        if (!solution_) throw std::logic_error("candidate enrichment has no LOD value");
        const Clock::time_point total_start = Clock::now();
        const Clock::time_point close_start = Clock::now();
        const ReferenceEpochRefinementResult close_result =
            hierarchy_.close_candidate_over_proposed_coarse();
        const double time_close = seconds_since(close_start);
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
        result.marked_c = flux.marked_elements;
        const Clock::time_point enrich_start = Clock::now();
        ReferenceEpochRefinementResult enrich_result;
        if (!flux.marked_elements.empty())
            enrich_result = hierarchy_.enrich_candidate(flux.marked_elements);
        result.time_candidate_enrich = seconds_since(enrich_start);
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
                const SingularRegionClassification proposed_regions =
                    classify_singular_regions(
                        hierarchy_.proposed_coarse_mesh(),
                        {Point2(0.0, 0.0)}, regions_->ell);
                return prospective_reference_level_gap(
                    hierarchy_, &proposed_regions.in_regular);
            }
            return prospective_reference_level_gap(
                hierarchy_, &regions_->in_regular);
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
        return result;
    }

    void commit_coarse_refinement() override {
        ComplexVector prolonged_warm_start;
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
            const RefineOutput transfer = build_nested_mesh_embedding(
                hierarchy_.coarse_mesh(), proposed);
            const ComplexVector proposed_nodal =
                transfer.P_node.cast<Complex>() * old_nodal;
            prolonged_warm_start.resize(proposed_free.size());
            for (int index = 0;
                 index < static_cast<int>(proposed_free.size()); ++index) {
                prolonged_warm_start(index) = proposed_nodal(proposed_free[index]);
            }
        }
        const ReferenceEpochRefinementResult committed =
            hierarchy_.commit_coarse_refinement();
        if (!committed.changed())
            throw std::runtime_error("coarse proposal was not committable");
        localization_warm_start_ = std::move(prolonged_warm_start);
        model_.reset();
        solution_.reset();
        regions_.reset();
        ++committed_H_steps_;
        ReferenceEpochMeshSnapshot snapshot;
        snapshot.epoch = epoch_index_;
        snapshot.H_step = committed_H_steps_;
        snapshot.stage = "committed";
        snapshot.anchor_action =
            ReferenceEpochDriverAction::CommitCoarseRefinement;
        snapshot.anchor_occurrence = commit_occurrence_++;
        snapshot.coarse = hierarchy_.coarse_mesh();
        snapshot.candidate = hierarchy_.candidate_mesh();
        snapshots_.push_back(std::move(snapshot));
    }

    void refresh_reference(const int minimum_post_refresh_level_gap) override {
        hierarchy_.deepen_candidate_over_proposed_coarse(
            minimum_post_refresh_level_gap);
        ReferenceEpochMeshSnapshot snapshot;
        snapshot.epoch = epoch_index_;
        snapshot.H_step = committed_H_steps_;
        snapshot.stage = "pre_switch";
        snapshot.anchor_action = ReferenceEpochDriverAction::RefreshReference;
        snapshot.anchor_occurrence = refresh_occurrence_++;
        snapshot.coarse = hierarchy_.coarse_mesh();
        snapshot.candidate = hierarchy_.candidate_mesh();
        snapshots_.push_back(std::move(snapshot));
        hierarchy_.refresh_reference_from_candidate();
        gram_factor_cache_.clear();
        reference_solution_.reset();
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
        return result;
    }

    const ReferenceEpochHierarchy &hierarchy() const { return hierarchy_; }
    const std::vector<ReferenceEpochMeshSnapshot> &snapshots() const {
        return snapshots_;
    }
    const std::optional<TriMesh> &latest_solved_coarse() const {
        return latest_solved_coarse_;
    }

private:
    void ensure_reference_solution(const ComplexVector &load) {
        if (!reference_solution_)
            reference_solution_ = solve_helmholtz_fem(model_->operators(), load);
    }

    ReferenceEpochPaperConfig config_;
    PaperCaseData data_;
    ReferenceEpochHierarchy hierarchy_;
    HelmholtzCorrectorPatchCache corrector_cache_;
    // Numeric constraints change on most patches after a coarse commit.
    // Keep only a small tail: a 512-entry pilot retained several GiB for
    // very few hits, while 32 bounds the optional cross-step memory cost.
    ReferenceDefectGramFactorCache gram_factor_cache_{32};
    ComplexVector localization_warm_start_;
    std::unique_ptr<HelmholtzLodModel> model_;
    std::optional<HelmholtzLodSolution> solution_;
    std::optional<ComplexVector> reference_solution_;
    std::optional<SingularRegionClassification> regions_;
    std::optional<TriMesh> latest_solved_coarse_;
    bool epoch_started_ = false;
    std::size_t epoch_index_ = 0;
    std::size_t committed_H_steps_ = 0;
    std::size_t begin_occurrence_ = 0;
    std::size_t commit_occurrence_ = 0;
    std::size_t refresh_occurrence_ = 0;
    std::vector<ReferenceEpochMeshSnapshot> snapshots_;
};

void write_iterations(
    const std::filesystem::path &path,
    const ReferenceEpochPaperConfig &config,
    const std::string &run_id,
    const ReferenceEpochDriverResult &result) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "schema_version,manuscript_sha256,code_commit,case,method,kappa,run_id,repeat_index,"
           "epoch,iteration,state,action,stop_reason,N_H,N_h,N_c,ell,kappa_H_max,"
           "eta_H,Theta_loc,U_practical,eta_eq_c,eta_dual_c,L_gap_c,dual_check_performed,"
           "minimum_reference_level_gap,level_gap_dual_trigger,"
           "relative_reference_energy,relative_exact_energy,relative_exact_L2,marked_H,marked_c,"
           "rebuilt_correctors,skipped_correctors,skipped_corrector_work_units,"
           "time_corrector,time_theta,"
           "time_gram_prepare_structure,time_gram_prepare_factorization,"
           "time_gram_action_rhs,time_gram_action_patch_solve,time_gram_action_scatter,"
           "gram_action_calls,gram_patch_factorizations,gram_factor_cache_hits,"
           "gram_factor_cache_misses,gram_parallel_threads,"
           "localization_iterations,localization_relative_residual,"
           "localization_used_warm_start,time_lod_solve,"
           "time_reference_riesz,time_candidate_flux,time_candidate_close,"
           "time_candidate_operator_assembly,time_candidate_prolongation,"
           "time_candidate_flux_reconstruction,time_candidate_flux_prepare,"
           "time_candidate_flux_patch_solve,time_candidate_flux_merge,"
           "time_candidate_flux_audit,candidate_flux_parallel_threads,"
           "time_candidate_enrich,time_candidate_nvb_refine,"
           "time_candidate_embedding_composition,time_candidate_parent_map_update,"
           "time_candidate_quasi_interpolation,time_candidate_embedding_validation,"
           "time_candidate_dual,time_candidate_dual_operator_assembly,"
           "time_candidate_dual_load_assembly,time_candidate_dual_prolongation,"
           "time_candidate_dual_solve,time_candidate_dual_prepare,"
           "time_candidate_dual_patch_solve,time_candidate_dual_reduction,"
           "candidate_dual_patch_factorizations,candidate_dual_parallel_threads,"
           "time_mesh,"
           "time_total_cumulative,peak_memory_mb\n";
    for (const ReferenceEpochDriverRecord &row : result.journal) {
        out << config.schema_version << ',' << config.manuscript_sha256 << ','
            << config.git_commit << ',' << to_string(config.case_id) << ','
            << config.method << ',' << number(config.wavenumber) << ',' << run_id
            << ',' << config.repeat_index << ',' << row.epoch << ',' << row.sequence
            << ',' << reference_epoch_driver_state_name(row.state_after) << ','
            << reference_epoch_driver_action_name(row.action) << ','
            << csv(row.action == ReferenceEpochDriverAction::Fail
                       || row.action == ReferenceEpochDriverAction::StopWorkLimit
                       ? row.detail : std::string{})
            << ',' << row.coarse_unknowns << ',' << row.reference_unknowns << ','
            << row.candidate_unknowns << ',' << row.ell << ','
            << number(row.kappa_H_max) << ',' << number(row.eta_H) << ','
            << number(row.theta_loc) << ',' << number(row.U_practical) << ','
            << number(row.eta_eq_c) << ',' << number(row.eta_dual_c) << ','
            << number(row.L_gap_c) << ','
            << (row.action == ReferenceEpochDriverAction::ComputeCandidateDual
                    ? "true" : "false")
            << ','
            << (row.minimum_reference_level_gap
                        == std::numeric_limits<int>::max()
                    ? "NA"
                    : std::to_string(row.minimum_reference_level_gap))
            << ',' << (row.level_gap_dual_trigger ? "true" : "false")
            << ',' << number(row.relative_reference_energy) << ','
            << number(row.relative_exact_energy) << ','
            << number(row.relative_exact_L2) << ',' << row.marked_H << ','
            << row.marked_c << ',' << row.rebuilt_correctors << ','
            << row.skipped_correctors << ','
            << row.skipped_corrector_work_units << ','
            << number(row.time_corrector) << ','
            << number(row.time_theta) << ','
            << number(row.time_gram_prepare_structure) << ','
            << number(row.time_gram_prepare_factorization) << ','
            << number(row.time_gram_action_rhs) << ','
            << number(row.time_gram_action_patch_solve) << ','
            << number(row.time_gram_action_scatter) << ','
            << row.gram_action_calls << ',' << row.gram_patch_factorizations << ','
            << row.gram_factor_cache_hits << ','
            << row.gram_factor_cache_misses << ','
            << row.gram_parallel_threads << ',' << row.localization_iterations << ','
            << number(row.localization_relative_residual) << ','
            << (row.localization_used_warm_start ? "true" : "false") << ','
            << number(row.time_lod_solve) << ','
            << number(row.time_reference_riesz) << ','
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
            << number(row.time_candidate_dual_operator_assembly) << ','
            << number(row.time_candidate_dual_load_assembly) << ','
            << number(row.time_candidate_dual_prolongation) << ','
            << number(row.time_candidate_dual_solve) << ','
            << number(row.time_candidate_dual_prepare) << ','
            << number(row.time_candidate_dual_patch_solve) << ','
            << number(row.time_candidate_dual_reduction) << ','
            << row.candidate_dual_patch_factorizations << ','
            << row.candidate_dual_parallel_threads << ','
            << number(row.time_mesh)
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
    const std::optional<TriMesh> &latest_solved_coarse) {
    {
        std::ofstream out(directory / "summary.csv");
        out << "schema_version,run_id,row_type,target,achieved,epoch,iteration,N_H,relative_energy,time_total_cumulative,status,stop_reason\n";
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
                    << number(hit->time_total_cumulative);
            } else {
                out << ",NA,NA,NA,NA,NA";
            }
            out << ',' << reference_epoch_driver_state_name(result.state)
                << ',' << csv(result.stop_reason) << '\n';
        }
        out << config.schema_version << ',' << run_id
            << ",run,NA,NA,NA,NA,NA,NA,NA,"
            << reference_epoch_driver_state_name(result.state) << ','
            << csv(result.stop_reason) << '\n';
    }
    {
        std::ofstream out(directory / "epoch_history.csv");
        out << "schema_version,epoch,iteration,action,N_H,N_h,N_c,L_gap_c\n";
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
        out << "schema_version,epoch,iteration,ell,rebuilt,reused,skipped,skipped_work_units,"
               "time_corrector,time_theta,time_gram_prepare_structure,"
               "time_gram_prepare_factorization,time_gram_action_rhs,"
               "time_gram_action_patch_solve,time_gram_action_scatter,"
               "gram_action_calls,gram_patch_factorizations,gram_parallel_threads,"
               "localization_iterations,localization_relative_residual,"
               "localization_used_warm_start\n";
        for (const auto &row : result.journal) {
            if (row.action == ReferenceEpochDriverAction::AcceptCorrector
                || row.action == ReferenceEpochDriverAction::IncreaseGlobalEll) {
                out << config.schema_version << ',' << row.epoch << ',' << row.sequence
                    << ',' << row.ell << ',' << row.rebuilt_correctors << ",0,"
                    << row.skipped_correctors << ',' << row.skipped_corrector_work_units
                    << ',' << number(row.time_corrector) << ','
                    << number(row.time_theta) << ','
                    << number(row.time_gram_prepare_structure) << ','
                    << number(row.time_gram_prepare_factorization) << ','
                    << number(row.time_gram_action_rhs) << ','
                    << number(row.time_gram_action_patch_solve) << ','
                    << number(row.time_gram_action_scatter) << ','
                    << row.gram_action_calls << ','
                    << row.gram_patch_factorizations << ','
                    << row.gram_parallel_threads << ','
                    << row.localization_iterations << ','
                    << number(row.localization_relative_residual) << ','
                    << (row.localization_used_warm_start ? "true" : "false")
                    << '\n';
            }
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
        out << "case,epoch,iteration,stage,mesh_role,filename,N_cells,N_dofs\n";
        const auto anchor_sequence = [&](const ReferenceEpochMeshSnapshot &snapshot) {
            std::size_t occurrence = 0;
            for (const ReferenceEpochDriverRecord &row : result.journal) {
                if (row.action != snapshot.anchor_action) continue;
                if (occurrence++ == snapshot.anchor_occurrence) return row.sequence;
            }
            throw std::runtime_error("mesh snapshot has no matching journal action");
        };
        for (const ReferenceEpochMeshSnapshot &snapshot : snapshots) {
            const bool selected = snapshot.epoch == 0
                || snapshot.stage == "epoch_start"
                || snapshot.stage == "pre_switch";
            if (!selected) continue;
            const std::size_t iteration = anchor_sequence(snapshot);
            const std::string experiment = config.case_id == PaperCase::R1
                ? "E1" : (config.case_id == PaperCase::S
                    ? "E2" : std::string(to_string(config.case_id)));
            const std::string prefix = "mesh_" + experiment + "_e"
                + [&] {
                    std::ostringstream value;
                    value << std::setfill('0') << std::setw(3) << snapshot.epoch;
                    return value.str();
                }();
            if (snapshot.reference) {
                const std::string name = prefix + "_reference.vtu";
                io::write_vtu(directory / name, *snapshot.reference);
                out << to_string(config.case_id) << ',' << snapshot.epoch << ','
                    << iteration << ',' << snapshot.stage << ",reference,"
                    << name << ',' << snapshot.reference->elems.size() << ','
                    << unknowns(*snapshot.reference) << '\n';
            }
            std::ostringstream iteration_text;
            iteration_text << std::setfill('0') << std::setw(3) << iteration;
            const std::string coarse_snapshot = prefix + "_i"
                + iteration_text.str() + "_coarse.vtu";
            const std::string candidate_snapshot = prefix + "_i"
                + iteration_text.str() + "_candidate.vtu";
            io::write_vtu(directory / coarse_snapshot, snapshot.coarse);
            io::write_vtu(directory / candidate_snapshot, snapshot.candidate);
            out << to_string(config.case_id) << ',' << snapshot.epoch << ','
                << iteration << ',' << snapshot.stage << ",coarse,"
                << coarse_snapshot << ',' << snapshot.coarse.elems.size() << ','
                << unknowns(snapshot.coarse) << '\n';
            out << to_string(config.case_id) << ',' << snapshot.epoch << ','
                << iteration << ',' << snapshot.stage << ",candidate,"
                << candidate_snapshot << ',' << snapshot.candidate.elems.size()
                << ',' << unknowns(snapshot.candidate) << '\n';
        }
        const std::size_t iteration = result.journal.empty()
            ? 0 : result.journal.back().sequence;
        const std::size_t epoch = result.journal.empty()
            ? 0 : result.journal.back().epoch;
        const auto row = [&](const char *role, const std::string &name,
                             const TriMesh &mesh) {
            out << to_string(config.case_id) << ',' << epoch << ',' << iteration
                << ",final," << role << ',' << name << ',' << mesh.elems.size()
                << ',' << unknowns(mesh) << '\n';
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
                    classify_singular_regions(
                        *latest_solved_coarse, {Point2(0.0, 0.0)}, solved->ell);
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
                    << solved->sequence << ",final_solved,hybrid_regions,"
                    << name << ',' << latest_solved_coarse->elems.size() << ','
                    << unknowns(*latest_solved_coarse) << '\n';
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
            << "  \"manuscript_sha256\": "
            << json_string(config.manuscript_sha256)
            << ",\n  \"code_commit\": " << json_string(config.git_commit)
            << ",\n  \"build_hash\": " << json_string(config.build_hash)
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
        throw std::runtime_error("manuscript baseline does not match schema-v5 config");
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
        backend.snapshots(), backend.latest_solved_coarse());
    if (check) {
        for (const char *file : {"iterations.csv", "run.json", "summary.csv",
                 "epoch_history.csv", "mesh_manifest.csv", "corrector_work.csv"}) {
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

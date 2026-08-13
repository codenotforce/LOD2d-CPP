#pragma once

#include "helmholtz/adaptive/reference_retraction.h"
#include "helmholtz/model.h"

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace lod2d::helmholtz::adaptive {

enum class PracticalDriverState {
    CoarseAdmissibility,
    LocalizationCheck,
    SolveAndEstimate,
    RefineCoarse,
    Converged,
    TrajectoryComplete,
    ReferenceRefreshRequired,
    WorkLimitReached,
    Failed,
};

enum class PracticalDriverAction {
    InitializeReferenceEpoch,
    RefineCoarseForAdmissibility,
    AcceptCoarse,
    IncreaseGlobalEll,
    AcceptLocalization,
    AcceptFixedEll,
    SolveStandardLod,
    RefineUniformLod,
    SolveUniformFem,
    RefineUniformFem,
    SolveAdaptiveFem,
    RefineAdaptiveFem,
    FormCoarseMarking,
    RefineCoarse,
    CompleteReferenceEpoch,
    RefreshReferenceEpoch,
    Complete,
    CompleteTrajectory,
    StopReferenceRefreshRequired,
    StopWorkLimit,
    Fail,
};

struct PracticalDriverProblem {
    TriMesh initial_mesh;
    ComplexFunction source;
    QuadraturePolicy quadrature;
    QuadratureContext quadrature_context;
};

struct PracticalWorkLimits {
    std::size_t maximum_iterations = 100;
    std::size_t maximum_H_steps = 10;
    std::size_t maximum_unknowns = 1000000;
    std::size_t maximum_coarse_elements = 100000;
    std::size_t maximum_ambient_elements = 1000000;
    double maximum_wall_seconds = 0.0;
};

enum class PracticalLocalizationPolicy {
    AdaptiveGlobalEll,
    FixedGlobalEll,
};

enum class PracticalStopPolicy {
    IndicatorTolerance,
    FixedWorkHorizon,
};

struct PracticalDriverConfig {
    int initial_coarse_level = 1;
    int reference_level = 5;
    std::uint64_t reference_epoch = 0;
    int ell0 = 2;
    int ell_max = 6;
    double wavenumber = 4.0;
    double boundary_beta = 1.0;
    // Schema-v2 compatibility/reporting value. The practical driver does not
    // refine from this prior; pre-asymptotic points are diagnosed a posteriori
    // from the reference-error trajectory.
    double c_H = 0.5;
    double theta_loc = 0.25;
    double C0_usr = 1.0;
    double C1_usr = 1.0;
    double theta_H = 0.5;
    double rho_star = 0.25;
    double tolerance_reference = 1e-2;
    PracticalStopPolicy stop_policy =
        PracticalStopPolicy::IndicatorTolerance;
    PracticalLocalizationPolicy localization_policy =
        PracticalLocalizationPolicy::AdaptiveGlobalEll;
    HelmholtzPetrovMode mode = HelmholtzPetrovMode::TwoSided;
    HelmholtzPatchSolverConfig patch_solver;
    KernelRieszSolver riesz_solver = KernelRieszSolver::SaddlePoint;
    LocalizationEigenConfig localization_eigen;
    PracticalWorkLimits limits;
    // Cumulative accepted H-refinement counts after which the current
    // reference epoch is closed and the ambient mesh is promoted to the next
    // reference mesh.  The coarse mesh and H_steps are inherited unchanged.
    std::vector<std::size_t> reference_refresh_H_steps;
};

struct PracticalIterationRecord {
    std::size_t sequence = 0;
    PracticalDriverState state_before =
        PracticalDriverState::CoarseAdmissibility;
    PracticalDriverState state_after =
        PracticalDriverState::CoarseAdmissibility;
    PracticalDriverAction action =
        PracticalDriverAction::InitializeReferenceEpoch;
    std::uint64_t reference_epoch = 0;
    int ell = 0;
    std::size_t coarse_nodes = 0;
    std::size_t coarse_dofs = 0;
    std::size_t reference_nodes = 0;
    std::size_t ambient_nodes = 0;
    std::size_t coarse_elements = 0;
    std::size_t ambient_elements = 0;
    std::size_t marked_H = 0;
    std::size_t ambient_refined_elements = 0;
    std::size_t rebuilt_correctors = 0;
    double kappa_H_max = 0.0;
    double rho_ambient = 0.0;
    double eta_H = 0.0;
    double theta_loc = 0.0;
    double U_practical = 0.0;
    std::optional<double> reference_energy_error;
    std::optional<double> reference_L2_error;
    // Reporting-only errors against a manufactured solution.  The absolute
    // errors use the weighted Helmholtz energy and L2 norms; the relative
    // values are normalized by the corresponding exact-solution norm.
    std::optional<double> exact_energy_error;
    std::optional<double> exact_L2_error;
    std::optional<double> relative_exact_energy_error;
    std::optional<double> relative_exact_L2_error;
    // Reporting-only candidate exported to WP5 after MARK/STOP has completed.
    // The driver never receives an evaluation reference solution.
    ComplexVector evaluation_candidate;
    double time_mesh_seconds = 0.0;
    double time_corrector_seconds = 0.0;
    double time_certificate_seconds = 0.0;
    double time_solve_seconds = 0.0;
    double time_estimator_seconds = 0.0;
    double time_total_cumulative_seconds = 0.0;
    std::string detail;
};

struct PracticalDriverResult {
    PracticalDriverState state = PracticalDriverState::Failed;
    std::string stop_reason;
    int ell = 0;
    std::size_t H_steps = 0;
    double eta_H = 0.0;
    double theta_loc = 0.0;
    double U_practical = 0.0;
    std::vector<int> final_marked_H;
    std::vector<double> final_element_eta_squared;
    std::vector<PracticalIterationRecord> journal;
};

// Reporting-only, one-way sink. It is invoked after an iteration action and
// state transition have been fixed. No value is returned to the driver, so an
// evaluation reference cannot influence MARK/STOP.
using PracticalEvaluationSink =
    std::function<void(std::size_t, const ReferenceEpochHierarchy &,
                       const ComplexVector &)>;

struct PracticalTargetHit {
    double target = 0.0;
    std::optional<std::size_t> journal_index;
};

enum class PracticalConvergenceRegime {
    InsufficientData,
    PreAsymptotic,
    ObservedStableDecay,
};

struct PracticalConvergenceDiagnostic {
    PracticalConvergenceRegime regime =
        PracticalConvergenceRegime::InsufficientData;
    std::size_t distinct_points = 0;
    std::optional<double> previous_log_slope;
    std::optional<double> last_log_slope;
    std::optional<double> last_error_ratio;
    std::optional<double> last_log_improvement;
    std::optional<double> window_geometric_mean_ratio;
    std::optional<double> window_relative_oscillation;
    bool plateau_observed = false;
};

struct PracticalPlateauDiagnosticConfig {
    // A three-point window is plateau-like when its geometric-mean reduction
    // is small and the whole window remains in a narrow band. Reporting-only;
    // a small intermediate increase is deliberately permitted.
    double minimum_geometric_mean_ratio = 0.9;
    double maximum_relative_oscillation = 0.15;
    std::size_t window_steps = 2;
};

// Post-processing only: extract the first empirical reference-energy hit for
// every target from one already completed trajectory.
std::vector<PracticalTargetHit> extract_practical_target_hits(
    const std::vector<PracticalIterationRecord> &journal,
    const std::vector<double> &targets);

// Post-processing only. Uses existing reference-error samples, so it never
// triggers another solve or influences MARK/STOP. Stable decay means that the
// last two error-vs-DOF log slopes are positive and differ by at most factor 2.
PracticalConvergenceDiagnostic diagnose_practical_convergence_regime(
    const std::vector<PracticalIterationRecord> &journal,
    PracticalPlateauDiagnosticConfig plateau = {});

class PracticalAdaptiveDriver {
public:
    PracticalAdaptiveDriver(
        PracticalDriverProblem problem,
        PracticalDriverConfig config,
        PracticalEvaluationSink evaluation_sink = {});

    PracticalDriverResult run();

    const ReferenceEpochHierarchy &hierarchy() const {
        return *hierarchy_;
    }
    PracticalDriverState state() const { return state_; }

private:
    PracticalDriverProblem problem_;
    PracticalDriverConfig config_;
    std::unique_ptr<ReferenceEpochHierarchy> hierarchy_;
    PracticalDriverState state_ =
        PracticalDriverState::CoarseAdmissibility;
    int ell_ = 0;
    std::size_t H_steps_ = 0;
    std::size_t next_reference_refresh_ = 0;
    std::vector<int> pending_marking_;
    ComplexVector reference_load_;
    std::unique_ptr<HelmholtzLodModel> model_;
    std::unique_ptr<ReferenceLocalizationCertificate> localization_;
    std::unique_ptr<ReferenceResidualRiesz> estimator_;
    ComplexVector localization_warm_start_;
    std::vector<PracticalIterationRecord> journal_;
    std::chrono::steady_clock::time_point start_;
    double eta_H_ = 0.0;
    double theta_loc_ = 0.0;
    double U_practical_ = 0.0;
    PracticalEvaluationSink evaluation_sink_;
    double evaluation_seconds_excluded_ = 0.0;

    void validate_config() const;
    void invalidate_discrete_cache();
    void append_record(PracticalIterationRecord record);
    bool work_limit_exceeded(std::string &reason) const;
};

const char *practical_driver_state_name(PracticalDriverState state);
const char *practical_driver_action_name(PracticalDriverAction action);
const char *practical_convergence_regime_name(PracticalConvergenceRegime regime);

} // namespace lod2d::helmholtz::adaptive

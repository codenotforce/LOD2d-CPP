#pragma once

#include "helmholtz/adaptive/certified_driver.h"
#include "helmholtz/adaptive/practical_driver.h"
#include "helmholtz/adaptive/reference_epoch_driver.h"
#include "helmholtz/quadrature.h"

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace lod2d::helmholtz::experiments {

inline constexpr int paper_schema_version = 1;
inline constexpr int practical_paper_schema_version = 4;
inline constexpr int reference_epoch_paper_schema_version = 6;

enum class PaperCase { R1, R2a, R2b, S };
enum class PaperMethod {
    Calod,
    Hlod,
    SlodPrior,
    SlodMatched,
    Ufem,
    Afem,
    HlodProxy
};

// Frozen run-level outcomes. Resource exhaustion is censored data, never a
// silently dropped or successful run.
enum class PaperRunStatus {
    Success,
    Interrupted,
    CensoredWorkLimit,
    CensoredMemoryLimit,
    CensoredTimeLimit,
    CensoredIterationLimit,
    LinearAlgebraFailure,
    CertificateFailure,
    Unavailable
};

// Every nullable numeric output is paired with one of these states. JSON null
// is used for the value itself; NaN, infinity, and numeric sentinels are banned.
enum class PaperValueStatus {
    Valid,
    NotApplicable,
    NotComputed,
    InvalidDenominator,
    EnclosureFailed
};

// Frozen claim labels from Section 1.2 of the paper plan.  In particular,
// "certified" by itself is intentionally not a valid output status.
enum class PaperCertificateStatus {
    ImplementationStudy,
    Conditional,
    AuditCertified,
    ContinuousCertified,
    EmpiricalReference
};

struct CaseDefinition {
    PaperCase id;
    std::string name;
    std::string domain;
    bool has_exact_solution = false;
    bool has_mixed_boundary = false;
    std::optional<double> gaussian_sigma;
};

struct MethodDefinition {
    PaperMethod id;
    std::string name;
    bool paper_comparator = true;
    bool diagnostic_only = false;
    bool adapts_H = false;
    bool adapts_h = false;
    bool adapts_ell = false;
};

struct TolerancePolicy {
    double linear_relative_residual = 1e-10;
    double eigen_relative_residual = 1e-10;
    double interpolation_right_inverse = 1e-9;
    double prolongation_composition = 1e-10;
};

// Algorithm-defining inputs consumed by NumericalCertifiedBackend.  Keeping
// them in the immutable paper configuration prevents two numerically different
// backend runs from sharing a run id.
struct NumericalBackendPolicy {
    int initial_coarse_level = 0;
    int initial_fine_level = 1;
    int initial_oversampling = 1;
    double boundary_beta = 1.0;
    HelmholtzPetrovMode petrov_mode = HelmholtzPetrovMode::TwoSided;
    HelmholtzPatchSolverConfig patch_solver;
    adaptive::CorrectorCertificateConfig certificate = [] {
        adaptive::CorrectorCertificateConfig value;
        value.q0 = 0.25;
        return value;
    }();
    adaptive::KernelRieszSolver kernel_riesz_solver =
        adaptive::KernelRieszSolver::SaddlePoint;
    double audit_doerfler_theta = 0.5;
    double audit_saturation_factor = 0.05;
    // Content digest of the complete constant registry, not a caller-chosen
    // label.  Reusing a human-readable name with different values must change
    // the immutable run identity.
    // SHA-256 digest of the frozen constant-set artifact.  The WP6 runner
    // must recompute and compare it before constructing the registry.
    std::string certificate_constant_set_hash;
};

struct PaperConfig {
    int schema_version = paper_schema_version;
    PaperCase case_id = PaperCase::R1;
    PaperMethod method_id = PaperMethod::Calod;
    double wavenumber = 8.0;
    // Single source for both driver and numerical-backend coarse Doerfler
    // marking; a second backend copy is intentionally forbidden.
    double theta_H = 0.5;
    adaptive::CertifiedErrorTarget error_target =
        adaptive::CertifiedErrorTarget::AuditSpace;
    adaptive::CertifiedEvidencePolicy evidence_policy =
        adaptive::CertifiedEvidencePolicy::RequireVerified;
    double mu0 = 0.5;
    double q0 = 0.25;
    double tau = 0.5;
    double theta_h = 0.5;
    double rho_aud = 0.05;
    double tolerance = 0.1;
    NumericalBackendPolicy numerical_backend;
    adaptive::CertifiedWorkLimits work_limits;
    std::string hlod_prior_corrector_space_id;
    int hlod_prior_oversampling = -1;
    int timing_repeats = 5;
    int repeat_index = 0;
    std::array<double, 4> relative_energy_targets{{0.1, 0.05, 0.02, 0.01}};
    QuadraturePolicy quadrature;
    TolerancePolicy tolerances;
    std::string git_commit;
    std::string build_hash;
};

const std::vector<CaseDefinition> &paper_case_registry();
const std::vector<MethodDefinition> &paper_method_registry();
const CaseDefinition &case_definition(PaperCase id);
const MethodDefinition &method_definition(PaperMethod id);

std::string_view to_string(PaperCase id);
std::string_view to_string(PaperMethod id);
std::string_view to_string(PaperRunStatus status);
std::string_view to_string(PaperValueStatus status);
std::string_view to_string(PaperCertificateStatus status);
PaperCase parse_paper_case(std::string_view text);
PaperMethod parse_paper_method(std::string_view text);
PaperRunStatus parse_paper_run_status(std::string_view text);
PaperValueStatus parse_paper_value_status(std::string_view text);
PaperCertificateStatus parse_paper_certificate_status(std::string_view text);

// The output schema accepts exactly this set of typed numeric metrics.
const std::vector<std::string_view> &paper_output_metric_registry();
void validate_paper_output_metric(std::string_view name);

void validate_paper_config(const PaperConfig &config);
// CALOD/HLOD consume the exact same values that are hashed in PaperConfig;
// method and theta_H are derived rather than duplicated.
adaptive::CertifiedDriverConfig make_certified_driver_config(
    const PaperConfig &config);
std::string canonical_json(const PaperConfig &config);
PaperConfig parse_paper_config(std::string_view json);

// Stable FNV-1a 64-bit hash of canonical_json(config), written as 16 lowercase hex digits.
std::string canonical_config_hash(const PaperConfig &config);
std::string make_run_id(const PaperConfig &config);

bool operator==(const TolerancePolicy &lhs, const TolerancePolicy &rhs);
bool operator==(const NumericalBackendPolicy &lhs,
                const NumericalBackendPolicy &rhs);
bool operator==(const PaperConfig &lhs, const PaperConfig &rhs);

// WP5 practical-paper contract.  It is deliberately separate from the
// legacy certified v1 contract above: theta_h/q_h and certificate evidence
// fields cannot enter a PALOD run identity through this type.
enum class PracticalPaperMethod { Palod, HlodFixed, Slod, Ufem, Afem };
enum class PracticalTrajectoryPolicy {
    PracticalIndicator,
    FixedWorkHorizon,
};

struct PracticalReferenceAdequacyConfig {
    bool enabled = false;
    int refinement_levels = 1;
    // The reference is adequate when ||u_fine-Iu_ref|| / ||u_fine|| is no
    // larger than this fraction of the terminal reported method error.
    double maximum_terminal_error_fraction = 0.25;
};

struct PracticalPaperConfig {
    int schema_version = practical_paper_schema_version;
    PaperCase case_id = PaperCase::R1;
    PracticalPaperMethod method_id = PracticalPaperMethod::Palod;
    double wavenumber = 16.0;
    // Case S only: u=a[(1-gamma)+gamma exp(ikx)].  The historical case uses
    // gamma=1; gamma=0 isolates the reentrant-corner regularity.
    double singular_oscillatory_fraction = 1.0;
    // Case S only.  Increasing the outer cut-off radius smooths its transition
    // while retaining a flat zero trace and normal derivative at r=1.
    double singular_cutoff_outer_radius = 0.5;
    // Case S only.  False preserves the historical C-infinity cut-off; true
    // uses a C2 quintic transition to suppress finite-range annular curvature.
    bool singular_quintic_cutoff = false;
    // Case S only: coefficient B of an additive boundary-compatible smooth
    // wave psi(x,y) exp(i k x).  Zero preserves the historical case.
    double smooth_wave_amplitude = 0.0;
    // Case S only.  "radial-cutoff" preserves historical data, while
    // "boundary-weight-gaussian" selects the revised manuscript benchmark
    // without an artificial radial transition annulus.
    std::string singular_solution_profile = "radial-cutoff";
    std::string reference_mesh = "uniform-nvb";
    int reference_level = 6;
    std::string ambient_mesh = "reference-shadow";
    std::uint64_t reference_epoch = 0;
    // Optional cumulative H-step schedule for explicit in-run reference
    // refreshes.  Empty preserves the original single-epoch v4 protocol.
    std::vector<std::size_t> reference_refresh_H_steps;
    // Optional local-level-driven PALOD epoch policy.  Both values are zero
    // when disabled.
    int reference_refresh_level_gap = 0;
    int maximum_reference_level = 0;
    // Optional a-priori stop for single-reference-epoch trajectories.  Zero
    // disables the gate; k-robust runs use three remaining NVB levels.
    int minimum_reference_level_gap = 0;
    int initial_coarse_level = 2;
    int ell0 = 2;
    int ell_max = 6;
    double boundary_beta = 1.0;
    double c_H = 0.5;
    double theta_loc = 0.25;
    double C0_usr = 1.0;
    double C1_usr = 1.0;
    double theta_H = 0.5;
    double rho_star = 0.25;
    PracticalTrajectoryPolicy trajectory_policy =
        PracticalTrajectoryPolicy::PracticalIndicator;
    // Absolute U_prac threshold, deliberately separate from the relative
    // evaluation-reference targets below. Ignored by FixedWorkHorizon.
    double practical_stop_tolerance = 1e-2;
    adaptive::PracticalPlateauDiagnosticConfig plateau_diagnostic;
    PracticalReferenceAdequacyConfig reference_adequacy;
    HelmholtzPetrovMode petrov_mode = HelmholtzPetrovMode::TwoSided;
    HelmholtzPatchSolverKind patch_solver_kind =
        HelmholtzPatchSolverKind::DirectSaddle;
    int patch_symbolic_cache_slots = 1;
    bool patch_reuse_identical_factorization = false;
    // Zero uses all available OpenMP workers.  Positive values cap only
    // concurrent high-memory patch solves and are part of the run identity.
    int maximum_patch_threads = 0;
    // SLOD-only policy.  With a direct-saddle base solver, switch to
    // direct-Schur once the unconstrained reference-space DoF reaches this
    // threshold. Zero disables the switch. The selected solver is recorded
    // per trajectory row.
    std::size_t slod_direct_schur_min_reference_dofs = 0;
    // SLOD-only performance policy for manufactured-solution experiments.
    // When true, evaluate the exact errors directly on the synchronized fine
    // mesh and skip the redundant fixed-reference solve/error pass.
    bool manufactured_exact_only_errors = false;
    adaptive::KernelRieszSolver kernel_riesz_solver =
        adaptive::KernelRieszSolver::SaddlePoint;
    adaptive::PracticalWorkLimits work_limits;
    int timing_repeats = 1;
    int repeat_index = 0;
    std::array<double, 4> relative_energy_targets{{0.1, 0.05, 0.02, 0.01}};
    QuadraturePolicy quadrature;
    TolerancePolicy tolerances;
    std::string git_commit;
    std::string build_hash;
    std::string manuscript_sha256;
};

std::string_view to_string(PracticalPaperMethod id);
PracticalPaperMethod parse_practical_paper_method(std::string_view text);
std::string_view to_string(PracticalTrajectoryPolicy policy);
PracticalTrajectoryPolicy parse_practical_trajectory_policy(
    std::string_view text);
// Frozen E1 protocol: c_prior = 1 in ell = ceil(c_prior * log2(kappa)).
int standard_lod_prior_ell(double wavenumber);
void validate_practical_paper_config(const PracticalPaperConfig &config);
adaptive::PracticalDriverConfig make_practical_driver_config(
    const PracticalPaperConfig &config);
std::string canonical_json(const PracticalPaperConfig &config);
PracticalPaperConfig parse_practical_paper_config(std::string_view json);
std::string canonical_config_hash(const PracticalPaperConfig &config);
std::string make_run_id(const PracticalPaperConfig &config);
bool operator==(const PracticalPaperConfig &lhs,
                const PracticalPaperConfig &rhs);

// Revised manuscript production contract.  It deliberately contains no
// ambient mesh, ambient ratio, or retraction field.
struct ReferenceEpochPaperConfig {
    int schema_version = reference_epoch_paper_schema_version;
    PaperCase case_id = PaperCase::R1;
    std::string method = "PALOD-reference-epoch";
    double wavenumber = 16.0;
    // Case S manufactured solution.  The legacy radial-cutoff profile is
    // retained for old runs; the revised paper uses boundary-weight-gaussian.
    double singular_oscillatory_fraction = 1.0;
    double singular_cutoff_outer_radius = 0.5;
    bool singular_quintic_cutoff = false;
    double smooth_wave_amplitude = 0.0;
    std::string singular_solution_profile = "radial-cutoff";
    int initial_coarse_level = 3;
    int initial_reference_level = 5;
    bool singularity_hybrid = false;
    // Implementation-study option for the standard (fixed-reference within
    // each epoch) PALOD.  It changes only candidate marking: Omega_F and its
    // complement satisfy independent Doerfler bulk constraints.  It does not
    // enable moving-reference promotion or coarse/reference equality.
    bool candidate_split_regional_marking = false;
    double candidate_regional_minimum_physical_radius = 0.0;
    // Case-S hybrid only.  At every matching check ell_S is the smallest
    // graph radius such that B_R(S) lies in Omega_S=N^{ell_S}(S); the exact
    // matching region is Omega_F=N^{ell_S+ell}(S).  Zero is required for
    // non-hybrid runs.
    double hybrid_minimum_physical_radius = 0.0;
    // Abort before corrector assembly when an active coarse patch expands to
    // more reference elements than this guard.  This is a reproducible pilot
    // safety limit; zero disables it and is required for non-hybrid runs.
    std::size_t hybrid_maximum_corrector_patch_fine_elements = 0;
    int ell0 = 1;
    int ell_max = 4;
    double theta_loc_usr = 0.5;
    double localization_eigen_relative_tolerance = 1e-7;
    int localization_eigen_maximum_iterations = 300;
    double C_rel_usr = 1.0;
    double theta_H = 0.5;
    double theta_c = 0.5;
    std::size_t candidate_update_stride = 1;
    int candidate_force_level_gap = 0;
    bool candidate_closure_cost_aware_marking = false;
    std::size_t candidate_closure_cost_pool_factor = 2;
    // Corrector patch solver policy is part of the reference-epoch run
    // identity.  DirectSchur with identical-support reuse is intended for
    // standard PALOD once reference patches become large; moving PALOD keeps
    // the direct saddle default unless a pilot demonstrates a benefit.
    HelmholtzPatchSolverKind patch_solver_kind =
        HelmholtzPatchSolverKind::DirectSaddle;
    int patch_symbolic_cache_slots = 1;
    bool patch_reuse_identical_factorization = false;
    int maximum_patch_threads = 0;
    // Global reference validation solve. UMFPACK is optional at build time;
    // sparse_lu preserves the historical Eigen backend.
    HelmholtzFemSolverKind reference_solver_kind =
        HelmholtzFemSolverKind::SparseLu;
    // Moving-reference implementation-study policy. One validates every
    // promoted reference; larger values retain exact-error observations but
    // perform the expensive reference solve only at periodic checkpoints.
    std::size_t reference_validation_stride = 1;
    double q_dual = 0.5;
    std::size_t m_dual = 3;
    double tau_ep = 0.5;
    int reference_refresh_level_gap = 0;
    int reference_refresh_target_gap = 0;
    std::size_t minimum_H_steps_per_epoch = 0;
    std::size_t minimum_solved_points_per_new_epoch = 0;
    double tolerance_reference = 1e-2;
    double continuity_constant = 8.0;
    double overlap_constant = 8.0;
    adaptive::ReferenceEpochDriverLimits work_limits;
    int repeat_index = 0;
    QuadraturePolicy quadrature;
    std::string git_commit;
    std::string build_hash;
    std::string manuscript_sha256;
};

void validate_reference_epoch_paper_config(
    const ReferenceEpochPaperConfig &config);
adaptive::ReferenceEpochDriverConfig make_reference_epoch_driver_config(
    const ReferenceEpochPaperConfig &config);
std::string canonical_json(const ReferenceEpochPaperConfig &config);
ReferenceEpochPaperConfig parse_reference_epoch_paper_config(
    std::string_view json);
std::string canonical_config_hash(const ReferenceEpochPaperConfig &config);
std::string make_run_id(const ReferenceEpochPaperConfig &config);
bool is_reference_epoch_paper_config(std::string_view json);

} // namespace lod2d::helmholtz::experiments

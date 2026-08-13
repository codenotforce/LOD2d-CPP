#pragma once

#include "helmholtz/adaptive/certified_driver.h"
#include "helmholtz/adaptive/practical_driver.h"
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
    std::string reference_mesh = "uniform-nvb";
    int reference_level = 6;
    std::string ambient_mesh = "reference-shadow";
    std::uint64_t reference_epoch = 0;
    // Optional cumulative H-step schedule for explicit in-run reference
    // refreshes.  Empty preserves the original single-epoch v4 protocol.
    std::vector<std::size_t> reference_refresh_H_steps;
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

} // namespace lod2d::helmholtz::experiments

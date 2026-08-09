#include "helmholtz/adaptive/numerical_backend.h"

#include "helmholtz/operators.h"

#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <utility>

namespace lod2d::helmholtz::adaptive {
namespace {

class FingerprintBuilder {
public:
    void add_string(std::string_view value) {
        add_u64(static_cast<std::uint64_t>(value.size()));
        for (unsigned char byte : value) {
            hash_ ^= byte;
            hash_ *= 1099511628211ULL;
        }
    }

    void add_u64(std::uint64_t value) {
        for (int byte = 0; byte < 8; ++byte) {
            hash_ ^= static_cast<unsigned char>(
                (value >> (8 * byte)) & 0xffU);
            hash_ *= 1099511628211ULL;
        }
    }

    void add_i64(std::int64_t value) {
        add_u64(std::bit_cast<std::uint64_t>(value));
    }

    void add_double(double value) {
        if (value == 0.0) value = 0.0;
        add_u64(std::bit_cast<std::uint64_t>(value));
    }

    void add_bool(bool value) { add_u64(value ? 1U : 0U); }

    std::string finish(std::string_view prefix) const {
        std::ostringstream encoded;
        encoded << prefix << ":fnv1a64:" << std::hex << std::setfill('0')
                << std::setw(16) << hash_;
        return encoded.str();
    }

private:
    std::uint64_t hash_ = 14695981039346656037ULL;
};

void add_mesh(FingerprintBuilder &builder, const TriMesh &mesh) {
    builder.add_u64(mesh.nodes.size());
    for (const Point2 &node : mesh.nodes) {
        builder.add_double(node.x());
        builder.add_double(node.y());
    }
    builder.add_u64(mesh.elems.size());
    for (const Triangle &element : mesh.elems) {
        for (int node : element) builder.add_i64(node);
    }
    builder.add_u64(mesh.dirichlet.size());
    for (int node : mesh.dirichlet) builder.add_i64(node);
    builder.add_u64(mesh.boundary_edges.size());
    for (const BoundaryEdge &edge : mesh.boundary_edges) {
        builder.add_i64(edge.nodes[0]);
        builder.add_i64(edge.nodes[1]);
        builder.add_u64(static_cast<std::uint64_t>(edge.tag));
    }
}

void add_int_vector(FingerprintBuilder &builder,
                    const std::vector<int> &values) {
    builder.add_u64(values.size());
    for (int value : values) builder.add_i64(value);
}

void add_u64_vector(FingerprintBuilder &builder,
                    const std::vector<std::uint64_t> &values) {
    builder.add_u64(values.size());
    for (std::uint64_t value : values) builder.add_u64(value);
}

void add_quadrature(FingerprintBuilder &builder,
                    const QuadraturePolicy &policy,
                    const QuadratureContext &context) {
    builder.add_i64(policy.base_triangle_order);
    builder.add_i64(policy.gaussian_triangle_order);
    builder.add_i64(policy.singular_triangle_order);
    builder.add_i64(policy.max_recursive_subdivisions);
    builder.add_u64(static_cast<std::uint64_t>(context.integrand_class));
    builder.add_double(context.feature_point.x());
    builder.add_double(context.feature_point.y());
    builder.add_double(context.feature_scale);
}

void add_patch_solver(FingerprintBuilder &builder,
                      const HelmholtzPatchSolverConfig &config) {
    builder.add_u64(static_cast<std::uint64_t>(config.kind));
    builder.add_i64(config.symbolic_cache_slots);
    builder.add_bool(config.reuse_identical_factorization);
    builder.add_i64(config.gmres.restart);
    builder.add_i64(config.gmres.max_iterations);
    builder.add_double(config.gmres.relative_tolerance);
    builder.add_double(config.gmres.absolute_tolerance);
    builder.add_bool(config.gmres.reorthogonalize);
    builder.add_u64(static_cast<std::uint64_t>(config.shifted.rule));
    builder.add_double(config.shifted.alpha);
    builder.add_double(config.shifted.absolute_epsilon);
    builder.add_u64(static_cast<std::uint64_t>(config.shifted.inverse));
    builder.add_i64(config.shifted.pre_smooth);
    builder.add_i64(config.shifted.post_smooth);
    builder.add_i64(config.shifted.coarse_max_dofs);
    builder.add_double(config.shifted.jacobi_weight);
    builder.add_bool(config.fallback_to_direct);
}

void add_constants(FingerprintBuilder &builder,
                   const CertificateConstantRegistry &constants) {
    builder.add_u64(constants.entries().size());
    for (const auto &[name, constant] : constants.entries()) {
        builder.add_string(name);
        builder.add_double(constant.value);
        builder.add_u64(static_cast<std::uint64_t>(constant.direction));
        builder.add_string(constant.source);
        builder.add_string(constant.derivation);
        builder.add_string(constant.mesh_class);
        builder.add_string(constant.patch_policy_hash);
        builder.add_bool(constant.verified);
    }
}

CertificateConstantRegistry floating_constants(
    const CertificateConstantRegistry &supplied) {
    CertificateConstantRegistry result;
    for (const auto &[name, entry] : supplied.entries()) {
        CertificateConstant constant = entry;
        constant.verified = false;
        constant.mesh_fingerprint.clear();
        constant.pde_fingerprint.clear();
        constant.operator_fingerprint.clear();
        if (constant.source.empty())
            constant.source = "caller-supplied floating-point constant";
        if (constant.derivation.empty())
            constant.derivation = "conditional numerical backend input";
        if (constant.mesh_class.empty())
            constant.mesh_class = "runtime adaptive hierarchy";
        result.set(std::move(constant));
    }
    return result;
}

void validate_config(const NumericalCertifiedBackendConfig &config) {
    if (config.problem_id.empty())
        throw std::invalid_argument("numerical backend problem_id is empty");
    if (config.source_id.empty())
        throw std::invalid_argument("numerical backend source_id is empty");
    if (!config.source)
        throw std::invalid_argument("numerical backend source is empty");
    if (config.initial_mesh.nodes.empty() || config.initial_mesh.elems.empty())
        throw std::invalid_argument("numerical backend initial mesh is empty");
    if (config.initial_coarse_level < 0
        || config.initial_fine_level < config.initial_coarse_level
        || config.initial_oversampling < 0) {
        throw std::invalid_argument(
            "numerical backend requires 0 <= H <= h and ell >= 0");
    }
    if (!(std::isfinite(config.wavenumber) && config.wavenumber > 0.0)
        || !(std::isfinite(config.boundary_beta)
             && config.boundary_beta >= 0.0)) {
        throw std::invalid_argument("numerical backend PDE data are invalid");
    }
    const auto valid_theta = [](double theta) {
        return std::isfinite(theta) && theta > 0.0 && theta <= 1.0;
    };
    if (!valid_theta(config.coarse_doerfler_theta)
        || !valid_theta(config.audit_doerfler_theta)) {
        throw std::invalid_argument(
            "numerical backend Doerfler parameters must lie in (0,1]");
    }
    if (!(std::isfinite(config.audit_saturation_factor)
          && config.audit_saturation_factor >= 0.0
          && config.audit_saturation_factor < 1.0)) {
        throw std::invalid_argument(
            "numerical backend saturation factor must lie in [0,1)");
    }
    validate_quadrature_policy(config.quadrature);
}

template <class Observation>
void force_conditional(Observation &observation,
                       std::string source,
                       const std::string &state_fingerprint) {
    observation.evidence.level = CertifiedEvidenceLevel::Conditional;
    observation.evidence.source = std::move(source);
    observation.evidence.hash = state_fingerprint;
}

} // namespace

struct NumericalCertifiedBackend::Impl {
    using Clock = std::chrono::steady_clock;

    explicit Impl(NumericalCertifiedBackendConfig supplied)
        : config(std::move(supplied)),
          constants(floating_constants(config.constants)),
          hierarchy(std::make_unique<AdaptiveMeshHierarchy>(
              config.initial_mesh,
              config.initial_coarse_level,
              config.initial_fine_level)),
          oversampling(config.initial_oversampling),
          start(Clock::now()) {
        validate_config(config);
        initial_problem_id = immutable_problem_fingerprint();
        full_rebuild();
    }

    std::string immutable_problem_fingerprint() const {
        FingerprintBuilder builder;
        builder.add_string("lod2d-numerical-certified-problem-v1");
        builder.add_string(config.problem_id);
        builder.add_string(config.source_id);
        add_mesh(builder, config.initial_mesh);
        builder.add_i64(config.initial_coarse_level);
        builder.add_i64(config.initial_fine_level);
        builder.add_i64(config.initial_oversampling);
        builder.add_double(config.wavenumber);
        builder.add_double(config.boundary_beta);
        builder.add_u64(static_cast<std::uint64_t>(config.petrov_mode));
        add_quadrature(builder, config.quadrature, config.quadrature_context);
        add_patch_solver(builder, config.patch_solver);
        add_constants(builder, constants);
        builder.add_i64(config.certificate.precision_bits);
        builder.add_double(config.certificate.cluster_relative_gap);
        builder.add_double(config.certificate.cluster_absolute_gap);
        builder.add_double(config.certificate.conjugation_tolerance);
        builder.add_double(config.certificate.q0);
        builder.add_u64(static_cast<std::uint64_t>(
            config.kernel_riesz_solver));
        builder.add_double(config.coarse_doerfler_theta);
        builder.add_double(config.audit_doerfler_theta);
        builder.add_double(config.audit_saturation_factor);
        return builder.finish("problem");
    }

    std::string state_fingerprint() const {
        FingerprintBuilder builder;
        builder.add_string("lod2d-numerical-certified-state-v1");
        builder.add_string(initial_problem_id);
        builder.add_i64(oversampling);
        add_mesh(builder, hierarchy->coarse_mesh());
        add_mesh(builder, hierarchy->fine_mesh());
        add_mesh(builder, hierarchy->cert_audit_mesh());
        add_int_vector(builder, hierarchy->coarse_levels());
        add_int_vector(builder, hierarchy->fine_element_levels());
        add_int_vector(builder, hierarchy->cert_audit_element_levels());
        add_u64_vector(builder, hierarchy->coarse_element_ids());
        add_u64_vector(builder, hierarchy->coarse_parent_ids());
        builder.add_u64(hierarchy->coarse_mesh_version());
        builder.add_u64(hierarchy->fine_mesh_version());
        builder.add_u64(hierarchy->cert_audit_mesh_version());
        builder.add_u64(hierarchy->interpolation_version());
        builder.add_u64(hierarchy->boundary_version());
        builder.add_u64(hierarchy->corrector_space_version());
        // The WP4 context binds the assembled fine/audit operators, all
        // prolongations, the PDE configuration, patch policy, and the actual
        // corrector entries.  Mesh/version-only fingerprints would miss a
        // changed model or stale corrector build on the same hierarchy.
        builder.add_string(context.mesh);
        builder.add_string(context.pde);
        builder.add_string(context.patch_policy);
        builder.add_string(context.operators);
        return builder.finish("backend");
    }

    std::string corrector_space_id() const {
        FingerprintBuilder builder;
        builder.add_string("lod2d-corrector-fine-space-v1");
        builder.add_string(initial_problem_id);
        add_mesh(builder, hierarchy->fine_mesh());
        add_int_vector(builder, hierarchy->fine_element_levels());
        return builder.finish("corrector-space");
    }

    void clear_observation_cache() {
        kernel_estimate.reset();
        certificate.reset();
        candidate_on_audit.resize(0);
        audit_load.resize(0);
    }

    void update_peak_memory() {
        const auto saturated_add = [](std::uint64_t first,
                                      std::uint64_t second) {
            return second > std::numeric_limits<std::uint64_t>::max() - first
                ? std::numeric_limits<std::uint64_t>::max()
                : first + second;
        };
        const auto saturated_multiply = [](std::uint64_t first,
                                           std::uint64_t second) {
            return first != 0
                    && second
                        > std::numeric_limits<std::uint64_t>::max() / first
                ? std::numeric_limits<std::uint64_t>::max()
                : first * second;
        };
        const std::uint64_t coarse_nodes =
            hierarchy->coarse_mesh().nodes.size();
        const std::uint64_t fine_nodes = hierarchy->fine_mesh().nodes.size();
        const std::uint64_t audit_nodes =
            hierarchy->cert_audit_mesh().nodes.size();
        const std::uint64_t coarse_elements =
            hierarchy->coarse_mesh().elems.size();
        const std::uint64_t total_nodes = saturated_add(
            coarse_nodes, saturated_add(fine_nodes, audit_nodes));

        // Sparse factorizations can fill densely, and WP4 retains several
        // dense midpoint/radius Grams plus per-coarse-element enclosures.
        // Use a deliberately conservative dimension-based upper budget rather
        // than under-reporting only visible sparse nonzeros.  The 64 MiB floor
        // covers allocator/factorization metadata on tiny runs; saturating
        // arithmetic fails closed on impractically large dimensions.
        std::uint64_t upper = 64ULL * 1024ULL * 1024ULL;
        upper = saturated_add(
            upper,
            saturated_multiply(
                saturated_multiply(total_nodes, total_nodes), 512));
        upper = saturated_add(
            upper,
            saturated_multiply(
                saturated_multiply(
                    saturated_multiply(
                        coarse_elements, coarse_nodes),
                    coarse_nodes),
                256));
        upper = saturated_add(
            upper,
            saturated_multiply(
                saturated_multiply(coarse_elements, fine_nodes),
                3 * sizeof(HelmholtzCorrectorEntry)));
        upper = saturated_add(
            upper,
            saturated_multiply(
                hierarchy->coarse_mesh().elems.size()
                    + hierarchy->fine_mesh().elems.size()
                    + hierarchy->cert_audit_mesh().elems.size(),
                8 * sizeof(Triangle)));
        peak_memory_bytes = std::max(peak_memory_bytes, upper);
    }

    void full_rebuild() {
        HelmholtzProblemConfig model_config;
        model_config.H = config.initial_coarse_level;
        model_config.h = config.initial_fine_level;
        model_config.ell = oversampling;
        model_config.wavenumber = config.wavenumber;
        model_config.boundary_beta = config.boundary_beta;
        model_config.mode = config.petrov_mode;
        model_config.progress = false;
        model_config.initial_mesh = config.initial_mesh;
        model_config.patch_solver = config.patch_solver;
        model_config.quadrature = config.quadrature;
        model_config.quadrature_context = config.quadrature_context;

        model = std::make_unique<HelmholtzLodModel>(
            HelmholtzLodModel::build_adaptive(
                model_config,
                hierarchy->coarse_mesh(),
                hierarchy->coarse_levels(),
                hierarchy->fine_mesh(),
                hierarchy->fine_element_levels()));
        audit = std::make_unique<CertificationAuditService>(
            hierarchy->cert_audit_mesh(), config.wavenumber,
            std::vector<double>{}, std::vector<double>{},
            config.boundary_beta);
        const KernelPatchPolicy policy = audit_kernel_patch_policy(*hierarchy);
        context = certificate_context_fingerprint(
            *hierarchy, *model, audit->operators(), policy);
        clear_observation_cache();
        ++rebuild_count;
        ++work_units;
        update_peak_memory();
    }

    CorrectorCertificateResult build_certificate(
        double eta_H,
        const AuditKernelResidualEvidence &eta_H_evidence) const {
        // Deliberately empty: ordinary FE matrices and floating solves have no
        // directed-rounding assembly enclosure.  WP4 records the missing
        // evidence and returns Conditional.
        const CertificateAssemblyEvidence no_verified_assembly_evidence;
        return build_corrector_certificates(
            *hierarchy, *model, audit->operators(), constants, eta_H,
            eta_H_evidence,
            no_verified_assembly_evidence, config.certificate);
    }

    NumericalCertifiedBackendConfig config;
    CertificateConstantRegistry constants;
    std::unique_ptr<AdaptiveMeshHierarchy> hierarchy;
    std::unique_ptr<HelmholtzLodModel> model;
    std::unique_ptr<CertificationAuditService> audit;
    CertificateContextFingerprint context;
    int oversampling = 0;
    Clock::time_point start;
    std::string initial_problem_id;

    std::optional<AuditKernelResidualEstimate> kernel_estimate;
    std::optional<CorrectorCertificateResult> certificate;
    ComplexVector candidate_on_audit;
    ComplexVector audit_load;

    std::uint64_t work_units = 0;
    std::uint64_t peak_memory_bytes = 0;
    std::size_t rebuild_count = 0;
    std::size_t observations = 0;
};

NumericalCertifiedBackend::NumericalCertifiedBackend(
    NumericalCertifiedBackendConfig config)
    : impl_(std::make_unique<Impl>(std::move(config))) {}

NumericalCertifiedBackend::~NumericalCertifiedBackend() = default;
NumericalCertifiedBackend::NumericalCertifiedBackend(
    NumericalCertifiedBackend &&) noexcept = default;
NumericalCertifiedBackend &NumericalCertifiedBackend::operator=(
    NumericalCertifiedBackend &&) noexcept = default;

CertifiedWorkSnapshot NumericalCertifiedBackend::work_snapshot() const {
    CertifiedWorkSnapshot result;
    result.coarse_dofs = impl_->hierarchy->coarse_mesh().nodes.size();
    result.fine_dofs = impl_->hierarchy->fine_mesh().nodes.size();
    result.audit_dofs = impl_->hierarchy->cert_audit_mesh().nodes.size();
    result.backend_work_units = impl_->work_units;
    result.peak_memory_bytes = impl_->peak_memory_bytes;
    result.elapsed_seconds = std::chrono::duration<double>(
        Impl::Clock::now() - impl_->start).count();
    result.initial_problem_id = impl_->initial_problem_id;
    result.state_fingerprint = impl_->state_fingerprint();
    result.corrector_space_id = impl_->corrector_space_id();
    result.oversampling = impl_->oversampling;
    return result;
}

CoarseAdmissibilityObservation
NumericalCertifiedBackend::inspect_coarse_admissibility() {
    ++impl_->observations;
    ++impl_->work_units;
    CoarseAdmissibilityObservation result =
        make_coarse_admissibility_observation(
            *impl_->hierarchy, *impl_->model, impl_->audit->operators(),
            impl_->constants);
    force_conditional(
        result, "floating-point coarse admissibility bound",
        impl_->state_fingerprint());
    return result;
}

void NumericalCertifiedBackend::refine_coarse(
    const std::vector<int> &marked_elements) {
    if (marked_elements.empty())
        throw std::invalid_argument("coarse refinement marking is empty");
    impl_->hierarchy->refine(marked_elements);
    impl_->full_rebuild();
}

CorrectorCertificationObservation
NumericalCertifiedBackend::inspect_corrector_certification() {
    ++impl_->observations;
    ++impl_->work_units;
    // Step 2 is independent of eta_H; the overall certificate remains
    // conditional until Step 3 supplies estimator-owned provenance.
    const AuditKernelResidualEvidence no_eta_H_evidence;
    impl_->certificate = impl_->build_certificate(
        0.0, no_eta_H_evidence);
    CorrectorCertificationObservation result =
        make_corrector_certification_observation(*impl_->certificate);
    force_conditional(
        result, "floating-point WP4 corrector diagnostics",
        impl_->state_fingerprint());
    return result;
}

void NumericalCertifiedBackend::refine_corrector_patches(
    const std::vector<int> &marked_coarse_sources) {
    if (marked_coarse_sources.empty())
        throw std::invalid_argument("corrector refinement marking is empty");
    impl_->hierarchy->refine_fine_in_coarse_patch(marked_coarse_sources);
    impl_->full_rebuild();
}

void NumericalCertifiedBackend::increase_oversampling() {
    if (impl_->oversampling == std::numeric_limits<int>::max())
        throw std::overflow_error("oversampling level overflowed");
    ++impl_->oversampling;
    impl_->full_rebuild();
}

CoarseErrorObservation
NumericalCertifiedBackend::solve_and_estimate_coarse_error() {
    ++impl_->observations;
    impl_->work_units += 3;
    const HelmholtzLodSolution solution =
        impl_->model->solve_source(impl_->config.source);
    impl_->candidate_on_audit =
        impl_->hierarchy->fine_to_cert_audit().cast<Complex>()
        * solution.fine_values;
    impl_->audit_load = assemble_helmholtz_load(
        impl_->hierarchy->cert_audit_mesh(), impl_->config.source,
        impl_->config.quadrature, impl_->config.quadrature_context);
    impl_->kernel_estimate = estimate_audit_kernel_residual(
        *impl_->hierarchy, impl_->audit->operators(), impl_->audit_load,
        impl_->candidate_on_audit, impl_->config.coarse_doerfler_theta,
        impl_->config.kernel_riesz_solver);
    impl_->certificate = impl_->build_certificate(
        impl_->kernel_estimate->eta,
        impl_->kernel_estimate->evidence);
    CoarseErrorObservation result = make_coarse_error_observation(
        *impl_->kernel_estimate, *impl_->certificate);
    std::string evidence_source = "floating-point WP3/WP4 coarse error diagnostics";
    if (!impl_->kernel_estimate->evidence.backend().empty())
        evidence_source += ":" + impl_->kernel_estimate->evidence.backend();
    std::string evidence_hash = impl_->state_fingerprint();
    if (!impl_->kernel_estimate->evidence.diagnostic_fingerprint().empty())
        evidence_hash += "|"
            + impl_->kernel_estimate->evidence.diagnostic_fingerprint();
    force_conditional(result, std::move(evidence_source), evidence_hash);
    impl_->update_peak_memory();
    return result;
}

AuditControlObservation NumericalCertifiedBackend::inspect_audit_control() {
    // This observation is a deterministic solve on the current certification
    // audit space.  It intentionally has no hidden dependency on the cached
    // Step-3 LOD solution, so a checkpoint resumed directly at AuditControl
    // can reproduce it after mutation replay.
    ++impl_->observations;
    impl_->work_units += 2;
    const CertificationAuditSolution audit_solution =
        impl_->audit->solve_source(
            impl_->config.source, impl_->config.quadrature,
            impl_->config.quadrature_context);
    const EmpiricalSaturationAuditEstimate estimate =
        estimate_empirical_saturation_audit_error(
            impl_->hierarchy->cert_audit_mesh(), audit_solution.values,
            impl_->config.source, impl_->config.wavenumber,
            impl_->hierarchy->cert_audit_parent_fine_elements(),
            static_cast<int>(impl_->hierarchy->fine_mesh().elems.size()),
            impl_->config.audit_saturation_factor,
            impl_->config.audit_doerfler_theta,
            std::vector<double>{}, std::vector<double>{},
            impl_->config.boundary_beta, impl_->config.quadrature,
            impl_->config.quadrature_context);

    AuditControlObservation result;
    result.evidence.valid = std::isfinite(estimate.audit_error_lower)
        && std::isfinite(estimate.audit_error_upper)
        && estimate.audit_error_lower >= 0.0
        && estimate.audit_error_upper >= estimate.audit_error_lower;
    if (!result.evidence.valid)
        result.evidence.invalid_reason =
            "empirical saturation audit interval is invalid";
    result.interval_kind = AuditIntervalKind::EmpiricalSaturation;
    result.audit_error_lower = estimate.audit_error_lower;
    result.audit_error_upper = estimate.audit_error_upper;
    result.marked_fine_elements = estimate.marked_fine_elements;
    result.refinement_available = true;
    force_conditional(
        result,
        "empirical two-level saturation audit (assumption not verified)",
        impl_->state_fingerprint());
    return result;
}

void NumericalCertifiedBackend::refine_audit(
    const std::vector<int> &marked_fine_elements) {
    if (marked_fine_elements.empty())
        throw std::invalid_argument("audit refinement marking is empty");
    impl_->hierarchy->refine_cert_audit_from_fine_elements(
        marked_fine_elements);
    impl_->full_rebuild();
}

std::size_t NumericalCertifiedBackend::full_rebuild_count() const {
    return impl_->rebuild_count;
}

std::size_t NumericalCertifiedBackend::observation_count() const {
    return impl_->observations;
}

const AdaptiveMeshHierarchy &NumericalCertifiedBackend::hierarchy() const {
    return *impl_->hierarchy;
}

const CertificateContextFingerprint &
NumericalCertifiedBackend::certificate_context() const {
    return impl_->context;
}

} // namespace lod2d::helmholtz::adaptive

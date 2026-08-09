#pragma once

#include "helmholtz/operators.h"

#include <cstddef>
#include <functional>
#include <memory>
#include <vector>

namespace lod2d::helmholtz::adaptive {

enum class ErrorReferenceRole {
    Exact,
    ExternalEstimator,
    TwoLevelSaturation
};

using ExternalErrorEstimator =
    std::function<HelmholtzError(const TriMesh &, const ComplexVector &)>;

// Explicit role object used by diagnostics and validation. The reference data
// stay private, so callers cannot reinterpret an exact or two-level quantity
// as another role without constructing a different object.
class ErrorReference {
public:
    static ErrorReference exact(
        ComplexFunction exact_value,
        ComplexGradientFunction exact_gradient,
        QuadraturePolicy quadrature = {},
        QuadratureContext quadrature_context = {});
    static ErrorReference external(ExternalErrorEstimator estimator);
    static ErrorReference two_level(ComplexVector finer_values);

    ErrorReferenceRole role() const { return role_; }
    HelmholtzError evaluate(
        const TriMesh &mesh,
        const HelmholtzOperators &operators,
        const ComplexVector &candidate) const;

private:
    ErrorReference() = default;

    ErrorReferenceRole role_ = ErrorReferenceRole::Exact;
    ComplexFunction exact_value_;
    ComplexGradientFunction exact_gradient_;
    QuadraturePolicy quadrature_;
    QuadratureContext quadrature_context_;
    ExternalErrorEstimator external_estimator_;
    ComplexVector finer_values_;
};

HelmholtzError compute_discrete_helmholtz_error(
    const TriMesh &mesh,
    const HelmholtzOperators &operators,
    const ComplexVector &reference,
    const ComplexVector &approximation);

HelmholtzError compute_local_discrete_helmholtz_error(
    const TriMesh &mesh,
    const HelmholtzOperators &operators,
    const ComplexVector &reference,
    const ComplexVector &approximation,
    const std::vector<int> &elements);

HelmholtzError compute_local_exact_helmholtz_error(
    const TriMesh &mesh,
    const HelmholtzOperators &operators,
    const ComplexVector &approximation,
    const std::vector<int> &elements,
    const ComplexFunction &exact_value,
    const ComplexGradientFunction &exact_gradient,
    const QuadraturePolicy &quadrature = {},
    const QuadratureContext &quadrature_context = {});

struct FineSpaceTimings {
    double operator_assembly_ms = 0.0;
    double factorization_ms = 0.0;
    double load_assembly_ms = 0.0;
    double solve_ms = 0.0;
    double error_evaluation_ms = 0.0;
    std::size_t factorization_builds = 0;
    std::size_t solve_count = 0;
    bool excluded_from_method_time = false;

    double total_ms() const;
    double included_method_time_ms() const;
};

struct CertificationAuditSolution {
    ComplexVector values;
};

// Algorithm-facing service. Its costs are included in CALOD method time and
// its solution may be used only for certification quantities.
class CertificationAuditService {
public:
    CertificationAuditService(
        TriMesh mesh,
        double wavenumber,
        std::vector<double> diffusion = {},
        std::vector<double> refractive_index = {},
        double boundary_beta = 1.0);
    ~CertificationAuditService();
    CertificationAuditService(CertificationAuditService &&) noexcept;
    CertificationAuditService &operator=(CertificationAuditService &&) noexcept;
    CertificationAuditService(const CertificationAuditService &) = delete;
    CertificationAuditService &operator=(const CertificationAuditService &) = delete;

    CertificationAuditSolution solve_load(const ComplexVector &load);
    CertificationAuditSolution solve_source(
        const ComplexFunction &source,
        const QuadraturePolicy &quadrature = {},
        const QuadratureContext &quadrature_context = {});

    const TriMesh &mesh() const;
    const HelmholtzOperators &operators() const;
    const FineSpaceTimings &timings() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

struct EvaluationReferenceError {
    HelmholtzError absolute;
    double relative_energy = 0.0;
    double relative_l2 = 0.0;
};

// Post-processing-only service. It intentionally has no reference-solution
// getter: MARK/STOP code can receive CertificationAuditService, but cannot
// obtain evaluation-reference values from this type.
class EvaluationReferenceService {
public:
    EvaluationReferenceService(
        TriMesh mesh,
        double wavenumber,
        std::vector<double> diffusion = {},
        std::vector<double> refractive_index = {},
        double boundary_beta = 1.0);
    ~EvaluationReferenceService();
    EvaluationReferenceService(EvaluationReferenceService &&) noexcept;
    EvaluationReferenceService &operator=(EvaluationReferenceService &&) noexcept;
    EvaluationReferenceService(const EvaluationReferenceService &) = delete;
    EvaluationReferenceService &operator=(const EvaluationReferenceService &) = delete;

    void prepare_load(const ComplexVector &load);
    void prepare_source(
        const ComplexFunction &source,
        const QuadraturePolicy &quadrature = {},
        const QuadratureContext &quadrature_context = {});
    bool ready() const;

    EvaluationReferenceError evaluate_candidate_on_reference(
        const ComplexVector &candidate_on_reference);
    EvaluationReferenceError evaluate_nested_candidate(
        const TriMesh &candidate_mesh,
        const ComplexVector &candidate_values);

    const FineSpaceTimings &timings() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace lod2d::helmholtz::adaptive

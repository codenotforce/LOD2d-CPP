#include "helmholtz/adaptive/error_control.h"

#include "helmholtz/adaptive/estimator.h"
#include "helmholtz/adaptive/hierarchy.h"

#include <Eigen/SparseLU>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <numeric>
#include <set>
#include <stdexcept>
#include <utility>

namespace lod2d::helmholtz::adaptive {
namespace {

using Clock = std::chrono::steady_clock;

double elapsed_ms(const Clock::time_point &start) {
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

struct TriangleGeometry {
    double area = 0.0;
    std::array<Eigen::Vector2d, 3> gradients;
};

TriangleGeometry triangle_geometry(const TriMesh &mesh, const Triangle &triangle) {
    const Point2 &first = mesh.nodes[triangle[0]];
    const Point2 &second = mesh.nodes[triangle[1]];
    const Point2 &third = mesh.nodes[triangle[2]];
    const double determinant =
        (second.x() - first.x()) * (third.y() - first.y())
        - (second.y() - first.y()) * (third.x() - first.x());
    if (std::abs(determinant) <= 1e-15)
        throw std::invalid_argument("error evaluation encountered a degenerate triangle");
    TriangleGeometry result;
    result.area = 0.5 * std::abs(determinant);
    result.gradients[0] = Eigen::Vector2d(
        second.y() - third.y(), third.x() - second.x()) / determinant;
    result.gradients[1] = Eigen::Vector2d(
        third.y() - first.y(), first.x() - third.x()) / determinant;
    result.gradients[2] = Eigen::Vector2d(
        first.y() - second.y(), second.x() - first.x()) / determinant;
    return result;
}

void validate_error_inputs(
    const TriMesh &mesh,
    const HelmholtzOperators &operators,
    const ComplexVector &first,
    const ComplexVector &second) {
    const int node_count = static_cast<int>(mesh.nodes.size());
    const int element_count = static_cast<int>(mesh.elems.size());
    if (first.size() != node_count || second.size() != node_count)
        throw std::invalid_argument("error vectors must match mesh nodes");
    if (operators.system.rows() != node_count || operators.system.cols() != node_count)
        throw std::invalid_argument("error operators do not match the mesh");
    if (!(operators.wavenumber > 0.0))
        throw std::invalid_argument("error operators require a positive wavenumber");
    if ((!operators.diffusion.empty()
         && static_cast<int>(operators.diffusion.size()) != element_count)
        || (!operators.refractive_index.empty()
            && static_cast<int>(operators.refractive_index.size()) != element_count)) {
        throw std::invalid_argument("error coefficients do not match mesh elements");
    }
}

std::vector<int> checked_elements(
    const TriMesh &mesh,
    const std::vector<int> &elements) {
    std::set<int> unique;
    for (int element : elements) {
        if (element < 0 || element >= static_cast<int>(mesh.elems.size()))
            throw std::out_of_range("local error element index is out of range");
        if (!unique.insert(element).second)
            throw std::invalid_argument("local error element list contains a duplicate");
    }
    return {unique.begin(), unique.end()};
}

std::vector<int> all_elements(const TriMesh &mesh) {
    std::vector<int> result(mesh.elems.size());
    for (int element = 0; element < static_cast<int>(result.size()); ++element)
        result[element] = element;
    return result;
}

double coefficient(const std::vector<double> &values, int element) {
    return values.empty() ? 1.0 : values[element];
}

struct SquaredError {
    double energy = 0.0;
    double l2 = 0.0;
};

HelmholtzError square_root_error(const SquaredError &squared) {
    return {
        std::sqrt(std::max(0.0, squared.energy)),
        std::sqrt(std::max(0.0, squared.l2))};
}

HelmholtzError checked_error(HelmholtzError error) {
    if (!std::isfinite(error.energy) || !std::isfinite(error.l2)
        || error.energy < 0.0 || error.l2 < 0.0) {
        throw std::runtime_error(
            "error reference returned a non-finite or negative error");
    }
    return error;
}

SquaredError discrete_squared_error(
    const TriMesh &mesh,
    const HelmholtzOperators &operators,
    const ComplexVector &difference,
    const std::vector<int> &elements) {
    SquaredError result;
    static constexpr double mass_pattern[3][3] = {
        {2.0, 1.0, 1.0},
        {1.0, 2.0, 1.0},
        {1.0, 1.0, 2.0}};
    for (int element : elements) {
        const Triangle &triangle = mesh.elems[element];
        const TriangleGeometry geometry = triangle_geometry(mesh, triangle);
        Eigen::Vector2cd gradient = Eigen::Vector2cd::Zero();
        for (int local = 0; local < 3; ++local)
            gradient += difference(triangle[local])
                * geometry.gradients[local].cast<Complex>();

        double l2_squared = 0.0;
        for (int first = 0; first < 3; ++first) {
            for (int second = 0; second < 3; ++second) {
                l2_squared += geometry.area * mass_pattern[first][second] / 12.0
                    * std::real(std::conj(difference(triangle[first]))
                                * difference(triangle[second]));
            }
        }
        const double gradient_squared = geometry.area * gradient.squaredNorm();
        result.l2 += l2_squared;
        result.energy += coefficient(operators.diffusion, element) * gradient_squared
            + operators.wavenumber * operators.wavenumber
                * coefficient(operators.refractive_index, element) * l2_squared;
    }
    return result;
}

class CachedFineSpaceFactorization {
public:
    CachedFineSpaceFactorization(
        TriMesh mesh,
        double wavenumber,
        std::vector<double> diffusion,
        std::vector<double> refractive_index,
        double boundary_beta,
        bool excluded_from_method_time)
        : mesh_(std::move(mesh)) {
        timings_.excluded_from_method_time = excluded_from_method_time;
        const auto assembly_start = Clock::now();
        operators_ = assemble_helmholtz_operators(
            mesh_, wavenumber, diffusion, refractive_index, boundary_beta);
        timings_.operator_assembly_ms = elapsed_ms(assembly_start);

        std::vector<char> is_dirichlet(mesh_.nodes.size(), false);
        for (int node : operators_.dirichlet_nodes) {
            if (node < 0 || node >= static_cast<int>(is_dirichlet.size()))
                throw std::invalid_argument("Dirichlet node index is out of range");
            is_dirichlet[node] = true;
        }
        std::vector<int> global_to_free(mesh_.nodes.size(), -1);
        for (int node = 0; node < static_cast<int>(mesh_.nodes.size()); ++node) {
            if (is_dirichlet[node]) continue;
            global_to_free[node] = static_cast<int>(free_nodes_.size());
            free_nodes_.push_back(node);
        }
        if (free_nodes_.empty())
            throw std::invalid_argument("fine-space solver has no unconstrained degrees of freedom");

        std::vector<ComplexTriplet> triplets;
        for (int global_column : free_nodes_) {
            const int local_column = global_to_free[global_column];
            for (ComplexSparseMatrix::InnerIterator it(
                     operators_.system, global_column); it; ++it) {
                const int local_row = global_to_free[it.row()];
                if (local_row >= 0)
                    triplets.emplace_back(local_row, local_column, it.value());
            }
        }
        ComplexSparseMatrix reduced(free_nodes_.size(), free_nodes_.size());
        reduced.setFromTriplets(triplets.begin(), triplets.end());
        reduced.makeCompressed();

        const auto factorization_start = Clock::now();
        solver_.analyzePattern(reduced);
        solver_.factorize(reduced);
        timings_.factorization_ms = elapsed_ms(factorization_start);
        if (solver_.info() != Eigen::Success)
            throw std::runtime_error("cached fine-space sparse LU factorization failed");
        timings_.factorization_builds = 1;
    }

    ComplexVector solve(const ComplexVector &load) {
        if (load.size() != static_cast<int>(mesh_.nodes.size()) || !load.allFinite())
            throw std::invalid_argument("cached fine-space load is invalid");
        ComplexVector reduced_load(free_nodes_.size());
        for (int local = 0; local < static_cast<int>(free_nodes_.size()); ++local)
            reduced_load(local) = load(free_nodes_[local]);
        const auto solve_start = Clock::now();
        const ComplexVector reduced_solution = solver_.solve(reduced_load);
        timings_.solve_ms += elapsed_ms(solve_start);
        if (solver_.info() != Eigen::Success || !reduced_solution.allFinite())
            throw std::runtime_error("cached fine-space sparse LU solve failed");
        ++timings_.solve_count;
        ComplexVector solution = ComplexVector::Zero(mesh_.nodes.size());
        for (int local = 0; local < static_cast<int>(free_nodes_.size()); ++local)
            solution(free_nodes_[local]) = reduced_solution(local);
        return solution;
    }

    ComplexVector solve_source(
        const ComplexFunction &source,
        const QuadraturePolicy &quadrature,
        const QuadratureContext &quadrature_context) {
        const auto load_start = Clock::now();
        ComplexVector load = assemble_helmholtz_load(
            mesh_, source, quadrature, quadrature_context);
        timings_.load_assembly_ms += elapsed_ms(load_start);
        return solve(load);
    }

    const TriMesh &mesh() const { return mesh_; }
    const HelmholtzOperators &operators() const { return operators_; }
    FineSpaceTimings &timings() { return timings_; }
    const FineSpaceTimings &timings() const { return timings_; }

private:
    TriMesh mesh_;
    HelmholtzOperators operators_;
    std::vector<int> free_nodes_;
    Eigen::SparseLU<ComplexSparseMatrix> solver_;
    FineSpaceTimings timings_;
};

} // namespace

HelmholtzError compute_local_discrete_helmholtz_error(
    const TriMesh &mesh,
    const HelmholtzOperators &operators,
    const ComplexVector &reference,
    const ComplexVector &approximation,
    const std::vector<int> &elements) {
    validate_error_inputs(mesh, operators, reference, approximation);
    const std::vector<int> selected = checked_elements(mesh, elements);
    return square_root_error(discrete_squared_error(
        mesh, operators, reference - approximation, selected));
}

HelmholtzError compute_discrete_helmholtz_error(
    const TriMesh &mesh,
    const HelmholtzOperators &operators,
    const ComplexVector &reference,
    const ComplexVector &approximation) {
    return compute_local_discrete_helmholtz_error(
        mesh, operators, reference, approximation, all_elements(mesh));
}

HelmholtzError compute_local_exact_helmholtz_error(
    const TriMesh &mesh,
    const HelmholtzOperators &operators,
    const ComplexVector &approximation,
    const std::vector<int> &elements,
    const ComplexFunction &exact_value,
    const ComplexGradientFunction &exact_gradient,
    const QuadraturePolicy &quadrature,
    const QuadratureContext &quadrature_context) {
    validate_error_inputs(mesh, operators, approximation, approximation);
    if (!exact_value || !exact_gradient)
        throw std::invalid_argument("exact error role requires value and gradient callbacks");
    const std::vector<int> selected = checked_elements(mesh, elements);
    SquaredError squared;
    for (int element : selected) {
        const Triangle &triangle = mesh.elems[element];
        const TriangleGeometry geometry = triangle_geometry(mesh, triangle);
        Eigen::Vector2cd discrete_gradient = Eigen::Vector2cd::Zero();
        for (int local = 0; local < 3; ++local)
            discrete_gradient += approximation(triangle[local])
                * geometry.gradients[local].cast<Complex>();
        for (const auto &point : triangle_quadrature_points(
                 mesh, element, quadrature, quadrature_context)) {
            Complex discrete_value = 0.0;
            for (int local = 0; local < 3; ++local)
                discrete_value += point.barycentric[local]
                    * approximation(triangle[local]);
            const Complex value_error = exact_value(point.point) - discrete_value;
            const Eigen::Vector2cd gradient_error =
                exact_gradient(point.point) - discrete_gradient;
            const double l2_contribution = point.weight * std::norm(value_error);
            squared.l2 += l2_contribution;
            squared.energy += coefficient(operators.diffusion, element)
                    * point.weight * gradient_error.squaredNorm()
                + operators.wavenumber * operators.wavenumber
                    * coefficient(operators.refractive_index, element)
                    * l2_contribution;
        }
    }
    return square_root_error(squared);
}

EmpiricalSaturationAuditEstimate estimate_empirical_saturation_audit_error(
    const TriMesh &audit_mesh,
    const ComplexVector &audit_solution,
    const ComplexFunction &source,
    double wavenumber,
    const std::vector<int> &audit_parent_fine_elements,
    int fine_element_count,
    double saturation_factor,
    double doerfler_theta,
    const std::vector<double> &audit_diffusion,
    const std::vector<double> &audit_refractive_index,
    double boundary_beta,
    const QuadraturePolicy &quadrature,
    const QuadratureContext &quadrature_context) {
    const int audit_node_count = static_cast<int>(audit_mesh.nodes.size());
    const int audit_element_count = static_cast<int>(audit_mesh.elems.size());
    if (audit_solution.size() != audit_node_count)
        throw std::invalid_argument("audit solution does not match the audit mesh");
    if (!source)
        throw std::invalid_argument("empirical audit estimate requires a source");
    if (!(wavenumber > 0.0) || !(boundary_beta >= 0.0))
        throw std::invalid_argument("empirical audit PDE parameters are invalid");
    if (!(saturation_factor >= 0.0 && saturation_factor < 1.0))
        throw std::invalid_argument("audit saturation factor must lie in [0,1)");
    if (!(doerfler_theta > 0.0 && doerfler_theta <= 1.0))
        throw std::invalid_argument("audit Doerfler parameter must lie in (0,1]");
    if (fine_element_count <= 0
        || static_cast<int>(audit_parent_fine_elements.size())
            != audit_element_count)
        throw std::invalid_argument("audit-to-fine parent map is incompatible");
    for (int parent : audit_parent_fine_elements) {
        if (parent < 0 || parent >= fine_element_count)
            throw std::out_of_range("audit-to-fine parent index is invalid");
    }
    const auto validate_coefficients = [&](const std::vector<double> &values,
                                           const char *name) {
        if (!values.empty()
            && static_cast<int>(values.size()) != audit_element_count)
            throw std::invalid_argument(std::string(name)
                                        + " does not match audit elements");
    };
    validate_coefficients(audit_diffusion, "audit diffusion");
    validate_coefficients(audit_refractive_index, "audit refractive index");

    const RefineOutput probe = refine_nvb(audit_mesh);
    const auto prolong_coefficients = [&](const std::vector<double> &values) {
        if (values.empty()) return std::vector<double>{};
        const Eigen::Map<const Eigen::VectorXd> parent_values(
            values.data(), static_cast<Eigen::Index>(values.size()));
        const Eigen::VectorXd child_values = probe.P_elem * parent_values;
        return std::vector<double>(child_values.data(),
                                   child_values.data() + child_values.size());
    };
    const std::vector<double> probe_diffusion =
        prolong_coefficients(audit_diffusion);
    const std::vector<double> probe_refractive_index =
        prolong_coefficients(audit_refractive_index);
    const HelmholtzOperators probe_operators = assemble_helmholtz_operators(
        probe.mesh, wavenumber, probe_diffusion, probe_refractive_index,
        boundary_beta);
    const ComplexVector probe_load = assemble_helmholtz_load(
        probe.mesh, source, quadrature, quadrature_context);
    const ComplexVector probe_solution =
        solve_helmholtz_fem(probe_operators, probe_load);
    const ComplexVector prolonged_audit =
        probe.P_node.cast<Complex>() * audit_solution;
    const HelmholtzError two_level = compute_discrete_helmholtz_error(
        probe.mesh, probe_operators, probe_solution, prolonged_audit);

    EmpiricalSaturationAuditEstimate result;
    result.saturation_factor = saturation_factor;
    result.two_level_energy_error = two_level.energy;
    result.audit_error_lower = two_level.energy / (1.0 + saturation_factor);
    result.audit_error_upper = two_level.energy / (1.0 - saturation_factor);
    result.audit_element_error_squared.assign(audit_element_count, 0.0);
    result.fine_element_error_squared.assign(fine_element_count, 0.0);

    std::vector<std::vector<int>> probe_children(audit_element_count);
    for (int parent = 0; parent < probe.P_elem.outerSize(); ++parent) {
        for (Eigen::SparseMatrix<double>::InnerIterator entry(
                 probe.P_elem, parent);
             entry; ++entry) {
            if (std::abs(entry.value()) > 0.0)
                probe_children[parent].push_back(entry.row());
        }
    }
    for (int audit_element = 0; audit_element < audit_element_count;
         ++audit_element) {
        if (probe_children[audit_element].empty())
            throw std::runtime_error("uniform audit probe lost a parent element");
        const HelmholtzError local = compute_local_discrete_helmholtz_error(
            probe.mesh, probe_operators, probe_solution, prolonged_audit,
            probe_children[audit_element]);
        const double squared = local.energy * local.energy;
        result.audit_element_error_squared[audit_element] = squared;
        result.fine_element_error_squared[
            audit_parent_fine_elements[audit_element]] += squared;
    }
    const double allocated = std::accumulate(
        result.audit_element_error_squared.begin(),
        result.audit_element_error_squared.end(), 0.0);
    const double expected = two_level.energy * two_level.energy;
    result.allocation_relative_error = std::abs(allocated - expected)
        / std::max(1.0, expected);
    if (allocated > 0.0) {
        result.marked_fine_elements = mark_doerfler(
            result.fine_element_error_squared, doerfler_theta);
    }
    return result;
}

ErrorReference ErrorReference::exact(
    ComplexFunction exact_value,
    ComplexGradientFunction exact_gradient,
    QuadraturePolicy quadrature,
    QuadratureContext quadrature_context) {
    if (!exact_value || !exact_gradient)
        throw std::invalid_argument("exact error reference requires value and gradient callbacks");
    validate_quadrature_policy(quadrature);
    ErrorReference result;
    result.role_ = ErrorReferenceRole::Exact;
    result.exact_value_ = std::move(exact_value);
    result.exact_gradient_ = std::move(exact_gradient);
    result.quadrature_ = quadrature;
    result.quadrature_context_ = quadrature_context;
    return result;
}

ErrorReference ErrorReference::external(ExternalErrorEstimator estimator) {
    if (!estimator)
        throw std::invalid_argument("external error reference requires an estimator callback");
    ErrorReference result;
    result.role_ = ErrorReferenceRole::ExternalEstimator;
    result.external_estimator_ = std::move(estimator);
    return result;
}

ErrorReference ErrorReference::two_level(ComplexVector finer_values) {
    if (finer_values.size() == 0 || !finer_values.allFinite())
        throw std::invalid_argument("two-level error reference requires finite finer values");
    ErrorReference result;
    result.role_ = ErrorReferenceRole::TwoLevelSaturation;
    result.finer_values_ = std::move(finer_values);
    return result;
}

HelmholtzError ErrorReference::evaluate(
    const TriMesh &mesh,
    const HelmholtzOperators &operators,
    const ComplexVector &candidate) const {
    HelmholtzError result;
    switch (role_) {
    case ErrorReferenceRole::Exact:
        result = compute_local_exact_helmholtz_error(
            mesh, operators, candidate, all_elements(mesh), exact_value_, exact_gradient_,
            quadrature_, quadrature_context_);
        break;
    case ErrorReferenceRole::ExternalEstimator:
        result = external_estimator_(mesh, candidate);
        break;
    case ErrorReferenceRole::TwoLevelSaturation:
        result = compute_discrete_helmholtz_error(
            mesh, operators, finer_values_, candidate);
        break;
    }
    return checked_error(result);
}

double FineSpaceTimings::total_ms() const {
    return operator_assembly_ms + factorization_ms + load_assembly_ms
        + solve_ms + error_evaluation_ms;
}

double FineSpaceTimings::included_method_time_ms() const {
    return excluded_from_method_time ? 0.0 : total_ms();
}

struct CertificationAuditService::Impl {
    Impl(
        TriMesh mesh,
        double wavenumber,
        std::vector<double> diffusion,
        std::vector<double> refractive_index,
        double boundary_beta)
        : cache(
            std::move(mesh), wavenumber, std::move(diffusion),
            std::move(refractive_index), boundary_beta, false) {}

    CachedFineSpaceFactorization cache;
};

CertificationAuditService::CertificationAuditService(
    TriMesh mesh,
    double wavenumber,
    std::vector<double> diffusion,
    std::vector<double> refractive_index,
    double boundary_beta)
    : impl_(std::make_unique<Impl>(
          std::move(mesh), wavenumber, std::move(diffusion),
          std::move(refractive_index), boundary_beta)) {}

CertificationAuditService::~CertificationAuditService() = default;
CertificationAuditService::CertificationAuditService(
    CertificationAuditService &&) noexcept = default;
CertificationAuditService &CertificationAuditService::operator=(
    CertificationAuditService &&) noexcept = default;

CertificationAuditSolution CertificationAuditService::solve_load(
    const ComplexVector &load) {
    return {impl_->cache.solve(load)};
}

CertificationAuditSolution CertificationAuditService::solve_source(
    const ComplexFunction &source,
    const QuadraturePolicy &quadrature,
    const QuadratureContext &quadrature_context) {
    return {impl_->cache.solve_source(source, quadrature, quadrature_context)};
}

const TriMesh &CertificationAuditService::mesh() const {
    return impl_->cache.mesh();
}

const HelmholtzOperators &CertificationAuditService::operators() const {
    return impl_->cache.operators();
}

const FineSpaceTimings &CertificationAuditService::timings() const {
    return impl_->cache.timings();
}

struct EvaluationReferenceService::Impl {
    Impl(
        TriMesh mesh,
        double wavenumber,
        std::vector<double> diffusion,
        std::vector<double> refractive_index,
        double boundary_beta)
        : cache(
            std::move(mesh), wavenumber, std::move(diffusion),
            std::move(refractive_index), boundary_beta, true) {}

    EvaluationReferenceError evaluate(const ComplexVector &candidate) {
        if (!has_reference)
            throw std::logic_error("evaluation reference has not been prepared");
        const auto start = Clock::now();
        EvaluationReferenceError result;
        result.absolute = compute_discrete_helmholtz_error(
            cache.mesh(), cache.operators(), reference_values, candidate);
        const ComplexVector zero = ComplexVector::Zero(reference_values.size());
        const HelmholtzError reference_norm = compute_discrete_helmholtz_error(
            cache.mesh(), cache.operators(), reference_values, zero);
        if (!(reference_norm.energy > 0.0) || !(reference_norm.l2 > 0.0))
            throw std::runtime_error("evaluation reference has a zero error norm denominator");
        result.relative_energy = result.absolute.energy / reference_norm.energy;
        result.relative_l2 = result.absolute.l2 / reference_norm.l2;
        cache.timings().error_evaluation_ms += elapsed_ms(start);
        return result;
    }

    CachedFineSpaceFactorization cache;
    ComplexVector reference_values;
    bool has_reference = false;
};

EvaluationReferenceService::EvaluationReferenceService(
    TriMesh mesh,
    double wavenumber,
    std::vector<double> diffusion,
    std::vector<double> refractive_index,
    double boundary_beta)
    : impl_(std::make_unique<Impl>(
          std::move(mesh), wavenumber, std::move(diffusion),
          std::move(refractive_index), boundary_beta)) {}

EvaluationReferenceService::~EvaluationReferenceService() = default;
EvaluationReferenceService::EvaluationReferenceService(
    EvaluationReferenceService &&) noexcept = default;
EvaluationReferenceService &EvaluationReferenceService::operator=(
    EvaluationReferenceService &&) noexcept = default;

void EvaluationReferenceService::prepare_load(const ComplexVector &load) {
    impl_->reference_values = impl_->cache.solve(load);
    impl_->has_reference = true;
}

void EvaluationReferenceService::prepare_source(
    const ComplexFunction &source,
    const QuadraturePolicy &quadrature,
    const QuadratureContext &quadrature_context) {
    impl_->reference_values = impl_->cache.solve_source(
        source, quadrature, quadrature_context);
    impl_->has_reference = true;
}

bool EvaluationReferenceService::ready() const {
    return impl_->has_reference;
}

EvaluationReferenceError EvaluationReferenceService::evaluate_candidate_on_reference(
    const ComplexVector &candidate_on_reference) {
    return impl_->evaluate(candidate_on_reference);
}

EvaluationReferenceError EvaluationReferenceService::evaluate_nested_candidate(
    const TriMesh &candidate_mesh,
    const ComplexVector &candidate_values) {
    if (candidate_values.size() != static_cast<int>(candidate_mesh.nodes.size()))
        throw std::invalid_argument("nested evaluation candidate does not match its mesh");
    const auto embedding_start = Clock::now();
    const RefineOutput embedding = build_nested_mesh_embedding(
        candidate_mesh, impl_->cache.mesh());
    const ComplexVector prolonged =
        embedding.P_node.cast<Complex>() * candidate_values;
    impl_->cache.timings().error_evaluation_ms += elapsed_ms(embedding_start);
    return impl_->evaluate(prolonged);
}

const FineSpaceTimings &EvaluationReferenceService::timings() const {
    return impl_->cache.timings();
}

} // namespace lod2d::helmholtz::adaptive

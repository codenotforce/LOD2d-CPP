#include "helmholtz/adaptive/kernel_residual.h"

#include "helmholtz/adaptive/error_control.h"
#include "helmholtz/adaptive/estimator.h"
#include "helmholtz/boundary.h"

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <Eigen/SparseLU>

#include <algorithm>
#include <atomic>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <limits>
#include <mutex>
#include <numeric>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lod2d::helmholtz::adaptive {
namespace {

constexpr double kSupportTolerance = 1e-13;
constexpr double kRankTolerance = 2e-12;

class FingerprintBuilder {
public:
    void add_u64(std::uint64_t value) {
        for (int byte = 0; byte < 8; ++byte) {
            hash_ ^= static_cast<unsigned char>((value >> (8 * byte)) & 0xffU);
            hash_ *= 1099511628211ULL;
        }
    }

    void add_i64(std::int64_t value) {
        add_u64(static_cast<std::uint64_t>(value));
    }

    void add_double(double value) {
        // Positive and negative zero are numerically indistinguishable in all
        // estimator operations, so give them one canonical representation.
        if (value == 0.0) value = 0.0;
        add_u64(std::bit_cast<std::uint64_t>(value));
    }

    void add_complex(const Complex &value) {
        add_double(value.real());
        add_double(value.imag());
    }

    void add_string(const std::string &value) {
        add_u64(value.size());
        for (unsigned char byte : value) {
            hash_ ^= byte;
            hash_ *= 1099511628211ULL;
        }
    }

    std::string finish() const {
        std::ostringstream encoded;
        encoded << "fnv1a64:" << std::hex << std::setw(16)
                << std::setfill('0') << hash_;
        return encoded.str();
    }

private:
    std::uint64_t hash_ = 14695981039346656037ULL;
};

void add_fingerprint_value(FingerprintBuilder &builder, double value) {
    builder.add_double(value);
}

void add_fingerprint_value(FingerprintBuilder &builder, const Complex &value) {
    builder.add_complex(value);
}

template <typename SparseMatrix>
void add_sparse_fingerprint(
    FingerprintBuilder &builder,
    const SparseMatrix &matrix) {
    builder.add_i64(matrix.rows());
    builder.add_i64(matrix.cols());
    builder.add_i64(matrix.nonZeros());
    for (int outer = 0; outer < matrix.outerSize(); ++outer) {
        for (typename SparseMatrix::InnerIterator it(matrix, outer); it; ++it) {
            builder.add_i64(it.row());
            builder.add_i64(it.col());
            add_fingerprint_value(builder, it.value());
        }
    }
}

void add_mesh_fingerprint(FingerprintBuilder &builder, const TriMesh &mesh) {
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
        builder.add_u64(static_cast<std::uint8_t>(edge.tag));
    }
}

template <typename Value>
void add_integral_vector_fingerprint(
    FingerprintBuilder &builder,
    const std::vector<Value> &values) {
    builder.add_u64(values.size());
    for (Value value : values) builder.add_u64(static_cast<std::uint64_t>(value));
}

void add_double_vector_fingerprint(
    FingerprintBuilder &builder,
    const std::vector<double> &values) {
    builder.add_u64(values.size());
    for (double value : values) builder.add_double(value);
}

void add_complex_vector_fingerprint(
    FingerprintBuilder &builder,
    const ComplexVector &values) {
    builder.add_i64(values.size());
    for (int index = 0; index < values.size(); ++index)
        builder.add_complex(values(index));
}

struct ConstraintReduction {
    std::vector<int> active_rows;
    std::vector<int> independent_rows;
    Eigen::MatrixXd matrix;
};

struct PatchConstructionData {
    std::vector<std::vector<int>> graph;
    std::vector<std::vector<int>> hat_support;
    std::vector<std::vector<int>> interpolation_support;
    std::vector<int> discrete_element_parent;
};

std::vector<int> sorted_unique(std::vector<int> values) {
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
    return values;
}

std::vector<std::vector<int>> coarse_element_graph(const TriMesh &mesh) {
    std::vector<std::vector<int>> incident(mesh.nodes.size());
    for (int element = 0; element < static_cast<int>(mesh.elems.size()); ++element) {
        for (int node : mesh.elems[element]) incident[node].push_back(element);
    }
    std::vector<std::vector<int>> graph(mesh.elems.size());
    for (const std::vector<int> &star : incident) {
        for (int first : star) {
            graph[first].insert(graph[first].end(), star.begin(), star.end());
        }
    }
    for (std::vector<int> &neighbors : graph) neighbors = sorted_unique(std::move(neighbors));
    return graph;
}

std::vector<std::vector<int>> coarse_hat_support(const TriMesh &mesh) {
    std::vector<std::vector<int>> result(mesh.nodes.size());
    for (int element = 0; element < static_cast<int>(mesh.elems.size()); ++element) {
        for (int node : mesh.elems[element]) result[node].push_back(element);
    }
    for (std::vector<int> &support : result) support = sorted_unique(std::move(support));
    return result;
}

std::vector<int> expand_elements(
    const std::vector<std::vector<int>> &graph,
    const std::vector<int> &initial,
    int layers) {
    if (layers < 0) throw std::invalid_argument("patch layer count must be nonnegative");
    std::vector<char> selected(graph.size(), false);
    std::vector<int> frontier;
    for (int element : initial) {
        if (element < 0 || element >= static_cast<int>(graph.size()))
            throw std::out_of_range("patch contains an out-of-range coarse element");
        if (!selected[element]) {
            selected[element] = true;
            frontier.push_back(element);
        }
    }
    for (int layer = 0; layer < layers && !frontier.empty(); ++layer) {
        std::vector<int> next;
        for (int element : frontier) {
            for (int neighbor : graph[element]) {
                if (selected[neighbor]) continue;
                selected[neighbor] = true;
                next.push_back(neighbor);
            }
        }
        frontier = std::move(next);
    }
    std::vector<int> result;
    for (int element = 0; element < static_cast<int>(selected.size()); ++element) {
        if (selected[element]) result.push_back(element);
    }
    return result;
}

int support_distance(
    const std::vector<std::vector<int>> &graph,
    const std::vector<int> &sources,
    const std::vector<int> &targets) {
    if (targets.empty()) return 0;
    std::vector<int> distance(graph.size(), -1);
    std::queue<int> queue;
    for (int source : sources) {
        distance[source] = 0;
        queue.push(source);
    }
    while (!queue.empty()) {
        const int current = queue.front();
        queue.pop();
        for (int neighbor : graph[current]) {
            if (distance[neighbor] >= 0) continue;
            distance[neighbor] = distance[current] + 1;
            queue.push(neighbor);
        }
    }
    int result = 0;
    for (int target : targets) {
        if (target < 0 || target >= static_cast<int>(distance.size())
            || distance[target] < 0) {
            throw std::runtime_error("interpolation support is disconnected from its coarse hat");
        }
        result = std::max(result, distance[target]);
    }
    return result;
}

bool contains_all(const std::vector<int> &container, const std::vector<int> &values) {
    for (int value : values) {
        if (!std::binary_search(container.begin(), container.end(), value)) return false;
    }
    return true;
}

std::string policy_hash(const KernelPatchPolicy &policy) {
    std::ostringstream description;
    description << "v1|" << policy.interpolation_support_layers << '|'
                << policy.patch_layers << '|' << policy.enlargement_layers << '|'
                << policy.adjacency << '|' << policy.audit_enlargement;
    std::uint64_t hash = 14695981039346656037ULL;
    for (unsigned char byte : description.str()) {
        hash ^= byte;
        hash *= 1099511628211ULL;
    }
    std::ostringstream encoded;
    encoded << "fnv1a64:" << std::hex << std::setw(16) << std::setfill('0') << hash;
    return encoded.str();
}

std::string audit_kernel_context_fingerprint(
    const AdaptiveMeshHierarchy &hierarchy,
    const HelmholtzOperators &operators,
    const KernelPatchPolicy &policy) {
    FingerprintBuilder builder;
    builder.add_string("audit-kernel-context-v1");
    add_mesh_fingerprint(builder, hierarchy.initial_mesh());
    add_mesh_fingerprint(builder, hierarchy.coarse_mesh());
    add_mesh_fingerprint(builder, hierarchy.fine_mesh());
    add_mesh_fingerprint(builder, hierarchy.cert_audit_mesh());
    add_integral_vector_fingerprint(builder, hierarchy.coarse_levels());
    add_integral_vector_fingerprint(builder, hierarchy.coarse_element_ids());
    add_integral_vector_fingerprint(builder, hierarchy.coarse_parent_ids());
    add_integral_vector_fingerprint(builder, hierarchy.fine_element_levels());
    add_integral_vector_fingerprint(
        builder, hierarchy.cert_audit_element_levels());
    builder.add_u64(hierarchy.coarse_mesh_version());
    builder.add_u64(hierarchy.fine_mesh_version());
    builder.add_u64(hierarchy.cert_audit_mesh_version());
    builder.add_u64(hierarchy.interpolation_version());
    builder.add_u64(hierarchy.boundary_version());
    builder.add_u64(hierarchy.corrector_space_version());
    add_sparse_fingerprint(builder, hierarchy.coarse_to_cert_audit());
    add_sparse_fingerprint(builder, hierarchy.fine_to_cert_audit());
    add_sparse_fingerprint(builder, hierarchy.coarse_elements_to_cert_audit());
    add_sparse_fingerprint(builder, hierarchy.cert_audit_quasi_interpolation());
    builder.add_string(policy.hash);

    builder.add_double(operators.wavenumber);
    builder.add_double(operators.boundary_beta);
    add_double_vector_fingerprint(builder, operators.diffusion);
    add_double_vector_fingerprint(builder, operators.refractive_index);
    add_integral_vector_fingerprint(builder, operators.dirichlet_nodes);
    add_sparse_fingerprint(builder, operators.stiffness);
    add_sparse_fingerprint(builder, operators.mass);
    add_sparse_fingerprint(builder, operators.boundary_mass);
    add_sparse_fingerprint(builder, operators.system);
    return builder.finish();
}

std::string audit_kernel_diagnostic_fingerprint(
    const std::string &context_fingerprint,
    KernelRieszSolver solver,
    double doerfler_theta,
    const ComplexVector &audit_load,
    const ComplexVector &candidate_on_audit,
    const AuditKernelResidualEstimate &estimate) {
    FingerprintBuilder builder;
    builder.add_string("audit-kernel-diagnostic-v1");
    builder.add_string(context_fingerprint);
    builder.add_u64(static_cast<std::uint64_t>(solver));
    builder.add_double(doerfler_theta);
    add_complex_vector_fingerprint(builder, audit_load);
    add_complex_vector_fingerprint(builder, candidate_on_audit);
    add_complex_vector_fingerprint(builder, estimate.global_residual);
    add_double_vector_fingerprint(builder, estimate.node_eta);
    add_double_vector_fingerprint(builder, estimate.node_eta_squared);
    add_double_vector_fingerprint(builder, estimate.element_eta_squared);
    add_integral_vector_fingerprint(builder, estimate.marked_elements);
    builder.add_double(estimate.eta);
    builder.add_double(estimate.allocation_relative_error);
    return builder.finish();
}

std::string audit_kernel_result_fingerprint(
    const AuditKernelResidualEstimate &estimate) {
    FingerprintBuilder builder;
    builder.add_string("audit-kernel-result-v1");
    builder.add_u64(static_cast<std::uint64_t>(estimate.solver));
    builder.add_i64(estimate.policy.interpolation_support_layers);
    builder.add_i64(estimate.policy.patch_layers);
    builder.add_i64(estimate.policy.enlargement_layers);
    builder.add_string(estimate.policy.adjacency);
    builder.add_string(estimate.policy.audit_enlargement);
    builder.add_string(estimate.policy.hash);
    add_complex_vector_fingerprint(builder, estimate.global_residual);
    builder.add_u64(estimate.patches.size());
    for (const AuditKernelPatch &patch : estimate.patches) {
        builder.add_i64(patch.coarse_node);
        add_integral_vector_fingerprint(builder, patch.coarse_hat_support);
        add_integral_vector_fingerprint(builder, patch.interpolation_support);
        add_integral_vector_fingerprint(builder, patch.coarse_elements);
        add_integral_vector_fingerprint(
            builder, patch.enlarged_coarse_elements);
        add_integral_vector_fingerprint(builder, patch.discrete_dofs);
        add_integral_vector_fingerprint(
            builder, patch.enlarged_discrete_elements);
        add_integral_vector_fingerprint(
            builder, patch.active_constraint_rows);
        add_integral_vector_fingerprint(
            builder, patch.independent_constraint_rows);
        builder.add_i64(patch.constraints.rows());
        builder.add_i64(patch.constraints.cols());
        for (Eigen::Index row = 0; row < patch.constraints.rows(); ++row) {
            for (Eigen::Index column = 0;
                 column < patch.constraints.cols(); ++column) {
                builder.add_double(patch.constraints(row, column));
            }
        }
    }
    builder.add_u64(estimate.local_results.size());
    for (const LocalKernelRieszResult &local : estimate.local_results) {
        builder.add_double(local.eta);
        builder.add_double(local.eta_squared);
        builder.add_double(local.constraint_relative_residual);
        builder.add_double(local.riesz_relative_residual);
        builder.add_double(local.energy_identity_relative_error);
        builder.add_i64(local.kernel_dimension);
        builder.add_i64(local.saddle_unknowns);
        add_complex_vector_fingerprint(builder, local.local_values);
        add_complex_vector_fingerprint(builder, local.lagrange_multipliers);
    }
    add_double_vector_fingerprint(builder, estimate.node_eta);
    add_double_vector_fingerprint(builder, estimate.node_eta_squared);
    add_double_vector_fingerprint(builder, estimate.element_eta_squared);
    add_integral_vector_fingerprint(builder, estimate.marked_elements);
    builder.add_double(estimate.eta);
    builder.add_double(estimate.allocation_relative_error);
    return builder.finish();
}

PatchConstructionData patch_construction_data(
    const AdaptiveMeshHierarchy &hierarchy) {
    const TriMesh &coarse = hierarchy.coarse_mesh();
    const TriMesh &audit = hierarchy.cert_audit_mesh();
    const Eigen::SparseMatrix<double> &interpolation =
        hierarchy.cert_audit_quasi_interpolation();
    if (interpolation.rows() != static_cast<int>(coarse.nodes.size())
        || interpolation.cols() != static_cast<int>(audit.nodes.size())) {
        throw std::invalid_argument("audit quasi-interpolation dimensions do not match the hierarchy");
    }

    PatchConstructionData result;
    result.graph = coarse_element_graph(coarse);
    result.hat_support = coarse_hat_support(coarse);
    result.discrete_element_parent = fine_element_parents(
        hierarchy.coarse_elements_to_cert_audit(),
        static_cast<int>(audit.elems.size()),
        static_cast<int>(coarse.elems.size()));

    std::vector<std::vector<int>> audit_node_parents(audit.nodes.size());
    for (int element = 0; element < static_cast<int>(audit.elems.size()); ++element) {
        const int parent = result.discrete_element_parent[element];
        for (int node : audit.elems[element]) audit_node_parents[node].push_back(parent);
    }
    for (std::vector<int> &parents : audit_node_parents)
        parents = sorted_unique(std::move(parents));

    result.interpolation_support.resize(coarse.nodes.size());
    for (int audit_node = 0; audit_node < interpolation.outerSize(); ++audit_node) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(interpolation, audit_node); it; ++it) {
            if (std::abs(it.value()) <= kSupportTolerance) continue;
            std::vector<int> &support = result.interpolation_support[it.row()];
            support.insert(
                support.end(),
                audit_node_parents[audit_node].begin(),
                audit_node_parents[audit_node].end());
        }
    }
    for (std::vector<int> &support : result.interpolation_support)
        support = sorted_unique(std::move(support));
    return result;
}

ConstraintReduction reduce_constraints(
    const Eigen::SparseMatrix<double> &interpolation,
    const std::vector<int> &discrete_dofs);

struct ReferenceEpochSpaceView {
    const TriMesh *coarse = nullptr;
    const TriMesh *discrete = nullptr;
    const Eigen::SparseMatrix<double> *interpolation = nullptr;
    const Eigen::SparseMatrix<double> *coarse_elements_to_discrete = nullptr;
};

ReferenceEpochSpaceView reference_epoch_space_view(
    const ReferenceEpochHierarchy &hierarchy,
    KernelRieszSpace space) {
    ReferenceEpochSpaceView result;
    result.coarse = &hierarchy.coarse_mesh();
    switch (space) {
    case KernelRieszSpace::ReferenceResidual:
        result.discrete = &hierarchy.reference_mesh();
        result.interpolation = &hierarchy.reference_quasi_interpolation();
        result.coarse_elements_to_discrete =
            &hierarchy.coarse_elements_to_reference();
        break;
    case KernelRieszSpace::AmbientDefect:
        result.discrete = &hierarchy.ambient_mesh();
        result.interpolation = &hierarchy.ambient_quasi_interpolation();
        result.coarse_elements_to_discrete =
            &hierarchy.coarse_elements_to_ambient();
        break;
    default:
        throw std::invalid_argument("unknown kernel Riesz discrete space");
    }
    return result;
}

PatchConstructionData patch_construction_data(
    const ReferenceEpochSpaceView &view) {
    if (view.coarse == nullptr || view.discrete == nullptr
        || view.interpolation == nullptr
        || view.coarse_elements_to_discrete == nullptr) {
        throw std::invalid_argument("kernel Riesz space view is incomplete");
    }
    if (view.interpolation->rows()
            != static_cast<int>(view.coarse->nodes.size())
        || view.interpolation->cols()
            != static_cast<int>(view.discrete->nodes.size())) {
        throw std::invalid_argument(
            "kernel Riesz interpolation dimensions do not match the selected space");
    }

    PatchConstructionData result;
    result.graph = coarse_element_graph(*view.coarse);
    result.hat_support = coarse_hat_support(*view.coarse);
    result.discrete_element_parent = fine_element_parents(
        *view.coarse_elements_to_discrete,
        static_cast<int>(view.discrete->elems.size()),
        static_cast<int>(view.coarse->elems.size()));

    std::vector<std::vector<int>> node_parents(view.discrete->nodes.size());
    for (int element = 0;
         element < static_cast<int>(view.discrete->elems.size()); ++element) {
        const int parent = result.discrete_element_parent[element];
        for (int node : view.discrete->elems[element])
            node_parents[node].push_back(parent);
    }
    for (std::vector<int> &parents : node_parents)
        parents = sorted_unique(std::move(parents));

    result.interpolation_support.resize(view.coarse->nodes.size());
    for (int node = 0; node < view.interpolation->outerSize(); ++node) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(
                 *view.interpolation, node); it; ++it) {
            if (std::abs(it.value()) <= kSupportTolerance) continue;
            std::vector<int> &support = result.interpolation_support[it.row()];
            support.insert(
                support.end(), node_parents[node].begin(),
                node_parents[node].end());
        }
    }
    for (std::vector<int> &support : result.interpolation_support)
        support = sorted_unique(std::move(support));
    return result;
}

KernelPatchPolicy kernel_patch_policy(const PatchConstructionData &data) {
    int propagation_layers = 0;
    for (int node = 0; node < static_cast<int>(data.hat_support.size()); ++node) {
        propagation_layers = std::max(
            propagation_layers,
            support_distance(
                data.graph, data.hat_support[node],
                data.interpolation_support[node]));
    }
    KernelPatchPolicy result;
    result.interpolation_support_layers = propagation_layers;
    result.patch_layers = propagation_layers + 1;
    result.enlargement_layers = 1;
    result.adjacency = "shared_coarse_vertex";
    result.audit_enlargement = "one_shared_coarse_vertex_layer";
    result.hash = policy_hash(result);
    return result;
}

bool same_policy(const KernelPatchPolicy &left,
                 const KernelPatchPolicy &right) {
    return left.interpolation_support_layers
            == right.interpolation_support_layers
        && left.patch_layers == right.patch_layers
        && left.enlargement_layers == right.enlargement_layers
        && left.adjacency == right.adjacency
        && left.audit_enlargement == right.audit_enlargement
        && left.hash == right.hash;
}

std::vector<KernelRieszPatch> build_kernel_riesz_patches_impl(
    const ReferenceEpochSpaceView &view,
    const PatchConstructionData &data,
    const KernelPatchPolicy &policy) {
    const KernelPatchPolicy expected = kernel_patch_policy(data);
    if (!same_policy(policy, expected)) {
        throw std::invalid_argument(
            "kernel Riesz patch policy does not match the selected space");
    }

    std::vector<int> global_incidence(view.discrete->nodes.size(), 0);
    for (const Triangle &triangle : view.discrete->elems) {
        for (int node : triangle) ++global_incidence[node];
    }
    std::vector<char> is_dirichlet(view.discrete->nodes.size(), false);
    for (int node : dirichlet_nodes(*view.discrete)) is_dirichlet[node] = true;

    std::vector<KernelRieszPatch> result;
    result.reserve(view.coarse->nodes.size());
    for (int coarse_node = 0;
         coarse_node < static_cast<int>(view.coarse->nodes.size());
         ++coarse_node) {
        KernelRieszPatch patch;
        patch.coarse_node = coarse_node;
        patch.coarse_hat_support = data.hat_support[coarse_node];
        patch.interpolation_support = data.interpolation_support[coarse_node];
        patch.coarse_elements = expand_elements(
            data.graph, patch.coarse_hat_support, policy.patch_layers);
        patch.enlarged_coarse_elements = expand_elements(
            data.graph, patch.coarse_elements, policy.enlargement_layers);
        if (!contains_all(patch.coarse_elements, patch.coarse_hat_support)
            || !contains_all(
                patch.coarse_elements, patch.interpolation_support)) {
            throw std::runtime_error(
                "D_z does not contain all required coarse supports");
        }

        std::vector<char> in_patch(view.coarse->elems.size(), false);
        std::vector<char> in_enlarged(view.coarse->elems.size(), false);
        for (int element : patch.coarse_elements) in_patch[element] = true;
        for (int element : patch.enlarged_coarse_elements)
            in_enlarged[element] = true;
        std::vector<int> local_incidence(view.discrete->nodes.size(), 0);
        for (int element = 0;
             element < static_cast<int>(view.discrete->elems.size()); ++element) {
            const int parent = data.discrete_element_parent[element];
            if (in_enlarged[parent])
                patch.enlarged_discrete_elements.push_back(element);
            if (!in_patch[parent]) continue;
            for (int node : view.discrete->elems[element])
                ++local_incidence[node];
        }
        for (int node = 0;
             node < static_cast<int>(view.discrete->nodes.size()); ++node) {
            if (!is_dirichlet[node] && local_incidence[node] > 0
                && local_incidence[node] == global_incidence[node]) {
                patch.discrete_dofs.push_back(node);
            }
        }

        const ConstraintReduction constraints = reduce_constraints(
            *view.interpolation, patch.discrete_dofs);
        patch.active_constraint_rows = constraints.active_rows;
        patch.independent_constraint_rows = constraints.independent_rows;
        patch.constraints = constraints.matrix;
        result.push_back(std::move(patch));
    }
    return result;
}

ConstraintReduction reduce_constraints(
    const Eigen::SparseMatrix<double> &interpolation,
    const std::vector<int> &discrete_dofs) {
    Eigen::MatrixXd all = Eigen::MatrixXd::Zero(
        interpolation.rows(), discrete_dofs.size());
    for (int local = 0; local < static_cast<int>(discrete_dofs.size()); ++local) {
        const int discrete = discrete_dofs[local];
        if (discrete < 0 || discrete >= interpolation.cols())
            throw std::out_of_range("local Riesz degree of freedom is out of range");
        for (Eigen::SparseMatrix<double>::InnerIterator it(
                 interpolation, discrete); it; ++it)
            all(it.row(), local) = it.value();
    }

    ConstraintReduction result;
    for (int row = 0; row < all.rows(); ++row) {
        if (all.row(row).norm() > kSupportTolerance) result.active_rows.push_back(row);
    }
    if (result.active_rows.empty()) {
        result.matrix.resize(0, discrete_dofs.size());
        return result;
    }

    Eigen::MatrixXd active(result.active_rows.size(), discrete_dofs.size());
    for (int row = 0; row < static_cast<int>(result.active_rows.size()); ++row)
        active.row(row) = all.row(result.active_rows[row]);
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(active.transpose());
    qr.setThreshold(kRankTolerance);
    const int rank = qr.rank();
    const Eigen::VectorXi permutation = qr.colsPermutation().indices();
    result.matrix.resize(rank, discrete_dofs.size());
    for (int row = 0; row < rank; ++row) {
        const int active_row = permutation(row);
        const int original_row = result.active_rows[active_row];
        result.independent_rows.push_back(original_row);
        result.matrix.row(row) = all.row(original_row);
    }
    return result;
}

Eigen::SparseMatrix<double> energy_matrix(const HelmholtzOperators &operators) {
    Eigen::SparseMatrix<double> result = operators.stiffness;
    result += operators.wavenumber * operators.wavenumber * operators.mass;
    result.makeCompressed();
    return result;
}

Eigen::SparseMatrix<double> restrict_sparse_matrix(
    const Eigen::SparseMatrix<double> &matrix,
    const std::vector<int> &dofs) {
    std::vector<int> local_index(matrix.rows(), -1);
    for (int local = 0; local < static_cast<int>(dofs.size()); ++local) {
        if (dofs[local] < 0 || dofs[local] >= matrix.rows())
            throw std::out_of_range("matrix restriction degree of freedom is out of range");
        local_index[dofs[local]] = local;
    }
    std::vector<Eigen::Triplet<double>> triplets;
    for (int local_column = 0; local_column < static_cast<int>(dofs.size()); ++local_column) {
        const int global_column = dofs[local_column];
        for (Eigen::SparseMatrix<double>::InnerIterator it(matrix, global_column); it; ++it) {
            const int local_row = local_index[it.row()];
            if (local_row >= 0)
                triplets.emplace_back(local_row, local_column, it.value());
        }
    }
    Eigen::SparseMatrix<double> result(dofs.size(), dofs.size());
    result.setFromTriplets(triplets.begin(), triplets.end());
    result.makeCompressed();
    return result;
}

Eigen::MatrixXd restrict_dense_matrix(
    const Eigen::SparseMatrix<double> &matrix,
    const std::vector<int> &dofs) {
    return Eigen::MatrixXd(restrict_sparse_matrix(matrix, dofs));
}

ComplexVector restrict_vector(
    const ComplexVector &values,
    const std::vector<int> &dofs) {
    ComplexVector result(dofs.size());
    for (int local = 0; local < static_cast<int>(dofs.size()); ++local)
        result(local) = values(dofs[local]);
    return result;
}

Eigen::MatrixXd kernel_basis(const Eigen::MatrixXd &constraints) {
    const int ambient_dimension = constraints.cols();
    if (ambient_dimension == 0) return Eigen::MatrixXd(0, 0);
    if (constraints.rows() == 0)
        return Eigen::MatrixXd::Identity(ambient_dimension, ambient_dimension);
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(constraints, Eigen::ComputeFullV);
    const Eigen::VectorXd singular_values = svd.singularValues();
    const double scale = singular_values.size() == 0 ? 0.0 : singular_values(0);
    const double threshold = kRankTolerance
        * std::max(constraints.rows(), constraints.cols()) * scale;
    int rank = 0;
    for (double singular_value : singular_values) {
        if (singular_value > threshold) ++rank;
    }
    return svd.matrixV().rightCols(ambient_dimension - rank);
}

LocalKernelRieszResult finalize_local_riesz(
    const Eigen::SparseMatrix<double> &energy,
    const Eigen::MatrixXd &constraints,
    const ComplexVector &rhs,
    const ComplexVector &local_values,
    const ComplexVector &lagrange_multipliers,
    const int kernel_dimension) {
    LocalKernelRieszResult result;
    const int unknowns = rhs.size();
    const int constraint_count = constraints.rows();
    result.kernel_dimension = kernel_dimension;
    result.saddle_unknowns = unknowns + constraint_count;
    result.local_values = local_values;
    result.lagrange_multipliers = lagrange_multipliers;

    const ComplexVector constraint_residual =
        constraints.cast<Complex>() * result.local_values;
    const ComplexVector stationarity =
        energy.cast<Complex>() * result.local_values
        + constraints.transpose().cast<Complex>() * result.lagrange_multipliers
        - rhs;
    result.constraint_relative_residual = constraint_residual.norm()
        / std::max(1.0, constraints.norm() * result.local_values.norm());
    result.riesz_relative_residual = stationarity.norm() / std::max(1.0, rhs.norm());

    const ComplexVector energy_times_solution =
        energy.cast<Complex>() * result.local_values;
    const double energy_squared =
        std::real(result.local_values.dot(energy_times_solution));
    const double residual_action = std::real(result.local_values.dot(rhs));
    if (energy_squared < -1e-11 * std::max(1.0, std::abs(residual_action)))
        throw std::runtime_error("kernel Riesz energy is negative");
    result.eta_squared = std::max(0.0, energy_squared);
    result.eta = std::sqrt(result.eta_squared);
    result.energy_identity_relative_error = std::abs(energy_squared - residual_action)
        / std::max({1.0, std::abs(energy_squared), std::abs(residual_action)});
    return result;
}

std::vector<LocalKernelRieszResult> solve_local_riesz_columns(
    const Eigen::SparseMatrix<double> &energy,
    const Eigen::MatrixXd &constraints,
    const ComplexMatrix &right_hand_sides,
    const KernelRieszSolver solver) {
    if (energy.rows() != energy.cols()
        || energy.rows() != right_hand_sides.rows()
        || constraints.cols() != right_hand_sides.rows()) {
        throw std::invalid_argument("local Riesz dimensions are inconsistent");
    }
    const int unknowns = right_hand_sides.rows();
    const int constraint_count = constraints.rows();
    const int column_count = right_hand_sides.cols();
    ComplexMatrix local_values = ComplexMatrix::Zero(unknowns, column_count);
    ComplexMatrix multipliers = ComplexMatrix::Zero(
        constraint_count, column_count);
    int kernel_dimension = unknowns - constraint_count;

    if (unknowns > 0 && solver == KernelRieszSolver::KernelBasisReference) {
        const Eigen::MatrixXd dense_energy(energy);
        const Eigen::MatrixXd basis = kernel_basis(constraints);
        kernel_dimension = basis.cols();
        if (basis.cols() > 0) {
            const Eigen::MatrixXd reduced_energy = basis.transpose() * dense_energy * basis;
            Eigen::LLT<Eigen::MatrixXd> factorization(reduced_energy);
            if (factorization.info() != Eigen::Success)
                throw std::runtime_error("kernel-basis Riesz factorization failed");
            const ComplexMatrix reduced_rhs =
                basis.transpose().cast<Complex>() * right_hand_sides;
            ComplexMatrix coefficients(basis.cols(), column_count);
            coefficients.real() = factorization.solve(reduced_rhs.real());
            coefficients.imag() = factorization.solve(reduced_rhs.imag());
            if (factorization.info() != Eigen::Success || !coefficients.allFinite())
                throw std::runtime_error("kernel-basis Riesz solve failed");
            local_values = basis.cast<Complex>() * coefficients;
        }
        if (constraint_count > 0) {
            const ComplexMatrix imbalance = right_hand_sides
                - energy.cast<Complex>() * local_values;
            Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(
                constraints.transpose());
            multipliers.real() = qr.solve(imbalance.real());
            multipliers.imag() = qr.solve(imbalance.imag());
        }
    } else if (unknowns > 0) {
        std::vector<ComplexTriplet> triplets;
        triplets.reserve(static_cast<std::size_t>(
            energy.nonZeros() + 2 * constraints.size()));
        for (int column = 0; column < energy.outerSize(); ++column) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(energy, column); it; ++it)
                triplets.emplace_back(it.row(), it.col(), it.value());
        }
        for (int row = 0; row < constraint_count; ++row) {
            for (int column = 0; column < unknowns; ++column) {
                if (constraints(row, column) == 0.0) continue;
                triplets.emplace_back(unknowns + row, column, constraints(row, column));
                triplets.emplace_back(column, unknowns + row, constraints(row, column));
            }
        }
        ComplexSparseMatrix saddle(unknowns + constraint_count, unknowns + constraint_count);
        saddle.setFromTriplets(triplets.begin(), triplets.end());
        saddle.makeCompressed();
        Eigen::SparseLU<ComplexSparseMatrix> factorization;
        factorization.analyzePattern(saddle);
        factorization.factorize(saddle);
        if (factorization.info() != Eigen::Success)
            throw std::runtime_error("kernel Riesz saddle factorization failed");
        ComplexMatrix saddle_rhs = ComplexMatrix::Zero(
            unknowns + constraint_count, column_count);
        saddle_rhs.topRows(unknowns) = right_hand_sides;
        const ComplexMatrix saddle_solution = factorization.solve(saddle_rhs);
        if (factorization.info() != Eigen::Success || !saddle_solution.allFinite())
            throw std::runtime_error("kernel Riesz saddle solve failed");
        local_values = saddle_solution.topRows(unknowns);
        multipliers = saddle_solution.bottomRows(constraint_count);
    }

    std::vector<LocalKernelRieszResult> result;
    result.reserve(column_count);
    for (int column = 0; column < column_count; ++column) {
        result.push_back(finalize_local_riesz(
            energy, constraints, right_hand_sides.col(column),
            local_values.col(column), multipliers.col(column),
            kernel_dimension));
    }
    return result;
}

LocalKernelRieszResult solve_local_riesz(
    const Eigen::SparseMatrix<double> &energy,
    const Eigen::MatrixXd &constraints,
    const ComplexVector &rhs,
    const KernelRieszSolver solver) {
    ComplexMatrix right_hand_sides(rhs.size(), 1);
    right_hand_sides.col(0) = rhs;
    return solve_local_riesz_columns(
        energy, constraints, right_hand_sides, solver).front();
}

double quantile(const std::vector<double> &sorted, double probability) {
    if (sorted.empty()) return 0.0;
    const double position = probability * static_cast<double>(sorted.size() - 1);
    const std::size_t lower = static_cast<std::size_t>(std::floor(position));
    const std::size_t upper = static_cast<std::size_t>(std::ceil(position));
    const double fraction = position - static_cast<double>(lower);
    return (1.0 - fraction) * sorted[lower] + fraction * sorted[upper];
}

Eigen::MatrixXd energy_orthonormal_basis(
    const Eigen::MatrixXd &basis,
    const Eigen::MatrixXd &energy) {
    if (basis.cols() == 0) return basis;
    const Eigen::MatrixXd gram = basis.transpose() * energy * basis;
    Eigen::LLT<Eigen::MatrixXd> factorization(gram);
    if (factorization.info() != Eigen::Success)
        throw std::runtime_error("kernel energy orthonormalization failed");
    const Eigen::MatrixXd inverse_transpose = factorization.matrixU().solve(
        Eigen::MatrixXd::Identity(gram.rows(), gram.cols()));
    return basis * inverse_transpose;
}

double spectral_norm(const Eigen::MatrixXd &matrix) {
    if (matrix.rows() == 0 || matrix.cols() == 0) return 0.0;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(matrix);
    return svd.singularValues()(0);
}

} // namespace

bool AuditKernelResidualEvidence::matches_eta(double eta) const noexcept {
    return !result_fingerprint_.empty() && std::isfinite(eta)
        && eta == eta_;
}

bool AuditKernelResidualEvidence::matches_result(
    const AuditKernelResidualEstimate &estimate) const {
    return matches_eta(estimate.eta)
        && !result_fingerprint_.empty()
        && result_fingerprint_ == audit_kernel_result_fingerprint(estimate);
}

const char *kernel_riesz_solver_name(KernelRieszSolver solver) {
    return solver == KernelRieszSolver::SaddlePoint
        ? "saddle_point" : "kernel_basis_reference";
}

KernelPatchPolicy audit_kernel_patch_policy(
    const AdaptiveMeshHierarchy &hierarchy) {
    const PatchConstructionData data = patch_construction_data(hierarchy);
    int propagation_layers = 0;
    for (int node = 0; node < static_cast<int>(data.hat_support.size()); ++node) {
        propagation_layers = std::max(
            propagation_layers,
            support_distance(
                data.graph,
                data.hat_support[node],
                data.interpolation_support[node]));
    }
    KernelPatchPolicy result;
    result.interpolation_support_layers = propagation_layers;
    result.patch_layers = propagation_layers + 1;
    result.enlargement_layers = 1;
    result.adjacency = "shared_coarse_vertex";
    result.audit_enlargement = "one_shared_coarse_vertex_layer";
    result.hash = policy_hash(result);
    return result;
}

std::string audit_kernel_residual_context_fingerprint(
    const AdaptiveMeshHierarchy &hierarchy,
    const HelmholtzOperators &audit_operators,
    const KernelPatchPolicy &policy) {
    const KernelPatchPolicy expected = audit_kernel_patch_policy(hierarchy);
    if (policy.interpolation_support_layers != expected.interpolation_support_layers
        || policy.patch_layers != expected.patch_layers
        || policy.enlargement_layers != expected.enlargement_layers
        || policy.adjacency != expected.adjacency
        || policy.audit_enlargement != expected.audit_enlargement
        || policy.hash != expected.hash) {
        throw std::invalid_argument(
            "audit-kernel patch policy does not match the hierarchy");
    }
    return audit_kernel_context_fingerprint(
        hierarchy, audit_operators, policy);
}

std::vector<AuditKernelPatch> build_audit_kernel_patches(
    const AdaptiveMeshHierarchy &hierarchy,
    const KernelPatchPolicy &policy) {
    const KernelPatchPolicy expected = audit_kernel_patch_policy(hierarchy);
    if (policy.interpolation_support_layers != expected.interpolation_support_layers
        || policy.patch_layers != expected.patch_layers
        || policy.enlargement_layers != expected.enlargement_layers
        || policy.adjacency != expected.adjacency
        || policy.audit_enlargement != expected.audit_enlargement
        || policy.hash != expected.hash) {
        throw std::invalid_argument("audit-kernel patch policy does not match the hierarchy");
    }

    const PatchConstructionData data = patch_construction_data(hierarchy);
    const TriMesh &audit = hierarchy.cert_audit_mesh();
    const Eigen::SparseMatrix<double> &interpolation =
        hierarchy.cert_audit_quasi_interpolation();
    std::vector<int> global_incidence(audit.nodes.size(), 0);
    for (const Triangle &triangle : audit.elems) {
        for (int node : triangle) ++global_incidence[node];
    }
    std::vector<char> is_dirichlet(audit.nodes.size(), false);
    for (int node : dirichlet_nodes(audit)) is_dirichlet[node] = true;

    std::vector<AuditKernelPatch> result;
    result.reserve(hierarchy.coarse_mesh().nodes.size());
    for (int coarse_node = 0;
         coarse_node < static_cast<int>(hierarchy.coarse_mesh().nodes.size());
         ++coarse_node) {
        AuditKernelPatch patch;
        patch.coarse_node = coarse_node;
        patch.coarse_hat_support = data.hat_support[coarse_node];
        patch.interpolation_support = data.interpolation_support[coarse_node];
        patch.coarse_elements = expand_elements(
            data.graph, patch.coarse_hat_support, policy.patch_layers);
        patch.enlarged_coarse_elements = expand_elements(
            data.graph, patch.coarse_elements, policy.enlargement_layers);
        if (!contains_all(patch.coarse_elements, patch.coarse_hat_support)
            || !contains_all(patch.coarse_elements, patch.interpolation_support)) {
            throw std::runtime_error(
                "D_z does not contain the coarse hat and interpolation-propagated support");
        }

        std::vector<char> in_patch(hierarchy.coarse_mesh().elems.size(), false);
        std::vector<char> in_enlarged(hierarchy.coarse_mesh().elems.size(), false);
        for (int element : patch.coarse_elements) in_patch[element] = true;
        for (int element : patch.enlarged_coarse_elements) in_enlarged[element] = true;
        std::vector<int> local_incidence(audit.nodes.size(), 0);
        for (int element = 0; element < static_cast<int>(audit.elems.size()); ++element) {
            const int parent = data.discrete_element_parent[element];
            if (in_enlarged[parent]) patch.enlarged_discrete_elements.push_back(element);
            if (!in_patch[parent]) continue;
            for (int node : audit.elems[element]) ++local_incidence[node];
        }
        for (int node = 0; node < static_cast<int>(audit.nodes.size()); ++node) {
            if (!is_dirichlet[node]
                && local_incidence[node] > 0
                && local_incidence[node] == global_incidence[node]) {
                patch.discrete_dofs.push_back(node);
            }
        }

        const ConstraintReduction constraints = reduce_constraints(
            interpolation, patch.discrete_dofs);
        patch.active_constraint_rows = constraints.active_rows;
        patch.independent_constraint_rows = constraints.independent_rows;
        patch.constraints = constraints.matrix;
        result.push_back(std::move(patch));
    }
    return result;
}

KernelPatchPolicy kernel_riesz_patch_policy(
    const ReferenceEpochHierarchy &hierarchy,
    KernelRieszSpace space) {
    const ReferenceEpochSpaceView view = reference_epoch_space_view(
        hierarchy, space);
    return kernel_patch_policy(patch_construction_data(view));
}

std::vector<KernelRieszPatch> build_kernel_riesz_patches(
    const ReferenceEpochHierarchy &hierarchy,
    KernelRieszSpace space,
    const KernelPatchPolicy &policy) {
    const ReferenceEpochSpaceView view = reference_epoch_space_view(
        hierarchy, space);
    const PatchConstructionData data = patch_construction_data(view);
    return build_kernel_riesz_patches_impl(view, data, policy);
}

ReferenceResidualRiesz compute_reference_residual_riesz(
    const ReferenceEpochHierarchy &hierarchy,
    const HelmholtzOperators &reference_operators,
    const ComplexVector &reference_load,
    const ComplexVector &candidate_on_reference,
    double doerfler_theta,
    KernelRieszSolver solver) {
    const int reference_nodes = static_cast<int>(
        hierarchy.reference_mesh().nodes.size());
    if (reference_operators.system.rows() != reference_nodes
        || reference_operators.system.cols() != reference_nodes
        || reference_operators.stiffness.rows() != reference_nodes
        || reference_operators.mass.rows() != reference_nodes
        || reference_load.size() != reference_nodes
        || candidate_on_reference.size() != reference_nodes) {
        throw std::invalid_argument(
            "reference residual Riesz inputs do not match reference_mesh");
    }
    if (!(reference_operators.wavenumber > 0.0)
        || !reference_load.allFinite()
        || !candidate_on_reference.allFinite()) {
        throw std::invalid_argument(
            "reference residual Riesz inputs are invalid");
    }
    if (!(doerfler_theta > 0.0 && doerfler_theta <= 1.0)) {
        throw std::invalid_argument(
            "reference residual Doerfler theta must lie in (0,1]");
    }

    ReferenceResidualRiesz result;
    result.policy = kernel_riesz_patch_policy(
        hierarchy, KernelRieszSpace::ReferenceResidual);
    result.solver = solver;
    result.patches = build_kernel_riesz_patches(
        hierarchy, KernelRieszSpace::ReferenceResidual, result.policy);
    result.global_residual = reference_load
        - reference_operators.system * candidate_on_reference;
    for (int node : reference_operators.dirichlet_nodes)
        result.global_residual(node) = Complex(0.0, 0.0);

    const Eigen::SparseMatrix<double> global_energy =
        energy_matrix(reference_operators);
    const int coarse_nodes = static_cast<int>(
        hierarchy.coarse_mesh().nodes.size());
    result.node_eta.assign(coarse_nodes, 0.0);
    result.node_eta_squared.assign(coarse_nodes, 0.0);
    result.local_results.reserve(result.patches.size());
    for (const KernelRieszPatch &patch : result.patches) {
        const Eigen::SparseMatrix<double> local_energy =
            restrict_sparse_matrix(global_energy, patch.discrete_dofs);
        const ComplexVector local_rhs = restrict_vector(
            result.global_residual, patch.discrete_dofs);
        LocalKernelRieszResult local = solve_local_riesz(
            local_energy, patch.constraints, local_rhs, solver);
        result.node_eta[patch.coarse_node] = local.eta;
        result.node_eta_squared[patch.coarse_node] = local.eta_squared;
        result.local_results.push_back(std::move(local));
    }

    result.element_eta_squared.assign(
        hierarchy.coarse_mesh().elems.size(), 0.0);
    for (const KernelRieszPatch &patch : result.patches) {
        if (patch.coarse_hat_support.empty()) {
            throw std::runtime_error(
                "reference residual coarse node has empty element incidence");
        }
        const double share = result.node_eta_squared[patch.coarse_node]
            / static_cast<double>(patch.coarse_hat_support.size());
        for (int element : patch.coarse_hat_support)
            result.element_eta_squared[element] += share;
    }
    const double local_sum = std::accumulate(
        result.local_results.begin(), result.local_results.end(), 0.0,
        [](double sum, const LocalKernelRieszResult &local) {
            return sum + local.eta_squared;
        });
    const double node_sum = std::accumulate(
        result.node_eta_squared.begin(), result.node_eta_squared.end(), 0.0);
    const double element_sum = std::accumulate(
        result.element_eta_squared.begin(),
        result.element_eta_squared.end(), 0.0);
    result.eta = std::sqrt(std::max(0.0, node_sum));
    result.local_square_sum_relative_error = std::abs(local_sum - node_sum)
        / std::max(1.0, node_sum);
    result.allocation_relative_error = std::abs(element_sum - node_sum)
        / std::max(1.0, node_sum);
    if (node_sum > 0.0) {
        result.marked_elements = mark_doerfler(
            result.element_eta_squared, doerfler_theta);
    }
    return result;
}

AmbientDefectRiesz compute_ambient_defect_riesz(
    const ReferenceEpochHierarchy &hierarchy,
    const HelmholtzOperators &ambient_operators,
    const ComplexSparseMatrix &defect_rhs,
    KernelRieszSolver solver,
    AmbientDefectDetail detail,
    int maximum_parallel_patch_solves) {
    const int ambient_nodes = static_cast<int>(
        hierarchy.ambient_mesh().nodes.size());
    const int input_columns = static_cast<int>(defect_rhs.cols());
    if (ambient_operators.system.rows() != ambient_nodes
        || ambient_operators.system.cols() != ambient_nodes
        || ambient_operators.stiffness.rows() != ambient_nodes
        || ambient_operators.mass.rows() != ambient_nodes
        || defect_rhs.rows() != ambient_nodes
        || input_columns <= 0) {
        throw std::invalid_argument(
            "ambient defect Riesz inputs do not match ambient_mesh");
    }
    if (!(ambient_operators.wavenumber > 0.0)
        || !defect_rhs.coeffs().allFinite()) {
        throw std::invalid_argument("ambient defect Riesz inputs are invalid");
    }

    AmbientDefectRiesz result;
    result.policy = kernel_riesz_patch_policy(
        hierarchy, KernelRieszSpace::AmbientDefect);
    result.solver = solver;
    std::vector<KernelRieszPatch> patches = build_kernel_riesz_patches(
        hierarchy, KernelRieszSpace::AmbientDefect, result.policy);
    result.patch_count = patches.size();
    result.local_details_stored =
        detail == AmbientDefectDetail::StoreLocalDetails;
    if (result.local_details_stored)
        result.local_results.resize(patches.size());
    result.gram = ComplexMatrix::Zero(input_columns, input_columns);
    result.column_eta_squared = Eigen::VectorXd::Zero(input_columns);

    const Eigen::SparseMatrix<double> global_energy =
        energy_matrix(ambient_operators);
    using RowMajorComplexSparse =
        Eigen::SparseMatrix<Complex, Eigen::RowMajor>;
    const RowMajorComplexSparse defect_by_row = defect_rhs;
    std::atomic<bool> failed{false};
    std::exception_ptr first_exception;
    std::mutex exception_mutex;

    struct PatchContribution {
        AmbientDefectLocalRiesz local;
    };
    const auto solve_patch = [&](int patch_index, PatchContribution &output) {
        if (failed.load(std::memory_order_relaxed)) return;
        try {
            const KernelRieszPatch &patch = patches[patch_index];
            std::vector<char> active(input_columns, false);
            for (int ambient_dof : patch.discrete_dofs) {
                for (RowMajorComplexSparse::InnerIterator it(
                         defect_by_row, ambient_dof); it; ++it) {
                    if (it.value() != Complex(0.0, 0.0))
                        active[it.col()] = true;
                }
            }
            for (int column = 0; column < input_columns; ++column) {
                if (active[column]) output.local.active_columns.push_back(column);
            }
            if (output.local.active_columns.empty()) return;

            std::vector<int> compact_column(input_columns, -1);
            for (int local_column = 0;
                 local_column
                     < static_cast<int>(output.local.active_columns.size());
                 ++local_column) {
                compact_column[output.local.active_columns[local_column]] =
                    local_column;
            }
            ComplexMatrix local_right_hand_sides = ComplexMatrix::Zero(
                patch.discrete_dofs.size(),
                output.local.active_columns.size());
            for (int local_row = 0;
                 local_row < static_cast<int>(patch.discrete_dofs.size());
                 ++local_row) {
                const int ambient_dof = patch.discrete_dofs[local_row];
                for (RowMajorComplexSparse::InnerIterator it(
                         defect_by_row, ambient_dof); it; ++it) {
                    const int local_column = compact_column[it.col()];
                    if (local_column >= 0)
                        local_right_hand_sides(local_row, local_column) =
                            it.value();
                }
            }
            const Eigen::SparseMatrix<double> local_energy =
                restrict_sparse_matrix(global_energy, patch.discrete_dofs);
            std::vector<LocalKernelRieszResult> columns =
                solve_local_riesz_columns(
                    local_energy, patch.constraints,
                    local_right_hand_sides, solver);
            ComplexMatrix representatives = ComplexMatrix::Zero(
                patch.discrete_dofs.size(),
                output.local.active_columns.size());
            output.local.columns.reserve(output.local.active_columns.size());
            for (int column = 0;
                 column < static_cast<int>(output.local.active_columns.size());
                 ++column) {
                LocalKernelRieszResult local = std::move(columns[column]);
                representatives.col(column) = local.local_values;
                output.local.columns.push_back(std::move(local));
            }
            output.local.gram = representatives.adjoint()
                * local_energy.cast<Complex>() * representatives;
            output.local.hermitian_relative_error =
                (output.local.gram - output.local.gram.adjoint()).norm()
                / std::max(1.0, output.local.gram.norm());
        } catch (...) {
            failed.store(true, std::memory_order_relaxed);
            std::lock_guard<std::mutex> lock(exception_mutex);
            if (!first_exception) first_exception = std::current_exception();
        }
    };

    if (maximum_parallel_patch_solves < 0) {
        throw std::invalid_argument(
            "maximum parallel patch solves must be nonnegative");
    }
    int maximum_threads = 1;
#ifdef _OPENMP
    maximum_threads = std::max(1, omp_get_max_threads());
#endif
    if (maximum_parallel_patch_solves > 0) {
        maximum_threads = std::min(
            maximum_threads, maximum_parallel_patch_solves);
    }
    maximum_threads = std::min(
        maximum_threads, static_cast<int>(std::max<std::size_t>(1, patches.size())));
    const std::size_t batch_capacity = std::max<std::size_t>(
        1, 2 * static_cast<std::size_t>(maximum_threads));
    for (std::size_t batch_begin = 0;
         batch_begin < patches.size();
         batch_begin += batch_capacity) {
        const std::size_t batch_size = std::min(
            batch_capacity, patches.size() - batch_begin);
        std::vector<PatchContribution> contributions(batch_size);
        const auto patch_solve_begin = std::chrono::steady_clock::now();
#ifdef _OPENMP
        #pragma omp parallel num_threads(maximum_threads)
        {
            #pragma omp single
            result.parallel_threads = omp_get_num_threads();
            #pragma omp for schedule(dynamic, 1)
            for (int local_index = 0;
                 local_index < static_cast<int>(batch_size);
                 ++local_index) {
                solve_patch(
                    static_cast<int>(batch_begin) + local_index,
                    contributions[local_index]);
            }
        }
#else
        for (int local_index = 0;
             local_index < static_cast<int>(batch_size);
             ++local_index) {
            solve_patch(
                static_cast<int>(batch_begin) + local_index,
                contributions[local_index]);
        }
#endif
        result.patch_solve_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - patch_solve_begin).count();
        if (first_exception) std::rethrow_exception(first_exception);

        // Reduce each bounded batch in the original patch order.  Thread
        // scheduling therefore cannot change a floating-point addition, and
        // completed local vectors are released before the next batch.
        const auto gram_reduction_begin = std::chrono::steady_clock::now();
        for (std::size_t local_index = 0;
             local_index < batch_size;
             ++local_index) {
            const std::size_t patch_index = batch_begin + local_index;
            AmbientDefectLocalRiesz &local = contributions[local_index].local;
            result.maximum_active_columns = std::max(
                result.maximum_active_columns, local.active_columns.size());
            if (!local.active_columns.empty()) {
                ++result.patch_factorizations;
                result.right_hand_side_solves += local.active_columns.size();
            }
            const int active_count = static_cast<int>(
                local.active_columns.size());
            // A column owns a disjoint global Gram column, so this reduction
            // is race-free.  Patches remain ordered and a barrier separates
            // them, preserving the exact per-entry floating-point addition
            // order while removing the large serial r_patch^2 scatter.
#ifdef _OPENMP
            #pragma omp parallel for num_threads(maximum_threads) \
                schedule(static) if(active_count >= 256)
#endif
            for (int local_column = 0;
                 local_column < active_count;
                 ++local_column) {
                const int global_column = local.active_columns[local_column];
                result.column_eta_squared(global_column) +=
                    local.columns[local_column].eta_squared;
                for (int local_row = 0;
                     local_row < static_cast<int>(local.active_columns.size());
                     ++local_row) {
                    const int global_row = local.active_columns[local_row];
                    result.gram(global_row, global_column) +=
                        local.gram(local_row, local_column);
                }
            }
            if (result.local_details_stored)
                result.local_results[patch_index] = std::move(local);
        }
        result.gram_reduction_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - gram_reduction_begin).count();
    }

    for (int column = 0; column < input_columns; ++column) {
        const double gram_diagonal = std::real(result.gram(column, column));
        result.local_square_sum_relative_error = std::max(
            result.local_square_sum_relative_error,
            std::abs(gram_diagonal - result.column_eta_squared(column))
                / std::max(1.0, result.column_eta_squared(column)));
    }
    if (result.local_details_stored)
        result.patches = std::move(patches);
    return result;
}

AuditKernelResidualEstimate estimate_audit_kernel_residual(
    const AdaptiveMeshHierarchy &hierarchy,
    const HelmholtzOperators &audit_operators,
    const ComplexVector &audit_load,
    const ComplexVector &candidate_on_audit,
    double doerfler_theta,
    KernelRieszSolver solver) {
    const int audit_nodes = static_cast<int>(hierarchy.cert_audit_mesh().nodes.size());
    if (audit_operators.system.rows() != audit_nodes
        || audit_operators.system.cols() != audit_nodes
        || audit_operators.stiffness.rows() != audit_nodes
        || audit_operators.mass.rows() != audit_nodes
        || audit_load.size() != audit_nodes
        || candidate_on_audit.size() != audit_nodes) {
        throw std::invalid_argument("audit-kernel residual inputs do not match the audit mesh");
    }
    if (!(audit_operators.wavenumber > 0.0)
        || !audit_load.allFinite() || !candidate_on_audit.allFinite()) {
        throw std::invalid_argument("audit-kernel residual inputs are invalid");
    }
    if (!(doerfler_theta > 0.0 && doerfler_theta <= 1.0))
        throw std::invalid_argument("audit-kernel Doerfler theta must lie in (0,1]");

    AuditKernelResidualEstimate result;
    result.policy = audit_kernel_patch_policy(hierarchy);
    result.solver = solver;
    result.patches = build_audit_kernel_patches(hierarchy, result.policy);
    result.global_residual = audit_load
        - audit_operators.system * candidate_on_audit;
    for (int node : audit_operators.dirichlet_nodes)
        result.global_residual(node) = Complex(0.0, 0.0);

    const Eigen::SparseMatrix<double> global_energy = energy_matrix(audit_operators);
    const int coarse_nodes = static_cast<int>(hierarchy.coarse_mesh().nodes.size());
    result.node_eta.assign(coarse_nodes, 0.0);
    result.node_eta_squared.assign(coarse_nodes, 0.0);
    result.local_results.reserve(result.patches.size());
    for (const AuditKernelPatch &patch : result.patches) {
        const Eigen::SparseMatrix<double> local_energy =
            restrict_sparse_matrix(global_energy, patch.discrete_dofs);
        const ComplexVector local_rhs = restrict_vector(
            result.global_residual, patch.discrete_dofs);
        LocalKernelRieszResult local = solve_local_riesz(
            local_energy, patch.constraints, local_rhs, solver);
        result.node_eta[patch.coarse_node] = local.eta;
        result.node_eta_squared[patch.coarse_node] = local.eta_squared;
        result.local_results.push_back(std::move(local));
    }

    const int coarse_elements = static_cast<int>(hierarchy.coarse_mesh().elems.size());
    result.element_eta_squared.assign(coarse_elements, 0.0);
    for (const AuditKernelPatch &patch : result.patches) {
        if (patch.coarse_hat_support.empty())
            throw std::runtime_error("coarse node has empty element incidence");
        const double share = result.node_eta_squared[patch.coarse_node]
            / static_cast<double>(patch.coarse_hat_support.size());
        for (int element : patch.coarse_hat_support)
            result.element_eta_squared[element] += share;
    }
    const double node_sum = std::accumulate(
        result.node_eta_squared.begin(), result.node_eta_squared.end(), 0.0);
    const double element_sum = std::accumulate(
        result.element_eta_squared.begin(), result.element_eta_squared.end(), 0.0);
    result.eta = std::sqrt(std::max(0.0, node_sum));
    result.allocation_relative_error = node_sum == 0.0
        ? 0.0
        : std::abs(element_sum - node_sum) / node_sum;
    if (node_sum > 0.0)
        result.marked_elements = mark_doerfler(result.element_eta_squared, doerfler_theta);

    result.evidence.level_ = KernelResidualEvidenceLevel::Diagnostic;
    result.evidence.source_ =
        "ordinary floating-point audit-kernel Riesz computation";
    result.evidence.backend_ =
        std::string("Eigen/") + kernel_riesz_solver_name(solver);
    result.evidence.context_fingerprint_ = audit_kernel_context_fingerprint(
        hierarchy, audit_operators, result.policy);
    result.evidence.eta_ = result.eta;
    result.evidence.result_fingerprint_ =
        audit_kernel_result_fingerprint(result);
    result.evidence.diagnostic_fingerprint_ =
        audit_kernel_diagnostic_fingerprint(
            result.evidence.context_fingerprint_, solver, doerfler_theta,
            audit_load, candidate_on_audit, result);
    result.evidence.failure_reason_ =
        "ordinary Eigen Riesz solves have no directed-rounding enclosure";
    return result;
}

LocalKernelEffectivityReport evaluate_local_kernel_effectivity(
    const AdaptiveMeshHierarchy &hierarchy,
    const HelmholtzOperators &audit_operators,
    const AuditKernelResidualEstimate &estimate,
    const ComplexVector &certification_solution,
    const ComplexVector &candidate_on_audit) {
    const int audit_nodes = static_cast<int>(hierarchy.cert_audit_mesh().nodes.size());
    if (certification_solution.size() != audit_nodes
        || candidate_on_audit.size() != audit_nodes
        || estimate.patches.size() != estimate.local_results.size()) {
        throw std::invalid_argument("local effectivity inputs are inconsistent");
    }
    LocalKernelEffectivityReport report;
    std::vector<double> ratios;
    report.local.reserve(estimate.patches.size());
    const double global_scale = std::max({
        1.0,
        discrete_energy_norm(audit_operators, certification_solution),
        discrete_energy_norm(audit_operators, candidate_on_audit)});
    const double zero_tolerance = 64.0 * std::numeric_limits<double>::epsilon()
        * global_scale;
    for (int index = 0; index < static_cast<int>(estimate.patches.size()); ++index) {
        LocalKernelEffectivity entry;
        entry.coarse_node = estimate.patches[index].coarse_node;
        entry.eta = estimate.local_results[index].eta;
        entry.enlarged_patch_error = compute_local_discrete_helmholtz_error(
            hierarchy.cert_audit_mesh(),
            audit_operators,
            certification_solution,
            candidate_on_audit,
            estimate.patches[index].enlarged_discrete_elements).energy;
        if (entry.enlarged_patch_error > zero_tolerance) {
            entry.ratio = entry.eta / entry.enlarged_patch_error;
            entry.valid = std::isfinite(entry.ratio);
            if (entry.valid) ratios.push_back(entry.ratio);
        } else {
            ++report.distribution.excluded_zero_error;
        }
        report.local.push_back(entry);
    }
    std::sort(ratios.begin(), ratios.end());
    report.distribution.count = ratios.size();
    if (!ratios.empty()) {
        report.distribution.minimum = ratios.front();
        report.distribution.first_quartile = quantile(ratios, 0.25);
        report.distribution.median = quantile(ratios, 0.5);
        report.distribution.third_quartile = quantile(ratios, 0.75);
        report.distribution.percentile_90 = quantile(ratios, 0.9);
        report.distribution.maximum = ratios.back();
    }
    return report;
}

ExplicitGlobalKernelComparison compare_with_explicit_global_kernel(
    const AdaptiveMeshHierarchy &hierarchy,
    const HelmholtzOperators &audit_operators,
    const AuditKernelResidualEstimate &estimate,
    int maximum_free_dofs) {
    if (maximum_free_dofs <= 0)
        throw std::invalid_argument("explicit global-kernel DOF limit must be positive");
    const int audit_nodes = static_cast<int>(hierarchy.cert_audit_mesh().nodes.size());
    if (estimate.global_residual.size() != audit_nodes
        || estimate.patches.size() != estimate.local_results.size()) {
        throw std::invalid_argument("explicit global-kernel comparison inputs are inconsistent");
    }

    std::vector<char> is_dirichlet(audit_nodes, false);
    for (int node : audit_operators.dirichlet_nodes) is_dirichlet[node] = true;
    std::vector<int> free_dofs;
    std::vector<int> free_index(audit_nodes, -1);
    for (int node = 0; node < audit_nodes; ++node) {
        if (is_dirichlet[node]) continue;
        free_index[node] = static_cast<int>(free_dofs.size());
        free_dofs.push_back(node);
    }
    if (static_cast<int>(free_dofs.size()) > maximum_free_dofs)
        throw std::invalid_argument("explicit global-kernel comparison exceeds its DOF limit");

    const Eigen::SparseMatrix<double> sparse_energy = energy_matrix(audit_operators);
    const Eigen::SparseMatrix<double> global_sparse_energy =
        restrict_sparse_matrix(sparse_energy, free_dofs);
    const Eigen::MatrixXd global_energy(global_sparse_energy);
    const ConstraintReduction global_constraints = reduce_constraints(
        hierarchy.cert_audit_quasi_interpolation(), free_dofs);
    const Eigen::MatrixXd global_basis = kernel_basis(global_constraints.matrix);
    const ComplexVector free_residual = restrict_vector(
        estimate.global_residual, free_dofs);
    const LocalKernelRieszResult global_riesz = solve_local_riesz(
        global_sparse_energy,
        global_constraints.matrix,
        free_residual,
        KernelRieszSolver::KernelBasisReference);

    ExplicitGlobalKernelComparison result;
    result.global_kernel_dimension = global_basis.cols();
    result.global_eta = global_riesz.eta;
    result.localized_eta = estimate.eta;
    if (global_basis.cols() == 0) {
        result.lower_bound = result.localized_eta;
        result.upper_bound = result.localized_eta;
        const double tolerance = 2e-11;
        result.lower_direction_holds = result.lower_bound <= result.global_eta + tolerance;
        result.upper_direction_holds = result.global_eta <= result.upper_bound + tolerance;
        return result;
    }

    std::vector<Eigen::MatrixXd> local_bases;
    int total_columns = 0;
    for (const AuditKernelPatch &patch : estimate.patches) {
        const Eigen::MatrixXd basis = kernel_basis(patch.constraints);
        if (basis.cols() == 0) {
            local_bases.emplace_back(patch.discrete_dofs.size(), 0);
            continue;
        }
        const Eigen::MatrixXd local_energy = restrict_dense_matrix(
            sparse_energy, patch.discrete_dofs);
        Eigen::MatrixXd normalized = energy_orthonormal_basis(basis, local_energy);
        total_columns += normalized.cols();
        local_bases.push_back(std::move(normalized));
    }
    result.local_kernel_columns = total_columns;
    if (total_columns == 0) {
        result.overlap_constant = std::numeric_limits<double>::infinity();
        result.stable_decomposition_constant = std::numeric_limits<double>::infinity();
        result.upper_bound = std::numeric_limits<double>::infinity();
        return result;
    }

    Eigen::MatrixXd synthesis = Eigen::MatrixXd::Zero(free_dofs.size(), total_columns);
    int first_column = 0;
    for (int patch_index = 0;
         patch_index < static_cast<int>(estimate.patches.size());
         ++patch_index) {
        const Eigen::MatrixXd &basis = local_bases[patch_index];
        for (int local_node = 0;
             local_node < static_cast<int>(estimate.patches[patch_index].discrete_dofs.size());
             ++local_node) {
            const int row = free_index[estimate.patches[patch_index].discrete_dofs[local_node]];
            if (row < 0)
                throw std::runtime_error("local audit kernel contains a Dirichlet degree of freedom");
            synthesis.block(row, first_column, 1, basis.cols()) = basis.row(local_node);
        }
        first_column += basis.cols();
    }

    Eigen::MatrixXd synthesis_gram = synthesis.transpose()
        * global_energy * synthesis;
    synthesis_gram = 0.5 * (synthesis_gram + synthesis_gram.transpose());
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> overlap_solver(synthesis_gram);
    if (overlap_solver.info() != Eigen::Success)
        throw std::runtime_error("explicit overlap eigenproblem failed");
    result.overlap_constant = std::sqrt(std::max(
        0.0, overlap_solver.eigenvalues().maxCoeff()));

    const Eigen::MatrixXd normalized_global = energy_orthonormal_basis(
        global_basis, global_energy);
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> decomposition(synthesis);
    const Eigen::MatrixXd coefficients = decomposition.solve(normalized_global);
    result.span_relative_residual = (synthesis * coefficients - normalized_global).norm()
        / std::max(1.0, normalized_global.norm());
    result.stable_decomposition_constant = spectral_norm(coefficients);
    result.lower_bound = result.localized_eta
        / std::max(result.overlap_constant, std::numeric_limits<double>::min());
    result.upper_bound = result.stable_decomposition_constant * result.localized_eta;
    const double tolerance = 2e-10 * std::max({
        1.0, result.global_eta, result.lower_bound,
        std::isfinite(result.upper_bound) ? result.upper_bound : 1.0});
    result.lower_direction_holds = result.lower_bound <= result.global_eta + tolerance;
    result.upper_direction_holds = result.span_relative_residual <= 2e-10
        && result.global_eta <= result.upper_bound + tolerance;
    return result;
}

} // namespace lod2d::helmholtz::adaptive

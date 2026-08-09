#include "helmholtz/patch_system.h"
#include "helmholtz/boundary.h"

#include <Eigen/QR>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace lod2d::helmholtz {
namespace {

struct PatchWorkspace {
    std::vector<int> node_count;
    std::vector<int> node_seen;
    std::vector<int> local_index;
    std::vector<int> local_seen;
    std::vector<int> coarse_row_index;
    std::vector<int> coarse_row_seen;
    int stamp = 0;

    void begin(std::size_t fine_nodes, std::size_t coarse_nodes) {
        bool resized = false;
        if (node_count.size() != fine_nodes) {
            node_count.assign(fine_nodes, 0);
            node_seen.assign(fine_nodes, 0);
            local_index.assign(fine_nodes, -1);
            local_seen.assign(fine_nodes, 0);
            resized = true;
        }
        if (coarse_row_index.size() != coarse_nodes) {
            coarse_row_index.assign(coarse_nodes, -1);
            coarse_row_seen.assign(coarse_nodes, 0);
            resized = true;
        }
        if (resized) {
            stamp = 0;
            std::fill(node_seen.begin(), node_seen.end(), 0);
            std::fill(local_seen.begin(), local_seen.end(), 0);
            std::fill(coarse_row_seen.begin(), coarse_row_seen.end(), 0);
        }
        if (stamp == std::numeric_limits<int>::max()) {
            std::fill(node_seen.begin(), node_seen.end(), 0);
            std::fill(local_seen.begin(), local_seen.end(), 0);
            std::fill(coarse_row_seen.begin(), coarse_row_seen.end(), 0);
            stamp = 0;
        }
        ++stamp;
    }
};

thread_local PatchWorkspace patch_workspace;

std::uint64_t edge_key(int a, int b) {
    const auto lo = static_cast<std::uint32_t>(std::min(a, b));
    const auto hi = static_cast<std::uint32_t>(std::max(a, b));
    return (static_cast<std::uint64_t>(lo) << 32U) | hi;
}

std::vector<std::vector<int>> fine_element_children(
    const Eigen::SparseMatrix<double> &prolongation,
    int coarse_element_count) {
    if (prolongation.cols() != coarse_element_count)
        throw std::invalid_argument("fine element prolongation has the wrong column count");
    std::vector<std::vector<int>> children(coarse_element_count);
    for (int coarse_element = 0; coarse_element < prolongation.outerSize(); ++coarse_element) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(
                 prolongation, coarse_element); it; ++it) {
            if (it.value() != 0.0) children[coarse_element].push_back(it.row());
        }
    }
    return children;
}

std::vector<int> node_incidence(const TriMesh &mesh) {
    std::vector<int> incidence(mesh.nodes.size(), 0);
    for (const Triangle &triangle : mesh.elems)
        for (int vertex : triangle) ++incidence[vertex];
    return incidence;
}

std::unordered_map<std::uint64_t, int> edge_counts(const TriMesh &mesh) {
    std::unordered_map<std::uint64_t, int> counts;
    counts.reserve(static_cast<std::size_t>(3 * mesh.elems.size()));
    for (const Triangle &triangle : mesh.elems) {
        ++counts[edge_key(triangle[0], triangle[1])];
        ++counts[edge_key(triangle[1], triangle[2])];
        ++counts[edge_key(triangle[2], triangle[0])];
    }
    return counts;
}

template <typename LocalIndex>
Eigen::SparseMatrix<double> restrict_matrix(
    const Eigen::SparseMatrix<double> &global,
    const std::vector<int> &local_vertices,
    LocalIndex local_index_of) {
    std::vector<Eigen::Triplet<double>> triplets;
    for (int local_col = 0; local_col < static_cast<int>(local_vertices.size()); ++local_col) {
        const int global_col = local_vertices[local_col];
        for (Eigen::SparseMatrix<double>::InnerIterator it(global, global_col); it; ++it) {
            const int local_row = local_index_of(it.row());
            if (local_row >= 0 && it.value() != 0.0)
                triplets.emplace_back(local_row, local_col, it.value());
        }
    }
    Eigen::SparseMatrix<double> local(
        static_cast<int>(local_vertices.size()),
        static_cast<int>(local_vertices.size()));
    local.setFromTriplets(triplets.begin(), triplets.end());
    local.makeCompressed();
    return local;
}


struct LocalPatchLevel {
    std::vector<int> elements;
    std::vector<int> vertices;
    std::vector<int> global_to_local;
};

LocalPatchLevel build_local_patch_level(
    const TriMesh &mesh,
    std::vector<int> patch_elements,
    const std::vector<int> &global_incidence,
    const std::vector<int> *forced_vertices = nullptr) {
    LocalPatchLevel level;
    level.elements = std::move(patch_elements);
    level.global_to_local.assign(mesh.nodes.size(), -1);
    if (forced_vertices) {
        level.vertices = *forced_vertices;
    } else {
        std::vector<int> counts(mesh.nodes.size(), 0);
        std::vector<int> touched;
        touched.reserve(3 * level.elements.size());
        for (int element : level.elements) {
            for (int vertex : mesh.elems[element]) {
                if (counts[vertex]++ == 0) touched.push_back(vertex);
            }
        }
        for (int vertex : touched) {
            if (counts[vertex] == global_incidence[vertex]
                && !std::binary_search(
                    mesh.dirichlet.begin(), mesh.dirichlet.end(), vertex))
                level.vertices.push_back(vertex);
        }
    }
    for (int local = 0; local < static_cast<int>(level.vertices.size()); ++local)
        level.global_to_local[level.vertices[local]] = local;
    return level;
}

std::vector<int> prolong_patch_elements(
    const Eigen::SparseMatrix<double> &element_prolongation,
    const std::vector<int> &coarse_elements) {
    std::vector<int> fine_elements;
    for (int coarse : coarse_elements) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(
                 element_prolongation, coarse); it; ++it) {
            if (it.value() != 0.0) fine_elements.push_back(it.row());
        }
    }
    return fine_elements;
}

Eigen::SparseMatrix<double> restrict_prolongation(
    const Eigen::SparseMatrix<double> &global,
    const LocalPatchLevel &coarse,
    const LocalPatchLevel &fine) {
    std::vector<Eigen::Triplet<double>> triplets;
    for (int local_col = 0; local_col < static_cast<int>(coarse.vertices.size());
         ++local_col) {
        const int global_col = coarse.vertices[local_col];
        for (Eigen::SparseMatrix<double>::InnerIterator it(global, global_col);
             it; ++it) {
            const int local_row = fine.global_to_local[it.row()];
            if (local_row >= 0 && it.value() != 0.0)
                triplets.emplace_back(local_row, local_col, it.value());
        }
    }
    Eigen::SparseMatrix<double> local(
        static_cast<int>(fine.vertices.size()),
        static_cast<int>(coarse.vertices.size()));
    local.setFromTriplets(triplets.begin(), triplets.end());
    local.makeCompressed();
    return local;
}
} // namespace

HelmholtzPatchAssembler::HelmholtzPatchAssembler(
    const TriMesh &coarse,
    const TriMesh &fine,
    const Eigen::SparseMatrix<double> &fine_element_prolongation,
    const Eigen::SparseMatrix<double> &fine_dg_prolongation,
    const Eigen::SparseMatrix<double> &quasi_interpolation,
    const Eigen::SparseMatrix<double> &patches,
    const std::vector<TriMesh> &hierarchy_meshes,
    const std::vector<Eigen::SparseMatrix<double>> &node_level_prolongations,
    const std::vector<Eigen::SparseMatrix<double>> &element_level_prolongations,
    const HelmholtzOperators &operators)
    : coarse_(coarse),
      fine_(fine),
      fine_dg_prolongation_(fine_dg_prolongation),
      quasi_interpolation_(quasi_interpolation),
      patches_(patches),
      operators_(operators),
      children_(fine_element_children(
          fine_element_prolongation, static_cast<int>(coarse.elems.size()))),
      fine_incidence_(node_incidence(fine)),
      fine_dirichlet_(fine.nodes.size(), false),
      fine_edge_counts_(edge_counts(fine)),
      hierarchy_meshes_(hierarchy_meshes),
      node_level_prolongations_(node_level_prolongations),
      element_level_prolongations_(element_level_prolongations) {
    for (int node : dirichlet_nodes(fine_)) {
        if (node >= 0 && node < static_cast<int>(fine_dirichlet_.size()))
            fine_dirichlet_[node] = true;
    }
    const int coarse_element_count = static_cast<int>(coarse_.elems.size());
    if (patches_.rows() != coarse_element_count
        || patches_.cols() != coarse_element_count)
        throw std::invalid_argument("Helmholtz patch matrix has the wrong dimensions");
    if (fine_dg_prolongation_.rows() != 3 * static_cast<int>(fine_.elems.size())
        || fine_dg_prolongation_.cols() != 3 * coarse_element_count)
        throw std::invalid_argument("Helmholtz DG prolongation has the wrong dimensions");
    if (quasi_interpolation_.rows() != static_cast<int>(coarse_.nodes.size())
        || quasi_interpolation_.cols() != static_cast<int>(fine_.nodes.size()))
        throw std::invalid_argument("Helmholtz quasi-interpolation has the wrong dimensions");
    if (operators_.element_blocks.size() != fine_.elems.size())
        throw std::invalid_argument("Helmholtz element block count must match fine elements");
    if (operators_.system.rows() != static_cast<int>(fine_.nodes.size()))
        throw std::invalid_argument("Helmholtz operator size must match fine nodes");
    if (!hierarchy_meshes_.empty()) {
        if (hierarchy_meshes_.size() != node_level_prolongations_.size() + 1
            || hierarchy_meshes_.size() != element_level_prolongations_.size() + 1)
            throw std::invalid_argument("inconsistent Helmholtz geometric hierarchy sizes");
        if (hierarchy_meshes_.front().nodes.size() != coarse_.nodes.size()
            || hierarchy_meshes_.front().elems.size() != coarse_.elems.size()
            || hierarchy_meshes_.back().nodes.size() != fine_.nodes.size()
            || hierarchy_meshes_.back().elems.size() != fine_.elems.size())
            throw std::invalid_argument("Helmholtz geometric hierarchy endpoints do not match");
        hierarchy_incidence_.reserve(hierarchy_meshes_.size());
        for (const TriMesh &mesh : hierarchy_meshes_)
            hierarchy_incidence_.push_back(node_incidence(mesh));
    } else if (!node_level_prolongations_.empty()
               || !element_level_prolongations_.empty()) {
        throw std::invalid_argument("Helmholtz hierarchy prolongations require meshes");
    }
}

std::size_t HelmholtzPatchAssembler::patch_cost(int target) const {
    if (target < 0 || target >= patch_count())
        throw std::out_of_range("Helmholtz patch target is out of range");
    std::size_t cost = 0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(patches_, target); it; ++it) {
        if (it.value() != 0.0) cost += children_[it.row()].size();
    }
    return cost;
}

HelmholtzPatchSystem HelmholtzPatchAssembler::assemble(int target) const {
    if (target < 0 || target >= patch_count())
        throw std::out_of_range("Helmholtz patch target is out of range");

    HelmholtzPatchSystem system;
    system.target_element = target;
    system.wavenumber = operators_.wavenumber;
    for (Eigen::SparseMatrix<double>::InnerIterator it(patches_, target); it; ++it) {
        if (it.value() == 0.0) continue;
        const auto &element_children = children_[it.row()];
        system.patch_elements.insert(
            system.patch_elements.end(), element_children.begin(), element_children.end());
    }
    if (system.patch_elements.empty())
        throw std::runtime_error("Helmholtz corrector patch is empty");

    patch_workspace.begin(fine_.nodes.size(), coarse_.nodes.size());
    const int stamp = patch_workspace.stamp;
    std::vector<int> touched_vertices;
    touched_vertices.reserve(3 * system.patch_elements.size());
    Point2 lower = fine_.nodes[fine_.elems[system.patch_elements.front()][0]];
    Point2 upper = lower;
    for (int element : system.patch_elements) {
        for (int vertex : fine_.elems[element]) {
            lower = lower.cwiseMin(fine_.nodes[vertex]);
            upper = upper.cwiseMax(fine_.nodes[vertex]);
            if (patch_workspace.node_seen[vertex] != stamp) {
                patch_workspace.node_seen[vertex] = stamp;
                patch_workspace.node_count[vertex] = 0;
                touched_vertices.push_back(vertex);
            }
            ++patch_workspace.node_count[vertex];
        }
    }
    system.diameter = (upper - lower).norm();

    system.local_vertices.reserve(touched_vertices.size());
    for (int vertex : touched_vertices) {
        if (patch_workspace.node_count[vertex] != fine_incidence_[vertex]) continue;
        if (fine_dirichlet_[vertex]) continue;
        patch_workspace.local_seen[vertex] = stamp;
        patch_workspace.local_index[vertex]
            = static_cast<int>(system.local_vertices.size());
        system.local_vertices.push_back(vertex);
    }
    if (system.local_vertices.empty())
        throw std::runtime_error("Helmholtz corrector patch has no unconstrained fine vertices");
    auto local_index_of = [&](int vertex) {
        return patch_workspace.local_seen[vertex] == stamp
            ? patch_workspace.local_index[vertex]
            : -1;
    };
    if (!hierarchy_meshes_.empty()) {
        std::vector<int> coarse_patch_elements;
        for (Eigen::SparseMatrix<double>::InnerIterator it(patches_, target);
             it; ++it) {
            if (it.value() != 0.0) coarse_patch_elements.push_back(it.row());
        }
        LocalPatchLevel previous = build_local_patch_level(
            hierarchy_meshes_.front(), std::move(coarse_patch_elements),
            hierarchy_incidence_.front());
        for (int level = 0;
             level < static_cast<int>(node_level_prolongations_.size());
             ++level) {
            std::vector<int> next_elements = prolong_patch_elements(
                element_level_prolongations_[level], previous.elements);
            const bool finest =
                level + 1 == static_cast<int>(node_level_prolongations_.size());
            LocalPatchLevel next = build_local_patch_level(
                hierarchy_meshes_[level + 1], std::move(next_elements),
                hierarchy_incidence_[level + 1],
                finest ? &system.local_vertices : nullptr);
            if (!previous.vertices.empty()) {
                if (next.vertices.empty())
                    throw std::runtime_error(
                        "nonempty coarse patch space has an empty fine space");
                system.geometric_prolongations.push_back(
                    restrict_prolongation(
                        node_level_prolongations_[level], previous, next));
            }
            previous = std::move(next);
        }
        if (previous.vertices != system.local_vertices)
            throw std::runtime_error(
                "geometric hierarchy finest patch DOFs do not match the operator");
    }


    static constexpr int local_edges[3][2] = {{0, 1}, {1, 2}, {2, 0}};
    for (int element : system.patch_elements) {
        const Triangle &triangle = fine_.elems[element];
        for (const auto &edge : local_edges) {
            const Edge physical_edge = canonical_edge(
                triangle[edge[0]], triangle[edge[1]]);
            if (fine_edge_counts_.at(edge_key(
                    physical_edge[0], physical_edge[1])) == 1
                && boundary_tag(fine_, physical_edge) == BoundaryTag::Robin) {
                system.touches_physical_boundary = true;
                break;
            }
        }
        if (system.touches_physical_boundary) break;
    }

    system.stiffness = restrict_matrix(
        operators_.stiffness, system.local_vertices, local_index_of);
    system.mass = restrict_matrix(
        operators_.mass, system.local_vertices, local_index_of);
    system.robin = restrict_matrix(
        operators_.boundary_mass, system.local_vertices, local_index_of);
    system.helmholtz = system.stiffness.cast<Complex>();
    system.helmholtz -= system.wavenumber * system.wavenumber
        * system.mass.cast<Complex>();
    system.helmholtz -= Complex(0.0, system.wavenumber)
        * system.robin.cast<Complex>();
    system.helmholtz.makeCompressed();

    const int local_size = static_cast<int>(system.local_vertices.size());
    system.rhs = ComplexMatrix::Zero(local_size, 3);
    for (int element : children_[target]) {
        const Triangle &triangle = fine_.elems[element];
        const Eigen::Matrix3cd &block = operators_.element_blocks[element];
        for (int i = 0; i < 3; ++i) {
            const int row = local_index_of(triangle[i]);
            if (row < 0) continue;
            for (int coarse_local = 0; coarse_local < 3; ++coarse_local) {
                Complex value = 0.0;
                for (int trial_local = 0; trial_local < 3; ++trial_local) {
                    value += block(i, trial_local)
                        * fine_dg_prolongation_.coeff(
                            3 * element + trial_local,
                            3 * target + coarse_local);
                }
                system.rhs(row, coarse_local) += value;
            }
        }
    }

    std::vector<int> active_constraint_rows;
    for (int local_col = 0; local_col < local_size; ++local_col) {
        const int global_col = system.local_vertices[local_col];
        for (Eigen::SparseMatrix<double>::InnerIterator it(
                 quasi_interpolation_, global_col); it; ++it) {
            const int coarse_row = it.row();
            if (patch_workspace.coarse_row_seen[coarse_row] == stamp) continue;
            patch_workspace.coarse_row_seen[coarse_row] = stamp;
            patch_workspace.coarse_row_index[coarse_row]
                = static_cast<int>(active_constraint_rows.size());
            active_constraint_rows.push_back(coarse_row);
        }
    }

    Eigen::MatrixXd active_constraints = Eigen::MatrixXd::Zero(
        active_constraint_rows.size(), local_size);
    for (int local_col = 0; local_col < local_size; ++local_col) {
        const int global_col = system.local_vertices[local_col];
        for (Eigen::SparseMatrix<double>::InnerIterator it(
                 quasi_interpolation_, global_col); it; ++it) {
            const int row = patch_workspace.coarse_row_index[it.row()];
            active_constraints(row, local_col) = it.value();
        }
    }

    if (active_constraints.rows() > 0) {
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(
            active_constraints.transpose());
        qr.setThreshold(1e-11);
        const int rank = qr.rank();
        system.constraints.resize(rank, local_size);
        const auto permutation = qr.colsPermutation().indices();
        for (int row = 0; row < rank; ++row)
            system.constraints.row(row) = active_constraints.row(permutation(row));
    } else {
        system.constraints.resize(0, local_size);
    }
    return system;
}

} // namespace lod2d::helmholtz

#include "helmholtz/patch_hierarchy.h"

#include <stdexcept>
#include <utility>

namespace lod2d::helmholtz {
namespace {

struct LocalPatchLevel {
    std::vector<int> elements;
    std::vector<int> vertices;
    std::vector<int> global_to_local;
};

std::vector<int> mesh_node_incidence(const TriMesh &mesh) {
    std::vector<int> incidence(mesh.nodes.size(), 0);
    for (const Triangle &triangle : mesh.elems)
        for (int vertex : triangle) ++incidence[vertex];
    return incidence;
}

LocalPatchLevel build_level(
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
        for (int vertex : touched)
            if (counts[vertex] == global_incidence[vertex])
                level.vertices.push_back(vertex);
    }
    for (int local = 0; local < static_cast<int>(level.vertices.size()); ++local)
        level.global_to_local[level.vertices[local]] = local;
    return level;
}

std::vector<int> prolong_elements(
    const Eigen::SparseMatrix<double> &prolongation,
    const std::vector<int> &coarse_elements) {
    std::vector<int> fine_elements;
    for (int coarse : coarse_elements) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(
                 prolongation, coarse);
             it; ++it) {
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
    for (int local_column = 0;
         local_column < static_cast<int>(coarse.vertices.size());
         ++local_column) {
        const int global_column = coarse.vertices[local_column];
        for (Eigen::SparseMatrix<double>::InnerIterator it(
                 global, global_column);
             it; ++it) {
            const int local_row = fine.global_to_local[it.row()];
            if (local_row >= 0 && it.value() != 0.0)
                triplets.emplace_back(
                    local_row, local_column, it.value());
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

HelmholtzDirichletPatchHierarchyBuilder::
HelmholtzDirichletPatchHierarchyBuilder(
    const std::vector<TriMesh> &meshes,
    const std::vector<Eigen::SparseMatrix<double>> &node_prolongations,
    const std::vector<Eigen::SparseMatrix<double>> &element_prolongations)
    : meshes_(meshes),
      node_prolongations_(node_prolongations),
      element_prolongations_(element_prolongations) {
    if (meshes_.empty()) {
        if (!node_prolongations_.empty() || !element_prolongations_.empty())
            throw std::invalid_argument(
                "patch hierarchy prolongations require meshes");
        return;
    }
    if (meshes_.size() != node_prolongations_.size() + 1
        || meshes_.size() != element_prolongations_.size() + 1)
        throw std::invalid_argument("inconsistent patch hierarchy sizes");
    node_incidence_.reserve(meshes_.size());
    for (const TriMesh &mesh : meshes_)
        node_incidence_.push_back(mesh_node_incidence(mesh));
}

std::vector<Eigen::SparseMatrix<double>>
HelmholtzDirichletPatchHierarchyBuilder::build(
    const std::vector<int> &coarse_patch_elements,
    const std::vector<int> &finest_patch_dofs) const {
    if (!available())
        throw std::runtime_error(
            "Dirichlet patch hierarchy requires nested NVB levels");
    if (coarse_patch_elements.empty() || finest_patch_dofs.empty())
        throw std::invalid_argument(
            "Dirichlet patch hierarchy inputs must be nonempty");

    LocalPatchLevel previous = build_level(
        meshes_.front(),
        coarse_patch_elements,
        node_incidence_.front());
    std::vector<Eigen::SparseMatrix<double>> result;
    result.reserve(node_prolongations_.size());
    for (int level = 0;
         level < static_cast<int>(node_prolongations_.size());
         ++level) {
        std::vector<int> next_elements = prolong_elements(
            element_prolongations_[level], previous.elements);
        const bool finest =
            level + 1 == static_cast<int>(node_prolongations_.size());
        LocalPatchLevel next = build_level(
            meshes_[level + 1],
            std::move(next_elements),
            node_incidence_[level + 1],
            finest ? &finest_patch_dofs : nullptr);
        if (!previous.vertices.empty()) {
            if (next.vertices.empty())
                throw std::runtime_error(
                    "nonempty coarse patch space has an empty fine space");
            result.push_back(restrict_prolongation(
                node_prolongations_[level], previous, next));
        }
        previous = std::move(next);
    }
    if (previous.vertices != finest_patch_dofs)
        throw std::runtime_error(
            "patch hierarchy finest DOFs do not match local operator");
    if (result.empty())
        throw std::runtime_error(
            "patch hierarchy has no nontrivial prolongation");
    return result;
}

} // namespace lod2d::helmholtz

#include "helmholtz/schwarz_patch.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace lod2d::helmholtz {
namespace {

std::uint64_t edge_key(int a, int b) {
    const auto lo = static_cast<std::uint32_t>(std::min(a, b));
    const auto hi = static_cast<std::uint32_t>(std::max(a, b));
    return (static_cast<std::uint64_t>(lo) << 32U) | hi;
}

std::array<int, 2> edge_vertices(std::uint64_t key) {
    return {
        static_cast<int>(key >> 32U),
        static_cast<int>(key & 0xffffffffU)};
}

const std::vector<TriMesh> &empty_meshes() {
    static const std::vector<TriMesh> value;
    return value;
}

const std::vector<Eigen::SparseMatrix<double>> &empty_prolongations() {
    static const std::vector<Eigen::SparseMatrix<double>> value;
    return value;
}

} // namespace

HelmholtzSchwarzPatchAssembler::HelmholtzSchwarzPatchAssembler(
    const TriMesh &fine_mesh,
    const Eigen::SparseMatrix<double> &fine_element_prolongation,
    const Eigen::SparseMatrix<double> &element_patches,
    int coarse_element_count,
    const HelmholtzOperators &operators)
    : HelmholtzSchwarzPatchAssembler(
          fine_mesh,
          fine_element_prolongation,
          element_patches,
          coarse_element_count,
          empty_meshes(),
          empty_prolongations(),
          empty_prolongations(),
          operators) {}

HelmholtzSchwarzPatchAssembler::HelmholtzSchwarzPatchAssembler(
    const TriMesh &fine_mesh,
    const Eigen::SparseMatrix<double> &fine_element_prolongation,
    const Eigen::SparseMatrix<double> &element_patches,
    int coarse_element_count,
    const std::vector<TriMesh> &hierarchy_meshes,
    const std::vector<Eigen::SparseMatrix<double>> &node_prolongations,
    const std::vector<Eigen::SparseMatrix<double>> &element_prolongations,
    const HelmholtzOperators &operators)
    : fine_mesh_(fine_mesh),
      element_patches_(element_patches),
      operators_(operators),
      children_(coarse_element_count),
      fine_incidence_(fine_mesh.nodes.size(), 0),
      hierarchy_builder_(
          std::make_unique<HelmholtzDirichletPatchHierarchyBuilder>(
              hierarchy_meshes, node_prolongations, element_prolongations)) {
    if (coarse_element_count <= 0)
        throw std::invalid_argument("Schwarz coarse element count must be positive");
    if (fine_element_prolongation.rows()
            != static_cast<int>(fine_mesh.elems.size())
        || fine_element_prolongation.cols() != coarse_element_count)
        throw std::invalid_argument(
            "Schwarz fine element prolongation has the wrong size");
    if (element_patches.rows() != coarse_element_count
        || element_patches.cols() != coarse_element_count)
        throw std::invalid_argument("Schwarz patch matrix has the wrong size");
    if (operators.system.rows() != static_cast<int>(fine_mesh.nodes.size())
        || operators.element_blocks.size() != fine_mesh.elems.size()
        || operators.refractive_index.size() != fine_mesh.elems.size())
        throw std::invalid_argument("Schwarz operator/mesh size mismatch");

    for (int coarse = 0; coarse < coarse_element_count; ++coarse) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(
                 fine_element_prolongation, coarse);
             it; ++it) {
            if (it.value() != 0.0) children_[coarse].push_back(it.row());
        }
        if (children_[coarse].empty())
            throw std::runtime_error("Schwarz coarse element has no fine children");
    }

    global_edge_counts_.reserve(3 * fine_mesh.elems.size());
    for (const Triangle &triangle : fine_mesh.elems) {
        for (int vertex : triangle) ++fine_incidence_[vertex];
        ++global_edge_counts_[edge_key(triangle[0], triangle[1])];
        ++global_edge_counts_[edge_key(triangle[1], triangle[2])];
        ++global_edge_counts_[edge_key(triangle[2], triangle[0])];
    }
}

HelmholtzSchwarzLocalSystem HelmholtzSchwarzPatchAssembler::assemble(
    int target,
    HelmholtzSchwarzArtificialBoundary boundary,
    double artificial_impedance_beta,
    bool assemble_mass,
    bool assemble_hierarchy) const {
    if (target < 0 || target >= patch_count())
        throw std::out_of_range("Schwarz patch target is out of range");
    if (boundary != HelmholtzSchwarzArtificialBoundary::HomogeneousDirichlet
        && boundary != HelmholtzSchwarzArtificialBoundary::Impedance)
        throw std::invalid_argument(
            "unknown Schwarz artificial boundary condition");
    if (!(artificial_impedance_beta >= 0.0)
        || !std::isfinite(artificial_impedance_beta))
        throw std::invalid_argument(
            "Schwarz artificial impedance beta must be finite and nonnegative");

    std::vector<int> patch_elements;
    std::vector<int> coarse_patch_elements;
    for (Eigen::SparseMatrix<double>::InnerIterator it(
             element_patches_, target);
         it; ++it) {
        if (it.value() == 0.0) continue;
        coarse_patch_elements.push_back(it.row());
        const std::vector<int> &children = children_[it.row()];
        patch_elements.insert(
            patch_elements.end(), children.begin(), children.end());
    }
    if (patch_elements.empty())
        throw std::runtime_error("Schwarz patch has no fine elements");

    std::vector<char> core_vertex(fine_mesh_.nodes.size(), 0);
    for (int element : children_[target]) {
        for (int vertex : fine_mesh_.elems[element]) core_vertex[vertex] = 1;
    }

    std::vector<int> patch_incidence(fine_mesh_.nodes.size(), 0);
    std::vector<int> touched_vertices;
    touched_vertices.reserve(3 * patch_elements.size());
    std::unordered_map<std::uint64_t, int> patch_edge_counts;
    patch_edge_counts.reserve(3 * patch_elements.size());
    for (int element : patch_elements) {
        const Triangle &triangle = fine_mesh_.elems[element];
        for (int vertex : triangle) {
            if (patch_incidence[vertex] == 0) touched_vertices.push_back(vertex);
            ++patch_incidence[vertex];
        }
        ++patch_edge_counts[edge_key(triangle[0], triangle[1])];
        ++patch_edge_counts[edge_key(triangle[1], triangle[2])];
        ++patch_edge_counts[edge_key(triangle[2], triangle[0])];
    }

    HelmholtzSchwarzLocalSystem system;
    system.global_dofs.reserve(touched_vertices.size());
    for (int vertex : touched_vertices) {
        if (boundary == HelmholtzSchwarzArtificialBoundary::Impedance
            || patch_incidence[vertex] == fine_incidence_[vertex])
            system.global_dofs.push_back(vertex);
    }
    if (system.global_dofs.empty())
        throw std::runtime_error("Schwarz patch has no local DOFs");

    std::vector<int> local_index(fine_mesh_.nodes.size(), -1);
    for (int local = 0;
         local < static_cast<int>(system.global_dofs.size()); ++local)
        local_index[system.global_dofs[local]] = local;
    for (int global : system.global_dofs) {
        if (core_vertex[global]) system.core_global_dofs.push_back(global);
    }

    std::vector<ComplexTriplet> triplets;
    std::vector<Eigen::Triplet<double>> mass_triplets;
    triplets.reserve(9 * patch_elements.size() + 4 * patch_edge_counts.size());
    if (assemble_mass) mass_triplets.reserve(9 * patch_elements.size());
    for (int element : patch_elements) {
        const Triangle &triangle = fine_mesh_.elems[element];
        const Eigen::Matrix3cd &block = operators_.element_blocks[element];
        Eigen::Matrix3d mass_block = Eigen::Matrix3d::Zero();
        if (assemble_mass) {
            static constexpr double mass_pattern[3][3] = {
                {2.0, 1.0, 1.0},
                {1.0, 2.0, 1.0},
                {1.0, 1.0, 2.0}};
            const Point2 first = fine_mesh_.nodes[triangle[0]];
            const Point2 second = fine_mesh_.nodes[triangle[1]];
            const Point2 third = fine_mesh_.nodes[triangle[2]];
            const double area = 0.5 * std::abs(
                (second.x() - first.x()) * (third.y() - first.y())
                - (second.y() - first.y()) * (third.x() - first.x()));
            for (int row = 0; row < 3; ++row)
                for (int column = 0; column < 3; ++column)
                    mass_block(row, column) =
                        operators_.refractive_index[element] * area
                        * mass_pattern[row][column] / 12.0;
        }
        for (int row = 0; row < 3; ++row) {
            const int local_row = local_index[triangle[row]];
            if (local_row < 0) continue;
            for (int column = 0; column < 3; ++column) {
                const int local_column = local_index[triangle[column]];
                if (local_column >= 0) {
                    triplets.emplace_back(
                        local_row, local_column, block(row, column));
                    if (assemble_mass)
                        mass_triplets.emplace_back(
                            local_row, local_column, mass_block(row, column));
                }
            }
        }
    }

    for (const auto &[key, count] : patch_edge_counts) {
        if (count != 1) continue;
        const std::array<int, 2> edge = edge_vertices(key);
        if (global_edge_counts_.at(key) == 1) {
            ++system.physical_boundary_edges;
            continue;
        }
        ++system.artificial_boundary_edges;
        if (boundary != HelmholtzSchwarzArtificialBoundary::Impedance)
            continue;

        const int first = local_index[edge[0]];
        const int second = local_index[edge[1]];
        if (first < 0 || second < 0)
            throw std::runtime_error(
                "impedance Schwarz boundary edge is missing a local DOF");
        const double length =
            (fine_mesh_.nodes[edge[0]] - fine_mesh_.nodes[edge[1]]).norm();
        const Complex diagonal = -Complex(0.0, operators_.wavenumber)
            * artificial_impedance_beta * length / 3.0;
        const Complex off_diagonal = 0.5 * diagonal;
        triplets.emplace_back(first, first, diagonal);
        triplets.emplace_back(second, second, diagonal);
        triplets.emplace_back(first, second, off_diagonal);
        triplets.emplace_back(second, first, off_diagonal);
    }

    const int local_size = static_cast<int>(system.global_dofs.size());
    system.matrix.resize(local_size, local_size);
    system.matrix.setFromTriplets(triplets.begin(), triplets.end());
    system.matrix.makeCompressed();
    if (assemble_mass) {
        system.mass.resize(local_size, local_size);
        system.mass.setFromTriplets(
            mass_triplets.begin(), mass_triplets.end());
        system.mass.makeCompressed();
    }
    if (assemble_hierarchy) {
        if (boundary != HelmholtzSchwarzArtificialBoundary::HomogeneousDirichlet)
            throw std::invalid_argument(
                "geometric Schwarz V-cycle currently requires Dirichlet "
                "artificial patch boundaries");
        system.geometric_prolongations = hierarchy_builder_->build(
            coarse_patch_elements, system.global_dofs);
    }
    return system;
}

} // namespace lod2d::helmholtz

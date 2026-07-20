#include "helmholtz/coarse_schwarz.h"

#include <algorithm>
#include <limits>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lod2d::helmholtz {
namespace {

std::vector<int> patch_coarse_vertices(
    const TriMesh &coarse_mesh,
    const Eigen::SparseMatrix<double> &element_patches,
    int target) {
    std::vector<char> seen(coarse_mesh.nodes.size(), 0);
    std::vector<int> vertices;
    for (Eigen::SparseMatrix<double>::InnerIterator it(element_patches, target);
         it; ++it) {
        if (it.value() == 0.0) continue;
        for (int vertex : coarse_mesh.elems[it.row()]) {
            if (seen[vertex]) continue;
            seen[vertex] = 1;
            vertices.push_back(vertex);
        }
    }
    std::sort(vertices.begin(), vertices.end());
    return vertices;
}

ComplexSparseMatrix principal_submatrix(
    const ComplexSparseMatrix &matrix,
    const std::vector<int> &indices) {
    std::vector<int> local_index(matrix.rows(), -1);
    for (int local = 0; local < static_cast<int>(indices.size()); ++local)
        local_index[indices[local]] = local;

    std::vector<ComplexTriplet> triplets;
    for (int local_col = 0; local_col < static_cast<int>(indices.size());
         ++local_col) {
        const int global_col = indices[local_col];
        for (ComplexSparseMatrix::InnerIterator it(matrix, global_col); it; ++it) {
            const int local_row = local_index[it.row()];
            if (local_row >= 0)
                triplets.emplace_back(local_row, local_col, it.value());
        }
    }

    ComplexSparseMatrix local(indices.size(), indices.size());
    local.setFromTriplets(triplets.begin(), triplets.end());
    local.makeCompressed();
    return local;
}

} // namespace

HelmholtzCoarseRasPreconditioner::HelmholtzCoarseRasPreconditioner(
    const ComplexSparseMatrix &coarse_operator,
    const TriMesh &coarse_mesh,
    const Eigen::SparseMatrix<double> &element_patches)
    : dimension_(coarse_operator.rows()) {
    if (dimension_ <= 0 || coarse_operator.cols() != dimension_)
        throw std::invalid_argument("coarse LOD Schwarz operator must be square");
    if (dimension_ != static_cast<int>(coarse_mesh.nodes.size()))
        throw std::invalid_argument("coarse LOD Schwarz operator/mesh size mismatch");
    const int element_count = static_cast<int>(coarse_mesh.elems.size());
    if (element_patches.rows() != element_count
        || element_patches.cols() != element_count)
        throw std::invalid_argument("coarse LOD Schwarz patch matrix size mismatch");

    Eigen::VectorXi multiplicity = Eigen::VectorXi::Zero(dimension_);
    diagnostics_.min_local_dofs = std::numeric_limits<int>::max();
    subdomains_.reserve(element_count);
    for (int target = 0; target < element_count; ++target) {
        Subdomain subdomain;
        subdomain.global_dofs =
            patch_coarse_vertices(coarse_mesh, element_patches, target);
        if (subdomain.global_dofs.empty())
            throw std::runtime_error("coarse LOD Schwarz patch has no DOFs");
        for (int dof : subdomain.global_dofs) ++multiplicity(dof);

        subdomain.solver =
            std::make_unique<Eigen::SparseLU<ComplexSparseMatrix>>();

        const int local_dofs = static_cast<int>(subdomain.global_dofs.size());
        diagnostics_.min_local_dofs =
            std::min(diagnostics_.min_local_dofs, local_dofs);
        diagnostics_.max_local_dofs =
            std::max(diagnostics_.max_local_dofs, local_dofs);
        subdomains_.push_back(std::move(subdomain));
    }

    std::vector<char> factorization_failed(element_count, 0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int target = 0; target < element_count; ++target) {
        Subdomain &subdomain = subdomains_[target];
        const ComplexSparseMatrix local =
            principal_submatrix(coarse_operator, subdomain.global_dofs);
        subdomain.solver->analyzePattern(local);
        subdomain.solver->factorize(local);
        if (subdomain.solver->info() != Eigen::Success)
            factorization_failed[target] = 1;
    }
    for (int target = 0; target < element_count; ++target) {
        if (factorization_failed[target])
            throw std::runtime_error(
                "coarse LOD Schwarz local factorization failed at patch "
                + std::to_string(target));
    }


    injection_weights_.resize(dimension_);
    Eigen::VectorXd partition_sum = Eigen::VectorXd::Zero(dimension_);
    for (int dof = 0; dof < dimension_; ++dof) {
        if (multiplicity(dof) <= 0)
            throw std::runtime_error("coarse LOD Schwarz patches do not cover every DOF");
        injection_weights_(dof) = 1.0 / static_cast<double>(multiplicity(dof));
    }
    for (const Subdomain &subdomain : subdomains_)
        for (int dof : subdomain.global_dofs)
            partition_sum(dof) += injection_weights_(dof);

    diagnostics_.subdomains = static_cast<int>(subdomains_.size());
    diagnostics_.partition_unity_error =
        (partition_sum - Eigen::VectorXd::Ones(dimension_)).lpNorm<Eigen::Infinity>();
}

ComplexVector HelmholtzCoarseRasPreconditioner::apply(
    const ComplexVector &right_hand_side) const {
    if (right_hand_side.size() != dimension_)
        throw std::invalid_argument("coarse LOD Schwarz right-hand side size mismatch");
    if (!right_hand_side.allFinite())
        throw std::invalid_argument("coarse LOD Schwarz right-hand side is non-finite");

    int worker_count = 1;
#ifdef _OPENMP
    worker_count = omp_get_max_threads();
#endif
    std::vector<ComplexVector> partial_results(
        worker_count, ComplexVector::Zero(dimension_));
    std::vector<char> solve_failed(subdomains_.size(), 0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int subdomain_index = 0;
         subdomain_index < static_cast<int>(subdomains_.size());
         ++subdomain_index) {
        const Subdomain &subdomain = subdomains_[subdomain_index];
        int worker = 0;
#ifdef _OPENMP
        worker = omp_get_thread_num();
#endif
        ComplexVector local_rhs(subdomain.global_dofs.size());
        for (int local = 0; local < local_rhs.size(); ++local)
            local_rhs(local) = right_hand_side(subdomain.global_dofs[local]);
        const ComplexVector local_solution = subdomain.solver->solve(local_rhs);
        if (subdomain.solver->info() != Eigen::Success
            || !local_solution.allFinite()) {
            solve_failed[subdomain_index] = 1;
            continue;
        }
        for (int local = 0; local < local_solution.size(); ++local) {
            const int global = subdomain.global_dofs[local];
            partial_results[worker](global) +=
                injection_weights_(global) * local_solution(local);
        }
    }
    if (std::find(solve_failed.begin(), solve_failed.end(), 1)
        != solve_failed.end())
        throw std::runtime_error("coarse LOD Schwarz local solve failed");

    ComplexVector result = ComplexVector::Zero(dimension_);
    for (const ComplexVector &partial : partial_results) result += partial;
    return result;
}

} // namespace lod2d::helmholtz

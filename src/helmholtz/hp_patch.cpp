#include "helmholtz/hp_patch.h"

#include <Eigen/QR>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace lod2d::helmholtz {
namespace {

std::vector<std::vector<int>> element_children(
    const Eigen::SparseMatrix<double> &prolongation,
    int coarse_count) {
    if (prolongation.cols() != coarse_count)
        throw std::invalid_argument(
            "hp patch element prolongation has wrong column count");
    std::vector<std::vector<int>> result(coarse_count);
    for (int coarse = 0; coarse < prolongation.outerSize(); ++coarse) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(
                 prolongation, coarse); it; ++it)
            if (it.value() != 0.0) result[coarse].push_back(it.row());
    }
    return result;
}

template <typename Scalar>
Eigen::SparseMatrix<Scalar> restrict_sparse(
    const Eigen::SparseMatrix<Scalar> &global,
    const std::vector<int> &local_dofs,
    const std::vector<int> &global_to_local) {
    std::vector<Eigen::Triplet<Scalar>> triplets;
    for (int local_col = 0;
         local_col < static_cast<int>(local_dofs.size()); ++local_col) {
        const int global_col = local_dofs[local_col];
        for (typename Eigen::SparseMatrix<Scalar>::InnerIterator it(
                 global, global_col); it; ++it) {
            const int local_row = global_to_local[it.row()];
            if (local_row >= 0 && it.value() != Scalar(0))
                triplets.emplace_back(local_row, local_col, it.value());
        }
    }
    Eigen::SparseMatrix<Scalar> result(
        static_cast<int>(local_dofs.size()),
        static_cast<int>(local_dofs.size()));
    result.setFromTriplets(triplets.begin(), triplets.end());
    result.makeCompressed();
    return result;
}

} // namespace

HelmholtzHpPatchAssembler::HelmholtzHpPatchAssembler(
    const TriMesh &coarse,
    const HpTriSpace &fine_space,
    const Eigen::SparseMatrix<double> &fine_element_prolongation,
    const Eigen::SparseMatrix<double> &patches,
    const HelmholtzHpInterpolation &interpolation,
    const HelmholtzHpOperators &operators)
    : coarse_(coarse),
      fine_space_(fine_space),
      patches_(patches),
      interpolation_(interpolation),
      operators_(operators),
      children_(element_children(
          fine_element_prolongation,
          static_cast<int>(coarse.elems.size()))) {
    const int coarse_elements = static_cast<int>(coarse.elems.size());
    if (patches.rows() != coarse_elements || patches.cols() != coarse_elements)
        throw std::invalid_argument("hp patch matrix has wrong dimensions");
    if (interpolation.coarse_injection.rows() != fine_space.dof_count()
        || interpolation.coarse_injection.cols()
            != static_cast<int>(coarse.nodes.size()))
        throw std::invalid_argument("hp coarse injection has wrong dimensions");
    if (interpolation.quasi_interpolation.rows()
            != static_cast<int>(coarse.nodes.size())
        || interpolation.quasi_interpolation.cols() != fine_space.dof_count())
        throw std::invalid_argument(
            "hp quasi-interpolation has wrong dimensions");
    if (operators.system.rows() != fine_space.dof_count()
        || operators.element_blocks.size()
            != fine_space.mesh().elems.size())
        throw std::invalid_argument("hp Helmholtz operators have wrong size");
}

HelmholtzPatchSystem HelmholtzHpPatchAssembler::assemble(int target) const {
    if (target < 0 || target >= patch_count())
        throw std::out_of_range("hp Helmholtz patch target is out of range");

    HelmholtzPatchSystem system;
    system.target_element = target;
    system.wavenumber = operators_.wavenumber;
    for (Eigen::SparseMatrix<double>::InnerIterator it(patches_, target);
         it; ++it) {
        if (it.value() == 0.0) continue;
        const auto &fine_children = children_[it.row()];
        system.patch_elements.insert(
            system.patch_elements.end(),
            fine_children.begin(), fine_children.end());
    }
    if (system.patch_elements.empty())
        throw std::runtime_error("hp Helmholtz corrector patch is empty");

    const int global_dofs = fine_space_.dof_count();
    std::vector<int> patch_incidence(global_dofs, 0);
    std::vector<int> touched;
    touched.reserve(
        system.patch_elements.size() * fine_space_.local_dof_count());
    Point2 lower = fine_space_.dof_points()[
        fine_space_.element_dofs()[system.patch_elements.front()][0]];
    Point2 upper = lower;
    for (int element : system.patch_elements) {
        for (int dof : fine_space_.element_dofs()[element]) {
            lower = lower.cwiseMin(fine_space_.dof_points()[dof]);
            upper = upper.cwiseMax(fine_space_.dof_points()[dof]);
            if (patch_incidence[dof]++ == 0) touched.push_back(dof);
        }
    }
    system.diameter = (upper - lower).norm();
    std::sort(touched.begin(), touched.end());
    for (int dof : touched) {
        if (patch_incidence[dof] == fine_space_.dof_incidence()[dof])
            system.local_vertices.push_back(dof);
    }
    if (system.local_vertices.empty())
        throw std::runtime_error(
            "hp Helmholtz patch has no unconstrained fine DOFs");

    std::vector<int> global_to_local(global_dofs, -1);
    for (int local = 0;
         local < static_cast<int>(system.local_vertices.size()); ++local)
        global_to_local[system.local_vertices[local]] = local;
    system.stiffness = restrict_sparse(
        operators_.stiffness, system.local_vertices, global_to_local);
    system.mass = restrict_sparse(
        operators_.mass, system.local_vertices, global_to_local);
    system.robin = restrict_sparse(
        operators_.boundary_mass, system.local_vertices, global_to_local);
    system.helmholtz = restrict_sparse(
        operators_.system, system.local_vertices, global_to_local);

    const int local_size = static_cast<int>(system.local_vertices.size());
    system.rhs = ComplexMatrix::Zero(local_size, 3);
    const Triangle &coarse_target = coarse_.elems[target];
    for (int element : children_[target]) {
        const auto &dofs = fine_space_.element_dofs()[element];
        Eigen::MatrixXd target_basis(dofs.size(), 3);
        for (int hp_local = 0;
             hp_local < static_cast<int>(dofs.size()); ++hp_local) {
            for (int coarse_local = 0; coarse_local < 3; ++coarse_local)
                target_basis(hp_local, coarse_local) =
                    interpolation_.coarse_injection.coeff(
                        dofs[hp_local], coarse_target[coarse_local]);
        }
        const ComplexMatrix local_rhs =
            operators_.element_blocks[element] * target_basis.cast<Complex>();
        for (int hp_local = 0;
             hp_local < static_cast<int>(dofs.size()); ++hp_local) {
            const int local_row = global_to_local[dofs[hp_local]];
            if (local_row >= 0)
                system.rhs.row(local_row) += local_rhs.row(hp_local);
        }
    }

    std::vector<int> active_rows;
    std::vector<int> row_position(coarse_.nodes.size(), -1);
    for (int global_col : system.local_vertices) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(
                 interpolation_.quasi_interpolation, global_col); it; ++it) {
            if (row_position[it.row()] >= 0) continue;
            row_position[it.row()] = static_cast<int>(active_rows.size());
            active_rows.push_back(it.row());
        }
    }
    Eigen::MatrixXd active = Eigen::MatrixXd::Zero(
        active_rows.size(), local_size);
    for (int local_col = 0; local_col < local_size; ++local_col) {
        const int global_col = system.local_vertices[local_col];
        for (Eigen::SparseMatrix<double>::InnerIterator it(
                 interpolation_.quasi_interpolation, global_col); it; ++it)
            active(row_position[it.row()], local_col) = it.value();
    }
    if (active.rows() > 0) {
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(active.transpose());
        qr.setThreshold(1e-11);
        const int rank = qr.rank();
        const auto permutation = qr.colsPermutation().indices();
        system.constraints.resize(rank, local_size);
        for (int row = 0; row < rank; ++row)
            system.constraints.row(row) = active.row(permutation(row));
    } else {
        system.constraints.resize(0, local_size);
    }

    for (int element : system.patch_elements) {
        const Triangle &triangle = fine_space_.mesh().elems[element];
        static constexpr int local_edges[3][2] =
            {{0, 1}, {1, 2}, {2, 0}};
        for (const auto &local_edge : local_edges) {
            const int edge = fine_space_.edge_index(
                triangle[local_edge[0]], triangle[local_edge[1]]);
            if (fine_space_.boundary_edges()[edge]) {
                system.touches_physical_boundary = true;
                break;
            }
        }
        if (system.touches_physical_boundary) break;
    }
    return system;
}

HelmholtzHpPatchSolveResult
HelmholtzHpPatchAssembler::solve_direct_saddle(int target) const {
    HelmholtzHpPatchSolveResult result;
    result.system = assemble(target);
    HelmholtzPatchSolverConfig config;
    config.kind = HelmholtzPatchSolverKind::DirectSaddle;
    config.fallback_to_direct = false;
    result.primal = solve_helmholtz_patch(result.system, config);
    result.adjoint_corrector = result.primal.corrector.conjugate();
    result.adjoint_multipliers = result.primal.multipliers.conjugate();
    return result;
}

} // namespace lod2d::helmholtz

#include "helmholtz/model.h"
#include "helmholtz/patch_system.h"
#include "helmholtz/schwarz_patch.h"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace lod2d::helmholtz;

namespace {

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

HelmholtzSchwarzPatchAssembler make_schwarz_assembler(
    const HelmholtzLodModel &model) {
    const HelmholtzProblemData &problem = model.problem();
    return HelmholtzSchwarzPatchAssembler(
        problem.fine,
        problem.fine_element_prolongation,
        problem.patches,
        static_cast<int>(problem.coarse.elems.size()),
        model.operators());
}

HelmholtzPatchAssembler make_corrector_assembler(
    const HelmholtzLodModel &model) {
    const HelmholtzProblemData &problem = model.problem();
    return HelmholtzPatchAssembler(
        problem.coarse,
        problem.fine,
        problem.fine_element_prolongation,
        problem.fine_dg_prolongation,
        problem.quasi_interpolation,
        problem.patches,
        problem.fine_hierarchy_meshes,
        problem.fine_node_level_prolongations,
        problem.fine_element_level_prolongations,
        model.operators());
}

ComplexSparseMatrix restrict_matrix(
    const ComplexSparseMatrix &global_matrix,
    const std::vector<int> &global_dofs) {
    std::vector<int> local_index(global_matrix.rows(), -1);
    for (int local = 0; local < static_cast<int>(global_dofs.size()); ++local)
        local_index[global_dofs[local]] = local;

    std::vector<ComplexTriplet> triplets;
    triplets.reserve(global_matrix.nonZeros());
    for (int column = 0; column < global_matrix.outerSize(); ++column) {
        for (ComplexSparseMatrix::InnerIterator it(global_matrix, column);
             it; ++it) {
            const int local_row = local_index[it.row()];
            const int local_column = local_index[it.col()];
            if (local_row >= 0 && local_column >= 0)
                triplets.emplace_back(local_row, local_column, it.value());
        }
    }
    ComplexSparseMatrix result(global_dofs.size(), global_dofs.size());
    result.setFromTriplets(triplets.begin(), triplets.end());
    result.makeCompressed();
    return result;
}

} // namespace

int main() {
    try {
        HelmholtzProblemConfig config;
        config.H = 3;
        config.h = 6;
        config.ell = 1;
        config.wavenumber = 2.0;
        HelmholtzLodModel model = HelmholtzLodModel::build(config);

        HelmholtzSchwarzPatchAssembler schwarz =
            make_schwarz_assembler(model);
        HelmholtzPatchAssembler corrector =
            make_corrector_assembler(model);
        require(schwarz.patch_count() == corrector.patch_count(),
                "Schwarz and corrector patch counts disagree");

        int artificial_edges = 0;
        for (int target = 0; target < schwarz.patch_count(); ++target) {
            const HelmholtzPatchSystem old_system =
                corrector.assemble(target);
            const HelmholtzSchwarzLocalSystem dirichlet = schwarz.assemble(
                target,
                HelmholtzSchwarzArtificialBoundary::HomogeneousDirichlet,
                1.0, true);
            require(dirichlet.global_dofs == old_system.local_vertices,
                    "lightweight Dirichlet patch changed the local DOF order");
            const double relative_matrix_error =
                (dirichlet.matrix - old_system.helmholtz).norm()
                / std::max(1.0, old_system.helmholtz.norm());
            require(relative_matrix_error < 1e-13,
                    "lightweight Dirichlet patch changed the local matrix");
            const double relative_mass_error =
                (dirichlet.mass - old_system.mass).norm()
                / std::max(1.0, old_system.mass.norm());
            require(relative_mass_error < 1e-13,
                    "lightweight Dirichlet patch changed the local mass matrix");

            const HelmholtzSchwarzLocalSystem impedance = schwarz.assemble(
                target, HelmholtzSchwarzArtificialBoundary::Impedance);
            require(impedance.global_dofs.size() >= dirichlet.global_dofs.size(),
                    "impedance patch unexpectedly removed boundary DOFs");
            require(impedance.matrix.rows()
                        == static_cast<int>(impedance.global_dofs.size()),
                    "impedance patch matrix has the wrong size");
            artificial_edges += impedance.artificial_boundary_edges;
        }
        require(artificial_edges > 0,
                "proper test patches did not contain artificial boundaries");

        HelmholtzProblemConfig full_config;
        full_config.H = 1;
        full_config.h = 3;
        full_config.ell = 20;
        full_config.wavenumber = 2.0;
        HelmholtzLodModel full_model = HelmholtzLodModel::build(full_config);
        HelmholtzSchwarzPatchAssembler full =
            make_schwarz_assembler(full_model);
        const HelmholtzSchwarzLocalSystem full_patch = full.assemble(
            0, HelmholtzSchwarzArtificialBoundary::Impedance);
        require(full_patch.artificial_boundary_edges == 0,
                "whole-domain patch retained an artificial boundary");
        require(full_patch.physical_boundary_edges > 0,
                "whole-domain patch lost the physical Robin boundary");
        require(full_patch.global_dofs.size()
                    == full_model.problem().fine.nodes.size(),
                "whole-domain impedance patch lost fine-grid DOFs");

        const ComplexSparseMatrix expected = restrict_matrix(
            full_model.operators().system, full_patch.global_dofs);
        const double full_error = (full_patch.matrix - expected).norm()
            / std::max(1.0, expected.norm());
        require(full_error < 1e-13,
                "whole-domain impedance patch differs from the global operator");

        std::cout << "Dirichlet matrix golden patches="
                  << schwarz.patch_count()
                  << " artificial_edges=" << artificial_edges
                  << " whole_domain_error=" << full_error << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_schwarz_patch failed: "
                  << error.what() << '\n';
        return 1;
    }
}

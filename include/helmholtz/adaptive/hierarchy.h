#pragma once

#include "mesh/refine.h"

#include <cstdint>
#include <vector>

namespace lod2d::helmholtz::adaptive {

struct NestedFineMesh {
    RefineOutput refinement;
    std::vector<int> element_levels;
};

class AdaptiveMeshHierarchy {
public:
    AdaptiveMeshHierarchy(
        const TriMesh &initial_mesh,
        int initial_coarse_level,
        int fine_level);

    const TriMesh &initial_mesh() const { return initial_mesh_; }
    const TriMesh &coarse_mesh() const { return coarse_mesh_; }
    const std::vector<int> &coarse_levels() const { return coarse_levels_; }
    const std::vector<std::uint64_t> &coarse_element_ids() const { return coarse_element_ids_; }
    const std::vector<std::uint64_t> &coarse_parent_ids() const { return coarse_parent_ids_; }

    // The initial uniform corrector level. Locally refined fine elements may
    // have larger entries in fine_element_levels().
    int fine_level() const { return fine_level_; }
    const TriMesh &fine_mesh() const { return fine_completion_.refinement.mesh; }
    const TriMesh &cert_audit_mesh() const { return cert_audit_mesh_; }
    const std::vector<int> &fine_element_levels() const {
        return fine_completion_.element_levels;
    }
    const std::vector<int> &cert_audit_element_levels() const {
        return cert_audit_element_levels_;
    }

    const Eigen::SparseMatrix<double> &coarse_to_fine() const {
        return fine_completion_.refinement.P_node;
    }
    const Eigen::SparseMatrix<double> &fine_to_cert_audit() const {
        return fine_to_cert_audit_.P_node;
    }
    const Eigen::SparseMatrix<double> &coarse_to_cert_audit() const {
        return coarse_to_cert_audit_;
    }
    const Eigen::SparseMatrix<double> &coarse_elements_to_fine() const {
        return fine_completion_.refinement.P_elem;
    }
    const Eigen::SparseMatrix<double> &fine_elements_to_cert_audit() const {
        return fine_to_cert_audit_.P_elem;
    }
    const Eigen::SparseMatrix<double> &coarse_elements_to_cert_audit() const {
        return coarse_elements_to_cert_audit_;
    }
    const Eigen::SparseMatrix<double> &coarse_dg_to_fine() const {
        return fine_completion_.refinement.P_dg;
    }
    const Eigen::SparseMatrix<double> &fine_dg_to_cert_audit() const {
        return fine_to_cert_audit_.P_dg;
    }
    const Eigen::SparseMatrix<double> &coarse_dg_to_cert_audit() const {
        return coarse_dg_to_cert_audit_;
    }
    const Eigen::SparseMatrix<double> &fine_quasi_interpolation() const {
        return fine_quasi_interpolation_;
    }
    const Eigen::SparseMatrix<double> &cert_audit_quasi_interpolation() const {
        return cert_audit_quasi_interpolation_;
    }

    const std::vector<int> &fine_parent_coarse_elements() const {
        return fine_parent_coarse_elements_;
    }
    const std::vector<int> &cert_audit_parent_fine_elements() const {
        return cert_audit_parent_fine_elements_;
    }

    std::uint64_t coarse_mesh_version() const { return coarse_mesh_version_; }
    std::uint64_t fine_mesh_version() const { return fine_mesh_version_; }
    std::uint64_t cert_audit_mesh_version() const { return cert_audit_mesh_version_; }
    std::uint64_t interpolation_version() const { return interpolation_version_; }
    std::uint64_t boundary_version() const { return boundary_version_; }
    std::uint64_t corrector_space_version() const { return corrector_space_version_; }

    NestedFineMesh build_nested_fine_mesh() const;
    void refine(const std::vector<int> &marked_elements);
    void refine_fine_elements(const std::vector<int> &marked_fine_elements);
    std::vector<int> fine_elements_in_coarse_patch(
        const std::vector<int> &coarse_patch_elements) const;
    void refine_fine_in_coarse_patch(
        const std::vector<int> &coarse_patch_elements);
    void refine_cert_audit_from_fine_elements(
        const std::vector<int> &marked_fine_elements);

    Eigen::SparseMatrix<double> fine_kernel_constraints(
        const std::vector<int> &fine_dofs) const;
    Eigen::SparseMatrix<double> cert_audit_kernel_constraints(
        const std::vector<int> &cert_audit_dofs) const;

private:
    TriMesh initial_mesh_;
    TriMesh coarse_mesh_;
    std::vector<int> coarse_levels_;
    std::vector<std::uint64_t> coarse_element_ids_;
    std::vector<std::uint64_t> coarse_parent_ids_;
    int fine_level_ = 0;
    std::uint64_t next_element_id_ = 1;
    NestedFineMesh fine_completion_;
    TriMesh cert_audit_mesh_;
    RefineOutput fine_to_cert_audit_;
    std::vector<int> cert_audit_element_levels_;
    Eigen::SparseMatrix<double> coarse_to_cert_audit_;
    Eigen::SparseMatrix<double> coarse_elements_to_cert_audit_;
    Eigen::SparseMatrix<double> coarse_dg_to_cert_audit_;
    Eigen::SparseMatrix<double> fine_quasi_interpolation_;
    Eigen::SparseMatrix<double> cert_audit_quasi_interpolation_;
    std::vector<int> fine_parent_coarse_elements_;
    std::vector<int> cert_audit_parent_fine_elements_;
    std::uint64_t coarse_mesh_version_ = 0;
    std::uint64_t fine_mesh_version_ = 0;
    std::uint64_t cert_audit_mesh_version_ = 0;
    std::uint64_t interpolation_version_ = 0;
    std::uint64_t boundary_version_ = 0;
    std::uint64_t corrector_space_version_ = 0;

    void refresh_fine_embedding();
    void refresh_coarse_to_cert_audit();
    void validate_current_embeddings() const;
};

NestedFineMesh complete_to_fine_level(
    const TriMesh &coarse_mesh,
    const std::vector<int> &coarse_levels,
    int fine_level);

// Build the exact P1/element/DG embedding for two conforming nested meshes.
// The routine is intentionally geometry based and serves as the full-rebuild
// correctness path after independent H, h, and audit refinements.
RefineOutput build_nested_mesh_embedding(
    const TriMesh &parent_mesh,
    const TriMesh &child_mesh);

std::vector<int> refinement_child_levels(
    const TriMesh &parent_mesh,
    const std::vector<int> &parent_levels,
    const RefineOutput &refinement);

std::vector<int> fine_element_parents(
    const Eigen::SparseMatrix<double> &fine_element_prolongation,
    int fine_element_count,
    int coarse_element_count);

Eigen::SparseMatrix<double> restrict_constraint_columns(
    const Eigen::SparseMatrix<double> &constraints,
    const std::vector<int> &ambient_dofs);

} // namespace lod2d::helmholtz::adaptive

#pragma once

#include "mesh/refine.h"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <optional>
#include <string>
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

enum class ReferenceEpochRefinementStatus {
    NoChange,
    Refined,
    ReferenceRefreshRequired,
};

struct ReferenceEpochRefinementResult {
    ReferenceEpochRefinementStatus status =
        ReferenceEpochRefinementStatus::NoChange;
    std::size_t previous_element_count = 0;
    std::size_t current_element_count = 0;
    std::string detail;
    // Incremental-refinement phase timings. They remain zero for operations
    // that do not use the corresponding path.
    double time_nvb_refine = 0.0;
    double time_embedding_composition = 0.0;
    double time_parent_map_update = 0.0;
    double time_proposed_embedding_update = 0.0;
    double time_candidate_quasi_interpolation = 0.0;
    double time_embedding_validation = 0.0;

    bool changed() const {
        return status == ReferenceEpochRefinementStatus::Refined;
    }
};

// Immutable result of previewing one marked NVB refinement of the committed
// coarse mesh.  The expensive NVB closure and reference embedding can be
// transferred into a pending proposal when the preview is selected.  Mesh
// version stamps prevent a preview from being adopted after its source
// hierarchy has changed.
struct ReferenceEpochCoarseRefinementPreview {
    struct CachedEmbedding {
        Eigen::SparseMatrix<double> P_node;
        Eigen::SparseMatrix<double> P_elem;
        Eigen::SparseMatrix<double> P_dg;
    };

    std::vector<int> marked_elements;
    RefineOutput refinement;
    std::vector<int> element_levels;
    std::optional<CachedEmbedding> to_reference;
    std::uint64_t coarse_mesh_version = 0;
    std::uint64_t reference_mesh_version = 0;
    double time_nvb_refine = 0.0;
    double time_reference_embedding_update = 0.0;

    bool reference_contained() const { return to_reference.has_value(); }
};

enum class ReferenceEpochRefinementGuardPoint {
    BeforeNvb,
    AfterNvb,
};

using ReferenceEpochRefinementGuard = std::function<void(
    ReferenceEpochRefinementGuardPoint, const TriMesh &)>;

struct ReferenceEpochCandidateDeepeningProbe {
    bool target_satisfied = false;
    std::size_t refinement_steps = 0;
    std::size_t final_element_count = 0;
    std::size_t final_node_count = 0;
    std::size_t protected_parent_spill = 0;
    int minimum_gap_margin = 0;
    int minimum_full_target_gap = 0;
};

struct AmbientRatioEnforcementResult {
    bool changed = false;
    std::size_t refinement_steps = 0;
    std::size_t refined_elements = 0;
    double maximum_ratio = 0.0;
};

// Practical-adaptive hierarchy for one fixed reference epoch:
//
//     V_H subset V_h^ref subset V_amb.
//
// Unlike AdaptiveMeshHierarchy, coarse refinement never expands the
// reference mesh.  A candidate H refinement is committed only after its
// embedding into the immutable reference mesh has been built successfully.
// The ambient shadow is advanced explicitly by enforce_ambient_ratio().
class ReferenceEpochHierarchy {
public:
    ReferenceEpochHierarchy(
        const TriMesh &initial_mesh,
        int initial_coarse_level,
        int reference_level,
        std::uint64_t initial_reference_epoch = 0);

    const TriMesh &initial_mesh() const { return initial_mesh_; }
    const TriMesh &coarse_mesh() const { return coarse_mesh_; }
    const TriMesh &reference_mesh() const {
        return reference_completion_.refinement.mesh;
    }
    const TriMesh &ambient_mesh() const {
        return reference_to_ambient_.mesh;
    }

    // Revised reference-epoch terminology.  The legacy ambient accessors
    // above remain available for old runners, but new code must treat this
    // mesh only as the candidate for the next reference epoch.
    const TriMesh &candidate_mesh() const { return ambient_mesh(); }
    const std::vector<int> &candidate_element_levels() const {
        return ambient_element_levels_;
    }
    const Eigen::SparseMatrix<double> &reference_to_candidate() const {
        return reference_to_ambient();
    }
    const Eigen::SparseMatrix<double> &coarse_to_candidate() const {
        return coarse_to_ambient();
    }
    const Eigen::SparseMatrix<double> &reference_elements_to_candidate() const {
        return reference_elements_to_ambient();
    }
    const Eigen::SparseMatrix<double> &coarse_elements_to_candidate() const {
        return coarse_elements_to_ambient();
    }
    const Eigen::SparseMatrix<double> &reference_dg_to_candidate() const {
        return reference_dg_to_ambient();
    }
    const Eigen::SparseMatrix<double> &coarse_dg_to_candidate() const {
        return coarse_dg_to_ambient();
    }
    const Eigen::SparseMatrix<double> &candidate_quasi_interpolation() const {
        return ambient_quasi_interpolation();
    }
    const std::vector<int> &candidate_parent_coarse_elements() const {
        return ambient_parent_coarse_elements();
    }
    const std::vector<int> &candidate_parent_reference_elements() const {
        return ambient_parent_reference_elements();
    }

    const std::vector<int> &coarse_levels() const { return coarse_levels_; }
    const std::vector<int> &reference_element_levels() const {
        return reference_completion_.element_levels;
    }
    const std::vector<int> &ambient_element_levels() const {
        return ambient_element_levels_;
    }
    int reference_level() const { return reference_level_; }

    const Eigen::SparseMatrix<double> &coarse_to_reference() const {
        return reference_completion_.refinement.P_node;
    }
    const Eigen::SparseMatrix<double> &reference_to_ambient() const {
        return reference_to_ambient_.P_node;
    }
    const Eigen::SparseMatrix<double> &coarse_to_ambient() const {
        return coarse_to_ambient_.P_node;
    }
    const Eigen::SparseMatrix<double> &coarse_elements_to_reference() const {
        return reference_completion_.refinement.P_elem;
    }
    const Eigen::SparseMatrix<double> &reference_elements_to_ambient() const {
        return reference_to_ambient_.P_elem;
    }
    const Eigen::SparseMatrix<double> &coarse_elements_to_ambient() const {
        return coarse_to_ambient_.P_elem;
    }
    const Eigen::SparseMatrix<double> &coarse_dg_to_reference() const {
        return reference_completion_.refinement.P_dg;
    }
    const Eigen::SparseMatrix<double> &reference_dg_to_ambient() const {
        return reference_to_ambient_.P_dg;
    }
    const Eigen::SparseMatrix<double> &coarse_dg_to_ambient() const {
        return coarse_to_ambient_.P_dg;
    }
    const Eigen::SparseMatrix<double> &reference_quasi_interpolation() const {
        return reference_quasi_interpolation_;
    }
    const Eigen::SparseMatrix<double> &ambient_quasi_interpolation() const {
        return ambient_quasi_interpolation_;
    }

    const std::vector<int> &reference_parent_coarse_elements() const {
        return reference_parent_coarse_elements_;
    }
    const std::vector<int> &ambient_parent_coarse_elements() const {
        return ambient_parent_coarse_elements_;
    }
    const std::vector<int> &ambient_parent_reference_elements() const {
        return ambient_parent_reference_elements_;
    }

    std::uint64_t reference_epoch() const { return reference_epoch_; }
    std::uint64_t coarse_mesh_version() const { return coarse_mesh_version_; }
    std::uint64_t reference_mesh_version() const {
        return reference_mesh_version_;
    }
    std::uint64_t ambient_mesh_version() const { return ambient_mesh_version_; }
    std::uint64_t candidate_mesh_version() const {
        return ambient_mesh_version_;
    }
    std::uint64_t interpolation_version() const {
        return interpolation_version_;
    }
    std::uint64_t boundary_version() const { return boundary_version_; }
    std::uint64_t corrector_space_version() const {
        return corrector_space_version_;
    }

    ReferenceEpochRefinementResult refine_coarse_preserving_reference(
        const std::vector<int> &marked_elements);
    AmbientRatioEnforcementResult enforce_ambient_ratio(double rho_star);
    bool reference_embedding_holds() const;
    double ambient_ratio() const;

    // Explicit epoch boundary.  This is the only operation that may change
    // reference_mesh() and advance reference_mesh_version().
    void refresh_reference_from_ambient();

    // Transactional hierarchy API used by the revised paper.  A proposed
    // coarse refinement is kept separate from the committed coarse mesh.
    // Candidate enrichment may then make the proposal geometrically
    // representable, while commit remains forbidden until the fixed
    // reference itself contains the proposal.
    void begin_reference_epoch();
    ReferenceEpochCoarseRefinementPreview preview_coarse_refinement(
        const std::vector<int> &marked_elements,
        const ReferenceEpochRefinementGuard &resource_guard = {}) const;
    ReferenceEpochRefinementResult propose_coarse_refinement(
        const std::vector<int> &marked_elements,
        const ReferenceEpochRefinementGuard &resource_guard = {});
    // Adopt a previously accepted preview without repeating its NVB closure
    // or reference-embedding construction.  The candidate embedding is built
    // against the current candidate at adoption time.
    ReferenceEpochRefinementResult propose_coarse_refinement(
        ReferenceEpochCoarseRefinementPreview preview,
        const ReferenceEpochRefinementGuard &resource_guard = {});
    // Open a pending transaction with T_H^+=T_H.  Algorithm 2 needs this
    // when the regular-region indicator mass vanishes but candidate
    // enrichment or an epoch decision must still be completed.
    ReferenceEpochRefinementResult propose_identity_coarse();
    // Refine an already pending proposal without changing the committed
    // coarse mesh.  Existing proposal-to-reference/candidate embeddings are
    // updated from their child-parent maps, so no global geometric embedding
    // search is required.
    ReferenceEpochRefinementResult refine_proposed_coarse(
        const std::vector<int> &marked_proposed_elements,
        const ReferenceEpochRefinementGuard &resource_guard = {});
    ReferenceEpochRefinementResult enrich_candidate(
        const std::vector<int> &marked_candidate_elements,
        const ReferenceEpochRefinementGuard &resource_guard = {});
    ReferenceEpochRefinementResult close_candidate_over_proposed_coarse(
        const ReferenceEpochRefinementGuard &resource_guard = {});
    ReferenceEpochRefinementResult deepen_candidate_over_proposed_coarse(
        int minimum_level_gap,
        const ReferenceEpochRefinementGuard &resource_guard = {});
    ReferenceEpochRefinementResult deepen_candidate_over_proposed_coarse(
        int minimum_level_gap,
        const std::vector<char> &included_proposed_elements,
        const ReferenceEpochRefinementGuard &resource_guard = {});
    ReferenceEpochRefinementResult deepen_candidate_over_proposed_coarse(
        const std::vector<int> &target_level_gaps,
        const ReferenceEpochRefinementGuard &resource_guard = {});
    // Run exactly the same marked-NVB depth loop on candidate mesh/level/
    // parent-map scratch data only.  No hierarchy matrices or state mutate.
    ReferenceEpochCandidateDeepeningProbe
    probe_candidate_deepening_over_proposed_coarse(
        const std::vector<int> &target_level_gaps,
        const std::vector<char> &protected_proposed_elements,
        const std::vector<char> &full_target_proposed_elements,
        const ReferenceEpochRefinementGuard &resource_guard = {}) const;
    bool candidate_contains_proposed_coarse() const;
    bool reference_contains_proposed_coarse() const;
    int minimum_proposed_reference_level_gap() const;
    // Minimum over reference children whose proposed coarse parent is selected.
    // The hybrid method selects only the regular region, thereby excluding
    // the deliberately matching Omega_F cells.  A selected zero gap is a real
    // reserve exhaustion event and is therefore retained in the minimum.
    int minimum_proposed_reference_level_gap(
        const std::vector<char> &included_proposed_elements) const;
    int minimum_proposed_candidate_level_gap() const;
    int minimum_proposed_candidate_level_gap(
        const std::vector<char> &included_proposed_elements) const;
    int minimum_proposed_reference_level_gap_margin(
        const std::vector<int> &target_level_gaps) const;
    int minimum_proposed_candidate_level_gap_margin(
        const std::vector<int> &target_level_gaps) const;
    bool has_proposed_coarse_refinement() const {
        return proposed_coarse_.has_value();
    }
    const TriMesh &proposed_coarse_mesh() const;
    const std::vector<int> &proposed_coarse_levels() const;
    const Eigen::SparseMatrix<double> &
    coarse_elements_to_proposed_coarse() const;
    const Eigen::SparseMatrix<double> &
    coarse_nodes_to_proposed_coarse() const;
    const Eigen::SparseMatrix<double> &
    proposed_coarse_elements_to_reference() const;
    const Eigen::SparseMatrix<double> &
    proposed_coarse_elements_to_candidate() const;
    std::vector<int> reference_parent_proposed_coarse_elements() const;
    std::vector<int> candidate_parent_proposed_coarse_elements() const;
    ReferenceEpochRefinementResult commit_coarse_refinement();
    void refresh_reference_from_candidate();

private:
    TriMesh initial_mesh_;
    TriMesh coarse_mesh_;
    std::vector<int> coarse_levels_;
    int reference_level_ = 0;
    NestedFineMesh reference_completion_;
    RefineOutput reference_to_ambient_;
    std::vector<int> ambient_element_levels_;
    RefineOutput coarse_to_ambient_;
    Eigen::SparseMatrix<double> reference_quasi_interpolation_;
    Eigen::SparseMatrix<double> ambient_quasi_interpolation_;
    std::vector<int> reference_parent_coarse_elements_;
    std::vector<int> ambient_parent_coarse_elements_;
    std::vector<int> ambient_parent_reference_elements_;
    std::uint64_t reference_epoch_ = 0;
    std::uint64_t coarse_mesh_version_ = 0;
    std::uint64_t reference_mesh_version_ = 0;
    std::uint64_t ambient_mesh_version_ = 0;
    std::uint64_t interpolation_version_ = 0;
    std::uint64_t boundary_version_ = 0;
    std::uint64_t corrector_space_version_ = 0;

    struct ProposedCoarseState {
        struct CachedEmbedding {
            Eigen::SparseMatrix<double> P_node;
            Eigen::SparseMatrix<double> P_elem;
            Eigen::SparseMatrix<double> P_dg;
        };
        RefineOutput refinement;
        std::vector<int> levels;
        std::optional<CachedEmbedding> to_reference;
        std::optional<CachedEmbedding> to_candidate;
        bool identity = false;
    };
    std::optional<ProposedCoarseState> proposed_coarse_;

    void refresh_embeddings();
    void refresh_embedding_metadata();
    void validate_current_embeddings() const;
    void validate_incremental_candidate_metadata() const;
    void replace_candidate(
        TriMesh mesh,
        std::vector<int> element_levels);
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

// Update a nested embedding after one NVB refinement of the parent mesh.
// The old child-to-parent element map restricts each geometric containment
// check to the few new descendants of the known old parent, avoiding a global
// parent-element search.  Sparse prolongations are rebuilt in linear time in
// the child mesh size because the parent column layout changed.
RefineOutput update_nested_mesh_embedding_after_parent_refinement(
    const TriMesh &old_parent_mesh,
    const RefineOutput &parent_refinement,
    const TriMesh &child_mesh,
    const std::vector<int> &old_child_parent_elements);

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

#pragma once

#include "helmholtz/adaptive/hierarchy.h"

#include <cstddef>
#include <vector>

namespace lod2d::helmholtz::adaptive {

struct SingularRegionClassification {
    // `l_s` is kept as the public spelling for compatibility, but denotes the
    // manuscript's ell_S in N^{ell_S}(S). Because N^0(S)=S, the seed coarse
    // cells belong to N^1(S); the stored value is therefore one plus their
    // zero-based cell-graph radius. It is independent of `corrector_ell`.
    int corrector_ell = 0;
    int l_s = 0;
    double minimum_physical_radius = 0.0;
    double covered_physical_radius = 0.0;
    std::vector<int> seed_elements;
    std::vector<int> omega_s_elements;
    std::vector<int> omega_f_elements;
    std::vector<int> regular_elements;
    std::vector<char> in_omega_s;
    std::vector<char> in_omega_f;
    std::vector<char> in_regular;
};

// A conforming fine mesh cannot in general have zero level gap on Omega_F
// (exact matching) and a fixed positive gap on every immediately adjacent
// regular cell.  A neutral conformity collar is left unrefined first; beyond
// it the attainable reserve grows by one per coarse element-graph layer and
// saturates at `target_gap`.
struct HybridGradedReserveProfile {
    int neutral_collar_layers = 0;
    int maximum_graph_distance = 0;
    std::size_t full_target_elements = 0;
    std::vector<int> target_level_gaps;
    std::vector<char> at_full_target;
};

HybridGradedReserveProfile make_hybrid_graded_reserve_profile(
    const TriMesh &coarse_mesh,
    const SingularRegionClassification &regions,
    int target_gap,
    int neutral_collar_layers = 1);

SingularRegionClassification classify_singular_regions(
    const TriMesh &coarse_mesh,
    const std::vector<Point2> &singular_points,
    int ell,
    double geometric_tolerance = 1e-12);

// Physical-size variant used by the production E2 hybrid method.  It chooses
// ell_S as the smallest manuscript neighborhood index for which every coarse
// element intersecting B_{minimum_physical_radius}(S) is contained in
// Omega_S=N^{ell_S}(S), and then sets
// Omega_F=N^ell(Omega_S)=N^{ell_S+ell}(S).  Thus the physical core is
// independent of the corrector oversampling width ell.
SingularRegionClassification classify_singular_regions_with_physical_radius(
    const TriMesh &coarse_mesh,
    const std::vector<Point2> &singular_points,
    int ell,
    double minimum_physical_radius,
    double geometric_tolerance = 1e-12);

struct HybridMatchingResult {
    SingularRegionClassification regions;
    std::size_t refinement_rounds = 0;
    std::size_t refined_coarse_elements = 0;
    bool reference_unchanged = true;
    bool matching_holds = false;
};

// Refine only the hybrid coarse mesh (plus conforming closure) until every
// coarse cell in Omega_F is exactly one reference cell.  The fixed reference
// and candidate meshes are never refined by this operation.
HybridMatchingResult restore_hybrid_reference_matching(
    ReferenceEpochHierarchy &hierarchy,
    const std::vector<Point2> &singular_points,
    int ell,
    std::size_t maximum_rounds = 64);

HybridMatchingResult restore_hybrid_reference_matching_with_physical_radius(
    ReferenceEpochHierarchy &hierarchy,
    const std::vector<Point2> &singular_points,
    int ell,
    double minimum_physical_radius,
    std::size_t maximum_rounds = 64);

bool hybrid_reference_matching_holds(
    const ReferenceEpochHierarchy &hierarchy,
    const SingularRegionClassification &regions);

// Dörfler marking relative to the regular-region indicator mass only.
std::vector<int> mark_hybrid_regular_region(
    const std::vector<double> &element_eta_squared,
    const SingularRegionClassification &regions,
    double theta);

// Independent regional Dörfler marking for a split F/R indicator.  Each
// nonzero regional mass receives its own theta-bulk marking, so a small F mass
// cannot be hidden by a much larger regular-region mass.
struct SplitRegionalDoerflerMarking {
    double omega_f_mass = 0.0;
    double regular_mass = 0.0;
    double marked_omega_f_mass = 0.0;
    double marked_regular_mass = 0.0;
    std::vector<int> marked_omega_f_elements;
    std::vector<int> marked_regular_elements;
    std::vector<int> marked_elements;
};

SplitRegionalDoerflerMarking mark_split_regional_doerfler(
    const std::vector<double> &element_eta_squared,
    const std::vector<char> &in_omega_f,
    double theta);

} // namespace lod2d::helmholtz::adaptive

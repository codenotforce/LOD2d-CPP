#pragma once

#include "helmholtz/adaptive/hierarchy.h"

#include <cstddef>
#include <vector>

namespace lod2d::helmholtz::adaptive {

struct SingularRegionClassification {
    // The corrector oversampling radius is only a lower bound.  The actual
    // singular-neighborhood radius l_s is enlarged, when necessary, until
    // N^{2 l_s}(z) covers the requested physical radius.
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

SingularRegionClassification classify_singular_regions(
    const TriMesh &coarse_mesh,
    const std::vector<Point2> &singular_points,
    int ell,
    double geometric_tolerance = 1e-12);

// Physical-size variant used by the production E2 hybrid method.  It chooses
// the smallest l_s >= ell for which every coarse element intersecting
// B_{minimum_physical_radius}(z) is contained in N^{2 l_s}(z).  Consequently
// the matching/transition region cannot collapse in physical size as H is
// locally refined.
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

} // namespace lod2d::helmholtz::adaptive

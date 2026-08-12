#pragma once

#include "helmholtz/operators.h"

#include <filesystem>
#include <optional>
#include <string>

namespace lod2d::helmholtz::experiments {

struct ReferenceSolutionCacheLookup {
    std::string key;
    std::filesystem::path path;
    std::optional<ComplexVector> solution;
};

// The identity string must bind the reference solver implementation and code
// revision. The numerical fingerprint additionally covers the complete mesh,
// assembled operator and load. FNV locates the cache file; an independent
// payload checksum and exact dimensions reject corruption on load.
std::string reference_solution_cache_key(
    const TriMesh &mesh,
    const HelmholtzOperators &operators,
    const ComplexVector &load,
    const std::string &identity);

ReferenceSolutionCacheLookup load_reference_solution_cache(
    const std::filesystem::path &directory,
    const std::string &key,
    Eigen::Index expected_size);

void store_reference_solution_cache(
    const std::filesystem::path &directory,
    const std::string &key,
    const ComplexVector &solution);

} // namespace lod2d::helmholtz::experiments

#pragma once

#include <filesystem>
#include <string_view>

namespace lod2d::helmholtz::experiments {

// WP7 entry point used by the existing paper benchmark when schema_version=5.
// Returns a process-style status code and writes the complete reference-epoch
// output contract below output_directory/run_id.
int run_reference_epoch_paper(
    std::string_view config_json,
    const std::filesystem::path &output_directory,
    const std::filesystem::path &manuscript_baseline,
    bool check,
    bool validate_only);

} // namespace lod2d::helmholtz::experiments

#include "helmholtz/experiments/paper_config.h"

#include <iostream>
#include <stdexcept>
#include <string>

using namespace lod2d::helmholtz::experiments;

namespace {

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

template <class Function>
void require_invalid(Function &&function, const char *message) {
    try {
        function();
    } catch (const std::invalid_argument &) {
        return;
    }
    throw std::runtime_error(message);
}

PaperConfig sample_config() {
    PaperConfig config;
    config.case_id = PaperCase::R2b;
    config.method_id = PaperMethod::SlodMatched;
    config.wavenumber = 16.0;
    config.theta_H = 0.7;
    config.repeat_index = 3;
    config.git_commit = "c8b5520d4f2b1e0f225a64e837dce5220daa1a2d";
    config.build_hash = "gcc-13-release-abc123";
    return config;
}

void verify_registries() {
    require(paper_case_registry().size() == 4, "paper case registry is incomplete");
    require(paper_method_registry().size() == 7, "paper method registry is incomplete");
    require(case_definition(PaperCase::R2a).gaussian_sigma == 1.0 / 32.0,
            "R2a sigma is not frozen to 2^-5");
    require(case_definition(PaperCase::R2b).gaussian_sigma == 1.0 / 64.0,
            "R2b sigma is not frozen to 2^-6");
    require(case_definition(PaperCase::S).has_mixed_boundary,
            "S must declare mixed boundary data");
    require(method_definition(PaperMethod::HlodProxy).diagnostic_only,
            "HLOD-proxy must remain diagnostic");
    require(!method_definition(PaperMethod::HlodProxy).paper_comparator,
            "HLOD-proxy must not enter the paper comparator matrix");
}

void verify_round_trip_and_hash() {
    const PaperConfig original = sample_config();
    const std::string encoded = canonical_json(original);
    const PaperConfig decoded = parse_paper_config(encoded);
    require(decoded == original, "paper config JSON round trip lost fields");
    require(canonical_json(decoded) == encoded, "canonical JSON changed after round trip");
    require(canonical_config_hash(decoded) == canonical_config_hash(original),
            "canonical hash changed after round trip");
    require(make_run_id(original) == make_run_id(decoded),
            "run ID is not deterministic");

    PaperConfig changed = original;
    changed.repeat_index = 4;
    require(canonical_config_hash(changed) != canonical_config_hash(original),
            "canonical hash ignores repeat index");
    require(make_run_id(changed) != make_run_id(original),
            "run ID ignores repeat index");
}

void verify_strict_validation() {
    const std::string encoded = canonical_json(sample_config());
    std::string unknown = encoded;
    unknown.insert(unknown.size() - 1, ",\"future_field\":1");
    require_invalid([&] { (void)parse_paper_config(unknown); },
                    "unknown config field was accepted");

    std::string wrong_version = encoded;
    const std::string needle = "\"schema_version\":1";
    wrong_version.replace(wrong_version.find(needle), needle.size(), "\"schema_version\":2");
    require_invalid([&] { (void)parse_paper_config(wrong_version); },
                    "unknown schema version was accepted");

    PaperConfig invalid = sample_config();
    invalid.wavenumber = 7.0;
    require_invalid([&] { (void)canonical_json(invalid); },
                    "non-protocol wavenumber was accepted");
}

} // namespace

int main() {
    try {
        verify_registries();
        verify_round_trip_and_hash();
        verify_strict_validation();
        std::cout << "Helmholtz paper configuration protocol passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_paper_config failed: " << error.what() << '\n';
        return 1;
    }
}

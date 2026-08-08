#pragma once

#include "helmholtz/quadrature.h"

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace lod2d::helmholtz::experiments {

inline constexpr int paper_schema_version = 1;

enum class PaperCase { R1, R2a, R2b, S };
enum class PaperMethod {
    Calod,
    Hlod,
    SlodPrior,
    SlodMatched,
    Ufem,
    Afem,
    HlodProxy
};

struct CaseDefinition {
    PaperCase id;
    std::string name;
    std::string domain;
    bool has_exact_solution = false;
    bool has_mixed_boundary = false;
    std::optional<double> gaussian_sigma;
};

struct MethodDefinition {
    PaperMethod id;
    std::string name;
    bool paper_comparator = true;
    bool diagnostic_only = false;
    bool adapts_H = false;
    bool adapts_h = false;
    bool adapts_ell = false;
};

struct TolerancePolicy {
    double linear_relative_residual = 1e-10;
    double eigen_relative_residual = 1e-10;
    double interpolation_right_inverse = 1e-9;
    double prolongation_composition = 1e-10;
};

struct PaperConfig {
    int schema_version = paper_schema_version;
    PaperCase case_id = PaperCase::R1;
    PaperMethod method_id = PaperMethod::Calod;
    double wavenumber = 8.0;
    double theta_H = 0.5;
    int timing_repeats = 5;
    int repeat_index = 0;
    std::array<double, 4> relative_energy_targets{{0.1, 0.05, 0.02, 0.01}};
    QuadraturePolicy quadrature;
    TolerancePolicy tolerances;
    std::string git_commit;
    std::string build_hash;
};

const std::vector<CaseDefinition> &paper_case_registry();
const std::vector<MethodDefinition> &paper_method_registry();
const CaseDefinition &case_definition(PaperCase id);
const MethodDefinition &method_definition(PaperMethod id);

std::string_view to_string(PaperCase id);
std::string_view to_string(PaperMethod id);
PaperCase parse_paper_case(std::string_view text);
PaperMethod parse_paper_method(std::string_view text);

void validate_paper_config(const PaperConfig &config);
std::string canonical_json(const PaperConfig &config);
PaperConfig parse_paper_config(std::string_view json);

// Stable FNV-1a 64-bit hash of canonical_json(config), written as 16 lowercase hex digits.
std::string canonical_config_hash(const PaperConfig &config);
std::string make_run_id(const PaperConfig &config);

bool operator==(const TolerancePolicy &lhs, const TolerancePolicy &rhs);
bool operator==(const PaperConfig &lhs, const PaperConfig &rhs);

} // namespace lod2d::helmholtz::experiments

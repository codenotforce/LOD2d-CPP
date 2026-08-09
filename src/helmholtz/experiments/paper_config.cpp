#include "helmholtz/experiments/paper_config.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <variant>

namespace lod2d::helmholtz::experiments {
namespace {

using JsonArray = std::vector<struct JsonValue>;
using JsonObject = std::map<std::string, struct JsonValue>;

struct JsonValue {
    std::variant<double, std::string, JsonArray, JsonObject> value;
};

class JsonParser {
public:
    explicit JsonParser(std::string_view input) : input_(input) {}

    JsonValue parse() {
        JsonValue value = parse_value();
        skip_space();
        if (position_ != input_.size()) fail("trailing content");
        return value;
    }

private:
    JsonValue parse_value() {
        skip_space();
        if (position_ >= input_.size()) fail("unexpected end of input");
        if (input_[position_] == '{') return JsonValue{parse_object()};
        if (input_[position_] == '[') return JsonValue{parse_array()};
        if (input_[position_] == '"') return JsonValue{parse_string()};
        return JsonValue{parse_number()};
    }

    JsonObject parse_object() {
        expect('{');
        JsonObject result;
        skip_space();
        if (consume('}')) return result;
        while (true) {
            skip_space();
            if (position_ >= input_.size() || input_[position_] != '"') {
                fail("object key must be a string");
            }
            std::string key = parse_string();
            expect(':');
            if (!result.emplace(key, parse_value()).second) {
                fail("duplicate object key: " + key);
            }
            skip_space();
            if (consume('}')) break;
            expect(',');
        }
        return result;
    }

    JsonArray parse_array() {
        expect('[');
        JsonArray result;
        skip_space();
        if (consume(']')) return result;
        while (true) {
            result.push_back(parse_value());
            skip_space();
            if (consume(']')) break;
            expect(',');
        }
        return result;
    }

    std::string parse_string() {
        expect('"');
        std::string result;
        while (position_ < input_.size()) {
            const char c = input_[position_++];
            if (c == '"') return result;
            if (c == '\\') {
                if (position_ >= input_.size()) fail("unfinished string escape");
                const char escaped = input_[position_++];
                switch (escaped) {
                case '"': result.push_back('"'); break;
                case '\\': result.push_back('\\'); break;
                case '/': result.push_back('/'); break;
                case 'b': result.push_back('\b'); break;
                case 'f': result.push_back('\f'); break;
                case 'n': result.push_back('\n'); break;
                case 'r': result.push_back('\r'); break;
                case 't': result.push_back('\t'); break;
                default: fail("unsupported string escape");
                }
            } else {
                if (static_cast<unsigned char>(c) < 0x20) fail("control character in string");
                result.push_back(c);
            }
        }
        fail("unterminated string");
    }

    double parse_number() {
        const std::size_t begin = position_;
        if (consume('-') && position_ >= input_.size()) fail("incomplete number");
        if (position_ >= input_.size()) fail("expected JSON number");

        if (input_[position_] == '0') {
            ++position_;
            if (position_ < input_.size()
                && input_[position_] >= '0' && input_[position_] <= '9') {
                fail("leading zero in number");
            }
        } else if (input_[position_] >= '1' && input_[position_] <= '9') {
            do {
                ++position_;
            } while (position_ < input_.size()
                     && input_[position_] >= '0' && input_[position_] <= '9');
        } else {
            fail("expected JSON number");
        }

        if (position_ < input_.size() && input_[position_] == '.') {
            ++position_;
            const std::size_t fraction_begin = position_;
            while (position_ < input_.size()
                   && input_[position_] >= '0' && input_[position_] <= '9') {
                ++position_;
            }
            if (fraction_begin == position_) fail("fraction requires a digit");
        }
        if (position_ < input_.size()
            && (input_[position_] == 'e' || input_[position_] == 'E')) {
            ++position_;
            if (position_ < input_.size()
                && (input_[position_] == '+' || input_[position_] == '-')) {
                ++position_;
            }
            const std::size_t exponent_begin = position_;
            while (position_ < input_.size()
                   && input_[position_] >= '0' && input_[position_] <= '9') {
                ++position_;
            }
            if (exponent_begin == position_) fail("exponent requires a digit");
        }

        double result = 0.0;
        const char *first = input_.data() + begin;
        const char *last = input_.data() + position_;
        const auto parsed = std::from_chars(first, last, result, std::chars_format::general);
        if (parsed.ec != std::errc{} || parsed.ptr != last || !std::isfinite(result)) {
            fail("invalid finite number");
        }
        return result;
    }

    void skip_space() {
        while (position_ < input_.size()
               && (input_[position_] == ' ' || input_[position_] == '\n'
                   || input_[position_] == '\r' || input_[position_] == '\t')) {
            ++position_;
        }
    }

    bool consume(char expected) {
        if (position_ < input_.size() && input_[position_] == expected) {
            ++position_;
            return true;
        }
        return false;
    }

    void expect(char expected) {
        skip_space();
        if (!consume(expected)) fail(std::string("expected '") + expected + "'");
    }

    [[noreturn]] void fail(const std::string &message) const {
        throw std::invalid_argument(
            "invalid paper config JSON at byte " + std::to_string(position_) + ": " + message);
    }

    std::string_view input_;
    std::size_t position_ = 0;
};

const JsonObject &as_object(const JsonValue &value, std::string_view field) {
    const auto *result = std::get_if<JsonObject>(&value.value);
    if (!result) throw std::invalid_argument(std::string(field) + " must be an object");
    return *result;
}

const JsonArray &as_array(const JsonValue &value, std::string_view field) {
    const auto *result = std::get_if<JsonArray>(&value.value);
    if (!result) throw std::invalid_argument(std::string(field) + " must be an array");
    return *result;
}

const std::string &as_string(const JsonValue &value, std::string_view field) {
    const auto *result = std::get_if<std::string>(&value.value);
    if (!result) throw std::invalid_argument(std::string(field) + " must be a string");
    return *result;
}

double as_number(const JsonValue &value, std::string_view field) {
    const auto *result = std::get_if<double>(&value.value);
    if (!result) throw std::invalid_argument(std::string(field) + " must be a number");
    return *result;
}

int as_integer(const JsonValue &value, std::string_view field) {
    const double number = as_number(value, field);
    if (number < static_cast<double>(std::numeric_limits<int>::min())
        || number > static_cast<double>(std::numeric_limits<int>::max())
        || std::floor(number) != number) {
        throw std::invalid_argument(std::string(field) + " must be an integer");
    }
    return static_cast<int>(number);
}

void require_keys(
    const JsonObject &object,
    std::initializer_list<std::string_view> expected,
    std::string_view field) {
    if (object.size() != expected.size()) {
        throw std::invalid_argument(std::string(field) + " has missing or unknown fields");
    }
    for (const std::string_view key : expected) {
        if (!object.contains(std::string(key))) {
            throw std::invalid_argument(
                std::string(field) + " is missing field '" + std::string(key) + "'");
        }
    }
}

const JsonValue &get(const JsonObject &object, std::string_view key) {
    const auto found = object.find(std::string(key));
    if (found == object.end()) throw std::invalid_argument("missing field: " + std::string(key));
    return found->second;
}

std::string json_string(std::string_view value) {
    std::string result = "\"";
    for (const char c : value) {
        switch (c) {
        case '"': result += "\\\""; break;
        case '\\': result += "\\\\"; break;
        case '\b': result += "\\b"; break;
        case '\f': result += "\\f"; break;
        case '\n': result += "\\n"; break;
        case '\r': result += "\\r"; break;
        case '\t': result += "\\t"; break;
        default:
            if (static_cast<unsigned char>(c) < 0x20)
                throw std::invalid_argument("paper metadata string contains a control character");
            result.push_back(c);
            break;
        }
    }
    result.push_back('"');
    return result;
}

std::string number(double value) {
    std::ostringstream stream;
    stream << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
    return stream.str();
}

std::string safe_component(std::string_view value) {
    std::string result;
    for (const unsigned char c : value) {
        result.push_back((c >= '0' && c <= '9') || (c >= 'a' && c <= 'z')
                             || (c >= 'A' && c <= 'Z') || c == '-'
                         ? static_cast<char>(c)
                         : '-');
    }
    return result;
}

} // namespace

const std::vector<CaseDefinition> &paper_case_registry() {
    static const std::vector<CaseDefinition> registry{
        {PaperCase::R1, "smooth-manufactured", "unit-square", true, false, std::nullopt},
        {PaperCase::R2a, "localized-gaussian-2^-5", "unit-square", false, false, 1.0 / 32.0},
        {PaperCase::R2b, "localized-gaussian-2^-6", "unit-square", false, false, 1.0 / 64.0},
        {PaperCase::S, "l-shape-singular", "l-shape", true, true, std::nullopt},
    };
    return registry;
}

const std::vector<MethodDefinition> &paper_method_registry() {
    static const std::vector<MethodDefinition> registry{
        {PaperMethod::Calod, "CALOD", true, false, true, true, true},
        {PaperMethod::Hlod, "HLOD", true, false, true, false, false},
        {PaperMethod::SlodPrior, "SLOD-prior", true, false, false, false, false},
        {PaperMethod::SlodMatched, "SLOD-matched", true, false, false, false, false},
        {PaperMethod::Ufem, "UFEM", true, false, false, false, false},
        {PaperMethod::Afem, "AFEM", true, false, true, false, false},
        {PaperMethod::HlodProxy, "HLOD-proxy", false, true, true, false, false},
    };
    return registry;
}

const CaseDefinition &case_definition(PaperCase id) {
    const auto &registry = paper_case_registry();
    const auto found = std::find_if(registry.begin(), registry.end(), [id](const auto &entry) {
        return entry.id == id;
    });
    if (found == registry.end()) throw std::invalid_argument("unknown paper case enum");
    return *found;
}

const MethodDefinition &method_definition(PaperMethod id) {
    const auto &registry = paper_method_registry();
    const auto found = std::find_if(registry.begin(), registry.end(), [id](const auto &entry) {
        return entry.id == id;
    });
    if (found == registry.end()) throw std::invalid_argument("unknown paper method enum");
    return *found;
}

std::string_view to_string(PaperCase id) {
    switch (id) {
    case PaperCase::R1: return "R1";
    case PaperCase::R2a: return "R2a";
    case PaperCase::R2b: return "R2b";
    case PaperCase::S: return "S";
    }
    throw std::invalid_argument("unknown paper case enum");
}

std::string_view to_string(PaperMethod id) {
    return method_definition(id).name;
}

std::string_view to_string(PaperRunStatus status) {
    switch (status) {
    case PaperRunStatus::Success: return "success";
    case PaperRunStatus::Interrupted: return "interrupted";
    case PaperRunStatus::CensoredWorkLimit: return "censored_work_limit";
    case PaperRunStatus::CensoredMemoryLimit: return "censored_memory_limit";
    case PaperRunStatus::CensoredTimeLimit: return "censored_time_limit";
    case PaperRunStatus::CensoredIterationLimit: return "censored_iteration_limit";
    case PaperRunStatus::LinearAlgebraFailure: return "linear_algebra_failure";
    case PaperRunStatus::CertificateFailure: return "certificate_failure";
    case PaperRunStatus::Unavailable: return "unavailable";
    }
    throw std::invalid_argument("unknown paper run status enum");
}

std::string_view to_string(PaperValueStatus status) {
    switch (status) {
    case PaperValueStatus::Valid: return "valid";
    case PaperValueStatus::NotApplicable: return "not_applicable";
    case PaperValueStatus::NotComputed: return "not_computed";
    case PaperValueStatus::InvalidDenominator: return "invalid_denominator";
    case PaperValueStatus::EnclosureFailed: return "enclosure_failed";
    }
    throw std::invalid_argument("unknown paper value status enum");
}

PaperRunStatus parse_paper_run_status(std::string_view text) {
    for (PaperRunStatus status : {
             PaperRunStatus::Success, PaperRunStatus::Interrupted,
             PaperRunStatus::CensoredWorkLimit, PaperRunStatus::CensoredMemoryLimit,
             PaperRunStatus::CensoredTimeLimit, PaperRunStatus::CensoredIterationLimit,
             PaperRunStatus::LinearAlgebraFailure, PaperRunStatus::CertificateFailure,
             PaperRunStatus::Unavailable}) {
        if (to_string(status) == text) return status;
    }
    throw std::invalid_argument("unknown paper run status: " + std::string(text));
}

PaperValueStatus parse_paper_value_status(std::string_view text) {
    for (PaperValueStatus status : {
             PaperValueStatus::Valid, PaperValueStatus::NotApplicable,
             PaperValueStatus::NotComputed, PaperValueStatus::InvalidDenominator,
             PaperValueStatus::EnclosureFailed}) {
        if (to_string(status) == text) return status;
    }
    throw std::invalid_argument("unknown paper value status: " + std::string(text));
}

PaperCase parse_paper_case(std::string_view text) {
    for (const auto &entry : paper_case_registry()) {
        if (to_string(entry.id) == text) return entry.id;
    }
    throw std::invalid_argument("unknown paper case: " + std::string(text));
}

PaperMethod parse_paper_method(std::string_view text) {
    for (const auto &entry : paper_method_registry()) {
        if (entry.name == text) return entry.id;
    }
    throw std::invalid_argument("unknown paper method: " + std::string(text));
}

void validate_paper_config(const PaperConfig &config) {
    if (config.schema_version != paper_schema_version) {
        throw std::invalid_argument(
            "unsupported paper config schema version: " + std::to_string(config.schema_version));
    }
    (void)case_definition(config.case_id);
    (void)method_definition(config.method_id);
    if (!(config.wavenumber == 8.0 || config.wavenumber == 16.0
          || config.wavenumber == 32.0)) {
        throw std::invalid_argument("paper wavenumber must be 8, 16, or 32");
    }
    if (!(config.theta_H == 0.3 || config.theta_H == 0.5 || config.theta_H == 0.7)) {
        throw std::invalid_argument("theta_H must be 0.3, 0.5, or 0.7");
    }
    if (config.timing_repeats != 5 || config.repeat_index < 0
        || config.repeat_index >= config.timing_repeats) {
        throw std::invalid_argument("formal timing requires repeat_index in [0,5)");
    }
    constexpr std::array<double, 4> targets{{0.1, 0.05, 0.02, 0.01}};
    if (config.relative_energy_targets != targets) {
        throw std::invalid_argument("relative energy targets differ from the frozen protocol");
    }
    if (config.quadrature.base_triangle_order <= 0
        || config.quadrature.gaussian_triangle_order <= 0
        || config.quadrature.singular_triangle_order <= 0
        || config.quadrature.max_recursive_subdivisions < 0) {
        throw std::invalid_argument("quadrature policy contains invalid orders");
    }
    const std::array<double, 4> tolerances{
        config.tolerances.linear_relative_residual,
        config.tolerances.eigen_relative_residual,
        config.tolerances.interpolation_right_inverse,
        config.tolerances.prolongation_composition};
    if (std::any_of(tolerances.begin(), tolerances.end(), [](double value) {
            return !std::isfinite(value) || !(value > 0.0);
        })) {
        throw std::invalid_argument("tolerances must be positive finite numbers");
    }
    if (config.git_commit.empty() || config.build_hash.empty()) {
        throw std::invalid_argument("git_commit and build_hash are required for immutable runs");
    }
    const auto contains_control = [](const std::string &value) {
        return std::any_of(value.begin(), value.end(), [](unsigned char character) {
            return character < 0x20;
        });
    };
    if (contains_control(config.git_commit) || contains_control(config.build_hash))
        throw std::invalid_argument("git_commit and build_hash cannot contain control characters");
}

std::string canonical_json(const PaperConfig &config) {
    validate_paper_config(config);
    std::ostringstream out;
    out << "{\"build_hash\":" << json_string(config.build_hash)
        << ",\"case\":" << json_string(to_string(config.case_id))
        << ",\"git_commit\":" << json_string(config.git_commit)
        << ",\"method\":" << json_string(to_string(config.method_id))
        << ",\"quadrature\":{\"base_triangle_order\":" << config.quadrature.base_triangle_order
        << ",\"gaussian_triangle_order\":" << config.quadrature.gaussian_triangle_order
        << ",\"max_recursive_subdivisions\":" << config.quadrature.max_recursive_subdivisions
        << ",\"singular_triangle_order\":" << config.quadrature.singular_triangle_order << "}"
        << ",\"relative_energy_targets\":[";
    for (std::size_t i = 0; i < config.relative_energy_targets.size(); ++i) {
        if (i) out << ',';
        out << number(config.relative_energy_targets[i]);
    }
    out << "]"
        << ",\"repeat_index\":" << config.repeat_index
        << ",\"schema_version\":" << config.schema_version
        << ",\"theta_H\":" << number(config.theta_H)
        << ",\"timing_repeats\":" << config.timing_repeats
        << ",\"tolerances\":{\"eigen_relative_residual\":"
        << number(config.tolerances.eigen_relative_residual)
        << ",\"interpolation_right_inverse\":"
        << number(config.tolerances.interpolation_right_inverse)
        << ",\"linear_relative_residual\":"
        << number(config.tolerances.linear_relative_residual)
        << ",\"prolongation_composition\":"
        << number(config.tolerances.prolongation_composition) << "}"
        << ",\"wavenumber\":" << number(config.wavenumber) << '}';
    return out.str();
}

PaperConfig parse_paper_config(std::string_view json) {
    const JsonValue parsed_root = JsonParser(json).parse();
    const JsonObject &root = as_object(parsed_root, "root");
    require_keys(root,
        {"build_hash", "case", "git_commit", "method", "quadrature",
         "relative_energy_targets", "repeat_index", "schema_version", "theta_H",
         "timing_repeats", "tolerances", "wavenumber"},
        "root");
    PaperConfig config;
    config.build_hash = as_string(get(root, "build_hash"), "build_hash");
    config.case_id = parse_paper_case(as_string(get(root, "case"), "case"));
    config.git_commit = as_string(get(root, "git_commit"), "git_commit");
    config.method_id = parse_paper_method(as_string(get(root, "method"), "method"));
    config.repeat_index = as_integer(get(root, "repeat_index"), "repeat_index");
    config.schema_version = as_integer(get(root, "schema_version"), "schema_version");
    config.theta_H = as_number(get(root, "theta_H"), "theta_H");
    config.timing_repeats = as_integer(get(root, "timing_repeats"), "timing_repeats");
    config.wavenumber = as_number(get(root, "wavenumber"), "wavenumber");

    const JsonArray &targets = as_array(get(root, "relative_energy_targets"), "relative_energy_targets");
    if (targets.size() != config.relative_energy_targets.size()) {
        throw std::invalid_argument("relative_energy_targets must contain four values");
    }
    for (std::size_t i = 0; i < targets.size(); ++i) {
        config.relative_energy_targets[i] = as_number(targets[i], "relative_energy_targets");
    }

    const JsonObject &quadrature = as_object(get(root, "quadrature"), "quadrature");
    require_keys(quadrature,
        {"base_triangle_order", "gaussian_triangle_order", "max_recursive_subdivisions",
         "singular_triangle_order"},
        "quadrature");
    config.quadrature.base_triangle_order = as_integer(
        get(quadrature, "base_triangle_order"), "base_triangle_order");
    config.quadrature.gaussian_triangle_order = as_integer(
        get(quadrature, "gaussian_triangle_order"), "gaussian_triangle_order");
    config.quadrature.max_recursive_subdivisions = as_integer(
        get(quadrature, "max_recursive_subdivisions"), "max_recursive_subdivisions");
    config.quadrature.singular_triangle_order = as_integer(
        get(quadrature, "singular_triangle_order"), "singular_triangle_order");

    const JsonObject &tolerances = as_object(get(root, "tolerances"), "tolerances");
    require_keys(tolerances,
        {"eigen_relative_residual", "interpolation_right_inverse", "linear_relative_residual",
         "prolongation_composition"},
        "tolerances");
    config.tolerances.eigen_relative_residual = as_number(
        get(tolerances, "eigen_relative_residual"), "eigen_relative_residual");
    config.tolerances.interpolation_right_inverse = as_number(
        get(tolerances, "interpolation_right_inverse"), "interpolation_right_inverse");
    config.tolerances.linear_relative_residual = as_number(
        get(tolerances, "linear_relative_residual"), "linear_relative_residual");
    config.tolerances.prolongation_composition = as_number(
        get(tolerances, "prolongation_composition"), "prolongation_composition");
    validate_paper_config(config);
    return config;
}

std::string canonical_config_hash(const PaperConfig &config) {
    const std::string canonical = canonical_json(config);
    std::uint64_t hash = 14695981039346656037ULL;
    for (const unsigned char byte : canonical) {
        hash ^= byte;
        hash *= 1099511628211ULL;
    }
    std::ostringstream out;
    out << std::hex << std::setfill('0') << std::setw(16) << hash;
    return out.str();
}

std::string make_run_id(const PaperConfig &config) {
    validate_paper_config(config);
    std::ostringstream out;
    out << to_string(config.case_id) << '_' << safe_component(to_string(config.method_id))
        << "_k" << static_cast<int>(config.wavenumber)
        << "_r" << config.repeat_index << '_' << canonical_config_hash(config);
    return out.str();
}

bool operator==(const TolerancePolicy &lhs, const TolerancePolicy &rhs) {
    return lhs.linear_relative_residual == rhs.linear_relative_residual
        && lhs.eigen_relative_residual == rhs.eigen_relative_residual
        && lhs.interpolation_right_inverse == rhs.interpolation_right_inverse
        && lhs.prolongation_composition == rhs.prolongation_composition;
}

bool operator==(const PaperConfig &lhs, const PaperConfig &rhs) {
    return lhs.schema_version == rhs.schema_version && lhs.case_id == rhs.case_id
        && lhs.method_id == rhs.method_id && lhs.wavenumber == rhs.wavenumber
        && lhs.theta_H == rhs.theta_H && lhs.timing_repeats == rhs.timing_repeats
        && lhs.repeat_index == rhs.repeat_index
        && lhs.relative_energy_targets == rhs.relative_energy_targets
        && lhs.quadrature == rhs.quadrature && lhs.tolerances == rhs.tolerances
        && lhs.git_commit == rhs.git_commit && lhs.build_hash == rhs.build_hash;
}

} // namespace lod2d::helmholtz::experiments

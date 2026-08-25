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

bool valid_sha256_digest(std::string_view value) {
    constexpr std::string_view prefix = "sha256:";
    return value.size() == prefix.size() + 64
        && value.starts_with(prefix)
        && std::all_of(
            value.begin() + static_cast<std::ptrdiff_t>(prefix.size()),
            value.end(),
            [](char character) {
                return (character >= '0' && character <= '9')
                    || (character >= 'a' && character <= 'f');
            });
}

struct JsonValue {
    std::variant<bool, double, std::string, JsonArray, JsonObject> value;
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
        if (input_.substr(position_, 4) == "true") {
            position_ += 4;
            return JsonValue{true};
        }
        if (input_.substr(position_, 5) == "false") {
            position_ += 5;
            return JsonValue{false};
        }
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

bool as_bool(const JsonValue &value, std::string_view field) {
    const auto *result = std::get_if<bool>(&value.value);
    if (!result) throw std::invalid_argument(std::string(field) + " must be a boolean");
    return *result;
}

constexpr std::uint64_t max_exact_json_integer = 9007199254740991ULL;

std::uint64_t as_uint64(const JsonValue &value, std::string_view field) {
    const double number = as_number(value, field);
    if (number < 0.0 || number > static_cast<double>(max_exact_json_integer)
        || std::floor(number) != number) {
        throw std::invalid_argument(
            std::string(field) + " must be a nonnegative JSON-safe integer");
    }
    return static_cast<std::uint64_t>(number);
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

adaptive::CertifiedErrorTarget parse_error_target(std::string_view text) {
    for (const auto value : {
             adaptive::CertifiedErrorTarget::AuditSpace,
             adaptive::CertifiedErrorTarget::Continuous}) {
        if (adaptive::to_string(value) == text) return value;
    }
    throw std::invalid_argument("unknown certified error target: " + std::string(text));
}

adaptive::CertifiedEvidencePolicy parse_evidence_policy(std::string_view text) {
    for (const auto value : {
             adaptive::CertifiedEvidencePolicy::RequireVerified,
             adaptive::CertifiedEvidencePolicy::AllowConditional}) {
        if (adaptive::to_string(value) == text) return value;
    }
    throw std::invalid_argument("unknown certified evidence policy: " + std::string(text));
}

std::string_view petrov_mode_name(HelmholtzPetrovMode value) {
    return value == HelmholtzPetrovMode::TwoSided
        ? std::string_view("two_sided")
        : std::string_view("corrected_test_only");
}

HelmholtzPetrovMode parse_petrov_mode(std::string_view text) {
    if (text == "two_sided") return HelmholtzPetrovMode::TwoSided;
    if (text == "corrected_test_only")
        return HelmholtzPetrovMode::CorrectedTestOnly;
    throw std::invalid_argument("unknown Helmholtz Petrov mode: " + std::string(text));
}

std::string_view patch_solver_kind_name(HelmholtzPatchSolverKind value) {
    switch (value) {
    case HelmholtzPatchSolverKind::DirectSaddle: return "direct_saddle";
    case HelmholtzPatchSolverKind::DirectSchur: return "direct_schur";
    case HelmholtzPatchSolverKind::ShiftedGmres: return "shifted_gmres";
    }
    throw std::invalid_argument("unknown Helmholtz patch solver kind");
}

HelmholtzPatchSolverKind parse_patch_solver_kind(std::string_view text) {
    if (text == "direct_saddle") return HelmholtzPatchSolverKind::DirectSaddle;
    if (text == "direct_schur") return HelmholtzPatchSolverKind::DirectSchur;
    if (text == "shifted_gmres") return HelmholtzPatchSolverKind::ShiftedGmres;
    throw std::invalid_argument("unknown Helmholtz patch solver kind: " + std::string(text));
}

std::string_view shift_rule_name(HelmholtzShiftRule value) {
    switch (value) {
    case HelmholtzShiftRule::KappaSquared: return "kappa_squared";
    case HelmholtzShiftRule::PatchScaled: return "patch_scaled";
    case HelmholtzShiftRule::Absolute: return "absolute";
    }
    throw std::invalid_argument("unknown Helmholtz shift rule");
}

HelmholtzShiftRule parse_shift_rule(std::string_view text) {
    if (text == "kappa_squared") return HelmholtzShiftRule::KappaSquared;
    if (text == "patch_scaled") return HelmholtzShiftRule::PatchScaled;
    if (text == "absolute") return HelmholtzShiftRule::Absolute;
    throw std::invalid_argument("unknown Helmholtz shift rule: " + std::string(text));
}

std::string_view shifted_inverse_name(HelmholtzShiftedInverseKind value) {
    switch (value) {
    case HelmholtzShiftedInverseKind::Identity: return "identity";
    case HelmholtzShiftedInverseKind::SparseLu: return "sparse_lu";
    case HelmholtzShiftedInverseKind::GeometricVcycle: return "geometric_vcycle";
    }
    throw std::invalid_argument("unknown shifted inverse kind");
}

HelmholtzShiftedInverseKind parse_shifted_inverse(std::string_view text) {
    if (text == "identity") return HelmholtzShiftedInverseKind::Identity;
    if (text == "sparse_lu") return HelmholtzShiftedInverseKind::SparseLu;
    if (text == "geometric_vcycle")
        return HelmholtzShiftedInverseKind::GeometricVcycle;
    throw std::invalid_argument("unknown shifted inverse kind: " + std::string(text));
}

std::string_view kernel_solver_name(adaptive::KernelRieszSolver value) {
    return value == adaptive::KernelRieszSolver::SaddlePoint
        ? std::string_view("saddle_point")
        : std::string_view("kernel_basis_reference");
}

adaptive::KernelRieszSolver parse_kernel_solver(std::string_view text) {
    if (text == "saddle_point") return adaptive::KernelRieszSolver::SaddlePoint;
    if (text == "kernel_basis_reference")
        return adaptive::KernelRieszSolver::KernelBasisReference;
    throw std::invalid_argument("unknown kernel Riesz solver: " + std::string(text));
}

bool equal_work_limits(
    const adaptive::CertifiedWorkLimits &lhs,
    const adaptive::CertifiedWorkLimits &rhs) {
    return lhs.max_state_transitions == rhs.max_state_transitions
        && lhs.max_coarse_refinements == rhs.max_coarse_refinements
        && lhs.max_corrector_refinements == rhs.max_corrector_refinements
        && lhs.max_oversampling_increments == rhs.max_oversampling_increments
        && lhs.max_audit_refinements == rhs.max_audit_refinements
        && lhs.max_coarse_dofs == rhs.max_coarse_dofs
        && lhs.max_fine_dofs == rhs.max_fine_dofs
        && lhs.max_audit_dofs == rhs.max_audit_dofs
        && lhs.max_backend_work_units == rhs.max_backend_work_units
        && lhs.max_peak_memory_bytes == rhs.max_peak_memory_bytes
        && lhs.max_elapsed_seconds == rhs.max_elapsed_seconds
        && lhs.max_oversampling == rhs.max_oversampling;
}

} // namespace

const std::vector<CaseDefinition> &paper_case_registry() {
    static const std::vector<CaseDefinition> registry{
        {PaperCase::R1, "localized-smooth-oscillatory", "unit-square", true, true, std::nullopt},
        {PaperCase::R2a, "localized-gaussian-2^-5", "unit-square", false, false, 1.0 / 32.0},
        {PaperCase::R2b, "localized-gaussian-2^-6", "unit-square", false, false, 1.0 / 64.0},
        {PaperCase::S, "l-shape-singular-interior-wave", "l-shape", true, true, std::nullopt},
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

std::string_view to_string(PaperCertificateStatus status) {
    switch (status) {
    case PaperCertificateStatus::ImplementationStudy: return "implementation-study";
    case PaperCertificateStatus::Conditional: return "conditional";
    case PaperCertificateStatus::AuditCertified: return "audit-certified";
    case PaperCertificateStatus::ContinuousCertified: return "continuous-certified";
    case PaperCertificateStatus::EmpiricalReference: return "empirical-reference";
    }
    throw std::invalid_argument("unknown paper certificate status enum");
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

PaperCertificateStatus parse_paper_certificate_status(std::string_view text) {
    for (PaperCertificateStatus status : {
             PaperCertificateStatus::ImplementationStudy,
             PaperCertificateStatus::Conditional,
             PaperCertificateStatus::AuditCertified,
             PaperCertificateStatus::ContinuousCertified,
             PaperCertificateStatus::EmpiricalReference}) {
        if (to_string(status) == text) return status;
    }
    throw std::invalid_argument(
        "unknown paper certificate status: " + std::string(text));
}

const std::vector<std::string_view> &paper_output_metric_registry() {
    static const std::vector<std::string_view> registry{
        "reference_energy_error",
        "reference_l2_error",
        "exact_energy_error",
        "exact_l2_error",
        "cert_audit_error",
        "eta_H",
        "theta_total_lower",
        "theta_total_upper",
        "theta_h_upper",
        "delta_total_lower",
        "delta_total_upper",
        "delta_h_upper",
        "delta_ell_lower",
        "delta_ell_upper",
        "q_total",
        "q_h",
        "q_ell",
        "mu_H",
        "inf_sup_lower",
        "stability_margin",
        "lod_error_lower",
        "lod_error_upper",
        "audit_error_lower",
        "audit_error_upper",
        "true_error_lower",
        "true_error_upper",
        "upper_effectivity",
        "certificate_gap",
        "patch_work_cumulative",
        "peak_memory_bytes",
    };
    return registry;
}

void validate_paper_output_metric(std::string_view name) {
    const auto &registry = paper_output_metric_registry();
    if (std::find(registry.begin(), registry.end(), name) == registry.end()) {
        throw std::invalid_argument("unknown paper output metric: " + std::string(name));
    }
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
    (void)adaptive::to_string(config.error_target);
    (void)adaptive::to_string(config.evidence_policy);
    if (!(std::isfinite(config.mu0) && config.mu0 > 0.0 && config.mu0 < 1.0))
        throw std::invalid_argument("mu0 must lie in (0,1)");
    if (!(std::isfinite(config.q0) && config.q0 > 0.0))
        throw std::invalid_argument("q0 must be positive and finite");
    if (!(std::isfinite(config.tau) && config.tau > 0.0 && config.tau < 1.0))
        throw std::invalid_argument("tau must lie in (0,1)");
    if (!(std::isfinite(config.theta_h) && config.theta_h > 0.0
          && config.theta_h <= 1.0))
        throw std::invalid_argument("theta_h must lie in (0,1]");
    if (!(std::isfinite(config.rho_aud) && config.rho_aud > 0.0
          && config.rho_aud < 1.0))
        throw std::invalid_argument("rho_aud must lie in (0,1)");
    if (!(std::isfinite(config.tolerance) && config.tolerance > 0.0))
        throw std::invalid_argument("tolerance must be positive and finite");
    const NumericalBackendPolicy &backend = config.numerical_backend;
    if (backend.initial_coarse_level < 0
        || backend.initial_fine_level < backend.initial_coarse_level
        || backend.initial_oversampling < 0) {
        throw std::invalid_argument(
            "numerical backend requires 0 <= initial_coarse_level <= initial_fine_level and initial_oversampling >= 0");
    }
    if (!(std::isfinite(backend.boundary_beta) && backend.boundary_beta >= 0.0))
        throw std::invalid_argument("numerical backend boundary_beta is invalid");
    (void)petrov_mode_name(backend.petrov_mode);
    (void)patch_solver_kind_name(backend.patch_solver.kind);
    (void)shift_rule_name(backend.patch_solver.shifted.rule);
    (void)shifted_inverse_name(backend.patch_solver.shifted.inverse);
    (void)kernel_solver_name(backend.kernel_riesz_solver);
    if (backend.patch_solver.symbolic_cache_slots <= 0
        || backend.patch_solver.maximum_parallel_solves < 0
        || backend.patch_solver.gmres.restart <= 0
        || backend.patch_solver.gmres.max_iterations <= 0
        || !(std::isfinite(backend.patch_solver.gmres.relative_tolerance)
             && backend.patch_solver.gmres.relative_tolerance > 0.0)
        || !(std::isfinite(backend.patch_solver.gmres.absolute_tolerance)
             && backend.patch_solver.gmres.absolute_tolerance >= 0.0)
        || !(std::isfinite(backend.patch_solver.shifted.alpha)
             && backend.patch_solver.shifted.alpha >= 0.0)
        || !(std::isfinite(backend.patch_solver.shifted.absolute_epsilon)
             && backend.patch_solver.shifted.absolute_epsilon >= 0.0)
        || backend.patch_solver.shifted.pre_smooth < 0
        || backend.patch_solver.shifted.post_smooth < 0
        || backend.patch_solver.shifted.coarse_max_dofs <= 0
        || !(std::isfinite(backend.patch_solver.shifted.jacobi_weight)
             && backend.patch_solver.shifted.jacobi_weight > 0.0)) {
        throw std::invalid_argument("numerical backend patch solver policy is invalid");
    }
    if (backend.certificate.precision_bits < 64
        || !(std::isfinite(backend.certificate.cluster_relative_gap)
             && backend.certificate.cluster_relative_gap >= 0.0)
        || !(std::isfinite(backend.certificate.cluster_absolute_gap)
             && backend.certificate.cluster_absolute_gap >= 0.0)
        || !(std::isfinite(backend.certificate.conjugation_tolerance)
             && backend.certificate.conjugation_tolerance >= 0.0)
        || backend.certificate.q0 != config.q0) {
        throw std::invalid_argument(
            "numerical backend certificate policy is invalid or its q0 differs from the driver q0");
    }
    if (!(std::isfinite(backend.audit_doerfler_theta)
          && backend.audit_doerfler_theta > 0.0
          && backend.audit_doerfler_theta <= 1.0)
        || !(std::isfinite(backend.audit_saturation_factor)
             && backend.audit_saturation_factor >= 0.0
             && backend.audit_saturation_factor < 1.0)) {
        throw std::invalid_argument("numerical backend audit policy is invalid");
    }
    if (!valid_sha256_digest(backend.certificate_constant_set_hash)) {
        throw std::invalid_argument(
            "certificate_constant_set_hash must be sha256:<64 lowercase hex digits>");
    }
    if (config.work_limits.max_state_transitions == 0)
        throw std::invalid_argument("max_state_transitions must be positive");
    if (!std::isfinite(config.work_limits.max_elapsed_seconds)
        || config.work_limits.max_elapsed_seconds < 0.0)
        throw std::invalid_argument("max_elapsed_seconds is invalid");
    if (config.work_limits.max_oversampling < -1)
        throw std::invalid_argument("max_oversampling is invalid");
    const std::array<std::uint64_t, 10> integer_limits{
        config.work_limits.max_state_transitions,
        config.work_limits.max_coarse_refinements,
        config.work_limits.max_corrector_refinements,
        config.work_limits.max_oversampling_increments,
        config.work_limits.max_audit_refinements,
        config.work_limits.max_coarse_dofs,
        config.work_limits.max_fine_dofs,
        config.work_limits.max_audit_dofs,
        config.work_limits.max_backend_work_units,
        config.work_limits.max_peak_memory_bytes,
    };
    if (std::any_of(integer_limits.begin(), integer_limits.end(), [](std::uint64_t value) {
            return value > max_exact_json_integer;
        })) {
        throw std::invalid_argument("work limits must be JSON-safe integers");
    }
    if (config.method_id == PaperMethod::Hlod
        && (config.hlod_prior_corrector_space_id.empty()
            || config.hlod_prior_oversampling < 0)) {
        throw std::invalid_argument(
            "HLOD requires a frozen corrector-space id and oversampling value");
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
    if (contains_control(config.git_commit) || contains_control(config.build_hash)
        || contains_control(config.numerical_backend.certificate_constant_set_hash)
        || contains_control(config.hlod_prior_corrector_space_id)) {
        throw std::invalid_argument(
            "provenance and HLOD prior identifiers cannot contain control characters");
    }
}

adaptive::CertifiedDriverConfig make_certified_driver_config(
    const PaperConfig &config) {
    validate_paper_config(config);
    adaptive::CertifiedDriverConfig result;
    if (config.method_id == PaperMethod::Calod) {
        result.method = adaptive::CertifiedMethod::Calod;
    } else if (config.method_id == PaperMethod::Hlod) {
        result.method = adaptive::CertifiedMethod::Hlod;
    } else {
        throw std::invalid_argument(
            "only CALOD and HLOD paper methods use CertifiedAdaptiveDriver");
    }
    result.error_target = config.error_target;
    result.evidence_policy = config.evidence_policy;
    result.mu0 = config.mu0;
    result.q0 = config.q0;
    result.tau = config.tau;
    result.theta_H = config.theta_H;
    result.theta_h = config.theta_h;
    result.rho_aud = config.rho_aud;
    result.tolerance = config.tolerance;
    result.limits = config.work_limits;
    result.hlod_prior_corrector_space_id = config.hlod_prior_corrector_space_id;
    result.hlod_prior_oversampling = config.hlod_prior_oversampling;
    return result;
}

std::string canonical_json(const PaperConfig &config) {
    validate_paper_config(config);
    std::ostringstream out;
    out << "{\"build_hash\":" << json_string(config.build_hash)
        << ",\"case\":" << json_string(to_string(config.case_id))
        << ",\"error_target\":"
        << json_string(adaptive::to_string(config.error_target))
        << ",\"evidence_policy\":"
        << json_string(adaptive::to_string(config.evidence_policy))
        << ",\"git_commit\":" << json_string(config.git_commit)
        << ",\"hlod_prior_corrector_space_id\":"
        << json_string(config.hlod_prior_corrector_space_id)
        << ",\"hlod_prior_oversampling\":" << config.hlod_prior_oversampling
        << ",\"method\":" << json_string(to_string(config.method_id))
        << ",\"mu0\":" << number(config.mu0)
        << ",\"numerical_backend\":{\"audit_doerfler_theta\":"
        << number(config.numerical_backend.audit_doerfler_theta)
        << ",\"audit_saturation_factor\":"
        << number(config.numerical_backend.audit_saturation_factor)
        << ",\"boundary_beta\":" << number(config.numerical_backend.boundary_beta)
        << ",\"certificate\":{\"cluster_absolute_gap\":"
        << number(config.numerical_backend.certificate.cluster_absolute_gap)
        << ",\"cluster_relative_gap\":"
        << number(config.numerical_backend.certificate.cluster_relative_gap)
        << ",\"conjugation_tolerance\":"
        << number(config.numerical_backend.certificate.conjugation_tolerance)
        << ",\"precision_bits\":"
        << config.numerical_backend.certificate.precision_bits
        << ",\"q0\":" << number(config.numerical_backend.certificate.q0) << "}"
        << ",\"certificate_constant_set_hash\":"
        << json_string(config.numerical_backend.certificate_constant_set_hash)
        << ",\"initial_coarse_level\":"
        << config.numerical_backend.initial_coarse_level
        << ",\"initial_fine_level\":" << config.numerical_backend.initial_fine_level
        << ",\"initial_oversampling\":"
        << config.numerical_backend.initial_oversampling
        << ",\"kernel_riesz_solver\":"
        << json_string(kernel_solver_name(config.numerical_backend.kernel_riesz_solver))
        << ",\"patch_solver\":{\"fallback_to_direct\":"
        << (config.numerical_backend.patch_solver.fallback_to_direct ? "true" : "false")
        << ",\"gmres\":{\"absolute_tolerance\":"
        << number(config.numerical_backend.patch_solver.gmres.absolute_tolerance)
        << ",\"max_iterations\":"
        << config.numerical_backend.patch_solver.gmres.max_iterations
        << ",\"relative_tolerance\":"
        << number(config.numerical_backend.patch_solver.gmres.relative_tolerance)
        << ",\"reorthogonalize\":"
        << (config.numerical_backend.patch_solver.gmres.reorthogonalize ? "true" : "false")
        << ",\"restart\":" << config.numerical_backend.patch_solver.gmres.restart << "}"
        << ",\"kind\":"
        << json_string(patch_solver_kind_name(config.numerical_backend.patch_solver.kind))
        << ",\"maximum_parallel_solves\":"
        << config.numerical_backend.patch_solver.maximum_parallel_solves
        << ",\"reuse_identical_factorization\":"
        << (config.numerical_backend.patch_solver.reuse_identical_factorization
                ? "true" : "false")
        << ",\"shifted\":{\"absolute_epsilon\":"
        << number(config.numerical_backend.patch_solver.shifted.absolute_epsilon)
        << ",\"alpha\":" << number(config.numerical_backend.patch_solver.shifted.alpha)
        << ",\"coarse_max_dofs\":"
        << config.numerical_backend.patch_solver.shifted.coarse_max_dofs
        << ",\"inverse\":"
        << json_string(shifted_inverse_name(
               config.numerical_backend.patch_solver.shifted.inverse))
        << ",\"jacobi_weight\":"
        << number(config.numerical_backend.patch_solver.shifted.jacobi_weight)
        << ",\"post_smooth\":"
        << config.numerical_backend.patch_solver.shifted.post_smooth
        << ",\"pre_smooth\":" << config.numerical_backend.patch_solver.shifted.pre_smooth
        << ",\"rule\":"
        << json_string(shift_rule_name(config.numerical_backend.patch_solver.shifted.rule))
        << "},\"symbolic_cache_slots\":"
        << config.numerical_backend.patch_solver.symbolic_cache_slots << "}"
        << ",\"petrov_mode\":"
        << json_string(petrov_mode_name(config.numerical_backend.petrov_mode)) << "}"
        << ",\"q0\":" << number(config.q0)
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
        << ",\"rho_aud\":" << number(config.rho_aud)
        << ",\"schema_version\":" << config.schema_version
        << ",\"tau\":" << number(config.tau)
        << ",\"theta_H\":" << number(config.theta_H)
        << ",\"theta_h\":" << number(config.theta_h)
        << ",\"timing_repeats\":" << config.timing_repeats
        << ",\"tolerance\":" << number(config.tolerance)
        << ",\"tolerances\":{\"eigen_relative_residual\":"
        << number(config.tolerances.eigen_relative_residual)
        << ",\"interpolation_right_inverse\":"
        << number(config.tolerances.interpolation_right_inverse)
        << ",\"linear_relative_residual\":"
        << number(config.tolerances.linear_relative_residual)
        << ",\"prolongation_composition\":"
        << number(config.tolerances.prolongation_composition) << "}"
        << ",\"wavenumber\":" << number(config.wavenumber)
        << ",\"work_limits\":{\"max_audit_dofs\":"
        << config.work_limits.max_audit_dofs
        << ",\"max_audit_refinements\":"
        << config.work_limits.max_audit_refinements
        << ",\"max_backend_work_units\":"
        << config.work_limits.max_backend_work_units
        << ",\"max_coarse_dofs\":" << config.work_limits.max_coarse_dofs
        << ",\"max_coarse_refinements\":"
        << config.work_limits.max_coarse_refinements
        << ",\"max_corrector_refinements\":"
        << config.work_limits.max_corrector_refinements
        << ",\"max_elapsed_seconds\":"
        << number(config.work_limits.max_elapsed_seconds)
        << ",\"max_fine_dofs\":" << config.work_limits.max_fine_dofs
        << ",\"max_oversampling\":" << config.work_limits.max_oversampling
        << ",\"max_oversampling_increments\":"
        << config.work_limits.max_oversampling_increments
        << ",\"max_peak_memory_bytes\":"
        << config.work_limits.max_peak_memory_bytes
        << ",\"max_state_transitions\":"
        << config.work_limits.max_state_transitions << "}}";
    return out.str();
}

PaperConfig parse_paper_config(std::string_view json) {
    const JsonValue parsed_root = JsonParser(json).parse();
    const JsonObject &root = as_object(parsed_root, "root");
    require_keys(root,
        {"build_hash", "case", "error_target", "evidence_policy", "git_commit",
         "hlod_prior_corrector_space_id", "hlod_prior_oversampling", "method",
         "mu0", "numerical_backend", "q0", "quadrature", "relative_energy_targets", "repeat_index",
         "rho_aud", "schema_version", "tau", "theta_H", "theta_h",
         "timing_repeats", "tolerance", "tolerances", "wavenumber", "work_limits"},
        "root");
    PaperConfig config;
    config.build_hash = as_string(get(root, "build_hash"), "build_hash");
    config.case_id = parse_paper_case(as_string(get(root, "case"), "case"));
    config.error_target = parse_error_target(
        as_string(get(root, "error_target"), "error_target"));
    config.evidence_policy = parse_evidence_policy(
        as_string(get(root, "evidence_policy"), "evidence_policy"));
    config.git_commit = as_string(get(root, "git_commit"), "git_commit");
    config.hlod_prior_corrector_space_id = as_string(
        get(root, "hlod_prior_corrector_space_id"),
        "hlod_prior_corrector_space_id");
    config.hlod_prior_oversampling = as_integer(
        get(root, "hlod_prior_oversampling"), "hlod_prior_oversampling");
    config.method_id = parse_paper_method(as_string(get(root, "method"), "method"));
    config.mu0 = as_number(get(root, "mu0"), "mu0");
    config.q0 = as_number(get(root, "q0"), "q0");
    config.repeat_index = as_integer(get(root, "repeat_index"), "repeat_index");
    config.rho_aud = as_number(get(root, "rho_aud"), "rho_aud");
    config.schema_version = as_integer(get(root, "schema_version"), "schema_version");
    config.tau = as_number(get(root, "tau"), "tau");
    config.theta_H = as_number(get(root, "theta_H"), "theta_H");
    config.theta_h = as_number(get(root, "theta_h"), "theta_h");
    config.timing_repeats = as_integer(get(root, "timing_repeats"), "timing_repeats");
    config.tolerance = as_number(get(root, "tolerance"), "tolerance");
    config.wavenumber = as_number(get(root, "wavenumber"), "wavenumber");

    const JsonObject &backend = as_object(
        get(root, "numerical_backend"), "numerical_backend");
    require_keys(backend,
        {"audit_doerfler_theta", "audit_saturation_factor", "boundary_beta",
         "certificate", "certificate_constant_set_hash", "initial_coarse_level",
         "initial_fine_level", "initial_oversampling", "kernel_riesz_solver",
         "patch_solver", "petrov_mode"},
        "numerical_backend");
    config.numerical_backend.audit_doerfler_theta = as_number(
        get(backend, "audit_doerfler_theta"), "audit_doerfler_theta");
    config.numerical_backend.audit_saturation_factor = as_number(
        get(backend, "audit_saturation_factor"), "audit_saturation_factor");
    config.numerical_backend.boundary_beta = as_number(
        get(backend, "boundary_beta"), "boundary_beta");
    config.numerical_backend.certificate_constant_set_hash = as_string(
        get(backend, "certificate_constant_set_hash"), "certificate_constant_set_hash");
    config.numerical_backend.initial_coarse_level = as_integer(
        get(backend, "initial_coarse_level"), "initial_coarse_level");
    config.numerical_backend.initial_fine_level = as_integer(
        get(backend, "initial_fine_level"), "initial_fine_level");
    config.numerical_backend.initial_oversampling = as_integer(
        get(backend, "initial_oversampling"), "initial_oversampling");
    config.numerical_backend.kernel_riesz_solver = parse_kernel_solver(
        as_string(get(backend, "kernel_riesz_solver"), "kernel_riesz_solver"));
    config.numerical_backend.petrov_mode = parse_petrov_mode(
        as_string(get(backend, "petrov_mode"), "petrov_mode"));

    const JsonObject &certificate = as_object(
        get(backend, "certificate"), "numerical_backend.certificate");
    require_keys(certificate,
        {"cluster_absolute_gap", "cluster_relative_gap", "conjugation_tolerance",
         "precision_bits", "q0"},
        "numerical_backend.certificate");
    config.numerical_backend.certificate.cluster_absolute_gap = as_number(
        get(certificate, "cluster_absolute_gap"), "cluster_absolute_gap");
    config.numerical_backend.certificate.cluster_relative_gap = as_number(
        get(certificate, "cluster_relative_gap"), "cluster_relative_gap");
    config.numerical_backend.certificate.conjugation_tolerance = as_number(
        get(certificate, "conjugation_tolerance"), "conjugation_tolerance");
    config.numerical_backend.certificate.precision_bits = as_integer(
        get(certificate, "precision_bits"), "precision_bits");
    config.numerical_backend.certificate.q0 = as_number(
        get(certificate, "q0"), "numerical_backend.certificate.q0");

    const JsonObject &patch_solver = as_object(
        get(backend, "patch_solver"), "numerical_backend.patch_solver");
    const bool has_maximum_parallel_solves =
        patch_solver.contains("maximum_parallel_solves");
    JsonObject patch_solver_contract = patch_solver;
    patch_solver_contract.erase("maximum_parallel_solves");
    require_keys(patch_solver_contract,
        {"fallback_to_direct", "gmres", "kind", "reuse_identical_factorization",
         "shifted", "symbolic_cache_slots"},
        "numerical_backend.patch_solver");
    config.numerical_backend.patch_solver.fallback_to_direct = as_bool(
        get(patch_solver, "fallback_to_direct"), "fallback_to_direct");
    config.numerical_backend.patch_solver.kind = parse_patch_solver_kind(
        as_string(get(patch_solver, "kind"), "kind"));
    config.numerical_backend.patch_solver.reuse_identical_factorization = as_bool(
        get(patch_solver, "reuse_identical_factorization"),
        "reuse_identical_factorization");
    config.numerical_backend.patch_solver.symbolic_cache_slots = as_integer(
        get(patch_solver, "symbolic_cache_slots"), "symbolic_cache_slots");
    if (has_maximum_parallel_solves) {
        config.numerical_backend.patch_solver.maximum_parallel_solves = as_integer(
            get(patch_solver, "maximum_parallel_solves"),
            "maximum_parallel_solves");
    }

    const JsonObject &gmres = as_object(get(patch_solver, "gmres"), "gmres");
    require_keys(gmres,
        {"absolute_tolerance", "max_iterations", "relative_tolerance",
         "reorthogonalize", "restart"},
        "gmres");
    config.numerical_backend.patch_solver.gmres.absolute_tolerance = as_number(
        get(gmres, "absolute_tolerance"), "absolute_tolerance");
    config.numerical_backend.patch_solver.gmres.max_iterations = as_integer(
        get(gmres, "max_iterations"), "max_iterations");
    config.numerical_backend.patch_solver.gmres.relative_tolerance = as_number(
        get(gmres, "relative_tolerance"), "relative_tolerance");
    config.numerical_backend.patch_solver.gmres.reorthogonalize = as_bool(
        get(gmres, "reorthogonalize"), "reorthogonalize");
    config.numerical_backend.patch_solver.gmres.restart = as_integer(
        get(gmres, "restart"), "restart");

    const JsonObject &shifted = as_object(get(patch_solver, "shifted"), "shifted");
    require_keys(shifted,
        {"absolute_epsilon", "alpha", "coarse_max_dofs", "inverse",
         "jacobi_weight", "post_smooth", "pre_smooth", "rule"},
        "shifted");
    config.numerical_backend.patch_solver.shifted.absolute_epsilon = as_number(
        get(shifted, "absolute_epsilon"), "absolute_epsilon");
    config.numerical_backend.patch_solver.shifted.alpha = as_number(
        get(shifted, "alpha"), "alpha");
    config.numerical_backend.patch_solver.shifted.coarse_max_dofs = as_integer(
        get(shifted, "coarse_max_dofs"), "coarse_max_dofs");
    config.numerical_backend.patch_solver.shifted.inverse = parse_shifted_inverse(
        as_string(get(shifted, "inverse"), "inverse"));
    config.numerical_backend.patch_solver.shifted.jacobi_weight = as_number(
        get(shifted, "jacobi_weight"), "jacobi_weight");
    config.numerical_backend.patch_solver.shifted.post_smooth = as_integer(
        get(shifted, "post_smooth"), "post_smooth");
    config.numerical_backend.patch_solver.shifted.pre_smooth = as_integer(
        get(shifted, "pre_smooth"), "pre_smooth");
    config.numerical_backend.patch_solver.shifted.rule = parse_shift_rule(
        as_string(get(shifted, "rule"), "rule"));

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

    const JsonObject &work_limits = as_object(get(root, "work_limits"), "work_limits");
    require_keys(work_limits,
        {"max_audit_dofs", "max_audit_refinements", "max_backend_work_units",
         "max_coarse_dofs", "max_coarse_refinements", "max_corrector_refinements",
         "max_elapsed_seconds", "max_fine_dofs", "max_oversampling",
         "max_oversampling_increments", "max_peak_memory_bytes",
         "max_state_transitions"},
        "work_limits");
    config.work_limits.max_audit_dofs = as_uint64(
        get(work_limits, "max_audit_dofs"), "max_audit_dofs");
    config.work_limits.max_audit_refinements = as_uint64(
        get(work_limits, "max_audit_refinements"), "max_audit_refinements");
    config.work_limits.max_backend_work_units = as_uint64(
        get(work_limits, "max_backend_work_units"), "max_backend_work_units");
    config.work_limits.max_coarse_dofs = as_uint64(
        get(work_limits, "max_coarse_dofs"), "max_coarse_dofs");
    config.work_limits.max_coarse_refinements = as_uint64(
        get(work_limits, "max_coarse_refinements"), "max_coarse_refinements");
    config.work_limits.max_corrector_refinements = as_uint64(
        get(work_limits, "max_corrector_refinements"), "max_corrector_refinements");
    config.work_limits.max_elapsed_seconds = as_number(
        get(work_limits, "max_elapsed_seconds"), "max_elapsed_seconds");
    config.work_limits.max_fine_dofs = as_uint64(
        get(work_limits, "max_fine_dofs"), "max_fine_dofs");
    config.work_limits.max_oversampling = as_integer(
        get(work_limits, "max_oversampling"), "max_oversampling");
    config.work_limits.max_oversampling_increments = as_uint64(
        get(work_limits, "max_oversampling_increments"),
        "max_oversampling_increments");
    config.work_limits.max_peak_memory_bytes = as_uint64(
        get(work_limits, "max_peak_memory_bytes"), "max_peak_memory_bytes");
    config.work_limits.max_state_transitions = as_uint64(
        get(work_limits, "max_state_transitions"), "max_state_transitions");
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

std::string_view to_string(const PracticalPaperMethod id) {
    switch (id) {
    case PracticalPaperMethod::Palod: return "PALOD";
    case PracticalPaperMethod::HlodFixed: return "HLOD-fixed";
    case PracticalPaperMethod::Slod: return "SLOD";
    case PracticalPaperMethod::Ufem: return "UFEM";
    case PracticalPaperMethod::Afem: return "AFEM";
    }
    throw std::invalid_argument("unknown practical paper method enum");
}

PracticalPaperMethod parse_practical_paper_method(const std::string_view text) {
    for (const PracticalPaperMethod method : {
             PracticalPaperMethod::Palod,
             PracticalPaperMethod::HlodFixed,
             PracticalPaperMethod::Slod,
             PracticalPaperMethod::Ufem,
             PracticalPaperMethod::Afem}) {
        if (to_string(method) == text) return method;
    }
    throw std::invalid_argument(
        "unknown practical paper method: " + std::string(text));
}

std::string_view to_string(const PracticalTrajectoryPolicy policy) {
    switch (policy) {
    case PracticalTrajectoryPolicy::PracticalIndicator:
        return "practical_indicator";
    case PracticalTrajectoryPolicy::FixedWorkHorizon:
        return "fixed_work_horizon";
    }
    throw std::invalid_argument("unknown practical trajectory policy enum");
}

PracticalTrajectoryPolicy parse_practical_trajectory_policy(
    const std::string_view text) {
    for (const PracticalTrajectoryPolicy policy : {
             PracticalTrajectoryPolicy::PracticalIndicator,
             PracticalTrajectoryPolicy::FixedWorkHorizon}) {
        if (to_string(policy) == text) return policy;
    }
    throw std::invalid_argument(
        "unknown practical trajectory policy: " + std::string(text));
}

int standard_lod_prior_ell(const double wavenumber) {
    if (!std::isfinite(wavenumber) || !(wavenumber > 0.0)) {
        throw std::invalid_argument(
            "standard LOD prior requires a positive finite wavenumber");
    }
    constexpr double prior_coefficient = 1.0;
    return std::max(
        1, static_cast<int>(std::ceil(
               prior_coefficient * std::log2(wavenumber))));
}

void validate_practical_paper_config(const PracticalPaperConfig &config) {
    if (config.schema_version != practical_paper_schema_version) {
        throw std::invalid_argument("unsupported practical paper schema version");
    }
    (void)case_definition(config.case_id);
    (void)to_string(config.method_id);
    (void)to_string(config.trajectory_policy);
    if (!(config.wavenumber == 2.0 || config.wavenumber == 4.0 ||
          config.wavenumber == 8.0 || config.wavenumber == 16.0 ||
          config.wavenumber == 32.0)) {
        throw std::invalid_argument(
            "practical paper wavenumber must be 2, 4, 8, 16, or 32");
    }
    if (!std::isfinite(config.singular_oscillatory_fraction)
        || config.singular_oscillatory_fraction < 0.0
        || config.singular_oscillatory_fraction > 1.0
        || (config.case_id != PaperCase::S
            && config.singular_oscillatory_fraction != 1.0)) {
        throw std::invalid_argument(
            "singular_oscillatory_fraction must lie in [0,1] and may differ "
            "from one only for case S");
    }
    if (!std::isfinite(config.singular_cutoff_outer_radius)
        || !(config.singular_cutoff_outer_radius > 0.25)
        || config.singular_cutoff_outer_radius > 1.0
        || (config.case_id != PaperCase::S
            && config.singular_cutoff_outer_radius != 0.5)) {
        throw std::invalid_argument(
            "singular_cutoff_outer_radius must lie in (0.25,1] and may "
            "differ from 0.5 only for case S");
    }
    if (config.case_id != PaperCase::S && config.singular_quintic_cutoff) {
        throw std::invalid_argument(
            "singular_quintic_cutoff may be enabled only for case S");
    }
    if (!std::isfinite(config.smooth_wave_amplitude)
        || config.smooth_wave_amplitude < 0.0
        || config.smooth_wave_amplitude > 1.0
        || (config.case_id != PaperCase::S
            && config.smooth_wave_amplitude != 0.0)) {
        throw std::invalid_argument(
            "smooth_wave_amplitude must lie in [0,1] and may be nonzero "
            "only for case S");
    }
    if (config.singular_solution_profile != "radial-cutoff"
        && config.singular_solution_profile != "boundary-weight-gaussian") {
        throw std::invalid_argument("unknown singular_solution_profile");
    }
    if (config.case_id != PaperCase::S
        && config.singular_solution_profile != "radial-cutoff") {
        throw std::invalid_argument(
            "singular_solution_profile may differ from radial-cutoff only for case S");
    }
    if (config.singular_solution_profile == "boundary-weight-gaussian"
        && (config.singular_oscillatory_fraction != 0.0
            || config.singular_cutoff_outer_radius != 0.5
            || config.singular_quintic_cutoff)) {
        throw std::invalid_argument(
            "boundary-weight-gaussian requires a nonoscillatory corner and default legacy cutoff fields");
    }
    if (config.reference_mesh != "uniform-nvb" ||
        config.ambient_mesh != "reference-shadow") {
        throw std::invalid_argument(
            "v4 requires uniform-nvb reference_mesh and reference-shadow ambient_mesh");
    }
    if (config.initial_coarse_level < 0 ||
        config.reference_level <= config.initial_coarse_level ||
        config.ell0 < 0 || config.ell_max < config.ell0
        || config.reference_refresh_level_gap < 0
        || config.maximum_reference_level < 0
        || config.minimum_reference_level_gap < 0
        || config.minimum_reference_level_gap
            >= config.reference_level - config.initial_coarse_level) {
        throw std::invalid_argument("practical mesh levels or ell limits are invalid");
    }
    if (config.maximum_patch_threads < 0
        || static_cast<std::uint64_t>(config.maximum_patch_threads)
            > max_exact_json_integer) {
        throw std::invalid_argument(
            "maximum_patch_threads must be a JSON-safe nonnegative integer");
    }
    if (config.patch_symbolic_cache_slots <= 0
        || static_cast<std::uint64_t>(config.patch_symbolic_cache_slots)
            > max_exact_json_integer
        || config.slod_direct_schur_min_reference_dofs
            > max_exact_json_integer) {
        throw std::invalid_argument(
            "patch cache and SLOD solver-switch parameters must be JSON-safe");
    }
    if (config.slod_direct_schur_min_reference_dofs > 0
        && (config.method_id != PracticalPaperMethod::Slod
            || config.patch_solver_kind
                != HelmholtzPatchSolverKind::DirectSaddle)) {
        throw std::invalid_argument(
            "the automatic direct-Schur threshold requires SLOD with a direct-saddle base solver");
    }
    if (config.manufactured_exact_only_errors
        && config.method_id != PracticalPaperMethod::Slod) {
        throw std::invalid_argument(
            "manufactured_exact_only_errors is supported only for SLOD");
    }
    std::size_t previous_refresh = 0;
    for (const std::size_t refresh : config.reference_refresh_H_steps) {
        if (refresh == 0 || refresh > config.work_limits.maximum_H_steps
            || refresh <= previous_refresh) {
            throw std::invalid_argument(
                "reference_refresh_H_steps must be strictly increasing and lie in "
                "[1, maximum_H_steps]");
        }
        previous_refresh = refresh;
    }
    if (!config.reference_refresh_H_steps.empty()
        && (config.trajectory_policy
                != PracticalTrajectoryPolicy::FixedWorkHorizon
            || (config.method_id != PracticalPaperMethod::Palod
                && config.method_id != PracticalPaperMethod::HlodFixed))) {
        throw std::invalid_argument(
            "scheduled reference refresh is supported only for fixed-horizon PALOD/HLOD-fixed");
    }
    const bool has_level_gap_refresh =
        config.reference_refresh_level_gap > 0
        || config.maximum_reference_level > 0;
    if (has_level_gap_refresh
        && (config.reference_refresh_level_gap <= 0
            || config.maximum_reference_level <= config.reference_level
            || config.reference_refresh_level_gap
                >= config.reference_level - config.initial_coarse_level
            || config.trajectory_policy
                != PracticalTrajectoryPolicy::FixedWorkHorizon
            || config.method_id != PracticalPaperMethod::Palod
            || !config.reference_refresh_H_steps.empty()
            || config.minimum_reference_level_gap > 0)) {
        throw std::invalid_argument(
            "level-gap reference refresh requires fixed-horizon PALOD, an initial gap above the trigger, a larger maximum reference level, and no other refresh/gap policy");
    }
    if (config.minimum_reference_level_gap > 0
        && (config.trajectory_policy
                != PracticalTrajectoryPolicy::FixedWorkHorizon
            || !config.reference_refresh_H_steps.empty()
            || (config.method_id != PracticalPaperMethod::Palod
                && config.method_id != PracticalPaperMethod::Afem))) {
        throw std::invalid_argument(
            "minimum_reference_level_gap is supported only for single-epoch fixed-horizon PALOD/AFEM");
    }
    if (config.method_id == PracticalPaperMethod::HlodFixed
        && config.ell0 != config.ell_max) {
        throw std::invalid_argument(
            "HLOD-fixed requires one frozen ell value (ell0 == ell_max)");
    }
    if (config.method_id == PracticalPaperMethod::Slod) {
        if (config.ell0 < 1 || config.ell0 != config.ell_max) {
            throw std::invalid_argument(
                "SLOD requires one positive frozen empirical oversampling depth");
        }
    }
    if ((config.method_id == PracticalPaperMethod::Ufem
         || config.method_id == PracticalPaperMethod::Afem)
        && (config.ell0 != 0 || config.ell_max != 0)) {
        throw std::invalid_argument(
            "conforming FEM baselines require ell0 == ell_max == 0");
    }
    const std::array<double, 9> positive{
        config.boundary_beta, config.c_H, config.theta_loc,
        config.C0_usr, config.C1_usr, config.theta_H, config.rho_star,
        config.practical_stop_tolerance,
        config.plateau_diagnostic.minimum_geometric_mean_ratio};
    if (std::any_of(positive.begin(), positive.end(), [](const double value) {
            return !std::isfinite(value) || !(value > 0.0);
        }) || config.theta_H > 1.0 || config.rho_star > 1.0
        || config.plateau_diagnostic.minimum_geometric_mean_ratio > 1.0
        || !std::isfinite(
            config.plateau_diagnostic.maximum_relative_oscillation)
        || config.plateau_diagnostic.maximum_relative_oscillation < 0.0
        || config.plateau_diagnostic.window_steps == 0
        || config.reference_adequacy.refinement_levels != 1
        || !std::isfinite(
            config.reference_adequacy.maximum_terminal_error_fraction)
        || !(config.reference_adequacy.maximum_terminal_error_fraction > 0.0)) {
        throw std::invalid_argument("practical decision parameters are invalid");
    }
    (void)petrov_mode_name(config.petrov_mode);
    (void)patch_solver_kind_name(config.patch_solver_kind);
    (void)kernel_solver_name(config.kernel_riesz_solver);
    if (config.work_limits.maximum_iterations == 0 ||
        config.work_limits.maximum_unknowns == 0 ||
        config.work_limits.maximum_coarse_elements == 0 ||
        config.work_limits.maximum_ambient_elements == 0 ||
        !std::isfinite(config.work_limits.maximum_wall_seconds) ||
        config.work_limits.maximum_wall_seconds < 0.0) {
        throw std::invalid_argument("practical work limits are invalid");
    }
    const std::array<std::size_t, 6> integer_limits{
        config.work_limits.maximum_iterations,
        config.work_limits.maximum_H_steps,
        config.work_limits.maximum_unknowns,
        config.work_limits.maximum_coarse_elements,
        config.work_limits.maximum_ambient_elements,
        config.plateau_diagnostic.window_steps};
    if (std::any_of(integer_limits.begin(), integer_limits.end(), [](const std::size_t value) {
            return value > max_exact_json_integer;
        })) {
        throw std::invalid_argument("practical work limits must be JSON-safe integers");
    }
    constexpr std::array<double, 4> targets{{0.1, 0.05, 0.02, 0.01}};
    if (config.relative_energy_targets != targets) {
        throw std::invalid_argument("relative energy targets differ from the frozen v2 protocol");
    }
    if (!((config.timing_repeats == 1 || config.timing_repeats == 3) &&
          config.repeat_index >= 0 && config.repeat_index < config.timing_repeats)) {
        throw std::invalid_argument("practical timing repeat policy is invalid");
    }
    if (config.quadrature.base_triangle_order <= 0 ||
        config.quadrature.gaussian_triangle_order <= 0 ||
        config.quadrature.singular_triangle_order <= 0 ||
        config.quadrature.max_recursive_subdivisions < 0) {
        throw std::invalid_argument("practical quadrature policy is invalid");
    }
    const std::array<double, 4> tolerances{
        config.tolerances.linear_relative_residual,
        config.tolerances.eigen_relative_residual,
        config.tolerances.interpolation_right_inverse,
        config.tolerances.prolongation_composition};
    if (std::any_of(tolerances.begin(), tolerances.end(), [](const double value) {
            return !std::isfinite(value) || !(value > 0.0);
        })) {
        throw std::invalid_argument("practical tolerances must be positive and finite");
    }
    const bool manuscript_hash_valid = config.manuscript_sha256.size() == 64 &&
        std::all_of(config.manuscript_sha256.begin(), config.manuscript_sha256.end(),
                    [](const char character) {
                        return (character >= '0' && character <= '9') ||
                               (character >= 'a' && character <= 'f');
                    });
    if (config.git_commit.empty() || config.build_hash.empty() ||
        !manuscript_hash_valid) {
        throw std::invalid_argument("practical provenance is incomplete or malformed");
    }
}

adaptive::PracticalDriverConfig make_practical_driver_config(
    const PracticalPaperConfig &config) {
    validate_practical_paper_config(config);
    if (config.method_id != PracticalPaperMethod::Palod
        && config.method_id != PracticalPaperMethod::HlodFixed) {
        throw std::invalid_argument(
            "the practical LOD driver executes only PALOD and HLOD-fixed");
    }
    adaptive::PracticalDriverConfig result;
    result.initial_coarse_level = config.initial_coarse_level;
    result.reference_level = config.reference_level;
    result.reference_epoch = config.reference_epoch;
    result.reference_refresh_H_steps = config.reference_refresh_H_steps;
    result.reference_refresh_level_gap =
        config.reference_refresh_level_gap;
    result.maximum_reference_level = config.maximum_reference_level;
    result.minimum_reference_level_gap =
        config.minimum_reference_level_gap;
    result.ell0 = config.ell0;
    result.ell_max = config.ell_max;
    result.wavenumber = config.wavenumber;
    result.boundary_beta = config.boundary_beta;
    result.c_H = config.c_H;
    result.theta_loc = config.theta_loc;
    result.C0_usr = config.C0_usr;
    result.C1_usr = config.C1_usr;
    result.theta_H = config.theta_H;
    result.rho_star = config.rho_star;
    result.tolerance_reference = config.practical_stop_tolerance;
    result.stop_policy = config.trajectory_policy
            == PracticalTrajectoryPolicy::FixedWorkHorizon
        ? adaptive::PracticalStopPolicy::FixedWorkHorizon
        : adaptive::PracticalStopPolicy::IndicatorTolerance;
    result.localization_policy = config.method_id
            == PracticalPaperMethod::HlodFixed
        ? adaptive::PracticalLocalizationPolicy::FixedGlobalEll
        : adaptive::PracticalLocalizationPolicy::AdaptiveGlobalEll;
    result.mode = config.petrov_mode;
    result.patch_solver.kind = config.patch_solver_kind;
    result.patch_solver.symbolic_cache_slots =
        config.patch_symbolic_cache_slots;
    result.patch_solver.reuse_identical_factorization =
        config.patch_reuse_identical_factorization;
    result.patch_solver.maximum_parallel_solves =
        config.maximum_patch_threads;
    result.patch_solver.gmres.relative_tolerance =
        config.tolerances.linear_relative_residual;
    result.riesz_solver = config.kernel_riesz_solver;
    result.localization_eigen.relative_tolerance =
        config.tolerances.eigen_relative_residual;
    result.limits = config.work_limits;
    return result;
}

std::string canonical_json(const PracticalPaperConfig &config) {
    validate_practical_paper_config(config);
    std::ostringstream out;
    out << "{\"C0_usr\":" << number(config.C0_usr)
        << ",\"C1_usr\":" << number(config.C1_usr)
        << ",\"ambient_mesh\":" << json_string(config.ambient_mesh)
        << ",\"boundary_beta\":" << number(config.boundary_beta)
        << ",\"build_hash\":" << json_string(config.build_hash)
        << ",\"c_H\":" << number(config.c_H)
        << ",\"case\":" << json_string(to_string(config.case_id))
        << ",\"ell0\":" << config.ell0
        << ",\"ell_max\":" << config.ell_max
        << ",\"git_commit\":" << json_string(config.git_commit)
        << ",\"initial_coarse_level\":" << config.initial_coarse_level
        << ",\"kernel_riesz_solver\":"
        << json_string(kernel_solver_name(config.kernel_riesz_solver))
        << ",\"manuscript_sha256\":" << json_string(config.manuscript_sha256)
        << ",\"maximum_patch_threads\":" << config.maximum_patch_threads
        << ",\"maximum_reference_level\":"
        << config.maximum_reference_level
        << ",\"minimum_reference_level_gap\":"
        << config.minimum_reference_level_gap
        << ",\"method\":" << json_string(to_string(config.method_id))
        << ",\"patch_solver_kind\":"
        << json_string(patch_solver_kind_name(config.patch_solver_kind))
        << ",\"patch_symbolic_cache_slots\":"
        << config.patch_symbolic_cache_slots
        << ",\"patch_reuse_identical_factorization\":"
        << (config.patch_reuse_identical_factorization ? "true" : "false")
        << ",\"petrov_mode\":" << json_string(petrov_mode_name(config.petrov_mode))
        << ",\"plateau_diagnostic\":{\"maximum_relative_oscillation\":"
        << number(config.plateau_diagnostic.maximum_relative_oscillation)
        << ",\"minimum_geometric_mean_ratio\":"
        << number(config.plateau_diagnostic.minimum_geometric_mean_ratio)
        << ",\"window_steps\":"
        << config.plateau_diagnostic.window_steps << "}"
        << ",\"practical_stop_tolerance\":"
        << number(config.practical_stop_tolerance)
        << ",\"quadrature\":{\"base_triangle_order\":"
        << config.quadrature.base_triangle_order
        << ",\"gaussian_triangle_order\":"
        << config.quadrature.gaussian_triangle_order
        << ",\"max_recursive_subdivisions\":"
        << config.quadrature.max_recursive_subdivisions
        << ",\"singular_triangle_order\":"
        << config.quadrature.singular_triangle_order << "}"
        << ",\"reference_epoch\":" << config.reference_epoch;
    if (!config.reference_refresh_H_steps.empty()) {
        out << ",\"reference_refresh_H_steps\":[";
        for (std::size_t index = 0;
             index < config.reference_refresh_H_steps.size(); ++index) {
            if (index) out << ',';
            out << config.reference_refresh_H_steps[index];
        }
        out << "]";
    }
    out << ",\"reference_adequacy\":{\"enabled\":"
        << (config.reference_adequacy.enabled ? "true" : "false")
        << ",\"maximum_terminal_error_fraction\":"
        << number(config.reference_adequacy.maximum_terminal_error_fraction)
        << ",\"refinement_levels\":"
        << config.reference_adequacy.refinement_levels << "}"
        << ",\"reference_level\":" << config.reference_level
        << ",\"reference_refresh_level_gap\":"
        << config.reference_refresh_level_gap
        << ",\"reference_mesh\":" << json_string(config.reference_mesh)
        << ",\"relative_energy_targets\":[";
    for (std::size_t index = 0; index < config.relative_energy_targets.size(); ++index) {
        if (index) out << ',';
        out << number(config.relative_energy_targets[index]);
    }
    out << "]"
        << ",\"repeat_index\":" << config.repeat_index
        << ",\"rho_star\":" << number(config.rho_star)
        << ",\"schema_version\":" << config.schema_version;
    if (config.singular_oscillatory_fraction != 1.0) {
        out << ",\"singular_oscillatory_fraction\":"
            << number(config.singular_oscillatory_fraction);
    }
    if (config.singular_cutoff_outer_radius != 0.5) {
        out << ",\"singular_cutoff_outer_radius\":"
            << number(config.singular_cutoff_outer_radius);
    }
    if (config.singular_quintic_cutoff) {
        out << ",\"singular_quintic_cutoff\":true";
    }
    if (config.smooth_wave_amplitude != 0.0) {
        out << ",\"smooth_wave_amplitude\":"
            << number(config.smooth_wave_amplitude);
    }
    if (config.singular_solution_profile != "radial-cutoff") {
        out << ",\"singular_solution_profile\":"
            << json_string(config.singular_solution_profile);
    }
    if (config.slod_direct_schur_min_reference_dofs > 0) {
        out << ",\"slod_direct_schur_min_reference_dofs\":"
            << config.slod_direct_schur_min_reference_dofs;
    }
    if (config.manufactured_exact_only_errors) {
        out << ",\"manufactured_exact_only_errors\":true";
    }
    out << ",\"theta_H\":" << number(config.theta_H)
        << ",\"theta_loc\":" << number(config.theta_loc)
        << ",\"timing_repeats\":" << config.timing_repeats
        << ",\"trajectory_policy\":"
        << json_string(to_string(config.trajectory_policy))
        << ",\"tolerances\":{\"eigen_relative_residual\":"
        << number(config.tolerances.eigen_relative_residual)
        << ",\"interpolation_right_inverse\":"
        << number(config.tolerances.interpolation_right_inverse)
        << ",\"linear_relative_residual\":"
        << number(config.tolerances.linear_relative_residual)
        << ",\"prolongation_composition\":"
        << number(config.tolerances.prolongation_composition) << "}"
        << ",\"wavenumber\":" << number(config.wavenumber)
        << ",\"work_limits\":{\"maximum_H_steps\":"
        << config.work_limits.maximum_H_steps
        << ",\"maximum_ambient_elements\":"
        << config.work_limits.maximum_ambient_elements
        << ",\"maximum_coarse_elements\":"
        << config.work_limits.maximum_coarse_elements
        << ",\"maximum_iterations\":"
        << config.work_limits.maximum_iterations
        << ",\"maximum_unknowns\":"
        << config.work_limits.maximum_unknowns
        << ",\"maximum_wall_seconds\":"
        << number(config.work_limits.maximum_wall_seconds) << "}}";
    return out.str();
}

PracticalPaperConfig parse_practical_paper_config(const std::string_view json) {
    const JsonValue parsed_root = JsonParser(json).parse();
    const JsonObject &root = as_object(parsed_root, "root");
    const bool has_minimum_reference_level_gap =
        root.contains("minimum_reference_level_gap");
    const bool has_reference_refresh_level_gap =
        root.contains("reference_refresh_level_gap");
    const bool has_maximum_reference_level =
        root.contains("maximum_reference_level");
    const bool has_maximum_patch_threads =
        root.contains("maximum_patch_threads");
    const bool has_patch_symbolic_cache_slots =
        root.contains("patch_symbolic_cache_slots");
    const bool has_patch_reuse_identical_factorization =
        root.contains("patch_reuse_identical_factorization");
    const bool has_slod_direct_schur_min_reference_dofs =
        root.contains("slod_direct_schur_min_reference_dofs");
    const bool has_manufactured_exact_only_errors =
        root.contains("manufactured_exact_only_errors");
    const bool has_singular_oscillatory_fraction =
        root.contains("singular_oscillatory_fraction");
    const bool has_singular_cutoff_outer_radius =
        root.contains("singular_cutoff_outer_radius");
    const bool has_singular_quintic_cutoff =
        root.contains("singular_quintic_cutoff");
    const bool has_smooth_wave_amplitude =
        root.contains("smooth_wave_amplitude");
    const bool has_singular_solution_profile =
        root.contains("singular_solution_profile");
    JsonObject contract_root = root;
    contract_root.erase("minimum_reference_level_gap");
    contract_root.erase("reference_refresh_level_gap");
    contract_root.erase("maximum_reference_level");
    contract_root.erase("maximum_patch_threads");
    contract_root.erase("patch_symbolic_cache_slots");
    contract_root.erase("patch_reuse_identical_factorization");
    contract_root.erase("slod_direct_schur_min_reference_dofs");
    contract_root.erase("manufactured_exact_only_errors");
    contract_root.erase("singular_oscillatory_fraction");
    contract_root.erase("singular_cutoff_outer_radius");
    contract_root.erase("singular_quintic_cutoff");
    contract_root.erase("smooth_wave_amplitude");
    contract_root.erase("singular_solution_profile");
    const bool has_refresh_schedule =
        root.contains("reference_refresh_H_steps");
    if (has_refresh_schedule) {
        require_keys(contract_root,
            {"C0_usr", "C1_usr", "ambient_mesh", "boundary_beta", "build_hash",
             "c_H", "case", "ell0", "ell_max", "git_commit", "initial_coarse_level",
             "kernel_riesz_solver", "manuscript_sha256", "method", "patch_solver_kind",
             "petrov_mode", "plateau_diagnostic", "practical_stop_tolerance",
             "quadrature", "reference_adequacy", "reference_epoch",
             "reference_refresh_H_steps", "reference_level", "reference_mesh",
             "relative_energy_targets", "repeat_index", "rho_star", "schema_version",
             "theta_H", "theta_loc", "timing_repeats", "tolerances",
             "trajectory_policy", "wavenumber", "work_limits"}, "root");
    } else {
        require_keys(contract_root,
            {"C0_usr", "C1_usr", "ambient_mesh", "boundary_beta", "build_hash",
             "c_H", "case", "ell0", "ell_max", "git_commit", "initial_coarse_level",
             "kernel_riesz_solver", "manuscript_sha256", "method", "patch_solver_kind",
             "petrov_mode", "plateau_diagnostic", "practical_stop_tolerance",
             "quadrature", "reference_adequacy", "reference_epoch", "reference_level",
             "reference_mesh", "relative_energy_targets", "repeat_index", "rho_star",
             "schema_version", "theta_H", "theta_loc", "timing_repeats", "tolerances",
             "trajectory_policy", "wavenumber", "work_limits"}, "root");
    }
    PracticalPaperConfig config;
    config.schema_version = as_integer(get(root, "schema_version"), "schema_version");
    config.case_id = parse_paper_case(as_string(get(root, "case"), "case"));
    config.method_id = parse_practical_paper_method(
        as_string(get(root, "method"), "method"));
    config.wavenumber = as_number(get(root, "wavenumber"), "wavenumber");
    if (has_singular_oscillatory_fraction) {
        config.singular_oscillatory_fraction = as_number(
            get(root, "singular_oscillatory_fraction"),
            "singular_oscillatory_fraction");
    }
    if (has_singular_cutoff_outer_radius) {
        config.singular_cutoff_outer_radius = as_number(
            get(root, "singular_cutoff_outer_radius"),
            "singular_cutoff_outer_radius");
    }
    if (has_singular_quintic_cutoff) {
        config.singular_quintic_cutoff = as_bool(
            get(root, "singular_quintic_cutoff"),
            "singular_quintic_cutoff");
    }
    if (has_smooth_wave_amplitude) {
        config.smooth_wave_amplitude = as_number(
            get(root, "smooth_wave_amplitude"),
            "smooth_wave_amplitude");
    }
    if (has_singular_solution_profile) {
        config.singular_solution_profile = as_string(
            get(root, "singular_solution_profile"),
            "singular_solution_profile");
    }
    config.reference_mesh = as_string(get(root, "reference_mesh"), "reference_mesh");
    config.reference_level = as_integer(get(root, "reference_level"), "reference_level");
    config.ambient_mesh = as_string(get(root, "ambient_mesh"), "ambient_mesh");
    config.reference_epoch = as_uint64(get(root, "reference_epoch"), "reference_epoch");
    if (has_minimum_reference_level_gap) {
        config.minimum_reference_level_gap = as_integer(
            get(root, "minimum_reference_level_gap"),
            "minimum_reference_level_gap");
    }
    if (has_reference_refresh_level_gap) {
        config.reference_refresh_level_gap = as_integer(
            get(root, "reference_refresh_level_gap"),
            "reference_refresh_level_gap");
    }
    if (has_maximum_reference_level) {
        config.maximum_reference_level = as_integer(
            get(root, "maximum_reference_level"),
            "maximum_reference_level");
    }
    if (has_maximum_patch_threads) {
        config.maximum_patch_threads = as_integer(
            get(root, "maximum_patch_threads"), "maximum_patch_threads");
    }
    if (has_patch_symbolic_cache_slots) {
        config.patch_symbolic_cache_slots = as_integer(
            get(root, "patch_symbolic_cache_slots"),
            "patch_symbolic_cache_slots");
    }
    if (has_patch_reuse_identical_factorization) {
        config.patch_reuse_identical_factorization = as_bool(
            get(root, "patch_reuse_identical_factorization"),
            "patch_reuse_identical_factorization");
    }
    if (has_slod_direct_schur_min_reference_dofs) {
        config.slod_direct_schur_min_reference_dofs =
            static_cast<std::size_t>(as_uint64(
                get(root, "slod_direct_schur_min_reference_dofs"),
                "slod_direct_schur_min_reference_dofs"));
    }
    if (has_manufactured_exact_only_errors) {
        config.manufactured_exact_only_errors = as_bool(
            get(root, "manufactured_exact_only_errors"),
            "manufactured_exact_only_errors");
    }
    if (has_refresh_schedule) {
        const JsonArray &refreshes = as_array(
            get(root, "reference_refresh_H_steps"),
            "reference_refresh_H_steps");
        config.reference_refresh_H_steps.reserve(refreshes.size());
        for (const JsonValue &refresh : refreshes) {
            config.reference_refresh_H_steps.push_back(
                static_cast<std::size_t>(as_uint64(
                    refresh, "reference_refresh_H_steps")));
        }
    }
    config.initial_coarse_level = as_integer(
        get(root, "initial_coarse_level"), "initial_coarse_level");
    config.ell0 = as_integer(get(root, "ell0"), "ell0");
    config.ell_max = as_integer(get(root, "ell_max"), "ell_max");
    config.boundary_beta = as_number(get(root, "boundary_beta"), "boundary_beta");
    config.c_H = as_number(get(root, "c_H"), "c_H");
    config.theta_loc = as_number(get(root, "theta_loc"), "theta_loc");
    config.C0_usr = as_number(get(root, "C0_usr"), "C0_usr");
    config.C1_usr = as_number(get(root, "C1_usr"), "C1_usr");
    config.theta_H = as_number(get(root, "theta_H"), "theta_H");
    config.rho_star = as_number(get(root, "rho_star"), "rho_star");
    config.trajectory_policy = parse_practical_trajectory_policy(
        as_string(get(root, "trajectory_policy"), "trajectory_policy"));
    config.practical_stop_tolerance = as_number(
        get(root, "practical_stop_tolerance"), "practical_stop_tolerance");
    config.petrov_mode = parse_petrov_mode(
        as_string(get(root, "petrov_mode"), "petrov_mode"));
    config.patch_solver_kind = parse_patch_solver_kind(
        as_string(get(root, "patch_solver_kind"), "patch_solver_kind"));
    config.kernel_riesz_solver = parse_kernel_solver(
        as_string(get(root, "kernel_riesz_solver"), "kernel_riesz_solver"));
    config.timing_repeats = as_integer(get(root, "timing_repeats"), "timing_repeats");
    config.repeat_index = as_integer(get(root, "repeat_index"), "repeat_index");
    config.git_commit = as_string(get(root, "git_commit"), "git_commit");
    config.build_hash = as_string(get(root, "build_hash"), "build_hash");
    config.manuscript_sha256 = as_string(
        get(root, "manuscript_sha256"), "manuscript_sha256");

    const JsonObject &reference_adequacy = as_object(
        get(root, "reference_adequacy"), "reference_adequacy");
    require_keys(reference_adequacy,
        {"enabled", "maximum_terminal_error_fraction", "refinement_levels"},
        "reference_adequacy");
    config.reference_adequacy.enabled = as_bool(
        get(reference_adequacy, "enabled"), "enabled");
    config.reference_adequacy.maximum_terminal_error_fraction = as_number(
        get(reference_adequacy, "maximum_terminal_error_fraction"),
        "maximum_terminal_error_fraction");
    config.reference_adequacy.refinement_levels = as_integer(
        get(reference_adequacy, "refinement_levels"), "refinement_levels");

    const JsonObject &plateau = as_object(
        get(root, "plateau_diagnostic"), "plateau_diagnostic");
    require_keys(plateau,
        {"maximum_relative_oscillation", "minimum_geometric_mean_ratio",
         "window_steps"},
        "plateau_diagnostic");
    config.plateau_diagnostic.window_steps =
        static_cast<std::size_t>(as_uint64(
            get(plateau, "window_steps"), "window_steps"));
    config.plateau_diagnostic.minimum_geometric_mean_ratio = as_number(
        get(plateau, "minimum_geometric_mean_ratio"),
        "minimum_geometric_mean_ratio");
    config.plateau_diagnostic.maximum_relative_oscillation = as_number(
        get(plateau, "maximum_relative_oscillation"),
        "maximum_relative_oscillation");

    const JsonArray &targets = as_array(
        get(root, "relative_energy_targets"), "relative_energy_targets");
    if (targets.size() != config.relative_energy_targets.size()) {
        throw std::invalid_argument("relative_energy_targets must have four entries");
    }
    for (std::size_t index = 0; index < targets.size(); ++index) {
        config.relative_energy_targets[index] =
            as_number(targets[index], "relative_energy_targets");
    }
    const JsonObject &quadrature = as_object(get(root, "quadrature"), "quadrature");
    require_keys(quadrature,
        {"base_triangle_order", "gaussian_triangle_order",
         "max_recursive_subdivisions", "singular_triangle_order"}, "quadrature");
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
        {"eigen_relative_residual", "interpolation_right_inverse",
         "linear_relative_residual", "prolongation_composition"}, "tolerances");
    config.tolerances.eigen_relative_residual = as_number(
        get(tolerances, "eigen_relative_residual"), "eigen_relative_residual");
    config.tolerances.interpolation_right_inverse = as_number(
        get(tolerances, "interpolation_right_inverse"), "interpolation_right_inverse");
    config.tolerances.linear_relative_residual = as_number(
        get(tolerances, "linear_relative_residual"), "linear_relative_residual");
    config.tolerances.prolongation_composition = as_number(
        get(tolerances, "prolongation_composition"), "prolongation_composition");
    const JsonObject &limits = as_object(get(root, "work_limits"), "work_limits");
    require_keys(limits,
        {"maximum_H_steps", "maximum_ambient_elements", "maximum_coarse_elements",
         "maximum_iterations", "maximum_unknowns", "maximum_wall_seconds"},
        "work_limits");
    config.work_limits.maximum_H_steps = static_cast<std::size_t>(as_uint64(
        get(limits, "maximum_H_steps"), "maximum_H_steps"));
    config.work_limits.maximum_ambient_elements = static_cast<std::size_t>(as_uint64(
        get(limits, "maximum_ambient_elements"), "maximum_ambient_elements"));
    config.work_limits.maximum_coarse_elements = static_cast<std::size_t>(as_uint64(
        get(limits, "maximum_coarse_elements"), "maximum_coarse_elements"));
    config.work_limits.maximum_iterations = static_cast<std::size_t>(as_uint64(
        get(limits, "maximum_iterations"), "maximum_iterations"));
    config.work_limits.maximum_unknowns = static_cast<std::size_t>(as_uint64(
        get(limits, "maximum_unknowns"), "maximum_unknowns"));
    config.work_limits.maximum_wall_seconds = as_number(
        get(limits, "maximum_wall_seconds"), "maximum_wall_seconds");
    validate_practical_paper_config(config);
    return config;
}

std::string canonical_config_hash(const PracticalPaperConfig &config) {
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

std::string make_run_id(const PracticalPaperConfig &config) {
    validate_practical_paper_config(config);
    std::ostringstream out;
    out << to_string(config.case_id) << '_'
        << safe_component(to_string(config.method_id))
        << "_k" << static_cast<int>(config.wavenumber)
        << "_r" << config.repeat_index << '_'
        << canonical_config_hash(config);
    return out.str();
}

bool operator==(const PracticalPaperConfig &lhs,
                const PracticalPaperConfig &rhs) {
    return lhs.schema_version == rhs.schema_version && lhs.case_id == rhs.case_id &&
        lhs.method_id == rhs.method_id && lhs.wavenumber == rhs.wavenumber &&
        lhs.singular_oscillatory_fraction ==
            rhs.singular_oscillatory_fraction &&
        lhs.singular_cutoff_outer_radius ==
            rhs.singular_cutoff_outer_radius &&
        lhs.singular_quintic_cutoff == rhs.singular_quintic_cutoff &&
        lhs.smooth_wave_amplitude == rhs.smooth_wave_amplitude &&
        lhs.singular_solution_profile == rhs.singular_solution_profile &&
        lhs.reference_mesh == rhs.reference_mesh &&
        lhs.reference_level == rhs.reference_level &&
        lhs.ambient_mesh == rhs.ambient_mesh &&
        lhs.reference_epoch == rhs.reference_epoch &&
        lhs.reference_refresh_H_steps == rhs.reference_refresh_H_steps &&
        lhs.reference_refresh_level_gap ==
            rhs.reference_refresh_level_gap &&
        lhs.maximum_reference_level == rhs.maximum_reference_level &&
        lhs.minimum_reference_level_gap ==
            rhs.minimum_reference_level_gap &&
        lhs.initial_coarse_level == rhs.initial_coarse_level &&
        lhs.ell0 == rhs.ell0 && lhs.ell_max == rhs.ell_max &&
        lhs.boundary_beta == rhs.boundary_beta && lhs.c_H == rhs.c_H &&
        lhs.theta_loc == rhs.theta_loc && lhs.C0_usr == rhs.C0_usr &&
        lhs.C1_usr == rhs.C1_usr && lhs.theta_H == rhs.theta_H &&
        lhs.rho_star == rhs.rho_star &&
        lhs.trajectory_policy == rhs.trajectory_policy &&
        lhs.practical_stop_tolerance == rhs.practical_stop_tolerance &&
        lhs.plateau_diagnostic.minimum_geometric_mean_ratio ==
            rhs.plateau_diagnostic.minimum_geometric_mean_ratio &&
        lhs.plateau_diagnostic.maximum_relative_oscillation ==
            rhs.plateau_diagnostic.maximum_relative_oscillation &&
        lhs.plateau_diagnostic.window_steps ==
            rhs.plateau_diagnostic.window_steps &&
        lhs.reference_adequacy.enabled == rhs.reference_adequacy.enabled &&
        lhs.reference_adequacy.refinement_levels ==
            rhs.reference_adequacy.refinement_levels &&
        lhs.reference_adequacy.maximum_terminal_error_fraction ==
            rhs.reference_adequacy.maximum_terminal_error_fraction &&
        lhs.petrov_mode == rhs.petrov_mode &&
        lhs.patch_solver_kind == rhs.patch_solver_kind &&
        lhs.patch_symbolic_cache_slots == rhs.patch_symbolic_cache_slots &&
        lhs.patch_reuse_identical_factorization ==
            rhs.patch_reuse_identical_factorization &&
        lhs.maximum_patch_threads == rhs.maximum_patch_threads &&
        lhs.slod_direct_schur_min_reference_dofs ==
            rhs.slod_direct_schur_min_reference_dofs &&
        lhs.manufactured_exact_only_errors ==
            rhs.manufactured_exact_only_errors &&
        lhs.kernel_riesz_solver == rhs.kernel_riesz_solver &&
        lhs.work_limits.maximum_iterations == rhs.work_limits.maximum_iterations &&
        lhs.work_limits.maximum_H_steps == rhs.work_limits.maximum_H_steps &&
        lhs.work_limits.maximum_unknowns == rhs.work_limits.maximum_unknowns &&
        lhs.work_limits.maximum_coarse_elements == rhs.work_limits.maximum_coarse_elements &&
        lhs.work_limits.maximum_ambient_elements == rhs.work_limits.maximum_ambient_elements &&
        lhs.work_limits.maximum_wall_seconds == rhs.work_limits.maximum_wall_seconds &&
        lhs.timing_repeats == rhs.timing_repeats &&
        lhs.repeat_index == rhs.repeat_index &&
        lhs.relative_energy_targets == rhs.relative_energy_targets &&
        lhs.quadrature == rhs.quadrature && lhs.tolerances == rhs.tolerances &&
        lhs.git_commit == rhs.git_commit && lhs.build_hash == rhs.build_hash &&
        lhs.manuscript_sha256 == rhs.manuscript_sha256;
}

bool operator==(const TolerancePolicy &lhs, const TolerancePolicy &rhs) {
    return lhs.linear_relative_residual == rhs.linear_relative_residual
        && lhs.eigen_relative_residual == rhs.eigen_relative_residual
        && lhs.interpolation_right_inverse == rhs.interpolation_right_inverse
        && lhs.prolongation_composition == rhs.prolongation_composition;
}

bool operator==(const NumericalBackendPolicy &lhs,
                const NumericalBackendPolicy &rhs) {
    return lhs.initial_coarse_level == rhs.initial_coarse_level
        && lhs.initial_fine_level == rhs.initial_fine_level
        && lhs.initial_oversampling == rhs.initial_oversampling
        && lhs.boundary_beta == rhs.boundary_beta
        && lhs.petrov_mode == rhs.petrov_mode
        && lhs.patch_solver.kind == rhs.patch_solver.kind
        && lhs.patch_solver.symbolic_cache_slots
            == rhs.patch_solver.symbolic_cache_slots
        && lhs.patch_solver.maximum_parallel_solves
            == rhs.patch_solver.maximum_parallel_solves
        && lhs.patch_solver.reuse_identical_factorization
            == rhs.patch_solver.reuse_identical_factorization
        && lhs.patch_solver.gmres.restart == rhs.patch_solver.gmres.restart
        && lhs.patch_solver.gmres.max_iterations
            == rhs.patch_solver.gmres.max_iterations
        && lhs.patch_solver.gmres.relative_tolerance
            == rhs.patch_solver.gmres.relative_tolerance
        && lhs.patch_solver.gmres.absolute_tolerance
            == rhs.patch_solver.gmres.absolute_tolerance
        && lhs.patch_solver.gmres.reorthogonalize
            == rhs.patch_solver.gmres.reorthogonalize
        && lhs.patch_solver.shifted.rule == rhs.patch_solver.shifted.rule
        && lhs.patch_solver.shifted.alpha == rhs.patch_solver.shifted.alpha
        && lhs.patch_solver.shifted.absolute_epsilon
            == rhs.patch_solver.shifted.absolute_epsilon
        && lhs.patch_solver.shifted.inverse == rhs.patch_solver.shifted.inverse
        && lhs.patch_solver.shifted.pre_smooth == rhs.patch_solver.shifted.pre_smooth
        && lhs.patch_solver.shifted.post_smooth == rhs.patch_solver.shifted.post_smooth
        && lhs.patch_solver.shifted.coarse_max_dofs
            == rhs.patch_solver.shifted.coarse_max_dofs
        && lhs.patch_solver.shifted.jacobi_weight
            == rhs.patch_solver.shifted.jacobi_weight
        && lhs.patch_solver.fallback_to_direct == rhs.patch_solver.fallback_to_direct
        && lhs.certificate.precision_bits == rhs.certificate.precision_bits
        && lhs.certificate.cluster_relative_gap
            == rhs.certificate.cluster_relative_gap
        && lhs.certificate.cluster_absolute_gap
            == rhs.certificate.cluster_absolute_gap
        && lhs.certificate.conjugation_tolerance
            == rhs.certificate.conjugation_tolerance
        && lhs.certificate.q0 == rhs.certificate.q0
        && lhs.kernel_riesz_solver == rhs.kernel_riesz_solver
        && lhs.audit_doerfler_theta == rhs.audit_doerfler_theta
        && lhs.audit_saturation_factor == rhs.audit_saturation_factor
        && lhs.certificate_constant_set_hash == rhs.certificate_constant_set_hash;
}

bool operator==(const PaperConfig &lhs, const PaperConfig &rhs) {
    return lhs.schema_version == rhs.schema_version && lhs.case_id == rhs.case_id
        && lhs.method_id == rhs.method_id && lhs.wavenumber == rhs.wavenumber
        && lhs.theta_H == rhs.theta_H && lhs.error_target == rhs.error_target
        && lhs.evidence_policy == rhs.evidence_policy && lhs.mu0 == rhs.mu0
        && lhs.q0 == rhs.q0 && lhs.tau == rhs.tau
        && lhs.theta_h == rhs.theta_h && lhs.rho_aud == rhs.rho_aud
        && lhs.tolerance == rhs.tolerance
        && lhs.numerical_backend == rhs.numerical_backend
        && equal_work_limits(lhs.work_limits, rhs.work_limits)
        && lhs.hlod_prior_corrector_space_id == rhs.hlod_prior_corrector_space_id
        && lhs.hlod_prior_oversampling == rhs.hlod_prior_oversampling
        && lhs.timing_repeats == rhs.timing_repeats
        && lhs.repeat_index == rhs.repeat_index
        && lhs.relative_energy_targets == rhs.relative_energy_targets
        && lhs.quadrature == rhs.quadrature && lhs.tolerances == rhs.tolerances
        && lhs.git_commit == rhs.git_commit && lhs.build_hash == rhs.build_hash;
}

HelmholtzFemSolverKind parse_reference_solver_kind(
    const std::string_view value) {
    if (value == "sparse_lu") return HelmholtzFemSolverKind::SparseLu;
    if (value == "umfpack") return HelmholtzFemSolverKind::Umfpack;
    throw std::invalid_argument("unknown reference_solver_kind");
}

void validate_reference_epoch_paper_config(
    const ReferenceEpochPaperConfig &config) {
    if (config.schema_version != reference_epoch_paper_schema_version)
        throw std::invalid_argument("unsupported reference-epoch schema version");
    (void)case_definition(config.case_id);
    if (config.method != "PALOD-reference-epoch"
        && config.method != "PALOD-hybrid-reference-epoch")
        throw std::invalid_argument("unsupported reference-epoch method");
    if (!(config.wavenumber > 0.0) || !std::isfinite(config.wavenumber)
        || !std::isfinite(config.singular_oscillatory_fraction)
        || config.singular_oscillatory_fraction < 0.0
        || config.singular_oscillatory_fraction > 1.0
        || !std::isfinite(config.singular_cutoff_outer_radius)
        || !(config.singular_cutoff_outer_radius > 0.25)
        || config.singular_cutoff_outer_radius > 1.0
        || !std::isfinite(config.smooth_wave_amplitude)
        || config.smooth_wave_amplitude < 0.0
        || config.smooth_wave_amplitude > 1.0
        || !std::isfinite(config.hybrid_minimum_physical_radius)
        || config.hybrid_minimum_physical_radius < 0.0
        || config.hybrid_minimum_physical_radius > 1.0
        || !std::isfinite(config.candidate_regional_minimum_physical_radius)
        || config.candidate_regional_minimum_physical_radius < 0.0
        || config.candidate_regional_minimum_physical_radius > 1.0
        || config.initial_coarse_level < 0
        || config.initial_reference_level <= config.initial_coarse_level
        || config.ell0 < 0 || config.ell_max < config.ell0
        || !(config.theta_loc_usr >= 0.0)
        || !(config.localization_eigen_relative_tolerance > 0.0)
        || config.localization_eigen_relative_tolerance > 1e-3
        || config.localization_eigen_maximum_iterations <= 0
        || !(config.C_rel_usr > 0.0)
        || !(config.theta_H > 0.0 && config.theta_H <= 1.0)
        || !(config.theta_c > 0.0 && config.theta_c <= 1.0)
        || config.candidate_update_stride == 0
        || config.candidate_force_level_gap < 0
        || config.candidate_closure_cost_pool_factor == 0
        || config.patch_symbolic_cache_slots <= 0
        || config.maximum_patch_threads < 0
        || config.reference_validation_stride == 0
        || !(config.q_dual > 0.0 && config.q_dual < 1.0)
        || config.m_dual == 0 || !(config.tau_ep > 0.0)
        || config.reference_refresh_level_gap < 0
        || config.reference_refresh_target_gap < 0
        || (config.reference_refresh_target_gap > 0
            && config.reference_refresh_target_gap
                <= config.reference_refresh_level_gap)
        || !(config.tolerance_reference >= 0.0)
        || !(config.continuity_constant > 0.0)
        || !(config.overlap_constant > 0.0)
        || config.work_limits.maximum_H_steps == 0
        || config.minimum_H_steps_per_epoch != 0
        || config.minimum_solved_points_per_new_epoch
            > config.work_limits.maximum_H_steps
        || config.work_limits.maximum_epochs == 0
        || config.work_limits.maximum_dual_checks == 0
        || !(config.work_limits.maximum_wall_seconds > 0.0)
        || config.quadrature.base_triangle_order <= 0
        || config.quadrature.gaussian_triangle_order <= 0
        || config.quadrature.singular_triangle_order <= 0
        || config.quadrature.max_recursive_subdivisions < 0
        || config.repeat_index < 0)
        throw std::invalid_argument("reference-epoch numerical inputs are invalid");
    if (!config.singularity_hybrid
        && config.reference_validation_stride != 1)
        throw std::invalid_argument(
            "reference_validation_stride may exceed one only for moving-reference runs");
    if (config.singular_solution_profile != "radial-cutoff"
        && config.singular_solution_profile != "boundary-weight-gaussian")
        throw std::invalid_argument("unknown singular_solution_profile");
    if (config.case_id != PaperCase::S
        && (config.singular_oscillatory_fraction != 1.0
            || config.singular_cutoff_outer_radius != 0.5
            || config.singular_quintic_cutoff
            || config.smooth_wave_amplitude != 0.0
            || config.singular_solution_profile != "radial-cutoff"))
        throw std::invalid_argument(
            "S manufactured-solution parameters require case S");
    if (config.singular_solution_profile == "boundary-weight-gaussian"
        && (config.case_id != PaperCase::S
            || config.singular_oscillatory_fraction != 0.0
            || config.singular_cutoff_outer_radius != 0.5
            || config.singular_quintic_cutoff))
        throw std::invalid_argument(
            "boundary-weight-gaussian requires case S, a nonoscillatory corner, and default legacy cutoff fields");
    if (config.singularity_hybrid !=
            (config.method == "PALOD-hybrid-reference-epoch"))
        throw std::invalid_argument("hybrid method and singularity_hybrid disagree");
    if (config.singularity_hybrid && config.case_id != PaperCase::S)
        throw std::invalid_argument("hybrid reference-epoch mode is restricted to case S");
    if (config.candidate_split_regional_marking
        && (config.case_id != PaperCase::S || config.singularity_hybrid))
        throw std::invalid_argument(
            "split regional candidate marking is restricted to standard case-S reference epochs");
    if (config.candidate_split_regional_marking
        != (config.candidate_regional_minimum_physical_radius > 0.0))
        throw std::invalid_argument(
            "candidate regional radius must be positive exactly for split regional marking");
    if (config.singularity_hybrid
        != (config.hybrid_minimum_physical_radius > 0.0))
        throw std::invalid_argument(
            "hybrid_minimum_physical_radius must be positive exactly for hybrid runs");
    if (config.singularity_hybrid
        != (config.hybrid_maximum_corrector_patch_fine_elements > 0))
        throw std::invalid_argument(
            "hybrid_maximum_corrector_patch_fine_elements must be positive exactly for hybrid runs");
    if (config.singularity_hybrid
        && (config.reference_refresh_level_gap != 0
            || config.reference_refresh_target_gap != 0
            || config.minimum_solved_points_per_new_epoch != 0
            || config.work_limits.maximum_epochs
                < config.work_limits.maximum_H_steps)) {
        throw std::invalid_argument(
            "moving-reference hybrid runs require zero reserve guards and maximum_epochs >= maximum_H_steps");
    }
    if (!valid_sha256_digest(config.manuscript_sha256))
        throw std::invalid_argument("manuscript_sha256 must be a sha256 digest");
}

adaptive::ReferenceEpochDriverConfig make_reference_epoch_driver_config(
    const ReferenceEpochPaperConfig &config) {
    validate_reference_epoch_paper_config(config);
    adaptive::ReferenceEpochDriverConfig result;
    result.moving_reference = config.singularity_hybrid;
    result.ell0 = config.ell0;
    result.ell_max = config.ell_max;
    result.tau_loc = config.theta_loc_usr;
    result.q_dual = config.q_dual;
    result.m_dual = config.m_dual;
    result.tau_ep = config.tau_ep;
    result.reference_refresh_level_gap =
        config.reference_refresh_level_gap;
    result.reference_refresh_target_gap =
        config.reference_refresh_target_gap;
    result.candidate_update_stride = config.candidate_update_stride;
    result.candidate_force_level_gap = config.candidate_force_level_gap;
    result.minimum_H_steps_per_epoch =
        config.minimum_H_steps_per_epoch;
    result.minimum_solved_points_per_new_epoch =
        config.minimum_solved_points_per_new_epoch;
    result.tolerance_reference = config.tolerance_reference;
    result.limits = config.work_limits;
    return result;
}

std::string canonical_json(const ReferenceEpochPaperConfig &config) {
    validate_reference_epoch_paper_config(config);
    std::ostringstream out;
    out << "{\"C_rel_usr\":" << number(config.C_rel_usr)
        << ",\"build_hash\":" << json_string(config.build_hash)
        << ",\"case\":" << json_string(to_string(config.case_id))
        << ",\"candidate_closure_cost_aware_marking\":"
        << (config.candidate_closure_cost_aware_marking ? "true" : "false")
        << ",\"candidate_closure_cost_pool_factor\":"
        << config.candidate_closure_cost_pool_factor
        << ",\"candidate_force_level_gap\":"
        << config.candidate_force_level_gap
        << ",\"candidate_regional_minimum_physical_radius\":"
        << number(config.candidate_regional_minimum_physical_radius)
        << ",\"candidate_split_regional_marking\":"
        << (config.candidate_split_regional_marking ? "true" : "false")
        << ",\"candidate_update_stride\":"
        << config.candidate_update_stride
        << ",\"continuity_constant\":" << number(config.continuity_constant)
        << ",\"ell0\":" << config.ell0
        << ",\"ell_max\":" << config.ell_max
        << ",\"git_commit\":" << json_string(config.git_commit)
        << ",\"hybrid_minimum_physical_radius\":"
        << number(config.hybrid_minimum_physical_radius)
        << ",\"hybrid_maximum_corrector_patch_fine_elements\":"
        << config.hybrid_maximum_corrector_patch_fine_elements
        << ",\"initial_coarse_level\":" << config.initial_coarse_level
        << ",\"initial_reference_level\":" << config.initial_reference_level
        << ",\"localization_eigen_maximum_iterations\":"
        << config.localization_eigen_maximum_iterations
        << ",\"localization_eigen_relative_tolerance\":"
        << number(config.localization_eigen_relative_tolerance)
        << ",\"m_dual\":" << config.m_dual
        << ",\"minimum_H_steps_per_epoch\":"
        << config.minimum_H_steps_per_epoch
        << ",\"minimum_solved_points_per_new_epoch\":"
        << config.minimum_solved_points_per_new_epoch
        << ",\"manuscript_sha256\":" << json_string(config.manuscript_sha256)
        << ",\"method\":" << json_string(config.method)
        << ",\"overlap_constant\":" << number(config.overlap_constant)
        << ",\"maximum_patch_threads\":"
        << config.maximum_patch_threads
        << ",\"patch_reuse_identical_factorization\":"
        << (config.patch_reuse_identical_factorization ? "true" : "false")
        << ",\"patch_solver_kind\":"
        << json_string(patch_solver_kind_name(config.patch_solver_kind))
        << ",\"patch_symbolic_cache_slots\":"
        << config.patch_symbolic_cache_slots
        << ",\"reference_solver_kind\":"
        << json_string(helmholtz_fem_solver_kind_name(
               config.reference_solver_kind))
        << ",\"reference_validation_stride\":"
        << config.reference_validation_stride
        << ",\"quadrature\":{\"base_triangle_order\":"
        << config.quadrature.base_triangle_order
        << ",\"gaussian_triangle_order\":"
        << config.quadrature.gaussian_triangle_order
        << ",\"max_recursive_subdivisions\":"
        << config.quadrature.max_recursive_subdivisions
        << ",\"singular_triangle_order\":"
        << config.quadrature.singular_triangle_order << "}"
        << ",\"q_dual\":" << number(config.q_dual)
        << ",\"reference_refresh_level_gap\":"
        << config.reference_refresh_level_gap
        << ",\"reference_refresh_target_gap\":"
        << config.reference_refresh_target_gap
        << ",\"repeat_index\":" << config.repeat_index
        << ",\"schema_version\":" << config.schema_version
        << ",\"singular_cutoff_outer_radius\":"
        << number(config.singular_cutoff_outer_radius)
        << ",\"singular_oscillatory_fraction\":"
        << number(config.singular_oscillatory_fraction)
        << ",\"singular_quintic_cutoff\":"
        << (config.singular_quintic_cutoff ? "true" : "false")
        << ",\"singular_solution_profile\":"
        << json_string(config.singular_solution_profile)
        << ",\"singularity_hybrid\":"
        << (config.singularity_hybrid ? "true" : "false")
        << ",\"smooth_wave_amplitude\":"
        << number(config.smooth_wave_amplitude)
        << ",\"tau_ep\":" << number(config.tau_ep)
        << ",\"theta_H\":" << number(config.theta_H)
        << ",\"theta_c\":" << number(config.theta_c)
        << ",\"theta_loc_usr\":" << number(config.theta_loc_usr)
        << ",\"tolerance_reference\":" << number(config.tolerance_reference)
        << ",\"wavenumber\":" << number(config.wavenumber)
        << ",\"work_limits\":{\"maximum_H_steps\":"
        << config.work_limits.maximum_H_steps
        << ",\"maximum_candidate_unknowns\":"
        << config.work_limits.maximum_candidate_unknowns
        << ",\"maximum_dual_checks\":"
        << config.work_limits.maximum_dual_checks
        << ",\"maximum_epochs\":" << config.work_limits.maximum_epochs
        << ",\"maximum_reference_unknowns\":"
        << config.work_limits.maximum_reference_unknowns
        << ",\"maximum_wall_seconds\":"
        << number(config.work_limits.maximum_wall_seconds) << "}}";
    return out.str();
}

ReferenceEpochPaperConfig parse_reference_epoch_paper_config(
    const std::string_view json) {
    const JsonValue parsed = JsonParser(json).parse();
    const JsonObject &root = as_object(parsed, "root");
    const bool has_patch_solver = root.contains("patch_solver_kind");
    const bool has_singular_solution_profile =
        root.contains("singular_solution_profile");
    const bool has_candidate_update_stride =
        root.contains("candidate_update_stride");
    const bool has_candidate_force_level_gap =
        root.contains("candidate_force_level_gap");
    const bool has_candidate_split_regional_marking =
        root.contains("candidate_split_regional_marking");
    const bool has_candidate_regional_minimum_physical_radius =
        root.contains("candidate_regional_minimum_physical_radius");
    if (has_candidate_split_regional_marking
        != has_candidate_regional_minimum_physical_radius) {
        throw std::invalid_argument(
            "candidate regional marking requires both policy and radius");
    }
    const bool has_reference_solver_kind =
        root.contains("reference_solver_kind");
    const bool has_reference_validation_stride =
        root.contains("reference_validation_stride");
    if (has_candidate_update_stride != has_candidate_force_level_gap) {
        throw std::invalid_argument(
            "candidate batching config requires both stride and force gap");
    }
    const bool has_candidate_closure_cost_aware_marking =
        root.contains("candidate_closure_cost_aware_marking");
    const bool has_candidate_closure_cost_pool_factor =
        root.contains("candidate_closure_cost_pool_factor");
    if (has_candidate_closure_cost_aware_marking
        != has_candidate_closure_cost_pool_factor) {
        throw std::invalid_argument(
            "candidate closure-cost config requires policy and pool factor");
    }
    JsonObject contract_root = root;
    contract_root.erase("singular_solution_profile");
    contract_root.erase("candidate_update_stride");
    contract_root.erase("candidate_force_level_gap");
    contract_root.erase("candidate_split_regional_marking");
    contract_root.erase("candidate_regional_minimum_physical_radius");
    contract_root.erase("candidate_closure_cost_aware_marking");
    contract_root.erase("candidate_closure_cost_pool_factor");
    contract_root.erase("reference_solver_kind");
    contract_root.erase("reference_validation_stride");
    if (has_patch_solver) {
        require_keys(contract_root,
        {"C_rel_usr", "build_hash", "case", "continuity_constant", "ell0",
         "ell_max", "git_commit",
         "hybrid_maximum_corrector_patch_fine_elements",
         "hybrid_minimum_physical_radius",
         "initial_coarse_level",
           "initial_reference_level", "localization_eigen_maximum_iterations",
           "localization_eigen_relative_tolerance", "m_dual", "manuscript_sha256", "method",
         "minimum_H_steps_per_epoch",
         "minimum_solved_points_per_new_epoch",
         "maximum_patch_threads", "overlap_constant",
         "patch_reuse_identical_factorization", "patch_solver_kind",
         "patch_symbolic_cache_slots", "quadrature", "q_dual", "reference_refresh_level_gap",
         "reference_refresh_target_gap",
         "repeat_index", "schema_version",
         "singular_cutoff_outer_radius", "singular_oscillatory_fraction",
         "singular_quintic_cutoff", "singularity_hybrid",
         "smooth_wave_amplitude", "tau_ep", "theta_H", "theta_c",
         "theta_loc_usr", "tolerance_reference", "wavenumber", "work_limits"},
        "root");
    } else {
        // Read pre-solver-policy schema-v6 files with the exact historical
        // key set.  Canonicalization writes the explicit defaults, so every
        // newly frozen run records its numerical backend.
        require_keys(contract_root,
        {"C_rel_usr", "build_hash", "case", "continuity_constant", "ell0",
         "ell_max", "git_commit",
         "hybrid_maximum_corrector_patch_fine_elements",
         "hybrid_minimum_physical_radius", "initial_coarse_level",
         "initial_reference_level", "localization_eigen_maximum_iterations",
         "localization_eigen_relative_tolerance", "m_dual", "manuscript_sha256",
         "method", "minimum_H_steps_per_epoch",
         "minimum_solved_points_per_new_epoch", "overlap_constant", "quadrature",
         "q_dual", "reference_refresh_level_gap", "reference_refresh_target_gap",
         "repeat_index", "schema_version", "singular_cutoff_outer_radius",
         "singular_oscillatory_fraction", "singular_quintic_cutoff",
         "singularity_hybrid", "smooth_wave_amplitude", "tau_ep", "theta_H",
         "theta_c", "theta_loc_usr", "tolerance_reference", "wavenumber",
         "work_limits"}, "root");
    }
    ReferenceEpochPaperConfig config;
    config.schema_version = as_integer(get(root, "schema_version"), "schema_version");
    config.case_id = parse_paper_case(as_string(get(root, "case"), "case"));
    config.method = as_string(get(root, "method"), "method");
    config.wavenumber = as_number(get(root, "wavenumber"), "wavenumber");
    config.hybrid_minimum_physical_radius = as_number(
        get(root, "hybrid_minimum_physical_radius"),
        "hybrid_minimum_physical_radius");
    config.hybrid_maximum_corrector_patch_fine_elements =
        static_cast<std::size_t>(as_uint64(
            get(root, "hybrid_maximum_corrector_patch_fine_elements"),
            "hybrid_maximum_corrector_patch_fine_elements"));
    config.singular_oscillatory_fraction = as_number(
        get(root, "singular_oscillatory_fraction"),
        "singular_oscillatory_fraction");
    config.singular_cutoff_outer_radius = as_number(
        get(root, "singular_cutoff_outer_radius"),
        "singular_cutoff_outer_radius");
    config.singular_quintic_cutoff = as_bool(
        get(root, "singular_quintic_cutoff"),
        "singular_quintic_cutoff");
    config.smooth_wave_amplitude = as_number(
        get(root, "smooth_wave_amplitude"), "smooth_wave_amplitude");
    if (has_singular_solution_profile) {
        config.singular_solution_profile = as_string(
            get(root, "singular_solution_profile"),
            "singular_solution_profile");
    }
    config.initial_coarse_level = as_integer(
        get(root, "initial_coarse_level"), "initial_coarse_level");
    config.initial_reference_level = as_integer(
        get(root, "initial_reference_level"), "initial_reference_level");
    config.localization_eigen_maximum_iterations = as_integer(
        get(root, "localization_eigen_maximum_iterations"),
        "localization_eigen_maximum_iterations");
    config.localization_eigen_relative_tolerance = as_number(
        get(root, "localization_eigen_relative_tolerance"),
        "localization_eigen_relative_tolerance");
    config.singularity_hybrid = as_bool(
        get(root, "singularity_hybrid"), "singularity_hybrid");
    if (has_candidate_split_regional_marking) {
        config.candidate_split_regional_marking = as_bool(
            get(root, "candidate_split_regional_marking"),
            "candidate_split_regional_marking");
        config.candidate_regional_minimum_physical_radius = as_number(
            get(root, "candidate_regional_minimum_physical_radius"),
            "candidate_regional_minimum_physical_radius");
    }
    config.ell0 = as_integer(get(root, "ell0"), "ell0");
    config.ell_max = as_integer(get(root, "ell_max"), "ell_max");
    config.theta_loc_usr = as_number(get(root, "theta_loc_usr"), "theta_loc_usr");
    config.C_rel_usr = as_number(get(root, "C_rel_usr"), "C_rel_usr");
    config.theta_H = as_number(get(root, "theta_H"), "theta_H");
    config.theta_c = as_number(get(root, "theta_c"), "theta_c");
    if (has_candidate_update_stride) {
        config.candidate_update_stride = static_cast<std::size_t>(as_uint64(
            get(root, "candidate_update_stride"),
            "candidate_update_stride"));
        config.candidate_force_level_gap = as_integer(
            get(root, "candidate_force_level_gap"),
            "candidate_force_level_gap");
    }
    if (has_candidate_closure_cost_aware_marking) {
        config.candidate_closure_cost_aware_marking = as_bool(
            get(root, "candidate_closure_cost_aware_marking"),
            "candidate_closure_cost_aware_marking");
        config.candidate_closure_cost_pool_factor =
            static_cast<std::size_t>(as_uint64(
                get(root, "candidate_closure_cost_pool_factor"),
                "candidate_closure_cost_pool_factor"));
    }
    if (has_patch_solver) {
        config.patch_solver_kind = parse_patch_solver_kind(
            as_string(get(root, "patch_solver_kind"), "patch_solver_kind"));
        config.patch_symbolic_cache_slots = as_integer(
            get(root, "patch_symbolic_cache_slots"),
            "patch_symbolic_cache_slots");
        config.patch_reuse_identical_factorization = as_bool(
            get(root, "patch_reuse_identical_factorization"),
            "patch_reuse_identical_factorization");
        config.maximum_patch_threads = as_integer(
            get(root, "maximum_patch_threads"), "maximum_patch_threads");
    }
    if (has_reference_solver_kind) {
        config.reference_solver_kind = parse_reference_solver_kind(
            as_string(get(root, "reference_solver_kind"),
                      "reference_solver_kind"));
    }
    if (has_reference_validation_stride) {
        config.reference_validation_stride = static_cast<std::size_t>(
            as_uint64(get(root, "reference_validation_stride"),
                      "reference_validation_stride"));
    }
    config.q_dual = as_number(get(root, "q_dual"), "q_dual");
    config.reference_refresh_level_gap = as_integer(
        get(root, "reference_refresh_level_gap"),
        "reference_refresh_level_gap");
    config.reference_refresh_target_gap = as_integer(
        get(root, "reference_refresh_target_gap"),
        "reference_refresh_target_gap");
    config.m_dual = static_cast<std::size_t>(as_uint64(get(root, "m_dual"), "m_dual"));
    config.minimum_H_steps_per_epoch = static_cast<std::size_t>(as_uint64(
        get(root, "minimum_H_steps_per_epoch"),
        "minimum_H_steps_per_epoch"));
    config.minimum_solved_points_per_new_epoch =
        static_cast<std::size_t>(as_uint64(
            get(root, "minimum_solved_points_per_new_epoch"),
            "minimum_solved_points_per_new_epoch"));
    config.tau_ep = as_number(get(root, "tau_ep"), "tau_ep");
    config.tolerance_reference = as_number(
        get(root, "tolerance_reference"), "tolerance_reference");
    config.continuity_constant = as_number(
        get(root, "continuity_constant"), "continuity_constant");
    config.overlap_constant = as_number(
        get(root, "overlap_constant"), "overlap_constant");
    const JsonObject &quadrature = as_object(
        get(root, "quadrature"), "quadrature");
    require_keys(quadrature,
        {"base_triangle_order", "gaussian_triangle_order",
         "max_recursive_subdivisions", "singular_triangle_order"},
        "quadrature");
    config.quadrature.base_triangle_order = as_integer(
        get(quadrature, "base_triangle_order"), "base_triangle_order");
    config.quadrature.gaussian_triangle_order = as_integer(
        get(quadrature, "gaussian_triangle_order"), "gaussian_triangle_order");
    config.quadrature.max_recursive_subdivisions = as_integer(
        get(quadrature, "max_recursive_subdivisions"),
        "max_recursive_subdivisions");
    config.quadrature.singular_triangle_order = as_integer(
        get(quadrature, "singular_triangle_order"), "singular_triangle_order");
    config.repeat_index = as_integer(get(root, "repeat_index"), "repeat_index");
    config.git_commit = as_string(get(root, "git_commit"), "git_commit");
    config.build_hash = as_string(get(root, "build_hash"), "build_hash");
    config.manuscript_sha256 = as_string(
        get(root, "manuscript_sha256"), "manuscript_sha256");
    const JsonObject &limits = as_object(get(root, "work_limits"), "work_limits");
    require_keys(limits,
        {"maximum_H_steps", "maximum_candidate_unknowns", "maximum_dual_checks",
         "maximum_epochs", "maximum_reference_unknowns", "maximum_wall_seconds"},
        "work_limits");
    config.work_limits.maximum_H_steps = static_cast<std::size_t>(as_uint64(
        get(limits, "maximum_H_steps"), "maximum_H_steps"));
    config.work_limits.maximum_candidate_unknowns = static_cast<std::size_t>(as_uint64(
        get(limits, "maximum_candidate_unknowns"), "maximum_candidate_unknowns"));
    config.work_limits.maximum_dual_checks = static_cast<std::size_t>(as_uint64(
        get(limits, "maximum_dual_checks"), "maximum_dual_checks"));
    config.work_limits.maximum_epochs = static_cast<std::size_t>(as_uint64(
        get(limits, "maximum_epochs"), "maximum_epochs"));
    config.work_limits.maximum_reference_unknowns = static_cast<std::size_t>(as_uint64(
        get(limits, "maximum_reference_unknowns"), "maximum_reference_unknowns"));
    config.work_limits.maximum_wall_seconds = as_number(
        get(limits, "maximum_wall_seconds"), "maximum_wall_seconds");
    validate_reference_epoch_paper_config(config);
    return config;
}

std::string canonical_config_hash(const ReferenceEpochPaperConfig &config) {
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

std::string make_run_id(const ReferenceEpochPaperConfig &config) {
    validate_reference_epoch_paper_config(config);
    std::ostringstream out;
    out << to_string(config.case_id) << '_' << safe_component(config.method)
        << "_k" << static_cast<int>(config.wavenumber)
        << "_r" << config.repeat_index << '_' << canonical_config_hash(config);
    return out.str();
}

bool is_reference_epoch_paper_config(const std::string_view json) {
    try {
        const JsonValue parsed = JsonParser(json).parse();
        const JsonObject &root = as_object(parsed, "root");
        return root.contains("schema_version")
            && as_integer(get(root, "schema_version"), "schema_version")
                == reference_epoch_paper_schema_version;
    } catch (...) {
        return false;
    }
}

} // namespace lod2d::helmholtz::experiments

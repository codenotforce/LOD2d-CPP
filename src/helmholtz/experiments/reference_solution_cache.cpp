#include "helmholtz/experiments/reference_solution_cache.h"

#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <system_error>
#include <type_traits>

namespace lod2d::helmholtz::experiments {
namespace {

constexpr std::uint64_t fnv_offset = 14695981039346656037ull;
constexpr std::uint64_t fnv_prime = 1099511628211ull;
constexpr std::array<char, 16> magic{{
    'L','O','D','2','D','R','E','F','C','A','C','H','E','0','0','1'}};

void hash_bytes(std::uint64_t &hash, const void *data, const std::size_t size) {
    const auto *bytes = static_cast<const unsigned char *>(data);
    for (std::size_t index = 0; index < size; ++index) {
        hash ^= bytes[index];
        hash *= fnv_prime;
    }
}

template <class Scalar>
void hash_scalar(std::uint64_t &hash, const Scalar &value) {
    static_assert(std::is_trivially_copyable_v<Scalar>);
    hash_bytes(hash, std::addressof(value), sizeof(value));
}

void hash_string(std::uint64_t &hash, const std::string &value) {
    hash_scalar(hash, value.size());
    hash_bytes(hash, value.data(), value.size());
}

template <class Scalar>
void hash_vector(std::uint64_t &hash, const std::vector<Scalar> &values) {
    hash_scalar(hash, values.size());
    for (const Scalar &value : values) hash_scalar(hash, value);
}

template <class Scalar>
void hash_sparse(
    std::uint64_t &hash,
    const Eigen::SparseMatrix<Scalar> &matrix) {
    hash_scalar(hash, matrix.rows());
    hash_scalar(hash, matrix.cols());
    hash_scalar(hash, matrix.nonZeros());
    for (int outer = 0; outer < matrix.outerSize(); ++outer) {
        for (typename Eigen::SparseMatrix<Scalar>::InnerIterator it(matrix, outer);
             it; ++it) {
            hash_scalar(hash, it.row());
            hash_scalar(hash, it.col());
            hash_scalar(hash, it.value());
        }
    }
}

void hash_mesh(std::uint64_t &hash, const TriMesh &mesh) {
    hash_scalar(hash, mesh.nodes.size());
    for (const Point2 &node : mesh.nodes) {
        hash_scalar(hash, node[0]);
        hash_scalar(hash, node[1]);
    }
    hash_vector(hash, mesh.elems);
    hash_vector(hash, mesh.dirichlet);
    hash_scalar(hash, mesh.boundary_edges.size());
    for (const BoundaryEdge &edge : mesh.boundary_edges) {
        hash_scalar(hash, edge.nodes);
        hash_scalar(hash, edge.tag);
    }
}

void hash_complex_vector(std::uint64_t &hash, const ComplexVector &values) {
    hash_scalar(hash, values.size());
    for (Eigen::Index index = 0; index < values.size(); ++index) {
        hash_scalar(hash, values[index].real());
        hash_scalar(hash, values[index].imag());
    }
}

std::string hex64(const std::uint64_t value) {
    std::ostringstream out;
    out << std::hex << std::setfill('0') << std::setw(16) << value;
    return out.str();
}

template <class Scalar>
bool read_exact(std::ifstream &input, Scalar &value) {
    input.read(reinterpret_cast<char *>(std::addressof(value)), sizeof(value));
    return input.good();
}

template <class Scalar>
void write_exact(std::ofstream &output, const Scalar &value) {
    output.write(
        reinterpret_cast<const char *>(std::addressof(value)), sizeof(value));
    if (!output) throw std::runtime_error("failed to write reference cache");
}

std::uint64_t solution_checksum(const ComplexVector &solution) {
    std::uint64_t hash = fnv_offset;
    hash_complex_vector(hash, solution);
    return hash;
}

std::filesystem::path cache_path(
    const std::filesystem::path &directory,
    const std::string &key) {
    return directory / ("reference_" + key + ".bin");
}

} // namespace

std::string reference_solution_cache_key(
    const TriMesh &mesh,
    const HelmholtzOperators &operators,
    const ComplexVector &load,
    const std::string &identity) {
    if (operators.system.rows() != static_cast<int>(mesh.nodes.size())
        || load.size() != static_cast<Eigen::Index>(mesh.nodes.size()))
        throw std::invalid_argument(
            "reference cache inputs have inconsistent dimensions");
    std::uint64_t hash = fnv_offset;
    hash_string(hash, "lod2d-reference-sparse-lu-v1");
    hash_string(hash, identity);
    hash_mesh(hash, mesh);
    hash_scalar(hash, operators.wavenumber);
    hash_scalar(hash, operators.boundary_beta);
    hash_vector(hash, operators.diffusion);
    hash_vector(hash, operators.refractive_index);
    hash_vector(hash, operators.dirichlet_nodes);
    hash_sparse(hash, operators.system);
    hash_complex_vector(hash, load);
    return hex64(hash);
}

ReferenceSolutionCacheLookup load_reference_solution_cache(
    const std::filesystem::path &directory,
    const std::string &key,
    const Eigen::Index expected_size) {
    if (key.size() != 16 || expected_size < 0)
        throw std::invalid_argument("invalid reference cache lookup");
    ReferenceSolutionCacheLookup result;
    result.key = key;
    result.path = cache_path(directory, key);
    std::ifstream input(result.path, std::ios::binary);
    if (!input) return result;

    std::array<char, magic.size()> file_magic{};
    input.read(file_magic.data(), static_cast<std::streamsize>(file_magic.size()));
    std::uint64_t count = 0;
    std::uint64_t checksum = 0;
    if (!input || file_magic != magic || !read_exact(input, count)
        || !read_exact(input, checksum)
        || count != static_cast<std::uint64_t>(expected_size)
        || count > static_cast<std::uint64_t>(
            std::numeric_limits<Eigen::Index>::max()))
        return result;

    ComplexVector solution(static_cast<Eigen::Index>(count));
    for (Eigen::Index index = 0; index < solution.size(); ++index) {
        double real = 0.0;
        double imaginary = 0.0;
        if (!read_exact(input, real) || !read_exact(input, imaginary))
            return result;
        solution[index] = Complex(real, imaginary);
    }
    char extra = 0;
    if (input.read(&extra, 1) || !input.eof() || !solution.allFinite()
        || solution_checksum(solution) != checksum)
        return result;
    result.solution = std::move(solution);
    return result;
}

void store_reference_solution_cache(
    const std::filesystem::path &directory,
    const std::string &key,
    const ComplexVector &solution) {
    if (key.size() != 16 || !solution.allFinite())
        throw std::invalid_argument("invalid reference cache store");
    std::filesystem::create_directories(directory);
    const std::filesystem::path destination = cache_path(directory, key);
    const std::filesystem::path temporary = destination.string() + ".tmp."
        + std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count());
    try {
        std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
        if (!output)
            throw std::runtime_error("cannot create reference cache file");
        output.write(magic.data(), static_cast<std::streamsize>(magic.size()));
        const std::uint64_t count = static_cast<std::uint64_t>(solution.size());
        const std::uint64_t checksum = solution_checksum(solution);
        write_exact(output, count);
        write_exact(output, checksum);
        for (Eigen::Index index = 0; index < solution.size(); ++index) {
            write_exact(output, solution[index].real());
            write_exact(output, solution[index].imag());
        }
        output.close();
        if (!output) throw std::runtime_error("failed to close reference cache");
        std::error_code error;
        std::filesystem::rename(temporary, destination, error);
        if (error) {
            // Another process may have won the race with identical content.
            if (!std::filesystem::is_regular_file(destination))
                throw std::filesystem::filesystem_error(
                    "cannot publish reference cache", temporary, destination, error);
            std::filesystem::remove(temporary, error);
        }
    } catch (...) {
        std::error_code ignored;
        std::filesystem::remove(temporary, ignored);
        throw;
    }
}

} // namespace lod2d::helmholtz::experiments

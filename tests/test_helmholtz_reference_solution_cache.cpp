#include "helmholtz/experiments/reference_solution_cache.h"
#include "helmholtz/model.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

using namespace lod2d;
using namespace lod2d::helmholtz;
using namespace lod2d::helmholtz::experiments;

namespace {

void require(const bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

std::filesystem::path fresh_directory() {
    const std::filesystem::path path =
        std::filesystem::temp_directory_path()
        / "lod2d-reference-solution-cache-test";
    std::error_code ignored;
    std::filesystem::remove_all(path, ignored);
    std::filesystem::create_directories(path);
    return path;
}

void verify_round_trip_and_fail_closed_identity() {
    const TriMesh mesh = make_helmholtz_unit_square_mesh();
    const HelmholtzOperators operators =
        assemble_helmholtz_operators(mesh, 4.0);
    ComplexVector load = ComplexVector::Zero(
        static_cast<Eigen::Index>(mesh.nodes.size()));
    for (Eigen::Index index = 0; index < load.size(); ++index)
        load[index] = Complex(0.25 * index, -0.5);
    const ComplexVector solution = solve_helmholtz_fem(operators, load);
    const std::filesystem::path directory = fresh_directory();

    const std::string key = reference_solution_cache_key(
        mesh, operators, load, "commit-a/build-a");
    require(!load_reference_solution_cache(directory, key, solution.size()).solution,
            "empty reference cache unexpectedly hit");
    store_reference_solution_cache(directory, key, solution);
    const ReferenceSolutionCacheLookup loaded =
        load_reference_solution_cache(directory, key, solution.size());
    require(loaded.solution
                && (*loaded.solution - solution).norm() == 0.0,
            "reference cache did not round-trip exactly");

    ComplexVector changed_load = load;
    changed_load[0] += Complex(1e-12, 0.0);
    require(reference_solution_cache_key(
                mesh, operators, changed_load, "commit-a/build-a") != key,
            "reference cache key ignored the load");
    require(reference_solution_cache_key(
                mesh, operators, load, "commit-b/build-a") != key,
            "reference cache key ignored the code identity");
    const HelmholtzOperators changed_operators =
        assemble_helmholtz_operators(mesh, 5.0);
    require(reference_solution_cache_key(
                mesh, changed_operators, load, "commit-a/build-a") != key,
            "reference cache key ignored the PDE");

    std::filesystem::remove_all(directory);
}

void verify_corruption_and_dimension_miss() {
    const TriMesh mesh = make_helmholtz_unit_square_mesh();
    const HelmholtzOperators operators =
        assemble_helmholtz_operators(mesh, 4.0);
    const ComplexVector load = ComplexVector::Ones(
        static_cast<Eigen::Index>(mesh.nodes.size()));
    const ComplexVector solution = solve_helmholtz_fem(operators, load);
    const std::filesystem::path directory = fresh_directory();
    const std::string key = reference_solution_cache_key(
        mesh, operators, load, "commit-a/build-a");
    store_reference_solution_cache(directory, key, solution);
    require(!load_reference_solution_cache(
                directory, key, solution.size() + 1).solution,
            "reference cache accepted the wrong vector dimension");

    const std::filesystem::path path = directory / ("reference_" + key + ".bin");
    const auto original_size = std::filesystem::file_size(path);
    std::filesystem::resize_file(path, original_size - 1);
    require(!load_reference_solution_cache(directory, key, solution.size()).solution,
            "reference cache accepted a truncated payload");
    std::filesystem::remove_all(directory);
}

} // namespace

int main() {
    try {
        verify_round_trip_and_fail_closed_identity();
        verify_corruption_and_dimension_miss();
        std::cout << "Helmholtz reference solution cache tests passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_reference_solution_cache failed: "
                  << error.what() << '\n';
        return 1;
    }
}

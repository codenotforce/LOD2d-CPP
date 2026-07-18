#include "helmholtz/adaptive/estimator.h"
#include "helmholtz/adaptive/hierarchy.h"
#include "helmholtz/model.h"

#include <algorithm>
#include <cmath>
#include <array>
#include <iostream>
#include <stdexcept>

using namespace lod2d;
using namespace lod2d::helmholtz;
using namespace lod2d::helmholtz::adaptive;

namespace {

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

ComplexFunction source_function() {
    return [](const Point2 &point) {
        const double dx = point.x() - 0.31;
        const double dy = point.y() - 0.63;
        return Complex(std::exp(-45.0 * (dx * dx + dy * dy)), 0.15 * point.x());
    };
}

using TriangleSignature = std::array<long long, 6>;

std::vector<TriangleSignature> canonical_triangles(const TriMesh &mesh) {
    std::vector<TriangleSignature> result;
    result.reserve(mesh.elems.size());
    for (const Triangle &triangle : mesh.elems) {
        std::array<std::pair<long long, long long>, 3> vertices;
        for (int i = 0; i < 3; ++i) {
            const Point2 &point = mesh.nodes[triangle[i]];
            vertices[i] = {
                std::llround(point.x() * 1e12),
                std::llround(point.y() * 1e12)};
        }
        std::sort(vertices.begin(), vertices.end());
        TriangleSignature signature{};
        for (int i = 0; i < 3; ++i) {
            signature[2 * i] = vertices[i].first;
            signature[2 * i + 1] = vertices[i].second;
        }
        result.push_back(signature);
    }
    std::sort(result.begin(), result.end());
    return result;
}

void verify_hierarchy_and_estimator() {
    AdaptiveMeshHierarchy hierarchy(make_helmholtz_unit_square_mesh(), 1, 4);
    const NestedFineMesh initial_fine = hierarchy.build_nested_fine_mesh();
    require(std::all_of(
                initial_fine.element_levels.begin(),
                initial_fine.element_levels.end(),
                [](int level) { return level == 4; }),
            "fixed fine completion has the wrong level");
    const TriMesh uniform_fine =
        refine_mesh_nvb(make_helmholtz_unit_square_mesh(), 4).mesh;
    require(canonical_triangles(initial_fine.refinement.mesh) == canonical_triangles(uniform_fine),
            "initial fixed fine completion differs from the uniform master mesh");
    require(initial_fine.refinement.mesh.nodes.size() == uniform_fine.nodes.size(),
            "initial fixed fine completion contains duplicate node coordinates");

    hierarchy.refine({0});
    require(*std::min_element(hierarchy.coarse_levels().begin(), hierarchy.coarse_levels().end()) == 1,
            "local NVB unexpectedly refined every coarse element");
    require(*std::max_element(hierarchy.coarse_levels().begin(), hierarchy.coarse_levels().end()) > 1,
            "local NVB did not refine the marked region");
    require(hierarchy.coarse_element_ids().size() == hierarchy.coarse_mesh().elems.size(),
            "stable coarse element IDs do not match the mesh");
    const NestedFineMesh adaptive_fine = hierarchy.build_nested_fine_mesh();
    require(canonical_triangles(adaptive_fine.refinement.mesh) == canonical_triangles(uniform_fine),
            "adaptive completion changed the fixed fine mesh geometry");

    HelmholtzProblemConfig config;
    config.H = 1;
    require(adaptive_fine.refinement.mesh.nodes.size() == uniform_fine.nodes.size(),
            "adaptive completion contains duplicate node coordinates");
    for (int step = 0; step < 3; ++step) {
        std::vector<int> marked;
        for (int t = 0; t < static_cast<int>(hierarchy.coarse_levels().size()); ++t) {
            if (hierarchy.coarse_levels()[t] < hierarchy.fine_level()
                && (t + step) % 3 == 0) {
                marked.push_back(t);
            }
        }
        if (marked.empty()) break;
        hierarchy.refine(marked);
        const NestedFineMesh completed = hierarchy.build_nested_fine_mesh();
        require(
            canonical_triangles(completed.refinement.mesh) == canonical_triangles(uniform_fine),
            "repeated adaptive completion changed the fixed fine mesh geometry");
        require(completed.refinement.mesh.nodes.size() == uniform_fine.nodes.size(),
                "repeated adaptive completion contains duplicate node coordinates");
    }
    config.h = 4;
    config.ell = 1;
    config.wavenumber = 3.0;
    HelmholtzLodModel model = HelmholtzLodModel::build_adaptive(
        config, hierarchy.coarse_mesh(), hierarchy.coarse_levels());
    require(model.problem().coarse.elems.size() == hierarchy.coarse_mesh().elems.size(),
            "adaptive model changed the coarse mesh");
    require(model.problem().fine_level == 4, "adaptive model lost the fixed fine level");

    const ComplexFunction source = source_function();
    const ComplexVector load = assemble_helmholtz_load(model.problem().fine, source);
    const HelmholtzLodSolution solution = model.solve_load(load);
    const ComplexVector reference = model.solve_fine_reference(load);
    const ComplexVector error = reference - solution.fine_values;
    require(discrete_energy_norm(model.operators(), error) > 0.0,
            "adaptive test unexpectedly has zero LOD error");

    const HelmholtzResidualContributions residual =
        assemble_helmholtz_residual_contributions(
            model.problem(), model.operators(), solution.fine_values, load, source);
    std::cout << "Residual reconstruction difference: "
              << residual.algebraic_relative_difference << '\n';
    require(residual.algebraic_relative_difference < 1e-10,
            "integrated residual disagrees with b-Au");
    const HelmholtzIndicatorSet indicators = build_helmholtz_indicators(
        model.problem(), residual);
    require(indicators.fine > 0.0 && indicators.mixed > 0.0 && indicators.macro > 0.0,
            "adaptive residual indicators must be positive");
    require(indicators.fine_squared.size() == model.problem().coarse.elems.size(),
            "coarse indicator count is incorrect");

    const std::vector<double> dual = build_local_dual_indicators(
        model.problem(), model.operators(), residual, 1);
    require(dual.size() == model.problem().coarse.elems.size(),
            "local dual indicator count is incorrect");
    require(*std::max_element(dual.begin(), dual.end()) > 0.0,
            "local dual indicators are all zero");

    std::vector<char> eligible(hierarchy.coarse_levels().size(), true);
    const std::vector<int> marked = mark_doerfler(
        indicators.mixed_squared, 0.5, eligible);
    require(!marked.empty(), "Doerfler marking returned an empty set");
}

void verify_doerfler_minimal_set() {
    const std::vector<int> marked = mark_doerfler({4.0, 3.0, 2.0, 1.0}, 0.6);
    require(marked == std::vector<int>({0, 1}),
            "Doerfler marking is not deterministic or minimal");
}

} // namespace

int main() {
    try {
        verify_hierarchy_and_estimator();
        verify_doerfler_minimal_set();
        std::cout << "Adaptive Helmholtz hierarchy and residual estimators passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_helmholtz_adaptive failed: " << error.what() << '\n';
        return 1;
    }
}

#include "helmholtz/adaptive/singularity_hybrid.h"

#include "helmholtz/adaptive/estimator.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <queue>
#include <stdexcept>

namespace lod2d::helmholtz::adaptive {
namespace {

std::vector<std::vector<int>> element_graph(const TriMesh &mesh) {
    std::vector<std::vector<int>> incidence(mesh.nodes.size());
    for (int element = 0; element < static_cast<int>(mesh.elems.size());
         ++element) {
        for (int node : mesh.elems[element]) incidence[node].push_back(element);
    }
    std::vector<std::vector<int>> graph(mesh.elems.size());
    for (const std::vector<int> &star : incidence) {
        for (int left : star)
            for (int right : star)
                if (left != right) graph[left].push_back(right);
    }
    for (std::vector<int> &neighbors : graph) {
        std::sort(neighbors.begin(), neighbors.end());
        neighbors.erase(
            std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
    }
    return graph;
}

bool triangle_contains_point(
    const TriMesh &mesh, int element, const Point2 &point, double tolerance) {
    const Triangle &triangle = mesh.elems[element];
    const Point2 &p0 = mesh.nodes[triangle[0]];
    Eigen::Matrix2d jacobian;
    jacobian.col(0) = mesh.nodes[triangle[1]] - p0;
    jacobian.col(1) = mesh.nodes[triangle[2]] - p0;
    if (std::abs(jacobian.determinant()) <= tolerance)
        throw std::invalid_argument("singular-region mesh has a degenerate cell");
    const Eigen::Vector2d local = jacobian.inverse() * (point - p0);
    const Eigen::Vector3d barycentric(
        1.0 - local.x() - local.y(), local.x(), local.y());
    return barycentric.minCoeff() >= -tolerance
        && barycentric.maxCoeff() <= 1.0 + tolerance;
}

std::vector<int> graph_neighborhood(
    const std::vector<std::vector<int>> &graph,
    const std::vector<int> &seeds,
    int layers) {
    std::vector<int> distance(graph.size(), std::numeric_limits<int>::max());
    std::queue<int> pending;
    for (int seed : seeds) {
        distance[seed] = 0;
        pending.push(seed);
    }
    while (!pending.empty()) {
        const int current = pending.front();
        pending.pop();
        if (distance[current] >= layers) continue;
        for (int neighbor : graph[current]) {
            if (distance[neighbor] <= distance[current] + 1) continue;
            distance[neighbor] = distance[current] + 1;
            pending.push(neighbor);
        }
    }
    std::vector<int> result;
    for (int element = 0; element < static_cast<int>(distance.size());
         ++element)
        if (distance[element] <= layers) result.push_back(element);
    return result;
}

double squared_distance_to_segment(
    const Point2 &point, const Point2 &left, const Point2 &right) {
    const Point2 direction = right - left;
    const double length_squared = direction.squaredNorm();
    if (!(length_squared > 0.0)) return (point - left).squaredNorm();
    const double parameter = std::clamp(
        (point - left).dot(direction) / length_squared, 0.0, 1.0);
    return (point - (left + parameter * direction)).squaredNorm();
}

double distance_to_triangle(
    const TriMesh &mesh, int element, const Point2 &point, double tolerance) {
    if (triangle_contains_point(mesh, element, point, tolerance)) return 0.0;
    const Triangle &triangle = mesh.elems[element];
    double squared = std::numeric_limits<double>::infinity();
    for (int edge = 0; edge < 3; ++edge) {
        squared = std::min(
            squared,
            squared_distance_to_segment(
                point, mesh.nodes[triangle[edge]],
                mesh.nodes[triangle[(edge + 1) % 3]]));
    }
    return std::sqrt(squared);
}

std::vector<int> graph_distances(
    const std::vector<std::vector<int>> &graph,
    const std::vector<int> &seeds) {
    std::vector<int> distance(graph.size(), std::numeric_limits<int>::max());
    std::queue<int> pending;
    for (int seed : seeds) {
        distance[seed] = 0;
        pending.push(seed);
    }
    while (!pending.empty()) {
        const int current = pending.front();
        pending.pop();
        for (int neighbor : graph[current]) {
            if (distance[neighbor] <= distance[current] + 1) continue;
            distance[neighbor] = distance[current] + 1;
            pending.push(neighbor);
        }
    }
    return distance;
}

SingularRegionClassification classify_impl(
    const TriMesh &mesh,
    const std::vector<Point2> &singular_points,
    int ell,
    double minimum_physical_radius,
    double geometric_tolerance) {
    if (singular_points.empty() || ell < 0
        || !(minimum_physical_radius >= 0.0)
        || !std::isfinite(minimum_physical_radius)
        || !(geometric_tolerance > 0.0)) {
        throw std::invalid_argument("singular-region inputs are invalid");
    }
    SingularRegionClassification result;
    result.corrector_ell = ell;
    result.minimum_physical_radius = minimum_physical_radius;
    for (int element = 0; element < static_cast<int>(mesh.elems.size());
         ++element) {
        if (std::any_of(
                singular_points.begin(), singular_points.end(),
                [&](const Point2 &point) {
                    return triangle_contains_point(
                        mesh, element, point, geometric_tolerance);
                })) {
            result.seed_elements.push_back(element);
        }
    }
    if (result.seed_elements.empty())
        throw std::invalid_argument("singular set does not meet the coarse mesh");

    const auto graph = element_graph(mesh);
    const std::vector<int> distance = graph_distances(
        graph, result.seed_elements);
    int physical_layers = 0;
    if (minimum_physical_radius > 0.0) {
        for (int element = 0; element < static_cast<int>(mesh.elems.size());
             ++element) {
            const bool intersects_fixed_ball = std::any_of(
                singular_points.begin(), singular_points.end(),
                [&](const Point2 &point) {
                    return distance_to_triangle(
                               mesh, element, point, geometric_tolerance)
                        <= minimum_physical_radius + geometric_tolerance;
                });
            if (intersects_fixed_ball) {
                if (distance[element] == std::numeric_limits<int>::max())
                    throw std::runtime_error(
                        "physical singular neighborhood is graph-disconnected");
                physical_layers = std::max(
                    physical_layers, (distance[element] + 1) / 2);
            }
        }
    }
    result.l_s = std::max(ell, physical_layers);
    result.omega_s_elements = graph_neighborhood(
        graph, result.seed_elements, result.l_s);
    result.omega_f_elements = graph_neighborhood(
        graph, result.seed_elements, 2 * result.l_s);
    result.in_omega_s.assign(mesh.elems.size(), false);
    result.in_omega_f.assign(mesh.elems.size(), false);
    result.in_regular.assign(mesh.elems.size(), true);
    for (int element : result.omega_s_elements)
        result.in_omega_s[element] = true;
    for (int element : result.omega_f_elements) {
        result.in_omega_f[element] = true;
        result.in_regular[element] = false;
    }
    for (int element = 0; element < static_cast<int>(mesh.elems.size());
         ++element)
        if (result.in_regular[element]) result.regular_elements.push_back(element);

    double covered = std::numeric_limits<double>::infinity();
    for (int element : result.regular_elements) {
        for (const Point2 &point : singular_points) {
            covered = std::min(
                covered,
                distance_to_triangle(
                    mesh, element, point, geometric_tolerance));
        }
    }
    if (!std::isfinite(covered)) {
        covered = 0.0;
        for (const Point2 &node : mesh.nodes)
            for (const Point2 &point : singular_points)
                covered = std::max(covered, (node - point).norm());
    }
    result.covered_physical_radius = covered;
    if (covered + geometric_tolerance < minimum_physical_radius)
        throw std::runtime_error(
            "hybrid graph neighborhood does not cover the physical radius");
    return result;
}

std::size_t reference_children(
    const Eigen::SparseMatrix<double> &embedding, int coarse_element) {
    std::size_t count = 0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(
             embedding, coarse_element); it; ++it) {
        if (std::abs(it.value()) > 1e-14) ++count;
    }
    return count;
}

} // namespace

SingularRegionClassification classify_singular_regions(
    const TriMesh &mesh,
    const std::vector<Point2> &singular_points,
    int ell,
    double geometric_tolerance) {
    return classify_impl(
        mesh, singular_points, ell, 0.0, geometric_tolerance);
}

SingularRegionClassification classify_singular_regions_with_physical_radius(
    const TriMesh &mesh,
    const std::vector<Point2> &singular_points,
    int ell,
    double minimum_physical_radius,
    double geometric_tolerance) {
    if (!(minimum_physical_radius > 0.0))
        throw std::invalid_argument(
            "physical singular-neighborhood radius must be positive");
    return classify_impl(
        mesh, singular_points, ell, minimum_physical_radius,
        geometric_tolerance);
}

bool hybrid_reference_matching_holds(
    const ReferenceEpochHierarchy &hierarchy,
    const SingularRegionClassification &regions) {
    if (regions.in_omega_f.size() != hierarchy.coarse_mesh().elems.size())
        return false;
    for (int element : regions.omega_f_elements) {
        if (reference_children(
                hierarchy.coarse_elements_to_reference(), element) != 1)
            return false;
    }
    return true;
}

HybridMatchingResult restore_hybrid_reference_matching(
    ReferenceEpochHierarchy &hierarchy,
    const std::vector<Point2> &singular_points,
    int ell,
    std::size_t maximum_rounds) {
    if (maximum_rounds == 0 || hierarchy.has_proposed_coarse_refinement())
        throw std::invalid_argument(
            "hybrid matching requires a positive limit and no pending proposal");
    const std::uint64_t reference_version = hierarchy.reference_mesh_version();
    const std::size_t reference_elements = hierarchy.reference_mesh().elems.size();
    HybridMatchingResult result;
    for (;;) {
        result.regions = classify_singular_regions(
            hierarchy.coarse_mesh(), singular_points, ell);
        std::vector<int> marked;
        for (int element : result.regions.omega_f_elements) {
            if (reference_children(
                    hierarchy.coarse_elements_to_reference(), element) > 1)
                marked.push_back(element);
        }
        if (marked.empty()) break;
        if (result.refinement_rounds >= maximum_rounds)
            throw std::runtime_error("hybrid matching refinement limit reached");
        const std::size_t before = hierarchy.coarse_mesh().elems.size();
        const ReferenceEpochRefinementResult proposed =
            hierarchy.propose_coarse_refinement(marked);
        if (!proposed.changed() || !hierarchy.reference_contains_proposed_coarse())
            throw std::runtime_error(
                "hybrid matching proposal is not contained in the fixed reference");
        const ReferenceEpochRefinementResult committed =
            hierarchy.commit_coarse_refinement();
        if (!committed.changed())
            throw std::runtime_error("hybrid matching coarse commit failed");
        result.refined_coarse_elements +=
            hierarchy.coarse_mesh().elems.size() - before;
        ++result.refinement_rounds;
    }
    result.regions = classify_singular_regions(
        hierarchy.coarse_mesh(), singular_points, ell);
    result.matching_holds = hybrid_reference_matching_holds(
        hierarchy, result.regions);
    result.reference_unchanged =
        hierarchy.reference_mesh_version() == reference_version
        && hierarchy.reference_mesh().elems.size() == reference_elements;
    if (!result.matching_holds || !result.reference_unchanged)
        throw std::runtime_error("hybrid matching invariant was not restored");
    return result;
}


HybridMatchingResult restore_hybrid_reference_matching_with_physical_radius(
    ReferenceEpochHierarchy &hierarchy,
    const std::vector<Point2> &singular_points,
    int ell,
    double minimum_physical_radius,
    std::size_t maximum_rounds) {
    if (maximum_rounds == 0 || hierarchy.has_proposed_coarse_refinement())
        throw std::invalid_argument(
            "hybrid matching requires a positive limit and no pending proposal");
    const std::uint64_t reference_version = hierarchy.reference_mesh_version();
    const std::size_t reference_elements = hierarchy.reference_mesh().elems.size();
    HybridMatchingResult result;
    for (;;) {
        result.regions = classify_singular_regions_with_physical_radius(
            hierarchy.coarse_mesh(), singular_points, ell,
            minimum_physical_radius);
        std::vector<int> marked;
        for (int element : result.regions.omega_f_elements) {
            if (reference_children(
                    hierarchy.coarse_elements_to_reference(), element) > 1)
                marked.push_back(element);
        }
        if (marked.empty()) break;
        if (result.refinement_rounds >= maximum_rounds)
            throw std::runtime_error("hybrid matching refinement limit reached");
        const std::size_t before = hierarchy.coarse_mesh().elems.size();
        const ReferenceEpochRefinementResult proposed =
            hierarchy.propose_coarse_refinement(marked);
        if (!proposed.changed() || !hierarchy.reference_contains_proposed_coarse())
            throw std::runtime_error(
                "hybrid matching proposal is not contained in the fixed reference");
        const ReferenceEpochRefinementResult committed =
            hierarchy.commit_coarse_refinement();
        if (!committed.changed())
            throw std::runtime_error("hybrid matching coarse commit failed");
        result.refined_coarse_elements +=
            hierarchy.coarse_mesh().elems.size() - before;
        ++result.refinement_rounds;
    }
    result.regions = classify_singular_regions_with_physical_radius(
        hierarchy.coarse_mesh(), singular_points, ell,
        minimum_physical_radius);
    result.matching_holds = hybrid_reference_matching_holds(
        hierarchy, result.regions);
    result.reference_unchanged =
        hierarchy.reference_mesh_version() == reference_version
        && hierarchy.reference_mesh().elems.size() == reference_elements;
    if (!result.matching_holds || !result.reference_unchanged)
        throw std::runtime_error("hybrid matching invariant was not restored");
    return result;
}

std::vector<int> mark_hybrid_regular_region(
    const std::vector<double> &element_eta_squared,
    const SingularRegionClassification &regions,
    double theta) {
    if (element_eta_squared.size() != regions.in_regular.size())
        throw std::invalid_argument(
            "hybrid marking indicators do not match the coarse mesh");
    std::vector<double> regular = element_eta_squared;
    for (int element = 0; element < static_cast<int>(regular.size()); ++element)
        if (!regions.in_regular[element]) regular[element] = 0.0;
    const double mass = std::accumulate(regular.begin(), regular.end(), 0.0);
    if (!(mass > 0.0)) return {};
    return mark_doerfler(regular, theta, regions.in_regular);
}

} // namespace lod2d::helmholtz::adaptive

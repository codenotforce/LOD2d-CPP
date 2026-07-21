#include "helmholtz/adaptive/hierarchy.h"
#include "helmholtz/model.h"
#include "helmholtz/operators.h"
#include "lod/patches.h"

#include <Eigen/Eigenvalues>
#include <algorithm>
#include <atomic>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace lod2d;
using namespace lod2d::helmholtz;
using namespace lod2d::helmholtz::adaptive;

namespace {

using Clock = std::chrono::steady_clock;

struct Options {
    double wavenumber = 2.0;
    int initial_level = 3;
    int fine_level = 12;
    int steps = 6;
    int ell = 3;
    std::string mark = "fraction";
    double mark_fraction = 0.25;
    Point2 seed = Point2(0.25, 0.25);
    double mass_threshold = 1e-12;
    int neighbor_layers = 0;
    std::string basis_selection = "all";
    std::string denominator_selection = "all";
    bool csv = false;
    bool check = false;
    bool only_final = false;
    std::string element_output;
    std::string mesh_output;
};

struct MeshDiagnostics {
    int level_min = 0;
    int level_max = 0;
    double H_min = 0.0;
    double H_max = 0.0;
    double grading_ratio = 1.0;
    double neighbor_ratio = 1.0;
    std::vector<double> diameters;
    std::vector<char> boundary;
    std::vector<char> transition;
};

struct NestingDiagnostics {
    double nodal_residual = 0.0;
    double element_residual = 0.0;
    double dg_residual = 0.0;
    double projection_residual = 0.0;
    double maximum() const {
        return std::max({nodal_residual, element_residual, dg_residual, projection_residual});
    }
};

struct LocalValue {
    double gradient = 0.0;
    double energy = 0.0;
    int mass_rank = 0;
    double mass_condition = 0.0;
    double stiffness_hermitian_defect = 0.0;
    double mass_hermitian_defect = 0.0;
    double eigen_residual = 0.0;
    double energy_identity_error = 0.0;
    int patch_elements = 0;
    double patch_diameter_over_H = 0.0;
};

struct Stats {
    double minimum = 0.0;
    double median = 0.0;
    double p90 = 0.0;
    double p99 = 0.0;
    double maximum = 0.0;
    int argmax = -1;
};

struct SummaryRow {
    int iteration = 0;
    std::string basis;
    std::string denominator;
    Stats stats;
    int argmax_level = 0;
    Point2 argmax_point = Point2::Zero();
    bool argmax_boundary = false;
    bool argmax_transition = false;
    LocalValue argmax_value;
    double global_H_scaled_max = 0.0;
};

struct RefinementCounts {
    int marked = 0;
    int closure_added = 0;
    int refined_total = 0;
};

struct EdgeHash {
    std::size_t operator()(const Edge &edge) const noexcept {
        return (static_cast<std::size_t>(static_cast<unsigned>(edge[0])) << 32U)
             ^ static_cast<std::size_t>(static_cast<unsigned>(edge[1]));
    }
};

Edge make_edge(int a, int b) {
    return a < b ? Edge{a, b} : Edge{b, a};
}

double elapsed_ms(const Clock::time_point &start) {
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

int parse_int(const std::string &text, const char *name) {
    std::size_t consumed = 0;
    const int value = std::stoi(text, &consumed);
    if (consumed != text.size())
        throw std::invalid_argument(std::string("invalid ") + name + ": " + text);
    return value;
}

double parse_double(const std::string &text, const char *name) {
    std::size_t consumed = 0;
    const double value = std::stod(text, &consumed);
    if (consumed != text.size() || !std::isfinite(value))
        throw std::invalid_argument(std::string("invalid ") + name + ": " + text);
    return value;
}

Options parse_options(int argc, char **argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        auto value = [&](const std::string &prefix) { return argument.substr(prefix.size()); };
        if (argument.rfind("--k=", 0) == 0)
            options.wavenumber = parse_double(value("--k="), "wavenumber");
        else if (argument.rfind("--initial-H=", 0) == 0)
            options.initial_level = parse_int(value("--initial-H="), "initial H");
        else if (argument.rfind("--fine-level=", 0) == 0)
            options.fine_level = parse_int(value("--fine-level="), "fine level");
        else if (argument.rfind("--steps=", 0) == 0)
            options.steps = parse_int(value("--steps="), "steps");
        else if (argument.rfind("--ell=", 0) == 0)
            options.ell = parse_int(value("--ell="), "ell");
        else if (argument.rfind("--mark=", 0) == 0)
            options.mark = value("--mark=");
        else if (argument.rfind("--mark-fraction=", 0) == 0)
            options.mark_fraction = parse_double(value("--mark-fraction="), "mark fraction");
        else if (argument.rfind("--seed-x=", 0) == 0)
            options.seed.x() = parse_double(value("--seed-x="), "seed x");
        else if (argument.rfind("--seed-y=", 0) == 0)
            options.seed.y() = parse_double(value("--seed-y="), "seed y");
        else if (argument.rfind("--mass-threshold=", 0) == 0)
            options.mass_threshold = parse_double(value("--mass-threshold="), "mass threshold");
        else if (argument.rfind("--neighbor-layers=", 0) == 0)
            options.neighbor_layers = parse_int(value("--neighbor-layers="), "neighbor layers");
        else if (argument.rfind("--basis=", 0) == 0)
            options.basis_selection = value("--basis=");
        else if (argument.rfind("--denominators=", 0) == 0)
            options.denominator_selection = value("--denominators=");
        else if (argument.rfind("--element-out=", 0) == 0)
            options.element_output = value("--element-out=");
        else if (argument.rfind("--mesh-out=", 0) == 0)
            options.mesh_output = value("--mesh-out=");
        else if (argument == "--format=csv") options.csv = true;
        else if (argument == "--check") options.check = true;
        else if (argument == "--only-final") options.only_final = true;
        else if (argument == "--help") {
            std::cout
                << "Usage: bench_helmholtz_local_inverse [--k=2] [--initial-H=3] "
                   "[--fine-level=12] [--steps=6] [--ell=3] "
                   "[--mark=fraction|single-chain|boundary-chain|argmax-element|argmax-patch] "
                   "[--mark-fraction=0.25] [--seed-x=0.25] [--seed-y=0.25] "
                   "[--mass-threshold=1e-12] [--neighbor-layers=0] [--element-out=PATH] "
                   "[--basis=all|trial] [--denominators=all|matched|element-matched] "
                   "[--mesh-out=PATH] [--format=csv] [--check] [--only-final]\n";
            std::exit(0);
        } else throw std::invalid_argument("unknown option: " + argument);
    }
    if (!(options.wavenumber > 0.0)) throw std::invalid_argument("k must be positive");
    if (options.initial_level < 0 || options.fine_level <= options.initial_level)
        throw std::invalid_argument("require 0 <= initial-H < fine-level");
    if (options.steps < 0 || options.ell < 0) throw std::invalid_argument("steps and ell must be nonnegative");
    if (!(options.mark_fraction > 0.0 && options.mark_fraction <= 1.0))
        throw std::invalid_argument("mark-fraction must lie in (0,1]");
    if (!(options.mass_threshold > 0.0 && options.mass_threshold < 1.0))
        throw std::invalid_argument("mass-threshold must lie in (0,1)");
    if (options.mark != "fraction" && options.mark != "single-chain"
        && options.mark != "boundary-chain" && options.mark != "argmax-element"
        && options.mark != "argmax-patch")
        throw std::invalid_argument(
            "mark must be fraction, single-chain, boundary-chain, argmax-element, or argmax-patch");
    if (options.neighbor_layers < 0)
        throw std::invalid_argument("neighbor-layers must be nonnegative");
    if (options.basis_selection != "all" && options.basis_selection != "trial")
        throw std::invalid_argument("basis must be all or trial");
    if (options.denominator_selection != "all"
        && options.denominator_selection != "matched"
        && options.denominator_selection != "element-matched")
        throw std::invalid_argument(
            "denominators must be all, matched, or element-matched");
    if (options.mark == "argmax-element" && options.denominator_selection == "matched")
        throw std::invalid_argument("argmax-element requires an element denominator");
    if (options.only_final
        && (options.mark == "argmax-element" || options.mark == "argmax-patch"))
        throw std::invalid_argument("only-final is incompatible with feedback argmax marking");
    return options;
}

double triangle_area(const TriMesh &mesh, int element) {
    const Triangle &tri = mesh.elems[element];
    const Point2 &a = mesh.nodes[tri[0]];
    const Point2 &b = mesh.nodes[tri[1]];
    const Point2 &c = mesh.nodes[tri[2]];
    return 0.5 * std::abs((b.x() - a.x()) * (c.y() - a.y())
                        - (b.y() - a.y()) * (c.x() - a.x()));
}

double element_diameter(const TriMesh &mesh, int element) {
    const Triangle &tri = mesh.elems[element];
    return std::max({
        (mesh.nodes[tri[0]] - mesh.nodes[tri[1]]).norm(),
        (mesh.nodes[tri[1]] - mesh.nodes[tri[2]]).norm(),
        (mesh.nodes[tri[2]] - mesh.nodes[tri[0]]).norm()});
}

Point2 centroid(const TriMesh &mesh, int element) {
    const Triangle &tri = mesh.elems[element];
    return (mesh.nodes[tri[0]] + mesh.nodes[tri[1]] + mesh.nodes[tri[2]]) / 3.0;
}

std::pair<Eigen::Matrix3d, Eigen::Matrix3d> local_stiffness_mass(
    const TriMesh &mesh,
    int element) {
    const Triangle &tri = mesh.elems[element];
    const Point2 &a = mesh.nodes[tri[0]];
    const Point2 &b = mesh.nodes[tri[1]];
    const Point2 &c = mesh.nodes[tri[2]];
    const double det = (b.x() - a.x()) * (c.y() - a.y())
                     - (b.y() - a.y()) * (c.x() - a.x());
    if (std::abs(det) <= 1e-15) throw std::runtime_error("degenerate fine triangle");
    const double area = 0.5 * std::abs(det);
    std::array<Point2, 3> gradients{
        Point2(b.y() - c.y(), c.x() - b.x()) / det,
        Point2(c.y() - a.y(), a.x() - c.x()) / det,
        Point2(a.y() - b.y(), b.x() - a.x()) / det};
    Eigen::Matrix3d stiffness;
    Eigen::Matrix3d mass;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            stiffness(i, j) = area * gradients[i].dot(gradients[j]);
            mass(i, j) = area * (i == j ? 2.0 : 1.0) / 12.0;
        }
    }
    return {stiffness, mass};
}

Eigen::SparseMatrix<double> cg_to_dg(const TriMesh &mesh) {
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(3 * mesh.elems.size());
    for (int t = 0; t < static_cast<int>(mesh.elems.size()); ++t)
        for (int i = 0; i < 3; ++i) triplets.emplace_back(3 * t + i, mesh.elems[t][i], 1.0);
    Eigen::SparseMatrix<double> result(
        3 * static_cast<int>(mesh.elems.size()), static_cast<int>(mesh.nodes.size()));
    result.setFromTriplets(triplets.begin(), triplets.end());
    return result;
}

std::vector<std::vector<int>> fine_children(const HelmholtzProblemData &problem) {
    std::vector<std::vector<int>> children(problem.coarse.elems.size());
    std::vector<int> parent_count(problem.fine.elems.size(), 0);
    for (int coarse = 0; coarse < problem.fine_element_prolongation.outerSize(); ++coarse) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(
                 problem.fine_element_prolongation, coarse); it; ++it) {
            if (std::abs(it.value()) <= 1e-14) continue;
            if (std::abs(it.value() - 1.0) > 1e-12)
                throw std::runtime_error("fine element prolongation is not a parent map");
            children[coarse].push_back(it.row());
            ++parent_count[it.row()];
        }
    }
    if (std::find_if(parent_count.begin(), parent_count.end(),
                     [](int count) { return count != 1; }) != parent_count.end())
        throw std::runtime_error("a fine element does not have exactly one coarse parent");
    return children;
}

NestingDiagnostics verify_nesting(const HelmholtzProblemData &problem) {
    NestingDiagnostics diagnostics;
    const Eigen::SparseMatrix<double> &P = problem.coarse_to_fine;
    Eigen::VectorXd coarse_one = Eigen::VectorXd::Ones(problem.coarse.nodes.size());
    Eigen::VectorXd coarse_x(problem.coarse.nodes.size());
    Eigen::VectorXd coarse_y(problem.coarse.nodes.size());
    Eigen::VectorXd fine_x(problem.fine.nodes.size());
    Eigen::VectorXd fine_y(problem.fine.nodes.size());
    for (int i = 0; i < coarse_x.size(); ++i) {
        coarse_x(i) = problem.coarse.nodes[i].x();
        coarse_y(i) = problem.coarse.nodes[i].y();
    }
    for (int i = 0; i < fine_x.size(); ++i) {
        fine_x(i) = problem.fine.nodes[i].x();
        fine_y(i) = problem.fine.nodes[i].y();
    }
    diagnostics.nodal_residual = std::max({
        (P * coarse_one - Eigen::VectorXd::Ones(problem.fine.nodes.size())).lpNorm<Eigen::Infinity>(),
        (P * coarse_x - fine_x).lpNorm<Eigen::Infinity>(),
        (P * coarse_y - fine_y).lpNorm<Eigen::Infinity>()});

    std::vector<double> accumulated(problem.coarse.elems.size(), 0.0);
    std::vector<int> parents(problem.fine.elems.size(), 0);
    for (int coarse = 0; coarse < problem.fine_element_prolongation.outerSize(); ++coarse) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(
                 problem.fine_element_prolongation, coarse); it; ++it) {
            if (std::abs(it.value()) <= 1e-14) continue;
            accumulated[coarse] += it.value() * triangle_area(problem.fine, it.row());
            ++parents[it.row()];
        }
    }
    for (int fine_parent_count : parents)
        diagnostics.element_residual = std::max(
            diagnostics.element_residual, std::abs(fine_parent_count - 1.0));
    for (int coarse = 0; coarse < static_cast<int>(accumulated.size()); ++coarse) {
        const double area = triangle_area(problem.coarse, coarse);
        diagnostics.element_residual = std::max(
            diagnostics.element_residual, std::abs(accumulated[coarse] - area) / area);
    }

    const Eigen::SparseMatrix<double> coarse_dg = cg_to_dg(problem.coarse);
    const Eigen::SparseMatrix<double> fine_dg = cg_to_dg(problem.fine);
    diagnostics.dg_residual = (problem.fine_dg_prolongation * coarse_dg
                             - fine_dg * problem.coarse_to_fine).norm();
    const Eigen::MatrixXd projection(problem.quasi_interpolation * problem.coarse_to_fine);
    diagnostics.projection_residual =
        (projection - Eigen::MatrixXd::Identity(projection.rows(), projection.cols())).norm();
    return diagnostics;
}

using CanonicalTriangle = std::array<long long, 6>;

std::vector<CanonicalTriangle> canonical_mesh(const TriMesh &mesh, int fine_level) {
    const double scale = std::ldexp(1.0, fine_level + 2);
    std::vector<CanonicalTriangle> result;
    result.reserve(mesh.elems.size());
    for (const Triangle &tri : mesh.elems) {
        std::array<std::pair<long long, long long>, 3> points;
        for (int i = 0; i < 3; ++i) {
            const Point2 &point = mesh.nodes[tri[i]];
            points[i] = {std::llround(scale * point.x()), std::llround(scale * point.y())};
        }
        std::sort(points.begin(), points.end());
        result.push_back({points[0].first, points[0].second,
                          points[1].first, points[1].second,
                          points[2].first, points[2].second});
    }
    std::sort(result.begin(), result.end());
    return result;
}

MeshDiagnostics mesh_diagnostics(
    const TriMesh &mesh,
    const std::vector<int> &levels) {
    MeshDiagnostics result;
    result.level_min = *std::min_element(levels.begin(), levels.end());
    result.level_max = *std::max_element(levels.begin(), levels.end());
    result.diameters.resize(mesh.elems.size());
    result.boundary.assign(mesh.elems.size(), false);
    result.transition.assign(mesh.elems.size(), false);
    result.H_min = std::numeric_limits<double>::infinity();
    result.H_max = 0.0;
    std::unordered_map<Edge, std::vector<int>, EdgeHash> owners;
    owners.reserve(3 * mesh.elems.size());
    for (int t = 0; t < static_cast<int>(mesh.elems.size()); ++t) {
        result.diameters[t] = element_diameter(mesh, t);
        result.H_min = std::min(result.H_min, result.diameters[t]);
        result.H_max = std::max(result.H_max, result.diameters[t]);
        const Triangle &tri = mesh.elems[t];
        owners[make_edge(tri[0], tri[1])].push_back(t);
        owners[make_edge(tri[1], tri[2])].push_back(t);
        owners[make_edge(tri[2], tri[0])].push_back(t);
    }
    result.grading_ratio = result.H_max / result.H_min;
    for (const auto &[edge, adjacent] : owners) {
        (void)edge;
        if (adjacent.size() == 1) result.boundary[adjacent[0]] = true;
        else if (adjacent.size() == 2) {
            const int a = adjacent[0];
            const int b = adjacent[1];
            const double ratio = std::max(
                result.diameters[a] / result.diameters[b],
                result.diameters[b] / result.diameters[a]);
            result.neighbor_ratio = std::max(result.neighbor_ratio, ratio);
            if (levels[a] != levels[b]) {
                result.transition[a] = true;
                result.transition[b] = true;
            }
        } else throw std::runtime_error("nonmanifold coarse edge");
    }
    return result;
}

std::vector<std::vector<int>> patches_for_layers(const TriMesh &coarse, int layers) {
    if (layers < 0) throw std::invalid_argument("patch layers must be nonnegative");
    std::vector<std::vector<int>> result(coarse.elems.size());
    if (layers == 0) {
        for (int target = 0; target < static_cast<int>(result.size()); ++target)
            result[target].push_back(target);
        return result;
    }
    const Eigen::SparseMatrix<double> matrix = build_patches(coarse, layers);
    for (int target = 0; target < matrix.outerSize(); ++target) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(matrix, target); it; ++it)
            if (it.value() != 0.0) result[target].push_back(it.row());
    }
    return result;
}

double patch_diameter(
    const TriMesh &mesh,
    const std::vector<int> &elements) {
    std::set<int> vertices;
    for (int element : elements)
        for (int vertex : mesh.elems[element]) vertices.insert(vertex);
    double diameter = 0.0;
    for (auto a = vertices.begin(); a != vertices.end(); ++a)
        for (auto b = std::next(a); b != vertices.end(); ++b)
            diameter = std::max(diameter, (mesh.nodes[*a] - mesh.nodes[*b]).norm());
    return diameter;
}

using BasisRows = std::vector<std::vector<std::pair<int, Complex>>>;

BasisRows basis_rows(const ComplexSparseMatrix &basis) {
    BasisRows rows(basis.rows());
    for (int column = 0; column < basis.outerSize(); ++column)
        for (ComplexSparseMatrix::InnerIterator it(basis, column); it; ++it)
            if (std::abs(it.value()) > 0.0) rows[it.row()].push_back({column, it.value()});
    return rows;
}

LocalValue local_value(
    const TriMesh &coarse,
    const TriMesh &fine,
    int coarse_element,
    const std::vector<std::vector<int>> &children,
    const std::vector<std::vector<int>> &patches,
    const BasisRows &rows,
    int basis_columns,
    bool patch_denominator,
    double wavenumber,
    double mass_threshold) {
    const std::vector<int> &denominator_coarse = patch_denominator
        ? patches[coarse_element]
        : std::vector<int>{coarse_element};
    std::vector<int> denominator_fine;
    for (int element : denominator_coarse)
        denominator_fine.insert(
            denominator_fine.end(), children[element].begin(), children[element].end());

    std::vector<char> active_marker(basis_columns, false);
    std::vector<int> active;
    for (int fine_element : denominator_fine) {
        for (int vertex : fine.elems[fine_element]) {
            for (const auto &[column, value] : rows[vertex]) {
                (void)value;
                if (!active_marker[column]) {
                    active_marker[column] = true;
                    active.push_back(column);
                }
            }
        }
    }
    if (active.empty()) return {};
    std::sort(active.begin(), active.end());
    std::vector<int> local_column(basis_columns, -1);
    for (int i = 0; i < static_cast<int>(active.size()); ++i) local_column[active[i]] = i;
    const int n = static_cast<int>(active.size());
    ComplexMatrix stiffness = ComplexMatrix::Zero(n, n);
    ComplexMatrix mass_element = ComplexMatrix::Zero(n, n);
    ComplexMatrix mass_denominator = ComplexMatrix::Zero(n, n);

    auto element_basis = [&](int fine_element) {
        ComplexMatrix B = ComplexMatrix::Zero(3, n);
        for (int i = 0; i < 3; ++i) {
            const int vertex = fine.elems[fine_element][i];
            for (const auto &[column, value] : rows[vertex]) {
                const int local = local_column[column];
                if (local >= 0) B(i, local) = value;
            }
        }
        return B;
    };
    for (int fine_element : children[coarse_element]) {
        const auto [Se, Me] = local_stiffness_mass(fine, fine_element);
        const ComplexMatrix B = element_basis(fine_element);
        stiffness.noalias() += B.adjoint() * Se.cast<Complex>() * B;
        mass_element.noalias() += B.adjoint() * Me.cast<Complex>() * B;
    }
    for (int fine_element : denominator_fine) {
        const auto [Se, Me] = local_stiffness_mass(fine, fine_element);
        (void)Se;
        const ComplexMatrix B = element_basis(fine_element);
        mass_denominator.noalias() += B.adjoint() * Me.cast<Complex>() * B;
    }

    LocalValue result;
    result.patch_elements = static_cast<int>(denominator_coarse.size());
    result.patch_diameter_over_H = patch_diameter(coarse, denominator_coarse)
                                 / element_diameter(coarse, coarse_element);
    result.stiffness_hermitian_defect = (stiffness - stiffness.adjoint()).norm()
        / std::max(stiffness.norm(), 1e-30);
    result.mass_hermitian_defect = (mass_denominator - mass_denominator.adjoint()).norm()
        / std::max(mass_denominator.norm(), 1e-30);
    stiffness = 0.5 * (stiffness + stiffness.adjoint());
    mass_element = 0.5 * (mass_element + mass_element.adjoint());
    mass_denominator = 0.5 * (mass_denominator + mass_denominator.adjoint());

    Eigen::SelfAdjointEigenSolver<ComplexMatrix> mass_solver(mass_denominator);
    if (mass_solver.info() != Eigen::Success) throw std::runtime_error("local mass eigensolve failed");
    const Eigen::VectorXd mass_values = mass_solver.eigenvalues();
    const double mass_max = mass_values.size() ? mass_values.maxCoeff() : 0.0;
    if (!(mass_max > 0.0)) return result;
    const double tolerance = mass_threshold * mass_max;
    for (int i = 0; i < mass_values.size(); ++i)
        if (mass_values(i) > tolerance) ++result.mass_rank;
    if (result.mass_rank == 0) return result;
    ComplexMatrix Z(n, result.mass_rank);
    double mass_min = mass_max;
    int retained = 0;
    for (int i = 0; i < mass_values.size(); ++i) {
        if (mass_values(i) <= tolerance) continue;
        mass_min = std::min(mass_min, mass_values(i));
        Z.col(retained++) = mass_solver.eigenvectors().col(i) / std::sqrt(mass_values(i));
    }
    result.mass_condition = mass_max / mass_min;

    auto largest = [&](const ComplexMatrix &numerator) {
        ComplexMatrix reduced = Z.adjoint() * numerator * Z;
        reduced = 0.5 * (reduced + reduced.adjoint());
        Eigen::SelfAdjointEigenSolver<ComplexMatrix> solver(reduced);
        if (solver.info() != Eigen::Success) throw std::runtime_error("local stiffness eigensolve failed");
        const int index = solver.eigenvalues().size() - 1;
        const double lambda = std::max(0.0, solver.eigenvalues()(index));
        const ComplexVector vector = solver.eigenvectors().col(index);
        const double residual = (reduced * vector - lambda * vector).norm()
                              / std::max(reduced.norm() * vector.norm(), 1e-30);
        return std::pair<double, double>{lambda, residual};
    };
    const auto [lambda_gradient, residual_gradient] = largest(stiffness);
    const ComplexMatrix energy = stiffness + wavenumber * wavenumber * mass_element;
    const auto [lambda_energy, residual_energy] = largest(energy);
    const double H = element_diameter(coarse, coarse_element);
    result.gradient = H * std::sqrt(lambda_gradient);
    result.energy = H * std::sqrt(lambda_energy);
    result.eigen_residual = std::max(residual_gradient, residual_energy);
    if (!patch_denominator) {
        const double identity_scale = std::max({
            result.energy * result.energy,
            result.gradient * result.gradient + wavenumber * wavenumber * H * H,
            1.0});
        result.energy_identity_error = std::abs(
            result.energy * result.energy - result.gradient * result.gradient
            - wavenumber * wavenumber * H * H) / identity_scale;
    }
    return result;
}

Stats summarize(const std::vector<LocalValue> &values) {
    std::vector<int> order(values.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        if (values[a].gradient != values[b].gradient)
            return values[a].gradient < values[b].gradient;
        return a < b;
    });
    auto quantile = [&](double q) {
        const int index = std::clamp(
            static_cast<int>(std::ceil(q * order.size())) - 1,
            0, static_cast<int>(order.size()) - 1);
        return values[order[index]].gradient;
    };
    Stats stats;
    stats.minimum = values[order.front()].gradient;
    stats.median = quantile(0.50);
    stats.p90 = quantile(0.90);
    stats.p99 = quantile(0.99);
    stats.maximum = values[order.back()].gradient;
    stats.argmax = order.back();
    return stats;
}

std::vector<int> choose_marked(
    const Options &options,
    const AdaptiveMeshHierarchy &hierarchy) {
    const auto &levels = hierarchy.coarse_levels();
    const int maximum = *std::max_element(levels.begin(), levels.end());
    Point2 seed = options.seed;
    if (options.mark == "boundary-chain") seed = Point2(0.0, 0.5);
    std::vector<int> candidates;
    for (int t = 0; t < static_cast<int>(levels.size()); ++t)
        if (levels[t] == maximum) candidates.push_back(t);
    std::sort(candidates.begin(), candidates.end(), [&](int a, int b) {
        const double da = (centroid(hierarchy.coarse_mesh(), a) - seed).squaredNorm();
        const double db = (centroid(hierarchy.coarse_mesh(), b) - seed).squaredNorm();
        if (std::abs(da - db) > 1e-15) return da < db;
        return hierarchy.coarse_element_ids()[a] < hierarchy.coarse_element_ids()[b];
    });
    int count = 1;
    if (options.mark == "fraction")
        count = std::max(1, static_cast<int>(std::ceil(options.mark_fraction * candidates.size())));
    candidates.resize(std::min(count, static_cast<int>(candidates.size())));
    return candidates;
}

std::vector<int> choose_argmax_marked(
    const AdaptiveMeshHierarchy &hierarchy,
    int target,
    int neighbor_layers) {
    if (target < 0 || target >= static_cast<int>(hierarchy.coarse_mesh().elems.size()))
        throw std::runtime_error("argmax feedback target is invalid");
    std::vector<int> marked;
    if (neighbor_layers == 0) {
        marked.push_back(target);
    } else {
        const Eigen::SparseMatrix<double> patches = build_patches(
            hierarchy.coarse_mesh(), neighbor_layers);
        for (Eigen::SparseMatrix<double>::InnerIterator it(patches, target); it; ++it)
            if (it.value() != 0.0) marked.push_back(it.row());
    }
    std::sort(marked.begin(), marked.end());
    marked.erase(std::unique(marked.begin(), marked.end()), marked.end());
    marked.erase(std::remove_if(marked.begin(), marked.end(), [&](int element) {
        return hierarchy.coarse_levels()[element] >= hierarchy.fine_level();
    }), marked.end());
    if (marked.empty())
        throw std::runtime_error("argmax feedback patch has no refinable coarse element");
    return marked;
}

RefinementCounts refine_hierarchy(
    AdaptiveMeshHierarchy &hierarchy,
    const std::vector<int> &marked) {
    const std::vector<std::uint64_t> old_ids = hierarchy.coarse_element_ids();
    std::set<std::uint64_t> old_set(old_ids.begin(), old_ids.end());
    hierarchy.refine(marked);
    std::set<std::uint64_t> new_set(
        hierarchy.coarse_element_ids().begin(), hierarchy.coarse_element_ids().end());
    int preserved = 0;
    for (std::uint64_t id : old_set) if (new_set.count(id)) ++preserved;
    RefinementCounts result;
    result.marked = static_cast<int>(marked.size());
    result.refined_total = static_cast<int>(old_set.size()) - preserved;
    result.closure_added = result.refined_total - result.marked;
    if (result.closure_added < 0) throw std::runtime_error("invalid NVB closure count");
    return result;
}

void write_mesh_header(std::ofstream &output) {
    output << "iteration,element,id,parent_id,level,n0,n1,n2,x0,y0,x1,y1,x2,y2\n";
}

void write_mesh(
    std::ofstream &output,
    int iteration,
    const AdaptiveMeshHierarchy &hierarchy) {
    if (!output) return;
    const TriMesh &mesh = hierarchy.coarse_mesh();
    for (int t = 0; t < static_cast<int>(mesh.elems.size()); ++t) {
        const Triangle &tri = mesh.elems[t];
        output << iteration << ',' << t << ',' << hierarchy.coarse_element_ids()[t] << ','
               << hierarchy.coarse_parent_ids()[t] << ',' << hierarchy.coarse_levels()[t] << ','
               << tri[0] << ',' << tri[1] << ',' << tri[2];
        for (int vertex : tri)
            output << ',' << mesh.nodes[vertex].x() << ',' << mesh.nodes[vertex].y();
        output << '\n';
    }
}

void write_element_header(std::ofstream &output) {
    output << "iteration,basis,denominator,element,id,level,H,boundary,transition,"
              "Q_gradient,Q_energy,mass_rank,mass_condition,stiffness_hermitian_defect,"
              "mass_hermitian_defect,eigen_residual,energy_identity_error,"
              "patch_elements,patch_diameter_over_H,x,y\n";
}

void print_summary_header() {
    std::cout
        << "iteration,k,ell,mark,neighbor_layers,basis_selection,denominator_selection,"
           "basis,denominator,coarse_nodes,coarse_elements,"
           "fine_nodes,fine_elements,L_min,L_max,level_gap,H_min,H_max,grading_ratio,"
           "neighbor_ratio,q_min,q_max,marked,closure_added,refined_total,nesting_residual,"
           "fixed_fine_mesh,petrov_residual,primal_residual,adjoint_residual,"
           "constraint_residual,trial_test_conjugacy,min,median,p90,p99,max,global_H_scaled_max,argmax,"
           "argmax_level,argmax_x,argmax_y,argmax_boundary,argmax_transition,"
           "argmax_energy,mass_rank,mass_condition,hermitian_defect,eigen_residual,"
           "energy_identity_error,patch_elements,patch_diameter_over_H,build_ms,inverse_ms\n";
}

void print_summary_row(
    const Options &options,
    const SummaryRow &row,
    const HelmholtzLodModel &model,
    const MeshDiagnostics &mesh,
    const NestingDiagnostics &nesting,
    const RefinementCounts &counts,
    double q_min,
    double q_max,
    bool fixed_fine_mesh,
    double petrov_residual,
    double conjugacy,
    double inverse_ms) {
    const auto &diagnostics = model.correctors().diagnostics;
    std::cout << std::setprecision(17)
              << row.iteration << ',' << options.wavenumber << ',' << options.ell << ','
              << options.mark << ',' << options.neighbor_layers << ','
              << options.basis_selection << ',' << options.denominator_selection << ','
              << row.basis << ',' << row.denominator << ','
              << model.problem().coarse.nodes.size() << ',' << model.problem().coarse.elems.size() << ','
              << model.problem().fine.nodes.size() << ',' << model.problem().fine.elems.size() << ','
              << mesh.level_min << ',' << mesh.level_max << ','
              << mesh.level_max - mesh.level_min << ',' << mesh.H_min << ',' << mesh.H_max << ','
              << mesh.grading_ratio << ',' << mesh.neighbor_ratio << ','
              << q_min << ',' << q_max << ','
              << counts.marked << ',' << counts.closure_added << ',' << counts.refined_total << ','
              << nesting.maximum() << ',' << (fixed_fine_mesh ? 1 : 0) << ','
              << petrov_residual << ',' << diagnostics.max_primal_residual << ','
              << diagnostics.max_adjoint_residual << ',' << diagnostics.max_constraint_residual << ','
              << conjugacy << ',' << row.stats.minimum << ',' << row.stats.median << ','
              << row.stats.p90 << ',' << row.stats.p99 << ',' << row.stats.maximum << ','
              << row.global_H_scaled_max << ','
              << row.stats.argmax << ',' << row.argmax_level << ',' << row.argmax_point.x() << ','
              << row.argmax_point.y() << ',' << row.argmax_boundary << ',' << row.argmax_transition << ','
              << row.argmax_value.energy << ',' << row.argmax_value.mass_rank << ','
              << row.argmax_value.mass_condition << ','
              << std::max(row.argmax_value.stiffness_hermitian_defect,
                          row.argmax_value.mass_hermitian_defect) << ','
              << row.argmax_value.eigen_residual << ',' << row.argmax_value.energy_identity_error << ','
              << row.argmax_value.patch_elements << ',' << row.argmax_value.patch_diameter_over_H << ','
              << model.build_timings().total_ms << ',' << inverse_ms << '\n';
}

} // namespace

int main(int argc, char **argv) {
    try {
        const Options options = parse_options(argc, argv);
        AdaptiveMeshHierarchy hierarchy(
            make_helmholtz_unit_square_mesh(), options.initial_level, options.fine_level);
        const TriMesh reference_fine = refine_mesh_nvb(
            make_helmholtz_unit_square_mesh(), options.fine_level).mesh;
        const auto reference_signature = canonical_mesh(reference_fine, options.fine_level);

        std::ofstream element_output;
        if (!options.element_output.empty()) {
            element_output.open(options.element_output);
            if (!element_output) throw std::runtime_error("cannot open element output");
            element_output << std::setprecision(17);
            write_element_header(element_output);
        }
        std::ofstream mesh_output;
        if (!options.mesh_output.empty()) {
            mesh_output.open(options.mesh_output);
            if (!mesh_output) throw std::runtime_error("cannot open mesh output");
            mesh_output << std::setprecision(17);
            write_mesh_header(mesh_output);
        }
        if (options.csv) print_summary_header();

        RefinementCounts incoming_counts;
        const double initial_grading = mesh_diagnostics(
            hierarchy.coarse_mesh(), hierarchy.coarse_levels()).grading_ratio;
        double final_grading = 1.0;
        std::map<std::pair<std::string, std::string>, std::vector<double>> maxima;
        for (int iteration = 0; iteration <= options.steps; ++iteration) {
            if (options.only_final && iteration < options.steps) {
                const std::vector<int> marked = choose_marked(options, hierarchy);
                const int maximum_level = *std::max_element(
                    hierarchy.coarse_levels().begin(), hierarchy.coarse_levels().end());
                for (int element : marked) {
                    if (hierarchy.coarse_levels()[element] != maximum_level)
                        throw std::runtime_error("marking selected a non-finest element");
                    if (hierarchy.coarse_levels()[element] >= options.fine_level)
                        throw std::runtime_error("coarse refinement reached the fixed fine level");
                }
                incoming_counts = refine_hierarchy(hierarchy, marked);
                continue;
            }
            HelmholtzProblemConfig config;
            config.h = options.fine_level;
            config.ell = options.ell;
            config.wavenumber = options.wavenumber;
            config.mode = HelmholtzPetrovMode::TwoSided;
            HelmholtzLodModel model = HelmholtzLodModel::build_adaptive(
                config, hierarchy.coarse_mesh(), hierarchy.coarse_levels());
            const NestingDiagnostics nesting = verify_nesting(model.problem());
            const bool fixed_fine_mesh = canonical_mesh(
                model.problem().fine, options.fine_level) == reference_signature;
            const MeshDiagnostics mesh = mesh_diagnostics(
                hierarchy.coarse_mesh(), hierarchy.coarse_levels());
            final_grading = mesh.grading_ratio;
            if (options.check && (nesting.maximum() > 1e-10 || !fixed_fine_mesh))
                throw std::runtime_error("V_H is not exactly nested in the fixed V_h");

            const HelmholtzLodSolution calibration = model.solve_source(
                [](const Point2 &point) { return Complex(1.0 + point.x(), -point.y()); });
            ComplexSparseMatrix conjugate_trial = model.corrected_trial_basis().conjugate();
            const double conjugacy = (model.corrected_test_basis() - conjugate_trial).norm()
                                   / std::max(model.corrected_test_basis().norm(), 1e-30);
            if (options.check && (calibration.petrov_residual > 1e-8 || conjugacy > 1e-10))
                throw std::runtime_error("Helmholtz basis or Petrov residual check failed");

            const auto children = fine_children(model.problem());
            const auto element_denominators = patches_for_layers(model.problem().coarse, 0);
            const auto patch1_denominators = patches_for_layers(model.problem().coarse, 1);
            const auto patch_ell_denominators = patches_for_layers(
                model.problem().coarse, options.ell);
            double q_min = std::numeric_limits<double>::infinity();
            double q_max = 0.0;
            for (int T = 0; T < static_cast<int>(children.size()); ++T) {
                double child_h = 0.0;
                for (int fine_element : children[T])
                    child_h = std::max(child_h, element_diameter(model.problem().fine, fine_element));
                const double q = child_h / mesh.diameters[T];
                q_min = std::min(q_min, q);
                q_max = std::max(q_max, q);
            }
            ComplexSparseMatrix coarse_basis;
            std::vector<std::pair<std::string, const ComplexSparseMatrix *>> bases;
            if (options.basis_selection == "all") {
                coarse_basis = model.problem().coarse_to_fine.cast<Complex>();
                bases.push_back({"coarse", &coarse_basis});
            }
            bases.push_back({"trial", &model.corrected_trial_basis()});
            if (options.basis_selection == "all")
                bases.push_back({"test", &model.corrected_test_basis()});
            const auto inverse_start = Clock::now();
            int feedback_element_argmax = -1;
            int feedback_patch_argmax = -1;
            for (const auto &[basis_name, basis_pointer] : bases) {
                const ComplexSparseMatrix &basis = *basis_pointer;
                const BasisRows rows = basis_rows(basis);
                std::vector<std::pair<std::string, const std::vector<std::vector<int>> *>> denominators;
                const std::string ell_denominator = "patch" + std::to_string(options.ell);
                if (options.denominator_selection == "all"
                    || options.denominator_selection == "element-matched")
                    denominators.push_back({"element", &element_denominators});
                if (options.denominator_selection == "all")
                    denominators.push_back({"patch1", &patch1_denominators});
                if (options.ell != 1 || options.denominator_selection != "all")
                    denominators.push_back({ell_denominator, &patch_ell_denominators});
                for (const auto &[denominator, denominator_patches] : denominators) {
                    const bool patch_denominator = denominator != "element";
                    std::vector<LocalValue> values(model.problem().coarse.elems.size());
                    double maximum_defect = 0.0;
                    double maximum_eigen_residual = 0.0;
                    double maximum_identity_error = 0.0;
                    std::atomic<bool> local_failed{false};
                    std::exception_ptr local_exception;
                    std::mutex local_exception_mutex;
#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic, 1)
#endif
                    for (int T = 0; T < static_cast<int>(values.size()); ++T) {
                        if (local_failed.load(std::memory_order_relaxed)) continue;
                        try {
                            values[T] = local_value(
                                model.problem().coarse, model.problem().fine, T,
                                children, *denominator_patches, rows, basis.cols(), patch_denominator,
                                options.wavenumber, options.mass_threshold);
                        } catch (...) {
                            local_failed.store(true, std::memory_order_relaxed);
                            std::lock_guard<std::mutex> lock(local_exception_mutex);
                            if (!local_exception) local_exception = std::current_exception();
                        }
                    }
                    if (local_exception) std::rethrow_exception(local_exception);
                    for (int T = 0; T < static_cast<int>(values.size()); ++T) {
                        maximum_defect = std::max(maximum_defect, std::max(
                            values[T].stiffness_hermitian_defect,
                            values[T].mass_hermitian_defect));
                        maximum_eigen_residual = std::max(
                            maximum_eigen_residual, values[T].eigen_residual);
                        maximum_identity_error = std::max(
                            maximum_identity_error, values[T].energy_identity_error);
                        if (element_output) {
                            const Point2 center = centroid(model.problem().coarse, T);
                            element_output
                                << iteration << ',' << basis_name << ',' << denominator << ',' << T << ','
                                << hierarchy.coarse_element_ids()[T] << ',' << hierarchy.coarse_levels()[T] << ','
                                << mesh.diameters[T] << ',' << static_cast<int>(mesh.boundary[T]) << ','
                                << static_cast<int>(mesh.transition[T]) << ',' << values[T].gradient << ','
                                << values[T].energy << ',' << values[T].mass_rank << ','
                                << values[T].mass_condition << ',' << values[T].stiffness_hermitian_defect << ','
                                << values[T].mass_hermitian_defect << ',' << values[T].eigen_residual << ','
                                << values[T].energy_identity_error << ',' << values[T].patch_elements << ','
                                << values[T].patch_diameter_over_H << ',' << center.x() << ',' << center.y() << '\n';
                        }
                    }
                    if (options.check && (maximum_defect > 1e-10
                                          || maximum_eigen_residual > 1e-8
                                          || (!patch_denominator && maximum_identity_error > 1e-8)))
                        throw std::runtime_error("local Hermitian eigenproblem check failed");
                    SummaryRow row;
                    row.iteration = iteration;
                    row.basis = basis_name;
                    row.denominator = denominator;
                    row.stats = summarize(values);
                    row.argmax_level = hierarchy.coarse_levels()[row.stats.argmax];
                    row.argmax_point = centroid(model.problem().coarse, row.stats.argmax);
                    row.argmax_boundary = mesh.boundary[row.stats.argmax];
                    row.argmax_transition = mesh.transition[row.stats.argmax];
                    row.argmax_value = values[row.stats.argmax];
                    if (basis_name == "trial" && denominator == "element")
                        feedback_element_argmax = row.stats.argmax;
                    if (basis_name == "trial" && denominator == ell_denominator)
                        feedback_patch_argmax = row.stats.argmax;
                    for (int T = 0; T < static_cast<int>(values.size()); ++T)
                        row.global_H_scaled_max = std::max(
                            row.global_H_scaled_max,
                            values[T].gradient * mesh.H_max / mesh.diameters[T]);
                    maxima[{basis_name, denominator}].push_back(row.stats.maximum);
                    if (options.csv) {
                        print_summary_row(
                            options, row, model, mesh, nesting, incoming_counts,
                            q_min, q_max,
                            fixed_fine_mesh, calibration.petrov_residual, conjugacy,
                            elapsed_ms(inverse_start));
                    } else {
                        std::cout << "iteration=" << iteration << " basis=" << basis_name
                                  << " denominator=" << denominator << " levels="
                                  << mesh.level_min << ':' << mesh.level_max
                                  << " grading=" << mesh.grading_ratio
                                  << " maxQ=" << row.stats.maximum
                                  << " nesting=" << nesting.maximum() << '\n';
                    }
                }
            }
            write_mesh(mesh_output, iteration, hierarchy);

            if (iteration == options.steps) break;
            std::vector<int> marked;
            const bool feedback = options.mark == "argmax-element" || options.mark == "argmax-patch";
            if (options.mark == "argmax-element") {
                marked = choose_argmax_marked(
                    hierarchy, feedback_element_argmax, options.neighbor_layers);
            } else if (options.mark == "argmax-patch") {
                marked = choose_argmax_marked(
                    hierarchy, feedback_patch_argmax, options.neighbor_layers);
            } else {
                marked = choose_marked(options, hierarchy);
            }
            for (int element : marked) {
                if (!feedback && hierarchy.coarse_levels()[element] != mesh.level_max)
                    throw std::runtime_error("marking selected a non-finest element");
                if (hierarchy.coarse_levels()[element] >= options.fine_level)
                    throw std::runtime_error("coarse refinement reached the fixed fine level");
            }
            incoming_counts = refine_hierarchy(hierarchy, marked);
        }

        if (options.check) {
            if (options.steps > 0 && !(final_grading > initial_grading))
                throw std::runtime_error("non-quasi-uniform grading did not increase");
            for (const auto &[key, values] : maxima) {
                (void)key;
                for (double value : values)
                    if (!(std::isfinite(value) && value > 0.0))
                        throw std::runtime_error("local inverse constant is not finite and positive");
            }
            const auto coarse_element = maxima.find({"coarse", "element"});
            if (coarse_element != maxima.end()) {
                const double coarse_min = *std::min_element(
                    coarse_element->second.begin(), coarse_element->second.end());
                const double coarse_max = *std::max_element(
                    coarse_element->second.begin(), coarse_element->second.end());
                if (coarse_max / coarse_min > 1.30)
                    throw std::runtime_error("C=0 local-H baseline is not stable");
            }
        }
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "bench_helmholtz_local_inverse failed: " << error.what() << '\n';
        return 1;
    }
}

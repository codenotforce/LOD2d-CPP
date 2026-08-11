#include "helmholtz/adaptive/reference_retraction.h"

#include "helmholtz/boundary.h"

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include <Eigen/SVD>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace lod2d::helmholtz::adaptive {
namespace {

constexpr double kRankTolerance = 2e-12;

struct PointKey {
    std::uint64_t x = 0;
    std::uint64_t y = 0;

    bool operator==(const PointKey &) const = default;
};

struct PointKeyHash {
    std::size_t operator()(const PointKey &key) const noexcept {
        const std::uint64_t mixed = key.x
            ^ (key.y + 0x9e3779b97f4a7c15ULL + (key.x << 6U)
               + (key.x >> 2U));
        return static_cast<std::size_t>(mixed);
    }
};

PointKey point_key(const Point2 &point) {
    double x = point.x();
    double y = point.y();
    if (x == 0.0) x = 0.0;
    if (y == 0.0) y = 0.0;
    return {std::bit_cast<std::uint64_t>(x),
            std::bit_cast<std::uint64_t>(y)};
}

Eigen::SparseMatrix<double> energy_matrix(
    const HelmholtzOperators &operators) {
    Eigen::SparseMatrix<double> result = operators.stiffness;
    result += operators.wavenumber * operators.wavenumber * operators.mass;
    result.makeCompressed();
    return result;
}

template <typename SparseMatrix>
double sparse_relative_difference(
    const SparseMatrix &left,
    const SparseMatrix &right) {
    if (left.rows() != right.rows() || left.cols() != right.cols())
        return std::numeric_limits<double>::infinity();
    return (left - right).norm()
        / std::max({1.0, left.norm(), right.norm()});
}

void validate_operator_contract(
    const TriMesh &mesh,
    const HelmholtzOperators &operators,
    const char *space_name) {
    const HelmholtzOperators rebuilt = assemble_helmholtz_operators(
        mesh,
        operators.wavenumber,
        operators.diffusion,
        operators.refractive_index,
        operators.boundary_beta);
    const double tolerance = 2e-13;
    if (rebuilt.dirichlet_nodes != operators.dirichlet_nodes
        || sparse_relative_difference(
               rebuilt.stiffness, operators.stiffness) > tolerance
        || sparse_relative_difference(rebuilt.mass, operators.mass) > tolerance
        || sparse_relative_difference(
               rebuilt.boundary_mass, operators.boundary_mass) > tolerance
        || sparse_relative_difference(rebuilt.system, operators.system)
            > tolerance) {
        throw std::invalid_argument(
            std::string(space_name)
            + " Helmholtz operators do not match the current mesh/PDE");
    }
}

Eigen::SparseMatrix<double> select_columns(
    const Eigen::SparseMatrix<double> &matrix,
    const std::vector<int> &columns) {
    std::vector<int> output_column(matrix.cols(), -1);
    for (int column = 0; column < static_cast<int>(columns.size()); ++column) {
        if (columns[column] < 0 || columns[column] >= matrix.cols())
            throw std::out_of_range("selected coarse basis node is out of range");
        if (output_column[columns[column]] >= 0)
            throw std::invalid_argument("selected coarse basis nodes are not unique");
        output_column[columns[column]] = column;
    }
    std::vector<Eigen::Triplet<double>> triplets;
    for (int input_column : columns) {
        const int output = output_column[input_column];
        for (Eigen::SparseMatrix<double>::InnerIterator it(
                 matrix, input_column); it; ++it) {
            triplets.emplace_back(it.row(), output, it.value());
        }
    }
    Eigen::SparseMatrix<double> result(matrix.rows(), columns.size());
    result.setFromTriplets(triplets.begin(), triplets.end());
    result.makeCompressed();
    return result;
}

Eigen::SparseMatrix<double> restrict_square(
    const Eigen::SparseMatrix<double> &matrix,
    const std::vector<int> &dofs) {
    std::vector<int> local_index(matrix.rows(), -1);
    for (int local = 0; local < static_cast<int>(dofs.size()); ++local) {
        if (dofs[local] < 0 || dofs[local] >= matrix.rows())
            throw std::out_of_range("restricted degree of freedom is out of range");
        local_index[dofs[local]] = local;
    }
    std::vector<Eigen::Triplet<double>> triplets;
    for (int local_column = 0;
         local_column < static_cast<int>(dofs.size()); ++local_column) {
        const int global_column = dofs[local_column];
        for (Eigen::SparseMatrix<double>::InnerIterator it(
                 matrix, global_column); it; ++it) {
            const int local_row = local_index[it.row()];
            if (local_row >= 0)
                triplets.emplace_back(local_row, local_column, it.value());
        }
    }
    Eigen::SparseMatrix<double> result(dofs.size(), dofs.size());
    result.setFromTriplets(triplets.begin(), triplets.end());
    result.makeCompressed();
    return result;
}

Eigen::MatrixXd constraint_columns(
    const Eigen::SparseMatrix<double> &interpolation,
    const std::vector<int> &dofs) {
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(
        interpolation.rows(), dofs.size());
    for (int local = 0; local < static_cast<int>(dofs.size()); ++local) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(
                 interpolation, dofs[local]); it; ++it) {
            result(it.row(), local) = it.value();
        }
    }
    return result;
}

Eigen::MatrixXd kernel_basis(const Eigen::MatrixXd &constraints) {
    const int dimension = constraints.cols();
    if (dimension == 0) return Eigen::MatrixXd(0, 0);
    if (constraints.rows() == 0)
        return Eigen::MatrixXd::Identity(dimension, dimension);
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        constraints, Eigen::ComputeFullV);
    const Eigen::VectorXd singular_values = svd.singularValues();
    const double scale = singular_values.size() == 0
        ? 0.0 : singular_values(0);
    const double threshold = kRankTolerance
        * std::max(constraints.rows(), constraints.cols()) * scale;
    int rank = 0;
    for (double singular_value : singular_values) {
        if (singular_value > threshold) ++rank;
    }
    return svd.matrixV().rightCols(dimension - rank);
}

Eigen::MatrixXd energy_orthonormal_basis(
    const Eigen::MatrixXd &basis,
    const Eigen::MatrixXd &energy) {
    if (basis.cols() == 0) return basis;
    const Eigen::MatrixXd gram = basis.transpose() * energy * basis;
    Eigen::LLT<Eigen::MatrixXd> factorization(gram);
    if (factorization.info() != Eigen::Success)
        throw std::runtime_error("kernel energy orthonormalization failed");
    const Eigen::MatrixXd inverse_upper = factorization.matrixU().solve(
        Eigen::MatrixXd::Identity(gram.rows(), gram.cols()));
    return basis * inverse_upper;
}

double spectral_norm(const Eigen::MatrixXd &matrix) {
    if (matrix.rows() == 0 || matrix.cols() == 0) return 0.0;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(matrix);
    return svd.singularValues()(0);
}

double triangle_area(const TriMesh &mesh, const Triangle &triangle) {
    const Point2 first = mesh.nodes[triangle[1]] - mesh.nodes[triangle[0]];
    const Point2 second = mesh.nodes[triangle[2]] - mesh.nodes[triangle[0]];
    return 0.5 * std::abs(
        first.x() * second.y() - first.y() * second.x());
}

ComplexMatrix whitened_matrix(
    const ComplexMatrix &numerator,
    const ComplexMatrix &denominator,
    ComplexMatrix &inverse_lower) {
    if (numerator.rows() != numerator.cols()
        || denominator.rows() != denominator.cols()
        || numerator.rows() != denominator.rows()) {
        throw std::invalid_argument("generalized spectrum dimensions disagree");
    }
    const ComplexMatrix hermitian_denominator =
        0.5 * (denominator + denominator.adjoint());
    Eigen::LLT<ComplexMatrix> factorization(hermitian_denominator);
    if (factorization.info() != Eigen::Success)
        throw std::runtime_error("coarse energy matrix is not positive definite");
    const ComplexMatrix lower = factorization.matrixL();
    inverse_lower = lower.triangularView<Eigen::Lower>().solve(
        ComplexMatrix::Identity(lower.rows(), lower.cols()));
    const ComplexMatrix hermitian_numerator =
        0.5 * (numerator + numerator.adjoint());
    ComplexMatrix result = inverse_lower * hermitian_numerator
        * inverse_lower.adjoint();
    return 0.5 * (result + result.adjoint());
}

LocalizationSpectrum largest_generalized_eigenvalue(
    const ComplexMatrix &numerator,
    const ComplexMatrix &denominator,
    const LocalizationEigenConfig &config) {
    if (config.maximum_iterations <= 0
        || !(config.relative_tolerance > 0.0)
        || config.dense_cross_check_max_dimension < 0) {
        throw std::invalid_argument("localization eigen configuration is invalid");
    }
    ComplexMatrix inverse_lower;
    const ComplexMatrix whitened = whitened_matrix(
        numerator, denominator, inverse_lower);
    const int dimension = whitened.rows();
    LocalizationSpectrum result;
    if (dimension == 0) {
        result.converged = true;
        return result;
    }

    ComplexVector iterate;
    if (config.warm_start.size() == dimension
        && config.warm_start.allFinite()
        && config.warm_start.norm() > 0.0) {
        const ComplexMatrix lower_inverse_adjoint =
            inverse_lower.adjoint();
        // y=L^*x and x=L^{-H}y.
        iterate = lower_inverse_adjoint.fullPivLu().solve(
            config.warm_start);
    } else {
        iterate = ComplexVector::Ones(dimension);
        for (int index = 0; index < dimension; ++index)
            iterate(index) += Complex(0.0, 0.125 * (index + 1));
    }
    iterate.normalize();
    const double matrix_scale = std::max(1.0, whitened.norm());
    for (int iteration = 1; iteration <= config.maximum_iterations; ++iteration) {
        const ComplexVector applied = whitened * iterate;
        const double norm = applied.norm();
        if (norm <= std::numeric_limits<double>::epsilon() * matrix_scale) {
            result.lambda_max = 0.0;
            result.relative_residual = 0.0;
            result.iterations = iteration;
            result.converged = true;
            break;
        }
        iterate = applied / norm;
        const ComplexVector next = whitened * iterate;
        result.lambda_max = std::max(
            0.0, std::real(iterate.dot(next)));
        result.relative_residual =
            (next - result.lambda_max * iterate).norm() / matrix_scale;
        result.iterations = iteration;
        if (result.relative_residual <= config.relative_tolerance) {
            result.converged = true;
            break;
        }
    }
    if (!result.converged)
        throw std::runtime_error("localization largest-eigenvalue iteration did not converge");

    result.dominant_vector = inverse_lower.adjoint() * iterate;
    const double energy_norm = std::sqrt(std::max(
        0.0, std::real(result.dominant_vector.dot(
            denominator * result.dominant_vector))));
    if (energy_norm > 0.0) result.dominant_vector /= energy_norm;

    if (dimension <= config.dense_cross_check_max_dimension) {
        Eigen::SelfAdjointEigenSolver<ComplexMatrix> dense_solver(whitened);
        if (dense_solver.info() != Eigen::Success)
            throw std::runtime_error("dense localization eigen cross-check failed");
        result.dense_cross_checked = true;
        result.dense_lambda_max = std::max(
            0.0, dense_solver.eigenvalues().maxCoeff());
        result.dense_relative_difference = std::abs(
            result.lambda_max - result.dense_lambda_max)
            / std::max(1.0, result.dense_lambda_max);
    }
    return result;
}

double dense_largest_generalized_eigenvalue(
    const ComplexMatrix &numerator,
    const ComplexMatrix &denominator) {
    ComplexMatrix inverse_lower;
    const ComplexMatrix whitened = whitened_matrix(
        numerator, denominator, inverse_lower);
    if (whitened.rows() == 0) return 0.0;
    Eigen::SelfAdjointEigenSolver<ComplexMatrix> solver(whitened);
    if (solver.info() != Eigen::Success)
        throw std::runtime_error("dense generalized eigenproblem failed");
    return std::max(0.0, solver.eigenvalues().maxCoeff());
}

} // namespace

ComplexVector ReferenceRetraction::apply(
    const ComplexVector &ambient_values) const {
    if (ambient_values.size() != retraction_.cols()
        || !ambient_values.allFinite()) {
        throw std::invalid_argument("reference retraction vector has the wrong size");
    }
    return retraction_.cast<Complex>() * ambient_values;
}

ComplexMatrix ReferenceRetraction::apply(
    const ComplexMatrix &ambient_values) const {
    if (ambient_values.rows() != retraction_.cols()
        || !ambient_values.allFinite()) {
        throw std::invalid_argument("reference retraction matrix has the wrong row count");
    }
    return retraction_.cast<Complex>() * ambient_values;
}

ReferenceRetraction build_reference_retraction(
    const ReferenceEpochHierarchy &hierarchy) {
    if (!hierarchy.reference_embedding_holds())
        throw std::invalid_argument("reference epoch embedding is invalid");
    const TriMesh &reference = hierarchy.reference_mesh();
    const TriMesh &ambient = hierarchy.ambient_mesh();
    const int reference_nodes = static_cast<int>(reference.nodes.size());
    const int ambient_nodes = static_cast<int>(ambient.nodes.size());

    std::unordered_map<PointKey, int, PointKeyHash> ambient_node;
    ambient_node.reserve(2 * ambient.nodes.size());
    for (int node = 0; node < ambient_nodes; ++node) {
        if (!ambient_node.emplace(point_key(ambient.nodes[node]), node).second)
            throw std::runtime_error("ambient mesh contains duplicate coordinates");
    }
    // Verify geometric nesting independently of the prolongation before
    // assembling the local projector.
    for (int node = 0; node < reference_nodes; ++node) {
        if (!ambient_node.contains(point_key(reference.nodes[node]))) {
            throw std::runtime_error(
                "reference node is missing from the nested ambient mesh");
        }
    }

    std::vector<char> reference_is_dirichlet(reference_nodes, false);
    for (int node : dirichlet_nodes(reference))
        reference_is_dirichlet[node] = true;
    std::vector<std::vector<int>> incident(reference_nodes);
    for (int element = 0;
         element < static_cast<int>(reference.elems.size()); ++element) {
        for (int node : reference.elems[element])
            incident[node].push_back(element);
    }
    std::vector<std::vector<int>> ambient_children(reference.elems.size());
    const std::vector<int> &ambient_parent =
        hierarchy.ambient_parent_reference_elements();
    for (int element = 0;
         element < static_cast<int>(ambient.elems.size()); ++element) {
        const int parent = ambient_parent[element];
        if (parent < 0 || parent >= static_cast<int>(reference.elems.size()))
            throw std::runtime_error("ambient/reference element parent is invalid");
        ambient_children[parent].push_back(element);
    }

    std::vector<Eigen::Triplet<double>> projector_triplets;
    for (int node = 0; node < reference_nodes; ++node) {
        if (reference_is_dirichlet[node]) continue;
        if (incident[node].empty())
            throw std::runtime_error("reference node has no incident element");
        const int reference_element = incident[node].front();
        const Triangle &reference_triangle =
            reference.elems[reference_element];
        int local_node = -1;
        for (int local = 0; local < 3; ++local) {
            if (reference_triangle[local] == node) local_node = local;
        }
        if (local_node < 0)
            throw std::runtime_error("reference node incidence is inconsistent");
        const double area = triangle_area(reference, reference_triangle);
        if (!(area > 0.0))
            throw std::runtime_error("reference projector element is degenerate");
        Eigen::Matrix3d local_mass;
        local_mass << 2.0, 1.0, 1.0,
                      1.0, 2.0, 1.0,
                      1.0, 1.0, 2.0;
        local_mass *= area / 12.0;
        const Eigen::Vector3d dual = local_mass.ldlt().solve(
            Eigen::Vector3d::Unit(local_node));

        std::unordered_map<int, double> row_weights;
        for (int ambient_element : ambient_children[reference_element]) {
            const Triangle &child = ambient.elems[ambient_element];
            const double child_area = triangle_area(ambient, child);
            Eigen::Vector3d dual_values;
            for (int local = 0; local < 3; ++local) {
                dual_values(local) = 0.0;
                for (int reference_local = 0; reference_local < 3;
                     ++reference_local) {
                    dual_values(local) += dual(reference_local)
                        * hierarchy.reference_to_ambient().coeff(
                            child[local],
                            reference_triangle[reference_local]);
                }
            }
            for (int local = 0; local < 3; ++local) {
                const double weight = child_area / 12.0
                    * (2.0 * dual_values(local)
                       + dual_values((local + 1) % 3)
                       + dual_values((local + 2) % 3));
                row_weights[child[local]] += weight;
            }
        }
        for (const auto &[ambient_node_index, weight] : row_weights) {
            if (std::abs(weight) > 1e-15)
                projector_triplets.emplace_back(node, ambient_node_index, weight);
        }
    }

    ReferenceRetraction result;
    result.projector_.resize(reference_nodes, ambient_nodes);
    result.projector_.setFromTriplets(
        projector_triplets.begin(), projector_triplets.end());
    result.projector_.makeCompressed();

    Eigen::SparseMatrix<double> space_identity(reference_nodes, reference_nodes);
    std::vector<Eigen::Triplet<double>> identity_triplets;
    for (int node = 0; node < reference_nodes; ++node) {
        if (!reference_is_dirichlet[node])
            identity_triplets.emplace_back(node, node, 1.0);
    }
    space_identity.setFromTriplets(
        identity_triplets.begin(), identity_triplets.end());
    Eigen::SparseMatrix<double> reference_kernel_projection = space_identity
        - hierarchy.coarse_to_reference()
            * hierarchy.reference_quasi_interpolation();
    reference_kernel_projection.makeCompressed();
    result.retraction_ = reference_kernel_projection * result.projector_;
    result.retraction_.prune(0.0, 1e-14);
    result.retraction_.makeCompressed();

    const Eigen::SparseMatrix<double> projector_identity =
        result.projector_ * hierarchy.reference_to_ambient() - space_identity;
    result.diagnostics_.projector_identity_relative_error =
        projector_identity.norm() / std::max(1.0, space_identity.norm());
    const Eigen::SparseMatrix<double> kernel_identity =
        result.retraction_ * hierarchy.reference_to_ambient()
        - reference_kernel_projection;
    result.diagnostics_.reference_kernel_identity_relative_error =
        kernel_identity.norm()
        / std::max(1.0, reference_kernel_projection.norm());
    const Eigen::SparseMatrix<double> kernel_constraint =
        hierarchy.reference_quasi_interpolation() * result.retraction_;
    result.diagnostics_.kernel_constraint_relative_error =
        kernel_constraint.norm() / std::max(1.0, result.retraction_.norm());
    for (int column = 0; column < result.retraction_.outerSize(); ++column) {
        int nonzeros = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(
                 result.retraction_, column); it; ++it) {
            ++nonzeros;
        }
        result.diagnostics_.maximum_column_nonzeros = std::max(
            result.diagnostics_.maximum_column_nonzeros, nonzeros);
    }
    return result;
}

ReferenceLocalizationCertificate compute_reference_localization_certificate(
    const ReferenceEpochHierarchy &hierarchy,
    const HelmholtzOperators &reference_operators,
    const HelmholtzOperators &ambient_operators,
    const ComplexSparseMatrix &localized_adjoint_basis,
    const std::vector<int> &coarse_basis_nodes,
    KernelRieszSolver riesz_solver,
    const LocalizationEigenConfig &eigen_config) {
    const int reference_nodes = static_cast<int>(
        hierarchy.reference_mesh().nodes.size());
    const int ambient_nodes = static_cast<int>(
        hierarchy.ambient_mesh().nodes.size());
    if (coarse_basis_nodes.empty()
        || localized_adjoint_basis.rows() != reference_nodes
        || localized_adjoint_basis.cols()
            != static_cast<int>(coarse_basis_nodes.size())
        || reference_operators.system.rows() != reference_nodes
        || reference_operators.system.cols() != reference_nodes
        || reference_operators.stiffness.rows() != reference_nodes
        || reference_operators.mass.rows() != reference_nodes
        || ambient_operators.system.rows() != ambient_nodes
        || ambient_operators.system.cols() != ambient_nodes
        || ambient_operators.stiffness.rows() != ambient_nodes
        || ambient_operators.mass.rows() != ambient_nodes
        || !localized_adjoint_basis.coeffs().allFinite()
        || reference_operators.wavenumber != ambient_operators.wavenumber) {
        throw std::invalid_argument(
            "reference localization certificate inputs are inconsistent");
    }
    std::unordered_set<int> dirichlet;
    for (int node : dirichlet_nodes(hierarchy.coarse_mesh()))
        dirichlet.insert(node);
    for (int node : coarse_basis_nodes) {
        if (dirichlet.contains(node)) {
            throw std::invalid_argument(
                "localization certificate basis includes a Dirichlet coarse node");
        }
    }
    validate_operator_contract(
        hierarchy.reference_mesh(), reference_operators, "reference");
    validate_operator_contract(
        hierarchy.ambient_mesh(), ambient_operators, "ambient");
    if (reference_operators.boundary_beta != ambient_operators.boundary_beta
        || reference_operators.diffusion.size()
            != hierarchy.reference_mesh().elems.size()
        || ambient_operators.diffusion.size()
            != hierarchy.ambient_mesh().elems.size()
        || reference_operators.refractive_index.size()
            != hierarchy.reference_mesh().elems.size()
        || ambient_operators.refractive_index.size()
            != hierarchy.ambient_mesh().elems.size()) {
        throw std::invalid_argument(
            "reference/ambient PDE coefficient contracts disagree");
    }
    const std::vector<int> &ambient_parent =
        hierarchy.ambient_parent_reference_elements();
    for (int element = 0;
         element < static_cast<int>(ambient_parent.size()); ++element) {
        const int parent = ambient_parent[element];
        if (ambient_operators.diffusion[element]
                != reference_operators.diffusion[parent]
            || ambient_operators.refractive_index[element]
                != reference_operators.refractive_index[parent]) {
            throw std::invalid_argument(
                "ambient PDE coefficients do not inherit the reference parent");
        }
    }

    ReferenceLocalizationCertificate result;
    result.retraction = build_reference_retraction(hierarchy);
    const ComplexMatrix localized(localized_adjoint_basis);
    const ComplexMatrix reference_action =
        reference_operators.system.adjoint() * localized;
    result.defect_rhs = result.retraction.matrix().transpose().cast<Complex>()
        * reference_action;
    result.ambient_riesz = compute_ambient_defect_riesz(
        hierarchy, ambient_operators, result.defect_rhs, riesz_solver);

    const Eigen::SparseMatrix<double> coarse_basis = select_columns(
        hierarchy.coarse_to_reference(), coarse_basis_nodes);
    const Eigen::SparseMatrix<double> reference_energy =
        energy_matrix(reference_operators);
    result.coarse_energy = ComplexMatrix(
        coarse_basis.transpose() * reference_energy * coarse_basis).cast<Complex>();
    result.spectrum = largest_generalized_eigenvalue(
        result.ambient_riesz.gram, result.coarse_energy, eigen_config);
    result.theta_loc = std::sqrt(std::max(0.0, result.spectrum.lambda_max));
    return result;
}

SmallMatrixLocalizationValidation
validate_reference_localization_certificate_small_matrix(
    const ReferenceEpochHierarchy &hierarchy,
    const HelmholtzOperators &reference_operators,
    const HelmholtzOperators &ambient_operators,
    const ComplexSparseMatrix &localized_adjoint_basis,
    const ComplexSparseMatrix &ideal_adjoint_basis,
    const std::vector<int> &coarse_basis_nodes,
    const ReferenceLocalizationCertificate &certificate,
    int maximum_free_dofs) {
    if (maximum_free_dofs <= 0)
        throw std::invalid_argument("small-matrix DOF limit must be positive");
    if (localized_adjoint_basis.rows() != ideal_adjoint_basis.rows()
        || localized_adjoint_basis.cols() != ideal_adjoint_basis.cols()
        || localized_adjoint_basis.cols()
            != static_cast<int>(coarse_basis_nodes.size())) {
        throw std::invalid_argument(
            "ideal/localized reference basis dimensions disagree");
    }

    SmallMatrixLocalizationValidation result;
    const Eigen::SparseMatrix<double> reference_energy =
        energy_matrix(reference_operators);
    const ComplexMatrix localized_dense(localized_adjoint_basis);
    const ComplexMatrix expected_defect =
        certificate.retraction.matrix().transpose().cast<Complex>()
        * (reference_operators.system.adjoint() * localized_dense);
    if (expected_defect.rows() != certificate.defect_rhs.rows()
        || expected_defect.cols() != certificate.defect_rhs.cols()
        || !expected_defect.isApprox(certificate.defect_rhs, 2e-12)) {
        throw std::invalid_argument(
            "localization certificate is stale for the localized basis");
    }
    const Eigen::SparseMatrix<double> expected_coarse_basis = select_columns(
        hierarchy.coarse_to_reference(), coarse_basis_nodes);
    const ComplexMatrix expected_coarse_energy = ComplexMatrix(
        expected_coarse_basis.transpose() * reference_energy
        * expected_coarse_basis).cast<Complex>();
    if (!expected_coarse_energy.isApprox(certificate.coarse_energy, 2e-12)) {
        throw std::invalid_argument(
            "localization certificate coarse energy is stale");
    }
    const ComplexMatrix difference = ComplexMatrix(ideal_adjoint_basis)
        - localized_dense;
    const ComplexMatrix direct_gram = difference.adjoint()
        * reference_energy.cast<Complex>() * difference;
    result.direct_delta = std::sqrt(dense_largest_generalized_eigenvalue(
        direct_gram, certificate.coarse_energy));

    const int ambient_nodes = static_cast<int>(
        hierarchy.ambient_mesh().nodes.size());
    std::vector<char> is_dirichlet(ambient_nodes, false);
    for (int node : ambient_operators.dirichlet_nodes)
        is_dirichlet[node] = true;
    std::vector<int> free_dofs;
    std::vector<int> free_index(ambient_nodes, -1);
    for (int node = 0; node < ambient_nodes; ++node) {
        if (is_dirichlet[node]) continue;
        free_index[node] = static_cast<int>(free_dofs.size());
        free_dofs.push_back(node);
    }
    if (static_cast<int>(free_dofs.size()) > maximum_free_dofs) {
        throw std::invalid_argument(
            "small-matrix validation exceeds its ambient DOF limit");
    }

    const Eigen::SparseMatrix<double> ambient_energy_sparse =
        energy_matrix(ambient_operators);
    const Eigen::MatrixXd ambient_energy(
        restrict_square(ambient_energy_sparse, free_dofs));
    std::vector<int> reference_free_dofs;
    std::vector<char> reference_is_dirichlet(
        hierarchy.reference_mesh().nodes.size(), false);
    for (int node : reference_operators.dirichlet_nodes)
        reference_is_dirichlet[node] = true;
    for (int node = 0;
         node < static_cast<int>(reference_is_dirichlet.size()); ++node) {
        if (!reference_is_dirichlet[node])
            reference_free_dofs.push_back(node);
    }
    const Eigen::MatrixXd reference_free_energy(
        restrict_square(energy_matrix(reference_operators),
                        reference_free_dofs));
    const Eigen::MatrixXd projector_dense(
        certificate.retraction.projector());
    const Eigen::MatrixXd retraction_dense(
        certificate.retraction.matrix());
    Eigen::MatrixXd projector_free(
        reference_free_dofs.size(), free_dofs.size());
    Eigen::MatrixXd retraction_free(
        reference_free_dofs.size(), free_dofs.size());
    for (int row = 0; row < static_cast<int>(reference_free_dofs.size()); ++row) {
        for (int column = 0; column < static_cast<int>(free_dofs.size()); ++column) {
            projector_free(row, column) = projector_dense(
                reference_free_dofs[row], free_dofs[column]);
            retraction_free(row, column) = retraction_dense(
                reference_free_dofs[row], free_dofs[column]);
        }
    }
    const ComplexMatrix ambient_energy_complex =
        ambient_energy.cast<Complex>();
    result.projector_stability_constant = std::sqrt(
        dense_largest_generalized_eigenvalue(
            (projector_free.transpose() * reference_free_energy
             * projector_free).cast<Complex>(),
            ambient_energy_complex));
    result.retraction_stability_constant = std::sqrt(
        dense_largest_generalized_eigenvalue(
            (retraction_free.transpose() * reference_free_energy
             * retraction_free).cast<Complex>(),
            ambient_energy_complex));
    Eigen::SparseMatrix<double> ambient_real_sparse =
        ambient_operators.stiffness
        - ambient_operators.wavenumber * ambient_operators.wavenumber
            * ambient_operators.mass;
    ambient_real_sparse.makeCompressed();
    const Eigen::MatrixXd ambient_real(
        restrict_square(ambient_real_sparse, free_dofs));
    const Eigen::MatrixXd global_constraints = constraint_columns(
        hierarchy.ambient_quasi_interpolation(), free_dofs);
    const Eigen::MatrixXd global_kernel = kernel_basis(global_constraints);
    result.ambient_kernel_dimension = global_kernel.cols();
    if (global_kernel.cols() == 0) {
        result.ambient_kernel_coercivity = 1.0;
        result.stable_decomposition_constant = 1.0;
        result.upper_certificate = certificate.theta_loc;
        result.one_sided_control_holds =
            result.direct_delta <= result.upper_certificate + 2e-10;
        return result;
    }

    const Eigen::MatrixXd kernel_energy = global_kernel.transpose()
        * ambient_energy * global_kernel;
    const Eigen::MatrixXd kernel_real = global_kernel.transpose()
        * ambient_real * global_kernel;
    Eigen::LLT<Eigen::MatrixXd> kernel_energy_factor(kernel_energy);
    if (kernel_energy_factor.info() != Eigen::Success)
        throw std::runtime_error("ambient kernel energy factorization failed");
    const Eigen::MatrixXd inverse_lower =
        kernel_energy_factor.matrixL().solve(
            Eigen::MatrixXd::Identity(
                kernel_energy.rows(), kernel_energy.cols()));
    Eigen::MatrixXd coercivity_matrix = inverse_lower * kernel_real
        * inverse_lower.transpose();
    coercivity_matrix = 0.5
        * (coercivity_matrix + coercivity_matrix.transpose());
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> coercivity_solver(
        coercivity_matrix);
    if (coercivity_solver.info() != Eigen::Success)
        throw std::runtime_error("ambient kernel coercivity eigenproblem failed");
    result.ambient_kernel_coercivity =
        coercivity_solver.eigenvalues().minCoeff();
    if (!(result.ambient_kernel_coercivity > 0.0))
        throw std::runtime_error("ambient kernel is not Helmholtz coercive");

    std::vector<Eigen::MatrixXd> local_bases;
    int total_columns = 0;
    for (const KernelRieszPatch &patch : certificate.ambient_riesz.patches) {
        Eigen::MatrixXd basis = kernel_basis(patch.constraints);
        if (basis.cols() > 0) {
            const Eigen::MatrixXd local_energy(
                restrict_square(ambient_energy_sparse, patch.discrete_dofs));
            basis = energy_orthonormal_basis(basis, local_energy);
            total_columns += basis.cols();
        }
        local_bases.push_back(std::move(basis));
    }
    result.local_kernel_columns = total_columns;
    if (total_columns == 0)
        throw std::runtime_error("ambient local kernels are all trivial");

    Eigen::MatrixXd synthesis = Eigen::MatrixXd::Zero(
        free_dofs.size(), total_columns);
    int first_column = 0;
    for (int patch_index = 0;
         patch_index < static_cast<int>(local_bases.size()); ++patch_index) {
        const Eigen::MatrixXd &basis = local_bases[patch_index];
        const KernelRieszPatch &patch =
            certificate.ambient_riesz.patches[patch_index];
        for (int local = 0;
             local < static_cast<int>(patch.discrete_dofs.size()); ++local) {
            const int row = free_index[patch.discrete_dofs[local]];
            if (row < 0)
                throw std::runtime_error(
                    "ambient local kernel contains a Dirichlet node");
            synthesis.block(row, first_column, 1, basis.cols()) =
                basis.row(local);
        }
        first_column += basis.cols();
    }
    const Eigen::MatrixXd normalized_global = energy_orthonormal_basis(
        global_kernel, ambient_energy);
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> decomposition(
        synthesis);
    const Eigen::MatrixXd coefficients = decomposition.solve(normalized_global);
    result.decomposition_relative_residual =
        (synthesis * coefficients - normalized_global).norm()
        / std::max(1.0, normalized_global.norm());
    result.stable_decomposition_constant = spectral_norm(coefficients);
    result.upper_certificate = result.stable_decomposition_constant
        / result.ambient_kernel_coercivity * certificate.theta_loc;
    const double tolerance = 2e-9 * std::max({
        1.0, result.direct_delta, result.upper_certificate});
    result.one_sided_control_holds =
        result.decomposition_relative_residual <= 2e-10
        && result.direct_delta <= result.upper_certificate + tolerance;
    return result;
}

} // namespace lod2d::helmholtz::adaptive

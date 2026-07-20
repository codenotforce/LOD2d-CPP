#include "helmholtz/corrector.h"

#include <stdexcept>
#include <vector>

namespace lod2d::helmholtz {

ComplexSparseMatrix build_helmholtz_corrector_matrix(
    const TriMesh &coarse,
    int fine_node_count,
    const std::vector<HelmholtzElementCorrector> &correctors) {
    if (correctors.size() != coarse.elems.size())
        throw std::invalid_argument("Helmholtz corrector count must match coarse elements");
    std::vector<ComplexTriplet> triplets;
    std::size_t nonzeros = 0;
    for (const auto &element : correctors) nonzeros += element.size();
    triplets.reserve(nonzeros);
    for (int element = 0; element < static_cast<int>(coarse.elems.size()); ++element) {
        for (const auto &entry : correctors[element]) {
            triplets.emplace_back(
                entry.row,
                coarse.elems[element][entry.local_coarse_vertex],
                entry.value);
        }
    }
    ComplexSparseMatrix matrix(fine_node_count, static_cast<int>(coarse.nodes.size()));
    matrix.setFromTriplets(triplets.begin(), triplets.end());
    return matrix;
}

ComplexSparseMatrix build_helmholtz_corrected_basis(
    const Eigen::SparseMatrix<double> &coarse_to_fine,
    const TriMesh &coarse,
    int fine_node_count,
    const std::vector<HelmholtzElementCorrector> &correctors) {
    ComplexSparseMatrix basis = coarse_to_fine.cast<Complex>();
    basis -= build_helmholtz_corrector_matrix(
        coarse, fine_node_count, correctors);
    basis.prune(Complex(0.0, 0.0), 1e-14);
    basis.makeCompressed();
    return basis;
}

} // namespace lod2d::helmholtz

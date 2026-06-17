#include "fem/assemble_dg.h"

namespace lod2d {

Eigen::SparseMatrix<double> assemble_dg(const TriMesh &mesh,
                                         const std::vector<double> &coeff) {
    // Placeholder — will be implemented in Phase C
    int ndof = 3 * static_cast<int>(mesh.elems.size());
    Eigen::SparseMatrix<double> S(ndof, ndof);
    return S;
}

} // namespace lod2d

#include "fem/assemble_dg.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace lod2d {

ElementStiffnessBlocks assemble_element_stiffness(const TriMesh &mesh,
                                                   const std::vector<double> &coeff) {
    const int nt = static_cast<int>(mesh.elems.size());
    if (coeff.size() != mesh.elems.size())
        throw std::invalid_argument("assemble_element_stiffness coefficient count must match element count");

    ElementStiffnessBlocks blocks(nt);
    for (int t = 0; t < nt; ++t) {
        int a = mesh.elems[t][0];
        int b = mesh.elems[t][1];
        int c = mesh.elems[t][2];

        const auto &pa = mesh.nodes[a];
        const auto &pb = mesh.nodes[b];
        const auto &pc = mesh.nodes[c];

        double ve1_x = pc.x() - pb.x(),  ve1_y = pc.y() - pb.y();
        double ve2_x = pa.x() - pc.x(),  ve2_y = pa.y() - pc.y();
        double ve3_x = pb.x() - pa.x(),  ve3_y = pb.y() - pa.y();

        double area2 = std::abs(ve3_x * ve2_y - ve3_y * ve2_x);
        double area  = 0.5 * area2;
        double inv_a2 = (area2 > 0) ? 1.0 / area2 : 0.0;

        double G[3][2] = {
            {-ve1_y * inv_a2,  ve1_x * inv_a2},
            {-ve2_y * inv_a2,  ve2_x * inv_a2},
            {-ve3_y * inv_a2,  ve3_x * inv_a2}
        };

        double cA = coeff[t] * area;
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                blocks[t](j, k) = cA * (G[j][0] * G[k][0] + G[j][1] * G[k][1]);
            }
        }
    }
    return blocks;
}

Eigen::SparseMatrix<double> assemble_dg_from_element_stiffness(
    const ElementStiffnessBlocks &element_stiffness) {
    const int nt = static_cast<int>(element_stiffness.size());
    const int ndof = 3 * nt;

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(9 * nt);

    for (int t = 0; t < nt; ++t) {
        int dg[3] = {3*t, 3*t + 1, 3*t + 2};
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                const double val = element_stiffness[t](j, k);
                if (val != 0.0)
                    triplets.emplace_back(dg[j], dg[k], val);
            }
        }
    }

    Eigen::SparseMatrix<double> S(ndof, ndof);
    S.setFromTriplets(triplets.begin(), triplets.end());
    return S;
}

Eigen::SparseMatrix<double> assemble_dg(const TriMesh &mesh,
                                         const std::vector<double> &coeff) {
    return assemble_dg_from_element_stiffness(assemble_element_stiffness(mesh, coeff));
}

} // namespace lod2d

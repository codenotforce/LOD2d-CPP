#include "fem/assemble_dg.h"
#include <Eigen/Sparse>
#include <vector>
#include <cmath>

namespace lod2d {

Eigen::SparseMatrix<double> assemble_dg(const TriMesh &mesh,
                                         const std::vector<double> &coeff) {
    int nt = static_cast<int>(mesh.elems.size());
    int ndof = 3 * nt;

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(9 * nt);

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
        int dg[3] = {3*t, 3*t + 1, 3*t + 2};

        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                double val = cA * (G[j][0] * G[k][0] + G[j][1] * G[k][1]);
                if (val != 0.0)
                    triplets.emplace_back(dg[j], dg[k], val);
            }
        }
    }

    Eigen::SparseMatrix<double> S(ndof, ndof);
    S.setFromTriplets(triplets.begin(), triplets.end());
    return S;
}

} // namespace lod2d

#include "lod/quasi_interp.h"
#include "mesh/refine.h"
#include <cmath>
#include <iostream>
#include <string>

using namespace lod2d;

namespace {

TriMesh unit_square() {
    TriMesh mesh;
    mesh.nodes = {{0,0}, {1,0}, {1,1}, {0,1}};
    mesh.elems = {{0,1,3}, {1,2,3}};
    mesh.dirichlet = {0,1,2,3};
    return mesh;
}

bool check(bool condition, const std::string &name, int &failed) {
    if (condition) {
        std::cout << "  [PASS] " << name << "\n";
        return true;
    }
    std::cout << "  [FAIL] " << name << "\n";
    ++failed;
    return false;
}

} // namespace

int main() {
    std::cout << "=== Quasi-Interpolation Tests ===\n";
    int failed = 0;

    TriMesh T0 = unit_square();
    auto coarse_out = refine_mesh(T0, 3);
    auto fine_out = refine_mesh(coarse_out.mesh, 2);
    const auto &coarse = coarse_out.mesh;
    const auto &fine = fine_out.mesh;

    const int Nhdg = static_cast<int>(fine.elems.size()) * 3;
    std::vector<Eigen::Triplet<double>> cg_t;
    cg_t.reserve(Nhdg);
    for (int e = 0; e < static_cast<int>(fine.elems.size()); ++e) {
        for (int i = 0; i < 3; ++i) cg_t.emplace_back(3 * e + i, fine.elems[e][i], 1.0);
    }
    Eigen::SparseMatrix<double> cg2dgh(Nhdg, static_cast<int>(fine.nodes.size()));
    cg2dgh.setFromTriplets(cg_t.begin(), cg_t.end());

    auto IH = build_quasi_interp(coarse, fine, fine_out.P_dg, cg2dgh,
                                  static_cast<int>(fine.nodes.size()),
                                  static_cast<int>(coarse.nodes.size()));

    check(IH.rows() == static_cast<int>(coarse.nodes.size()), "IH row count", failed);
    check(IH.cols() == static_cast<int>(fine.nodes.size()), "IH column count", failed);
    check(IH.nonZeros() > 0, "IH nonzero", failed);

    bool finite = true;
    for (int col = 0; col < IH.outerSize(); ++col) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(IH, col); it; ++it) {
            if (!std::isfinite(it.value())) finite = false;
        }
    }
    check(finite, "IH finite values", failed);

    if (failed == 0) {
        std::cout << "\nAll quasi-interpolation tests passed!\n";
        return 0;
    }
    std::cout << "\nQuasi-interpolation tests failed: " << failed << "\n";
    return 1;
}

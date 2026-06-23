#include "lod/patches.h"
#include "mesh/refine.h"
#include <cmath>
#include <iostream>
#include <set>
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
    std::cout << "=== Patch Construction Tests ===\n";
    int failed = 0;

    TriMesh T0 = unit_square();
    auto coarse_out = refine_mesh(T0, 4);
    const auto &mesh = coarse_out.mesh;
    const int nt = static_cast<int>(mesh.elems.size());

    for (int ell : {0, 1, 2}) {
        auto patch = build_patches(mesh, ell);
        check(patch.rows() == nt && patch.cols() == nt, "patch size ell=" + std::to_string(ell), failed);
        bool diagonal_ok = true;
        bool finite_ok = true;
        bool monotone_ok = true;
        for (int col = 0; col < patch.outerSize(); ++col) {
            if (patch.coeff(col, col) <= 0.0) diagonal_ok = false;
            for (Eigen::SparseMatrix<double>::InnerIterator it(patch, col); it; ++it) {
                if (!std::isfinite(it.value()) || it.value() <= 0.0) finite_ok = false;
                if (ell == 0 && it.row() != col) monotone_ok = false;
            }
        }
        check(diagonal_ok, "patch contains seed ell=" + std::to_string(ell), failed);
        check(finite_ok, "patch finite levels ell=" + std::to_string(ell), failed);
        check(monotone_ok, "patch ell=0 identity behavior", failed);
    }

    if (failed == 0) {
        std::cout << "\nAll patch tests passed!\n";
        return 0;
    }
    std::cout << "\nPatch tests failed: " << failed << "\n";
    return 1;
}

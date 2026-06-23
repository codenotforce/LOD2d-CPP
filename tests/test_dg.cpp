#include "fem/assemble_dg.h"
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

bool near(double a, double b, double tol = 1e-12) {
    return std::abs(a - b) <= tol;
}

} // namespace

int main() {
    std::cout << "=== DG Assembly Tests ===\n";
    int failed = 0;

    {
        TriMesh simple;
        simple.nodes = {{0,0}, {1,0}, {0,1}};
        simple.elems = {{0,1,2}};
        auto S = assemble_dg(simple, std::vector<double>{1.0});
        double expected[9] = {1.0, -0.5, -0.5, -0.5, 0.5, 0.0, -0.5, 0.0, 0.5};
        bool ok = true;
        for (int j = 0; j < 3; ++j) {
            for (int i = 0; i < 3; ++i) {
                if (!near(S.coeff(i, j), expected[j * 3 + i])) ok = false;
            }
        }
        check(ok, "A=1 right triangle values", failed);
    }

    {
        TriMesh mesh = refine_mesh(unit_square(), 5).mesh;
        std::vector<double> coeff(mesh.elems.size(), 1.0);
        auto S = assemble_dg(mesh, coeff);
        const int ndg = 3 * static_cast<int>(mesh.elems.size());
        check(S.rows() == ndg && S.cols() == ndg, "refined DG matrix size", failed);
        check(S.nonZeros() > 0, "refined DG nonzero", failed);
        bool symmetric = true;
        for (int col = 0; col < S.outerSize(); ++col) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(S, col); it; ++it) {
                if (!near(it.value(), S.coeff(it.col(), it.row()), 1e-10)) symmetric = false;
            }
        }
        check(symmetric, "refined DG symmetry", failed);
    }

    if (failed == 0) {
        std::cout << "\nAll DG tests passed!\n";
        return 0;
    }
    std::cout << "\nDG tests failed: " << failed << "\n";
    return 1;
}

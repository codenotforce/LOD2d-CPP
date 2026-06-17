#include "mesh/refine.h"
#include <iostream>
#include <chrono>
#include <cmath>

using namespace lod2d;

int main() {
    std::cout << "=== LOD2d C++ — Mesh Refinement Test ===\n\n";

    // Unit square: two triangles
    TriMesh mesh;
    mesh.nodes = {{0,0}, {1,0}, {1,1}, {0,1}};
    mesh.elems = {{0,1,3}, {1,2,3}};
    mesh.dirichlet = {0,1,2,3};

    int nref = 6;
    std::cout << "Initial: " << mesh.nodes.size() << " nodes, "
              << mesh.elems.size() << " elements\n";
    std::cout << "Refining " << nref << " levels...\n";

    auto t0 = std::chrono::high_resolution_clock::now();
    RefineOutput out = refine_mesh(mesh, nref);
    auto t1 = std::chrono::high_resolution_clock::now();

    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    const auto &fine = out.mesh;

    std::cout << "Result:  " << fine.nodes.size() << " nodes, "
              << fine.elems.size() << " elements\n";
    std::cout << "Time:    " << ms << " ms\n\n";

    // Verify expected counts
    int expected_elems = 2;
    for (int i = 0; i < nref; ++i) expected_elems *= 4;
    int expected_nodes_approx = expected_elems / 2 + 85;  // rough

    bool elems_ok = (fine.elems.size() == static_cast<size_t>(expected_elems));
    std::cout << "Element count check: " << (elems_ok ? "PASS" : "FAIL")
              << " (expected " << expected_elems << ", got " << fine.elems.size() << ")\n";

    return elems_ok ? 0 : 1;
}

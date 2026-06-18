/// Mesh refinement timing benchmark.
#include "mesh/refine.h"
#include <chrono>
#include <iostream>

using namespace lod2d;
namespace chr = std::chrono;

int main() {
    TriMesh T0;
    T0.nodes={{0,0},{1,0},{1,1},{0,1}};
    T0.elems={{0,1,3},{1,2,3}};
    T0.dirichlet={0,1,2,3};

    std::cout << "=== Mesh refinement benchmark ===\n";
    for (int level = 1; level <= 8; ++level) {
        auto t0 = chr::high_resolution_clock::now();
        auto out = refine_mesh(T0, level);
        auto t1 = chr::high_resolution_clock::now();
        double ms = chr::duration<double, std::milli>(t1 - t0).count();
        std::cout << "level=" << level
                  << " nodes=" << out.mesh.nodes.size()
                  << " elems=" << out.mesh.elems.size()
                  << " P_node_nnz=" << out.P_node.nonZeros()
                  << " P_elem_nnz=" << out.P_elem.nonZeros()
                  << " P_dg_nnz=" << out.P_dg.nonZeros()
                  << " time_ms=" << ms << "\n";
    }
    return 0;
}

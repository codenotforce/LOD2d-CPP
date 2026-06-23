/// Mesh refinement timing benchmark.
#include "mesh/refine.h"
#include <chrono>
#include <iostream>
#include <string>

using namespace lod2d;
namespace chr = std::chrono;

namespace {

template <class RefineFn>
void run_series(const std::string &name, const TriMesh &T0, RefineFn refine_fn) {
    std::cout << "-- " << name << " --\n";
    for (int level = 1; level <= 8; ++level) {
        auto t0 = chr::high_resolution_clock::now();
        auto out = refine_fn(T0, level);
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
}

} // namespace

int main() {
    TriMesh T0;
    T0.nodes = {{0,0},{1,0},{1,1},{0,1}};
    T0.elems = {{0,1,3},{1,2,3}};
    T0.dirichlet = {0,1,2,3};

    std::cout << "=== Mesh refinement benchmark ===\n";
    run_series("longest-edge bisection", T0, [](const TriMesh &mesh, int level) {
        return refine_mesh(mesh, level);
    });
    run_series("newest-vertex bisection", T0, [](const TriMesh &mesh, int level) {
        return refine_mesh_nvb(mesh, level);
    });
    return 0;
}

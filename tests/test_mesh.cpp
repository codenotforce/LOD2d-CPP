#include "mesh/refine.h"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace lod2d;

int main() {
    std::cout << "=== Mesh Refinement Tests ===\n";

    // ---- Test 1: Edge enumeration ----
    {
        TriMesh mesh;
        mesh.nodes = {{0,0}, {1,0}, {1,1}, {0,1}};
        mesh.elems = {{0,1,3}, {1,2,3}};

        auto [edges, is_bdy] = compute_edges(mesh);
        assert(edges.size() == 5);       // 5 unique edges in 2‑tri square
        int bdy_count = 0;
        for (bool b : is_bdy) if (b) ++bdy_count;
        assert(bdy_count == 4);          // 4 boundary edges
        std::cout << "  [PASS] Edge enumeration\n";
    }

    // ---- Test 2: Area computation ----
    {
        TriMesh mesh;
        mesh.nodes = {{0,0}, {1,0}, {1,1}, {0,1}};
        mesh.elems = {{0,1,3}, {1,2,3}};

        auto areas = compute_area(mesh);
        assert(areas.size() == 2);
        assert(std::abs(areas[0] - 0.5) < 1e-12);
        assert(std::abs(areas[1] - 0.5) < 1e-12);
        std::cout << "  [PASS] Area calculation\n";
    }

    // ---- Test 3: Single-level refinement ----
    {
        TriMesh mesh;
        mesh.nodes = {{0,0}, {1,0}, {1,1}, {0,1}};
        mesh.elems = {{0,1,3}, {1,2,3}};
        mesh.dirichlet = {0,1,2,3};

        RefineOutput out = refine(mesh);
        const auto &fine = out.mesh;

        assert(fine.elems.size() == 8);          // 2 → 8
        assert(fine.nodes.size() == 9);          // 4 corners + 5 edge midpoints
        assert(out.P_node.rows() == 9);
        assert(out.P_node.cols() == 4);
        assert(out.P_elem.rows() == 8);
        assert(out.P_elem.cols() == 2);
        assert(out.P_dg.rows() == 24);           // 8 elems × 3 DG dofs
        assert(out.P_dg.cols() == 6);            // 2 elems × 3 DG dofs

        // Check one prolongation entry: new node 4 (midpoint of (0,1))
        // should have weight 0.5 from nodes 0 and 1
        assert(std::abs(out.P_node.coeff(4, 0) - 0.5) < 1e-12);
        assert(std::abs(out.P_node.coeff(4, 1) - 0.5) < 1e-12);

        // Dirichlet: all 4 corners + all boundary midpoints should be Dirichlet
        assert(fine.dirichlet.size() >= 8);

        std::cout << "  [PASS] Single-level refinement\n";
    }

    // ---- Test 4: Multi-level refinement ----
    {
        TriMesh mesh;
        mesh.nodes = {{0,0}, {1,0}, {1,1}, {0,1}};
        mesh.elems = {{0,1,3}, {1,2,3}};
        mesh.dirichlet = {0,1,2,3};

        for (int nref : {2, 4, 6}) {
            RefineOutput out = refine_mesh(mesh, nref);
            int expected_elems = 2;
            for (int i = 0; i < nref; ++i) expected_elems *= 4;
            assert(out.mesh.elems.size() == static_cast<size_t>(expected_elems));

            // P_node should accumulate correctly
            assert(out.P_node.cols() == static_cast<int>(mesh.nodes.size()));
            assert(out.P_elem.cols() == static_cast<int>(mesh.elems.size()));
        }
        std::cout << "  [PASS] Multi-level refinement (2,4,6 levels)\n";
    }

    // ---- Test 5: Refinement preserves geometry (area unchanged) ----
    {
        TriMesh mesh;
        mesh.nodes = {{0,0}, {1,0}, {1,1}, {0,1}};
        mesh.elems = {{0,1,3}, {1,2,3}};

        double orig_area = 0;
        for (double a : compute_area(mesh)) orig_area += a;

        RefineOutput out = refine_mesh(mesh, 4);
        double fine_area = 0;
        for (double a : compute_area(out.mesh)) fine_area += a;

        assert(std::abs(orig_area - fine_area) < 1e-10);
        std::cout << "  [PASS] Area preserved after 4 refinements\n";
    }

    std::cout << "\nAll tests passed!\n";
    return 0;
}

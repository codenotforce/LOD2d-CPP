add_executable(bench_refine bench_refine.cpp)
target_link_libraries(bench_refine PRIVATE lod2d_core)

add_executable(bench_H4h8 bench_H4h8.cpp)
target_link_libraries(bench_H4h8 PRIVATE lod2d_core OpenMP::OpenMP_CXX)

add_executable(bench_profile bench_profile.cpp)
target_link_libraries(bench_profile PRIVATE lod2d_core OpenMP::OpenMP_CXX)

add_executable(bench_H4h9 bench_H4h9.cpp)
target_link_libraries(bench_H4h9 PRIVATE lod2d_core OpenMP::OpenMP_CXX)

add_executable(bench_reuse_rhs bench_reuse_rhs.cpp)
target_link_libraries(bench_reuse_rhs PRIVATE lod2d_core OpenMP::OpenMP_CXX)

add_executable(bench_inverse_inequality bench_inverse_inequality.cpp)
target_link_libraries(bench_inverse_inequality PRIVATE lod2d_core OpenMP::OpenMP_CXX)

add_executable(bench_saddle_h3h10 bench_saddle_h3h10.cpp)
target_link_libraries(bench_saddle_h3h10 PRIVATE lod2d_core OpenMP::OpenMP_CXX)

add_executable(bench_helmholtz_k bench_helmholtz_k.cpp)
target_link_libraries(bench_helmholtz_k PRIVATE lod2d_core)
if(LOD2D_BUILD_TESTS)
    add_test(NAME helmholtz_pollution_manufactured_smoke
        COMMAND bench_helmholtz_k
            --k=2 --H=2 --h=6 --ell=1 --source=manufactured
            --mode=two-sided --solver=schur --threads=2
            --check --format=csv)
endif()

add_executable(bench_helmholtz_manufactured bench_helmholtz_manufactured.cpp)
target_link_libraries(bench_helmholtz_manufactured PRIVATE lod2d_core)

add_executable(bench_helmholtz_visualization bench_helmholtz_visualization.cpp)
target_link_libraries(
    bench_helmholtz_visualization PRIVATE lod2d_core lod2d_visualization_io)
if(LOD2D_BUILD_TESTS)
    add_test(NAME helmholtz_visualization_export_smoke
        COMMAND bench_helmholtz_visualization
            --k=2 --H=1 --h=3 --ell=1
            --output-dir=${CMAKE_CURRENT_BINARY_DIR}/visualization-smoke)
endif()

add_executable(bench_helmholtz_H_convergence bench_helmholtz_H_convergence.cpp)
target_link_libraries(bench_helmholtz_H_convergence PRIVATE lod2d_core)
if(LOD2D_BUILD_TESTS)
    add_test(NAME helmholtz_H_convergence_smoke
        COMMAND bench_helmholtz_H_convergence
            --k=2 --h=6 --H-levels=2,3 --ell=1
            --solver=schur --threads=2 --check --format=csv)
endif()

add_executable(bench_helmholtz_hp_convergence bench_helmholtz_hp_convergence.cpp)
target_link_libraries(bench_helmholtz_hp_convergence PRIVATE lod2d_core)

add_executable(bench_helmholtz_local_inverse bench_helmholtz_local_inverse.cpp)
target_link_libraries(bench_helmholtz_local_inverse PRIVATE lod2d_core)
if(LOD2D_BUILD_TESTS)
    add_test(NAME helmholtz_local_inverse_smoke
        COMMAND bench_helmholtz_local_inverse
            --k=2 --initial-H=2 --fine-level=7 --steps=2 --ell=2
            --mark=single-chain --check --format=csv)
    add_test(NAME helmholtz_local_inverse_server_smoke
        COMMAND bench_helmholtz_local_inverse
            --k=2 --initial-H=2 --fine-level=7 --steps=1 --ell=2
            --mark=argmax-patch --basis=trial --denominators=matched
            --check --format=csv)
endif()

add_executable(bench_helmholtz_adaptive bench_helmholtz_adaptive.cpp)
target_link_libraries(bench_helmholtz_adaptive PRIVATE lod2d_core)

add_executable(bench_helmholtz_certified bench_helmholtz_certified.cpp)
target_link_libraries(bench_helmholtz_certified PRIVATE lod2d_core)
if(LOD2D_BUILD_TESTS)
    add_test(NAME helmholtz_certified_backend_strict_smoke
        COMMAND bench_helmholtz_certified --evidence=strict --check)
    add_test(NAME helmholtz_certified_backend_conditional_smoke
        COMMAND bench_helmholtz_certified --evidence=conditional --check)
    set_tests_properties(
        helmholtz_certified_backend_strict_smoke
        helmholtz_certified_backend_conditional_smoke
        PROPERTIES TIMEOUT 60 LABELS "smoke;helmholtz;conditional")
endif()

add_executable(bench_helmholtz_patch_solver bench_helmholtz_patch_solver.cpp)
target_link_libraries(bench_helmholtz_patch_solver PRIVATE lod2d_core)

add_executable(bench_helmholtz_coarse_schwarz bench_helmholtz_coarse_schwarz.cpp)
target_link_libraries(bench_helmholtz_coarse_schwarz PRIVATE lod2d_core)

add_executable(bench_helmholtz_two_level_schwarz bench_helmholtz_two_level_schwarz.cpp)
target_link_libraries(bench_helmholtz_two_level_schwarz PRIVATE lod2d_core)

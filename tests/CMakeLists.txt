add_executable(test_mesh test_mesh.cpp)
target_link_libraries(test_mesh PRIVATE lod2d_core)
add_test(NAME mesh_refinement COMMAND test_mesh)

add_executable(test_golden test_golden.cpp)
target_link_libraries(test_golden PRIVATE lod2d_core)

add_executable(test_dg test_dg.cpp)
target_link_libraries(test_dg PRIVATE lod2d_core)
add_test(NAME dg_assembly COMMAND test_dg)

add_executable(test_nvb test_nvb.cpp)
target_link_libraries(test_nvb PRIVATE lod2d_core)
add_test(NAME newest_vertex_bisection COMMAND test_nvb)

add_executable(test_qi test_qi.cpp)
target_link_libraries(test_qi PRIVATE lod2d_core)
add_test(NAME quasi_interp COMMAND test_qi)

add_executable(test_patch test_patch.cpp)
target_link_libraries(test_patch PRIVATE lod2d_core)
add_test(NAME patches COMMAND test_patch)

add_executable(test_corr test_corr.cpp)
target_link_libraries(test_corr PRIVATE lod2d_core)

add_executable(test_full test_full.cpp)
target_link_libraries(test_full PRIVATE lod2d_core)
target_link_libraries(test_full PRIVATE OpenMP::OpenMP_CXX)

option(LOD2D_ENABLE_MATLAB_GOLDEN_TESTS "Register legacy MATLAB red-refinement golden tests with CTest" OFF)
if(LOD2D_ENABLE_MATLAB_GOLDEN_TESTS)
    add_test(NAME golden_mesh COMMAND test_golden)
    add_test(NAME corrector COMMAND test_corr)
    add_test(NAME full_lod COMMAND test_full)
endif()

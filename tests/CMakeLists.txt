add_executable(test_mesh test_mesh.cpp)
target_link_libraries(test_mesh PRIVATE lod2d_core)
add_test(NAME mesh_refinement COMMAND test_mesh)

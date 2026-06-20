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

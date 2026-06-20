cmake_minimum_required(VERSION 3.16)
project(LOD2d VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---- Build options ----
option(LOD2D_USE_OPENMP "Enable OpenMP parallelism" ON)
option(LOD2D_BUILD_TESTS "Build test suite" ON)
option(LOD2D_BUILD_BENCHMARKS "Build benchmarks" ON)

# ---- Dependencies ----
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(OpenMP)
find_package(SuiteSparse QUIET)
if(NOT SuiteSparse_FOUND)
    # Fallback: find CHOLMOD directly
    find_library(CHOLMOD_LIB cholmod)
    find_path(CHOLMOD_INC suitesparse/cholmod.h PATH_SUFFIXES suitesparse)
    if(CHOLMOD_LIB AND CHOLMOD_INC)
        add_library(CHOLMOD::CHOLMOD UNKNOWN IMPORTED)
        set_target_properties(CHOLMOD::CHOLMOD PROPERTIES
            IMPORTED_LOCATION ${CHOLMOD_LIB}
            INTERFACE_INCLUDE_DIRECTORIES ${CHOLMOD_INC})
    endif()
endif()
find_package(TBB QUIET)

# ---- Compiler flags ----
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-O3 -march=native -DNDEBUG)
    if(LOD2D_USE_OPENMP AND OpenMP_FOUND)
        add_compile_options(-fopenmp)
    endif()
else()
    add_compile_options(-g -O0 -Wall -Wextra -Wpedantic)
endif()

# ---- Library: lod2d_core ----
add_library(lod2d_core STATIC
    src/mesh/refine.cpp
    src/mesh/edges.cpp
    src/fem/assemble_dg.cpp
    src/lod/quasi_interp.cpp
    src/lod/patches.cpp
    src/lod/corrector.cpp
    src/lod/lod_model.cpp
    src/solver/coarse_solve.cpp
    src/solver/lod_reuse.cpp
)
target_include_directories(lod2d_core PUBLIC include)
target_link_libraries(lod2d_core PUBLIC Eigen3::Eigen)
if(SuiteSparse_FOUND)
    target_link_libraries(lod2d_core PUBLIC
        SuiteSparse::CHOLMOD SuiteSparse::AMD SuiteSparse::COLAMD)
elseif(CHOLMOD_LIB)
    target_link_libraries(lod2d_core PUBLIC CHOLMOD::CHOLMOD)
endif()
if(LOD2D_USE_OPENMP AND OpenMP_FOUND)
    target_link_libraries(lod2d_core PUBLIC OpenMP::OpenMP_CXX)
endif()
if(TBB_FOUND)
    target_link_libraries(lod2d_core PUBLIC TBB::tbb)
endif()

# ---- Main executable ----
add_executable(lod2d src/main.cpp)
target_link_libraries(lod2d PRIVATE lod2d_core)

# ---- Tests ----
if(LOD2D_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# ---- Benchmarks ----
if(LOD2D_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

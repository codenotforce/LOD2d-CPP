cmake_minimum_required(VERSION 3.16)
project(LOD2d VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---- Build options ----
option(LOD2D_USE_OPENMP "Enable OpenMP parallelism" ON)
option(LOD2D_BUILD_TESTS "Build test suite" ON)
option(LOD2D_BUILD_BENCHMARKS "Build benchmarks" ON)
option(LOD_ENABLE_VERIFIED_CERTIFICATES
    "Enable MPFR/MPFI directed-rounding certificate verification" OFF)

# ---- Dependencies ----
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(OpenMP)
find_package(SuiteSparse QUIET)
find_library(UMFPACK_LIBRARY NAMES umfpack)
find_path(UMFPACK_INCLUDE_DIR umfpack.h PATH_SUFFIXES suitesparse)
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

if(LOD_ENABLE_VERIFIED_CERTIFICATES)
    find_path(MPFR_INCLUDE_DIR mpfr.h)
    find_path(MPFI_INCLUDE_DIR mpfi.h)
    find_library(MPFR_LIBRARY NAMES mpfr)
    find_library(MPFI_LIBRARY NAMES mpfi)
    find_library(GMP_LIBRARY NAMES gmp)
    if(NOT MPFR_INCLUDE_DIR OR NOT MPFI_INCLUDE_DIR OR NOT MPFR_LIBRARY OR
       NOT MPFI_LIBRARY OR NOT GMP_LIBRARY)
        message(FATAL_ERROR
            "LOD_ENABLE_VERIFIED_CERTIFICATES requires MPFR, MPFI, and GMP")
    endif()
endif()

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
    src/fem/hp_triangle.cpp
    src/lod/quasi_interp.cpp
    src/helmholtz/manufactured.cpp
    src/helmholtz/boundary.cpp
    src/helmholtz/quadrature.cpp
    src/lod/patches.cpp
    src/lod/corrector.cpp
    src/lod/lod_model.cpp
    src/helmholtz/patch_system.cpp
    src/helmholtz/patch_solver.cpp
    src/helmholtz/corrector_pipeline.cpp
    src/helmholtz/coarse_schwarz.cpp
    src/helmholtz/schwarz_patch.cpp
    src/helmholtz/patch_hierarchy.cpp
    src/helmholtz/schwarz_local_solver.cpp
    src/helmholtz/two_level_schwarz.cpp
    src/helmholtz/shifted_laplacian.cpp
    src/helmholtz/patch_multigrid.cpp
    src/helmholtz/operators.cpp
    src/helmholtz/hp_operators.cpp
    src/helmholtz/hp_fem.cpp
    src/helmholtz/hp_interpolation.cpp
    src/helmholtz/hp_patch.cpp
    src/helmholtz/hp_model.cpp
    src/helmholtz/corrector.cpp
    src/helmholtz/model.cpp
    src/helmholtz/adaptive/hierarchy.cpp
    src/helmholtz/adaptive/error_control.cpp
    src/helmholtz/adaptive/kernel_residual.cpp
    src/helmholtz/adaptive/candidate_flux.cpp
    src/helmholtz/adaptive/candidate_dual.cpp
    src/helmholtz/adaptive/reference_epoch_driver.cpp
    src/helmholtz/adaptive/singularity_hybrid.cpp
    src/helmholtz/adaptive/reference_retraction.cpp
    src/helmholtz/adaptive/practical_driver.cpp
    src/helmholtz/adaptive/verified_spectrum.cpp
    src/helmholtz/adaptive/certificates.cpp
    src/helmholtz/adaptive/certified_driver.cpp
    src/helmholtz/adaptive/numerical_backend.cpp
    src/helmholtz/adaptive/estimator.cpp
    src/helmholtz/adaptive/reliability.cpp
    src/helmholtz/adaptive/driver.cpp
    src/helmholtz/experiments/paper_config.cpp
    src/helmholtz/experiments/reference_solution_cache.cpp
    src/helmholtz/benchmarks/paper_cases.cpp
    src/solver/coarse_solve.cpp
    src/solver/right_gmres.cpp
    src/solver/lod_reuse.cpp
)
target_include_directories(lod2d_core PUBLIC include)
target_link_libraries(lod2d_core PUBLIC Eigen3::Eigen)
if(UMFPACK_LIBRARY AND UMFPACK_INCLUDE_DIR)
    target_compile_definitions(lod2d_core PUBLIC LOD2D_HAVE_UMFPACK=1)
    target_include_directories(lod2d_core PUBLIC ${UMFPACK_INCLUDE_DIR})
    target_link_libraries(lod2d_core PUBLIC ${UMFPACK_LIBRARY})
endif()
if(LOD_ENABLE_VERIFIED_CERTIFICATES)
    target_compile_definitions(lod2d_core PUBLIC LOD_ENABLE_VERIFIED_CERTIFICATES=1)
    target_include_directories(lod2d_core PRIVATE
        ${MPFR_INCLUDE_DIR} ${MPFI_INCLUDE_DIR})
    target_link_libraries(lod2d_core PUBLIC
        ${MPFI_LIBRARY} ${MPFR_LIBRARY} ${GMP_LIBRARY})
endif()
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

# ---- Optional-use visualization I/O (not linked into lod2d_core) ----
add_library(lod2d_visualization_io STATIC
    src/io/vtk_writer.cpp
)
target_include_directories(lod2d_visualization_io PUBLIC include)
target_link_libraries(lod2d_visualization_io PUBLIC Eigen3::Eigen)

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

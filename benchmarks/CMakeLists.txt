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

add_executable(bench_helmholtz_adaptive_paper
    bench_helmholtz_adaptive_paper.cpp
    reference_epoch_paper_runner.cpp)
target_link_libraries(
    bench_helmholtz_adaptive_paper PRIVATE lod2d_core lod2d_visualization_io)
if(LOD2D_BUILD_TESTS)
    add_test(NAME helmholtz_adaptive_paper_v4_smoke
        COMMAND bench_helmholtz_adaptive_paper
            --config=${CMAKE_SOURCE_DIR}/experiments/helmholtz_adaptive_paper/configs/R1-palod-smoke-v4.json
            --output-dir=${CMAKE_CURRENT_BINARY_DIR}/adaptive-paper-v4-smoke
            --manuscript-baseline=${CMAKE_SOURCE_DIR}/experiments/helmholtz_adaptive_paper/MANUSCRIPT_BASELINE_LEGACY_V4.sha256
            --check)
    add_test(NAME helmholtz_reference_epoch_paper_v5_smoke
        COMMAND bench_helmholtz_adaptive_paper
            --config=${CMAKE_SOURCE_DIR}/experiments/helmholtz_adaptive_paper/configs/R1-reference-epoch-smoke-v5.json
            --output-dir=${CMAKE_CURRENT_BINARY_DIR}/reference-epoch-paper-v5-smoke
            --manuscript-baseline=${CMAKE_SOURCE_DIR}/experiments/helmholtz_adaptive_paper/MANUSCRIPT_BASELINE.sha256
            --check)
    add_test(NAME helmholtz_hlod_fixed_paper_v4_smoke
        COMMAND bench_helmholtz_adaptive_paper
            --config=${CMAKE_SOURCE_DIR}/experiments/helmholtz_adaptive_paper/configs/R1-hlod-fixed-smoke-v4.json
            --output-dir=${CMAKE_CURRENT_BINARY_DIR}/hlod-fixed-paper-v4-smoke
            --manuscript-baseline=${CMAKE_SOURCE_DIR}/experiments/helmholtz_adaptive_paper/MANUSCRIPT_BASELINE_LEGACY_V4.sha256
            --check)
    add_test(NAME helmholtz_palod_fixed_horizon_v4_smoke
        COMMAND bench_helmholtz_adaptive_paper
            --config=${CMAKE_SOURCE_DIR}/experiments/helmholtz_adaptive_paper/configs/R1-palod-fixed-horizon-smoke-v4.json
            --output-dir=${CMAKE_CURRENT_BINARY_DIR}/palod-fixed-horizon-v4-smoke
            --manuscript-baseline=${CMAKE_SOURCE_DIR}/experiments/helmholtz_adaptive_paper/MANUSCRIPT_BASELINE_LEGACY_V4.sha256
            --check)
    add_test(NAME helmholtz_ufem_paper_v4_smoke
        COMMAND bench_helmholtz_adaptive_paper
            --config=${CMAKE_SOURCE_DIR}/experiments/helmholtz_adaptive_paper/configs/R1-ufem-smoke-v4.json
            --output-dir=${CMAKE_CURRENT_BINARY_DIR}/ufem-paper-v4-smoke
            --manuscript-baseline=${CMAKE_SOURCE_DIR}/experiments/helmholtz_adaptive_paper/MANUSCRIPT_BASELINE_LEGACY_V4.sha256
            --check)
    add_test(NAME helmholtz_slod_paper_v4_smoke
        COMMAND bench_helmholtz_adaptive_paper
            --config=${CMAKE_SOURCE_DIR}/experiments/helmholtz_adaptive_paper/configs/R1-slod-smoke-v4.json
            --output-dir=${CMAKE_CURRENT_BINARY_DIR}/slod-paper-v4-smoke
            --manuscript-baseline=${CMAKE_SOURCE_DIR}/experiments/helmholtz_adaptive_paper/MANUSCRIPT_BASELINE_LEGACY_V4.sha256
            --check)
    add_test(NAME helmholtz_afem_paper_v4_smoke
        COMMAND bench_helmholtz_adaptive_paper
            --config=${CMAKE_SOURCE_DIR}/experiments/helmholtz_adaptive_paper/configs/R1-afem-smoke-v4.json
            --output-dir=${CMAKE_CURRENT_BINARY_DIR}/afem-paper-v4-smoke
            --manuscript-baseline=${CMAKE_SOURCE_DIR}/experiments/helmholtz_adaptive_paper/MANUSCRIPT_BASELINE_LEGACY_V4.sha256
            --check)
    add_test(NAME helmholtz_afem_exact_only_reference_exhaustion_v4_smoke
        COMMAND bench_helmholtz_adaptive_paper
            --config=${CMAKE_SOURCE_DIR}/experiments/helmholtz_adaptive_paper/configs/S-afem-exact-only-reference-exhaustion-smoke-v4.json
            --output-dir=${CMAKE_CURRENT_BINARY_DIR}/afem-exact-only-reference-exhaustion-v4-smoke
            --manuscript-baseline=${CMAKE_SOURCE_DIR}/experiments/helmholtz_adaptive_paper/MANUSCRIPT_BASELINE_LEGACY_V4.sha256
            --check)
    set_tests_properties(
        helmholtz_adaptive_paper_v4_smoke
        helmholtz_hlod_fixed_paper_v4_smoke
        helmholtz_palod_fixed_horizon_v4_smoke
        helmholtz_ufem_paper_v4_smoke
        helmholtz_slod_paper_v4_smoke
        helmholtz_afem_paper_v4_smoke
        helmholtz_afem_exact_only_reference_exhaustion_v4_smoke
        PROPERTIES TIMEOUT 60 LABELS "smoke;helmholtz;practical-paper-v4")
    set_tests_properties(helmholtz_reference_epoch_paper_v5_smoke
        PROPERTIES TIMEOUT 300 LABELS "smoke;helmholtz;reference-epoch-paper-v5")
    set_tests_properties(
        helmholtz_slod_paper_v4_smoke PROPERTIES TIMEOUT 120)
endif()

add_executable(bench_helmholtz_reference_adequacy
    bench_helmholtz_reference_adequacy.cpp)
target_link_libraries(bench_helmholtz_reference_adequacy PRIVATE lod2d_core)
if(LOD2D_BUILD_TESTS)
    add_test(NAME helmholtz_reference_adequacy_v4_smoke
        COMMAND bench_helmholtz_reference_adequacy
            --config=${CMAKE_SOURCE_DIR}/experiments/helmholtz_adaptive_paper/configs/R1-reference-adequacy-smoke-v4.json
            --source-iterations=${CMAKE_SOURCE_DIR}/experiments/helmholtz_adaptive_paper/fixtures/R1-reference-audit-source-v3.csv
            --output-dir=${CMAKE_CURRENT_BINARY_DIR}/reference-adequacy-v4-smoke
            --reference-cache-dir=${CMAKE_CURRENT_BINARY_DIR}/reference-adequacy-v4-cache
            --manuscript-baseline=${CMAKE_SOURCE_DIR}/experiments/helmholtz_adaptive_paper/MANUSCRIPT_BASELINE_LEGACY_V4.sha256
            --check)
    set_tests_properties(helmholtz_reference_adequacy_v4_smoke
        PROPERTIES TIMEOUT 60 LABELS "smoke;helmholtz;practical-paper-v4")
endif()

add_executable(bench_helmholtz_e0_calibration
    bench_helmholtz_e0_calibration.cpp)
target_link_libraries(bench_helmholtz_e0_calibration PRIVATE lod2d_core)
if(LOD2D_BUILD_TESTS)
    add_test(NAME helmholtz_e0_R1_k16_calibration
        COMMAND bench_helmholtz_e0_calibration
            --output-dir=${CMAKE_CURRENT_BINARY_DIR}/e0-R1-k16-calibration
            --check)
    set_tests_properties(helmholtz_e0_R1_k16_calibration
        PROPERTIES TIMEOUT 420 LABELS "calibration;helmholtz;practical-paper-v2")
endif()

add_executable(bench_helmholtz_reference_epoch_e0
    bench_helmholtz_reference_epoch_e0.cpp)
target_link_libraries(bench_helmholtz_reference_epoch_e0 PRIVATE lod2d_core)
if(LOD2D_BUILD_TESTS)
    add_test(NAME helmholtz_reference_epoch_e0
        COMMAND bench_helmholtz_reference_epoch_e0
            --output-dir=${CMAKE_CURRENT_BINARY_DIR}/reference-epoch-e0
            --check)
    set_tests_properties(helmholtz_reference_epoch_e0
        PROPERTIES TIMEOUT 420 LABELS "calibration;helmholtz;reference-epoch-v5")
endif()

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

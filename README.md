# LOD2d-C++

High-performance C++20 implementation of the Localized Orthogonal
Decomposition (LOD) method for 2D elliptic diffusion problems.

This project is a MATLAB-to-C++ port of the LOD2d workflow.  The mesh layer now
provides both conforming longest-edge bisection and newest-vertex bisection
(NVB), so it can serve later adaptive LOD experiments.  The old MATLAB
red-refinement golden executables are still buildable as references, but they
are not part of default CTest because the mesh topology and numbering
intentionally differ.

## Dependencies

| Library | Purpose |
|---------|---------|
| Eigen 3 | Dense and sparse linear algebra |
| SuiteSparse / CHOLMOD | Experimental sparse Cholesky backend |
| OpenMP | Parallel corrector loop |
| CMake >= 3.16 | Build system |

Ubuntu / WSL setup:

```bash
sudo apt update
sudo apt install -y build-essential cmake g++ libeigen3-dev libsuitesparse-dev libtbb-dev
```

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

Build options:

| Option | Default | Description |
|--------|---------|-------------|
| `LOD2D_USE_OPENMP` | `ON` | Enable OpenMP parallelism |
| `LOD2D_BUILD_TESTS` | `ON` | Build test suite |
| `LOD2D_BUILD_BENCHMARKS` | `ON` | Build benchmarks |

## Run

Run focused tests from the repository root:

```bash
ctest --test-dir build --output-on-failure
./build/tests/test_mesh
./build/tests/test_dg
./build/tests/test_nvb
./build/tests/test_patch
./build/tests/test_qi
```

Run benchmarks:

```bash
./build/benchmarks/bench_refine
./build/benchmarks/bench_saddle_h3h10 --H=3 --h=10 --ell=3 --threads=8 --skip-reference
./build/benchmarks/bench_H4h9 --solver=eigen
./build/benchmarks/bench_H4h9 --solver=cholmod
./build/benchmarks/bench_H4h9 --solver=cholmod_cached --skip-reference
./build/benchmarks/bench_profile --solver=auto --skip-reference
./build/benchmarks/bench_reuse_rhs --solver=auto --rhs=5
./build/benchmarks/bench_inverse_inequality --sweep-H --H-min=2 --H-max=4 --h-minus-H=5 --ell=2 --coeff=unit --solver=eigen
./build/benchmarks/bench_inverse_inequality --basis=coarse --sweep-H --H-min=2 --H-max=5 --h-minus-H=5 --space=all --coeff=unit --solver=eigen
./build/benchmarks/bench_inverse_inequality --sweep-h --H=3 --h-min=6 --h-max=9 --ell=2 --basis=lod --coeff=unit --solver=eigen
./build/benchmarks/bench_inverse_inequality --H=3 --h=9 --ell=2 --numerator=corrector --coeff=unit --solver=eigen
```

`--solver=eigen` is the default for the small tested corrector sizes.
`--solver=auto` selects CHOLMOD for h>=10 benchmarks and Eigen otherwise.
`--solver=cholmod_cached` is a correctness-checked experimental path for
thread-local CHOLMOD symbolic factor reuse; it is not the default because the
current patch order does not produce enough repeated patterns to offset cache
overhead.


## High-Level API

For repeated solves with fixed `A`, mesh, `H/h`, and `ell` but different right-hand sides, use `LodModel` instead of rebuilding the correctors in each driver:

```cpp
#include "lod/lod_model.h"

using namespace lod2d;

LodProblemConfig config;
config.H = 4;
config.h = 10;
config.ell = 2;
config.solver = CorrectorSolver::Cholmod;
config.initial_mesh = make_unit_square_mesh();

LodModel model = LodModel::build(config, Ah);
LodReuseSolution sol = model.solve_from_coarse_values(f_coarse);
```

`LodModel` owns the reusable multiscale system and cached coarse factorization. By default it releases setup-only `P_elem` and `P_dg` after construction; set `keep_setup_matrices=true` only when inspecting those matrices.
Lower-level helpers such as `build_lod_problem_data`, `build_lod_operators`,
`build_lod_correctors`, and `build_lod_basis` remain available for benchmarks
that need phase timings or memory profiling.


## Server Inverse-Inequality Runs

For high-memory servers, use the helper script to probe the continuous-corrector
limit by fixing `H` and refining only `h`:

```bash
git clone https://github.com/codenotforce/LOD2d-CPP.git
cd LOD2d-CPP
sudo apt update
sudo apt install -y build-essential cmake g++ libeigen3-dev libsuitesparse-dev libtbb-dev
THREADS=8 MODE=h H=3 ELL=2 H_MIN=6 H_MAX=12 COEFF=unit SOLVER=auto ./scripts/run_inverse_server.sh
THREADS=8 MODE=H H_FIXED=10 ELL=2 H_MIN=2 H_MAX=5 NUMERATOR=corrector COEFF=unit SOLVER=auto ./scripts/run_inverse_server.sh
```

The script builds `bench_inverse_inequality`, runs each `h` separately, and
stores logs plus `/usr/bin/time -v` peak-memory output under
`results/inverse_inequality/`.

On a 377 GiB EPYC server, `H=3,h=11` and `h=12` should be feasible with the
current full-global benchmark. `h=13` is experimental and may be very slow;
`h>=14` is not recommended without a local/streaming implementation.
## Project Status

| Module | Status | Test coverage |
|--------|--------|---------------|
| Mesh refinement and prolongation | LEB + NVB | Conformity, area, prolongation, and MATLAB-derived NVB golden tests |
| DG stiffness assembly | Complete | `test_dg`: 10/10 pass |
| Quasi-interpolation | Complete | LEB structural smoke tests |
| Patch construction | Complete | Golden patch tests |
| Element corrector | Complete | `test_corr`: Eigen, CHOLMOD, and cached CHOLMOD pass |
| Full LOD pipeline | Complete | `test_full`: 3/3 pass |
| H4/h9 benchmark | Complete | Error check vs MATLAB reference |

## Recent Optimization Work

The current corrector path avoids repeated global sparse work and unnecessary
local allocations inside each element corrector.  Shared helpers now keep the
benchmark and test pipelines consistent:

1. `compute_all_correctors` centralizes the OpenMP corrector loop and solver
   selection.
2. `build_multiscale_basis` centralizes `G = P_node - C_ell` assembly from
   compact corrector entries.
3. `ElementStiffnessBlocks`: stores each fine element's local 3x3 stiffness
   block while assembling DG stiffness, so correctors no longer repeatedly scan
   `Shdg` columns to recover the same block.
4. `FineElementChildren`: stores the `coarse element -> fine children` mapping
   from `P_elem`, so correctors no longer compute `P0 * patch(:, k)` and scan
   all fine elements for every coarse element.
5. Multi-RHS local solves: the Eigen corrector now calls `llt.solve(RHS)` once
   instead of solving each RHS column separately after the same factorization.
6. Thread-local scratch buffers: large `Nh`-sized marker arrays are reused per
   OpenMP thread instead of being allocated and cleared for every corrector.
7. Local dense `IHp` and triplet `CTk` construction avoid many tiny sparse
   insertions in the hot loop.
8. `cholmod_cached` keeps a bounded thread-local CHOLMOD symbolic factor cache
   for reproducible experiments.  The cache is intentionally capped at one
   pattern per thread after an unlimited cache reached roughly 11.6 GB RSS and
   was killed on the 12 GB WSL test machine.
9. `LodReusableSystem` caches `G`, `G0`, and the coarse LOD factorization for
   repeated solves with the same `A`, mesh, `H/h`, and `ell` but different RHS
   values.

Historical timing tables and solver experiments are recorded in
`DEVELOPMENT.md`; README intentionally contains no benchmark result snapshots.

## Important Migration Notes

- The active mesh layer supports conforming longest-edge bisection and
  newest-vertex bisection.  Edge midpoint creation is hash-based and shared by
  adjacent elements, so local marking does not introduce hanging nodes.
- NVB uses the MATLAB convention from `lod.bisect`: local vertex 0 is the
  newest vertex and local edge (1,2) is the reference edge.
- Legacy MATLAB red-refinement golden files remain useful only for historical
  comparison; regenerate them before using `test_corr` or `test_full` as
  correctness gates under the bisection meshes.
- `P_dg` rows are grouped by sub-triangle type, not by parent element.
- Eigen keeps explicit zero triplets, so assembly skips zero values.
- Release builds are required for meaningful performance comparisons.
- The Windows tree is used as the GitHub upload mirror; WSL remains the tested
  build environment.

## License

Research and educational use.  Based on the LOD code from
*An Introduction to the Localized Orthogonal Decomposition Method*
by A. Malqvist and D. Peterseim.

## Helmholtz Petrov-Galerkin LOD

See `HELMHOLTZ_GUIDE.md` for the mathematical conventions, build and test
commands, benchmark options, output fields, server workflow, and reuse rules.

The first seven stages of `HELMHOLTZ_LOD_PLAN.md` are implemented in the
`lod2d::helmholtz` namespace:

- complex P1 Helmholtz FEM with homogeneous impedance Robin data;
- globally nested NVB coarse/fine meshes;
- complex-linear `I_H = E_H * Pi_H^dg` using the existing real sparse matrix;
- primal and adjoint localized correctors solved by complex sparse LU;
- two-sided Petrov-Galerkin and corrected-test-only model modes;
- cached coarse factorization for repeated right-hand sides;
- parallel local correctors with thread-local assembly workspaces and exact
  sparse-pattern symbolic-analysis reuse.

The assembled fine operator is

```text
A_h(k) = K_h - k^2 M_h - i k B_h,
```

with the first argument linear and the second argument conjugate-linear.
Physical-domain portions of a corrector patch retain the Robin term; artificial
patch boundaries are homogeneous Dirichlet boundaries.  The corrector and
coarse systems are not solved with LLT, CG, or PCG.

Minimal API example:

```cpp
#include "helmholtz/model.h"

using namespace lod2d::helmholtz;

HelmholtzProblemConfig config;
config.H = 1;
config.h = 4;
config.ell = 1;
config.wavenumber = 2.0;
config.mode = HelmholtzPetrovMode::TwoSided;

HelmholtzLodModel model = HelmholtzLodModel::build(config);
auto solution = model.solve_source([](const lod2d::Point2 &x) {
    return Complex(1.0 + x.x(), -x.y());
});
```

Build and run the correctness gates with:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 8
ctest --test-dir build -R helmholtz_ --output-on-failure
```

The default unit-square mesh uses NVB-compatible local vertex ordering: both
initial triangles use the shared diagonal as their reference edge.  Custom
initial meshes must provide a compatible reference-edge labeling.

### Wave-Number Scan

Run one wave number directly:

```bash
./build/benchmarks/bench_helmholtz_k \
  --k=16 --H=auto --h=auto --kH-target=1 \
  --fine-gap=8 --ell=auto --mode=two-sided
```

Run independent processes for the full sequence so memory is released after
each wave number:

```bash
THREADS=32 K_VALUES="4 8 16 32 64" FINE_GAP=8 KH_TARGET=1 \
  bash scripts/run_helmholtz_k_scan.sh
```

The resumable scan compares standard coarse P1 FEM and Petrov-Galerkin LOD
against the same fine P1 reference solution. It writes errors, `kH`, `kh`,
corrector residuals, an energy-scaled inf-sup value for small coarse systems,
phase timings, metadata, and `/usr/bin/time -v` logs to
`results/helmholtz_k_scan/`. Set `RESUME=0` to rerun completed points.

Measured result tables belong in `DEVELOPMENT.md`, not README.

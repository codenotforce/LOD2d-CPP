# LOD2d-C++

High-performance C++20 implementation of the Localized Orthogonal
Decomposition (LOD) method for 2D elliptic diffusion problems.

This project is a MATLAB-to-C++ port of the LOD2d workflow.  The current
implementation preserves the MATLAB reference ordering and numerical results,
then optimizes the expensive element-corrector phase with Eigen sparse linear
algebra and OpenMP.

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
./build/tests/test_dg
./build/tests/test_corr --solver=both
./build/tests/test_full
```

Run benchmarks:

```bash
./build/benchmarks/bench_refine
./build/benchmarks/bench_H4h8 --solver=eigen
./build/benchmarks/bench_H4h8 --solver=cholmod
./build/benchmarks/bench_H4h9 --solver=eigen
./build/benchmarks/bench_H4h9 --solver=cholmod
./build/benchmarks/bench_profile
```

`--solver=eigen` is the default and fastest path for the tested LOD corrector
sizes.  `--solver=cholmod` is kept as an experimental correctness-checked
backend.

## Project Status

| Module | Status | Test coverage |
|--------|--------|---------------|
| Mesh refinement and prolongation | Complete | Golden tests vs MATLAB |
| DG stiffness assembly | Complete | `test_dg`: 10/10 pass |
| Quasi-interpolation | Complete | Golden positions vs MATLAB |
| Patch construction | Complete | Golden patch tests |
| Element corrector | Complete | `test_corr`: Eigen and CHOLMOD pass |
| Full LOD pipeline | Complete | `test_full`: 3/3 pass |
| H4/h8 benchmark | Complete | Error check vs MATLAB reference |
| H4/h9 benchmark | Complete | Error check vs MATLAB reference |

## Recent Optimization Work

The current corrector path avoids repeated global sparse work and unnecessary
local allocations inside each element corrector:

1. `ElementStiffnessBlocks`: stores each fine element's local 3x3 stiffness
   block while assembling DG stiffness, so correctors no longer repeatedly scan
   `Shdg` columns to recover the same block.
2. `FineElementChildren`: stores the `coarse element -> fine children` mapping
   from `P_elem`, so correctors no longer compute `P0 * patch(:, k)` and scan
   all fine elements for every coarse element.
3. Multi-RHS local solves: the Eigen corrector now calls `llt.solve(RHS)` once
   instead of solving each RHS column separately after the same factorization.
4. Thread-local scratch buffers: large `Nh`-sized marker arrays are reused per
   OpenMP thread instead of being allocated and cleared for every corrector.
5. Local dense `IHp` and triplet `CTk` construction avoid many tiny sparse
   insertions in the hot loop.

### Benchmarks (H=4, h=9, ℓ=2, 7-run median, OMP_NUM_THREADS=16)

| Version | Median | Range | vs Serial |
|---------|--------|-------|-----------|
| **C++ 16t** | **9.79 s** | 9.48–10.01 s | **5.0×** |
| MATLAB parallel (4w) | 23.57 s | 22.90–24.07 s | 2.1× |
| MATLAB serial | 49.39 s | 48.19–52.17 s | 1.0× |

Errors identical across all versions:
```text
Energy: 0.0257411   L2: 0.00175159   FE-L2: 0.0157131
```

### Benchmarks (H=4, h=8, ℓ=2, 20-run median, OMP_NUM_THREADS=16)

| Version | Median | Range |
|---------|--------|-------|
| **C++ 16t** | **1.73 s** | 1.45–1.99 s |
| MATLAB parallel (4w) | 3.37 s | (hot cache) |
| MATLAB serial | 27.0 s | — |

Errors:
```text
Energy: 0.0247816   L2: 0.00166544   FE-L2: 0.0145842
```

CHOLMOD is now stable and exact against golden corrector data, but remains
slower than Eigen for H4/h8 correctors, so it is not the default.

## Important Migration Notes

- MATLAB sparse construction/order matters.  Edge ordering and refined triangle
  ordering are intentionally matched to MATLAB golden data.
- `P_dg` rows are grouped by sub-triangle type, not by parent element.
- Eigen keeps explicit zero triplets, so assembly skips zero values.
- Release builds are required for meaningful performance comparisons.
- The Windows tree is used as the GitHub upload mirror; WSL remains the tested
  build environment.

## License

Research and educational use.  Based on the LOD code from
*An Introduction to the Localized Orthogonal Decomposition Method*
by A. Malqvist and D. Peterseim.

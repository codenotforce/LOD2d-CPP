# LOD2d-C++

High‑performance C++20 implementation of the Localized Orthogonal
Decomposition (LOD) method for 2D elliptic diffusion problems.
Ported from [LOD2d-MATLAB](https://github.com/codenotforce/LOD2d-MATLAB).

[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue)](https://en.cppreference.com/w/cpp/20)
[![CMake](https://img.shields.io/badge/CMake-3.16%2B-green)](https://cmake.org)
[![License](https://img.shields.io/badge/license-research%20%2B%20educational-lightgrey)](LICENSE)

## Overview

LOD2d solves the linear elliptic diffusion problem

```
−∇·(A∇u) = f   on Ω ⊂ ℝ²,   u = 0 on ∂Ω
```

with a highly oscillatory coefficient A, using the Localized Orthogonal
Decomposition method on triangular meshes.  The key computational kernel
— element correctors — is parallelised with OpenMP.

### Why C++?

| | MATLAB (optimised) | C++ (target) |
|---|---|---|
| Corrector speedup | 2.8× (`parfor` + RCM) | **>5×** |
| Memory | GC‑managed, copies of large sparse matrices | Explicit, zero‑copy shared |
| Parallel efficiency | ~93% (Threads `parfor`) | **~98%** (OpenMP) |
| Deployment | Requires MATLAB + PCT | Single static binary |

## Dependencies

| Library | Version | Purpose |
|---------|---------|---------|
| Eigen 3 | ≥3.3 | Dense/sparse linear algebra, matrix containers |
| SuiteSparse | ≥5.10 | Sparse Cholesky (`CHOLMOD`) for patch solves |
| OpenMP | ≥4.5 | Thread‑level parallelism over coarse elements |
| TBB (optional) | ≥2020 | Alternative task scheduler |

### Install (Ubuntu / WSL)

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
| `LOD2D_BUILD_BENCHMARKS` | `OFF` | Build benchmarks |

## Quick Start

```bash
# Run tests
./build/tests/test_mesh

# Run the main LOD solver (once implemented)
./build/lod2d
```

## Project Structure

```
lod2d-cpp/
├── CMakeLists.txt              Build configuration
├── include/
│   ├── mesh/                   Mesh data types, refinement, edges, area
│   ├── fem/                    DG stiffness assembly
│   ├── lod/                    Quasi‑interpolation, patches, correctors
│   ├── solver/                 Coarse LOD system solver
│   └── utils/                  I/O, timing, helpers
├── src/                        Implementation files
│   ├── mesh/
│   ├── fem/
│   ├── lod/
│   ├── solver/
│   ├── io/
│   └── main.cpp                Entry point
├── tests/                      Unit tests
├── benchmarks/                 Performance benchmarks
├── scripts/                    Helper scripts
├── data/                       Test data & gold references
└── README.md
```

## Module Status

| Module | Status | File |
|--------|--------|------|
| Mesh types | ✅ Done | `include/mesh/types.h` |
| Edge enumeration | ✅ Done | `src/mesh/edges.cpp` |
| Area computation | ✅ Done | `src/mesh/edges.cpp` |
| Red refinement + prolongation | ✅ Done | `src/mesh/refine.cpp` |
| DG assembly | ⬜ Stub | `src/fem/assemble_dg.cpp` |
| Quasi‑interpolation | ⬜ Stub | `src/lod/quasi_interp.cpp` |
| Patch construction | ⬜ Stub | `src/lod/patches.cpp` |
| Corrector solver | ⬜ Stub | `src/lod/corrector.cpp` |
| Coarse LOD solve | ⬜ Stub | `src/solver/coarse_solve.cpp` |

## Test Results (Phase A)

```
=== Mesh Refinement Tests ===
  [PASS] Edge enumeration
  [PASS] Area calculation
  [PASS] Single-level refinement
  [PASS] Multi-level refinement (2,4,6 levels)
  [PASS] Area preserved after 4 refinements

6-level refinement: 2 → 8192 elements in 5.2 ms
```

## Migration Plan

See the [Migration Plan](https://github.com/codenotforce/LOD2d-MATLAB/blob/main/MIGRATION_PLAN.md)
in the MATLAB reference repository for the full roadmap.

## License

Research and educational use.  Based on the LOD code from
*An Introduction to the Localized Orthogonal Decomposition Method*
by A. Målqvist and D. Peterseim.

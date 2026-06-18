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
| Corrector speedup | 2.8× (`parfor` + RCM) | **8.3× (serial), >30× (OpenMP est.)** |
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

| Module | Status | File | Tests |
|--------|--------|------|-------|
| Mesh types | ✅ Phase A | `include/mesh/types.h` | 5/5 |
| Edge enumeration | ✅ Phase A | `src/mesh/edges.cpp` | 5/5 |
| Area computation | ✅ Phase A | `src/mesh/edges.cpp` | 5/5 |
| Red refinement + prolongation | ✅ Phase B | `src/mesh/refine.cpp` | 36/36 golden |
| DG assembly | ✅ Phase C | `src/fem/assemble_dg.cpp` | 10/10 golden |
| Quasi‑interpolation | ✅ Phase D | `src/lod/quasi_interp.cpp` | 3793/3793 golden |
| Patch construction | ✅ Phase E | `src/lod/patches.cpp` | 6/6 golden |
| Corrector solver | ✅ Phase F | `src/lod/corrector.cpp` | 3/3 golden |
| **Full LOD pipeline** | ✅ Phase G | `tests/test_full.cpp` | **3/3 golden** |

## Test Results

### Phase A — Mesh Refinement
```
All tests passed! (5/5)
6-level refinement: 2 → 8192 elements in ~5 ms
```

### Phase B — Golden‑data (vs MATLAB)
```
node/elem/counts · coords · connectivity · Dirichlet · area
P_node · P_elem · P_dg — all exact match across 4 levels
36/36 PASS
```

### Phase C — DG Assembly
```
Coefficient‑aware element stiffness · global sparse assembly
Exact match to MATLAB output (values, nnz, dimensions)
10/10 PASS (3 levels × 3 checks + 1 deterministic test)
```

### Phase D — Quasi‑Interpolation
```
Fine DG mass · inverse coarse DG mass · L2 projection · averaging
All 3793 golden positions match to 1.7e-16
2/2 PASS (dimensions + values)
```

### Phase E — Patch Construction
```
Vertex‑to‑element incidence · adjacency graph · ℓ‑step BFS expansion
6/6 PASS — mesh and nnz exact match with MATLAB across 6 configs (max err 0)
```

### Phase F — Corrector Solver
```
Full corrector: Sph assembly · rhsp · IHp · Sph\RHS · mu
3/3 PASS — CTk matrices exact match to 2.5e-16
```

### Phase G — Full LOD Pipeline (End‑to‑End)
```
H=3, h=5, ℓ=2 — complete LOD solve with reference comparison
uH max diff: 7.5e-15  |  uHms max diff: 7.7e-15  |  uh max diff: 8.9e-16
Energy/L²/FE‑L² errors match MATLAB to machine precision
3/3 PASS  |  C++ 57 ms vs MATLAB 470 ms = 8.3× speedup
```

## Migration Plan

See the [Migration Plan](https://github.com/codenotforce/LOD2d-MATLAB/blob/main/MIGRATION_PLAN.md)
in the MATLAB reference repository for the full roadmap.

## License

Research and educational use.  Based on the LOD code from
*An Introduction to the Localized Orthogonal Decomposition Method*
by A. Målqvist and D. Peterseim.

# LOD2d-C++

`LOD2d-C++` is a C++20 implementation of the two-dimensional Localized
Orthogonal Decomposition (LOD) method. It contains two related workflows:

- elliptic diffusion LOD, migrated from the original MATLAB implementation;
- complex-valued Petrov-Galerkin LOD for the Helmholtz equation.

The project provides conforming longest-edge bisection (LEB), newest-vertex
bisection (NVB), reusable correctors for repeated right-hand sides, benchmark
drivers, and correctness tests against MATLAB-derived or direct-solver
references.

## Quick Start

### Requirements

- CMake 3.16 or newer
- a C++20 compiler
- Eigen 3
- OpenMP, recommended for parallel corrector computation
- SuiteSparse/CHOLMOD, optional experimental sparse backend
- TBB, optional Eigen backend dependency

Ubuntu and WSL:

```bash
sudo apt update
sudo apt install -y build-essential cmake g++ \
  libeigen3-dev libsuitesparse-dev libtbb-dev
```

### Build

```bash
git clone https://github.com/codenotforce/LOD2d-CPP.git
cd LOD2d-CPP
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j "$(nproc)"
```

Available CMake options:

| Option | Default | Purpose |
|---|---:|---|
| `LOD2D_USE_OPENMP` | `ON` | Enable OpenMP parallelism |
| `LOD2D_BUILD_TESTS` | `ON` | Build the CTest suite |
| `LOD2D_BUILD_BENCHMARKS` | `ON` | Build benchmark executables |

### Validate

```bash
ctest --test-dir build --output-on-failure
```

Use a Release build for timings. Debug builds are intended for diagnostics and
correctness checks.

## First Runs

Elliptic LOD correctness and repeated-right-hand-side examples:

```bash
./build/tests/test_corr --solver=eigen
./build/tests/test_full
./build/benchmarks/bench_reuse_rhs --solver=auto --rhs=5
```

Helmholtz manufactured-solution and wave-number examples:

```bash
./build/benchmarks/bench_helmholtz_manufactured \
  --k=4 --h=10 --H-levels=5 --ell-levels=3 \
  --fem-levels=10 --check

./build/benchmarks/bench_helmholtz_k \
  --k=16 --H=auto --h=auto --kH-target=1 \
  --fine-gap=8 --ell=auto --mode=two-sided
```

Benchmark parameter conventions and reproducibility requirements are defined
in [BENCHMARK_GUIDE.md](BENCHMARK_GUIDE.md).

## Main Capabilities

### Elliptic LOD

- real-valued P1 finite elements and DG assembly;
- conforming LEB and NVB refinement with prolongation operators;
- quasi-interpolation of the form `E_H * Pi_H^dg`;
- localized element correctors with Eigen and experimental CHOLMOD backends;
- reusable `LodModel` for fixed coefficient, mesh, `H`, `h`, and `ell` with
  multiple right-hand sides;
- inverse-inequality and performance benchmark drivers.

### Helmholtz LOD

- complex P1 finite elements with homogeneous impedance Robin boundaries;
- globally nested NVB coarse and fine spaces;
- primal and adjoint localized correctors;
- two-sided Petrov-Galerkin LOD;
- direct saddle, explicit Schur, and experimental shifted-GMRES patch solvers;
- experimental two-level hybrid Schwarz solvers;
- stage-1 adaptive coarse-mesh infrastructure and estimator calibration.

The Helmholtz sesquilinear form is linear in its first argument and
conjugate-linear in its second. Artificial patch boundaries are homogeneous
Dirichlet; physical boundary portions retain the Robin term.

## Reusable APIs

For elliptic problems with fixed coefficient and discretization but changing
right-hand sides:

```cpp
#include "lod/lod_model.h"

lod2d::LodProblemConfig config;
config.H = 4;
config.h = 10;
config.ell = 2;

lod2d::LodModel model = lod2d::LodModel::build(config, coefficient);
auto solution = model.solve_from_coarse_values(f_coarse);
```

For Helmholtz problems:

```cpp
#include "helmholtz/model.h"

lod2d::helmholtz::HelmholtzProblemConfig config;
config.H = 5;
config.h = 10;
config.ell = 3;
config.wavenumber = 4.0;

auto model = lod2d::helmholtz::HelmholtzLodModel::build(config);
auto solution = model.solve_source(source);
```

Both models retain correctors and the factorized coarse operator. Changing
only the source does not recompute the correctors. Changing the coefficient,
wave number, mesh, `H`, `h`, `ell`, interpolation, or boundary model requires
a rebuild.

## Documentation

Each document has one role:

| Document | Audience and purpose |
|---|---|
| [HELMHOLTZ_GUIDE.md](HELMHOLTZ_GUIDE.md) | User guide for building, running, and configuring Helmholtz LOD |
| [HELMHOLTZ_ADAPTIVE_GUIDE.md](HELMHOLTZ_ADAPTIVE_GUIDE.md) | User guide for the implemented stage-1 adaptive experiment |
| [BENCHMARK_GUIDE.md](BENCHMARK_GUIDE.md) | Benchmark authoring, timing, memory, validation, and output conventions |
| [DEVELOPMENT.md](DEVELOPMENT.md) | Chronological implementation decisions, measured results, rejected experiments, and remaining engineering work |

Measured timing tables belong only in `DEVELOPMENT.md` and result directories.
User-facing command instructions belong in the corresponding guide.

## Current Status

| Area | Status |
|---|---|
| Elliptic LOD pipeline | Implemented and regression-tested |
| Elliptic repeated-RHS reuse | Implemented |
| Helmholtz Petrov-Galerkin foundation | Implemented |
| Helmholtz patch DirectSaddle/DirectSchur | Implemented; DirectSaddle remains the default |
| Shifted-GMRES and geometric V-cycle | Correct experimental paths; not runtime defaults |
| Fine-space two-level Schwarz | Experimental through stage S4e |
| Adaptive Helmholtz LOD | Stage 1 implemented; later stages paused |
| PML and three-dimensional extensions | Not implemented |

## Repository Layout

```text
include/       Public headers
src/           Library and executable sources
tests/         CTest and golden-data validation
benchmarks/    Performance and numerical experiment drivers
scripts/       Reproducible multi-case runners
results/       Checked-in reference measurements
data/          Input data used by selected experiments
```

## License

Research and educational use. The elliptic implementation is based on the LOD
workflow described in *An Introduction to the Localized Orthogonal
Decomposition Method* by A. Malqvist and D. Peterseim.

# LOD2d-C++

`LOD2d-C++` is a C++20 implementation of two-dimensional Localized
Orthogonal Decomposition (LOD). The repository contains two related workflows:

- elliptic diffusion LOD, migrated from the original MATLAB implementation;
- complex-valued Petrov-Galerkin LOD for the Helmholtz equation.

It provides conforming longest-edge bisection (LEB), newest-vertex bisection
(NVB), reusable correctors for repeated right-hand sides, benchmark drivers,
and regression tests against MATLAB-derived or direct-solver references.

> **Certified-adaptive status.** The Helmholtz paper workflow has a versioned
> experiment contract, paper cases, a validated three-mesh hierarchy, the
> audit-kernel estimator, a fail-closed certificate framework, a checkpointable
> four-stage controller, and a real floating-point numerical backend. That
> backend is intentionally `conditional`: it cannot emit a certified result.
> Verified constants, assembly/corrector and eta evidence, fully directed
> interval propagation, and the frozen six-method paper runner remain open.
> This revision is not a completed paper reproduction.

## Quick Start

### Requirements

- CMake 3.16 or newer
- a C++20 compiler
- Eigen 3
- OpenMP, recommended for parallel corrector computation
- SuiteSparse/CHOLMOD, optional experimental sparse backend
- TBB, optional Eigen backend dependency
- MPFR, MPFI, and GMP, optional, for the directed-rounding spectrum kernel

Ubuntu and WSL:

```bash
sudo apt update
sudo apt install -y build-essential cmake g++ \
  libeigen3-dev libsuitesparse-dev libtbb-dev
```

For the directed-rounding spectrum build, also install:

```bash
sudo apt install -y libmpfr-dev libmpfi-dev libgmp-dev
```

### Build And Validate

```bash
git clone https://github.com/codenotforce/LOD2d-CPP.git
cd LOD2d-CPP
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j "$(nproc)"
ctest --test-dir build --output-on-failure
```

Available CMake options:

| Option | Default | Purpose |
|---|---:|---|
| `LOD2D_USE_OPENMP` | `ON` | Enable OpenMP parallelism |
| `LOD2D_BUILD_TESTS` | `ON` | Build the CTest suite |
| `LOD2D_BUILD_BENCHMARKS` | `ON` | Build benchmark executables |
| `LOD_ENABLE_VERIFIED_CERTIFICATES` | `OFF` | Use MPFR/MPFI directed rounding for certificate verification |

Enable that kernel in a separate build tree (the overall WP4 chain remains
conditional until its other verified inputs and interval operations exist):

```bash
cmake -S . -B build-verified -DCMAKE_BUILD_TYPE=Release \
  -DLOD_ENABLE_VERIFIED_CERTIFICATES=ON
cmake --build build-verified -j "$(nproc)"
ctest --test-dir build-verified --output-on-failure
```

Use a Release build for timings. Debug and sanitizer builds are intended for
diagnostics and correctness checks.

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

./build/benchmarks/bench_helmholtz_H_convergence \
  --k=8 --h=11 --H-levels=4,5,6,7,8 --ell=3 \
  --solver=schur --threads=8 \
  --export-dir=results/helmholtz_H_convergence/snapshots \
  --summary-out=results/helmholtz_H_convergence/summary.csv \
  --check

./build/benchmarks/bench_helmholtz_k \
  --k=16 --H=auto --h=auto --kH-target=1 \
  --fine-gap=8 --ell=auto --source=manufactured \
  --mode=two-sided --check

K_VALUES="4 8 16 32" FINE_GAP=6 \
  bash scripts/run_helmholtz_pollution.sh
```

The following exercises the real WP5 backend as an implementation smoke. It
does not use formal paper parameters and must not be reported as a certified or
WP6 paper run:

```bash
./build/benchmarks/bench_helmholtz_certified \
  --evidence=strict --check
./build/benchmarks/bench_helmholtz_certified \
  --evidence=conditional --check
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

- complex P1 finite elements with impedance Robin boundaries;
- paper-case support for mixed Dirichlet/impedance boundary partitions;
- continuous P1/P2/P3 fine spaces for direct-saddle hp-corrector experiments;
- globally nested NVB spaces with primal and adjoint localized correctors;
- two-sided Petrov-Galerkin LOD;
- direct saddle, explicit Schur, and experimental shifted-GMRES patch solvers;
- experimental two-level hybrid Schwarz solvers;
- a legacy fixed-fine-space H-only calibration driver retained as
  `HLOD-proxy` for diagnostics;
- certified-adaptive WP0-WP5 library infrastructure described below.

The Helmholtz sesquilinear form is linear in its first argument and
conjugate-linear in its second. Artificial patch boundaries are homogeneous
Dirichlet; physical impedance portions retain the Robin term.

### Certified-Adaptive Helmholtz Foundation

| Work package | Implemented scope | Current boundary |
|---|---|---|
| WP0 | Frozen schemas, result statuses, backend/driver parameters, canonical run IDs, and strict parsing | Unified runner consumption is WP6 |
| WP1 | Paper cases, explicit boundary-edge tags, coefficients, sources, and quadrature validation | No production case-matrix runner |
| WP2 | Independently refinable coarse, corrector-fine, and certification-audit meshes with production invariant checks and capacity expansion | Full rebuild is the correctness path |
| WP3 | Audit-kernel Riesz solves, `eta_H`, allocation, real-LOD diagnostics, and immutable evidence fingerprints | The Eigen producer is Diagnostic; no verified `eta_H` producer exists |
| WP4 | Constant registry, spectrum hooks, corrector/stability/error formulas, and context-bound fail-closed evidence | Matrix/scalar directed interval propagation, verified inputs, and independent adjoint fallback remain open |
| WP5 | Four-stage controller, resource limits, checkpoints, no-fallback policy, and a real full-rebuild numerical backend | Backend results are conditional; formal verified runner is not complete |

The legacy `bench_helmholtz_adaptive` executable is not the certified paper
method and must not be reported as CALOD or as the frozen HLOD comparator. See
[HELMHOLTZ_ADAPTIVE_GUIDE.md](HELMHOLTZ_ADAPTIVE_GUIDE.md) for the distinction
and [HELMHOLTZ_ADAPTIVE_LOD_PLAN.md](HELMHOLTZ_ADAPTIVE_LOD_PLAN.md) for the
paper protocol and remaining gates.

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
only the source does not recompute the correctors. Changing a coefficient,
wave number, mesh, `H`, `h`, `ell`, interpolation, or boundary model requires
a rebuild.

## Documentation Map

### User And Contributor Guides

| Document | Purpose |
|---|---|
| [HELMHOLTZ_GUIDE.md](HELMHOLTZ_GUIDE.md) | Stable Helmholtz build, run, solver, and reproducibility guide |
| [HELMHOLTZ_ADAPTIVE_GUIDE.md](HELMHOLTZ_ADAPTIVE_GUIDE.md) | Legacy proxy and certified-adaptive WP0-WP5 guide |
| [BENCHMARK_GUIDE.md](BENCHMARK_GUIDE.md) | Benchmark implementation, validation, timing, memory, and output rules |
| [DEVELOPMENT.md](DEVELOPMENT.md) | Chronological decisions, measured results, defects, and rejected experiments |

### Paper Protocol And Operations

| Document | Purpose |
|---|---|
| [HELMHOLTZ_ADAPTIVE_LOD_PLAN.md](HELMHOLTZ_ADAPTIVE_LOD_PLAN.md) | Active paper-reproduction work packages, scientific gates, and result contract |
| [experiments/helmholtz_adaptive_paper/README.md](experiments/helmholtz_adaptive_paper/README.md) | Versioned experiment schemas, IDs, timing ownership, and current execution boundary |
| [HELMHOLTZ_HP_SERVER_RUNBOOK.md](HELMHOLTZ_HP_SERVER_RUNBOOK.md) | hp-LOD validation and server runs |
| [HELMHOLTZ_POLLUTION_SERVER_RUNBOOK.md](HELMHOLTZ_POLLUTION_SERVER_RUNBOOK.md) | Staged pollution scans, memory gates, and recovery |

Other `*_PLAN.md` files are retained as scoped design records. Their dated
status statements are historical unless the active adaptive plan explicitly
adopts them.

## Current Status

| Area | Status |
|---|---|
| Elliptic LOD pipeline | Implemented and regression-tested |
| Elliptic repeated-RHS reuse | Implemented |
| Helmholtz Petrov-Galerkin foundation | Implemented |
| Helmholtz continuous hp fine-space experiment | P1/P2/P3 implemented; deep server matrix pending |
| Helmholtz patch DirectSaddle/DirectSchur | Implemented; DirectSaddle remains the public default |
| Shifted-GMRES and geometric V-cycle | Correct experimental paths; not runtime defaults |
| Fine-space two-level Schwarz | Experimental through stage S4e |
| Legacy adaptive H-only proxy | Implemented for regression and calibration only |
| Certified adaptive Helmholtz LOD | Conditional WP0-WP5 implementation smoke available; verified certificate chain and WP6 paper runner pending |
| PML and three-dimensional extensions | Not implemented |

## Repository Layout

```text
include/       Public headers
src/           Library and executable sources
tests/         CTest and golden-data validation
benchmarks/    Performance and numerical experiment drivers
experiments/   Versioned experiment contracts and frozen inputs
scripts/       Reproducible multi-case runners
results/       Checked-in reference measurements
data/          Input data used by selected experiments
tools/         Validation and data-processing helpers
```

## License

Research and educational use. The elliptic implementation is based on the LOD
workflow described in *An Introduction to the Localized Orthogonal
Decomposition Method* by A. Malqvist and D. Peterseim.

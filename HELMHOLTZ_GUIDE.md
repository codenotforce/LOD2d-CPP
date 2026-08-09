# Helmholtz Petrov-Galerkin LOD User Guide

## Purpose

This guide explains how to run the implemented Helmholtz workflows and choose
their solver options. It does not contain performance history or future-task
definitions:

- measured results and engineering decisions: [DEVELOPMENT.md](DEVELOPMENT.md);
- legacy adaptive proxy and certified-adaptive WP0-WP5 boundary:
  [HELMHOLTZ_ADAPTIVE_GUIDE.md](HELMHOLTZ_ADAPTIVE_GUIDE.md).

## Mathematical Convention

The implemented model is

```math
-\nabla\cdot(A\nabla u)-k^2 n u=f \quad\text{in }\Omega,
\qquad
A\nabla u\cdot\nu-i k\beta u=0 \quad\text{on }\partial\Omega.
```

The finite-element form is linear in its first argument and conjugate-linear
in its second. The default unit-square problem uses `A=n=beta=1`.

Coarse and fine spaces come from one globally nested NVB hierarchy. On a
corrector patch, physical-domain boundary edges retain the Robin term and
artificial patch edges use homogeneous Dirichlet data. The complex Helmholtz
and saddle systems must not be solved with LLT, CG, or PCG.

## Build And Test

Install dependencies and perform the initial build as described in
[README.md](README.md). To rebuild only Helmholtz targets:

```bash
cmake --build build --target \
  test_helmholtz_fem test_helmholtz_interp \
  test_helmholtz_corrector test_helmholtz_model -j 8

ctest --test-dir build -R helmholtz_ --output-on-failure
```

## Wave-Number Experiment

Run one manufactured or Gaussian-source comparison with:

```bash
OMP_NUM_THREADS=8 ./build/benchmarks/bench_helmholtz_k \
  --k=16 --H=auto --h=auto --kH-target=1 \
  --fine-gap=8 --ell=auto --mode=two-sided
```

| Option | Meaning |
|---|---|
| `--k` | Wave number |
| `--H=auto` | First global NVB level satisfying the requested `kH` target |
| `--h=auto` | Set the fine level to `H + fine-gap` |
| `--fine-gap` | Additional global NVB levels used by correctors and reference FEM |
| `--ell=auto` | Use the benchmark's logarithmic oversampling policy |
| `--mode=two-sided` | Trial `(I-C)V_H`, test `(I-C*)V_H` |
| `--mode=test-only` | Trial `V_H`, test `(I-C*)V_H` |
| `--format=csv` | Emit one machine-readable row |

`fine-gap=8` is the accuracy-oriented default. A smaller gap is suitable for
smoke tests but does not establish fine-grid convergence.

### Reading The Output

The discrete energy norm is

```math
\|v\|_{1,k}^2=v^*(K+k^2M)v.
```

A wave-number robustness claim requires all of the following:

1. `kH` remains approximately fixed across cases.
2. `kh` is sufficiently small and increasing `fine-gap` does not materially
   change the LOD error.
3. Corrector, constraint, and Petrov residuals remain near solver tolerance.
4. LOD and standard FEM use the same source and fine reference.
5. Error, time, memory, and the actual mesh levels are reported together.

Do not infer a pollution rate from one source family. A Gaussian source changes
the solution with `k`; use the manufactured solution when an exact-error study
is required.

## Corrector Patch Solvers

`HelmholtzProblemConfig::patch_solver.kind` selects:

| Solver | Role |
|---|---|
| `DirectSaddle` | Full constrained saddle solve; default and gold standard |
| `DirectSchur` | Explicit Schur elimination; fastest tested large-patch experiment |
| `ShiftedGmres` | Right-preconditioned GMRES; correctness/research path |

For

```math
A_\omega=K_\omega-k^2M_\omega-i kR_\omega
```

the Schur paths compute

```math
Y=A_\omega^{-1}F,\qquad Z=A_\omega^{-1}B^*,
\qquad (BZ)\Lambda=BY,\qquad X=Y-Z\Lambda.
```

The shifted preconditioner is

```math
P_{\omega,\varepsilon}
=K_\omega-(k^2+i\varepsilon)M_\omega-i kR_\omega,
\qquad \varepsilon=\alpha k^2.
```

`ShiftedGmres` checks the original residual, not only a preconditioned
residual. Its inverse can be exact shifted SparseLU, identity for an
unpreconditioned baseline, or a fixed geometric V-cycle.

Compare all patch solvers on the same data with:

```bash
cmake --build build --target bench_helmholtz_patch_solver -j 8
OMP_NUM_THREADS=8 ./build/benchmarks/bench_helmholtz_patch_solver \
  --H=5 --h=10 --ell=3 --k=4 --threads=8 \
  --solver=all --inverse=vcycle --alpha=0.2 --tol=1e-10 \
  --pre-smooth=2 --post-smooth=2 --coarse-max=200 --omega=0.6
```

Use `fallback_to_direct=true` only in an explicitly labelled fallback
experiment. A Schur path may fail when unrestricted `A_omega` is locally
singular even though the constrained saddle problem is well posed.

Current policy: keep `DirectSaddle` as the public default. `DirectSchur` is the
preferred performance experiment for large patches. Shifted-LU, unpreconditioned
GMRES, and the current V-cycle are not runtime defaults.

## Experimental Two-Level Schwarz

`HelmholtzTwoLevelSchwarzPreconditioner` combines the factored LOD coarse
solve with overlapping fine-space local solves. Run the common comparison:

```bash
cmake --build build --target bench_helmholtz_two_level_schwarz -j 8
./build/benchmarks/bench_helmholtz_two_level_schwarz \
  --H=7 --h=12 --ell=3 --k=8 --threads=8 \
  --solver=all --boundary=dirichlet --extension=weighted \
  --local-solver=direct
```

Important switches:

| Option | Values and role |
|---|---|
| `--solver` | `identity`, `local`, `additive`, `hybrid`, or `all` |
| `--boundary` | Artificial `dirichlet` or experimental `impedance` |
| `--extension` | Inverse-multiplicity `weighted` or Boolean `restricted` |
| `--local-solver` | Exact `direct` or experimental `shifted-gmres` |
| `--local-inverse` | Shifted `lu` or fixed `vcycle` |
| `--factorization-reuse` | `none` or strict `identical` local-matrix reuse |

`factorization-reuse=identical` is valid only with the direct local solver.
It hashes and then exactly compares compressed local matrices before sharing
one factorization and using a multi-right-hand-side solve. It remains explicit
because grouping overhead can dominate small cases.

A tolerance-based inner GMRES makes the preconditioner variable, so the outer
method uses FGMRES. The direct local path is fixed and uses ordinary
right-preconditioned GMRES.

Current policy: Dirichlet, weighted overlap, direct local SparseLU is the
strongest tested action. ORAS, RestrictedCore, shifted local GMRES, and the
V-cycle remain research switches.

## Repeated Right-Hand Sides

For fixed mesh, `k`, coefficient, Robin data, interpolation, and `ell`, build
the model once:

```cpp
HelmholtzLodModel model = HelmholtzLodModel::build(config);
auto first = model.solve_source(source_one);
auto second = model.solve_source(source_two);
```

The correctors, Petrov-Galerkin bases, coarse operator, and coarse
factorization are reused. Rebuild after changing any parameter that changes
the operator or multiscale space.

## Multi-Case Server Runs

Run each wave number in a separate process so all patch data is released:

```bash
THREADS=32 JOBS=16 K_VALUES="4 8 16 32 64" \
FINE_GAP=8 KH_TARGET=1 MODE=two-sided \
bash scripts/run_helmholtz_k_scan.sh
```

The script writes `scan.csv`, per-case logs, `/usr/bin/time` files, and
metadata under `results/helmholtz_k_scan/`. Use `RESUME=1` to preserve
completed rows and `CONTINUE_ON_ERROR=1` only when later cases should continue
after a failure.

For two-level Schwarz scale studies:

```bash
THREADS=8 JOBS=8 CASES="5:10:4 7:12:8 9:14:16" \
ELL_VALUES=3 FACTORIZATION_REUSE=identical \
bash scripts/run_helmholtz_schwarz_scale.sh
```

## Reproducibility Rules

- Use a clean Release build and record the Git revision.
- Record actual `H`, `h`, `ell`, `kH`, `kh`, solver configuration, thread
  count, and OpenMP placement.
- Preserve CSV, metadata, command logs, and peak-RSS reports.
- Validate every experimental solver against DirectSaddle or fine SparseLU.
- Treat true residual and constraint residual as correctness gates.
- Put new measured results in `DEVELOPMENT.md`, not in this guide.

## Adaptive Workflow

Use [HELMHOLTZ_ADAPTIVE_GUIDE.md](HELMHOLTZ_ADAPTIVE_GUIDE.md) for both the
legacy fixed-fine-space H-only calibration driver and the separate certified
WP0-WP5 library foundation. The legacy executable is `HLOD-proxy`; it is not a
formal CALOD or frozen-HLOD paper runner. Uniform estimator reliability and a
completed paper experiment matrix are not currently claimed.

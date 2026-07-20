# Helmholtz Adaptive LOD Stage-1 User Guide

## Purpose And Status

Stage 1 is implemented as a correctness-first calibration workflow:

- the fine finite-element space is fixed;
- only the coarse NVB mesh is adaptively refined;
- every iteration fully rebuilds the LOD model;
- strong-residual candidates are compared with exact or fine-reference errors.

Incremental corrector reuse, L-shaped domains, joint `H/h/ell` adaptation, and
PML are not part of stage 1. This guide documents only the implemented
calibration workflow. Measured calibration histories are recorded in
[DEVELOPMENT.md](DEVELOPMENT.md).

## Implemented Invariants

- One NVB tree supplies all coarse meshes and the fixed master fine mesh.
- Coarse elements retain stable IDs and explicit parents.
- Every adaptive coarse mesh is completed to the same master fine level.
- Element and nodal prolongation identities are checked.
- NVB closure may not create duplicate global midpoint nodes.
- The reference fine solution and finite-element space remain unchanged
  throughout one adaptive run.

## Build And Test

```bash
cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug
cmake --build build-debug --target \
  test_helmholtz_adaptive test_helmholtz_reliability -j 8
ctest --test-dir build-debug -R 'helmholtz_(adaptive|reliability)' \
  --output-on-failure
```

## Run A Calibration

```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release --target bench_helmholtz_adaptive -j 8

./build-release/benchmarks/bench_helmholtz_adaptive \
  --k=4 --H=5 --h=10 --ell=3 --iterations=3 \
  --theta=0.5 --estimator=fine --q-limit=0.5 --threads=8 \
  --source=manufactured
```

Key options:

| Option | Meaning |
|---|---|
| `--H` | Initial coarse NVB level |
| `--h` | Fixed master fine level |
| `--ell` | Corrector oversampling layers |
| `--iterations` | Maximum adaptive iterations |
| `--theta` | Dorfler marking fraction |
| `--estimator` | `fine`, `mixed`, or `macro` aggregation |
| `--q-limit` | Coarse/fine scale-separation gate |
| `--source` | `manufactured` or `gaussian` |
| `--no-dual` | Skip expensive local energy-Riesz calibration |
| `--format=csv` | Emit a machine-readable history |
| `--mesh-out=PATH` | Write final mesh IDs, geometry, levels, and indicators |

The driver stops before a refinement whose NVB closure would violate the
configured scale-separation limit. Dense inf-sup diagnostics are omitted for
large coarse systems.

## Manufactured Solution

The exact validation problem uses

```math
\phi(t)=16t^2(1-t)^2,
\qquad
u(x,y)=\phi(x)\phi(y)e^{ikx},
```

and

```math
f=-e^{ikx}\left[
\phi''(x)\phi(y)+\phi(x)\phi''(y)
+2ik\phi'(x)\phi(y)\right].
```

Because `phi=phi'=0` at the endpoints, the homogeneous impedance condition is
satisfied. This source reports both the exact LOD error and the fixed fine-FEM
discretization floor.

## Residual Identity

The estimator assembly reconstructs the integrated nodal residual from broken
volume terms, interior flux jumps, and Robin boundary terms. It must agree
with

```math
r_h=b_h-A_hu_{\mathrm{LOD}}.
```

This identity is a correctness gate for signs, edge ownership, and complex
conjugation. It is not by itself a reliability theorem.

## Interpreting The Indicators

The strong-residual candidates target the continuous error
`u-u_LOD`. They are not estimators of the purely discrete difference
`u_h-u_LOD`.

When the coarse and fine meshes coincide, `u_h-u_LOD` and the algebraic
residual vanish, while `u-u_h` and squared strong volume/jump residuals need
not vanish. Therefore:

- do not divide a strong-residual indicator by a vanishing discrete error and
  interpret the resulting effectivity as a failure;
- use the manufactured exact solution for continuous-error calibration;
- use the local algebraic energy-dual residual as a discrete calibration
  quantity;
- require parameter sweeps and holdout validation before claiming uniform
  reliability or efficiency.

A future estimator specifically targeting `u_h-u_LOD` should use an algebraic
dual-residual localization or a hierarchical two-level construction.

## Output Checklist

For every adaptive iteration retain:

- coarse element and DOF counts;
- exact and/or fine-reference errors;
- all requested indicators and effectivities;
- marked elements and NVB closure growth;
- `q_max` and effective scale separation;
- residual identity and solver residuals;
- inf-sup value when available;
- solve, estimate, mark, refine, and total timings;
- final mesh output and run metadata.

## Scientific Gate

Stage 1 validates implementation and estimator behavior on the unit square.
It does not yet establish a theorem or an adaptive optimality result. Continue
to stage 2 only after the candidate estimator has a stable reliability
envelope and acceptable marking correlation on training and holdout cases.

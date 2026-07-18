# Helmholtz Petrov-Galerkin LOD Guide

## Scope

The current implementation solves the two-dimensional complex Helmholtz problem

```math
-\nabla\cdot(A\nabla u)-k^2 n u=f\quad\text{in }\Omega,
\qquad
A\nabla u\cdot\nu-i k\beta u=0\quad\text{on }\partial\Omega,
```

with P1 finite elements, globally nested NVB meshes, and localized
Petrov-Galerkin LOD. The default experiment uses `A=n=beta=1` on the unit
square. Stage-1 adaptive coarse refinement is available as a calibration
workflow; PML is not implemented.

## Build And Validate

Install the dependencies on Ubuntu:

```bash
sudo apt update
sudo apt install -y build-essential cmake g++ libeigen3-dev libsuitesparse-dev libtbb-dev
```

Configure a Release build and run every registered test:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 8
ctest --test-dir build --output-on-failure
```

The Helmholtz-specific tests check FEM convergence, complex interpolation,
physical/artificial patch boundaries, primal/adjoint correctors, and both
Petrov-Galerkin modes.

## Run One Wave Number

```bash
OMP_NUM_THREADS=8 ./build/benchmarks/bench_helmholtz_k \
  --k=16 \
  --H=auto \
  --h=auto \
  --kH-target=1 \
  --fine-gap=8 \
  --ell=auto \
  --mode=two-sided
```

Important options:

| Option | Meaning |
|---|---|
| `--k` | Wave number |
| `--H=auto` | First global NVB level satisfying `kH <= kH-target` |
| `--h=auto` | Use `h-level = H-level + fine-gap` |
| `--fine-gap` | Number of additional global NVB levels for correctors/reference |
| `--ell=auto` | Use `ceil(log2(k))` patch layers |
| `--mode=two-sided` | Trial `(I-C)V_H`, test `(I-C*)V_H` |
| `--mode=test-only` | Trial `V_H`, test `(I-C*)V_H` |
| `--stability-max-dofs` | Maximum coarse size for dense energy-scaled inf-sup SVD |
| `--format=csv` | Print one machine-readable result row |

`fine-gap=8` is the default accuracy-oriented scan. Smaller values such as 4
are useful for smoke/performance tests but should not be treated as the final
fine-grid convergence study.

## Server Scan

Always invoke the helper through `bash`; this also works when a Windows Git
checkout does not preserve the executable bit.

```bash
THREADS=32 JOBS=16 \
K_VALUES="4 8 16 32 64" \
FINE_GAP=8 KH_TARGET=1 MODE=two-sided \
bash scripts/run_helmholtz_k_scan.sh
```

The script runs every wave number in a separate process, so the previous mesh,
correctors, and sparse factors are released before the next case. Results are
written to `results/helmholtz_k_scan/`:

- `scan.csv`: one row per successful wave number;
- `k*.log`: command, CSV row, diagnostics, and errors;
- `k*.time`: `/usr/bin/time -v` resource report;
- `metadata.txt`: Git revision/status, compiler, kernel, OpenMP placement, and
  all scan parameters.

Resume a partially completed scan without repeating existing CSV rows:

```bash
RESUME=1 THREADS=32 K_VALUES="4 8 16 32 64" \
  bash scripts/run_helmholtz_k_scan.sh
```

Set `CONTINUE_ON_ERROR=1` when one failed wave number should not stop later
cases.

## Reading The Output

The main error columns compare standard coarse FEM and LOD with the same fine
P1 reference solution in

```math
\|v\|_{1,k}^2=v^*(K+k^2M)v.
```

For a wave-number robustness study, verify all of the following:

1. `kH` is approximately constant across rows.
2. `kh` is sufficiently small and a second run with a larger `fine-gap` does
   not materially change the reported LOD error.
3. Corrector, constraint, and Petrov residuals remain near solver tolerance.
4. The LOD relative error or `LOD energy/(H||f||)` remains bounded as `k`
   increases.
5. Standard FEM and LOD use exactly the same source and fine reference.

The fixed Gaussian-source scan is a reproducibility baseline. Because the
solution itself changes with `k`, a nonmonotone standard FEM error is possible;
do not infer a pollution slope from one source family alone.

## Corrector Performance

The corrector is parallelized over coarse elements. Each thread reuses stamped
node/constraint workspaces and the most recent matching `SparseLU` symbolic
pattern. Three local coarse basis right-hand sides are solved as one dense
block, and the adjoint corrector is obtained by conjugating the primal
corrector under the current real-coefficient assumptions.

Control parallelism with `OMP_NUM_THREADS` or the scan script's `THREADS`
variable. Report the thread count whenever publishing timing results.

## Reusing Different Right-Hand Sides

For fixed mesh, `k`, coefficients, boundary data, interpolation, and `ell`, the
correctors and coarse factorization do not depend on `f`:

```cpp
HelmholtzLodModel model = HelmholtzLodModel::build(config);
auto u1 = model.solve_source(source_one);
auto u2 = model.solve_source(source_two);
```

Changing `k`, the mesh, coefficients, Robin data, interpolation, or `ell`
requires rebuilding the model.

## Reproducibility Checklist

- Use a clean Release build and record the Git revision.
- Keep the NVB-compatible initial reference-edge labeling.
- Record `k`, actual `H/h`, `ell`, mode, thread count, and OpenMP placement.
- Preserve `scan.csv`, `metadata.txt`, and every `.time` file.
- Confirm the fine-grid result with a larger `fine-gap` at the largest `k`.
- Do not use LLT, CG, or PCG for the complex Helmholtz or saddle systems.
- Keep performance conclusions in `DEVELOPMENT.md`; README contains commands
  and stable user-facing behavior only.

## Adaptive Coarse-Mesh Calibration

The fixed-fine-space adaptive Helmholtz implementation is documented in
`HELMHOLTZ_ADAPTIVE_GUIDE.md`. It includes the mathematical status of the
strong-residual candidates, fixed-master-mesh invariants, benchmark fields,
validation commands, a manufactured exact solution, and pilot tables. The
strong-residual candidates are calibrated against `u-u_LOD`, not
`u_h-u_LOD`. Treat the current adaptive benchmark as a calibration experiment;
uniform reliability and efficiency are not yet established.

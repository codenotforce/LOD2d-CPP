# Helmholtz Adaptive LOD Stage 1

## Status

The stage-1 engineering baseline is implemented. The strong-residual candidates
must be interpreted against the continuous error `u-u_LOD`, not against the
purely discrete difference `u_h-u_LOD`. A manufactured exact solution is now
available for that distinction. A full reliability and efficiency theorem is
still outside the current numerical validation.

Implemented components:

- one NVB hierarchy with stable coarse-element IDs and explicit parent IDs;
- completion of every adaptive coarse mesh to one fixed master fine level;
- exact checks of element and nodal prolongation identities;
- an adaptive Helmholtz model entry point sharing the uniform model pipeline;
- one-pass broken element, interior-flux-jump, and Robin residual assembly;
- fine, mixed, and macro coarse-element aggregations;
- local energy-Riesz dual indicators for calibration;
- deterministic minimal Dorfler marking;
- reliability-envelope evaluation, constrained fitting, and holdout utilities;
- a full-rebuild SOLVE-ESTIMATE-MARK-REFINE driver;
- CSV/human benchmark output, mesh snapshots, timings, residual diagnostics,
  scale-separation diagnostics, and a dense energy inf-sup diagnostic.

Incremental corrector reuse is intentionally deferred. The full-rebuild path is
the numerical gold standard that an incremental implementation must match.

## Build And Test

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j 8
ctest --test-dir build --output-on-failure
```

The stage-1 tests are:

```bash
./build/tests/test_helmholtz_adaptive
./build/tests/test_helmholtz_reliability
```

They check fixed-fine-mesh geometry, adaptive NVB closure, prolongation,
complex residual reconstruction, all three aggregations, local Riesz solves,
Dorfler minimality, decay fitting, reliability-domain rejection, constrained
training envelopes, Spearman correlation, and marking overlap. The NVB tests
also reject coincident nodes with different global indices, which would change
the fixed fine finite element space during coarse-mesh adaptation.

## Run The Calibration Benchmark

```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release --target bench_helmholtz_adaptive -j 8

./build-release/benchmarks/bench_helmholtz_adaptive \
  --k=4 --H=5 --h=10 --ell=3 --iterations=3 \
  --theta=0.5 --estimator=fine --q-limit=0.5 --threads=8 \
  --source=manufactured
```

Available production-candidate names are `fine`, `mixed`, and `macro`.
`--source=manufactured` uses `u=phi(x)phi(y)exp(i*k*x)`, where
`phi(t)=16*t^2*(1-t)^2`, and reports both `u-u_LOD` and `u-u_h` errors.
Use `--no-dual` for larger timing runs. Use `--format=csv` for tables and
`--mesh-out=mesh.csv` for the final element IDs, levels, coordinates, and
indicators.

The driver stops before a marked refinement whose NVB closure would violate
the configured coarse/fine scale-separation limit. The dense inf-sup SVD is
reported only up to 512 coarse unknowns; larger systems report `nan`.

## Residual Identity

The implementation reconstructs the nodal residual from the integrated broken
volume, flux-jump, and Robin terms and compares it with

```text
r_h = b_h - A_h u_LOD.
```

The Debug test gives a relative difference of approximately `2.53e-17`,
well below the `1e-10` gate. This also caught and fixed a complex-flux bug:
Eigen's complex `dot` conjugates its left operand, so the real normal must be
placed on the left.

## Error Target And Scientific Gate

For `T_H = T_h`, the discrete difference `u_h-u_LOD` and algebraic residual
vanish, but the continuous error `u-u_h` generally does not. Likewise, the sum
of squared strong volume and jump residuals generally does not vanish. Galerkin
orthogonality is obtained only after the signed terms are tested and cancel;
squaring the terms separately removes that cancellation.

Consequences:

- `eta_fine`, `eta_mix`, and `eta_macro` are not estimators of
  `u_h-u_LOD`;
- their discrete effectivity `eta/||u_h-u_LOD||` can grow without bound;
- this does not invalidate them as candidates for the continuous error
  `u-u_LOD`, whose limiting value is the nonzero fine discretization error;
- the plan's original requirement that a strong residual indicator reach
  machine precision at `T_H=T_h` is mathematically incompatible with these
  definitions;
- the local algebraic energy-dual residual is still a valid calibration
  quantity because it acts on `V_h` and does vanish with `r_h`.

For a discrete-error estimator, the next candidate remains an algebraic
dual-residual localization or a hierarchical two-level estimator. For the
continuous-error target, the manufactured solution provides the appropriate
calibration; broader parameter sweeps are still required before claiming
uniform reliability or efficiency.

## Manufactured-Solution Validation

The exact pair is

```text
phi(t) = 16 t^2 (1-t)^2,
u(x,y) = phi(x) phi(y) exp(i k x),
f = -exp(i k x) [phi''(x)phi(y) + phi(x)phi''(y)
                  + 2 i k phi'(x)phi(y)].
```

Because `phi=phi'=0` at both endpoints, this solution satisfies the
homogeneous Robin condition. Direct differentiation gives
`-Delta u-k^2 u=f`.

Release run: `k=4,H=5,h=10,ell=3`, `theta=0.5`, seven iterations.

| iteration | coarse elements | exact LOD energy error | exact fine energy error | fine eta | exact effectivity |
|---:|---:|---:|---:|---:|---:|
| 0 | 64 | 0.529410 | 0.208842 | 1.26604 | 2.39142 |
| 1 | 84 | 0.417764 | 0.208842 | 1.25424 | 3.00226 |
| 2 | 116 | 0.363453 | 0.208842 | 1.24709 | 3.43122 |
| 3 | 176 | 0.280067 | 0.208842 | 1.24106 | 4.43129 |
| 4 | 282 | 0.240944 | 0.208842 | 1.23662 | 5.13238 |
| 5 | 404 | 0.228634 | 0.208842 | 1.23151 | 5.38640 |
| 6 | 648 | 0.214176 | 0.208842 | 1.22683 | 5.72816 |

The invariant fine error confirms that every adaptive coarse mesh uses the same
master fine finite element space. The exact LOD error decreases monotonically
toward that fine-discretization floor. In contrast, the discrete effectivity
grows from `2.61` to `25.91`, illustrating why the reference and exact
solutions answer different estimator questions.

## Gaussian Pilot Result

Release pilot: `k=4`, initial level 5, fine level 10, `ell=3`,
`theta=0.5`, Gaussian source, four threads.

| iteration | coarse elements | relative energy error | fine eta | effectivity | inf-sup |
|---:|---:|---:|---:|---:|---:|
| 0 | 64 | 0.122562 | 0.0103049 | 2.27068 | 0.318864 |
| 1 | 70 | 0.155091 | 0.0175353 | 7.04001 | 0.558892 |
| 2 | 88 | 0.108361 | 0.0170622 | 12.1898 | 0.664904 |

Increasing `ell` from 3 to 6 produced essentially the same trajectory.
The inf-sup value improves rather than collapses, so the temporary error
increase is not caused by Petrov-Galerkin instability. It is consistent with
non-nested multiscale spaces after local coarse refinement, while the rapidly
growing effectivity exposes the strong-residual mismatch described above.

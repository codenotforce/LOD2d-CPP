# Helmholtz Adaptive LOD Guide

## Scope And Status

The repository currently contains two different adaptive Helmholtz paths. They
serve different purposes and must not be reported under the same method name.

| Path | Entry point | Intended use | Paper status |
|---|---|---|---|
| Legacy H-only calibration | `bench_helmholtz_adaptive` | Regression, diagnostics, and strong-residual calibration on a fixed fine space | `HLOD-proxy`; excluded from the six-method comparator matrix |
| Conditional numerical smoke | `bench_helmholtz_certified` | Real hierarchy/LOD/WP3/WP4/audit wiring and strict fail-closed validation | Implementation smoke only; never a certified or WP6 paper run |
| Certified-adaptive foundation | `helmholtz/adaptive/*.h` plus CTest targets | WP0-WP5 contracts, estimators, certificates, state machine, and numerical backend | Formal verified CALOD/HLOD remains unavailable until G3 and WP6 pass |

The active scientific protocol and remaining work are defined in
[HELMHOLTZ_ADAPTIVE_LOD_PLAN.md](HELMHOLTZ_ADAPTIVE_LOD_PLAN.md). Versioned
experiment contracts are documented in
[experiments/helmholtz_adaptive_paper/README.md](experiments/helmholtz_adaptive_paper/README.md).

## Implemented Work Packages

### WP0: Experiment Contract

- strict v1 input and output schemas;
- canonical JSON configuration hashing and immutable run IDs;
- hashing of every driver/resource and numerical-backend decision parameter;
- five explicit certificate-status labels and a closed metric registry;
- explicit provenance, timing ownership, nullable numeric values, censored
  resource states, and method labels;
- rejection of unknown fields and non-RFC-8259 numeric values.

### WP1: Paper Problems

- registered paper cases and manufactured/source definitions;
- mixed Dirichlet/impedance boundary tagging;
- coefficient, geometry, and quadrature checks;
- exact-solution checks remain validation data and are not exposed to adaptive
  MARK/STOP decisions.

### WP2: Three-Mesh Hierarchy And Error Roles

- independently refinable coarse mesh `T_H`, corrector-fine mesh `T_h`, and
  certification-audit mesh `T_aud`;
- exact P1, element, and DG embeddings for conforming nested meshes;
- parent maps, version counters, local fine refinement, and kernel-constraint
  restriction;
- automatic fine/audit capacity expansion when coarse refinement catches the
  old fine mesh, followed by production composition/right-inverse checks;
- separate certification-audit and post-processing evaluation-reference
  services, with their timings owned by different accounting domains.

The full-rebuild path is the correctness reference after an `H`, `h`, or audit
mutation. Incremental reuse must reproduce it before it can become a default.

### WP3: Coarse Audit-Kernel Estimator

- residual Riesz solves in the constrained audit kernel;
- elementwise source allocation and the global `eta_H` quantity;
- primal/adjoint, constraint, stationarity, and energy-identity diagnostics;
- allocation-sum and conjugation checks;
- read-only provenance tokens binding the hierarchy, operators, inputs, eta,
  and element allocation.

`eta_H` is not automatically verified just because the linear solve succeeds.
The current Eigen producer always emits Diagnostic evidence; no verified
`eta_H` producer exists yet.

### WP4: Corrector And Stability Certificates

- directional certificate-constant registry and derived overlap constants;
- diagnostic matrix-radius propagation and verified generalized-spectrum
  kernel interfaces;
- total, fine-discretization, localization, stability, and LOD-error bounds;
- elementwise `eta_h` allocation and a fail-closed conjugation gate;
- MPFR/MPFI directed-rounding backend behind a CMake option.

The complete WP4 certificate is not yet verified. It lacks theorem constants,
verified FE/corrector assembly and `eta_H` producers, fully directed matrix and
scalar interval propagation, and an independent adjoint-corrector fallback.
The implementation therefore forces the result to `conditional`, even when
the standalone MPFR/MPFI spectrum kernel is enabled.

### WP5: Certified Decision State Machine

The checkpointable driver enforces this order:

```text
coarse admissibility
  -> corrector certification
  -> coarse error control
  -> audit control
  -> done / work limit / failure
```

It supports independent coarse, corrector-fine, oversampling, and audit
mutations; Dorfler marking; pending coarse marks while the audit space is
refined; structured terminal codes; resource limits; and deterministic
checkpoint/resume. CALOD never silently falls back to H-only adaptation.
The HLOD path checks its frozen prior corrector-space identifier and
oversampling value before every decision.

`NumericalCertifiedBackend` connects the controller to real meshes, LOD
correctors and solves, WP3, WP4, and an empirical two-level audit diagnostic.
All its observations remain Conditional; `RequireVerified` stops with
`UnverifiedEvidence`.

## Build And Test

### Standard Eigen Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j "$(nproc)"
ctest --test-dir build --output-on-failure
```

The WP0-WP5 regression targets are:

```bash
cmake --build build --target \
  test_helmholtz_adaptive \
  test_helmholtz_error_control \
  test_helmholtz_kernel_residual \
  test_verified_spectrum \
  test_helmholtz_certificates \
  test_helmholtz_certified_driver \
  test_helmholtz_numerical_backend \
  test_helmholtz_paper_config \
  test_helmholtz_mixed_boundary \
  test_helmholtz_paper_cases \
  test_helmholtz_quadrature -j "$(nproc)"

ctest --test-dir build \
  -R 'helmholtz_(adaptive|error_control|kernel_residual|certificates|certified_driver|numerical_backend|paper_config|mixed_boundary|paper_cases|quadrature)|verified_spectrum' \
  --output-on-failure
```

### Verified-Certificate Build

Install `libmpfr-dev`, `libmpfi-dev`, and `libgmp-dev`, then use a separate
build tree:

```bash
cmake -S . -B build-verified -DCMAKE_BUILD_TYPE=Release \
  -DLOD_ENABLE_VERIFIED_CERTIFICATES=ON
cmake --build build-verified -j "$(nproc)"
ctest --test-dir build-verified --output-on-failure
```

The option enables the standalone directed-rounding spectrum code. It does not
manufacture missing evidence or fix the remaining ordinary-double matrix and
scalar propagation, and therefore does not promote the WP4 chain to verified.

## Run The Conditional Implementation Smoke

```bash
cmake --build build --target bench_helmholtz_certified -j "$(nproc)"

./build/benchmarks/bench_helmholtz_certified \
  --evidence=strict --check
./build/benchmarks/bench_helmholtz_certified \
  --evidence=conditional --check
```

The strict command must stop with `UnverifiedEvidence`. The conditional command
must traverse the real coarse/corrector/error/audit path. The executable uses a
small R1 configuration and prints `runner_scope=implementation-smoke` and
`wp6_runner=false`; it neither consumes the frozen paper matrix nor emits the
common paper output schema.

## Run The Legacy HLOD Proxy

Build and run the fixed-fine-space calibration executable with:

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
| `--iterations` | Maximum H-only adaptive iterations |
| `--theta` | Dorfler marking fraction |
| `--estimator` | `fine`, `mixed`, or `macro` strong-residual aggregation |
| `--q-limit` | Legacy coarse/fine scale-separation gate |
| `--source` | `manufactured` or `gaussian` |
| `--no-dual` | Skip the local energy-Riesz calibration diagnostic |
| `--format=csv` | Emit a machine-readable history |
| `--mesh-out=PATH` | Write final mesh IDs, geometry, levels, and indicators |

This driver fixes `h` and `ell`, rebuilds after each coarse refinement, and
uses legacy strong-residual candidates. Its output is useful for regression
and calibration only. Label it `HLOD-proxy`; do not label it CALOD or the frozen
HLOD paper comparator.

### Manufactured Calibration Case

The legacy exact validation problem uses

```math
\phi(t)=16t^2(1-t)^2,
\qquad
u(x,y)=\phi(x)\phi(y)e^{ikx},
```

with

```math
f=-e^{ikx}\left[
\phi''(x)\phi(y)+\phi(x)\phi''(y)
+2ik\phi'(x)\phi(y)\right].
```

Because `phi=phi'=0` at the endpoints, the homogeneous impedance condition is
satisfied. Exact error and the fixed fine-FEM discretization floor are
diagnostics only.

### Residual Identity

The legacy estimator reconstructs the integrated nodal residual from broken
volume terms, interior flux jumps, and impedance boundary terms. It must agree
with

```math
r_h=b_h-A_hu_{\mathrm{LOD}}.
```

This is a sign, edge-ownership, and conjugation gate. It is not a reliability
theorem and is not the WP3 audit-kernel `eta_H` definition.

## Evidence And Reporting Rules

- An exact solution or evaluation-reference solution is post-processing data;
  it may not influence MARK or STOP.
- Certification-audit work belongs to CALOD method time. Evaluation-reference
  work is reported separately.
- Verified and conditional evidence are different run claims. Missing or
  invalid evidence causes a structured failure under `RequireVerified`.
- A resource limit is a censored terminal state, not convergence and not a
  generic numerical failure.
- HLOD and CALOD have no silent fallback path to `HLOD-proxy`.
- Preserve configuration, provenance hashes, checkpoint compatibility data,
  transition history, indicators, bounds, work counters, timings, and the
  final structured stop code.

## What Is Not Yet Available

WP6 must connect the immutable `PaperConfig` and all paper cases/comparators to
the common output schema, frozen manifests, matched targets, and production
CLI. Separately, G3 requires the missing verified numerical evidence and
directed interval chain. Until both are complete:

- there is no supported command for a formal CALOD or frozen-HLOD paper run;
- the six-method comparison matrix has not been executed by this code path;
- no paper table or figure should be claimed from WP0-WP5 unit tests;
- exact/reference quantities must remain isolated from algorithmic decisions.

Measured calibration histories and implementation decisions belong in
[DEVELOPMENT.md](DEVELOPMENT.md), not in this user guide.

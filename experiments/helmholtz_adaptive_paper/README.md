# Helmholtz Adaptive Paper Experiment Protocol

This directory owns the versioned contracts for the Helmholtz LOD paper
experiments. Schema v1 is the legacy certified-controller contract. Schema v2
is the practical paper contract and is consumed by
`bench_helmholtz_adaptive_paper`.

## Files

| File | Role |
|---|---|
| `schema-v1.json` | Per-run input configuration contract |
| `output-schema-v1.json` | Common run-output envelope, terminal states, timings, and nullable values |
| `schema-v2.json` | Practical PALOD input contract with fixed reference/ambient epoch semantics |
| `output-schema-v2.json` | Practical `run.json` contract and artifact manifest |
| `configs/R1-palod-smoke-v2.json` | Development-only end-to-end smoke configuration |

The C++ registry in `helmholtz/experiments/paper_config.h` is the executable
counterpart. The v1 and v2 types and parsers are separate. In particular, v2
rejects legacy `theta_h`, `q_h`, and corrector-fine-refinement fields rather
than silently retaining them in the PALOD run identity.

## Identity And Provenance

Practical run IDs have the form:

```text
case_method_kN_rN_hash
```

The final component is the lowercase 16-digit FNV-1a 64-bit hash of the
canonical JSON configuration. The v2 identity covers reference/ambient mesh
policies, reference epoch, all practical decisions, resource limits, solver
selection, quadrature/tolerances, Git/build provenance, and the frozen
manuscript SHA-256.

Every formal run must retain the canonical configuration, Git revision, build
hash, compiler and dependency metadata, host/thread information, and all
method parameters. A resumed run must pass checkpoint schema, configuration
fingerprint, problem identity, and state fingerprint checks.
The v2 runner compares `manuscript_sha256` with the versioned
`MANUSCRIPT_BASELINE.sha256` artifact before doing numerical work.

## Numeric And Terminal-State Rules

Every nullable numeric result is represented as a typed pair of `value` and
`value_status`. Non-valid values use JSON `null`, never zero, NaN, infinity, or
a string sentinel.

Resource limits are explicit censored terminal states. They are distinct from
convergence and from numerical or evidence failures. CALOD and frozen HLOD use
structured stop codes and never silently fall back to the legacy H-only proxy.

## Timing Ownership

- Certification-audit assembly, solves, and certificate checks are part of
  CALOD method time.
- Evaluation-reference construction and error evaluation are recorded in
  `evaluation_reference_time_ms` and excluded from method time.
- Exact solutions and evaluation-reference values are post-processing data and
  must not be accessible to MARK or STOP logic.

## Method Names

`HLOD-proxy` is intentionally a legacy diagnostic method only. It is not a v2
paper comparator and must not be relabelled as PALOD or HLOD-fixed.

WP5 executes PALOD only. `HLOD-fixed`, `SLOD`, `UFEM`, and `AFEM` are reserved
v2 method names, but the runner rejects them until their real backends share
the same reference/error/timing contract. It never relabels an old proxy.

## Current Implementation Boundary

| Scope | Status |
|---|---|
| v1 certified input/output schemas and strict C++ parsing | Preserved as legacy |
| v2 practical input/output schemas and strict C++ parsing | Implemented in WP5 |
| Paper case registry, mixed boundaries, sources, and quadrature | Implemented in WP1 |
| Fixed reference epoch and ambient shadow hierarchy | Implemented in WP1 |
| Reference/ambient kernel Riesz problems | Implemented in WP2 |
| Ambient-to-reference retraction and localization certificate | Implemented in WP3 |
| Independent practical driver state machine | Implemented in WP4 |
| PALOD paper runner and five-file artifact contract | Implemented in WP5 |
| HLOD-fixed/SLOD/UFEM/AFEM paper backends | Pending production experiment work |
| Six-method production matrix and paper tables/figures | Not run |

Do not publish paper-result claims from the development smoke configuration.
Formal manifests still require frozen E0 calibration and completion of the
acceptance gates in
[../../HELMHOLTZ_ADAPTIVE_LOD_PLAN.md](../../HELMHOLTZ_ADAPTIVE_LOD_PLAN.md).

## Validation

The contract and paper-case regressions are included in the normal CTest suite.
Targeted validation can be run with:

```bash
cmake --build build --target \
  test_helmholtz_paper_config \
  test_helmholtz_paper_cases \
  test_helmholtz_mixed_boundary \
  test_helmholtz_quadrature -j "$(nproc)"

ctest --test-dir build \
  -R 'helmholtz_(paper_config|paper_cases|mixed_boundary|quadrature)' \
  --output-on-failure
```

The v2 end-to-end smoke is registered as
`helmholtz_adaptive_paper_v2_smoke`. A direct run has the form:

```bash
benchmarks/bench_helmholtz_adaptive_paper \
  --config=experiments/helmholtz_adaptive_paper/configs/R1-palod-smoke-v2.json \
  --output-dir=results/wp5-smoke \
  --manuscript-baseline=experiments/helmholtz_adaptive_paper/MANUSCRIPT_BASELINE.sha256
```

Each run directory contains `iterations.csv`, `summary.csv`, `run.json`,
`ell_history.csv`, and `final_mesh.vtu`. Reference solves and post-run
reference-error evaluation are timed separately and excluded from method time.

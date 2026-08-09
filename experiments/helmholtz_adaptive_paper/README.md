# Helmholtz Adaptive Paper Experiment Protocol

This directory owns the versioned contracts for the certified-adaptive
Helmholtz LOD paper experiments. It defines how a run is described and
reported; it does not yet provide the WP6 production runner or frozen formal
run matrix.

## Files

| File | Role |
|---|---|
| `schema-v1.json` | Per-run input configuration contract |
| `output-schema-v1.json` | Common run-output envelope, terminal states, timings, and nullable values |

The C++ registry in `helmholtz/experiments/paper_config.h` is the executable
counterpart. It rejects unknown fields, non-RFC-8259 numbers, invalid enum
values, and schema versions other than v1.

## Identity And Provenance

Run IDs have the form:

```text
case_method_kN_rN_hash
```

The final component is the lowercase 16-digit FNV-1a 64-bit hash of the
canonical JSON configuration. The canonical configuration includes the Git
commit, build hash, repeat index, all controller/resource limits, and every
numerical-backend discretization, solver, certificate, audit, and constant-set
content hash. Changing any of these fields creates a different immutable run
ID.

Every formal run must retain the canonical configuration, Git revision, build
hash, compiler and dependency metadata, host/thread information, and all
method parameters. A resumed run must pass checkpoint schema, configuration
fingerprint, problem identity, and state fingerprint checks.
`certificate_constant_set_hash` is a lowercase SHA-256 content digest; the WP6
runner must recompute it from the frozen constant artifact before building the
registry. The schema alone cannot attest the artifact contents.

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

`HLOD-proxy` is intentionally a diagnostic method only. It is not one of the
six paper comparators and must not be relabelled as CALOD or frozen HLOD.

Formal method output is not valid merely because it conforms to the JSON
schema. It must also come from the frozen WP6+ runner, satisfy the configured
evidence policy, and preserve the complete transition and provenance record.

## Current Implementation Boundary

| Scope | Status |
|---|---|
| v1 input/output schemas and strict C++ parsing | Implemented in WP0 |
| Paper case registry, mixed boundaries, sources, and quadrature | Implemented in WP1 |
| Three-mesh hierarchy and error-reference isolation | Implemented in WP2 |
| Audit-kernel `eta_H` | Implemented in WP3 |
| Fail-closed conditional certificate machinery | Partially implemented in WP4; verified chain remains open |
| Checkpointable certified driver state machine | Controller implemented in WP5 |
| Real floating-point numerical backend | Implemented for conditional integration smoke only |
| Verified numerical backend, frozen manifests, and formal runner | Pending G3/WP6 |
| Six-method production matrix and paper tables/figures | Not run |

Do not create formal production manifests or publish paper-result claims from
this directory until the WP6 backend and acceptance gates in
[../../HELMHOLTZ_ADAPTIVE_LOD_PLAN.md](../../HELMHOLTZ_ADAPTIVE_LOD_PLAN.md)
are complete.

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

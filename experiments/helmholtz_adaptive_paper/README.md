# Helmholtz Adaptive Paper Experiment Protocol

This directory owns the versioned contracts for the Helmholtz LOD paper
experiments. Schema v1 is the legacy certified-controller contract. Schema v2
is retained for the completed resource pilots. Schema v3 records the first
fixed-horizon calibration. Schema v4 is the active practical paper contract consumed by
`bench_helmholtz_adaptive_paper`.

## Files

| File | Role |
|---|---|
| `schema-v1.json` | Per-run input configuration contract |
| `output-schema-v1.json` | Common run-output envelope, terminal states, timings, and nullable values |
| `schema-v2.json` | Practical PALOD input contract with fixed reference/ambient epoch semantics |
| `output-schema-v2.json` | Practical `run.json` contract and artifact manifest |
| `schema-v3.json` | Historical first fixed-horizon calibration contract |
| `output-schema-v3.json` | Historical first plateau-diagnostic output |
| `schema-v4.json` | Active three-point plateau, reference-audit and trajectory-completion contract |
| `output-schema-v4.json` | Active run output with normal fixed-horizon completion semantics |
| `reference-audit-output-schema-v1.json` | Level-to-level reference adequacy audit artifact |
| `configs/R1-palod-smoke-v2.json` | Development-only end-to-end smoke configuration |
| `configs/R2a-palod-k16-resource-pilot-v2.json` | Non-paper R2a resource pilot |
| `configs/S-palod-k16-resource-pilot-v2.json` | Non-paper S resource pilot |
| `configs/R2a-palod-k16-extended-calibration-v3.json` | Non-paper fixed-horizon R2a calibration |
| `configs/S-palod-k16-extended-calibration-v3.json` | Non-paper fixed-horizon S calibration |
| `configs/R2a-palod-k16-reference-audit-v4.json` | R2a level 10-to-11 audit without rerunning PALOD |
| `configs/R2a-palod-k16-epoch1-calibration-v4.json` | Non-paper R2a epoch-1 calibration |
| `configs/S-palod-k16-step6-calibration-v4.json` | Non-paper S six-step calibration |
| `calibration/reference-epoch-v4-2026-08-12.md` | Audited epoch-2 through epoch-5 gate record; never paper data |
| `configs/R2a-palod-k16-epoch6-level16-calibration-v4.json` | Next server-only R2a calibration |
| `configs/S-palod-k16-epoch3-level13-step6-calibration-v4.json` | Next server-only S calibration |
| `configs/R2a-palod-k16-epoch6-level16-reference-audit-v4.json` | R2a level 16-to-17 audit after the next calibration |
| `configs/S-palod-k16-epoch3-level13-reference-audit-v4.json` | S level 13-to-14 audit after the next calibration |

The C++ registry in `helmholtz/experiments/paper_config.h` is the executable
counterpart. The legacy v1 and active practical types are separate. In particular, v4
rejects legacy `theta_h`, `q_h`, and corrector-fine-refinement fields rather
than silently retaining them in the PALOD run identity.

## Identity And Provenance

Practical run IDs have the form:

```text
case_method_kN_rN_hash
```

The final component is the lowercase 16-digit FNV-1a 64-bit hash of the
canonical JSON configuration. The v4 identity covers reference/ambient mesh
policies, reference epoch, all practical decisions, resource limits, solver
selection, quadrature/tolerances, Git/build provenance, and the frozen
manuscript SHA-256. It also covers the trajectory policy, independent absolute
stopping tolerance, and reporting-only plateau policy.

Every formal run must retain the canonical configuration, Git revision, build
hash, compiler and dependency metadata, host/thread information, and all
method parameters. A resumed run must pass checkpoint schema, configuration
fingerprint, problem identity, and state fingerprint checks.
The v4 runner compares `manuscript_sha256` with the versioned
`MANUSCRIPT_BASELINE.sha256` artifact before doing numerical work.

## Numeric And Terminal-State Rules

Every nullable numeric result is represented as a typed pair of `value` and
`value_status`. Non-valid values use JSON `null`, never zero, NaN, infinity, or
a string sentinel.

Resource limits are explicit censored terminal states. They are distinct from
convergence and from numerical or evidence failures. CALOD and frozen HLOD use
structured stop codes and never silently fall back to the legacy H-only proxy.

Schema v4 retains the v3 separation between two quantities that v2 accidentally conflated:

- `practical_stop_tolerance` is an absolute threshold for `U_prac` and is used
  only with `trajectory_policy=practical_indicator`;
- `relative_energy_targets` are post-processing targets and never enter
  MARK/STOP;
- `fixed_work_horizon` ignores the practical threshold and is used both for
  calibration and for pre-frozen formal trajectory acquisition. Reference-error ratios, logarithmic
  improvements, and plateau flags remain reporting-only.

A plateau uses the final three points: the two-step geometric-mean ratio must
be at least 0.9, while the whole window may oscillate by at most 15%. This
admits a small intermediate increase without confusing strong decay with a
plateau. Fixed horizons end as `TrajectoryComplete`/`success`; actual resource
limits remain censored. After every completed calibration trajectory and before
an epoch can be frozen, `bench_helmholtz_reference_adequacy` compares the
reference with exactly one further uniform NVB level. This gate is independent
of whether a plateau was observed; it never reruns or mutates the adaptive
trajectory.

## Timing Ownership

- Certification-audit assembly, solves, and certificate checks are part of
  CALOD method time.
- Evaluation-reference construction and error evaluation are recorded in
  `evaluation_reference_time_ms` and excluded from method time.
- Exact solutions and evaluation-reference values are post-processing data and
  must not be accessible to MARK or STOP logic.

## Method Names

`HLOD-proxy` is intentionally a legacy diagnostic method only. It is not a v4
paper comparator and must not be relabelled as PALOD or HLOD-fixed.

The v4 runner executes real `PALOD`, `HLOD-fixed`, `SLOD`, `UFEM`, and
`AFEM` paths.
`HLOD-fixed` reuses the reference-kernel estimator and H marking with one
frozen global ell, while skipping the PALOD localization certificate and
ambient-shadow work. `UFEM` performs uniform conforming P1 refinements and
solves, then prolongs every candidate to the same fixed reference space for
post-processing. `SLOD` rebuilds a real LOD solution along a uniform coarse
sequence with the frozen prior localization rule. `AFEM` uses volume,
interior-flux-jump, and impedance-boundary residual contributions for local
P1 refinement. The runner never relabels an old proxy.

## Current Implementation Boundary

| Scope | Status |
|---|---|
| v1 certified input/output schemas and strict C++ parsing | Preserved as legacy |
| v2 practical input/output schemas and strict C++ parsing | Implemented in WP5 |
| v3 independent stop/calibration policy and plateau reporting | Implemented after resource pilot |
| v4 oscillation-tolerant plateau, explicit trajectory completion and reference audit | Implemented after extended calibration |
| Paper case registry, mixed boundaries, sources, and quadrature | Implemented in WP1 |
| Fixed reference epoch and ambient shadow hierarchy | Implemented in WP1 |
| Reference/ambient kernel Riesz problems | Implemented in WP2 |
| Ambient-to-reference retraction and localization certificate | Implemented in WP3 |
| Independent practical driver state machine | Implemented in WP4 |
| PALOD paper runner and five-file artifact contract | Implemented in WP5 |
| HLOD-fixed paper backend without localization/certificate cost | Implemented |
| UFEM uniform conforming P1 trajectory backend | Implemented |
| SLOD/AFEM paper backends | Implemented and smoke-tested |
| Epoch-2 adequacy gates | Both failed; explicit deeper epochs required |
| R2a local epoch-3/4/5 calibrations and adjacent audits | Completed; all audits failed the 0.25 gate |
| Five-method R2a/S production matrix and paper tables/figures | Blocked; not run |

Reference solutions are cached on disk across method processes using a key
that binds the assembled reference problem, load, solver format, Git revision,
and executable build identity. Candidate errors are evaluated as each
trajectory point is produced and the candidate vector is then released; the
journal does not retain all reference-sized vectors. Both operations remain
post-processing and are excluded from method time. Cache hits and keys are
recorded in `run.json`.

Do not publish paper-result claims from the development smoke configuration.
Formal manifests still require frozen E0 calibration and completion of the
acceptance gates in
[../../HELMHOLTZ_ADAPTIVE_LOD_PLAN.md](../../HELMHOLTZ_ADAPTIVE_LOD_PLAN.md).

The case-S manufactured-solution comparison at `kappa=16`, including the
accepted run IDs, PALOD target crossings, epoch starts, and reproduction
command, is recorded in
[`S_K16_METHOD_COMPARISON_2026-08-13.md`](S_K16_METHOD_COMPARISON_2026-08-13.md).
The bounded PALOD/SLOD thread-scaling and stage-profile study, including the
post-NVB-optimization bottlenecks and recommended next optimizations, is in
[`S_K16_PALOD_SLOD_BOTTLENECK_2026-08-13.md`](S_K16_PALOD_SLOD_BOTTLENECK_2026-08-13.md).
The common-commit deep manufactured-solution convergence configurations,
server acceptance rules, old-prefix validation, and full-trajectory plotting
command are recorded in
[`S_K16_DEEP_CONVERGENCE_SERVER_2026-08-13.md`](S_K16_DEEP_CONVERGENCE_SERVER_2026-08-13.md).

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

All five v4 methods have end-to-end smoke tests, and PALOD has a separate
fixed-horizon calibration-policy smoke:
`helmholtz_adaptive_paper_v4_smoke`,
`helmholtz_hlod_fixed_paper_v4_smoke`,
`helmholtz_palod_fixed_horizon_v4_smoke`,
`helmholtz_slod_paper_v4_smoke`, `helmholtz_ufem_paper_v4_smoke`, and
`helmholtz_afem_paper_v4_smoke`. The independent audit gate is
`helmholtz_reference_adequacy_v4_smoke`. A direct PALOD run has the form:

```bash
benchmarks/bench_helmholtz_adaptive_paper \
  --config=experiments/helmholtz_adaptive_paper/configs/R1-palod-smoke-v4.json \
  --output-dir=results/wp5-smoke \
  --reference-cache-dir=results/wp5-reference-cache \
  --manuscript-baseline=experiments/helmholtz_adaptive_paper/MANUSCRIPT_BASELINE.sha256
```

Each run directory contains `iterations.csv`, `summary.csv`, `run.json`,
`ell_history.csv`, and `final_mesh.vtu`. Reference solves and post-run
reference-error evaluation are timed separately and excluded from method time.
`iterations.csv` reports both `N_H` (coarse nodes) and `DoF_H` (unconstrained
coarse degrees of freedom). For the manufactured cases R1 and S it also reports
absolute and relative weighted-energy/L2 errors against the exact solution.
These exact values use the same quadrature policy as the run and remain a
one-way post-processing quantity: they cannot affect MARK or STOP. R2a/R2b
leave the exact-error columns empty.

The two `S-palod-k16-exploratory-e1-*-step2-v4.json` configurations are short,
explicitly non-paper trajectories for checking this output before the reference
adequacy gate is frozen. Plot one or more completed epochs with:

```bash
python3 tools/visualization/plot_helmholtz_adaptive_epochs.py \
  --run-dir=results/.../epoch-0-run-id \
  --run-dir=results/.../epoch-1-run-id \
  --output=results/.../S-PALOD-error-vs-DoF.png
```

The figure marks the first evaluated point in each supplied epoch and writes a
provenance CSV beside the image. A recalibration run may restart from the same
coarse mesh, so coincident epoch-start markers are intentional and must not be
interpreted as a continued adaptive trajectory.
For resource pilots and formal server runs, follow
[../../HELMHOLTZ_ADAPTIVE_PAPER_SERVER_RUNBOOK.md](../../HELMHOLTZ_ADAPTIVE_PAPER_SERVER_RUNBOOK.md).

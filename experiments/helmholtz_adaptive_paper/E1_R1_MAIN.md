# E1 R1 main experiment

E1 uses the localized smooth manufactured case at `kappa=16`.  All four
methods start from coarse level 2, so the first plotted mesh is deliberately
coarse.  The frozen E0 values are
`theta_loc_usr=1.961323606684523` and
`C_rel_usr=1.6358263950741063`.

Production configurations:

- `E1-R1-palod-reference-epoch-k16-H2-h12-gap4-step15-v5.json`:
  reference-epoch PALOD with `ell0=2`, at most 15 committed H refinements and
  four epochs; the proactive refresh guard is a local level gap of four after
  at least three H commits, candidate is deepened to gap six before promotion,
  `ell` is inherited, and a new epoch requires four remaining solved points;
- `E1-R1-hlod-fixed-k16-H2-h15-ell3-step15-v4.json`: fixed oversampling
  `ell=3`, reference level 15;
- `E1-R1-ufem-k16-H2-level18-step16-v4.json`: uniform P1 FEM through level 18;
- `E1-R1-afem-k16-H2-level18-step40-v4.json`: adaptive P1 FEM with a
  level-18 ceiling.

The v5 PALOD path records real epoch-0 reference/coarse/candidate VTU meshes.
`plot_reference_epoch_meshes.py` reads those VTUs and aligns coarse and
candidate panels by the same iteration.  `plot_reference_epoch_e1.py` reads
the mixed v4/v5 outputs and creates exact/reference energy-error plots against
DoF and cumulative wall time.  AFEM has a manufactured exact error but no
fixed-reference error in schema v4, so it is intentionally absent only from
the reference-error panels.

The matrix-free localization Gram operator is prepared once per corrector
check.  Its output records structure preparation, patch SparseLU
factorization, RHS formation, parallel patch solve, deterministic scatter,
action count, factorization count, thread count, Ritz iterations/residual and
warm-start use.  Patch actions use OpenMP `schedule(dynamic,1)` and fixed
patch-order reduction.  Multiple Ritz basis vectors use a block action, and a
dominant vector is prolonged through each committed coarse refinement.
The production matrix-free path now uses a four-vector locally optimal block
Rayleigh--Ritz iteration with retained search directions and reuses the
already-computed block action between iterations;
the dense path remains limited to small cross-checks and fallback dimensions.

Candidate RT2 reconstruction uses OpenMP over independent vertex patches,
caches element quadrature/Piola basis/local mass/P2 data, and merges patch
contributions in fixed vertex order. Candidate dual patch Riesz solves are
also parallel. Structural non-containment and the level-gap guard decide a
reference refresh before candidate dual construction, because the dual value
cannot change either forced decision. Detailed candidate stage timings and
thread/factorization counts are written to `iterations.csv`.

An eight-H-step WSL probe of the new production settings completed in
2:36.06 with 5,128,604 KiB peak RSS and no swap. It reached `N_H=242` and
crossed the former localization failure at `N_H=161`: the block solver used
5 iterations there (relative residual `5.18e-5`) and 3 iterations at
`N_H=242` (`8.37e-5`). The gap-four refresh recorded `SkipCandidateDual`;
ordinary interval checks at gaps seven and five still computed the dual.
The new timing split also shows that at `N_H=161,N_c=20722`, 17.94 s of the
22.68 s candidate total is hierarchy enrichment/embedding work, while the
RT2 reconstruction is 4.47 s. Incremental hierarchy/embedding updates are
therefore the next optimization target, not another reduction of `theta_c`.

That hierarchy optimization is now implemented. Candidate NVB refinement
composes its step nodal/element/DG prolongations with the existing
reference-to-candidate and coarse-to-candidate maps, propagates both parent
arrays through the step parent map, preserves the fixed reference
quasi-interpolation, and rebuilds only the candidate quasi-interpolation.
Two consecutive local-refinement tests agree with a full geometric embedding
rebuild to `1e-12`. At the same `N_H=161,N_c=20722` point, candidate
enrichment fell from 17.94 s to 0.227 s (about 79x); the otherwise identical
eight-step trajectory fell from 2:36.06 to 2:06.03 (19.2%). The CSV now
separates NVB refinement, embedding composition, parent-map propagation,
candidate quasi-interpolation, and validation.

## Superseded completed PALOD trajectory

The 15-H-step production run completed locally as
`R1_PALOD-reference-epoch_k16_r0_09fa933db6d58f2f`. It contains 15 solved
points in four epochs and exits with the intended structured
`WorkLimitReached: maximum_H_steps reached`. The exact relative energy error
after the three refresh boundaries changes from `0.274717` to `0.126200`,
from `0.0539417` to `0.0380823`, and from `0.0304663` to `0.0240618`.
Consequently there is neither the superseded h10 reference floor near 0.14
nor an error rebound at a refresh. The 10% and 5% targets are reached at
`N_H=161` and `N_H=479`; 2% and 1% are outside this fixed 15-step budget.

The run used four OpenMP threads on the 12 GiB WSL machine, took 9:06.04,
peaked at 9,000,792 KiB RSS, and used no swap. In addition to incremental
embedding, the prepared saddle Gram path now discards local energy and dense
constraint copies after SparseLU factorization and discards the row-major
defect RHS after patch preparation. These objects are not read by subsequent
Gram actions; kernel-residual, reference-epoch driver, and hierarchy tests
all pass with the release enabled. Cumulative recorded costs include 235.95 s
for localization checks (129.69 s Gram factorization and 55.12 s patch
actions), 82.57 s for candidate RT2 reconstruction, 27.99 s for reference
Riesz estimation, and only 5.43 s for all candidate hierarchy enrichments.

Generated implementation-study artifacts are:

- `figures/paper/E1-R1-k16-H2-h12-gap4-main-20260823.{pdf,png,json}`;
- `figures/paper/E1-R1-k16-H2-h12-gap4-meshes-20260823/` for the real epoch-0
  reference/coarse/candidate mesh evolution.

## Superseded epoch semantics and corrected production run

The completed `09fa...` trajectory above is now retained only as diagnostic
evidence. It reset `ell` to `ell0=2` after each refresh, started later epochs
with local reference/coarse gaps of only four and two, and opened a final
epoch with one solved point. Those choices are not suitable for the revised
paper algorithm or an epoch-wise convergence fit.

The corrected production contract now:

- inherits the current `ell` across a reference refresh;
- keeps the forced refresh threshold at four, but deepens candidate locally
  before promotion until the proposed-coarse/candidate gap is at least six;
- keeps three committed H refinements per epoch and requires at least four
  remaining solve observations before a new epoch may be opened;
- reports PALOD convergence exponents separately for every epoch and never
  fits across a refresh.

Proposed-to-reference/candidate nodal, element, and DG embeddings are cached
in the coarse-refinement transaction and updated by candidate step
prolongations. Containment, level-gap evaluation, candidate closure, refresh,
and commit reuse them. In the corrected eight-step probe, epoch 1 starts at
gap six with inherited `ell=3`; total mesh time is 1.45 s, lazy-decision time
is about 0.00056 s, and the last commit is 0.076 s. The no-cross-step-Gram-
cache baseline takes 2:32.36, peaks at 5,572,280 KiB, and uses no swap.

A trajectory-owned numerical-fingerprint LRU was also tested for unchanged
reference-defect SparseLU factors. A 512-entry cache is rejected: only 33 of
715 patch lookups hit while peak RSS rises to 11,593,028 KiB. The bounded
32-entry version records 21 hits/694 misses, reduces Gram factorization from
60.70 s to 57.88 s and wall time to 2:28.43, but raises peak RSS to
6,394,072 KiB. Production therefore uses only the small cache and records
`gram_factor_cache_hits` and `gram_factor_cache_misses`; this is a modest,
not headline, optimization. RT2 cross-H-step factor reuse remains pending.

The corrected 15-step run is too large for the 12 GiB WSL machine: inherited
`ell=3` reached about 11.11 GB RSS at 4:49 with only about 0.38 GB available,
so it was stopped before swap/OOM and produced no formal output. Run the
corrected main configuration on the 366 GB server.

The fixed-LOD production value `ell=3` was selected by a frozen h12
adequacy pilot against `ell=4`. At H-step 11, the relative exact energy
errors were `0.0623348` and `0.0620395`, respectively (a 0.48% difference),
while `ell=3` used 14.18 s instead of 30.70 s and 354 MB instead of 972 MB.
Thus `ell=4` is retained as audit evidence, not used for the E1 production
curve. No `ell=2` production run is claimed without an equivalent adequacy
check.

The fixed-LOD implementation groups patches with identical local systems,
reuses their factorization, uses a bounded trajectory cache, and schedules
independent groups with OpenMP. The reference residual estimator likewise
groups identical constrained Riesz systems, uses multi-RHS solves, and uses
an SPD factorization plus a small dense Schur complement when applicable.

## Superseded diagnostic single-run results

These are implementation-study timings from the older `H2/h10/12-step`
PALOD configuration on the local 12 GiB WSL machine. They are retained for
regression and bottleneck comparison, not as the final PALOD paper curve and
not as three-repeat paper medians.

| method | terminal state | plotted terminal DoF | exact relative energy error | wall time | peak RSS |
|---|---|---:|---:|---:|---:|
| reference-epoch PALOD | `WorkLimitReached` (planned H-step limit) | 662 | 0.0955932 | 2:36.63 | 5.76 GiB |
| fixed LOD, `ell=3` | `TrajectoryComplete` | 3660 | 0.0241365 | 3:47.04 | 1.44 GiB |
| UFEM | `TrajectoryComplete` | 263169 | 0.00700805 | 2:59.64 | 1.85 GiB |
| AFEM | `TrajectoryComplete` | 61901 | 0.0136403 | 0:19.71 | 0.304 GiB |

In that superseded run, PALOD finishes with a structured work-limit row after the final refinement;
the last solved point is therefore DoF 662, while the terminal mesh has
`N_H=892`. Its first reference refresh changes the reference-error target,
so the plot deliberately breaks that curve between epochs. Stars mark
reference refreshes and crosses mark changes of `ell`.

The old measured PALOD phase totals motivated the current optimizations:
candidate flux reconstruction (75.1 s) and localization/certificate work
(`time_theta` 48.8 s, including 25.6 s of Gram factorizations and 15.8 s of
patch actions). LOD solution time itself is negligible. For fixed LOD,
corrector/model construction remains dominant (167.4/181.3 s), followed by
the grouped reference estimator (49.0 s); global operator assembly and
coarse factorization are not current bottlenecks.

Corrected 366 GB server command:

```bash
/usr/bin/time -v env OMP_NUM_THREADS=16 \
  build-release/benchmarks/bench_helmholtz_adaptive_paper \
  --config=experiments/helmholtz_adaptive_paper/configs/E1-R1-palod-reference-epoch-k16-H2-h12-gap4-step15-v5.json \
  --output-dir=results/E1-R1-k16-H2-h12-gap6-ell-inherit-main \
  --manuscript-baseline=experiments/helmholtz_adaptive_paper/MANUSCRIPT_BASELINE.sha256 \
  --check
```

Structured `WorkLimitReached` with `maximum_H_steps reached` is the intended
end of a fixed-work E1 trajectory; it is not a numerical failure.

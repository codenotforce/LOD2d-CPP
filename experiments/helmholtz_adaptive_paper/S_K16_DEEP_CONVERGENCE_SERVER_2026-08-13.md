# Case-S kappa=16 deep convergence server run

## Purpose

This fixed-horizon run extends the manufactured-solution error-versus-DoF
comparison far enough to test whether the finite-range slopes settle toward
the expected `N^{-1/2}` adaptive/PALOD and `N^{-1/3}` uniform/SLOD guides. It
does not tune marking, localization, or stopping parameters from the observed
errors. Exact errors are post-processing values and never enter MARK or STOP.

The four methods start from the same level-5 coarse mesh. SLOD keeps `ell=2`
and the same five-level h/H gap by refining H and h together. PALOD uses
scheduled reference refreshes after H steps 3, 6, and 9; every refresh
inherits the terminal coarse mesh exactly.

Thus PALOD performs 12 coarse adaptive refinements in four reference epochs:
each epoch owns three new H-refinement steps.  An epoch evaluates its inherited
start mesh and the three refined meshes, so the complete journal contains 16
PALOD solution observations on 13 distinct coarse meshes.  The repeated DoF at
each of the three epoch boundaries is intentional: it records the effect of
refreshing the reference space without changing H.  SLOD/UFEM/AFEM do not use
this PALOD reference-epoch schedule. In particular, S-AFEM now runs in
`manufactured-exact-only` mode: its current adaptive mesh is not required to be
contained in a fixed reference mesh, and its plotted exact error is integrated
directly on that mesh.

## Status of the existing comparison data

The accepted 2026-08-13 trajectories remain useful for:

- selecting safe server horizons and memory limits;
- documenting the first common PALOD-error crossing;
- checking that an extended run reproduces every existing exact-error prefix;
- demonstrating the pre-asymptotic advantage in the already resolved range.

They are not the final frozen source for a deep paper figure. The methods were
run from different intermediate commits, some templates used `WORKTREE`, the
old SLOD/AFEM timings include now-removed mesh bottlenecks, and the last-three
rate estimates are too short to establish asymptotic behavior. The deep runs
therefore rerun the complete prefix from one clean commit and one binary. On
success, the new complete trajectories replace, rather than concatenate with,
the old points in the final figure. The old results are retained as an
independent prefix regression.

For case S the plotted error is against the analytic manufactured solution,
so reference adequacy is not needed to define the plotted error. PALOD's
reference hierarchy is nevertheless part of the algorithm and all epoch
starts must remain visible. This distinction does not unblock reference-error
figures for R2a or the general formal E1 matrix.

## Frozen deep configurations

| Method | Configuration | Horizon | Intended extra range |
|---|---|---:|---|
| PALOD | `configs/S-palod-k16-deep-convergence-step12-v4.json` | 12 H steps, 4 epochs | four additional adaptive points |
| SLOD | `configs/S-slod-k16-deep-convergence-ell2-fixed-ratio-step9-v4.json` | 9 synchronized H/h steps | two additional uniform-LOD points |
| UFEM | `configs/S-ufem-k16-deep-convergence-level16-step11-v4.json` | 11 uniform steps | three additional P1 points |
| AFEM | `configs/S-afem-k16-deep-convergence-level18-step28-v4.json` | 28 adaptive steps | fourteen additional adaptive points |

The first server attempt reached `maximum_coarse_elements=300000` before the
fixed 28-step horizon. It was a clean structured stop, not a solver failure:
`/usr/bin/time` therefore reported exit status zero, while the orchestration
script correctly withheld the `.done` marker. The v2 server limits raise the
coarse-element, unknown, and ambient caps to five million and the wall cap to
12 hours. These remain hard safety limits and do not change marking or the
numerical trajectory. New runs report the exact limit name in `stop_reason`.

These are paper-candidate exact-error trajectories, not certified error-bound
runs. A work limit, allocation failure, or non-`TrajectoryComplete` terminal
state fails closed and must not be plotted as a completed curve. A reference-
containment failure remains invalid for PALOD/SLOD but is not an AFEM condition
for manufactured case S.

The templates can be parsed without starting an experiment by passing
`--validate-only` together with `--config` and `--manuscript-baseline` to
`bench_helmholtz_adaptive_paper`.

The interrupted AFEM directory cannot be resumed. Rerun the AFEM template from
its initial mesh in a new `RESULT_DIR`; do not concatenate its partial records
with the replacement. For a final four-method figure, rerun all four methods
from the new commit in one fresh result directory so that every `run.json`
records one commit and binary. Running only AFEM is acceptable as a diagnostic
repair, but then the mixed-commit set is not a frozen paper dataset.

## Server command

Use at least 64 GiB available memory; 96--128 GiB is preferred for the final
two SLOD levels. Run the four configurations serially:

```bash
CONFIGS='experiments/helmholtz_adaptive_paper/configs/S-palod-k16-deep-convergence-step12-v4.json experiments/helmholtz_adaptive_paper/configs/S-slod-k16-deep-convergence-ell2-fixed-ratio-step9-v4.json experiments/helmholtz_adaptive_paper/configs/S-ufem-k16-deep-convergence-level16-step11-v4.json experiments/helmholtz_adaptive_paper/configs/S-afem-k16-deep-convergence-level18-step28-v4.json' \
MODE=custom PATCH_THREADS=4 JOBS=16 MIN_AVAILABLE_GIB=64 \
RESULT_DIR="$PWD/results/S-k16-deep-convergence-server" \
  scripts/run_helmholtz_adaptive_paper_server.sh
```

`PATCH_THREADS=4` is a conservative launch value, not a frozen experimental
parameter.  The DoF/error figure is thread-count independent and may use any
successful deep run.  For an error-versus-time figure, record the actual
thread count and do not mix timings from different thread counts or binaries.

## Independent acceptance

Require all four `.done` files, then verify the archive contract and resource
logs:

```bash
cd results/S-k16-deep-convergence-server
sha256sum -c SHA256SUMS
grep -H -E 'Elapsed|Percent of CPU|Maximum resident|Swaps|Exit status' logs/*.time
grep -H -E 'state=|convergence_regime=|reference_cache=' logs/*.stdout
```

For every `run.json`, require `status=success`,
`driver_state=TrajectoryComplete`, and
`stop_reason="fixed H-step trajectory complete"`. The final CSV action must
be `CompleteTrajectory`; swap must be zero. Additionally require
`error_evaluation_mode="manufactured-exact-only"` for AFEM, with blank
`reference_*_error` columns and populated `relative_exact_*_error` columns.
Inspect PALOD epoch starts and the AFEM final VTU near the re-entrant corner and
mixed-boundary junctions.

## Prefix validation and consolidated plot

For each method, verify the old accepted run against the new complete run:

```bash
python3 tools/visualization/compare_helmholtz_trajectory_prefix.py \
  --baseline=results/<old-run-directory> \
  --extended=results/<new-run-directory> \
  --output=results/S-k16-deep-convergence-server/prefix-<method>.json
```

Then generate the complete curves, all PALOD epoch markers, fitted legend
rates, and both theoretical slope families:

```bash
python3 tools/visualization/plot_helmholtz_method_comparison.py \
  --palod=results/<new-PALOD-run> \
  --slod=results/<new-SLOD-run> \
  --ufem=results/<new-UFEM-run> \
  --afem=results/<new-AFEM-run> \
  --full-trajectories \
  --output=results/S-k16-deep-convergence-server/S-k16-deep-error-vs-DoF.png
```

With `--full-trajectories`, the legend reports the last-four-distinct-DoF fit;
the accompanying JSON retains the whole-range, last-three, and legend fits.
For an asymptotic statement, require that the four-point window spans a factor
of four or more and compare neighboring rolling windows. A visual crossing of
a guide or one steep window is not sufficient evidence of the limiting rate.

## Five-epoch PALOD and deeper-UFEM supplement

The accepted four-epoch PALOD run may be extended by one complete epoch with
`configs/S-palod-k16-deep-convergence-step15-v4.json`.  Its refresh schedule is
`[3,6,9,12]`, so epoch 4 inherits exactly the H mesh at step 12 and contributes
three new adaptive H refinements.  The practical runner does not splice an old
journal into a new canonical run: the supplement intentionally reruns the full
five-epoch trajectory from level 5 and the first 12 H steps must reproduce the
accepted prefix.

Case S contains the corner factor `r^(2/3)`, so quasi-uniform P1 theory predicts
an eventual energy-error rate `h^(2/3) = N^(-1/3)`.  The current kappa=16 UFEM
tail is still close to `N^(-1/2)`.  To distinguish a finite-range smooth-wave
regime from an implementation error, use
`configs/S-ufem-k16-deep-convergence-level18-step13-v4.json`.  It adds two
uniform levels, yielding a new three-point terminal slope window.  This UFEM
diagnostic is separate from the extra PALOD epoch, but should be run from the
same commit and binary.

Run both supplements serially:

```bash
CONFIGS='experiments/helmholtz_adaptive_paper/configs/S-palod-k16-deep-convergence-step15-v4.json experiments/helmholtz_adaptive_paper/configs/S-ufem-k16-deep-convergence-level18-step13-v4.json' \
MODE=custom PATCH_THREADS=4 JOBS=16 MIN_AVAILABLE_GIB=32 \
RESULT_DIR="$PWD/results/S-k16-five-epoch-and-ufem-level18" \
  scripts/run_helmholtz_adaptive_paper_server.sh
```

Require two `.done` files, no swap, and `TrajectoryComplete`.  Compare the new
PALOD run against the accepted step-12 run with the prefix checker.  Do not
replace the S deep figure merely because the final UFEM slope moves toward the
expected guide; retain the measured rolling slopes and state the observed
finite-range regime.

## Corner-dominant manufactured-solution variant

Multiplying the historical exact solution by a constant cannot change a
relative-error slope. The final diagnostic variant instead freezes the
additive decomposition

`u = a(r,theta) + B psi(x,y) exp(i kappa x)`,

where `a=chi(r) r^(2/3) sin(2 theta/3)` and
`psi=(729/16) x^2(1-x^2)^2 y^2(1-y^2)^2`.  The normalization gives
`max psi=1`.  The double zeros of `psi` on the outer square make both its
value and normal derivative vanish there, while its `x^2 y^2` factor makes it
zero on the two reentrant Dirichlet rays.  Thus the smooth oscillatory term
satisfies the same homogeneous boundary contract as the singular term.

Historical case S remains exactly its old multiplicative plane-wave case.
The additive development experiment uses singular multiplicative fraction
zero, `B=0.05`, and a C2 quintic cut-off on `(0.25,1.0)`. The wide,
endpoint-flat polynomial transition keeps the singular core unchanged while
reducing the smooth annular curvature that can dominate a finite-range global
error. It still has zero value and zero gradient on the outer boundary.
Thus the variant keeps the same L-shaped domain, reentrant exponent, wave
number, PDE operator, boundary contract, marking parameters, and relative-
error norm. The corresponding source is derived from
`-Delta u-kappa^2 u`; the `kappa^2` terms of the additive plane wave cancel
analytically. It is not an independently tuned load.

The parameters are recorded as `singular_oscillatory_fraction` and
`singular_cutoff_outer_radius`, together with the
`singular_quintic_cutoff` selector and `smooth_wave_amplitude`, in the
canonical configuration and run ID.
They may differ from their historical defaults only for case S. These runs are
a separate manufactured-solution experiment and must not be concatenated
with, or silently substituted for, the historical trajectories.

The first 24-step attempt reached `epoch=5,H_step=15` and then failed because
the cold localization power iteration did not meet its tolerance in 1,000
iterations at dimension 17,448. The preceding solve on the identical inherited
coarse space had already computed a valid dominant vector, but the epoch
refresh discarded it. The driver now preserves this coarse-coordinate warm
start across a reference-only refresh; it is still cleared whenever H changes
and the vector dimension is no longer compatible.

The revised 366 GiB comparison deliberately uses the requested smaller
horizons:

- PALOD: 18 adaptive H steps in six inherited-coarse-grid epochs, with
  refreshes after steps 3, 6, 9, 12, and 15;
- AFEM: at most 28 adaptive steps, with a 4-million-unknown fail-safe;
- SLOD: 10 synchronized H/h steps at fixed `ell=2`;
- UFEM: level 20, i.e. 15 uniform H-refinement steps from level 5 and about
  3.2 million final unknowns by extrapolation from the accepted level-16 run.

All four runs are serial. The configurations allow 16 simultaneous patch
solves. The order runs the lower-memory adaptive methods
first and leaves SLOD/UFEM, whose final sparse factorizations dominate memory,
until last. Start the frozen suite with:

```bash
MODE=s-corner-wave-366g PATCH_THREADS=16 JOBS=16 MIN_AVAILABLE_GIB=128 \
RESULT_DIR="$PWD/results/S-corner-wave-k16-six-epoch" \
  scripts/run_helmholtz_adaptive_paper_server.sh
```

The accepted level-16 UFEM run used 1.32 GiB at about 197,000 unknowns. Four
further NVB levels give about 16 times as many unknowns; level 20 therefore
keeps a substantial margin below the failed level-22 ambition. The accepted
SLOD step-9 run used about 45 GiB, so SLOD step 10 remains the expected memory
leader of this revised suite. The older step-24 PALOD, level-22 UFEM, and
step-40 AFEM templates remain archived, but the named server mode no longer
selects them.

The script checks available memory before every method, validates the build
and output contract, skips completed `.done` cases on restart, and records
`/usr/bin/time -v`. If a run is interrupted, rerun the identical command: all
completed methods are skipped. Do not lower `MIN_AVAILABLE_GIB` while another
large job is resident.

The decisive diagnostic is the terminal rolling UFEM exponent. A transition
toward `N^(-1/3)` while the nonzero wave is retained supports the controlled-
crossover design. At kappa=16 the energy contribution of a pointwise wave
amplitude scales roughly like `B*kappa`, so `B=0.05` is not negligible in the
energy norm. If the deep `B=0.05` run remains near `N^(-1/2)`, do not silently
reduce B again: either report the observed crossover or introduce an explicitly
energy-normalized wave coefficient in a new protocol.

A local level-14 corner-only (`B=0`) contract run of the quintic variant completed in about
8 seconds of method time with roughly 305 MiB peak RSS.  Its last four
point-to-point exact-energy exponents decreased monotonically from `0.510` to
`0.472`, `0.455`, and `0.432`.  This is the expected turn away from the
historical `0.5` window toward `1/3`, and is the prerequisite for launching the
additive-wave screening; it is not itself an asymptotic paper result.

The same level-14 UFEM screen was then run with the additive polynomial wave.
Because one NVB call bisects every element once, adjacent levels have a visible
odd/even effect; the robust diagnostic therefore compares points two H steps
apart. For `B=0.1`, the last four two-step exponents were `0.753`, `0.686`,
`0.649`, and `0.562`. Reducing only B to `0.05` gave `0.712`, `0.631`, `0.599`,
and `0.517`. The latter moves the crossover earlier while retaining a
nonzero kappa=16 wave, so `B=0.05` is the frozen deep-run candidate. It has not
yet reached the `1/3` asymptote at level 14.

Corner-only (`B=0`) small adaptive contract runs provide a preliminary rate check. PALOD
with 9 H steps in 3 reference epochs completed in 177.7 s at 1,584 MiB peak
RSS; after retaining the last observation at repeated epoch-boundary DoFs, its
last-four and last-three fits were `N^(-0.688)` and `N^(-0.634)`, while the
terminal pair gave `N^(-0.546)`.  AFEM with 12 steps completed in 0.52 s at
17.7 MiB; its last-four and last-three fits were `N^(-0.587)` and
`N^(-0.548)`, with terminal pair `N^(-0.487)`.  Both terminal trends are
consistent with the adaptive `N^(-1/2)` target, but the short PALOD window is
not yet a stable asymptotic fit and must not replace the deep trajectory. The
final additive `B=0.05` PALOD/AFEM rates must be recomputed rather than copied
from these corner-only runs.

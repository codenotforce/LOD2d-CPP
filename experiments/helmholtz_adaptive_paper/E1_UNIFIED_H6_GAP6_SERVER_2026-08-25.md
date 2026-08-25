# E1 unified-H6/reduced-gap main experiment (2026-08-25)

This rerun keeps the revised-paper R1 manufactured solution and Algorithm 1
unchanged.  It changes only the experimental horizon and structural reserve so
that all comparison methods start from the same H6 mesh and the PALOD reference
does not grow to the unnecessary target gap 9 used by the superseded H4 run.

## Frozen design

All five methods start at `initial_coarse_level=6`.  Adaptive methods use
`theta_H=0.1`; this lengthens the DoF trajectory without changing any
estimator.  Exact manufactured-solution errors remain validation-only.

- reference-epoch PALOD: H6/h12, `ell0=2`, inherited `ell`,
  `theta_c=0.3`, proactive trigger/target gaps 2/6, DirectSchur with its exact
  small-patch saddle fallback, 36 H-step budget and `tol_ref=0.005`;
- fixed LOD: H6/h18, fixed `ell=3`, `theta_H=0.1`, 36 H-steps;
- standard uniform LOD: H6/h10, fixed gap 4 and `ell=3`, ten synchronous
  refinements ending at H16/h20;
- UFEM: H6 through uniform level 20, fourteen refinements;
- AFEM: H6, `theta_H=0.1`, adaptive-level cap 20 and 160-step safety budget.

The latest paper still lists PALOD, fixed LOD, UFEM, and AFEM for E1.  SLOD is
an additional standard-LOD comparator and must be labelled as such; it does not
replace fixed LOD.

## Optimization evidence

The old H4/target-gap9 PALOD main reached 436451 reference unknowns, used
315573716 KiB peak RSS and 1:56:00 wall time.  The new local gates used the same
mathematics and current code:

- trigger/target 3/6, ten-step budget: 23.12 s, 4902564 KiB, zero swap;
- trigger/target 2/6, twelve-step budget: 22.57 s, 4144700 KiB, zero swap.

The latter produced a five-point epoch followed by a four-point epoch.  Their
exact-error/DoF exponents were 3.38 and 1.51, so the reduced reserve did not
cause a rate loss.  Trigger 2 is retained because trigger 3 admitted a
two-point intermediate epoch.

For SLOD, identical H6/h10 five-step probes gave:

| ell | wall | peak RSS | tail-four exact exponent |
|---:|---:|---:|---:|
| 3 | 12.68 s | 766296 KiB | 0.7544 |
| 4 | 24.38 s | 1250572 KiB | 0.7552 |

The maximum pointwise relative difference of the exact energy errors was
0.137%.  Hence ell 3 is frozen for the main SLOD curve.  The implementation
already uses exact-only validation for manufactured SLOD, incremental uniform
hierarchy composition, identical-patch factor reuse, and automatic DirectSchur
after 24000 reference DoFs.

AFEM and UFEM were also run concurrently as a non-timing local gate.  Their
tail-four exponents were 0.532 and 0.518.  Production timing remains serial so
that CPU contention cannot bias wall time or peak RSS.

## Server preparation

```bash
ssh shuihan
cd /home/sutai/code/LOD2d-CPP
git switch codex/palod-streaming-gram
git pull --ff-only origin codex/palod-streaming-gram
git status --short
git rev-parse HEAD
free -h
df -h "$PWD"
```

## Short server gate

This runs PALOD trigger/target 2/6 and the two SLOD oversampling probes from
short to long:

```bash
MODE=e1-revised-h6-pilot VALIDATE=1 JOBS=16 PATCH_THREADS=16 \
RESULT_DIR="$PWD/results/E1-H6-gap6-pilot-$(git rev-parse --short HEAD)-$(date +%Y%m%d-%H%M%S)" \
  bash scripts/run_helmholtz_adaptive_paper_server.sh
```

The gate requires zero swap, a valid three-or-more-point post-initial PALOD
epoch, exact-error decay, corrector accounting, ell inheritance and a fixed
epoch-0 reference mesh.

## Candidate-optimized medium gate

Before applying candidate batching or closure-cost-aware marking to the long
PALOD trajectory, run their PALOD-only 16-step gate.  This uses stride 2,
forces a full RT2 sweep before every dual/refresh/termination decision, retains
the original global Dörfler mass with `theta_c=0.3`, and uses an NVB closure-cost
candidate pool of size `2m`.

Do not run it concurrently with E2 moving-reference PALOD.  The server mode
requires at least 96 GiB available memory and zero configured experiment swap.

```bash
cd /home/sutai/code/LOD2d-CPP
RESULT_DIR="$PWD/results/E1-candidate-optimized-gate-$(git rev-parse --short HEAD)-$(date +%Y%m%d-%H%M%S)"
MODE=e1-candidate-optimized-gate VALIDATE=1 JOBS=16 PATCH_THREADS=16 \
RESULT_DIR="$RESULT_DIR" \
  bash scripts/run_helmholtz_adaptive_paper_server.sh
```

Promotion to the 36-step main configuration requires all of:

- zero swap and peak RSS within the server gate;
- at least one post-initial epoch with three or more solved points;
- no cross-refresh rate fit;
- final exact error within 1% of the stride-1 policy at a comparable coarse
  DoF;
- candidate closure additions at least 20% below the stride-1 policy;
- no failed localization, reference stability, nestedness or Dörfler audit.

The implementation and ten-step ablation evidence are recorded in
`E1_CANDIDATE_CLOSURE_ABLATION_2026-08-25.md`.  Promotion-time candidate rebuild
and active-region gap reserve are not part of this gate.

## Main experiment

The production runner is deliberately serial and ordered as AFEM, UFEM, SLOD,
fixed LOD, PALOD.  AFEM and UFEM were already exercised concurrently during
the local non-timing gate; serial production is necessary for interpretable
paper timings and avoids concurrent access to the common build tree.

```bash
cd /home/sutai/code/LOD2d-CPP
SESSION=e1-h6-gap6-main
RESULT_DIR="$PWD/results/E1-H6-gap6-main-$(git rev-parse --short HEAD)-$(date +%Y%m%d-%H%M%S)"

tmux new-session -d -s "$SESSION" \
  "cd '$PWD' && systemd-run --user --scope -p MemoryMax=340G -p MemorySwapMax=0 env MODE=e1-revised-h6-main VALIDATE=1 JOBS=16 PATCH_THREADS=16 RESULT_DIR='$RESULT_DIR' bash scripts/run_helmholtz_adaptive_paper_server.sh 2>&1 | tee '$RESULT_DIR.launch.log'"
```

The mode applies a 300 GiB available-memory gate before every method and a
24-hour external timeout to each configuration.  It does not delete earlier
results.

## Monitoring

```bash
ROOT=/absolute/path/to/the/result-directory
watch -n 5 'free -h; pgrep -af bench_helmholtz_adaptive_paper'
find "$ROOT/runs" -name '*.done' -print
grep -H -E 'Elapsed|Percent of CPU|Maximum resident|Swaps|Exit status' \
  "$ROOT"/logs/*.time
tail -F "$ROOT"/logs/*.stdout
```

For PALOD progress:

```bash
RUN=$(find "$ROOT/runs" -path '*PALOD*' -name iterations.csv -print -quit)
tail -n 8 "$RUN"
```

## Post-processing

Do not plot cumulative wall time.  Plot exact error versus DoF for all five
methods and fit PALOD only within each epoch:

```bash
python3 tools/visualization/plot_reference_epoch_e1.py \
  --experiment E1 \
  --palod /path/to/PALOD/run-dir \
  --fixed-lod /path/to/fixed-LOD/run-dir \
  --slod /path/to/standard-LOD/run-dir \
  --ufem /path/to/UFEM/run-dir \
  --afem /path/to/AFEM/run-dir \
  --output figures/paper/E1-unified-H6-gap6-main
```

The no-rate-loss gate is an exact-error/DoF exponent of at least 0.4 on the
last PALOD epoch containing at least three solved points.  No fit may cross a
reference refresh.  Pre-asymptotic points remain visible.

After downloading the VTU payload, render deduplicated mesh states with:

```bash
python3 tools/visualization/plot_reference_epoch_meshes.py \
  --run-dir /path/to/PALOD/run-dir \
  --output-dir figures/paper/E1-unified-H6-gap6-main-meshes \
  --epochs 0,1,2 --all-checkpoints --deduplicate-identical \
  --experiment-label E1
```

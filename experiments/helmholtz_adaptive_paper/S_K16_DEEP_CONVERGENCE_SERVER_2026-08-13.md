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

These are paper-candidate exact-error trajectories, not certified error-bound
runs. A work limit, reference-containment failure, allocation failure, or
non-`TrajectoryComplete` terminal state fails closed and must not be plotted as
a completed curve.

The templates can be parsed without starting an experiment by passing
`--validate-only` together with `--config` and `--manuscript-baseline` to
`bench_helmholtz_adaptive_paper`.

## Server command

Use at least 64 GiB available memory; 96--128 GiB is preferred for the final
two SLOD levels. Run the four configurations serially:

```bash
CONFIGS='experiments/helmholtz_adaptive_paper/configs/S-palod-k16-deep-convergence-step12-v4.json experiments/helmholtz_adaptive_paper/configs/S-slod-k16-deep-convergence-ell2-fixed-ratio-step9-v4.json experiments/helmholtz_adaptive_paper/configs/S-ufem-k16-deep-convergence-level16-step11-v4.json experiments/helmholtz_adaptive_paper/configs/S-afem-k16-deep-convergence-level18-step28-v4.json' \
MODE=custom PATCH_THREADS=4 JOBS=16 MIN_AVAILABLE_GIB=64 \
RESULT_DIR="$PWD/results/S-k16-deep-convergence-server" \
  scripts/run_helmholtz_adaptive_paper_server.sh
```

Four threads preserve the existing server resource contract. Because the
PALOD Riesz stage is now parallel, separately repeat the bounded performance
configs at 4/8/16 threads before producing a new error-versus-time figure.
Do not mix timings from different thread counts or old binaries. The DoF/error
figure itself may use the successful four-thread deep run.

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
be `CompleteTrajectory`; swap must be zero. Inspect PALOD epoch starts and the
AFEM final VTU near the re-entrant corner and mixed-boundary junctions.

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

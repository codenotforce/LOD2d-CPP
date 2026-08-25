# E2 cut-off-free L-shaped experiment (2026-08-25)

This runbook reproduces the revised E2 manufactured-solution experiment from
`helmholtz_lod_certified_amsart_revised.tex`.  The frozen manuscript SHA-256 is

```text
b7fa79f0f7c90d5a9dcb2d2938635e0ed3ccb13224f4adc4dcd9782954fb7874
```

## Exact solution and scientific check

The former radial cut-off is not used.  With

```text
b(x,y) = (1-x^2)^2 (1-y^2)^2,
s(r,theta) = r^(2/3) sin(2 theta/3),
(x_osc,y_osc) = (-1/2,1/2), alpha_osc = 25, epsilon_osc = 0.25,
```

the manufactured solution is

```text
u = b s + epsilon_osc C_osc x y b
      exp(-alpha_osc ((x-x_osc)^2+(y-y_osc)^2))
      exp(i kappa (x-x_osc)).
```

The forcing is evaluated analytically.  Unit tests check the mixed homogeneous
boundary data, finite-difference gradients, the PDE residual away from the
corner, and the absence of a radial transition layer.  A valid adaptive mesh
must show grading at the reentrant corner and around the Gaussian oscillation
center, not a ring at the old cut-off radius.

## Frozen comparison

All methods start at coarse level H6 and use `theta_H=0.2`.  This is deliberately
less aggressive than the former 0.3 marking, giving more points without losing
the observed asymptotic slope.

| Method | Reference policy | ell | Candidate marking | Main horizon |
|---|---|---:|---:|---:|
| AFEM | none; exact manufactured error | n/a | n/a | up to local level 24 |
| Moving PALOD | H6/h10, immediate promotion, `R_star=0.125` | 2 | global Dörfler, `theta_c=0.2` | resource-capped, at most 24 H-steps |
| Standard LOD (SLOD) | uniform gap 4 | 2 | uniform H/h refinement | resource-limited, at most 10 H-steps |
| Standard PALOD | H6/h12, refresh gap 2, target gap 6 | inherited, starts at 2 | global Dörfler, `theta_c=0.2` | 24 H-steps |

The SLOD horizon is also bounded by 8 million ambient elements and 2 million
coarse elements.  Thus `maximum_H_steps=10` is a ceiling, not a request to build
an infeasible tenth uniform level.

The production quadrature uses triangle orders 12/16/24 and at most six
recursive subdivisions.  Against 16/20/32 with eight subdivisions, the local
SLOD pilot reduced wall time by about 30%, while the largest relative exact
energy-error change was approximately `5.1e-6`.  Further verified reductions:

- SLOD ell 3 -> 2: about 25% lower total wall time and 37% lower peak memory;
- fixed-LOD ell 3 -> 2: about 2x faster in the six-step probe;
- PALOD `theta_c` 0.3 -> 0.2: lower candidate work, with no observed tail-rate
  degradation;
- moving PALOD `R_star=0.125`: 8 steps in 19.1 s and 1.30 GiB locally, so the
  manuscript radius is retained.

Standard PALOD remains dominated by candidate RT2 flux reconstruction.  It is
placed last.  The production code already uses patch OpenMP, summary audit,
direct Schur patch solves, reduced quadrature and ell 2; no mathematical
shortcut is enabled.

## Resource-safe moving-PALOD continuation

The original 36-step moving-PALOD horizon is superseded.  Its reference mesh
reached about 2.08 million elements at H-step 33, consumed about 351 GiB, and
spent most late-stage time in a serial global reference factorization.  The
optimized continuation therefore uses:

- a pre-factorization guard at 900,000 reference and 1,100,000 candidate
  unconstrained unknowns;
- atomic `iterations.partial.csv` and `progress.json` checkpoints after every
  solved point, reference promotion, and terminal transition;
- deterministic OpenMP load quadrature followed by an element-order scatter;
- UMFPACK for the periodic global reference solve;
- `reference_validation_stride=2`: every exact PALOD error is retained, while
  the independent reference solution/error is evaluated every second H-step;
- closure-cost-aware global candidate Dörfler marking with pool factor 3;
- the original direct-saddle corrector backend, because the local 8-step
  direct-Schur ablation was slower and used more memory.

On the 12-GiB WSL gate, deterministic parallel load assembly reduced the
selected 8-step trajectory from 20.15 s to 17.06 s (about 15%), with zero swap
and a final exact relative energy error change below `2e-16`.  UMFPACK roughly
halved numeric factorization time at this small scale; its intended benefit is
the much larger late reference systems.

Run the moving method alone in three gates; do not launch the next gate until
the preceding tail fit is accepted:

```bash
MODE=e2-cutofffree-revised-optimized-gate \
RESULT_DIR="$PWD/results/E2-cutofffree-optimized-gate-$(git rev-parse --short HEAD)" \
JOBS=16 PATCH_THREADS=16 \
bash scripts/run_helmholtz_adaptive_paper_server.sh

MODE=e2-cutofffree-revised-optimized-medium \
RESULT_DIR="$PWD/results/E2-cutofffree-optimized-medium-$(git rev-parse --short HEAD)" \
JOBS=16 PATCH_THREADS=16 \
bash scripts/run_helmholtz_adaptive_paper_server.sh

MODE=e2-cutofffree-revised-optimized-main \
RESULT_DIR="$PWD/results/E2-cutofffree-optimized-main-$(git rev-parse --short HEAD)" \
JOBS=16 PATCH_THREADS=16 \
bash scripts/run_helmholtz_adaptive_paper_server.sh
```

For the 16-step gate, fit the final 6 and 8 exact-error points against `N_H`.
Proceed to 24 steps only when both fits decay, the final-six Pearson
correlation is at most `-0.98`, no late point increases by more than 2%, and
swap use is zero.  The 24-step paper fit is accepted when the final 8 and 12
point exponents are stable to within 0.15 and the final-8 exponent is at least
0.45.  A configured unknown limit is a valid censored stop, not permission to
extrapolate beyond the available data.

```bash
python3 tools/analysis/analyze_e2_tail.py \
  --run-dir results/E2-cutofffree-optimized-medium-$(git rev-parse --short HEAD) \
  --gate medium \
  --output results/E2-cutofffree-optimized-medium-tail.json

python3 tools/analysis/analyze_e2_tail.py \
  --run-dir results/E2-cutofffree-optimized-main-$(git rev-parse --short HEAD) \
  --gate main \
  --output results/E2-cutofffree-optimized-main-tail.json
```

## Server checkout and pilot

Do not update the shared checkout while another benchmark from that checkout is
running.  After the existing E1 job finishes:

```bash
cd /home/sutai/code/LOD2d-CPP
git fetch origin
git switch codex/palod-streaming-gram
git pull --ff-only origin codex/palod-streaming-gram
git rev-parse HEAD
```

Run the four pilot/configuration gates in shortest-to-longest order:

```bash
tmux new-session -d -s e2-cutofffree-pilot \
  "cd /home/sutai/code/LOD2d-CPP && \
   MODE=e2-cutofffree-revised-pilot \
   RESULT_DIR=/home/sutai/code/LOD2d-CPP/results/E2-cutofffree-pilot-$(git rev-parse --short HEAD)-20260825 \
   JOBS=16 PATCH_THREADS=16 \
   bash scripts/run_helmholtz_adaptive_paper_server.sh \
   > results/E2-cutofffree-pilot-launch.log 2>&1"
```

Monitor with

```bash
tmux capture-pane -pt e2-cutofffree-pilot -S -80
pgrep -af bench_helmholtz_adaptive_paper
find results/E2-cutofffree-pilot-*-20260825/runs -name '*.done' -print
grep -H -E 'Elapsed|Percent of CPU|Maximum resident|Swaps|Exit status' \
  results/E2-cutofffree-pilot-*-20260825/logs/*.time
```

The pilot is accepted only if every `.done` file exists, swap use is zero, the
moving-PALOD exact error has net decay, `ell_S >= ell`, and the final AFEM and
moving-PALOD meshes contain no annular refinement artifact.

## Main run

Use a fresh result directory after the pilot passes:

```bash
tmux new-session -d -s e2-cutofffree-main \
  "cd /home/sutai/code/LOD2d-CPP && \
   MODE=e2-cutofffree-revised-main \
   RESULT_DIR=/home/sutai/code/LOD2d-CPP/results/E2-cutofffree-main-$(git rev-parse --short HEAD)-20260825 \
   JOBS=16 PATCH_THREADS=16 \
   bash scripts/run_helmholtz_adaptive_paper_server.sh \
   > results/E2-cutofffree-main-launch.log 2>&1"
```

Runs are intentionally serial.  Although AFEM and moving PALOD can coexist in
memory, concurrent execution would contaminate method wall times and CPU/RSS
diagnostics.  Patch-level parallelism remains enabled within each LOD method.

## Plotting after the main run

Resolve the single run directory below each config stem, then generate the
error-versus-DoF figure (there is deliberately no cumulative-wall-time panel):

```bash
ROOT=results/E2-cutofffree-main-$(git rev-parse --short HEAD)-20260825
AFEM=$(find "$ROOT/runs/E2-S-afem-cutofffree-k16-H6-level24-thetaH02-step120-v4" -mindepth 1 -maxdepth 1 -type d | head -1)
MOVING=$(find "$ROOT/runs/E2-S-palod-moving-cutofffree-k16-H6-h10-radius0125-thetaH02-thetaC02-step36-v6" -mindepth 1 -maxdepth 1 -type d | head -1)
SLOD=$(find "$ROOT/runs/E2-S-slod-cutofffree-k16-H6-h10-gap4-ell2-step10-v4" -mindepth 1 -maxdepth 1 -type d | head -1)
STANDARD=$(find "$ROOT/runs/E2-S-palod-standard-cutofffree-k16-H6-h12-trigger2-target6-schur-thetaH02-thetaC02-step24-v6" -mindepth 1 -maxdepth 1 -type d | head -1)

python3 tools/visualization/plot_reference_epoch_e1.py \
  --experiment E2 \
  --hybrid "$MOVING" \
  --standard "$STANDARD" \
  --slod "$SLOD" \
  --afem "$AFEM" \
  --initial-coarse-level 6 \
  --output figures/paper/E2-cutofffree-main
```

Render distinct coarse/reference/candidate states for the first three moving
epochs and suppress identical meshes:

```bash
python3 tools/visualization/plot_reference_epoch_meshes.py \
  --run-dir "$MOVING" \
  --output-dir figures/paper/E2-cutofffree-main-moving-meshes \
  --epochs 0,1,2 \
  --all-checkpoints \
  --deduplicate-identical \
  --experiment-label E2
```

Render the final AFEM mesh separately:

```bash
python3 tools/visualization/plot_vtu_mesh.py \
  --mesh "$AFEM/final_mesh.vtu" \
  --output figures/paper/E2-cutofffree-main-AFEM-final \
  --title 'E2 AFEM final mesh'
```

Do not infer the claimed convergence rate from a cross-epoch PALOD fit.  Report
standard-PALOD fits epoch by epoch, and use a tail fit for moving PALOD, SLOD and
AFEM.  Retain all pre-asymptotic points in the figure.

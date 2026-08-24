# E2 unified-H6/global-candidate-mark experiment (2026-08-25)

This experiment replaces the earlier H3 comparison at commit `08f7e65`.
Production code and frozen templates start at commit `6528b81` on
`codex/palod-streaming-gram`.

## Frozen design

All four methods start from `initial_coarse_level=6` and use
`theta_H=0.1`.  This smaller Dörfler fraction produces a longer DoF trajectory
without changing the estimator or the exact-error evaluation.

- Moving PALOD: `h12` (initial level gap 6), `R_*=0.0625`, `ell0=2`,
  `ell_max=4`, and `theta_c=0.3`.  Candidate selection is one global Dörfler
  marking.  The reported F/R masses only partition that global marked set and
  are not two independent bulk constraints.
- Standard PALOD: `h12`, refresh trigger gap 3, refresh target gap 6,
  `ell0=2`, `ell_max=6`, and `theta_c=0.3`.  Correctors use DirectSchur with
  identical-support factorization reuse; patches for which Schur is not
  profitable retain the exact automatic DirectSaddle fallback.
- AFEM: H6, `theta_H=0.1`, at most 60 adaptive steps, evaluation reference
  level 20.
- Fixed LOD: H6, fixed `ell=3`, `theta_H=0.1`, reference level 16, at most
  24 adaptive steps.

The moving and standard main trajectories contain at most 24 H-steps.  The
server script orders the methods as AFEM, moving PALOD, fixed LOD, and standard
PALOD.

## Optimization evidence gate

The three-run solver probe is
`MODE=e2-unified-solver-probe`.  On the 2026-08-25 server run, DirectSchur and
DirectSaddle generated identical three-point standard-PALOD errors.  Schur
reduced corrector time from 0.828 s to 0.756 s and LOD-build time from 1.225 s
to 1.166 s.  Total wall time was still dominated by shared validation and was
12.06 s versus 12.02 s.  Schur is therefore retained only with the existing
small-patch saddle fallback.

The eight-step gate is `MODE=e2-unified-pilot`.  Both runs completed without
swap or failed states.  Moving PALOD used 56.75 s and 7.72 GiB peak RSS;
standard PALOD used 47.19 s and 12.63 GiB peak RSS.  Exact relative energy
errors were strictly decreasing:

- moving: `0.53407 -> 0.49131 -> 0.43453 -> 0.36840 -> 0.33710 -> 0.30682 -> 0.27030 -> 0.23670`;
- standard: `0.56909 -> 0.55144 -> 0.51657 -> 0.49077 -> 0.47146 -> 0.43114 -> 0.37776 -> 0.33971`.

The standard eight-step method time fell from 824.0 s in the old H3/gap9 run
to 40.3 s in H6/gap6, while its reference size fell from 80,952 to 17,574
unknowns.  This is the performance gate for the main experiment.  The short
pilot log-log slopes are diagnostic only because its DoF range is narrow; the
paper rate must be fitted on the completed main trajectory.

## Server commands

```bash
ssh shuihan
cd /home/sutai/code/LOD2d-CPP
git switch codex/palod-streaming-gram
git pull --ff-only origin codex/palod-streaming-gram
```

Solver probe:

```bash
MODE=e2-unified-solver-probe VALIDATE=1 JOBS=16 PATCH_THREADS=16 \
RESULT_DIR="$PWD/results/E2-unified-solver-probe-$(git rev-parse --short HEAD)-$(date +%Y%m%d-%H%M%S)" \
  bash scripts/run_helmholtz_adaptive_paper_server.sh
```

Eight-step gate:

```bash
MODE=e2-unified-pilot VALIDATE=1 JOBS=16 PATCH_THREADS=16 \
RESULT_DIR="$PWD/results/E2-unified-pilot-$(git rev-parse --short HEAD)-$(date +%Y%m%d-%H%M%S)" \
  bash scripts/run_helmholtz_adaptive_paper_server.sh
```

Main run:

```bash
tmux new -s e2-unified-main
MODE=e2-unified-main VALIDATE=1 JOBS=16 PATCH_THREADS=16 \
RESULT_DIR="$PWD/results/E2-unified-main-$(git rev-parse --short HEAD)-$(date +%Y%m%d-%H%M%S)" \
  bash scripts/run_helmholtz_adaptive_paper_server.sh
```

Monitor:

```bash
ROOT="$PWD/results/E2-unified-main-6528b81-server-20260825"
find "$ROOT/runs" -name '*.done' -print
tail -F "$ROOT"/logs/*.stdout
grep -H -E 'Elapsed|Percent of CPU|Maximum resident|Swaps|Exit status' \
  "$ROOT"/logs/*.time
```

Do not accept the main result solely because every process exits zero.  Verify
strict exact-error decay, global candidate bulk mass, moving promotion count,
reference/candidate resource growth, and fitted DoF rates before producing the
paper plot.

## Main result audit

The frozen main run is
`results/E2-unified-main-6528b81-server-20260825` on `shuihan`.  Its complete
`SHA256SUMS` passed.  A 53 MiB transport archive was also frozen on the server
with SHA-256
`c4d2a66bcd09874fb05a304328046e873bf64ea7fb87fdd03fc7c56c33016e58`.
The complete 232 MiB directory is retained on the server, in the WSL mirror,
and in the Windows workspace.  Both local copies pass the server-generated
`SHA256SUMS`; each contains 174 VTU mesh files.

All four jobs exited zero without swap:

| method | solved points | terminal DoF | terminal exact energy error | wall time | peak RSS |
| --- | ---: | ---: | ---: | ---: | ---: |
| AFEM | 61 | 1907 | 0.10379 | 25.61 s | 13 MiB |
| Moving PALOD | 24 | 2123 | 0.05343 | 7:36.30 | 35.76 GiB |
| fixed LOD, ell=3 | 25 | 1646 | 0.06852 | 25:53.48 | 26.40 GiB |
| standard PALOD | 22 | 1244 | 0.08796 | 8:30.89 | 64.07 GiB |

The standard run stopped with `WorkLimitReached` because the remaining budget
could not supply the configured minimum number of points in a new epoch.  This
is the intended protection against a one-point terminal epoch, not a failure.
It used epochs 0--8 and retained `ell=2` throughout.  Moving PALOD also retained
`ell=2`; every nonterminal candidate was promoted, giving epochs 0--23.

Tail-four fitted exact-error exponents `p` in `error ~ DoF^-p` are:

- Moving PALOD: 1.051;
- standard PALOD: 0.983 (fitted inside terminal epoch 8);
- fixed LOD: 0.843;
- AFEM: 0.633.

These finite-trajectory exponents show no rate degradation, but they are steeper
than the asymptotic low-regularity rates and must not be reported as measured
theoretical orders without a deeper range.  At common errors from 0.3 to 0.1,
Moving PALOD is not the least-DoF method.  For example, at error 0.1 the first
recorded DoFs are 1124 (fixed LOD), 1137 (standard PALOD), and 1349 (Moving
PALOD); AFEM does not reach 0.1 in its budget.  Moving PALOD is the only run that
reaches 0.06, at DoF 1944, because the other trajectories terminate earlier.

The four configurations share uniform base level H6 and `theta_H=0.1`, but the
algorithm-specific singular conformity collar is applied before the first
Moving-PALOD solve.  Its first recorded coarse size is therefore 520 rather
than 208.  This distinction must be stated when claiming a unified initial
configuration.  With `R_*=0.0625` and one global `theta_c=0.3` Dörfler mark, all
23 candidate selections satisfy a global mass ratio in
`[0.300007, 0.301289]`.  None of the directly selected cells lies in
`Omega_F`; the regular-region indicator dominates after halving the singular
radius.  Matching closure still produces the visible annular candidate update.

The dominant measured costs are:

- Moving PALOD: reference promotion/mesh work 152.3 s, candidate enrichment
  102.7 s (85.96 s RT2 reconstruction), localization certificate 68.0 s, and
  reference Riesz solves 52.2 s;
- standard PALOD: candidate enrichment 160.3 s (148.5 s RT2 reconstruction),
  localization certificate 76.0 s, reference Riesz solves 58.4 s, and LOD
  builds 36.8 s;
- fixed LOD: 769.0 s in model construction, including 641.2 s in correctors,
  plus 552.7 s in estimator work.

The H6/gap6 standard configuration reduced eight-step method time from 824.0 s
in the superseded H3/gap9 run to 40.3 s (20.4x), primarily by controlling
reference/patch size.  The DirectSchur probe itself supplied a smaller 4.8%
LOD-build and 8.7% corrector improvement at three steps; the large total gain
must not be attributed to the solver alone.

Final local figures:

- `figures/paper/E2-unified-H6-global-mark-6528b81-20260825.{png,pdf,json}`;
- `figures/paper/E2-unified-H6-global-mark-6528b81-epochs0-1-meshes-20260825/`;
- `figures/paper/E2-unified-H6-global-mark-6528b81-epochs0-2-unique-meshes-20260825/`;
- `figures/paper/E2-unified-H6-global-mark-6528b81-last-epoch23-unique-meshes-20260825/`;
- `figures/paper/E2-unified-H6-global-mark-6528b81-AFEM-final-mesh-20260825.{png,pdf}`.

The two-epoch mesh audit confirms that each reference mesh is bitwise unchanged
inside its epoch and that the epoch-0 candidate SHA-256 equals the epoch-1
reference SHA-256 after promotion.

The AFEM run already contains `final_mesh.vtu`; no recomputation was needed to
render its H-step-60 mesh (3788 cells).  Moving-PALOD epoch 23 has no mesh
change after initialization: its candidate and final snapshots are bitwise
identical to the 138494-cell reference, so the deduplicated last-epoch figure
draws that geometry once together with the 4242-cell coarse mesh.

The three-epoch `unique-meshes` audit was regenerated from the complete server
VTU payload.  It writes one page for each of epochs 0, 1, and 2 and renders the
coarse, reference, and candidate groups.  Repeated reference checkpoints and
an epoch-start candidate that is bitwise identical to the reference are listed
in the JSON audit but are not drawn twice.  Reproduce it with:

```bash
python3 tools/visualization/plot_reference_epoch_meshes.py \
  --run-dir "$RUN_DIR" \
  --output-dir figures/paper/E2-unified-H6-global-mark-6528b81-epochs0-2-unique-meshes-20260825 \
  --epochs 0,1,2 --all-checkpoints --deduplicate-identical \
  --experiment-label E2
```

## WSL SCP MTU fix

If SSH commands and approximately 1 KiB files work but larger SCP/SFTP/rsync
transfers create a zero-byte destination and stall, inspect the actual route and
interface MTU:

```bash
ip route get 47.97.230.226
ip -d link show <route-interface>
ping -c 2 -M do -s 1472 47.97.230.226
```

On 2026-08-25, mirrored WSL networking routed the server through `eth4`, which
inherited MTU 9000 from the Windows `Meta` virtual adapter.  The public path
supported only MTU 1500.  SSH control traffic and a 1 KiB file passed, while an
8 KiB random file stalled immediately after the SCP file header.  This is a
PMTUD black hole, not a VTU, filesystem-permission, or server-storage problem.

The immediate WSL fix, run from Windows PowerShell, is:

```powershell
wsl -d Ubuntu-22.04 -u root -- ip link set dev eth4 mtu 1500
```

Re-run `ip route get` after a WSL/network restart because the route interface
name can change.  The persistent host-side fix requires an elevated PowerShell
or Command Prompt and should target the adapter reported by
`Get-NetIPInterface`:

```powershell
netsh interface ipv4 set subinterface "Meta" mtu=1500 store=persistent
```

Validate the repair with a non-compressible file before downloading a large
result directory.  A successful transfer must match the server hash; file size
alone is insufficient.

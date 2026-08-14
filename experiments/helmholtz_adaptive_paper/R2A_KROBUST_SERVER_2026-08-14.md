# R2a wave-number robustness server experiment

## Question and comparison

This development experiment tests whether the single-reference-epoch PALOD
relative energy-error curves remain comparable for `kappa = 2, 4, 8, 16, 32`,
and contrasts them with conforming residual AFEM. It is intended to expose a
longer AFEM pre-asymptotic/pollution range as the wave number grows. It does
not assume that the coarse meshes satisfy a theoretical `kappa H` resolution
constant. Reference errors are evaluation-only and never enter marking or
stopping.

All methods start from uniform NVB level 4, use `theta_H=0.5`, and run in one
fixed reference epoch. PALOD uses the same `theta_loc`, user constants, and
`ell` bounds at every wave number. AFEM uses the same residual estimator and
marking parameter. Thus a wave-number trend is not produced by retuning an
algorithm for each curve.

## Frozen development matrix

| kappa | reference level | audit level | PALOD cap | AFEM cap |
|---:|---:|---:|---:|---:|
| 2 | 18 | 19 | 8 H steps | 8 H steps |
| 4 | 18 | 19 | 8 H steps | 8 H steps |
| 8 | 18 | 19 | 8 H steps | 8 H steps |
| 16 | 18 | 19 | 8 H steps | 8 H steps |
| 32 | 18 | 19 | 8 H steps | 8 H steps |

Every trajectory has `minimum_reference_level_gap=3` as a hard reference-
safety backup. The common eight-step work horizon is expected to stop first:
starting at level 4, even a maximally deep refinement chain ends at level 12,
leaving at least six NVB levels to the fixed level-18 reference. Eight
refinements give nine solution points and two adjacent rolling three-point
rate comparisons without driving the terminal error into the reference-error
floor. A local kappa=2 level-12 smoke reached a 0.6769 reference-difference/
terminal-error fraction, so the earlier low-wave-number levels were rejected
before server launch rather than being treated as valid candidates.

The tracked templates are:

- `configs/R2a-{palod,afem}-k{2,4,8,16,32}-krobust-server-level18-step8-v4.json`
- `configs/R2a-{palod,afem}-k{2,4,8,16,32}-krobust-server-level18-reference-audit-v4.json`

These are development evidence, not yet a manuscript claim. They may be
promoted only after all trajectory and reference-audit gates below pass and
the manuscript baseline is updated intentionally.

## Resource policy

Run the ten trajectories serially. Start with 16 patch workers and 16 build
jobs. PALOD patch work can use those workers; the current AFEM sparse direct
solve may remain close to one CPU core, which is expected and does not justify
running multiple memory-heavy trajectories concurrently. Use 32 patch workers
only after the kappa=16 PALOD log shows low CPU utilization and ample memory;
do not mix thread counts in an error-versus-time comparison.

The server has about 366 GiB RAM. The launch gate below requires 256 GiB
available before each new trajectory, leaves roughly 110 GiB to the OS and
other users, and runs with swap forbidden by the acceptance criteria. The
Every reference audit reaches level 19; PALOD at level 18 and those audit
solves are the likely memory maxima.

## Build and run all trajectories

From a clean checkout of the commit containing this document:

```bash
git fetch origin
git switch codex/palod-streaming-gram
git pull --ff-only origin codex/palod-streaming-gram

MODE=r2a-krobust PATCH_THREADS=16 JOBS=16 MIN_AVAILABLE_GIB=256 \
RESULT_DIR="$PWD/results/R2a-krobust-k2-32-server" \
  scripts/run_helmholtz_adaptive_paper_server.sh
```

Do not reuse a result directory from another commit. The runner builds once,
injects the actual commit and binary hash into every runtime config, executes
PALOD then AFEM for each wave number, and shares the immutable reference cache.
It writes a `.done` marker only after the run manifest and required artifacts
pass validation.

## Trajectory acceptance

```bash
cd results/R2a-krobust-k2-32-server
find runs -name '*.done' -print
sha256sum -c SHA256SUMS
grep -H -E 'Elapsed|Percent of CPU|Maximum resident|Swaps|Exit status' logs/*.time
grep -H -E 'state=|convergence_regime=|reference_cache=' logs/*.stdout
```

Require ten `.done` files, exit status zero, no swap, and
`driver_state=TrajectoryComplete`. Accept either `fixed H-step trajectory
complete` or `minimum reference/coarse level gap reached`; reject every
`WorkLimitReached`, solver failure, or missing relative reference-energy-error
column. AFEM should normally hit the shared reference cache after the PALOD
run at the same wave number.

## Reference-adequacy audits

Run ten level-plus-one audits after the trajectories complete: one PALOD and
one AFEM source per wave number. The adjacent-reference solution difference is
common, but the denominator is each method's own terminal error; checking both
prevents a smaller AFEM terminal error from inheriting an invalid PALOD ratio.

```bash
cd ~/code/LOD2d-CPP
ROOT="$PWD"
TRAJ="$ROOT/results/R2a-krobust-k2-32-server"
AUDIT="$ROOT/results/R2a-krobust-k2-32-reference-audits"

for k in 2 4 8 16 32; do
  for method in palod afem; do
    source=$(find "$TRAJ/runs/R2a-${method}-k${k}-krobust-server-level18-step8-v4" \
      -name iterations.csv -print -quit)
    SOURCE_ITERATIONS="$source" \
    TEMPLATE="$ROOT/experiments/helmholtz_adaptive_paper/configs/R2a-${method}-k${k}-krobust-server-level18-reference-audit-v4.json" \
    RESULT_DIR="$AUDIT/${method}-k${k}" \
    REFERENCE_CACHE_DIR="$TRAJ/reference-cache" \
    JOBS=16 MIN_AVAILABLE_GIB=256 \
      scripts/run_helmholtz_reference_adequacy_server.sh
  done
done
```

Inspect every `reference_adequacy.json`. The frozen adequacy gate is
`relative_reference_difference / terminal_reference_energy_error <= 0.25`.

If either method fails for a wave number, do not compare that method's terminal
portion. Increase that wave number's reference level, rerun both matching
trajectories from the initial coarse mesh, and repeat both audits. Do not tune
PALOD/AFEM marking or localization from the observed errors.

## Analysis after successful audits

Plot relative reference energy error against coarse DoF, with one PALOD and
one AFEM curve per wave number. Also report `kappa H_max`, PALOD `ell` changes,
and rolling log-log slopes. Compare methods at common error targets or
interpolated common DoF; do not compare raw iteration indices. A k-robust
PALOD conclusion requires overlap of the five PALOD curves over a resolved DoF
range and no systematic k-dependent growth after the pre-asymptotic prefix.
An AFEM non-robust conclusion requires a reproducible growth of its
pre-asymptotic range or DoF-to-target with kappa. If the data do not show those
patterns, report the observed result rather than extending reference meshes
until the expected claim appears.

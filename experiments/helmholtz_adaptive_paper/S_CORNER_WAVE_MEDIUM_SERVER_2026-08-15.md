# Case-S corner-wave four-method medium server run

This is the aligned-level rerun after the PALOD and SLOD performance changes.
It compares PALOD, SLOD, uniform FEM and adaptive FEM for the same manufactured
corner-wave solution (`kappa=16`, singular amplitude 1, smooth wave amplitude
0.05). It is a medium validation trajectory, not yet the final deepest paper
run.

## Frozen work horizons

| Method | Configuration | Horizon | Purpose |
|---|---|---:|---|
| UFEM | `configs/S-corner-wave-ufem-k16-H6-level20-step14-v4.json` | H=6 start, 14 uniform steps, terminal level 20 | uniform low-regularity comparison |
| AFEM | `configs/S-corner-wave-afem-k16-H6-level20-step14-v4.json` | H=6 start, 14 adaptive rounds | manufactured-exact AFEM through the level-20 work horizon |
| SLOD | `configs/S-corner-wave-slod-k16-H6-ell2-gap4-h20-step10-v4.json` | H=6/h=10 start, 10 synchronized H/h refinements, fixed gap 4 | terminal H=16/h=20; memory-safe direct-saddle path |
| PALOD | `configs/S-corner-wave-palod-k16-H6-h12-to-h20-gap4-step10-v4.json` | H=6/h=12 start; expected about 10 H steps | actual local level-gap refresh at `h-H <= 4`, stopping when the reference reaches level 20 |

The server executes the methods in the table order: UFEM, AFEM, SLOD, then
PALOD.  Old medium JSON files are retained only to reproduce the earlier
shorter run and are not selected by `MODE=s-corner-wave-medium`.

All four methods start from NVB level 6.  UFEM and the manufactured-exact AFEM
therefore use 14 refinement rounds to reach the level-20 work horizon.  SLOD
starts from H=6/h=10 and preserves the four-level difference through H=16/h=20.
PALOD starts from H=6/h=12.  It does not infer a level from the H-step count.
After every accepted local H refinement, the driver computes the actual NVB
level difference between every reference child and its parent coarse element.
When the minimum difference reaches four, the ambient shadow is locally
advanced by one additional level and promoted; the inherited coarse mesh and
ell are retained.  Once the promoted reference has reached level 20 and the
gap again reaches four, the trajectory completes.  Ten H steps are expected
for this case; `maximum_H_steps=20` is only a fail-safe ceiling.

SLOD remains on `direct_saddle`; automatic direct-Schur switching is disabled.
The previous level-19 run reached about 152 GiB with late direct Schur, whose
factorization growth makes a level-20 run unsafe on a 366-GB machine.  The
direct-saddle path is slower but has the safer observed memory scaling.  This
level-20 SLOD run is still the dominant memory risk and must be monitored.

## Pull, build and run

Run from a clean server checkout. Replace `<commit>` in the result directory
with the short commit printed by `git rev-parse --short HEAD`.

```bash
cd ~/code/LOD2d-CPP
git fetch origin
git switch codex/palod-streaming-gram
git pull --ff-only origin codex/palod-streaming-gram
git status --short
git rev-parse HEAD

MODE=s-corner-wave-medium \
PATCH_THREADS=16 JOBS=16 MIN_AVAILABLE_GIB=300 \
RESULT_DIR="$PWD/results/S-corner-wave-H6-to-level20-<commit>" \
bash scripts/run_helmholtz_adaptive_paper_server.sh
```

The script runs UFEM, AFEM, SLOD and PALOD sequentially, refuses a dirty tracked checkout,
validates each manifest, keeps a shared reference cache, enables live progress
and writes `/usr/bin/time -v` resource logs. A rerun with the same result
directory skips `.done` cases.

For PALOD, the returned `iterations.csv` and `ell_history.csv` contain the
model and localization-certificate stage profile for every accepted or failed
`ell` attempt. `run.json` aggregates the same columns under
`timing.stage_totals_seconds`; no log parsing is required to recover the
overall time split.

## Monitor without disturbing the run

```bash
RESULT=results/S-corner-wave-H6-to-level20-<commit>
tail -f "$RESULT"/logs/*.stdout
find "$RESULT"/runs -name '*.done' -print
grep -H -E 'LOD2D_(PALOD|SLOD|MODEL)_PROGRESS|LOD2D_MODEL_STAGES|state=' \
  "$RESULT"/logs/*.stdout | tail -n 80
```

`iteration` is the chronological driver-record number. Use `H_step` and
`reference_epoch` for the mathematical adaptive trajectory. For SLOD, inspect
`patch_solver_used` rather than assuming that the requested base solver was
used.

## One feedback bundle

After all four `.done` files appear, run:

```bash
RESULT=results/S-corner-wave-H6-to-level20-<commit>

find "$RESULT"/runs -name '*.done' -print
grep -H -E 'Elapsed|Percent of CPU|Maximum resident|Swaps|Exit status' \
  "$RESULT"/logs/*.time
grep -H -E 'state=|convergence_regime=|reference_cache=' \
  "$RESULT"/logs/*.stdout

python3 - "$RESULT" <<'PY'
import csv, json, pathlib, sys
root = pathlib.Path(sys.argv[1])
for path in sorted(root.glob('runs/*/*/iterations.csv')):
    rows = list(csv.DictReader(path.open()))
    solves = [r for r in rows if r['action'] in {
        'SolvePracticalLod', 'SolveStandardLod',
        'SolveUniformFem', 'SolveAdaptiveFem'}]
    if not solves:
        continue
    last = solves[-1]
    print(json.dumps({
        'run': path.parent.name,
        'method': last['method'],
        'H_step': last['H_step'],
        'epoch': last['reference_epoch'],
        'DoF_H': last['DoF_H'],
        'DoF_ref': last.get('DoF_ref', ''),
        'relative_error': last['relative_exact_energy_error'],
        'patch_solver': last.get('patch_solver_used', ''),
        'auto_schur': last.get('slod_auto_direct_schur', ''),
        'corrector_s': last.get('time_corrector', ''),
        'coarse_factor_s': last.get('time_coarse_factorization', ''),
        'factorization_reuses': last.get('corrector_factorization_reuses', ''),
        'peak_memory_mb': last['peak_memory_mb'],
        'method_seconds': last['time_total_cumulative'],
    }, sort_keys=True))
PY

tar -czf "${RESULT}.tar.gz" \
  "$RESULT"/server-build-identity.txt \
  "$RESULT"/runtime-configs "$RESULT"/logs "$RESULT"/runs \
  "$RESULT"/SHA256SUMS
```

Push the result directory or archive to a result branch and report its commit.
The next decision is based on one bundle: retain or disable the late SLOD
Schur switch, select the deepest safe SLOD horizon, and then generate the
four-method error-versus-DoF and stage-time figures from the same binary.

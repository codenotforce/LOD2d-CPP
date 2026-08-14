# Case-S corner-wave four-method medium server run

This is the single feedback run after the PALOD and SLOD performance changes.
It compares PALOD, SLOD, uniform FEM and adaptive FEM for the same manufactured
corner-wave solution (`kappa=16`, singular amplitude 1, smooth wave amplitude
0.05). It is a medium validation trajectory, not yet the final deepest paper
run.

## Frozen work horizons

| Method | Configuration | Horizon | Purpose |
|---|---|---:|---|
| UFEM | `configs/S-corner-wave-ufem-k16-medium-level20-step15-v4.json` | 15 uniform steps, terminal level 20 | deeper uniform low-regularity comparison |
| AFEM | `configs/S-corner-wave-afem-k16-medium-step28-v4.json` | 28 adaptive steps | adaptive low-regularity comparison |
| SLOD | `configs/S-corner-wave-slod-k16-medium-ell2-gap4-step10-v4.json` | 10 synchronized H/h refinements, fixed level gap 4 | one-slot factorization reuse; late Schur probe |
| PALOD | `configs/S-corner-wave-palod-k16-medium-step15-v4.json` | 15 H steps, refresh after 3/6/9/12 | five epochs, optimized localization path |

The server executes the methods in the table order: UFEM, AFEM, SLOD, then
PALOD.  Old medium JSON files are retained only to reproduce the earlier
shorter run and are not selected by `MODE=s-corner-wave-medium`.

The SLOD base solver is `direct_saddle`. It switches to `direct_schur` only
when the current reference space reaches 200,000 unconstrained DoFs. Every
row records `patch_solver_used` and `slod_auto_direct_schur`, so the server
data will show whether the late switch helped. With fixed `ell=2`, this is an
experimental probe rather than an assumption that Schur must be faster.

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
PATCH_THREADS=16 JOBS=16 MIN_AVAILABLE_GIB=32 \
RESULT_DIR="$PWD/results/S-corner-wave-medium-<commit>" \
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
RESULT=results/S-corner-wave-medium-<commit>
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
RESULT=results/S-corner-wave-medium-<commit>

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

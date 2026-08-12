# Exploratory S/PALOD manufactured-error run (2026-08-13)

This record is exploratory implementation evidence, not paper data.  The
R2a/S reference-adequacy gate was still open when the runs were authorized.
Both runs used the frozen E0 practical parameters, κ=16, two H refinements,
`fixed_work_horizon`, and four patch threads.

## Independent completion checks

The result root is
`results/exploratory-e1-S-exact-errors-2026-08-13`.  Both configuration-level
`.done` files exist.  The two `run.json` files independently report
`status=success`, `driver_state=TrajectoryComplete`, and
`stop_reason="fixed H-step trajectory complete"`; their final evaluated CSV
action is `CompleteTrajectory`.  The runtime build identity is commit
`8f3ef38304f49637a6d7ebb3241c0401dc749900`, binary SHA-256
`18f27d85c41ab72c90427e7e4ddb8933af35a3a41cb8bdf9bb0e03847ae2dc5e`,
and `patch_threads=4`.  A direct `sha256sum -c SHA256SUMS` passed for every
listed payload.

| epoch | reference level | run id | wall time | peak RSS | swaps |
|---:|---:|---|---:|---:|---:|
| 0 | 10 | `S_PALOD_k16_r0_91f566a0af9d84e6` | 19.49 s | 1,351,112 kB | 0 |
| 1 | 11 | `S_PALOD_k16_r0_e5802c171ad52f08` | 40.13 s | 1,825,720 kB | 0 |

`SwapFree` was 6,291,456 kB both before and after the batch.

## Error versus unconstrained coarse DoF

All quantities below are relative errors.  The exact columns compare the
PALOD candidate directly with the singular manufactured solution; the
reference columns compare it with the epoch's FEM evaluation reference.

| epoch | `N_H` | `DoF_H` | reference energy | exact energy | reference L2 | exact L2 |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 225 | 208 | 0.578242 | 0.641097 | 0.381926 | 0.436795 |
| 0 | 252 | 235 | 0.444799 | 0.533577 | 0.232352 | 0.305055 |
| 0 | 288 | 271 | 0.368955 | 0.478111 | 0.179293 | 0.261188 |
| 1 | 225 | 208 | 0.600955 | 0.635721 | 0.397424 | 0.426557 |
| 1 | 253 | 236 | 0.477625 | 0.526802 | 0.252307 | 0.290190 |
| 1 | 291 | 274 | 0.401514 | 0.463293 | 0.193564 | 0.236916 |

Both exact-error trajectories decrease.  At the final epoch-1 point, the
reference energy error is about 13.3% below the manufactured-exact energy
error, whereas the epoch-0 final discrepancy is about 22.8%.  This is direct
evidence that a deeper evaluation reference moves the reported reference error
toward the exact error, but it is not a replacement for the adequacy audit.
The two calibrations restart from the same initial coarse mesh; epoch-start
markers therefore coincide at `DoF_H=208` and do not imply a continuous
cross-epoch adaptive trajectory.

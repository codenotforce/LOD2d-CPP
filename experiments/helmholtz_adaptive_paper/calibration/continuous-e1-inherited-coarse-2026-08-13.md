# Continuous reference-epoch E1 pilot (2026-08-13)

This pilot corrects the semantics of the earlier pair of independent
exploratory runs.  A paper reference refresh promotes the current ambient mesh
to the next reference mesh while retaining the terminal coarse mesh.  It does
not reconstruct the initial hierarchy.

## Code and configuration

- continuous-driver implementation commit:
  `d309038345bf23ba6ec9481a5b3335d61b693f22`
- wrapper execution checkout: `6c9ae2f`
- configuration:
  `configs/S-palod-k16-continuous-e1-epoch0-to1-step4-v4.json`
- case/method: S/PALOD, kappa=16
- frozen E0 practical parameters; four patch threads
- fixed horizon: four cumulative H refinements
- scheduled refresh: after cumulative H-step 2
- audited wrapper run id: `S_PALOD_k16_r0_1c54bb4141b7fe83`

The driver records `CompleteReferenceEpoch` before promotion and
`RefreshReferenceEpoch` after promotion.  It rejects a refresh if coarse node,
DoF, element, or mesh-version identity changes.  The evaluation reference is
rebuilt independently for each epoch and remains one-way post-processing.

## Inheritance evidence

| boundary record | epoch | `N_H` | `DoF_H` | `N_ref` | `N_amb` | ell |
|---|---:|---:|---:|---:|---:|---:|
| terminal solve | 0 | 288 | 271 | 3201 | 4209 | 2 |
| refresh complete | 1 | 288 | 271 | 4209 | 4209 | 2 |
| first epoch-1 solve | 1 | 288 | 271 | 4209 | 4209 | 2 |

Thus epoch 1 starts on the exact terminal coarse grid of epoch 0, while its
reference mesh is the previous ambient mesh.

## Evaluated trajectory

| epoch | `N_H` | `DoF_H` | reference energy | relative exact energy | reference L2 | relative exact L2 |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 225 | 208 | 0.578242 | 0.641097 | 0.381926 | 0.436795 |
| 0 | 252 | 235 | 0.444799 | 0.533577 | 0.232352 | 0.305055 |
| 0 | 288 | 271 | 0.368955 | 0.478111 | 0.179293 | 0.261188 |
| 1 | 288 | 271 | 0.413134 | 0.456329 | 0.195247 | 0.228681 |
| 1 | 351 | 334 | 0.298751 | 0.358690 | 0.133022 | 0.168676 |
| 1 | 442 | 425 | 0.212844 | 0.292275 | 0.085746 | 0.125864 |

The fail-closed wrapper reports `status=success`,
`driver_state=TrajectoryComplete`, and the required fixed-horizon stop reason.
Its final CSV action is `CompleteTrajectory`; the `.done` sentinel exists only
after those checks.  Method time was about 63.00 s, `/usr/bin/time` wall time
64.08 s, process peak RSS 4,142,568 kB, and swaps zero.  Independent
`sha256sum -c SHA256SUMS` verification passed for every listed file, including
`run.json`, all CSV/VTK outputs, `.stdout`, `.time`, `.done`, runtime config,
reference caches, and build identity.  The audited root is
`results/continuous-e1-S-PALOD-k16-wrapper-2026-08-13`.

Reference adequacy is not enabled, so this remains a controlled E1 pilot rather
than paper data.

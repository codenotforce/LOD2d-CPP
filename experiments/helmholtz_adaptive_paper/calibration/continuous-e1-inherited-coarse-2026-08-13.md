# Continuous reference-epoch E1 pilot (2026-08-13)

This pilot corrects the semantics of the earlier pair of independent
exploratory runs.  A paper reference refresh promotes the current ambient mesh
to the next reference mesh while retaining the terminal coarse mesh.  It does
not reconstruct the initial hierarchy.

## Code and configuration

- implementation baseline: `d309038345bf23ba6ec9481a5b3335d61b693f22`
- configuration:
  `configs/S-palod-k16-continuous-e1-epoch0-to1-step4-v4.json`
- case/method: S/PALOD, kappa=16
- frozen E0 practical parameters; four patch threads
- fixed horizon: four cumulative H refinements
- scheduled refresh: after cumulative H-step 2
- run id: `S_PALOD_k16_r0_039f918ff7b6fdda`

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

The direct smoke reports `status=success`,
`driver_state=TrajectoryComplete`, and the required fixed-horizon stop reason.
Method time was about 59.22 s and process peak memory about 4046 MiB.  Because
this direct invocation did not produce the server wrapper's `.time`, `.done`,
and checksum artifacts, those resource/completion checks must be repeated with
the wrapper before treating this as a frozen result.  Reference adequacy is
also not enabled, so this remains a controlled E1 pilot rather than paper data.

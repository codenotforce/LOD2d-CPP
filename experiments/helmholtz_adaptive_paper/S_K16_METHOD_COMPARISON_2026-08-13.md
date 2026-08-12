# Case S, kappa=16 method comparison (2026-08-13)

## Scope and interpretation

This record freezes the manufactured-solution comparison requested for case S
at `kappa=16`. PALOD, SLOD, uniform P1 FEM (`UFEM`), and adaptive P1 FEM
(`AFEM`) all start from coarse level 5 with 104 unconstrained coarse degrees of
freedom. They use the same case registry, impedance coefficient, manuscript
baseline, and quadrature policy.

The plotted quantity is the relative weighted-energy error against the analytic
manufactured solution. It is post-processing data and never enters MARK or
STOP. Comparator curves are displayed only through their first evaluated point
at or below the terminal PALOD error. The full fixed-horizon run artifacts are
retained. The manifests label these runs `implementation-study`; promoting the
figure to a formal paper result remains a manuscript-level decision.

## Frozen trajectories

- PALOD: `S-palod-k16-comparison-3epoch-step8-v4.json`; three continuous
  reference epochs, with each new epoch inheriting the preceding terminal
  coarse mesh.
- SLOD: `S-slod-k16-comparison-ell2-fixed-ratio-step7-v4.json`; `ell=2`, and
  both H and h are uniformly refined by one level at every step, preserving the
  initial five-level fine/coarse gap. The completed run used `direct_saddle`.
- UFEM: `S-ufem-k16-comparison-level13-step8-v4.json`; uniform conforming P1
  refinements.
- AFEM: `S-afem-k16-comparison-level15-step14-v4.json`; residual Dörfler
  marking, reference level 15, and a 14-step fixed work horizon. The exact
  target does not stop the online trajectory.

## Accepted runs

| Method | Run ID | Result root | Acceptance |
|---|---|---|---|
| PALOD | `S_PALOD_k16_r0_3f07823497504606` | `results/S-k16-four-method-comparison-2026-08-13` | `success/TrajectoryComplete`, `.done`, SHA256 clean |
| SLOD | `S_SLOD_k16_r0_47640db67a48178a` | `results/S-k16-SLOD-ell2-fixed-ratio-comparison-2026-08-13` | `success/TrajectoryComplete`, `.done`, SHA256 clean, exit 0, swap 0 |
| UFEM | `S_UFEM_k16_r0_e77c9f13f580747f` | `results/S-k16-four-method-comparison-2026-08-13` | `success/TrajectoryComplete`, `.done`, SHA256 clean |
| AFEM | `S_AFEM_k16_r0_060a8c5e15a6cace` | `results/S-k16-AFEM-level15-comparison-2026-08-13` | `success/TrajectoryComplete`, `.done`, SHA256 clean, exit 0, swap 0 |

SLOD used 6809.65 s of method time and 1:54:42 wall time, with a peak RSS of
about 5.39 GiB. AFEM used 2463.52 s of method time and 44:49 wall time, with a
peak RSS of about 454 MiB and a reference-cache hit. For both methods the
dominant cost is mesh/refinement hierarchy work rather than the final linear
solve. Raw result directories and generated figures remain untracked.

## PALOD target and first crossings

The frozen comparison target is the terminal PALOD exact relative
weighted-energy error

`0.13779351073909829`.

| Method | Epoch/index | Iteration | N_H | DoF_H | Exact relative weighted-energy error | DoF / PALOD DoF |
|---|---:|---:|---:|---:|---:|---:|
| PALOD | 2 | 36 | 808 | 790 | 0.13779351073909829 | 1.000 |
| SLOD | 6 | 13 | 6273 | 6208 | 0.13560788885142366 | 7.858 |
| UFEM | 0 | 17 | 24833 | 24704 | 0.12131756410742389 | 31.271 |
| AFEM | 0 | 27 | 2080 | 2053 | 0.12524231704425701 | 2.599 |

For SLOD the reported epoch field is the synchronized fine-grid evaluation
index, not a PALOD reference-refresh decision.

The three PALOD epoch-start markers are:

| PALOD epoch | Iteration | N_H | DoF_H | Exact relative weighted-energy error |
|---:|---:|---:|---:|---:|
| 0 | 4 | 113 | 104 | 0.85259625657919824 |
| 1 | 17 | 206 | 195 | 0.46218741735638830 |
| 2 | 30 | 460 | 447 | 0.20821022694250013 |

The plotted `N^{-1/2}` segment is a reference slope, not a fitted convergence
claim. The generated CSV and JSON beside the figure preserve all plotted points
and the first-crossing records.

## Reproduction

```bash
python3 tools/visualization/plot_helmholtz_method_comparison.py \
  --palod=results/S-k16-four-method-comparison-2026-08-13/runs/S-palod-k16-comparison-3epoch-step8-v4/S_PALOD_k16_r0_3f07823497504606 \
  --slod=results/S-k16-SLOD-ell2-fixed-ratio-comparison-2026-08-13/runs/S-slod-k16-comparison-ell2-fixed-ratio-step7-v4/S_SLOD_k16_r0_47640db67a48178a \
  --ufem=results/S-k16-four-method-comparison-2026-08-13/runs/S-ufem-k16-comparison-level13-step8-v4/S_UFEM_k16_r0_e77c9f13f580747f \
  --afem=results/S-k16-AFEM-level15-comparison-2026-08-13/runs/S-afem-k16-comparison-level15-step14-v4/S_AFEM_k16_r0_060a8c5e15a6cace \
  --output=results/S-k16-method-comparison-2026-08-13/S-k16-error-vs-DoF.png
```

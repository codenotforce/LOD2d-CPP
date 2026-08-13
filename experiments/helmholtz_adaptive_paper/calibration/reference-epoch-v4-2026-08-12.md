# Schema-v4 reference-epoch gate record (2026-08-12)

This is a calibration and adequacy decision record, not paper data. The
underlying result directories remain outside the commit; each directory has a
portable, relative-path `SHA256SUMS` manifest. Reference errors from different
epochs are not joined into one convergence curve.

## Uploaded S epoch-2 calibration

The uploaded result bundle
`results/S-palod-k16-level12-calibration-v4-server` passed all 14 payload
hashes after remapping the original server-root prefix. Its run
`S_PALOD_k16_r0_6225f02ebe6c1f0b` has `status=success`,
`driver_state=TrajectoryComplete`, the exact fixed-horizon stop reason,
`reference_epoch=2`, and `reference_level=12`. `/usr/bin/time` recorded zero
process swaps and 17.31 GiB peak RSS. The host's pre-existing swap occupancy
was unchanged before and after the run.

The uploaded `server-build-identity.txt` records `patch_threads=8`, not the
frozen value 4. The numerical trajectory is retained as calibration evidence,
but it is not evidence for the four-thread resource protocol. New
calibration/custom runs are guarded against non-frozen patch thread counts.

| `N_H` | reference energy error | `U_prac` | `eta_H` | `Theta_loc` | `ell` |
|---:|---:|---:|---:|---:|---:|
| 225 | 0.6156762099 | 4.794521138 | 9.672941218 | 0.0723858540 | 2 |
| 253 | 0.4969073272 | 4.366342336 | 8.815473639 | 0.0663620617 | 2 |
| 291 | 0.4250116577 | 3.808531101 | 7.689490921 | 0.0661282540 | 2 |
| 354 | 0.3195878615 | 2.831885924 | 5.697764504 | 0.0951096507 | 2 |
| 443 | 0.2385495181 | 2.134613599 | 4.276311022 | 0.1312769189 | 2 |
| 564 | 0.1490872025 | 1.358488539 | 2.711117254 | 0.1633192020 | 2 |
| 720 | 0.1012683802 | 0.910037259 | 1.813457845 | 0.1758010174 | 2 |

The last two ratios are 0.624973816 and 0.679256022; their geometric mean is
0.651549866 and the last-three-point relative variation is 1.355617002.
Therefore the trajectory is in stable decay, not a plateau. The common target
0.1 is not reached (the terminal error is about 1.27% above it).

The final mesh has 720 points and 1369 triangles. About 99.16% of
`sum(eta_H_T^2)` lies in `r<0.5`, whereas the outer `r>0.75` region remains at
the coarse diameter 0.1767767. The strongest refinement is in the cut-off
transition annulus `0.25<r<0.5`; all six elements incident to the re-entrant
corner remain at diameter 0.125. This supports localization near the singular
manufactured-data support, but it does not yet pass the G7 corner-recovery
claim.

## Independent reference audits

Each row solves only two adjacent reference FEM problems and consumes the
listed completed PALOD trajectory as immutable input.

| case/epoch | reference levels | terminal error | adjacent-reference difference | terminal fraction | gate |
|---|---:|---:|---:|---:|---|
| R2a/2 | 12 to 13 | 0.0050982358 | 0.0434043966 | 8.513611 | fail |
| S/2 | 12 to 13 | 0.1012683802 | 0.1166080460 | 1.151475 | fail |
| R2a/3 | 13 to 14 | 0.0062762928 | 0.0289718425 | 4.616076 | fail |
| R2a/4 | 14 to 15 | 0.0072442005 | 0.0212730523 | 2.936563 | fail |
| R2a/5 | 15 to 16 | 0.0058738093 | 0.0143569347 | 2.444229 | fail |

All fractions exceed the frozen maximum 0.25, so neither case may enter E1.
R2a epoch-2 was reconstructed at code baseline `1061268` because the original
local result directory was absent; its terminal value exactly reproduces the
recorded value. The reconstruction used four patch threads, 1.25 GiB peak RSS,
and zero swaps.

## R2a deeper local calibrations

| epoch/level | run id | terminal error | window geometric ratio | regime | method seconds | peak MiB | swaps |
|---|---|---:|---:|---|---:|---:|---:|
| 3/13 | `R2a_PALOD_k16_r0_bb9f4efb02525125` | 0.0062762928 | 0.675498 | stable decay | 116.13 | 1555.88 | 0 |
| 4/14 | `R2a_PALOD_k16_r0_8408e4a4538484d3` | 0.0072442005 | 0.697018 | pre-asymptotic | 352.28 | 3248.39 | 0 |
| 5/15 | `R2a_PALOD_k16_r0_6ab2fd848cf533f4` | 0.0058738093 | 0.620047 | pre-asymptotic | 1472.40 | 8400.91 | 0 |

All three runs completed the same six-H-step horizon with the frozen E0
parameters and four patch threads. Level 15 took 24:49 wall time and 8.40 GiB
peak RSS. A level-16 calibration is therefore moved to the server instead of
risking the 23 GiB WSL memory boundary.

## Next gate

Run, sequentially on the server with `PATCH_THREADS=4`:

1. R2a epoch 6 / reference level 16, six H-steps;
2. S epoch 3 / reference level 13, six H-steps;
3. R2a level 16 to 17 and S level 13 to 14 adequacy audits.

Formal E1 configuration generation and all ten E1 trajectories remain
blocked until both final terminal fractions are at most 0.25.

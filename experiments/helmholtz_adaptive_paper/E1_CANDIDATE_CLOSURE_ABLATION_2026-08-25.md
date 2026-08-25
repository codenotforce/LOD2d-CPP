# E1 candidate closure/batching ablation (2026-08-25)

This gate evaluates two implementation-level changes before another large E1
run.  The H-step budget, `theta_H=0.1`, `theta_c=0.3`, H6/h12 start, reference
refresh policy and all estimators are otherwise unchanged.

- **closure-aware marking** keeps the original global Dörfler mass condition,
  but ranks the leading `2 m` standard candidates by indicator mass divided by
  an independently estimated NVB closure cost;
- **candidate batching** performs the RT2 accuracy sweep every second H-step,
  while still forcing a sweep before a structural/level-gap refresh, numerical
  dual check, tolerance decision or termination.  Intermediate steps perform
  only hierarchy containment closure.

The closure cost is a deterministic independent-seed proxy, not the exact
marginal cost after all preceding marks.  The final selected set is nevertheless
checked against the original, unmodified Dörfler mass.

## Local ten-step gate

All rows used the same current executable and manuscript baseline.  Times are
the runner's cumulative timers; peak memory is the runner's RSS sample.

| policy | RT2 sweeps | final candidate DoF | final reference DoF | direct marks | candidate added elements | closure-added elements | candidate time (s) | method time (s) | total time (s) | peak MiB | final exact energy error |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| baseline | 9 | 9069 | 8271 | 2471 | 3930 | 1459 | 7.665 | 23.603 | 24.212 | 3903.0 | 0.07586447 |
| stride 2 | 7 | 8591 | 8264 | 1647 | 2466 | 819 | 5.604 | 21.589 | 22.252 | 3863.7 | 0.07586619 |
| closure aware | 9 | 8970 | 8255 | 2437 | 3584 | 1147 | 7.023 | 22.735 | 23.364 | 4189.6 | 0.07588170 |
| stride 2 + closure aware | 7 | 8559 | 8255 | 1647 | 2288 | 641 | 5.133 | 20.767 | 21.380 | 4123.8 | 0.07588170 |

Relative to the baseline, the combined gate reduced candidate time by 33.0%,
method time by 12.0%, candidate closure additions by 56.1%, and final candidate
DoF by 5.6%.  The final exact error changed by about 0.023%.  Peak RSS increased
by about 5.7%, so production memory gates must remain enabled.

The candidate containment pass added no cells in these ten-step trajectories.
Consequently, changing the full-sweep order to containment-first would add
complexity without a measured E1 benefit: after a nontrivial containment change
the candidate solution and RT2 indicator would have to be recomputed on the new
mesh.  The batching path already performs containment-only updates safely.

`candidate_force_level_gap=4` must not be used with the present global minimum
gap: the local gate showed that it forces an RT2 sweep at almost every H-step and
neutralizes batching.  The tested configuration uses stride 2 and the existing
refresh/dual/termination force conditions.  A value of 3 merely makes the
pre-refresh intent explicit for the current trigger and does not change the
refresh decision.

## Remaining bottleneck and scope

For the combined gate, the largest phases were:

| phase | time (s) |
|---|---:|
| localization certificate (`theta`) | 9.345 |
| Gram patch solves | 4.753 |
| Gram factorization | 2.484 |
| candidate RT2 reconstruction | 4.740 |
| LOD build | 3.043 |
| reference Riesz | 1.764 |

The existing implementation already prepares Gram patch structures in
parallel, factors patch-local systems in parallel, applies block Gram actions,
prolongs the Ritz subspace across H-steps, and uses bounded cross-step factor
caches.  The ten-step gate needed only 19 Ritz iterations; enlarging thread
counts or the low-hit factor LRU is therefore not the next safe optimization.

Promotion-time candidate rebuild and active-region gap reserve are deliberately
not enabled by this gate.  They change which mesh becomes the next reference
and require new promotion-time indicator, nestedness, stability, localization
and (for numerical refresh) dual checks.  They need a separate algorithmic
ablation before they can replace the current paper trajectory.  Until that
ablation exists, the validated combined policy is an optional performance
configuration rather than a silent change to the frozen main run.

## Reproduction configs

- `E1-R1-palod-k16-H6-h12-gap6-schur-thetaH01-thetaC03-stride1-step10-ablation-v6.json`
- `E1-R1-palod-k16-H6-h12-gap6-schur-thetaH01-thetaC03-stride2-step10-pilot-v6.json`
- `E1-R1-palod-k16-H6-h12-gap6-schur-thetaH01-thetaC03-closurecost-step10-ablation-v6.json`
- `E1-R1-palod-k16-H6-h12-gap6-schur-thetaH01-thetaC03-stride2-closurecost-step10-ablation-v6.json`

No large experiment should be launched from these ten-step files.  The next
gate is a bounded medium trajectory with server memory monitoring; only after
its epoch-wise error and reference-floor audit passes should the combined
policy be copied into a production configuration.

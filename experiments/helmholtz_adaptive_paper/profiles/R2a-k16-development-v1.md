# R2a / kappa=16 development profile v1

This profile is an implementation diagnostic, not a paper result.  It uses
`R2a-palod-k16-development-profile-v2.json`, the E0 diagnostic admissibility
value `c_H=2.8567113959936523`, and `maximum_H_steps=0`.  The production
admissibility constant remains to be frozen separately.

Environment: WSL Ubuntu 22.04, GCC 11.4 Release, 16 hardware threads.  Both
runs have ID `R2a_PALOD_k16_r0_dfbd38712f97653a` and use the same numerical
configuration and fixed reference epoch.

| quantity (seconds) | repeated factorization | batched patch RHS | ratio |
|---|---:|---:|---:|
| method total | 10.216926386 | 1.509227582 | 6.77x |
| corrector rebuild | 0.033001048 | 0.030629280 | 1.08x |
| localization certificate | 10.106927037 | 1.404376975 | 7.20x |
| coarse solve | 0.000082328 | 0.000079313 | 1.04x |
| reference-kernel estimator | 0.032844600 | 0.029046000 | 1.13x |
| peak memory (MiB) | 84.91796875 | 87.36328125 | 0.97x |

The optimized ambient-defect Riesz path factorizes each patch saddle matrix
once and solves all coarse-input right-hand sides together.  Structural tests
check one factorization per patch.  The before/after runs agree on
`Theta_loc=0.18169031406509437`, the 23-element Doerfler marking, reference
errors, and the practical indicator to floating-point roundoff.

Conclusion: localization certification, not corrector construction, was the
measured bottleneck on this representative reduced trajectory.  The next
performance experiment should freeze a feasible production `c_H`, permit one
real H refinement, and only then evaluate corrector/patch caching against the
full-rebuild path.

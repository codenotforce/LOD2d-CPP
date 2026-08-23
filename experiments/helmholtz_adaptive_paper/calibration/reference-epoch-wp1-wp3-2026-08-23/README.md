# Reference-epoch WP1--WP5 calibration (2026-08-23)

This is an implementation-study record, not a rigorous certified result and
not a paper production trajectory.

## Frozen inputs

- manuscript:
  `../LOD_paper/helmholtz_lod_certified_amsart_revised.tex`
- manuscript SHA-256:
  `94b0c1469312ce006f3b76d08b30f920115d274f442e3912ca660ccf919bd3f9`
- code base: worktree based on `9aeac45`
- build: WSL Ubuntu 22.04, GCC 11.4, Release, OpenMP enabled
- algebraic cases: localized smooth R1 and mixed-boundary S, `kappa=16`

## WP1 hierarchy transaction

`test_helmholtz_reference_epoch_hierarchy` checks that:

1. candidate is initialized from the reference mesh;
2. proposing a coarse refinement does not mutate the committed coarse mesh;
3. candidate closure contains the proposed coarse mesh;
4. candidate enrichment leaves the reference geometry and version unchanged;
5. a proposal beyond the reference returns `ReferenceRefreshRequired`;
6. promotion of candidate to reference advances the epoch and then permits the
   pending coarse commit.

## WP2 reference corrector sweep

| ell | Theta_loc | direct defect | theorem lower | theorem upper |
|---:|---:|---:|---:|---:|
| 1 | 2.50414 | 4.55780 | 0.0391272 | 160.265 |
| 2 | 1.41851 | 3.09622 | 0.0221642 | 90.7844 |
| 3 | 4.61105e-08 | 0 | 7.20476e-10 | 2.95107e-06 |

The theorem bracket uses deliberately conservative offline test constants
`C_a=C_ov=C_sd=8` and `c_W=0.125`; it validates the implementation direction
but does not freeze production constants.  `G_loc` is applied matrix-free;
the small calibration dimension additionally receives a dense cross-check.

## WP3 candidate RT2/P2 residual audit

| case | eta_eq | residual dual | compatibility | P2 divergence | normal jump | boundary flux | patches | marked |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| R1 | 9.51245 | 5.99362 | 0 | 1.13909e-13 | 2.75559e-14 | 1.17175e-14 | 9 | 4 |
| S | 8.91630 | 0.806316 | 0 | 2.84226e-14 | 6.21783e-15 | 1.00418e-15 | 21 | 3 |

Both global reconstructions satisfy the compatibility, element-divergence,
and Robin boundary-flux identities to the linear algebra tolerance.  In both
cases the global indicator controls the independently solved discrete energy
Riesz norm.  The RT2 basis uses a contravariant Piola map from the reference
triangle.  The old RT0/P0 API remains an explicitly non-production smoke path.
The active-region mode is separately checked and is labeled as a marking
heuristic without a global reliability claim.

## WP4 candidate dual-Riesz gap audit

The candidate kernel path produced

```text
eta_dual_c=0.961709
L_c=0.0150267
candidate_error=1.13069
L_gap_c=0.0150267
reference_candidate_gap=1.13069
```

The direct candidate Galerkin solve is isolated in the small-grid test.  The
production entry point only receives the already available LOD field on the
candidate mesh.  Both theorem inequalities hold, and practical output is
named `L_gap_practical_c` instead of being promoted to `L_gap_c`.

## WP5 reference-epoch state machine

`test_helmholtz_reference_epoch_driver` verifies:

1. corrector failure returns only to `CorrectorCheck` after global `ell++`;
2. lazy dual checks may be skipped, while termination forces one check;
3. a structural hierarchy trigger refreshes reference before coarse commit;
4. a numerical dual gap refreshes instead of terminating;
5. resource and `ell_max` limits exit as `WorkLimitReached`;
6. the backend interface contains no candidate Helmholtz solve operation.

## Reproduction

```bash
cmake -S . -B build-wp3-wp5 -DCMAKE_BUILD_TYPE=Release -DLOD2D_USE_OPENMP=ON
cmake --build build-wp3-wp5 -j 8 --target \
  test_helmholtz_reference_epoch_hierarchy \
  test_helmholtz_reference_corrector_certificate \
  test_helmholtz_candidate_flux \
  test_helmholtz_candidate_dual \
  test_helmholtz_reference_epoch_driver
ctest --test-dir build-wp3-wp5 --output-on-failure \
  -R 'helmholtz_(reference_epoch_hierarchy|reference_corrector_certificate|candidate_flux|candidate_dual|reference_epoch_driver)$'
```

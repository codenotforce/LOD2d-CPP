# Case-S kappa=16 PALOD/SLOD Bottleneck Study

This note records a local performance study. The configurations and outputs
are development evidence only: they are not calibration data, reference-gate
evidence, or paper results.

## Scope and machine

- case: `S`, manufactured solution, `kappa=16`;
- coarse mesh: level 5;
- evaluation reference: level 10, served from the normal reference cache;
- CPU visible to WSL: Intel Core i9-14900HX, 16 cores / 32 hardware threads;
- SLOD: P1, `ell=2`, five synchronized H/h refinements, so both H and h are
  refined once per trajectory step;
- PALOD: three H refinements with its normal practical localization logic;
- every reported run ended with `status=success` and
  `driver_state=TrajectoryComplete` and used no swap.

The short SLOD case ends at 98,817 reference nodes and takes about 35--61 s
wall time, which is long enough to expose scaling while remaining suitable for
interactive profiling. All thread-count runs generated identical DoF/error
trajectories. Direct-saddle and direct-Schur errors agreed to
`4.44e-14` or better.

## Historical full-comparison baseline

The earlier seven-step SLOD comparison spent 6,644.9 s of its 6,809.7 s
method time in mesh/embedding work (97.6%). That measurement predates the NVB
lineage-based incremental embedding change in commit `83fe84d`; it describes
the removed global geometry-search bottleneck and must not be used to predict
the current SLOD cost.

The earlier three-epoch PALOD comparison is still representative of its cost
split: 234.9 s of 277.0 s method time (84.8%) was certificate work, versus
7.9 s (2.8%) for correctors.

## SLOD thread scaling after the NVB embedding optimization

The table uses the direct-Schur backend and excludes post-processing reference
evaluation from method time.

| Patch threads | Method (s) | Mesh update (s) | Correctors (s) | Wall (s) | Peak RSS (MiB) | Method speedup |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 46.587 | 4.832 | 30.218 | 60.89 | 1,180 | 1.00x |
| 2 | 32.733 | 4.887 | 16.251 | 44.53 | 1,196 | 1.42x |
| 4 | 26.110 | 5.035 | 9.218 | 37.80 | 1,209 | 1.78x |
| 8 | 23.511 | 5.014 | 6.404 | 35.02 | 1,222 | 1.98x |
| 16 | 23.140 | 5.154 | 5.220 | 32.65 | 1,297 | 2.01x |

Corrector construction scales by about 5.8x from 1 to 16 threads, but whole
method time scales by only 2.0x. From 8 to 16 threads the method improves by
only 1.6%, while peak memory increases by about 75 MiB. Eight threads was the
best local setting in this measurement; server runs may choose another value
and must record it in the build identity.

At four threads, optional stage probes give the following cumulative SLOD
model-build split over the six evaluated meshes:

| Stage | Time (s) | Share of profiled model build |
|---|---:|---:|
| coarse/reference embedding and quasi-interpolation | 8.923 | 42.0% |
| local correctors | 9.553 | 45.0% |
| multiscale basis and factorization | 2.335 | 11.0% |
| operator assembly | 0.424 | 2.0% |

There is another 4.98 s in synchronized hierarchy refinement outside the
model build. Consequently, after the lineage optimization the remaining SLOD
bottleneck is shared between local correctors and rebuilding the
coarse-to-reference embedding/interpolation inside every model construction.
The latter still calls `build_nested_mesh_embedding(coarse_mesh, fine_mesh)`
and repeats a full geometric reconstruction even though the hierarchy already
owns exact NVB prolongation matrices.

For this P1, `ell=2` case, direct-saddle is marginally faster than
direct-Schur at four threads:

| Solver | Method (s) | Correctors (s) | Peak RSS (MiB) |
|---|---:|---:|---:|
| direct-saddle | 25.993 | 8.847 | 1,217 |
| direct-Schur | 26.307 | 9.553 | 1,209 |

The difference is small, but it does not support making direct-Schur the
default for this P1 SLOD trajectory. Direct-Schur remains useful for larger hp
constraint systems where its reduced system offsets the Schur construction
cost.

## SLOD asymptotic work model and 2026-08-15 optimization

Let `M` be the number of coarse patches, `p` the largest local patch DoF,
`c` the number of patch constraints and `r` the number of local right-hand
sides. Ignoring sparsity-dependent constants, the current stages have the
following leading costs:

| Stage | Work model | Parallel status |
|---|---|---|
| synchronized H/h refinement and lineage update | `O(N_ref)` | serial hierarchy mutation |
| local patch assembly | `O(M p^2)` upper model | parallel over patches |
| direct-saddle factor/solve | `O(M (p+c)^3 + M (p+c)^2 r)` | parallel over patches |
| direct-Schur factor/solve | `O(M[p^3+p^2(c+r)+c^3])` | parallel over patches |
| basis packing and sparse assembly | `O(M p r)` | parallel patch packing, serial deterministic merge |
| coarse Galerkin operator | sparse products depending on basis overlap | Eigen sparse kernels |
| 2-D coarse direct factorization | approximately `O(N_H^(3/2))` work and `O(N_H log N_H)` memory | largely library-controlled |

For this synchronized P1 SLOD trajectory, fixed `ell=2` and a fixed fine/coarse
ratio keep `p`, `c` and `r` bounded while `M` grows. The measured maximum was
`p=728`, `c=34`. Therefore the total patch work grows mainly through more
independent patches; increasing global reference DoFs alone does not create a
larger Schur complement advantage. An automatic solver switch is nevertheless
available through `slod_direct_schur_min_reference_dofs`, but it is a measured
policy input and is never silently applied.

The corrector cache now has an explicit bounded symbolic/factorization reuse
policy. `patch_symbolic_cache_slots=1` keeps one exact local factorization per
worker and `patch_reuse_identical_factorization=true` reuses it only after the
existing matrix-value identity check succeeds. A five-step, 16-thread local
comparison gave:

| Policy | Method (s) | Model build (s) | Correctors (s) | Factorization reuses | Final solver | Final relative energy error |
|---|---:|---:|---:|---:|---|---:|
| direct-saddle, no reuse | 17.420 | 10.950 | 7.284 | 0 | direct-saddle | 0.20786163464178672 |
| direct-saddle, one-slot reuse | 16.538 | 9.774 | 6.117 | 2,659 | direct-saddle | 0.20786163464178672 |
| auto Schur at `DoF_ref>=24000`, one-slot reuse | 18.027 | 11.539 | 8.072 | 3,011 | direct-Schur | 0.20786163464178639 |

Thus bounded reuse improves this trajectory by about 5.1% without changing
the numerical result, whereas the early Schur switch is about 9.0% slower than
the optimized saddle path. The medium server experiment delays the Schur
probe until `DoF_ref>=200000` and records the actual solver at every point; its
purpose is to decide from a genuinely larger server point whether that policy
should remain in the paper configuration.

`iterations.csv` and `LOD2D_MODEL_STAGES` now expose mesh interpolation,
operator assembly, corrector wall time, summed patch assembly/solve/pack work,
basis assembly, coarse-operator assembly and coarse factorization separately.
They also record patch dimensions, worker count, symbolic analyses/reuses,
factorization reuses, requested/actual patch solver and whether the automatic
Schur threshold fired.

## PALOD thread scaling and certificate profile

| Patch threads | Method (s) | Correctors (s) | Certificate (s) | Wall (s) | Peak RSS (MiB) |
|---:|---:|---:|---:|---:|---:|
| 1 | 19.441 | 4.034 | 14.096 | 22.33 | 734 |
| 4 | 15.842 | 1.471 | 13.034 | 16.17 | 774 |
| 8 | 16.269 | 1.104 | 13.797 | 16.60 | 814 |

Four threads is the best tested PALOD setting. Increasing to eight threads
only accelerates the already-small corrector component and cannot compensate
for the serial certificate cost.

The four-thread localization probe accumulated the following times over five
certificate evaluations (the first mesh required an `ell=2` attempt followed
by `ell=3`):

| Certificate stage | Time (s) | Share of certificate |
|---|---:|---:|
| ambient defect Riesz patch solves | 12.921 | 96.5% |
| reference retraction | 0.274 | 2.0% |
| defect right-hand side | 0.159 | 1.2% |
| generalized spectrum | 0.032 | 0.2% |
| coarse energy | 0.002 | <0.1% |

The final generalized eigensolve used 694 iterations, but its measured time
was only 0.022 s; it is not the bottleneck. `compute_ambient_defect_riesz`
currently loops over ambient patches serially, assembles each dense local
energy/RHS block, solves the local Riesz columns, and accumulates a dense Gram
matrix. This stage is about 79% of the complete profiled PALOD method time.

## Implemented optimizations and measured result

The hierarchy-owned nodal/element/DG prolongations and quasi-interpolation are
now consumed directly by SLOD model construction.  The retained algebraic
dimension, element-level, DG-consistency and right-inverse checks replace the
redundant global geometry search.  At four threads, cumulative in-model
mesh/interpolation time fell from 9.194 s to 0.227 s (97.5%), and method time
fell from 25.993 s to 17.137 s (34.1%).  The eight-thread method time is
14.711 s.  The exact-error trajectory is unchanged.

The ambient-defect Riesz patches are now solved in parallel.  Each worker
writes to a fixed patch-index slot and the Gram matrix and column norms are
reduced afterwards in the original patch order.  One- and four-thread test
results are bitwise identical.  The bounded PALOD method times changed from
19.441/15.842/16.269 s at 1/4/8 threads before the optimization to
20.113/7.610/6.697 s afterwards; 16 threads takes 5.207 s.  At four threads,
certificate time fell from 13.034 s to 4.712 s.  The small one-thread
difference is timing noise/parallel bookkeeping and not a numerical change.

The old server choice of four patch threads was made while certificate work
was serial.  Error-versus-DoF data are thread-count independent.  A published
time comparison must record and use a common thread count and binary for all
methods, but the project does not freeze one global thread count in advance.

The follow-up streaming implementation keeps the defect RHS sparse, solves
only the columns active on each patch, reduces compact local Gram matrices in
the original patch order in bounded batches, and discards local vectors after
each batch.  Full patch/local-solution data are retained only by an explicit
diagnostic mode.  On the same local four-thread three-step configuration,
peak RSS fell from about 768 MiB to 211 MiB and wall time from 11.29 s to
10.31 s.  At the final bounded certificate, 27,787 active RHS columns were
solved instead of the dense 40,170 patch-column combinations.  Thread count
is no longer treated as a frozen DoF/error parameter; it must merely be
recorded and held common within any published timing comparison.

### Deep-PALOD spectrum and reduction optimization (2026-08-14)

The deep trajectory exposed two costs that the earlier small probe could not
see.  For a coarse dimension `n`, explicit Cholesky whitening formed dense
`n x n` inverse/whitened matrices and incurred cubic work before the power
iteration.  In addition, the compact patch Gram blocks were scattered into
the global Gram matrix serially.  The production path now:

- retains the sparse coarse energy operator and, for `n >= 1025`, applies a
  sparse-factorized generalized power iteration directly in coarse
  coordinates; it does not form an inverse Cholesky factor or whitened dense
  matrix;
- zero-extends the dominant coarse vector after an H refinement and retains
  it unchanged across a reference-epoch refresh;
- parallelizes the deterministic compact-Gram column scatter for large
  active patch blocks; and
- records eigen iterations, residual, sparse/dense route, warm-start use and
  patch-thread count in both live progress and `iterations.csv`.

The numerator Gram is still dense, so this is not a fully matrix-free
localization operator.  It removes the cubic whitening allocation/work and
the avoidable cold starts; replacing the dense numerator by a patch-operator
Lanczos action is the next optimization boundary if much deeper PALOD remains
too expensive.

On the local 16-thread nine-step S corner-dominant run, the optimized trajectory
completed in 149.34 method seconds with 1484.8 MiB recorded peak memory.  The
older nearby trajectory used 177.68 seconds and 1583.7 MiB, despite ending on
a slightly smaller coarse mesh, so these figures are a conservative
before/after indicator rather than a bitwise-identical timing pair.  At the
last three localization checks the sparse route was active and the spectrum
stage took 0.022, 0.515 and 0.578 seconds.  The final certificate spent 41.52
of 43.84 seconds in ambient Riesz construction: local patch assembly/solves
and dense-Gram production, rather than the generalized spectrum, are now the
dominant deep-PALOD cost.

The post-change six-step stage probe separates this further. At its last
certificate, ambient Riesz took 7.64 seconds: 6.67 seconds were local patch
assembly/factorization/solves plus compact local-Gram formation, while only
0.169 seconds were the deterministic global Gram scatter. Thus additional
scatter parallelism is complete enough for this range; the next optimization
must reduce or reuse the per-patch algebra rather than only parallelize the
final accumulation.

## Remaining optimization boundary

SLOD is now dominated by local correctors and the external synchronized NVB
hierarchy update.  PALOD may still benefit from caching invariant ambient
patch topology and local energy factorizations across an `ell` retry on the
same hierarchy version.  Beyond that, the dense numerator Gram is the next
memory boundary; a patch-operator generalized Lanczos implementation should
be considered before increasing the PALOD horizon far beyond the current
six-epoch server run.

## Reproduction inputs and optional probes

- `configs/S-slod-k16-bottleneck-ell2-step5-direct-schur-v4.json`
- `configs/S-slod-k16-bottleneck-ell2-step5-direct-saddle-v4.json`
- `configs/S-palod-k16-bottleneck-step3-v4.json`
- set `LOD2D_PROFILE_MODEL_STAGES=1` for SLOD model-build stage lines;
- set `LOD2D_PROFILE_LOCALIZATION_STAGES=1` for PALOD localization-certificate
  stage lines.

The local output root used for this study was
`results/S-k16-bottleneck-study-2026-08-13/`; it remains untracked and must not
be copied into formal E1 directories.

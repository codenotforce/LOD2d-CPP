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
only 1.6%, while peak memory increases by about 75 MiB. Eight threads is the
best local exploratory setting on this host. Formal server configurations
remain at their frozen `PATCH_THREADS=4` unless the resource contract is
deliberately revised.

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
was serial.  It must therefore be re-measured on the server before publishing
new time comparisons.  Error-versus-DoF data are thread-count independent;
formal time data must use one newly frozen thread count and common binary.

## Remaining optimization boundary

SLOD is now dominated by local correctors and the external synchronized NVB
hierarchy update.  PALOD may still benefit from caching invariant ambient
patch topology and local energy factorizations across an `ell` retry on the
same hierarchy version, but that is a separate cache-lifetime change and is
not required for the present deep-convergence runs.

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

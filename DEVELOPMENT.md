# LOD2d-C++ Development Log

## Document Role

This is the chronological engineering record. It is the only documentation
file that stores measured timing tables, memory measurements, rejected
experiments, migration defects, and performance conclusions.

- Public quick start and capability summary: [README.md](README.md)
- Helmholtz commands and stable options: [HELMHOLTZ_GUIDE.md](HELMHOLTZ_GUIDE.md)
- Benchmark implementation rules: [BENCHMARK_GUIDE.md](BENCHMARK_GUIDE.md)
- Future work and acceptance gates: the three `*_PLAN.md` files

Do not copy command catalogs or future task lists into this log. A command is
shown here only when it is needed to reproduce a recorded result.

## Current Engineering Summary

| Area | Current decision |
|---|---|
| Elliptic corrector | Eigen LLT remains the default; CHOLMOD and saddle GMRES are explicit experiments |
| Repeated elliptic RHS | Reuse `LodModel`, correctors, basis, and coarse factorization |
| Helmholtz corrector | DirectSaddle is the default; DirectSchur is the leading large-patch experiment |
| Shifted patch GMRES | Correct, but no runtime crossover over direct methods |
| Two-level Schwarz | S4e complete; useful iteration/memory behavior, still experimental |
| Adaptive Helmholtz | Stage-1 calibration implemented; no production estimator frozen |

## Navigation

- Elliptic migration baseline and correctness decisions: sections below
  through **Modular API Layers**.
- Inverse-inequality experiments: **Inverse Inequality Verification**.
- Elliptic performance work and rejected optimizations: **Latest Performance
  Changes** and **Failed or Rejected Experiments**.
- Helmholtz foundation and wave-number work: sections 13-16.
- Adaptive Helmholtz stage 1: section 17.
- Helmholtz patch solver: section 18.
- Coarse and fine-space Schwarz studies: sections 19-27.

## Elliptic Migration Baseline

This section records the established MATLAB-to-C++ elliptic baseline. It is
historical measurement data, not the public quick start.


Environment used for this recorded elliptic baseline:

- WSL2 Ubuntu 22.04
- g++ 11.4
- CMake 3.16+
- Eigen 3.4
- SuiteSparse / CHOLMOD 5.10
- OpenMP, 16 hardware threads available

Validated commands:

```bash
cmake --build build -j 8
./build/tests/test_dg
./build/tests/test_corr --solver=both
./build/tests/test_full
OMP_NUM_THREADS=16 ./build/benchmarks/bench_H4h9 --solver=eigen
OMP_NUM_THREADS=16 ./build/benchmarks/bench_H4h9 --solver=cholmod
./build/benchmarks/bench_refine
```

Latest test results:

- `test_dg`: 10 PASS, 0 FAIL
- `test_corr --solver=both`: 6 PASS, 0 FAIL
- `test_full`: 3 PASS, 0 FAIL

### H=4, h=9, ell=2 (7-run median, OMP_NUM_THREADS=16)

| Version | Median | Range |
|---------|--------|-------|
| **C++ 16t** | **9.79 s** | 9.48-10.01 s |
| MATLAB parallel (4w) | 23.57 s | 22.90-24.07 s |
| MATLAB serial | 49.39 s | 48.19-52.17 s |

Speedup: C++ 5.0x vs MATLAB serial, 2.4x vs MATLAB parallel.

Eigen remains the default corrector solver.

## Correctness Decisions

### Mesh Refinement

The mesh layer supports two conforming bisection rules:

- `refine(mesh)` and `refine_marked(mesh, marked_elements)` use geometric longest-edge bisection (LEB). Marked elements contribute their longest edge to a global split-edge set; all elements incident to a split edge are subdivided consistently.
- `refine_nvb(mesh)`, `bisect_newest_vertex(mesh, marked_elements)`, and `refine_mesh_nvb(mesh, nref)` implement newest-vertex bisection (NVB). The first local vertex is the newest vertex, so local edge (1,2) is the reference edge, matching MATLAB `lod.bisect`.
- NVB local refinement performs a recursive reference-edge closure. The C++ loop only removes a marked edge after that complete edge no longer exists in the current mesh, which also handles initially incompatible newest-vertex labels without hanging nodes.
- Both rules generate `P_node`, `P_elem`, and `P_dg` in the refinement pass. `P_dg` stores affine/barycentric interpolation from each parent DG triangle to child DG vertices.
- `tests/golden_nvb.txt` is derived from the MATLAB NVB implementation in `D:/code/femcode/LOD2d_MATLAB/src/+lod/bisect.m` and is checked by `test_nvb`.
- The old MATLAB red-refinement golden files no longer match the active bisection meshes. `test_mesh`, `test_dg`, `test_patch`, `test_qi`, and `test_nvb` are the default structural/golden tests; `test_corr` and `test_full` are buildable historical red-refinement references unless new bisection golden data is exported.

### DG Assembly

`assemble_dg` now has a split implementation:

- `assemble_element_stiffness(mesh, coeff)` computes and returns all local 3x3
  element stiffness blocks.
- `assemble_dg_from_element_stiffness(blocks)` builds the global DG sparse
  matrix from those blocks.
- `assemble_dg(mesh, coeff)` remains the compatible one-call wrapper.

The coefficient vector length is checked against the element count to avoid
silent out-of-bounds reads in Release builds.

### Saddle GMRES corrector experiment

Previous PCG experiments are not mathematically appropriate for the full corrector saddle system because the matrix is symmetric indefinite, not SPD. The experimental solver `CorrectorSolver::SaddleGmres` instead solves

```text
[Sph  IHp^T] [q     ] = [rhsp]
[IHp   0   ] [lambda]   [0   ]
```

with left-preconditioned GMRES. The preconditioner is the exact block Schur complement inverse using an LLT factorization of `Sph` and a dense LDLT factorization of `IHp*Sph^{-1}*IHp^T`. This is intentionally a correctness-first experiment: it tests the saddle formulation directly, but it still pays for essentially the same local factorization work as the Schur-eliminated Eigen path.

Benchmark command on the MATLAB-compatible red uniform mesh:

```bash
/usr/bin/time -v ./build/benchmarks/bench_saddle_h3h10 --H=3 --h=10 --ell=3 --threads=8 --skip-reference
```

Result on WSL 12 GiB:

| Path | Corrector time | Accuracy vs Eigen Schur | Notes |
|---|---:|---:|---|
| Eigen Schur elimination | 232.283 s | baseline | solves `Sph` for `[IHp^T,rhsp]` and forms the constrained correction by dense Schur solve |
| Saddle GMRES + exact block Schur | 256.436 s | max corrector diff 2.49e-14; max `uHms` diff 4.50e-15 | correct but about 10.4% slower |

Other timings: red mesh 4.18 s, operators 5.22 s, coarse solve 13.79/3.29 s for Eigen/Saddle-built bases, total wall 610.2 s for both corrector paths. Peak RSS was 11.65 GB, so the full reference solve was skipped in this WSL run to avoid OOM risk.

Conclusion: the saddle formulation is correct, but this exact block-Schur GMRES path is not a speed improvement for H=3,h=10,ell=3. It should remain an explicit experiment rather than a default solver. The next useful saddle experiment would need an actually cheaper approximate `Sph^{-1}` or Schur preconditioner; otherwise GMRES mainly adds Krylov overhead around the same factorization.

### Corrector Solver

The corrector still solves the same MATLAB saddle-point formulation:

1. Build coarse patch DOFs.
2. Build fine patch DOFs.
3. Assemble `Sph`.
4. Assemble `rhsp`.
5. Build `IHp`.
6. Solve `Sph \ [IHp', rhsp]`.
7. Compute `mu`.
8. Store the sparse corrector `CTk`.

`CorrectorSolver` supports:

- `EigenLLT`: default, fastest in current benchmarks.
- `Cholmod`: experimental, useful for h=10 corrector benchmarks but higher
  memory than Eigen.
- `CholmodCached`: explicit experiment that reuses a bounded thread-local
  CHOLMOD symbolic factor when the exact local sparsity pattern repeats.
- `SaddleGmres`: explicit saddle-system experiment with exact block Schur preconditioning; correct on H=3,h=10,ell=3 but slower than Eigen Schur.


## Modular API Layers

The LOD setup is now split into three reusable layers:

1. `LodProblemData` owns mesh-derived data: coarse/fine meshes, node incidence
   counts, DG index maps, and refinement prolongation matrices `P_node`,
   `P_elem`, and `P_dg`. Build it with `build_lod_problem_data(initial, H, h)`.
2. `LodOperators` owns coefficient-dependent setup data: element stiffness
   blocks, CG stiffness/mass matrices, patches, interpolation rows, and
   fine-element children. Build it with `build_lod_operators(problem, Ah, ell)`.
3. `LodModel` is the user-facing API for repeated RHS solves. It builds the
   problem data, operators, correctors, multiscale basis, and `LodReusableSystem`
   once, then exposes `solve_from_coarse_values` and `solve_from_fine_values`. It releases setup-only `P_elem` and `P_dg` by default to keep repeated RHS models memory-stable; set `keep_setup_matrices=true` when those matrices must remain inspectable.

Use `LodModel` in examples and application-style benchmarks where only the
right-hand side changes. Use the lower-level `build_lod_*` functions in
profilers that need to time mesh, operator, corrector, and `G` assembly phases
separately. This keeps the public path compact without hiding performance
bottlenecks from `bench_profile`.
## Inverse Inequality Verification

A new benchmark, `bench_inverse_inequality`, numerically checks whether

```text
sup_{v in V_H} H_T ||grad (1-C)v||_{L2(T)} / ||(1-C)v||_{L2(T)}
```

is stable with respect to the coarse mesh size `H_T`.

### Numerical formulation

For the multiscale basis `G = (1-C)P_H` and coefficient vector `a`, write
`w = G a`. On each coarse element `T`, the benchmark assembles local fine-scale
matrices only over the fine children of `T`:

```text
A_T = G_T' S_T G_T
M_T = G_T' M_T G_T
Q_T = H_T * sqrt(lambda_max(A_T, M_T)).
```

`S_T` is the unweighted geometric stiffness matrix for `||grad w||_{L2(T)}`.
This is intentional: the corrector may be built with a heterogeneous coefficient
`A`, but the quantity being tested here is the plain gradient/L2 inverse ratio.
The generalized eigenproblem is solved on the positive mass subspace, dropping
mass eigenvalues below `1e-12 * max(eig(M_T))` to avoid zero-mass directions.


### C=0 sanity check

The benchmark now supports `--basis=coarse`, which sets `G = P_node` and skips
corrector construction. This is the `C=0` case and tests the numerical method
itself rather than the LOD space.

Command:

```bash
./build/benchmarks/bench_inverse_inequality --basis=coarse --sweep-H --H-min=2 --H-max=5 --h-minus-H=5 --ell=2 --coeff=unit --solver=eigen --threads=8 --space=all
```

Results:

| H | h | basis | min | median | p90 | p99 | max |
|---:|---:|---|---:|---:|---:|---:|---:|
| 2 | 7 | coarse | 8.48528 | 8.48528 | 8.48528 | 8.48528 | 8.48528 |
| 3 | 8 | coarse | 8.48528 | 8.48528 | 8.48528 | 8.48528 | 8.48528 |
| 4 | 9 | coarse | 8.48528 | 8.48528 | 8.48528 | 8.48528 | 8.48528 |
| 5 | 10 | coarse | 8.48528 | 8.48528 | 8.48528 | 8.48528 | 8.48528 |

This confirms that the local `S_T/M_T` assembly, element diameter scaling, and
generalized eigenvalue calculation are behaving correctly. With the default
`--space=free`, boundary elements can have `Q_T=0` because all active coarse
basis functions on that element are Dirichlet nodes; `--space=all` removes that
boundary artifact for the sanity check.

### Interpretation caveat

For the fully discrete corrector `C_{H,h,T}^ell : V_H -> W_h^ell(T)`, every
function `(1-C_{H,h}^ell)v_H` is still a fine-grid finite element function.
Therefore a discrete inverse inequality must hold with a constant that can in
principle depend on the fine-scale polynomial space and on the ratio `H/h`.
The benchmark is useful because it measures whether the scaled local constant is
stable in the tested regime, but by itself it does not prove the corresponding
continuous corrector estimate.

A more meaningful numerical probe for the non-discrete corrector is a nested
fine-grid convergence experiment: fix `H`, `ell`, and `A`; compute the corrector
on increasingly fine meshes `h = H+m`; measure `Q_{H,h}`; and check whether
`Q_{H,h}` converges to a finite limit as `h -> 0`. If the limit is stable and
then remains stable under a separate `H` sweep, that is better evidence for the
continuous corrector version. This is still not a proof, but it tests the right
limit instead of only one discrete finite element space.
### Baseline H sweep

Command:

```bash
./build/benchmarks/bench_inverse_inequality --sweep-H --H-min=2 --H-max=4 --h-minus-H=5 --ell=2 --coeff=unit --solver=eigen --threads=8
```

Results on WSL, unit coefficient, free coarse-node space:

| H | h | ell | coarse elems | fine elems | min | median | p90 | p99 | max |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 7 | 2 | 32 | 32768 | 17.2464 | 32.2145 | 54.2416 | 54.2416 | 54.2416 |
| 3 | 8 | 2 | 128 | 131072 | 17.2456 | 72.3493 | 99.6182 | 101.896 | 101.896 |
| 4 | 9 | 2 | 512 | 524288 | 17.2456 | 94.0825 | 101.895 | 101.896 | 101.896 |

Interpretation: the supremum `max_T Q_T` stabilizes from `H=3` to `H=4` at
about `101.9`; the `H=2` mesh is too coarse to show the asymptotic value. This
supports mesh-size independence of the tested inverse ratio for the unit
coefficient in the observed range, although the constant is much larger than a
plain coarse P1 inverse constant because `(1-C)V_H` contains fine-scale corrector
structure inside each coarse element.

### Oversampling sensitivity

Command family:

```bash
./build/benchmarks/bench_inverse_inequality --H=4 --h=9 --ell=1 --coeff=unit --solver=eigen --threads=8
./build/benchmarks/bench_inverse_inequality --H=4 --h=9 --ell=2 --coeff=unit --solver=eigen --threads=8
./build/benchmarks/bench_inverse_inequality --H=4 --h=9 --ell=3 --coeff=unit --solver=eigen --threads=8
```

Results:

| H | h | ell | min | median | p90 | p99 | max |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 4 | 9 | 1 | 7.62808 | 35.0801 | 35.128 | 35.128 | 35.128 |
| 4 | 9 | 2 | 17.2456 | 94.0825 | 101.895 | 101.896 | 101.896 |
| 4 | 9 | 3 | 34.6497 | 100.348 | 103.405 | 103.413 | 103.42 |

Interpretation: increasing `ell` from 1 to 2 sharply increases the local inverse
constant, then `ell=3` changes the supremum only mildly. This suggests the
single-element inverse ratio sees additional fine-scale corrector structure as
oversampling grows, but appears to plateau for this unit-coefficient case.


### Nested h-refinement probe for the continuous corrector

To probe the non-discrete corrector, fix `H`, `ell`, and `A`, then refine only
the corrector mesh. This checks whether the discrete constants `Q_{H,h}` appear
to converge as `h -> infinity` before using them in an `H` sweep.

Commands:

```bash
./build/benchmarks/bench_inverse_inequality --sweep-h --H=3 --h-min=6 --h-max=9 --ell=2 --basis=lod --coeff=unit --solver=eigen --threads=8
./build/benchmarks/bench_inverse_inequality --H=3 --h=10 --ell=2 --basis=lod --coeff=unit --solver=auto --threads=8
```

Results for unit coefficient, free coarse-node space:

| H | h | ell | fine elems | median | p90 | p99 | max |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 6 | 2 | 8192 | 45.1931 | 49.3655 | 49.8567 | 49.8567 |
| 3 | 7 | 2 | 32768 | 62.2412 | 78.296 | 78.6532 | 78.6532 |
| 3 | 8 | 2 | 131072 | 72.3493 | 99.6182 | 101.896 | 101.896 |
| 3 | 9 | 2 | 524288 | 76.4449 | 109.619 | 113.477 | 113.477 |
| 3 | 10 | 2 | 2097152 | 77.6661 | 112.885 | 117.342 | 117.342 |

Interpretation: for fixed `H=3`, `ell=2`, and `A=1`, the local inverse constant
continues to increase as the corrector mesh is refined, but the increment from
`h=9` to `h=10` is much smaller than from earlier refinements. This is evidence
that the discrete constants may be approaching a finite continuous-corrector
limit near this value. The next useful check is to run the same nested
`h`-refinement for another `H`, then compare the apparent limiting constants.

This experiment is more relevant to the continuous corrector than a single
fully discrete run, because it tests the behavior of the sequence
`C_{H,h}^ell` as the fine discretization is refined while `H` and `ell` are
held fixed.


### Oversampling-limit probe for fixed H and h

The benchmark now supports `--sweep-ell`, which fixes `H`, `h`, and `A` while
increasing the oversampling radius `ell`. This probes whether the local inverse
constant approaches a plateau as the localized corrector approaches the global
corrector on the fixed fine grid.

Local smoke-test command:

```bash
./build/benchmarks/bench_inverse_inequality --sweep-ell --H=3 --h=8 --ell-min=1 --ell-max=4 --basis=lod --coeff=unit --solver=eigen --threads=4
```

Results for unit coefficient, free coarse-node space:

| H | h | ell | median | p90 | p99 | max |
|---:|---:|---:|---:|---:|---:|---:|
| 3 | 8 | 1 | 30.3253 | 35.1118 | 35.128 | 35.128 |
| 3 | 8 | 2 | 72.3493 | 99.6183 | 101.896 | 101.896 |
| 3 | 8 | 3 | 81.2125 | 101.579 | 103.42 | 103.42 |
| 3 | 8 | 4 | 82.4034 | 101.243 | 103.332 | 103.333 |

Interpretation: for this fixed `H=3,h=8` case, the maximum constant jumps from
`ell=1` to `ell=2`, then essentially plateaus for `ell=3,4`. This suggests a
finite oversampling-limit constant on the fixed fine grid.

Recommended server command for the same experiment at the resolved `h=12` level:

```bash
THREADS=8 MODE=ell H=3 H_FIXED=12 ELL_MIN=1 ELL_MAX=5 COEFF=unit SOLVER=auto bash scripts/run_inverse_server.sh
```

Use the server logs to compare `max Q_T` and peak RSS across `ell`. Very large
`ell` may become expensive because local patches grow and CHOLMOD factorization
cost increases nonlinearly.


### Corrector cache for repeated inverse experiments

`bench_inverse_inequality` now stores compact corrector entries on disk. This is
useful because the correctors depend on `H`, `h`, `ell`, `A`, and the solver,
but not on whether the inverse-inequality numerator is `grad (1-C)v` or
`grad C v`.

Default cache directory:

```text
results/corrector_cache
```

Sanity test on WSL:

```bash
rm -rf /tmp/lod_ct_cache_test
./build/benchmarks/bench_inverse_inequality --H=3 --h=9 --ell=2 --basis=lod --numerator=lod --coeff=unit --solver=eigen --threads=4 --cache-dir=/tmp/lod_ct_cache_test
./build/benchmarks/bench_inverse_inequality --H=3 --h=9 --ell=2 --basis=lod --numerator=corrector --coeff=unit --solver=eigen --threads=4 --cache-dir=/tmp/lod_ct_cache_test
./build/benchmarks/bench_inverse_inequality --H=4 --h=9 --ell=2 --basis=lod --numerator=corrector --coeff=unit --solver=eigen --threads=4 --cache-dir=/tmp/lod_ct_cache_test
```

Observed behavior:

| Case | Cache | Setup | Max Q |
|---|---|---:|---:|
| H=3,h=9,ell=2,numerator=lod | miss/write | 13.65 s | 113.477 |
| H=3,h=9,ell=2,numerator=corrector | hit | 2.74 s | 13292 |
| H=4,h=9,ell=2,numerator=corrector | hit after prior write | 3.18 s | 13350.6 |

The
umerator=corrector` quotient is much larger than the original LOD-basis
inverse quotient. The benchmark computes the generalized quotient on the
positive local mass subspace of `(1-C)V_H`; exact zero-denominator directions
are projected out, matching the earlier inverse-inequality experiments.

### Fixed h and ell, H-refinement probe

The server helper supports `MODE=H`, which fixes the fine resolution `h` and
oversampling radius `ell`, then sweeps the coarse level `H`. This is the direct
numerical test for whether the already-resolved right-hand side constant is
independent of `H`:

```bash
THREADS=8 MODE=H H_FIXED=12 ELL=3 H_MIN=2 H_MAX=4 COEFF=unit SOLVER=auto bash scripts/run_inverse_server.sh
```

Interpretation guide:

- Keep `H_FIXED` large enough that each tested `H` has a resolved corrector.
- Choose `ELL` from a plateau observed in `MODE=ell` first; for the unit test
  case, `ell=3` is a reasonable first choice after the `H=3,h=8` smoke test.
- Compare the `max Q_T` column across `H`. If it stays in the same range instead
  of growing systematically as `H` is refined, this supports an `H`-independent
  inverse constant for the numerically resolved ideal/localized corrector.

For a cheaper first pass, run:

```bash
THREADS=8 MODE=H H_FIXED=10 ELL=3 H_MIN=2 H_MAX=4 COEFF=unit SOLVER=auto bash scripts/run_inverse_server.sh
```

Then repeat with `H_FIXED=12` on the server if the trend is unclear.
### High-memory server plan

The WSL 12 GiB run killed `H=3,h=11` at about 11.28 GB RSS before completion.
The EPYC server reported 377 GiB total memory and 366 GiB available, so the full
benchmark can be pushed further without immediately rewriting the algorithm.

Recommended first server run:

```bash
THREADS=32 H=3 ELL=2 H_MIN=6 H_MAX=12 COEFF=unit SOLVER=auto ./scripts/run_inverse_server.sh
```

Expected feasibility on that server:

| Case | Fine elements | Recommendation |
|---|---:|---|
| H=3,h=11 | 8,388,608 | Should fit; use 32 threads first |
| H=3,h=12 | 33,554,432 | Should fit in memory; expect a long run |
| H=3,h=13 | 134,217,728 | Experimental; try only after h=12 succeeds |
| H=3,h>=14 | >=536,870,912 | Not recommended with the full-global implementation |

The script runs each `h` as a separate process so a failed high-h case does not
lose lower-h results. Each log includes `/usr/bin/time -v`, especially maximum
resident set size.
### Next checks

The current evidence is positive for `A=1`, but not a proof. The next numerical
checks should use:

- `--coeff=file:benchmarks/data_H4h10.txt` for the MATLAB-exported benchmark coefficient,
- `--coeff=checkerboard:1000` for high-contrast stress testing,
- a larger sweep such as `H=5,h=10` when memory/time allow.
## Latest Performance Changes

### 1. Element stiffness block cache

Previous corrector code recovered each fine element's 3x3 stiffness block by
iterating over the global sparse `Shdg` columns.  This happened in both `Sph`
and `rhsp` assembly for every corrector.

New code reuses `ElementStiffnessBlocks`, computed once during DG assembly.
This removes repeated sparse column scans and keeps the global `Shdg` matrix
available for the rest of the pipeline.

### 2. Fine element children cache

Previous corrector code computed:

```text
P0 * patch(:, k)
```

for every coarse element `k`, then scanned all fine elements to identify the
local patch. This repeatedly scanned the complete fine-element set once per
coarse corrector.

New code builds `FineElementChildren` once from `P_elem`, mapping:

```text
coarse element -> fine child elements
```

The corrector now expands only the coarse elements present in `patch(:, k)`.
The target RHS set still correctly uses only entries with `patch(:, k) > 1`.
Diagnostics confirmed that this is not the same as the full patch:

```text
k=0   fine_patch=2304  fine_target=256
k=63  fine_patch=9472  fine_target=256
```

### 3. CHOLMOD wrapper repair

The previous wrapper used `cholmod_l_*` and borrowed Eigen's internal sparse
arrays.  That was unsafe because of integer-width and ownership issues.

The new wrapper:

- uses the 32-bit `cholmod_*` API,
- copies Eigen sparse data into CHOLMOD-owned triplets,
- marks the triplet and sparse matrix as lower-triangular SPD storage,
- checks analyze/factorize/solve failures,
- returns dense multi-RHS solutions.

CHOLMOD now passes golden corrector tests.  Plain CHOLMOD remains the h>=10
benchmark option, while cached CHOLMOD is kept as an explicit experiment only.

### 4. Benchmark build fix

`bench_refine.cpp` previously contained only a comment while CMake built it as
an executable.  This caused full builds to fail at link time.  It is now a real
mesh-refinement timing benchmark.

### 5. Multi-RHS Eigen solve

The previous Eigen corrector branch factorized `Sph` once, but solved the
dense RHS block column by column:

```cpp
for (int jj = 0; jj < nd+d+1; ++jj)
    X.col(jj) = llt.solve(RHS.col(jj));
```

For H4/h9, the local patch matrices and the number of RHS columns are large
enough that this dominated runtime.  The corrector now uses Eigen's dense
multi-RHS solve directly:

```cpp
X = llt.solve(RHS);
```

This was the largest single speedup and reduced the H4/h9 corrector phase from
the previously observed 90 s range to single-digit seconds on the same
16-thread benchmark.

### 6. Corrector scratch-buffer reuse

The corrector previously allocated and cleared `Nh`-sized vectors for every
coarse element.  At H4/h9, that means 512 correctors repeatedly touching arrays
with 263169 entries, under OpenMP.

The current implementation uses thread-local marker arrays with integer stamps
for:

- fine vertex counts inside the current patch,
- global fine DOF to local patch DOF lookup.

Only vertices touched by the current patch are scanned.  This keeps the
mathematics unchanged and removes a large amount of allocator and memory
bandwidth overhead.

### 7. Local sparse insertion cleanup

`IHp` is now built as a dense local block because it is consumed immediately in
dense products.  `CTk` is assembled through triplets instead of repeated sparse
`insert` calls.  These changes are smaller than the multi-RHS solve but help
stabilize the OpenMP benchmark timings.


### 8. Bounded CHOLMOD symbolic factor cache

`solve_cholmod_cached` keeps a thread-local `cholmod_common` and a bounded map
from exact local sparse pattern to `cholmod_factor`.  Numeric factorization is
still refreshed for every corrector; only CHOLMOD's symbolic analysis is reused.

The cache is intentionally capped at one pattern per OpenMP thread.  An
unbounded cache retained too many large factors at H=4,h=10 and was killed by
WSL after reaching about 11.6 GB RSS on a 12 GB configuration.  The bounded
version is correct, but the measured h=10 profile was slower than plain CHOLMOD
because the current dynamic patch order has few immediate exact pattern hits.
Keep it behind `--solver=cholmod_cached` until patch grouping or a better reuse
policy proves a net speedup.

### 9. Thread-local local triplet buffers

The corrector now reuses the local `sph_t` and `rhs_t` triplet buffers per
OpenMP thread.  This avoids repeated large vector allocations while preserving
the same Eigen `SparseMatrix::setFromTriplets` path and numerical behavior.

Measured H=4,h=10 results were mixed because total time is noisy, but the
corrector phase was slightly lower in profile runs.  This is kept as a small,
low-risk allocation cleanup rather than a headline solver speedup.


### 10. Reusable LOD system for multiple RHS values

For fixed `A`, mesh, `H/h`, `ell`, patches, and interpolation, the element
correctors do not depend on the right-hand side.  `LodReusableSystem` stores the
assembled multiscale basis `G`, the free coarse basis `G0`, the mass matrix, the
coarse prolongation, and the factorization of `SHLOD0 = G0' * Sh * G0`.

Use this when solving many problems with the same coefficient and mesh but
different forcing terms.  The expensive corrector stage is paid once; each new
RHS only computes `G0' * Mh * f_fine`, solves the already-factorized coarse
system, and evaluates `uHms = G * uH`.

Validation benchmark:

```bash
./build/benchmarks/bench_reuse_rhs --solver=auto --rhs=5
```

Observed H=4,h=10 results on WSL:

- reusable setup: about 46-50 s,
- correctors inside setup: about 32-34 s,
- repeated RHS solves: about 50-100 ms per RHS,
- peak RSS: about 6.6-7.2 GB.

### 11. LodModel high-level API

`LodProblemData` and `LodOperators` remove the repeated setup code that had
spread across benchmark drivers. `LodModel` composes those pieces with
`LodReusableSystem`, giving users a single build step for fixed coefficient and
mesh cases followed by cheap repeated RHS solves.

Validation after this refactor:

```bash
cmake --build build --target lod2d_core bench_reuse_rhs test_corr test_full -j 8
./build/tests/test_corr --solver=both
./build/tests/test_full
./build/benchmarks/bench_reuse_rhs --solver=auto --rhs=2
```

Latest `bench_reuse_rhs --solver=auto --rhs=2` result on WSL H=4,h=10:

- reusable setup: 48.15 s,
- repeated RHS total: 189.30 ms,
- repeated RHS average: 94.65 ms.

### 12. Corrector pipeline helper extraction

The benchmark and full-test pipelines now use two shared helpers instead of
copying the same loops in every driver:

- `compute_all_correctors(...)` owns the OpenMP corrector loop and solver
  dispatch for all coarse elements.
- `build_multiscale_basis(...)` owns direct assembly of `G = P_node - C_ell`
  from compact corrector entries.

This is an organizational refactor, not a mathematical change.  It reduces the
chance that future optimizations, scheduling experiments, or bug fixes are
applied to one benchmark but forgotten in another.
## Failed or Rejected Experiments

| Experiment | Result |
|------------|--------|
| CHOLMOD as default corrector solver | Correct but slower than Eigen in the earlier small-case comparison |
| CHOLMOD for H4/h9 | Faster corrector phase than Eigen in one run, but slower total runtime |
| Unbounded CHOLMOD factor cache | Reached about 11.6 GB RSS and was killed on the 12 GB WSL machine |
| Bounded CHOLMOD factor cache | Correct and memory-safe, but h=10 profile was slower than plain CHOLMOD with current patch order |
| Precomputed corrector patch plans as default | Reduced h=10 corrector time slightly in one profile, but increased setup/RSS and did not improve total time robustly |
| Patch signature grouped corrector scheduling | Sorting by patch column signature caused pathological h=10 profile runtime and was reverted |
| Simultaneous `G`/`G0` triplet construction | Avoided one `G` scan in theory, but increased triplet pressure and caused pathological h=10 profile runtime; reverted |
| Direct CHOLMOD construction from local triplets | Correct, but was not consistently faster than the existing Eigen sparse path in h=10 profiles |
| `IHp` sparse iterator replacement | Previously broke golden data; kept `coeff()` extraction |
| `symrcm` on each `Sph` | Increased overhead/fill for tested patch matrices |
| Precomputing full sparse submatrices | Too much serial precompute/broadcast overhead |
| MATLAB GPU sparse path | Sparse GPU indexing does not support the needed submatrix access |
| PCG + incomplete Cholesky | Multi-RHS local solves favored direct Cholesky |

## Research Roadmap

This log does not duplicate detailed future task lists. The remaining public
research directions are:

- PML and additional Helmholtz boundary models;
- a production policy for shifted-GMRES patch solves;
- adaptive LOD beyond the current stage-1 calibration gate.

## 13. Helmholtz Petrov-Galerkin Foundation

The Helmholtz implementation is deliberately separate from the real SPD
elliptic corrector path while sharing mesh, NVB, patch, and quasi-interpolation
utilities.

### Module ownership

- `include/helmholtz/types.h`: complex scalar, vector, sparse matrix, and
  callback aliases.
- `include/helmholtz/operators.h`, `src/helmholtz/operators.cpp`: volume
  stiffness/mass, physical Robin boundary mass, complex load quadrature,
  direct fine solve, and error norms.
- `include/helmholtz/corrector.h`, `src/helmholtz/corrector.cpp`: local patch
  constraints, complex saddle solve, primal/adjoint correctors, diagnostics,
  and corrected basis assembly.
- `include/helmholtz/model.h`, `src/helmholtz/model.cpp`: global NVB hierarchy,
  Petrov-Galerkin spaces, coarse factorization, reconstruction, and repeated
  right-hand-side solves.

### Algebraic conventions

For coefficient vectors `x` and `y`, the variational form is represented as
`y.adjoint() * A * x`. The fine operator is complex symmetric for the current
real coefficients, but is not Hermitian. Generic complex code must use
`adjoint()`; `transpose()` is only valid where complex symmetry is explicitly
being tested.

The trial and test bases are

```text
Psi_trial = (I - C) P,
Psi_test  = (I - C*) P,
A_LOD     = Psi_test^* A_h Psi_trial.
```

For the present real mesh, coefficients, and interpolation matrix,
`C* v = conjugate(C conjugate(v))`. The code still stores primal and adjoint
correctors separately so a future general complex/PML operator can replace
this fast path without changing the model API.

### Patch boundary invariant

A fine vertex is a local patch unknown only when all of its incident domain
elements belong to the patch. This retains physical Robin boundary vertices
when their full domain star is present and removes vertices on the artificial
patch boundary. Robin edge blocks are stored in the fine element operator, so
a physical boundary edge is included exactly once whenever its element belongs
to the patch.

### NVB compatibility

NVB interprets local vertex 0 as the newest vertex and local edge `(1,2)` as
the reference edge. The Helmholtz unit-square initial mesh is therefore
`{0,1,3}` and `{2,3,1}`: both triangles reference the shared diagonal. The
older ordering `{0,1,3}`, `{1,2,3}` is not a compatible initial labeling for
repeated global NVB and must not be used for this hierarchy.

### Solver and reuse policy

The correctness baseline uses `Eigen::SparseLU` for local complex saddle
systems, the coarse Petrov-Galerkin system, and the fine reference system.
Local corrector failures are captured inside the OpenMP region and rethrown
after it, so parallel construction preserves clear solver diagnostics.
The built model retains the coarse factorization, so changing only `f` does not
recompute correctors or refactor the coarse operator. Changing `k`, mesh,
coefficients, boundary data, interpolation, or `ell` requires rebuilding.

### Validation on 2026-06-29

```text
Manufactured FEM energy rate : 0.985374
Manufactured FEM L2 rate     : 1.97936
Max primal residual          : 3.10527e-16
Max adjoint residual         : 3.10527e-16
Max I_H constraint residual  : 1.70725e-16
CTest                         : 9/9 passed
```

Both the two-sided and corrected-test-only Petrov-Galerkin systems agree with
independent dense coarse solves on the small verification mesh.

### Manufactured global-NVB convergence on 2026-07-19

`bench_helmholtz_manufactured` provides a reproducible convergence gate for
the homogeneous impedance Robin problem

```text
-Delta u - k^2 u = f,
partial_n u - i k u = 0 on the boundary.
```

The run used `k=4`, the two-sided Petrov-Galerkin formulation, a fixed fine
level `h=10`, coarse levels `H=3,4,5,6`, and oversampling levels `ell=2,3`.
The manufactured solution was

```text
u(x,y) = phi(x) phi(y) exp(i k x),
phi(t) = 16 t^2 (1-t)^2.
```

Both `phi` and `phi'` vanish at the endpoints, so the exact solution satisfies
the homogeneous Robin condition on all four sides. The source is evaluated by
the shared `make_polynomial_plane_wave_solution` implementation. The separate
fine-P1 calibration used global NVB levels 6, 8, and 10:

| NVB level | nodes / elements | measured `h_max` | energy error | energy rate | L2 error | L2 rate | relative linear residual |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 6 | 81 / 128 | 1.767767e-1 | 8.320692e-1 | - | 4.245765e-2 | - | 1.03e-15 |
| 8 | 289 / 512 | 8.838835e-2 | 4.176412e-1 | 0.9944 | 1.121713e-2 | 1.9203 | 4.27e-15 |
| 10 | 1,089 / 2,048 | 4.419417e-2 | 2.088422e-1 | 0.9999 | 2.843028e-3 | 1.9802 | 1.77e-14 |

The expected P1 rates are recovered: first order in the Helmholtz energy norm
and approximately second order in `L2`. On the fixed level-10 fine grid, the
LOD results were:

| `H` level | measured `H_max` | `ell` | exact energy error | local energy rate | exact L2 error | local L2 rate | LOD-to-fine energy | LOD-to-fine L2 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 5.000000e-1 | 2 | 1.448194e+0 | - | 1.259591e-1 | - | 1.428401e+0 | 1.242684e-1 |
| 4 | 3.535534e-1 | 2 | 8.549633e-1 | 1.5206 | 5.531563e-2 | 2.3744 | 8.259847e-1 | 5.376757e-2 |
| 5 | 2.500000e-1 | 2 | 5.303743e-1 | 1.3777 | 2.290905e-2 | 2.5435 | 4.856312e-1 | 2.140134e-2 |
| 6 | 1.767767e-1 | 2 | 3.304008e-1 | 1.3656 | 8.754737e-3 | 2.7756 | 2.550294e-1 | 7.244207e-3 |
| 3 | 5.000000e-1 | 3 | 1.449053e+0 | - | 1.260970e-1 | - | 1.429264e+0 | 1.244059e-1 |
| 4 | 3.535534e-1 | 3 | 8.550645e-1 | 1.5220 | 5.534252e-2 | 2.3761 | 8.260938e-1 | 5.379733e-2 |
| 5 | 2.500000e-1 | 3 | 5.294103e-1 | 1.3833 | 2.288520e-2 | 2.5479 | 4.845916e-1 | 2.138612e-2 |
| 6 | 1.767767e-1 | 3 | 3.284277e-1 | 1.3776 | 8.666529e-3 | 2.8018 | 2.524960e-1 | 7.168535e-3 |

Across the full `H=3` to `H=6` range, the exact-solution energy/L2 rates are
approximately `1.42/2.56` for `ell=2` and `1.43/2.58` for `ell=3`. The two
oversampling choices are nearly indistinguishable on this smooth test, so
localization is not the dominant error. Petrov residuals were below
`1.1e-14`, primal corrector residuals below `2.2e-15`, and interpolation
constraint residuals below `6.1e-16`.

NVB level numbers must not be reported as literal powers of the geometric
meshwidth. Each global NVB sweep bisects element area, and this initial mesh
satisfies approximately

```text
max_element_diameter(level j) = 2^((1-j)/2).
```

Consequently the requested project parameters `H=3,...,6`, `h=10` have actual
diameters shown above; a literal geometric `h approximately 2^-10` would need
about NVB level 21 and roughly 4.2 million triangles. That high-memory case is
not part of this WSL validation.

The automated `--check` gate verifies homogeneous Robin data, FEM rates,
algebraic residuals, and decreasing LOD-to-exact and LOD-to-fine errors for
each `ell`. Results and metadata are stored in
`results/helmholtz_manufactured/`. The archived run passed, and all six
registered `helmholtz_` CTest targets also passed.

## 14. Helmholtz Wave-Number Scan

`bench_helmholtz_k` implements the fifth Helmholtz milestone. For every wave
number it builds an independent globally refined NVB hierarchy, chooses the
first coarse level satisfying `kH <= target`, sets `ell = ceil(log2(k))`, and
compares:

- standard coarse P1 Galerkin FEM;
- two-sided Petrov-Galerkin LOD;
- the fine P1 reference solution used by both error calculations.

The fixed test source is

```text
f(x,y) = exp(-40 * ((x - 0.35)^2 + (y - 0.55)^2)).
```

The reported energy norm is `sqrt(v^* (K + k^2 M) v)`. The optional inf-sup
value uses exact trial/test energy Gram matrices and a dense SVD, so it is only
computed below `--stability-max-dofs`; larger cases report
an` instead of an
unscaled or misleading surrogate.

### WSL validation scan on 2026-06-29

This local validation used `kH = 1`, `kh = 0.25`, `fine-gap = 4`, and
`ell = ceil(log2(k))`. It validates the experiment and error pipeline, but the
finer server scan should use the script default `fine-gap = 8` before drawing a
final pollution conclusion.

| k | coarse/fine elements | FEM energy rel. | LOD energy rel. | LOD L2 rel. | correctors | total |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | 64 / 1,024 | 2.115e-1 | 7.298e-2 | 1.934e-2 | 0.072 s | 0.089 s |
| 8 | 256 / 4,096 | 1.585e-1 | 1.771e-2 | 3.636e-3 | 0.871 s | 1.363 s |
| 16 | 1,024 / 16,384 | 1.561e-1 | 8.811e-3 | 1.545e-3 | 10.318 s | 10.833 s |
| 32 | 4,096 / 65,536 | 3.449e-2 | 1.140e-2 | 1.544e-3 | 119.906 s | 124.193 s |

All four points had Petrov residuals below `4.3e-16`, primal corrector
residuals below `1.4e-15`, and constraint residuals below `3.8e-16`. The
`k=32` process peaked at about 2.11 GB RSS. Corrector construction occupied
roughly 96.5% of its total runtime and is the clear target for milestone 6.

These data show that the LOD error remains bounded over the tested range, but
they do not by themselves establish a clean pollution slope: the fixed-source
standard FEM error is not monotone in `k`. The full `k=4,...,64`,
`fine-gap=8` server run and a manufactured oscillatory solution are useful
follow-up checks before making a stronger numerical claim.

The full scan is run as separate processes:

```bash
K_VALUES="4 8 16 32 64" FINE_GAP=8 KH_TARGET=1 \
  bash scripts/run_helmholtz_k_scan.sh
```

Each process writes one CSV row and one `/usr/bin/time -v` log. This avoids
retaining the previous wave number's mesh, correctors, and sparse factors.

## 15. Helmholtz Corrector Performance

Milestone 6 keeps the complex SparseLU gold-standard solve but reduces its
assembly and scheduling overhead:

- sort patches by estimated fine-element cost before OpenMP dynamic scheduling;
- reuse stamped thread-local node, element, and constraint workspaces;
- assemble only active interpolation rows in each local saddle system;
- solve the three local coarse-basis right-hand sides as one block;
- cache symbolic analysis when a thread encounters the exact same sparse pattern;
- capture worker exceptions and rethrow them after the parallel region.

The following Release measurements used `fine-gap=4`, `kH=1`, and eight OpenMP
threads. Errors and residuals were unchanged from the serial milestone-5 run.

| k | milestone-5 correctors | optimized correctors | optimized total | speedup |
|---:|---:|---:|---:|---:|
| 16 | 10.318 s | 1.653 s | 2.119 s | 6.24x |
| 32 | 119.906 s | 18.551 s | 22.790 s | 6.46x |

For `k=32`, peak RSS was 2.14 GB versus approximately 2.11 GB before the
optimization. The optimized errors remained `3.4485076e-2` (coarse FEM energy),
`1.1401850e-2` (LOD energy), and `1.5435537e-3` (LOD L2), with corrector
residuals near `1e-15`. Symbolic reuse is useful in serial runs but contributes
less with many workers because each thread owns its cache; parallelism and
reduced local assembly are the main gains.

## 16. Helmholtz Reproduction Infrastructure

This milestone added configuration-aware resume keys, explicit OpenMP
placement, per-wave-number process isolation, metadata, CSV output, and
`/usr/bin/time` resource logs. These facts are retained here because they
explain the recorded measurements. Current commands and options are maintained
only in [HELMHOLTZ_GUIDE.md](HELMHOLTZ_GUIDE.md).

## 17. Adaptive Helmholtz Stage-1 Baseline

The full-rebuild adaptive baseline now supports a fixed master fine mesh,
locally refined NVB coarse meshes, residual reconstruction, three strong
residual aggregations, local energy-Riesz calibration, deterministic Dorfler
marking, reliability-envelope utilities, and per-iteration inf-sup diagnostics.
All 11 CTest targets pass. The integrated complex residual agrees with
`b-Au` to about `2.53e-17` on the Debug test.

A Release pilot with `k=4`, level `H=5`, `h=10`, `ell=3`, and four threads
produced relative energy errors `0.122562, 0.155091, 0.108361` on
`64, 70, 88` coarse elements. The corresponding inf-sup values were
`0.318864, 0.558892, 0.664904`; therefore the nonmonotone transition is not an
inf-sup collapse. Repeating with `ell=6` was essentially unchanged.

The candidate fine-indicator effectivity grew from `2.27` to `12.19` when it
was divided by `||u_h-u_LOD||`. This only shows that a strong residual is not
an estimator of the purely discrete difference, which vanishes when
`u_LOD` reaches `u_h`.

The benchmark now supports `--source=manufactured` with
`u=phi(x)phi(y)exp(i*k*x)` and `phi(t)=16*t^2(1-t)^2`. For
`k=4,H=5,h=10,ell=3`, seven adaptive iterations reduced the exact LOD energy
error from `0.529410` to `0.214176`; the fixed fine-grid exact error remained
`0.208842` in every iteration. The fine estimator decreased from `1.26604`
to `1.22683`, and its continuous-error effectivity increased from `2.39142`
to `5.72816`. Wider parameter studies are required before making a uniform
reliability or efficiency claim.

The manufactured test exposed a generic NVB defect: recursive closure could
create separate global node indices at the same shared-edge midpoint. The
duplicate degrees of freedom changed the nominally fixed fine space after each
adaptive step. The NVB core now merges coincident nodes and rebuilds `P_node`;
regression tests check unique coordinates and fixed master-fine node counts.

For an estimator of `u_h-u_LOD`, an algebraic dual-residual or hierarchical
two-level estimator remains the appropriate next candidate.

## 18. Helmholtz Patch Right-Preconditioned GMRES (M1-M4)

Implemented on 2026-07-19. Adaptive Helmholtz work remains paused after
stage 1; this change is confined to the local corrector solver.

### Architecture

- `HelmholtzPatchAssembler` owns patch topology, free-DOF selection, physical
  Robin/artificial Dirichlet boundaries, restricted `K/M/R`, independent
  quasi-interpolation constraints, and the three element right-hand sides.
- `solve_helmholtz_patch` owns `DirectSaddle`, `DirectSchur`, and
  `ShiftedGmres`. Corrector scheduling and compact output are separate.
- `solve_right_preconditioned_gmres` is a complex restarted implementation
  with modified Gram-Schmidt, optional reorthogonalization, and explicit true
  residual checks after every update.
- `HelmholtzProblemConfig::patch_solver` exposes the experiment without
  changing the default `DirectSaddle` behavior.
- `HelmholtzPatchVcycle` builds Galerkin coarse shifted operators from the
  uniform NVB hierarchy, keeps the patch boundary semantics on every level,
  uses fixed weighted-Jacobi smoothing, and applies SparseLU only on the
  selected coarse patch grid.

`DirectSchur` computes `Y=A^-1 F`, `Z=A^-1 B*`, solves
`(B Z) Lambda=B Y`, and returns `X=Y-Z Lambda`. `ShiftedGmres` replaces each
`A^-1` application by right-preconditioned GMRES with
`P=K-(k^2+i epsilon)M-i k R`. SparseLU remains the exact shifted-inverse
reference. M4 can instead apply a fixed geometric V-cycle. Uniform global NVB
data supplies the hierarchy; adaptive/non-nested data rejects this path
explicitly.

### Correctness

- All 16 registered Debug CTests pass, including both Schwarz tests.
- The standalone nonnormal complex test covers restarted GMRES, zero RHS,
  and one-iteration convergence with an exact right preconditioner.
- On `H=1,h=4,ell=1,k=1.5`, DirectSchur agrees with DirectSaddle below
  `1e-10`; both exact-shifted and V-cycle ShiftedGmres agree below `1e-8`
  with no fallback.
- Primal, adjoint, interpolation-constraint, Schur, and true GMRES residuals
  are all checked. Failures throw unless explicit fallback is enabled.
- The user-added global-NVB manufactured benchmark still runs through the
  default DirectSaddle path; a small `k=2,h=6,H=2,3` regression completed
  with Petrov and corrector residuals near machine precision.

### Release measurements

Both runs used `epsilon=0.2 k^2`, GMRES tolerance `1e-10`, and no fallback.

| case | DirectSaddle corrector | DirectSchur corrector | ShiftedGmres corrector | GMRES avg/max | GMRES basis difference |
|---|---:|---:|---:|---:|---:|
| `H=1,h=5,ell=1,k=4,1t` | 0.550 ms | 0.439 ms | 2.579 ms | 7.625 / 8 | 2.05e-11 |
| `H=2,h=7,ell=2,k=4,4t` | 2.706 ms | 1.383 ms | 8.230 ms | 7.500 / 8 | 1.20e-11 |

DirectSchur is 1.25-1.96 times faster than DirectSaddle in these small M3
tests. Exact-shifted GMRES is roughly 5-6 times slower because every Krylov
step applies a sparse triangular solve and adds orthogonalization and
true-residual work. M4 below tests whether a geometric V-cycle changes that
conclusion. `DirectSaddle` remains the default while DirectSchur still needs
larger wave-number scans for robustness against local `A_omega` resonances.

### M4 V-cycle scale study

The following Release runs used `H=5`, `ell=3`, `k=4`, eight OpenMP
threads, `epsilon=0.2 k^2`, two pre/post weighted-Jacobi steps,
`omega=0.6`, coarse threshold 200, GMRES tolerance `1e-10`, and no
fallback. Here `h` is the global NVB sweep count, not a mesh width.

| h | global fine triangles | largest patch DOF | DirectSaddle | DirectSchur | GMRES+shifted LU | GMRES+V-cycle | V-cycle avg/max it. | basis rel. |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 | 2,048 | 999 | 68.0 ms | 79.2 ms | 577.6 ms | 1,306.5 ms | 9.34 / 10 | 1.11e-11 |
| 12 | 8,192 | 3,919 | 887.2 ms | 541.7 ms | 4,319.7 ms | 5,518.0 ms | 9.66 / 10 | 1.28e-11 |
| 14 | 32,768 | 15,519 | 19,657.5 ms | 3,703.9 ms | 25,251.1 ms | 31,575.9 ms | 9.75 / 10 | 1.44e-11 |

For `h=10`, the V-cycle run used about 100 MiB peak RSS versus 123 MiB for
the exact shifted-LU comparison process. One V-cycle had worst observed
relative residual about `0.142`; GMRES true residual stayed below `1e-10`,
and the primal corrector residual stayed below `7.2e-11`.

The V-cycle is numerically correct and its iteration count is nearly
mesh-independent over this range. It is nevertheless slower than exact
shifted LU by `2.26x, 1.28x, 1.25x`, and much slower than DirectSchur. Each
patch has `m_omega+3` right-hand sides: direct Schur pays one factorization
and cheap triangular solves, whereas the current path performs about ten
V-cycles plus Krylov orthogonalization for every right-hand side. Therefore
`DirectSaddle` remains the public default and `DirectSchur` is the fastest
tested experimental path. The next performance work should prioritize
block/multi-RHS Krylov and reference-patch hierarchy reuse before further smoother tuning.

### Unpreconditioned patch GMRES baseline

`--inverse=none` applies the identity right preconditioner and therefore runs
GMRES directly on each local Helmholtz block. The `H=5,h=10,ell=3,k=4`,
eight-thread Release case used 2,240 right-hand sides and tolerance `1e-10`.

| restart / max-iters | status | corrector | avg/max iteration | basis rel. | primal residual |
|---:|---|---:|---:|---:|---:|
| 30 / 500 | failed | 0.26 s before abort | at least one RHS reached 500 | n/a | n/a |
| 100 / 2,000 | converged | 38.73 s | 259.76 / 385 | 4.66e-11 | 9.99e-11 |
| 500 / 2,000 | converged | 66.89 s | 160.44 / 200 | 4.96e-11 | 9.27e-11 |

The converged solutions are correct, but direct GMRES is not competitive:
shifted-LU and V-cycle preconditioning reduced the same problem to about
`7.35` and `9.34` average iterations. Restart 100 was faster than restart 500
despite more iterations because orthogonalizing a much larger Krylov basis is
expensive. Restart 30 stagnated, so a low iteration cap can hide the severity
of the unpreconditioned problem rather than produce an approximate corrector.

## 19. Coarse LOD Weighted RAS Pilot (M6/S0-S1)

Implemented on 2026-07-20 as an experimental solver for the assembled
Petrov-Galerkin LOD coefficient system. It is separate from the constrained
corrector patch solver and does not change any default.

`HelmholtzCoarseRasPreconditioner` maps every element patch to the union of
its coarse vertices, extracts the principal block
`R_T A_LOD R_T^*`, factors every block once with SparseLU, and injects local
solutions with inverse-overlap weights. The weights form a discrete partition
of unity. Independent local factorizations and applications use OpenMP; each
thread accumulates into a private global vector before deterministic reduction.

The dedicated test checks patch coverage, proper local subdomains, the
partition-of-unity identity, true GMRES residual, and agreement with a global
SparseLU reference. All 14 registered Debug CTests pass.

Release measurements used right GMRES with restart 100, tolerance `1e-10`,
`ell=3`, and eight threads. The first two rows are medians of three runs.

| case | coarse DOF | identity / Jacobi / RAS it. | SparseLU | RAS setup | identity / Jacobi / RAS solve |
|---|---:|---:|---:|---:|---:|
| `H=5,h=10,k=4` | 41 | 30 / 26 / 9 | 0.076 ms | 0.851 ms | 0.831 / 0.665 / 0.514 ms |
| `H=9,h=14,k=16` | 545 | 191 / 182 / 36 | 23.58 ms | 40.34 ms | 143.20 / 127.25 / 143.14 ms |
| `H=11,h=16,k=32` | 2113 | 450 / 461 / 145 | 369.04 ms | 211.02 ms | 1475.55 / 1492.89 / 2185.89 ms |

RAS substantially reduces iteration counts but does not reduce end-to-end
solve time at the tested sizes. At 545 DOFs it only matches identity GMRES
solve time and loses to Jacobi after setup; at 2113 DOFs it is slower even
before setup. Global SparseLU remains decisively faster. The reconstructed
fine-vector differences from SparseLU stayed below `5.0e-10`, so this is a
performance rejection rather than a correctness failure.

Keep this one-level method as a reproducible coarse-smoother candidate. The
fine-space two-level additive/hybrid implementation and its discrete
Helmholtz energy-error benchmark are recorded in Section 20.

## 20. Fine-Space Two-Level LOD Hybrid Schwarz (M6/S2a-S2b)

Implemented on 2026-07-20 for homogeneous Dirichlet boundaries and extended
the same day with an explicit artificial-impedance experiment. The
implementation is experimental and does not change the fine
FEM or LOD default solvers.

`HelmholtzTwoLevelSchwarzPreconditioner` exposes independent coarse, local,
additive, and hybrid actions. The coarse action reuses the model's factored
Petrov-Galerkin LOD system. The local action uses unconstrained Helmholtz
blocks from a dedicated `HelmholtzSchwarzPatchAssembler`, exact local SparseLU,
and inverse-overlap
partition-of-unity weights. It deliberately ignores the corrector constraint
matrix, Schur complement, and element right-hand sides.

For a residual `r`, the hybrid action is evaluated in factored form:

```text
coarse = B0 * r
projected = r - A * coarse
local = Bloc * projected
result = coarse + local - B0 * (A * local)
```

This is `B0 + (I-B0*A) Bloc (I-A*B0)`. The dedicated Schwarz assembler
selects complete-star vertices for homogeneous Dirichlet patches. In impedance
mode it retains artificial-boundary vertices and adds their boundary mass term.
Physical-domain vertices remain present and retain the original Robin term.

The dedicated Debug test checks full fine-DOF coverage, partition of unity,
the Petrov coarse equation, the factorized hybrid formula, true residual, and
the discrete Helmholtz energy error against fine-grid SparseLU. All 15
registered CTests pass. On the small gold case, additive/hybrid required 15/9
iterations, with true residuals `4.78e-11`/`7.33e-12` and energy errors
`3.94e-11`/`4.57e-12`.

Eight-thread Release scale points used direct saddle correctors, exact local
LU, restart 100, and tolerance `1e-10`:

| case | fine DOF | fine SparseLU | setup | identity it./solve | hybrid it./solve | energy rel. |
|---|---:|---:|---:|---:|---:|---:|
| `H=5,h=10,k=4,ell=3` | 1089 | 3.19 ms | 55.1 ms | 307 / 90.3 ms | 7 / 15.7 ms | 2.46e-12 |
| `H=7,h=12,k=8,ell=3` | 4225 | 18.34 ms | 398.3 ms | 1022 / 1067.5 ms | 8 / 112.1 ms | 7.56e-13 |
| `H=9,h=14,k=16,ell=3` | 16641 | 166.8 ms | 2370.7 ms | not run | 9 / 619.8 ms | 1.10e-12 |

Hybrid iterations are promisingly stable, but exact local setup and applying
all local solves every outer iteration remain more expensive than global fine
SparseLU at these sizes. The `H=9` process peaked near 1.91 GiB RSS.

At `H=7,h=12,k=8`, increasing `ell=2,3,4` changed hybrid iterations from
`9,8,7`, while setup rose from `171,398,854 ms` and solve time from
`66,112,170 ms`. Larger patches are therefore not a performance win with the
current exact-local-LU implementation.

S2b now assembles patch volume blocks directly and adds
`-i*k*beta*M_boundary` only on artificial edges. A 16-patch matrix golden
test matches the legacy Dirichlet blocks below `1e-13`; a whole-domain patch
matches the global operator to `3.93e-17`. Small-case impedance
additive/hybrid converge in 19/15 iterations with fine-SparseLU energy errors
below `1e-10`.

Release comparisons rejected impedance as the current default. At
`H=5,h=10,k=4,ell=3`, Dirichlet used 7 iterations and 15.7 ms, while
impedance beta 1 used 21 and 55.0 ms; beta 4 used 14 and 32.0 ms. At
`H=7,h=12,k=8,ell=3`, the corresponding values were 8/141.5 ms,
24/371.6 ms, and 15/204.6 ms. The beta scan 0.25, 0.5, 1, 2, 4 showed monotone
movement toward Dirichlet rather than an optimized Robin window.

The lightweight assembler nevertheless reduced the rerun Dirichlet setup from
55.1 to 42.8 ms at H=5 and from 398.3 to 289.1 ms at H=7. S3 should retain
Dirichlet exact-local LU as the gold path, then evaluate shifted-Laplacian
V-cycles, reusable local hierarchies/factorizations, and a mathematically
consistent ORAS restriction/extension separately.

## 21. Fine-Space Schwarz Shifted Local Solver (M6/S3a)

Implemented on 2026-07-20 as an explicit alternative to exact local SparseLU.
`HelmholtzSchwarzLocalSolver` owns one unconstrained patch solve and supports
`SparseLu` or right-preconditioned `ShiftedGmres` with

```math
P_j=A_j-i\alpha k^2M_j.
```

The Schwarz assembler forms `M_j` from patch elements, not a global principal
submatrix. It does so only for shifted local solvers; the default direct path
neither stores global element mass blocks nor assembles local mass matrices.
The `alpha=0` golden converges in one inner iteration. On a 44-DOF patch,
`alpha=0.2` used six iterations with true residual `8.52e-14` and matched
the direct local solution below `1e-10`.

Tolerance-based inner GMRES makes the Schwarz action variable. The outer
shifted-local path therefore calls the explicit FGMRES API. The Krylov core
already stored each preconditioned vector, so this required an API-level
clarification plus a stateful variable-preconditioner regression test, not a
second algorithm implementation. Direct local Schwarz continues to use
ordinary right GMRES.

Eight-thread Release results with Dirichlet artificial boundaries, hybrid
outer iteration, and tolerance `1e-10` were:

| case | local solver | alpha | setup | outer it. / solve | avg/max inner it. |
|---|---|---:|---:|---:|---:|
| `H=5,h=10,k=4,ell=3` | direct | - | 43.2 ms | 7 / 19.5 ms | 0 / 0 |
| same | shifted | 0 | 47.4 ms | 7 / 49.0 ms | 1 / 1 |
| same | shifted | 0.05 | 53.2 ms | 7 / 84.8 ms | 4.95 / 6 |
| same | shifted | 0.2 | 48.1 ms | 7 / 103.2 ms | 6.61 / 8 |
| `H=7,h=12,k=8,ell=3` | direct | - | 269.3 ms | 8 / 119.9 ms | 0 / 0 |
| same | shifted | 0 | 323.9 ms | 8 / 298.2 ms | 1 / 1 |
| same | shifted | 0.05 | 308.6 ms | 8 / 686.0 ms | 5.22 / 6 |
| same | shifted | 0.2 | 308.5 ms | 8 / 691.6 ms | 7.11 / 8 |

All shifted runs retained the same outer iteration count and fine-SparseLU
accuracy as direct local solves. The cost model explains the loss: the H=7
hybrid run applies 256 local solvers eight times, so even one Krylov step means
2,048 local triangular solves plus matrix products and orthogonalization.
Exact shifted SparseLU is therefore a correctness reference, not an
optimization.

S3a is complete and does not change defaults. S3b should test a fixed,
linear geometric V-cycle and reuse patch hierarchy/workspace. Only that path
can remove the per-patch shifted factorization; if it remains slower, further
inner tolerance tuning is not justified. A true ORAS restriction/extension
remains a separate experiment.


## 22. Fine-Space Schwarz Geometric V-Cycle (M6/S3b)

Implemented on 2026-07-20 as the second explicit shifted-local experiment.
`HelmholtzDirichletPatchHierarchyBuilder` restricts the model's nested NVB
node prolongations to complete-star Dirichlet patch spaces. It validates that
the final local ordering exactly matches the Schwarz block and constructs the
hierarchy only when `shifted_inverse=GeometricVcycle`; direct local LU and
shifted-LU setup remain unchanged.

For each patch, the fixed linear V-cycle approximates

```math
P_j^{-1},\qquad
P_j=A_j-\mathrm{i}\alpha k^2M_j,
```

using Galerkin coarse operators, damped Jacobi pre/post smoothing, and
SparseLU only on the selected coarsest local level. Right-preconditioned
local GMRES still solves the unshifted equation `A_j z_j=r_j` and stops on
its true residual. Artificial impedance boundaries are rejected for this
hierarchy because the current level construction represents homogeneous
Dirichlet patch spaces only.

The focused tests check hierarchy dimensions, fixed-V-cycle linearity,
agreement with a direct local solve, local true residual, the impedance guard,
and an end-to-end hybrid solve against fine SparseLU. On the 44-DOF local
golden, shifted-LU and V-cycle used 6 and 10 inner iterations, respectively,
with residuals below `9e-14`. The hybrid test retained 9 outer iterations
and `4.57e-12` relative discrete energy error. All 16 Debug CTests pass.

Eight-thread Release results used Dirichlet patches, `alpha=0.2`, local and
outer tolerance `1e-10`, and two pre/post Jacobi steps:

| case | local inverse | setup | outer it. / solve | avg/max inner it. | V-cycle levels / coarse DOF |
|---|---|---:|---:|---:|---:|
| `H=5,h=10,k=4,ell=3` | direct LU | 50.1 ms | 7 / 19.1 ms | 0 / 0 | - |
| same | shifted LU | 49.7 ms | 7 / 105.0 ms | 6.61 / 8 | - |
| same | V-cycle | 42.0 ms | 7 / 236.7 ms | 8.62 / 9 | 4 / 127 |
| `H=7,h=12,k=8,ell=3` | direct LU | 289.1 ms | 8 / 124.6 ms | 0 / 0 | - |
| same | V-cycle | 215.8 ms | 8 / 1601.4 ms | 9.30 / 11 | 4 / 127 |

All paths retained the same outer iteration count, true residual, and
fine-SparseLU accuracy. The V-cycle removes most full-size shifted
factorizations and cuts setup by about 16% at H=5 and 25% at H=7, but every
hybrid application performs all patch solves and each local Krylov step now
contains a multi-level sparse cycle. Solve time is therefore about 12.4 to
12.9 times direct LU in the tested cases.

S3b is complete as a correct modular experiment and is rejected as a runtime
default. Do not spend the next stage on looser inner tolerances or isolated
Jacobi tuning. A plausible next stage must amortize local work across the many
right-hand sides and translated patch types, for example block Krylov,
reference-patch hierarchy/operator reuse, or a separately derived ORAS
restriction/extension.


## 23. Fine-Space Schwarz Robustness Baseline (M6/S4a)

Implemented on 2026-07-20 without changing any solver default.
`bench_helmholtz_two_level_schwarz` now accepts
`--source=gaussian|manufactured`. Manufactured mode reuses
`make_polynomial_plane_wave_solution` and reports the continuous energy and
L2 errors, their exact-solution-normalized values, the LOD Petrov residual,
corrector/constraint residuals, and the Petrov residual of every outer Krylov
solution. It also prints physical coarse/fine meshwidths, `kH`, and a
`steady_clock` total; the latter is the trusted total when operating-system
wall time is disturbed by a clock adjustment.

The executable `scripts/run_helmholtz_schwarz_scale.sh` builds the Release
target, accepts resumable `H:h:k` cases and ell lists, records one log and
resource file per case, and writes a machine-readable `summary.csv`.
Defaults are the fixed-`kH=1` sequence
`5:10:4 7:12:8 9:14:16`, `ell=3`, manufactured source, Dirichlet
artificial boundaries, direct local LU, and hybrid outer GMRES.

The eight-thread S4a study covered `k=4,8,16`, fine gaps 4, 5, and 6.
For every gap, coarse patches increased from 64 to 1,024 and the hybrid
iteration count stayed in the narrow range 7 to 9:

| fine gap | k=4 | k=8 | k=16 |
|---:|---:|---:|---:|
| 4 | 7 | 8 | 9 |
| 5 | 8 | 8 | 9 |
| 6 | 8 | 8 | 9 |

All runs converged with fine-system true residual below `1e-10`; outer
Petrov residuals stayed below `1.6e-12`, and the hybrid/fine-SparseLU
relative discrete energy difference stayed below `3.6e-12`. Corrector and
constraint residuals remained below `2.2e-15` and `4.8e-16`.
The normalized fine-FEM manufactured energy error decreased along every
fixed-gap sequence:

| fine gap | k=4 | k=8 | k=16 |
|---:|---:|---:|---:|
| 4 | 9.99e-2 | 6.08e-2 | 4.92e-2 |
| 5 | 6.87e-2 | 3.93e-2 | 2.85e-2 |
| 6 | 5.00e-2 | 3.02e-2 | 2.37e-2 |

At the representative `H=7,h=12,k=8` point, `ell=2,3,4` required
9, 8, and 7 outer iterations. Larger oversampling reduced one iteration per
layer but increased setup from 138.5 to 289.5 to 531.9 ms and solve from
69.4 to 117.2 to 233.8 ms.

This establishes iteration robustness, not runtime competitiveness. The
largest local run `H=9,h=15,k=16,ell=3` used 33,025 fine DOFs, 1,024
patches, and about 4.54 GB peak RSS. Local setup plus hybrid solve was about
5.73 s, versus 0.447 s for the fine SparseLU reference; complete benchmark
time was 15.53 s. The default therefore remains global fine SparseLU for a
single fine solve and direct local LU inside the experimental hybrid action.

S4a established the fixed-`kH` hybrid baseline. S4b now compares identity,
local, additive, and hybrid actions on the same manufactured dataset; the
results and remaining S4c work are recorded below. M5 remains open because
no corrector iterative policy has crossed the DirectSaddle runtime threshold.

## 24. Fine-Space Schwarz Solver Comparison (M6/S4b)

Completed on 2026-07-20 without changing solver defaults. The resumable scale
driver now stores `source`, `solver`, `boundary`, `local_solver`, and
`local_inverse` in every CSV row and uses the same fields in its resume key.
Different outer actions can therefore share one result directory without
silently skipping or mixing cases.

Using the S4a manufactured problem, fixed `kH=1`, `ell=3`, Dirichlet
artificial boundaries, direct local LU, and eight threads gave:

| case | identity | local | additive | hybrid |
|---|---:|---:|---:|---:|
| `H=5,h=10,k=4` | 466 / 66.93 ms | 11 / 21.36 ms | 15 / 42.72 ms | 8 / 27.23 ms |
| `H=7,h=12,k=8` | 1285 / 679.04 ms | 20 / 219.00 ms | 18 / 234.41 ms | 8 / 124.27 ms |
| `H=9,h=14,k=16` | 2750 / 7533.72 ms | 51 / 2827.82 ms | 22 / 1397.20 ms | 9 / 627.71 ms |

Each entry is outer iterations / outer solve time. All 12 runs converged;
true residuals were below `1e-10`, Petrov residuals below `1.3e-10`, and
relative discrete energy errors against fine SparseLU below `2.5e-10`.
Manufactured-solution errors agree across all actions.

The comparison separates the roles of the two levels. Local correction alone
eventually degrades, additive remains useful, and the multiplicative hybrid
coarse/local composition is the only tested action staying at 8-9 iterations.
It is nevertheless not the runtime default: at `H=9,h=14,k=16`, hybrid setup
plus solve is about 2.22 s versus 0.151 s for fine SparseLU. The identity rows
also pay Schwarz setup because the benchmark constructs common diagnostics;
use `outer_solve_ms` for Krylov-only comparisons and setup plus solve for an
algorithm-level comparison.

Raw S4b logs and the configuration-aware CSV are in
`results/helmholtz_schwarz_s4b/`. The local higher-wave-number S4c extension
is recorded in the next section; its equal-accuracy server continuation is
deferred. M5 remains open because neither iterative path has crossed the
relevant direct-solver runtime.

## 25. Local High-Wavenumber Schwarz Study (M6/S4c)

Completed locally on 2026-07-20 after deferring the large-memory server run.
All cases use the manufactured source, fixed `kH=1`, Dirichlet artificial
boundaries, direct local LU, eight threads, restart 50, and outer tolerance
`1e-10`.

At `H=11,h=15,k=32,ell=3` (33,025 fine DOFs and 4,096 patches),
identity/local/additive/hybrid required `2869/385/29/9` iterations. Their
outer solve times were `16.93/39.18/3.64/1.45 s`. All true residuals were
below `1e-10`, and relative discrete energy errors against fine SparseLU
were below `1.1e-10`. One-level local is therefore not merely less robust:
its expensive all-patch actions make it about 2.3 times slower than identity.

The local memory limit still permits a controlled `k=64` probe by reducing
the fine gap from four to three. The hybrid oversampling comparison is:

| case | ell | iterations | setup | solve | peak RSS | LOD exact energy rel. |
|---|---:|---:|---:|---:|---:|---:|
| `k=32,H=11,h=15` | 2 | 9 | 1.31 s | 0.52 s | 1.66 GiB | 5.35e-2 |
| same | 3 | 9 | 3.43 s | 1.45 s | 3.64 GiB | 5.27e-2 |
| `k=64,H=13,h=16` | 2 | 8 | 3.02 s | 1.32 s | 3.91 GiB | 8.94e-2 |
| same | 3 | 9 | 6.87 s | 4.65 s | 6.93 GiB | 8.71e-2 |

The `k=64,ell=3` additive action took 34 iterations. Thus hybrid iteration
robustness extends through `k=64` in the tested range, while `ell=2` is much
cheaper than `ell=3`. This is empirical and does not remove the theoretical
need for oversampling that grows logarithmically with wavenumber.

The `k=64` fine-FEM normalized exact energy error is about `8.69e-2`, so the
gap-three case is an iteration/memory probe rather than an equal-accuracy
runtime comparison. Even the faster ell-two hybrid setup plus solve is
about 4.34 s versus 1.73 s for fine SparseLU. A gap-four `k=64,ell=3` run
was not attempted because extrapolated storage exceeds the 12 GiB WSL
budget. Raw data are in `results/helmholtz_schwarz_s4c_local/`.

S4c-local is complete. S4d now makes both ORAS extension choices explicit;
the following section records their validation and rejection as runtime
defaults. Block/multi-RHS and reference-patch reuse move to S4e.

## 26. Explicit ORAS Extensions (M6/S4d)

Completed on 2026-07-20 without changing the default action. The local
preconditioner is now represented as

```text
B_local = sum_j R_j^T D_j B_j^{-1} R_j,
sum_j R_j^T D_j R_j = I.
```

`WeightedOverlap` uses inverse multiplicity on every overlapping copy.
`RestrictedCore` assigns each fine DOF to one deterministic coarse-element
core, solves on the full overlapping patch, and injects only owner-core
values. With an impedance artificial boundary, the first path is weighted
ORAS and the second is a Boolean-core RAS/ORAS variant. This corrects the
earlier overly narrow statement that an impedance local solve was not ORAS
unless it used Boolean restriction.

The patch assembler now exposes core DOFs. The preconditioner stores
per-subdomain extension weights and reports min/max owned DOFs in addition
to local patch sizes. Both extensions satisfy partition of unity to roundoff.
The scale driver records `extension`, `impedance_beta`, and owned sizes in
the CSV identity and includes them in collision-free log names.

The Debug golden case gave:

| action | additive | hybrid |
|---|---:|---:|
| Dirichlet weighted | 15 | 9 |
| weighted ORAS, beta=1 | 19 | 15 |
| Boolean-core ORAS, beta=1 | 30 | 23 |

All converged solutions agree with fine SparseLU within the established
energy and true-residual tolerances. Full Debug CTest passes 16/16.

Release results with eight threads and `ell=3` were:

| case | method | iterations | solve |
|---|---|---:|---:|
| `k=8,H=7,h=12` | Dirichlet weighted hybrid | 8 | 0.108 s |
| same | weighted ORAS hybrid, beta=4 | 15 | 0.247 s |
| same | weighted ORAS additive, beta=4 | 22 | 0.321 s |
| same | RestrictedCore Dirichlet hybrid | 21 | 0.302 s |
| `k=16,H=9,h=14` | Dirichlet weighted hybrid | 9 | 0.674 s |
| same | weighted ORAS hybrid, beta=4 | 16 | 1.349 s |
| same | weighted ORAS additive, beta=4 | 25 | 1.693 s |
| same | RestrictedCore Dirichlet hybrid | 41 | 2.592 s |

At `k=8`, RestrictedCore ORAS with beta `1,2` failed at 4,000 iterations;
beta `4,8,16,32,64,128` required `276,74,37,25,22,21` iterations.
The beta-to-infinity limit matches the 21-step RestrictedCore Dirichlet
case. At `k=16`, RestrictedCore beta 16/64 required 74/42 iterations,
while RestrictedCore Dirichlet required 41.

The Boolean core owns only 9-25 DOFs per subdomain versus 509-1,561 local
patch DOFs. This sharp extension loses the smooth overlap exchange supplied
by multiplicity weights. Weighted ORAS is better than Boolean-core ORAS but
still loses to Dirichlet weighted hybrid on both scales. Both remain explicit
research switches; neither changes the runtime default.

Raw logs are in `results/helmholtz_schwarz_s4d_oras/`. S4e moves to
block/multi-RHS local solves and translated-patch symbolic/hierarchy reuse.

## 27. Identical-Matrix Schwarz Factorization Reuse (M6/S4e)

Completed on 2026-07-21 without changing the default. A Schwarz application
has one right-hand side per patch, so block GMRES inside one patch has no
natural batch. On uniform meshes with constant coefficients, however, many
translated patches have exactly the same local matrix. They are grouped as

```text
G_q = {j : A_j == A_q}
F_q = [R_j r] for j in G_q
X_q = SparseLU(A_q).solve(F_q)
```

The sparse hash covers dimensions, compressed outer/inner indices, and every
real and imaginary coefficient. A hash hit is followed by full exact matrix
comparison. This is a correctness boundary, not a tolerance-based geometric
guess. One `SparseLU` is owned per group, one dense multi-RHS solve handles
all group columns, and groups remain OpenMP-parallel. Consequently a shared
Eigen factorization is never called concurrently.

`HelmholtzSchwarzFactorizationReuse::IdenticalMatrix` is an explicit direct-LU
switch; `None` remains the default, and shifted GMRES/V-cycle rejects reuse.
Diagnostics expose solver group count, saved factorizations, and maximum group
size. Tests compare block columns against scalar SparseLU, compare reused and
unreused local actions, and validate the final fine-SparseLU solution.

Eight-thread Release measurements are in
`results/helmholtz_schwarz_s4e_reuse/summary.csv`:

| case | mode | groups | setup | hybrid solve | iterations | peak RSS |
|---|---|---:|---:|---:|---:|---:|
| `k=8,H=5,h=10` | none | 64 | 37.0 ms | 17.1 ms | 8 | 116,540 KiB |
| same | identical | 36 | 41.3 ms | 16.0 ms | 8 | 92,768 KiB |
| `k=8,H=7,h=12` | none median | 256 | 237.8 ms | 103.1 ms | 8 | 500,000 KiB |
| same | identical median | 75 | 194.6 ms | 87.7 ms | 8 | 289,552 KiB |
| `k=16,H=9,h=14` | none | 1024 | 1.520 s | 0.781 s | 9 | 2,010,232 KiB |
| same | identical | 166 | 1.086 s | 0.429 s | 9 | 823,180 KiB |

The `H=7` row is the median of three runs. True residuals and energy errors
are unchanged to roundoff. At `H=9`, 858 factorizations are removed; setup,
solve, and peak memory fall by about 29%, 45%, and 59%, respectively.
`H=3,h=7` remains slightly slower because grouping overhead dominates tiny
LU factors, so automatic/default reuse would be premature. The complete
Schwarz method still does not cross fine SparseLU wall time; M5 remains open.

## 28. Helmholtz Local Inverse Inequality on Graded NVB Meshes

Implemented and run locally on 2026-07-21. The experiment follows
`HELMHOLTZ_LOCAL_INVERSE_NONQUASI_PLAN.md` and deliberately distinguishes the
same-element quotient from the more robust one-layer-patch denominator. It
does not turn the numerical observations into a proof.

### Implementation and nesting boundary

`bench_helmholtz_local_inverse` starts from the compatible unit-square NVB
mesh and, at every iteration, marks only a deterministic subset of the
currently finest leaves. The resulting triangulation is globally conforming
and shape regular but increasingly non-quasi-uniform. The benchmark computes
the coarse P1, corrected trial, and corrected test spaces for

```text
Q_T = H_T sup ||grad G a||_T / ||G a||_T
```

and for the one-layer-patch denominator

```text
Q_T_patch = H_T sup ||grad G a||_T / ||G a||_{omega_T^1}.
```

The complex local matrices are `G^* S_T G` and `G^* M_T G`; they use
`adjoint()`, are checked for Hermitian defects, and are reduced to the positive
mass eigenspace before the maximum eigenvalue is computed. Trial and test are
kept as separate output rows even though they are conjugates for the current
real coefficients.

The fixed-master-fine-space invariant is checked independently on every
coarse mesh. The reported nesting residual is the maximum of:

1. prolongation errors for the constant, `x`, and `y` coarse P1 functions;
2. unique fine-element parent mapping and parent/child area conservation;
3. CG/DG prolongation consistency;
4. `||I_H P_H-I||`.

The fine mesh is also compared with the canonical global-NVB mesh at the
requested fixed level. Across all three primary mesh families the largest
nesting residual was `7.70e-15`, every canonical fine-mesh comparison passed,
and the maximum neighboring-element diameter ratio remained `sqrt(2)`. Thus
the experiment verifies `V_H subset V_h` algebraically rather than inferring
it from matching coordinates.

### Fixed-fine-space grading scan

The primary run used `k=2`, initial coarse level 3, fixed fine level 12,
`ell=3`, six local refinement steps, and a 25% deterministic fraction of the
finest leaves nearest `(0.25,0.25)`. `Gamma=H_max/H_min` increased from 1 to 8
while the local resolution ratio `q_max` increased from `0.0442` to `0.3536`.

| step | level range | Gamma | q_max | coarse P1 | coarse with wrong global H | trial element | trial patch1 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 3:3 | 1.000 | 0.0442 | 8.4853 | 8.4853 | 64.2240 | 3.4018 |
| 1 | 3:4 | 1.414 | 0.0625 | 8.4853 | 12.0000 | 99.6996 | 3.9427 |
| 2 | 3:5 | 2.000 | 0.0884 | 8.4853 | 16.9706 | 97.5129 | 6.0012 |
| 3 | 3:6 | 2.828 | 0.1250 | 8.4853 | 24.0000 | 101.7924 | 10.9635 |
| 4 | 3:7 | 4.000 | 0.1768 | 8.4853 | 33.9411 | 101.6788 | 12.7162 |
| 5 | 3:8 | 5.657 | 0.2500 | 8.4853 | 48.0000 | 101.4298 | 12.3763 |
| 6 | 3:9 | 8.000 | 0.3536 | 8.4853 | 67.8823 | 101.4366 | 12.4135 |

The `C=0` value remains exactly `8.485281` when each element uses its own
`H_T`. The deliberately wrong global-`H_max` control grows by the full grading
ratio. This validates the local diameter, child map, and generalized
eigenvalue pipeline. The trial/test element maximum reaches an apparent
fixed-`h` plateau near 101; its argmax stays on a level-3 element rather than
inside the finest core.

The final values for `fraction`, `single-chain`, and `boundary-chain` were:

| family | Gamma | trial element | trial patch1 |
|---|---:|---:|---:|
| fraction | 8 | 101.4366 | 12.4135 |
| single-chain | 8 | 101.4366 | 12.4135 |
| boundary-chain | 8 | 103.3774 | 9.8501 |

The fraction `theta=0.10,0.25,0.50` endpoint patch values were
`12.4135,12.4135,12.7699`. Trial and test values agree to roundoff. Petrov,
corrector, and constraint residuals over the primary data were at most
`1.25e-14`, `5.91e-15`, and `1.18e-15`; the maximum local Hermitian defect and
eigen residual were `5.16e-16` and `1.41e-13`.

### Fine-space, mass-rank, oversampling, and wave-number probes

Fixing the final fraction-family coarse mesh (levels 3:9) and refining only
the master fine space gave:

| fine level | q_max | trial element | trial patch1 | max nesting residual |
|---:|---:|---:|---:|---:|
| 11 | 0.5000 | 97.3016 | 11.9796 | 2.03e-15 |
| 12 | 0.3536 | 101.4366 | 12.4135 | 7.37e-15 |
| 13 | 0.2500 | 124.3089 | 12.4487 | 4.78e-15 |
| 14 | 0.1768 | 132.5017 | 12.5692 | 3.76e-14 |

The patch-denominator value is resolved to a narrow range, whereas the
same-element value is not converged. Its local mass condition number reaches
roughly `1e11-1e12`. Changing the positive-mass threshold from `1e-10` to
`1e-12` to `1e-14` changes the final element value from `82.88` to `101.44`
to `116.85`; the corresponding patch values are `10.62`, `12.41`, and
`12.44`. The `1e-14` run is retained as a diagnostic rather than a passing
gate because its energy-identity roundoff exceeds the registered tolerance.

For `ell=2,3,4`, the final element values were `93.50,101.44,97.28` and the
patch values were `9.85,12.41,11.21`; no systematic oversampling divergence
was observed. Resolution-compatible wave-number endpoints were:

| k | initial/final levels | k H_max | Gamma | initial/final element | initial/final patch1 |
|---:|---:|---:|---:|---:|---:|
| 1 | 3:3 -> 3:9 | 0.5 | 8 | 65.18 / 101.75 | 3.39 / 12.58 |
| 2 | 3:3 -> 3:9 | 1.0 | 8 | 64.22 / 101.44 | 3.40 / 12.41 |
| 4 | 5:5 -> 5:9 | 1.0 | 4 | 71.81 / 76.77 | 8.91 / 12.89 |

Fine level 14 completed in 1:58 with peak RSS `5,971,068 KiB`. The planned
fine-level-15 endpoint was killed by signal 9 after reaching `11,645,504 KiB`;
it is a server continuation, not a numerical failure.

### Feedback refinement at the observed maximum

The initial feedback grading experiment selected the next marked set from
the inverse-constant maximizer rather than from a prescribed geometric chain.
For iteration `j`, its two feedback rules were

```text
T*_j = argmax_T Q_T             (argmax-element),
T*_j = argmax_T Q_T_patch1      (legacy argmax-patch run),
M_j  = {T*_j} or N_1(T*_j),
T_{j+1} = NVB-close(refine(T_j, M_j)).
```

Ties are broken by the stable element index. `N_1(T*_j)` is the one-layer
vertex patch; the
=0` runs mark only `T*_j`. Unlike the earlier geometric
families, the feedback target is allowed to be at any current level. This is
essential because the observed maximum need not lie on a finest element.

All runs used `k=2`, initial level 3, fixed fine level 12, `ell=3`, and six
feedback steps. The endpoint comparison is:

| feedback quotient | neighbor layers | elements | level range | Gamma | trial element | trial patch1 |
|---|---:|---:|---:|---:|---:|---:|
| element | 0 | 28 | 3:5 | 2.000 | 83.8179 | 4.1160 |
| element | 1 | 105 | 5:7 | 2.000 | 76.9891 | 12.3634 |
| patch1 | 0 | 25 | 3:8 | 5.657 | 88.0053 | 8.6512 |
| patch1 | 1 | 70 | 3:7 | 4.000 | 94.0259 | 12.3047 |

The `argmax-patch,n=0` level ranges evolve as
`3:3, 3:4, 3:5, 3:6, 3:6, 3:7, 3:8`; it therefore gives the strongest
feedback-driven non-quasi-uniform grid with the fewest extra elements. Its
patch maxima are `3.4018, 3.8849, 9.3743, 8.7185, 6.3599, 8.9497, 8.6512`:
after the initial jump they fluctuate but do not grow monotonically with
`Gamma`. Marking a one-layer neighborhood is substantially less local. In
the element-feedback case it eventually raises the minimum level from 3 to 5,
so `Gamma` falls from its intermediate value 4 back to 2. This option is useful
as a robustness comparison, but not as the primary stress test for global
non-quasi-uniformity.

Every feedback row retained the canonical fixed fine mesh. The maximum
nesting residual was `7.99e-15`, the maximum Petrov residual was `1.28e-14`,
and the maximum constraint residual was `1.11e-15`. Hence feedback selection
and NVB closure preserve the required exact discrete inclusion
`V_H subset V_h` to roundoff.

### Oversampling-matched denominator and deep fine space

The one-layer denominator above does not match an `ell=3` localized
corrector. The benchmark now reports all three diagnostics

```text
Q_T,0   = H_T sup ||grad v||_T / ||v||_T,
Q_T,1   = H_T sup ||grad v||_T / ||v||_{omega_T^1},
Q_T,ell = H_T sup ||grad v||_T / ||v||_{omega_T^ell}.
```

`argmax-patch` now means `argmax_T Q_T,ell`; `patch1` remains in the output
only as a narrower-denominator diagnostic. Thus the denominator region uses
the same coarse-element patch construction and radius as the corrector
oversampling region.

The master fine mesh is global and fixed across coarse iterations. For a
coarse element `T`, the measured separation is

```text
q_T = max_{t subset T} h_t / H_T,    q_max = max_T q_T.
```

Under this NVB level convention a level difference of two approximately
halves the diameter, so selecting `L_h` is based on the finest coarse level:
`q_max` must be small for the level-`L_max` elements, not merely for the
initial coarse mesh. To isolate the `h` effect, a single fixed level-`3:5`
coarse mesh (21 elements) was evaluated at `L_h=11,12,13`:

| fine level | q_max | trial element | trial patch1 | trial patch3 | nesting residual |
|---:|---:|---:|---:|---:|---:|
| 11 | 0.125000 | 94.1403 | 5.0627 | 3.5289 | 1.89e-15 |
| 12 | 0.088388 | 94.2433 | 4.8820 | 3.5181 | 4.99e-15 |
| 13 | 0.062500 | 184.1021 | 4.9155 | 3.5124 | 5.62e-15 |

The oversampling-matched `patch3` maximum changes by only `0.47%` from
`L_h=11` to `L_h=13` and by `0.16%` over the last refinement. It is therefore
resolved at `h/H_T <= 1/16` on this mesh. In contrast, the same-element value
doubles at `L_h=13`, reinforcing the earlier warning that the strict
same-element quotient is sensitive to fine-scale near-null mass directions.

An `L_h=13`, three-step `argmax-patch`, target-only feedback run also passed.
Its final grid had levels `3:5`, `q_max=0.0625`, and trial maxima
`184.5663`, `11.1407`, and `3.5414` for element, patch1, and patch3
denominators. The patch3 trajectory was
`3.1806, 3.5243, 3.1863, 3.5414`; it stayed bounded while driving the mesh.
Every state used the canonical fixed fine mesh and the largest nesting
residual was `3.82e-15`.

The attempted `L_h=14` feedback run was terminated by signal 9 at
`11,668,680 KiB` peak RSS after 4:18. It is retained as a documented memory
limit. Since the successful `L_h=13` fixed-grid point already has
`q_max=1/16` and the patch3 values have reached a tight plateau, this failure
does not obstruct the local `h`-resolution conclusion.

### Conclusion

The strict same-element statement is **not numerically established** by this
experiment. Its value depends strongly on fine-space resolution and on how
near-null local mass directions are retained. The current data therefore do
not justify a grading-independent same-element inverse constant.

The oversampling-matched patch denominator is especially stable under the
fixed-grid deep-`h` scan and remains near `3.5` there. The earlier one-layer
form remains useful as a stricter diagnostic, but it should not be described
as support-matched when `ell>1`. The correct statement is therefore: the
matched patch form has encouraging resolved behavior, but the present local
data are still insufficient to claim a uniform grading-independent theorem.

`helmholtz_local_inverse_smoke` and its reduced server-mode counterpart are
registered in CTest. The full Release suite passes `18/18`. Raw summaries,
per-element spectra, meshes, timing logs, and
all auxiliary `h`, threshold, `ell`, and `k` scans are in
`results/helmholtz_local_inverse/`.

### EPYC server execution and performance path

The local-inverse postprocessing now supports OpenMP dynamic scheduling over
coarse elements and reduced server modes. `--basis=trial` avoids recomputing
the conjugate test and coarse diagnostic spectra after a full calibration;
`--denominators=matched` computes only `omega_T^ell`, while
`element-matched` retains the strict element diagnostic as well. Defaults
remain `all/all` for backward compatibility.

On the local WSL `L_h=12` fixed-grid calibration, full analysis took
`7.86 s` with one thread and `2.84 s` with eight threads. With eight threads,
trial/matched mode took `2.40 s`; the inverse-analysis portion fell from
about `154 ms` to `27 ms`, and the matched quotient was identical. A test of
the server script itself gave `7.44,2.81,1.94 s` for 1,4,8 threads.

`scripts/run_helmholtz_local_inverse_server.sh` is the reproducible entry
point for the AMD EPYC 9554 / 377 GiB run. It performs a `8,16,32,64` thread
pilot, fixed-coarse-grid `L_h=14,15,16` scan, and `patch3`-argmax feedback run
at `L_h=16`. OpenMP threads are pinned to physical-core places with spread
placement. Cases run sequentially and write success-only `.done` markers so
long jobs can resume safely. The complete commands, acceptance criteria,
memory escalation rule, and result-return list are in
`HELMHOLTZ_LOCAL_INVERSE_SERVER_RUNBOOK.md`.

The first server `L_h=16` feedback attempt exposed a configuration issue:
the command used `--denominators=element-matched`, so the known ill-conditioned
same-element mass problem could trip the hard eigenproblem check before the
support-matched `patch3` experiment completed. Deep-`h` server `hscan` and
`feedback` modes now use `--denominators=matched`; the small full calibration
continues to test `all/all`. Failure messages now identify the iteration,
basis, denominator, Hermitian defect, eigen residual, energy-identity error,
and maximum retained mass condition number instead of reporting a generic
local-eigenproblem failure.

### EPYC server results

The complete AMD EPYC 9554 / 377 GiB run was returned on 2026-07-22 in
`results/helmholtz_local_inverse_server/`. Every planned case has a `.done`
marker and exit status zero. For the fixed level-`3:8` coarse grid with
`Gamma=H_max/H_min=4 sqrt(2)`, define

```text
C_inv,3(h) = max_T H_T sup_{0 != v in V_H,3^ms}
             ||grad v||_T / ||v||_{omega_T^3},
q_max      = max_T max_{t subset T} h_t/H_T.
```

The deep fine-space scan is:

| fine level | q_max | same-element | patch3 | nesting residual | wall time | peak RSS |
|---:|---:|---:|---:|---:|---:|---:|
| 14 | 0.125000 | 132.0481 | 3.2019028413 | 3.80e-14 | 28.60 s | 18.2 GiB |
| 15 | 0.088388 | 143.3774 | 3.2017193578 | 2.96e-14 | 2:58 | 78.4 GiB |
| 16 | 0.062500 | 162.1273 | 3.2013850546 | 1.38e-13 | 16:08 | 180.8 GiB |

Thus

```text
|C_inv,3(L15)-C_inv,3(L14)| / C_inv,3(L14) = 5.73e-5,
|C_inv,3(L16)-C_inv,3(L15)| / C_inv,3(L15) = 1.044e-4,
```

and the total `L_h=14 -> 16` change is only `1.617e-4` (`0.01617%`).
The preregistered `2%` fine-space plateau gate is passed by more than two
orders of magnitude. At `q_max=1/16`, the resolved numerical statement is

```text
||grad v||_{L2(T)} <= 3.202 H_T^{-1} ||v||_{L2(omega_T^3)}
```

for every tested coarse element and discrete trial LOD function. The number
`3.202` is an experimental upper bound for this mesh family, not a proved
uniform constant. In contrast, the same-element diagnostic increases by
`22.78%` from level 14 to 16 and its retained mass condition reaches
`5.57e11`; it remains unsuitable as the main inequality.

The successful `L_h=16` matched-denominator feedback trajectory is

```text
C_inv,3^(j) = 3.178899, 3.522896, 3.184477, 3.539848,
              3.184809, 3.539944, 3.189363,
```

so `3.1788 <= C_inv,3^(j) <= 3.5400` for `j=0,...,6`. The maximizer moves
among symmetry-related physical-boundary elements. Consequently this feedback
family ends at levels `3:5` and `Gamma=2`; it supports boundedness under the
tested feedback rule but is not the strong global-grading stress test. The
fixed `3:8` scan above supplies that evidence.

Every server row has `fixed_fine_mesh=1`; the largest nesting, primal
corrector, and constraint residuals are `1.38e-13`, `1.32e-13`, and
`6.70e-15`. The neighboring diameter ratio remains `sqrt(2)`, and patch3
Hermitian/eigen residuals remain at roundoff. The thread pilot gives
`1.14,0.84,0.57,0.59 s` for `8,16,32,64` threads, confirming 32 threads as
the best tested setting. The feedback run takes `1:00:22` and peaks at
`237.1 GiB`; because level 16 is already converged and exceeds the runbook's
`120 GiB` escalation gate, level 17 should not be attempted on this server.

The final numerical conclusion is conditional but sharp: the data strongly
support an `h`-resolved oversampling-matched patch inverse inequality with
`C_inv,3` of order `3.2-3.54` on the tested graded NVB families. They do not
establish a theorem uniform over all admissible graded meshes, and they do not
support replacing `omega_T^3` by `T` in the denominator.

## 29. Continuous hp Helmholtz Fine Space and Single-Patch Corrector

The first three implementation steps of
`HELMHOLTZ_HP_CORRECTOR_CONVERGENCE_PLAN.md` were completed on 2026-07-24.
The implementation deliberately stops before the global hp
Petrov-Galerkin LOD model.

The new `HpTriSpace` supports continuous triangular `P1`, `P2`, and `P3`
spaces on globally conforming NVB meshes. Global numbering contains vertex
DOFs, consistently oriented shared-edge DOFs, and element-interior DOFs.
A nodal Vandermonde basis is used only for these low orders. Tensor
Gauss-Legendre quadrature under a Duffy map integrates triangle terms, and
one-dimensional Gauss quadrature integrates physical Robin edges.

`HelmholtzHpOperators` assembles

```text
A_hp = K_hp - kappa^2 M_hp - i kappa R_hp.
```

It also provides fine-hp load assembly, sparse-LU solution, and error
integration. For the homogeneous-Robin manufactured solution

```text
u(x,y) = phi(x) phi(y) exp(i kappa x),
phi(t) = 16 t^2 (1-t)^2,
```

the small three-grid regression observed energy rates
`0.743, 1.972, 2.738` for `p=1,2,3`. The `p=1` sequence is still
pre-asymptotic; these numbers are regression evidence, not the planned
final convergence study.

The hp interpolation path constructs the coarse P1 nodal injection
`P_Hhp` and the mixed-moment quasi-interpolation
`I_H=E_H Pi_H^dg`. It checks `I_H P_Hhp=I` during construction. At `p=1`,
the hp stiffness, mass, Robin matrix, injection, and quasi-interpolation
agree with the existing P1 implementation to the test tolerances.

`HelmholtzHpPatchAssembler` then restricts the global hp operators to one
patch. A DOF is retained exactly when all globally incident fine elements
belong to the patch. This removes every high-order artificial-boundary DOF
while retaining physical Robin-boundary DOFs. The existing
`HelmholtzPatchSystem` and `DirectSaddle` solver are reused. Tests cover an
interior target, a physical-boundary target, and a corner target; normalized
primal, constraint, and adjoint residuals are below `2e-10`.

Validation:

```text
cmake --build build --target test_helmholtz_hp_fem test_helmholtz_hp_patch -j 8
./build/tests/test_helmholtz_hp_fem
./build/tests/test_helmholtz_hp_patch
ctest --test-dir build --output-on-failure
```

All 20 registered CTest cases passed at the HP3 boundary.

## 30. hp Petrov-Galerkin LOD and H-Convergence Calibration

HP4-HP6 were advanced on 2026-07-24. `HelmholtzHpLodModel` now owns the
continuous `P1/P2/P3` fine space, hp interpolation, direct-saddle element
correctors, two-sided trial/test bases, coarse Petrov-Galerkin factorization,
and fine reference factorization. Repeated right-hand sides reuse all these
objects. The `p=1` golden test compares the complete hp path with the existing
P1 model; operator, basis, coarse-system, and solution differences are below
`3e-9`.

The fine-hp calibration at `k=4`, `L_h=4,6,8,10,12` gives:

| p | final energy order | final L2 order |
|---:|---:|---:|
| 1 | 1.0001 | 1.9951 |
| 2 | 1.9943 | 2.9895 |
| 3 | 2.9988 | 4.0028 |

Thus the hp basis, quadrature, Robin assembly, load, sparse solve, and exact
error integration pass the prerequisite calibration.

At `L_H=4`, `L_h=8`, the neighboring-localization differences
`||u^ell-u^(ell-1)||_(1,k)` fall from about `0.135-0.177` at `ell=2` to
`0.0070-0.0082` at `ell=3`, for `p=1,2,3`. At `ell=4` the patch already
covers the full domain and the difference is zero. Primal, adjoint,
constraint, and Petrov residuals remain between roughly `1e-16` and
`2e-14`. The frozen local rule used below is therefore `ell=3`.

The fixed `L_h=8`, `L_H=2,...,7` scan shows that P1 eventually approaches
its fine-space floor. P2/P3 delay that floor. Last-three-point log fits are:

| master fine | p | energy slope | energy R2 | L2 slope | L2 R2 |
|---:|---:|---:|---:|---:|---:|
| 8 | 2 | 1.9232 | 0.9980 | 3.2691 | 0.9978 |
| 8 | 3 | 1.7904 | 0.9997 | 3.0303 | 0.9997 |
| 10 | 2 | 1.5993 | 0.9975 | 2.7403 | 0.9981 |
| 10 | 3 | 1.5972 | 0.9975 | 2.7374 | 0.9982 |

The `L_h=10` fits use `L_H=4,5,6`. These data do not validate the
preregistered `O(H)` energy and `O(H^2)` L2 hypothesis. They show a stable
higher pre-asymptotic slope on the tested levels, while also showing that
raising only the fine corrector degree from P2 to P3 does not raise the
observed H order. No fit interval was selected after seeing the result.

The coupled `L_h-L_H=4` scan reproduces the same qualitative behavior. The
full preregistered `L_h=12`, `L_H=2,4,6,8` direct-saddle matrix is not a
suitable local WSL job because large low-H patches and P3 sparse saddle
factorizations dominate. It remains a server run:

```text
FULL=1 JOBS=8 ./scripts/run_helmholtz_hp_convergence_server.sh
```

Raw local CSV files and `/usr/bin/time -v` logs are in
`results/helmholtz_hp/`. Until the server matrix passes the fine-floor,
localization-floor, master-depth, and regression gates, the correct statement
is that HP4 is complete and HP5-HP6 passed local correctness/calibration, not
that the final asymptotic hypothesis has been numerically verified.
## 31. hp Corrector Parallelism and Long-Run Recovery

The initial HP6 server runner exposed only `JOBS`, which controls compilation;
the element-corrector loop itself was serial. A representative
`p=3,L_H=4,L_h=10,ell=3` direct-saddle case took 98.16 s at about one CPU
core, so a multi-level P1/P2/P3 matrix could run for many hours.

`HelmholtzHpLodModel` now accepts an explicit `corrector_threads` count.
Targets are ordered by estimated patch size and dispatched with dynamic
OpenMP scheduling. Each worker owns its sparse saddle solves, diagnostics,
and corrector triplet buffers; only completed buffers are merged. This keeps
Eigen sparse objects thread-local and avoids serialized insertion into the
global corrector matrix. The global fine-reference LU is also lazy and is
skipped when the benchmark supplies its reusable `FineReference`.

Final local measurements for the representative case are:

| patch threads | wall time | CPU | peak RSS |
|---:|---:|---:|---:|
| 1 | 98.16 s | 101% | 425 MiB |
| 8 | 34.70 s | 755% | 2.50 GiB |
| 16 | 34.97 s | lower scaling benefit | 4.68 GiB |

Eight threads give a 2.83x wall-time speedup locally. Sixteen threads do not
improve this machine and nearly double memory again, so the server default is
8 rather than the hardware thread count. The server should pilot 8/16/32 on
one case before selecting a larger value.

The convergence benchmark now supports `--threads`, `--progress`, and
`--stream`. Server runners flush one CSV row after every completed grid level,
write patch progress to `.time`, split heavy LOD scans into one recoverable
case per degree, and skip cases with an existing `.done` marker. Fine Pp
reference solutions are reused across every fixed-h H scan.

Correctness checks include serial-versus-four-thread correctors, adjoint
correctors, coarse operators, and final LOD solutions. Differences are below
`2e-11`. The final representative row retained the previous errors and gave
Petrov, corrector, and constraint residuals of `2.17e-14`, `1.75e-14`, and
`1.16e-15`. All 21 CTest cases pass, and both server scripts pass `bash -n`
and an end-to-end manufactured-solution smoke run.
## 32. hp DirectSchur Default

The existing `DirectSchur` patch solver was exposed through
`HelmholtzHpProblemConfig`, `HelmholtzHpPatchAssembler`, and
`bench_helmholtz_hp_convergence --solver=saddle|schur`. The server runners
also accept `HP_SOLVER`; their default is now `schur`, while the public model
configuration retains `DirectSaddle` as its compatibility default.

For one patch system, DirectSchur computes

```text
Z = A_patch^{-1} B^*,  Y = A_patch^{-1} F,
S = B Z,               S Lambda = B Y,
Q = Y - Z Lambda.
```

This is mathematically the same constrained corrector as the direct saddle
solve. In this hp experiment it is much cheaper because the constraints
cause substantial fill in the monolithic saddle LU, whereas the Helmholtz
A block remains sparse and the Schur complement is small.

Representative `p=3,L_H=4,L_h=10,ell=3` results:

| solver | threads | wall time | peak RSS |
|---|---:|---:|---:|
| DirectSaddle | 1 | 98.16 s | 425 MiB |
| DirectSaddle | 8 | 34.70 s | 2.50 GiB |
| DirectSchur | 1 | 7.41 s | 168 MiB |
| DirectSchur | 8 | 1.97 s | 628 MiB |
| DirectSchur | 16 | 2.09 s | 1.07 GiB |

Thus eight-thread DirectSchur is 17.6x faster than the parallel saddle path
and 49.8x faster than the original serial baseline. P1/P2/P3 tests at
`L_H=4,L_h=10,ell=3` retained the same manufactured, fine-reference, and LOD
errors. Maximum Schur residuals were below `4.6e-16`, minimum reciprocal
condition estimates were about `6.6e-2`, no fallback occurred, and final
corrector/Petrov residuals remained below `3e-14`.

The model now aggregates `max_schur_residual`, `min_schur_rcond`, and
`direct_fallback_count`. Streaming CSV output records these values and
`--check` rejects a large Schur residual or any fallback. This keeps the
faster route subject to the same correctness gate rather than treating a
successful factorization as sufficient evidence.

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

The `numerator=corrector` quotient is much larger than the original LOD-basis
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
computed below `--stability-max-dofs`; larger cases report `nan` instead of an
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

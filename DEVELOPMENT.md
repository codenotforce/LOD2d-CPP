# LOD2d-C++ Development Log

This file records the MATLAB-to-C++ migration decisions, correctness traps, and
performance experiments.

## Current Baseline

Environment used for the latest validation:

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
./build/benchmarks/bench_H4h8 --solver=eigen
./build/benchmarks/bench_H4h8 --solver=cholmod
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

### H=4, h=8, ell=2 (20-run median, OMP_NUM_THREADS=16)

| Version | Median | Range |
|---------|--------|-------|
| **C++ 16t** | **1.73 s** | 1.45-1.99 s |

Errors identical to MATLAB across all versions.

Eigen remains the default corrector solver.

## Correctness Decisions

### Mesh Refinement

The C++ refinement must match MATLAB sparse ordering exactly:

- MATLAB `find(sparse(row, col, value))` returns column-major order.
- Edge numbering therefore uses `(col, row)` ordering rather than the natural
  C++ row-major pair order.
- Refined triangles are grouped by sub-triangle type:
  `[all sub1; all sub2; all sub3; all sub4]`.
- `P_dg` is built as stacked sub-triangle prolongation blocks, matching MATLAB's
  `[kron(E, sub1); kron(E, sub2); ...]`.
- Explicit zero triplets are skipped because Eigen stores them as nonzeros.

### DG Assembly

`assemble_dg` now has a split implementation:

- `assemble_element_stiffness(mesh, coeff)` computes and returns all local 3x3
  element stiffness blocks.
- `assemble_dg_from_element_stiffness(blocks)` builds the global DG sparse
  matrix from those blocks.
- `assemble_dg(mesh, coeff)` remains the compatible one-call wrapper.

The coefficient vector length is checked against the element count to avoid
silent out-of-bounds reads in Release builds.

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
local patch.  For H4/h8 this means 512 correctors scanning 131072 fine elements.

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
| CHOLMOD as default corrector solver | Correct but slower than Eigen for H4/h8 |
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

## Remaining Opportunities

- Profile `C_ell + G` assembly and coarse/reference solves after corrector
  improvements; they are now a larger share of runtime.
- Consider a sparse column extraction helper for `G0`, `Sh_free`, and RHS
  assembly in benchmark/full pipeline code.
- Add CTest working-directory settings so `ctest --test-dir build` can find the
  golden files without running from the repository root.
- Keep Windows `D:\code\femcode\LOD2d_C++` as a mirror of the WSL-tested source
  to avoid line-ending and build-cache confusion.

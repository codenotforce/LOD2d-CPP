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

### H=4, h=9, ℓ=2 (7-run median, OMP_NUM_THREADS=16)

| Version | Median | Range |
|---------|--------|-------|
| **C++ 16t** | **9.79 s** | 9.48–10.01 s |
| MATLAB parallel (4w) | 23.57 s | 22.90–24.07 s |
| MATLAB serial | 49.39 s | 48.19–52.17 s |

Speedup: C++ 5.0× vs MATLAB serial, 2.4× vs MATLAB parallel.

### H=4, h=8, ℓ=2 (20-run median, OMP_NUM_THREADS=16)

| Version | Median | Range |
|---------|--------|-------|
| **C++ 16t** | **1.73 s** | 1.45–1.99 s |

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
- `Cholmod`: experimental, correctness checked but slower for current patch
  sizes.

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

CHOLMOD now passes golden corrector tests, but remains slower than Eigen for
the current LOD patch sizes.

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

## Failed or Rejected Experiments

| Experiment | Result |
|------------|--------|
| CHOLMOD as default corrector solver | Correct but slower than Eigen for H4/h8 |
| CHOLMOD for H4/h9 | Faster corrector phase than Eigen in one run, but slower total runtime |
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

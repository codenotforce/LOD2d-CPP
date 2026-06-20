# Benchmark and Corrector Coding Guide

This guide defines how LOD2d-C++ benchmarks should be written.  The goal is to
make every benchmark reproducible, comparable, memory-aware, and safe to clone
and run from a clean checkout.

## Scope

Use this guide for:

- `benchmarks/bench_H*` full or LOD-only benchmark drivers.
- `benchmarks/bench_profile.cpp` phase timing drivers.
- Any future benchmark that computes correctors, assembles `G`, solves the
  coarse LOD system, or optionally computes a reference solution.

Tests may use smaller or more direct code when they need to compare golden data,
but benchmark logic should follow this guide.

## Command-Line Options

Every full benchmark should support these options:

```bash
--solver=eigen|cholmod|cholmod_cached|auto
--threads=auto|env|N
--skip-reference
```

Recommended defaults:

- `--solver=auto`
- `--threads=auto`
- reference enabled unless `--skip-reference` is present

### Solver Policy

Use one enum from `CorrectorSolver` throughout the benchmark.

```cpp
if (value == "eigen") solver = CorrectorSolver::EigenLLT;
else if (value == "cholmod") solver = CorrectorSolver::Cholmod;
else if (value == "cholmod_cached") solver = CorrectorSolver::CholmodCached;
else if (value == "auto") solver = (h >= 10) ? CorrectorSolver::Cholmod
                                             : CorrectorSolver::EigenLLT;
```

Current empirical policy:

- `h <= 9`: Eigen is the default because it is stable and has lower overhead.
- `h >= 10`: plain CHOLMOD is faster for the corrector phase, but uses more memory.
- `cholmod_cached`: explicit experiment only; do not select it from `auto`
  unless a fresh benchmark proves it faster on the same machine.

Do not silently change this policy in one benchmark only.  If the policy changes,
update all benchmark drivers and this guide.

### Thread Policy

`--threads=auto` should cap h>=10 CHOLMOD runs to 8 OpenMP threads when the
environment requests more.  Current measurements on WSL show that 8 threads are
faster and much more memory-stable than 12 or 16 threads for h=10 CHOLMOD.

```cpp
#ifdef _OPENMP
void apply_thread_option(const Options &opt, int h) {
    if (opt.threads > 0) {
        omp_set_num_threads(opt.threads);
        return;
    }
    if (opt.threads == 0) return;  // --threads=env
    if (h >= 10 &&
        (opt.solver == CorrectorSolver::Cholmod || opt.solver == CorrectorSolver::CholmodCached) &&
        omp_get_max_threads() > 8)
        omp_set_num_threads(8);
}
#else
void apply_thread_option(const Options &, int) {}
#endif
```

Use `--threads=env` only when intentionally testing the raw `OMP_NUM_THREADS`
setting.

## Required Phase Order

Benchmarks should follow this phase order and print timings in this order:

1. Mesh refinement
2. Operator setup
3. Correctors
4. `C_ell + G`
5. Coarse LOD solve
6. Reference solve and error computation, unless skipped

Profile output should use this stable text shape:

```text
1. Mesh: ... ms -> Coarse:... Fine:...
2. Operators: ... ms (...)
3. Correctors: ... ms (... threads, ... elem/s)
4. C_ell + G: ... ms
5. Coarse solve: ... ms (SHLOD0: ...x...)
6. Reference: ... ms (DOFs: ...)
Total: ... ms
```

When reference is skipped:

```text
6. Reference: skipped
Total: ... ms
```

Full benchmark output should always show the selected solver and whether
reference was skipped.

## Corrector Calling Convention

Hot benchmark paths must use compact corrector output:

```cpp
std::vector<CorrectorEntries> CT(NTH);
#pragma omp parallel for schedule(dynamic)
for (int k = 0; k < NTH; ++k) {
    CT[k] = compute_corrector_entries(...);
}
```

Avoid returning or storing `SparseMatrix(Nh, 3)` for each corrector in benchmark
hot paths.  The sparse wrapper `compute_corrector(...)` remains for tests and
golden-data comparison.

### Required Corrector Caches

Build these once and pass them to every corrector:

```cpp
auto element_stiffness = assemble_element_stiffness(fine, Ah);
auto interpolation_rows = build_interpolation_rows(IH, NH);
auto fine_element_children = build_fine_element_children(f_out.P_elem, NTH);
```

Always pass these pointers when available:

```cpp
&element_stiffness
&fine_element_children
&interpolation_rows
```

Do not recompute `P0 * patch(:, k)` or scan all fine elements inside every
corrector when `fine_element_children` is available.

## Memory Lifetime Rules

Memory lifetime is part of benchmark correctness for h>=10.  The benchmark
should release large temporary data as soon as the next phases no longer need it.

### Operator Setup

Avoid long-lived DG matrices when direct CG matrices are enough.

Preferred:

```cpp
auto element_stiffness = assemble_element_stiffness(fine, Ah);
Eigen::SparseMatrix<double> Sh = assemble_cg_from_element_stiffness(fine, element_stiffness);
Eigen::SparseMatrix<double> Mh = assemble_cg_mass(fine, areas);
```

Avoid storing `Shdg` in benchmark drivers unless the code path actually needs it.
When `element_stiffness` is passed to `compute_corrector_entries`, `Shdg` is not
needed by the corrector hot path.

### Temporary DG-to-CG Map

`cg2dgh` is large for h=10.  Build it only inside a local scope for
`build_quasi_interp`:

```cpp
Eigen::SparseMatrix<double> IH;
{
    std::vector<Eigen::Triplet<double>> cg_t;
    cg_t.reserve(3 * NTh_f);
    ...
    Eigen::SparseMatrix<double> cg2dgh(3 * NTh_f, Nh);
    cg2dgh.setFromTriplets(cg_t.begin(), cg_t.end());
    IH = build_quasi_interp(coarse, fine, f_out.P_dg, cg2dgh, Nh, NH);
}
```

Do not keep `cg2dgh` alive after `IH` is built.

### Release After Cache Construction

After building `fine_element_children`, release `P_elem`:

```cpp
Eigen::SparseMatrix<double> empty_elem;
f_out.P_elem.swap(empty_elem);
```

After building `interpolation_rows`, release `IH`:

```cpp
Eigen::SparseMatrix<double> empty_ih;
IH.swap(empty_ih);
```

After correctors are computed, release corrector-only data:

```cpp
ElementStiffnessBlocks().swap(element_stiffness);
FineElementChildren().swap(fine_element_children);
InterpolationRows().swap(interpolation_rows);
Eigen::SparseMatrix<double> empty_pdg;
f_out.P_dg.swap(empty_pdg);
```

After assembling `G`, release compact correctors:

```cpp
std::vector<CorrectorEntries>().swap(CT);
```

Use this pattern only when the released data is truly no longer needed.  For
example, do not release `element_stiffness` before any later code that computes
energy norms from local element blocks.

## Building G

Do not build `cell_mat` and multiply by `cg2dgH` in benchmark hot paths.
Assemble `G = P_node - C_ell` directly from compact corrector entries:

```cpp
size_t corrector_nnz = 0;
for (const auto &entries : CT) corrector_nnz += entries.size();

std::vector<Eigen::Triplet<double>> g_t;
g_t.reserve(static_cast<size_t>(f_out.P_node.nonZeros()) + corrector_nnz);

for (int c = 0; c < f_out.P_node.outerSize(); ++c)
    for (Eigen::SparseMatrix<double>::InnerIterator it(f_out.P_node, c); it; ++it)
        g_t.emplace_back(it.row(), it.col(), it.value());

for (int k = 0; k < NTH; ++k)
    for (const auto &entry : CT[k])
        g_t.emplace_back(entry.row, coarse.elems[k][entry.col], -entry.value);

Eigen::SparseMatrix<double> G(Nh, NH);
G.setFromTriplets(g_t.begin(), g_t.end());
```

Always reserve with `P_node.nonZeros() + corrector_nnz`.

## Free-DOF Mapping

Use dense integer maps instead of `unordered_map` in performance-sensitive code.

Preferred:

```cpp
std::vector<int> dofH_map(NH, -1);
for (int i = 0; i < NH; ++i) {
    if (!is_dirH[i]) {
        dofH_map[i] = static_cast<int>(dofH.size());
        dofH.push_back(i);
    }
}
```

Then use:

```cpp
if (dofH_map[it.col()] >= 0) ...
```

Avoid repeated `unordered_map::count` inside sparse iterators.

## Reference Solve Policy

Reference solve is for validation, not for LOD performance timing.  It can be
very expensive for h>=10.

Benchmarks should support `--skip-reference` and clearly separate:

- LOD-only timing: mesh, operators, correctors, `G`, coarse solve.
- Full validation timing: LOD-only plus reference solve and error computation.

Do not compare LOD-only timing directly with a full MATLAB run that includes
reference computation unless the report clearly says so.

## Error Reporting

When reference is enabled, print all three errors:

```text
Errors: E=... L2=... FE=...
```

The same benchmark should be able to skip errors when `--skip-reference` is set.

## CHOLMOD Notes

CHOLMOD is correctness-checked but should remain an explicit or auto-selected
benchmark backend.  Do not make CHOLMOD the unconditional default.

Observed behavior:

- h=8 and h=9: Eigen is usually simpler and competitive.
- h=10: plain CHOLMOD is faster for correctors, but memory is higher.
- h=10 CHOLMOD: 8 OpenMP threads were faster and lower-memory than 12 or 16 on
  the current WSL machine.
- `cholmod_cached` reuses a bounded thread-local symbolic factor cache.  It is
  correct but currently experimental: an unbounded cache OOM-killed at about
  11.6 GB RSS, and the bounded one-pattern cache was slower than plain CHOLMOD
  in the current h=10 profile.

When changing CHOLMOD options, always compare:

```bash
/usr/bin/time -v ./build/benchmarks/bench_profile --solver=cholmod --threads=N --skip-reference
/usr/bin/time -v ./build/benchmarks/bench_profile --solver=cholmod_cached --threads=N --skip-reference
```

Record both wall time and `Maximum resident set size`.


## Reusing Correctors Across RHS Values

When `A`, mesh, `H/h`, `ell`, patches, and interpolation are unchanged, the LOD
correctors and multiscale basis are independent of the right-hand side.  Do not
recompute correctors for a batch of different forcing terms.

Preferred workflow:

```cpp
LodReusableSystem system(G, Sh, Mh, P_node, NH, coarse.dirichlet);
for (const Eigen::VectorXd &f : rhs_values) {
    LodReuseSolution sol = system.solve_from_coarse_values(f);
}
```

Use `solve_from_fine_values` when the RHS is already represented on fine nodes.
The benchmark entry is:

```bash
./build/benchmarks/bench_reuse_rhs --solver=auto --rhs=5
```

Report setup time separately from repeated RHS time.  The setup includes the
corrector computation and coarse factorization; repeated RHS timings should only
measure RHS projection, coarse back-substitution, and `G * uH`.
## Validation Checklist

Before committing benchmark changes, run at least:

```bash
cmake --build build --target test_corr test_full bench_H4h9 bench_profile -j 8
./build/tests/test_corr --solver=both
./build/tests/test_full
./build/benchmarks/bench_profile --solver=auto --skip-reference
```

For h>=10 or memory-sensitive changes, also run:

```bash
/usr/bin/time -v ./build/benchmarks/bench_profile --solver=auto --skip-reference
```

If changing solver policy, run both Eigen and CHOLMOD on the same benchmark and
record timing plus peak RSS.

## Documentation Checklist

When benchmark behavior changes, update:

- `README.md` for user-facing commands and current benchmark results.
- `DEVELOPMENT.md` for experiment history and rejected ideas.
- This guide for coding-policy changes.

Keep benchmark text ASCII-only unless there is a specific reason not to.  This
avoids encoding corruption in Windows/WSL round trips.

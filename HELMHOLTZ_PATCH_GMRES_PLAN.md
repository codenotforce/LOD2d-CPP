# Helmholtz Patch Solver And Two-Level Schwarz Plan

> Status date: 2026-07-21
>
> Corrector solver: M1-M4 complete; M5 performance policy open.
>
> Two-level Schwarz: S0-S4e complete as an experimental branch.
>
> Historical note: adaptive Helmholtz work was paused after stage 1 when this
> status was recorded. Current certified-adaptive work is tracked in
> [HELMHOLTZ_ADAPTIVE_LOD_PLAN.md](HELMHOLTZ_ADAPTIVE_LOD_PLAN.md).

## 1. Document Role

This plan defines the mathematics, implementation stages, acceptance gates,
and current execution status for two separate research tracks:

1. solving the constrained Helmholtz corrector patch problem;
2. using LOD patches in an experimental two-level Schwarz preconditioner.

User commands belong in [HELMHOLTZ_GUIDE.md](HELMHOLTZ_GUIDE.md). Measured
timings, failed parameter choices, and engineering decisions belong in
[DEVELOPMENT.md](DEVELOPMENT.md).

## 2. Corrector Patch Problem

### 2.1 Space And Boundary Conditions

For a coarse element `T`, let `omega=N^ell(T)` be its oversampling patch.
The localized fine-scale space is

```math
W_h(\omega)=\{w_h\in V_h:\operatorname{supp}w_h\subset\omega,
\ I_Hw_h=0\}.
```

Artificial patch boundaries use homogeneous Dirichlet data. Patch edges on
the physical boundary retain the impedance Robin term.

With fine patch basis functions, define

```math
A_\omega=K_\omega-k^2M_\omega-i kR_\omega.
```

The interpolation-kernel constraint is represented by `B_omega X=0`. For the
three local coarse basis right-hand sides collected in `F_T`, the discrete
corrector problem is

```math
\begin{pmatrix}
A_\omega & B_\omega^*\\
B_\omega & 0
\end{pmatrix}
\begin{pmatrix}X_T\\\Lambda_T\end{pmatrix}
=
\begin{pmatrix}F_T\\0\end{pmatrix}.
```

### 2.2 DirectSaddle Gold Standard

`DirectSaddle` factors the full constrained system. It remains the public
default because the constrained problem can be well posed even when the
unrestricted block `A_omega` is singular or close to a local resonance.

Acceptance gates:

- equation, constraint, primal, and adjoint residuals meet tolerance;
- internal, physical-boundary, and corner patches are covered;
- global Petrov-Galerkin solutions agree with direct references;
- no fallback or unconverged iterate is accepted silently.

### 2.3 DirectSchur

When `A_omega` is safely invertible, compute

```math
Y=A_\omega^{-1}F_T,\qquad
Z=A_\omega^{-1}B_\omega^*,
```

then

```math
S_\omega\Lambda_T=B_\omega Y,\qquad
S_\omega=B_\omega Z,\qquad
X_T=Y-Z\Lambda_T.
```

`DirectSchur` shares one sparse factorization across all constraint and
corrector right-hand sides. It must report factorization status, Schur
conditioning, Schur residual, equation residual, and constraint residual.
A detected singular or unreliable `A_omega` is an error unless an explicit
DirectSaddle fallback experiment was requested.

### 2.4 Shifted Right-Preconditioned GMRES

Use the shifted operator

```math
P_{\omega,\varepsilon}
=K_\omega-(k^2+i\varepsilon)M_\omega-i kR_\omega,
\qquad \varepsilon=\alpha k^2,
```

and solve

```math
A_\omega P_{\omega,\varepsilon}^{-1}y=g,
\qquad x=P_{\omega,\varepsilon}^{-1}y.
```

The solver must stop on the true residual

```math
\|g-A_\omega x\|_2/\|g\|_2,
```

not only a preconditioned residual. The shifted inverse has three controlled
forms:

- exact shifted SparseLU, used as the algebraic reference;
- identity, used only for the unpreconditioned GMRES baseline;
- fixed geometric V-cycle on the uniform NVB hierarchy.

A fixed linear V-cycle uses ordinary right GMRES. A changing or
tolerance-controlled inner inverse requires right FGMRES.

## 3. Corrector Solver Execution Status

| Stage | Deliverable | Status | Decision |
|---|---|---|---|
| M0 | Freeze DirectSaddle golden behavior | Complete | Retained as gold standard |
| M1 | Separate patch assembly, solver, scheduling, and output | Complete | No mathematical change |
| M2 | Implement and validate DirectSchur | Complete | Fastest tested large-patch experiment |
| M3 | Complex restarted right-GMRES with exact shifted inverse | Complete | Correct but slower than direct methods |
| M4 | Uniform-NVB geometric V-cycle | Complete | Mesh-robust iterations, but excessive multi-RHS cost |
| M5 | Establish crossover and automatic solver policy | Open | No automatic/default change justified |

### 3.1 Current Corrector Conclusions

- `DirectSaddle` is the most robust path and remains the default.
- `DirectSchur` is the primary performance candidate for large patches.
- Unpreconditioned GMRES converges too slowly to be competitive.
- Exact shifted-LU GMRES adds Krylov overhead around another factorization.
- The V-cycle reduces factorization memory and gives nearly mesh-independent
  iteration counts, but repeated solves for `m_omega+3` right-hand sides make
  it slower than DirectSchur in current two-dimensional tests.
- Correctors are independent of the external source `f`. With fixed operator,
  mesh, `H`, `h`, `ell`, interpolation, and boundary model, all correctors and
  coarse factors are reusable across right-hand sides.

### 3.2 Remaining M5 Work

M5 is complete only when at least one iterative or hybrid patch policy is
stably faster than the appropriate direct reference on multiple medium/large
cases while preserving every correctness gate.

Work order:

1. scan DirectSchur for local resonances over wider `k`, patch, and boundary
   classes;
2. report setup, factorization, multi-RHS solve, Schur, total time, and RSS
   separately;
3. evaluate block/recycling methods only where multiple right-hand sides share
   useful Krylov information;
4. evaluate symbolic-only or hierarchy reuse for heterogeneous coefficients;
5. define an automatic threshold only after a measured crossover exists;
6. prefer DirectSaddle near uncertain or failed thresholds.

## 4. Module Boundaries

| Module | Responsibility |
|---|---|
| `patch_system` | Patch topology, boundary classification, local matrices, constraints, and right-hand sides |
| `patch_solver` | DirectSaddle, DirectSchur, ShiftedGmres, recovery, and diagnostics |
| `shifted_laplacian` | Shift rule and shifted operator construction |
| `patch_hierarchy` | Restricted uniform-NVB hierarchy |
| `patch_multigrid` | Fixed V-cycle action |
| `right_gmres` | Complex right-GMRES and FGMRES algorithms |
| `corrector_pipeline` | Parallel scheduling and compact corrector output |

Solver modules must not rebuild mesh topology or reinterpret patch boundary
conditions.

## 5. Corrector Test Matrix

Required tests:

- complex nonnormal restarted GMRES, zero RHS, exact one-step preconditioner,
  restart, breakdown, iteration limit, and NaN/Inf rejection;
- exact identity `P_epsilon=A-i*epsilon*M`;
- DirectSchur and ShiftedGmres against DirectSaddle;
- physical Robin and artificial Dirichlet edge classification;
- V-cycle linearity, Galerkin operators, and level dimensions;
- primal/adjoint correctors, trial/test bases, Petrov matrix, final solution,
  and manufactured errors;
- failure paths with explicit patch ID and no silent fallback.

Numerical equivalence gates are

```math
\frac{\|X_{\mathrm{candidate}}-X_{\mathrm{saddle}}\|_F}
{\max(1,\|X_{\mathrm{saddle}}\|_F)}\le 10^{-8}
```

and matching true, equation, and constraint residual tolerances.

## 6. Cost Model

For `r=m_omega+3` right-hand sides, a direct Schur solve costs approximately

```math
C_{\mathrm{direct}}=C_{\mathrm{factor}}+rC_{\mathrm{tri}}.
```

An iterative shifted solve costs

```math
C_{\mathrm{iter}}=C_{\mathrm{setup}}+rqC_{\mathrm{inverse}},
```

plus Krylov orthogonalization. Iteration robustness alone is not a speedup;
the comparison must include every right-hand side, setup, and peak memory.

## 7. Reuse And Failure Rules

Numerical factors or hierarchies may be reused only when all operator-defining
data match: `k`, shift, coefficients, Robin parameter, physical boundary mask,
patch geometry, fine hierarchy, and numerical matrix values.

Heterogeneous coefficients may share topology or symbolic structure but not
stale numerical coarse operators. GMRES failure, non-finite data, local
factorization failure, singular Schur complement, or residual failure throws
by default. Explicit fallback must be counted and timed.

## 8. Two-Level Schwarz Extension

This is separate from the constrained corrector solver. Corrector matrices
cannot be inserted directly as Schwarz local inverses.

Let `B0` denote the factored Petrov-Galerkin LOD coarse action and let

```math
B_{\mathrm{loc}}r
=\sum_j R_j^TD_jA_j^{-1}R_jr,
\qquad \sum_jR_j^TD_jR_j=I.
```

The implemented actions are

```math
B_{\mathrm{add}}=B_0+B_{\mathrm{loc}},
```

and

```math
B_{\mathrm{hyb}}
=B_0+(I-B_0A_h)B_{\mathrm{loc}}(I-A_hB_0).
```

`WeightedOverlap` uses inverse multiplicity in `D_j`; `RestrictedCore` uses a
Boolean owner core. Artificial impedance plus weighted overlap is weighted
ORAS.

## 9. Schwarz Execution Status

| Stage | Deliverable | Status | Decision |
|---|---|---|---|
| S0-S1 | Coarse assembled weighted RAS pilot | Complete | Correct, not faster than coarse SparseLU |
| S2a | Fine-space local/additive/hybrid with LOD coarse action | Complete | Hybrid is strongest tested composition |
| S2b | Dedicated local assembler and artificial impedance | Complete | Dirichlet remains default |
| S3a | Shifted local GMRES with exact inverse | Complete | Correct, slower than direct local LU |
| S3b | Fixed local geometric V-cycle | Complete | Setup reduction does not offset solve cost |
| S4a | Fixed-`kH` robustness baseline | Complete | Hybrid iterations remain nearly constant in tested range |
| S4b | Identity/local/additive/hybrid comparison | Complete | Hybrid retained |
| S4c | Local high-wave-number study | Complete | Server equal-accuracy extension deferred |
| S4d | Weighted and RestrictedCore ORAS variants | Complete | Both remain research switches |
| S4e | Strict identical-matrix factor reuse and multi-RHS solve | Complete | Reduces medium/large setup, solve, and RSS |

### 9.1 Current Schwarz Conclusions

- Recommended research configuration: Dirichlet artificial boundaries,
  weighted overlap, direct local SparseLU, and optional strict identical-matrix
  reuse.
- Impedance ORAS and Boolean RestrictedCore did not improve the tested cases.
- Hybrid iteration counts show useful wave-number robustness, but the complete
  fine-space Schwarz method has not stably beaten fine SparseLU in total time.
- `factorization_reuse=identical` is explicit rather than default because
  grouping overhead can dominate small local systems.

### 9.2 Next Schwarz Decision

The S4 program is complete. Before adding more solver variants, choose one of
the following based on target problem size:

1. stop the two-dimensional performance branch and retain it as a scalability
   reference for larger or three-dimensional problems; or
2. continue with heterogeneous-coefficient symbolic/hierarchy reuse and an
   equal-accuracy server study.

No default solver changes until total runtime and correctness both improve.

## 10. Documentation And Reproduction

- Commands and stable options: [HELMHOLTZ_GUIDE.md](HELMHOLTZ_GUIDE.md).
- Shared benchmark rules: [BENCHMARK_GUIDE.md](BENCHMARK_GUIDE.md).
- Numerical evidence and rejected ideas: [DEVELOPMENT.md](DEVELOPMENT.md).
- Adaptive research stages: [HELMHOLTZ_ADAPTIVE_LOD_PLAN.md](HELMHOLTZ_ADAPTIVE_LOD_PLAN.md).

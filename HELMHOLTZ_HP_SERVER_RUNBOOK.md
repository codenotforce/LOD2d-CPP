# Helmholtz hp LOD Server Runbook

This runbook builds and validates the continuous `P1/P2/P3` fine-space
Helmholtz LOD experiment, then runs fixed-master-fine scans on a Linux
server. The current correctness reference is the direct saddle patch solver.

## 1. Hardware and storage

The `deep` cases are intended for the AMD EPYC server with 377 GiB RAM.
Use a local Linux filesystem, not an NFS-mounted build directory. Keep at
least 20 GiB free for build products and result logs.

`JOBS` controls compilation only. `PATCH_THREADS` controls independent
corrector patch solves. Each worker owns a sparse saddle factorization and
triplet buffers, so increasing `PATCH_THREADS` trades memory for wall time.
Start with `PATCH_THREADS=8`; only increase it after checking peak RSS.

## 2. Clone or update

For a new checkout:

```bash
git clone https://github.com/codenotforce/LOD2d-CPP.git
cd LOD2d-CPP
```

For an existing clean checkout:

```bash
cd ~/code/LOD2d-CPP
git fetch origin main
git switch main
git pull --ff-only origin main
```

Check the revision:

```bash
git status --short
git log --oneline -3
```

Do not use `git reset --hard` when the server checkout contains unreturned
results.

## 3. Install dependencies

Ubuntu:

```bash
sudo apt update
sudo apt install -y build-essential cmake g++ \
  libeigen3-dev libsuitesparse-dev libtbb-dev time
```

## 4. Release build

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DLOD2D_BUILD_TESTS=ON \
  -DLOD2D_BUILD_BENCHMARKS=ON

cmake --build build --target \
  lod2d_core \
  test_helmholtz_hp_fem \
  test_helmholtz_hp_patch \
  test_helmholtz_hp_model \
  bench_helmholtz_hp_convergence \
  -j 16
```

## 5. Correctness gate

Run the hp-specific tests first:

```bash
./build/tests/test_helmholtz_hp_fem
./build/tests/test_helmholtz_hp_patch
./build/tests/test_helmholtz_hp_model
```

Then run the complete suite:

```bash
ctest --test-dir build --output-on-failure
```

All tests must pass before starting a long scan. The hp model test checks
`p=1` equivalence with the original P1 implementation and repeated-RHS
corrector reuse. It also compares serial and four-thread corrector`nassembly, coarse operators, and final LOD solutions.

## 6. Script permissions and smoke test

```bash
chmod +x scripts/run_helmholtz_hp_convergence_server.sh
MODE=smoke JOBS=16 PATCH_THREADS=8 HP_SOLVER=schur \
  ./scripts/run_helmholtz_hp_convergence_server.sh
```

Expected files:

```text
results/helmholtz_hp/smoke_H4_h8.csv
results/helmholtz_hp/smoke_H4_h8.time
results/helmholtz_hp/smoke_H4_h8.done
```

The `.done` file is created only after a successful benchmark. Rows are
flushed into `.csv.tmp` as each grid level completes, then the file is
renamed to `.csv` after success. The `.time` file receives patch progress
such as `correctors=16/128`. Heavy LOD scans are split by polynomial degree,
so rerunning the same command skips every finished `p` case and resumes at
the first missing `.done` file.

Monitor a running case with:

```bash
tail -f results/helmholtz_hp/*.time
watch -n 5 'free -h; ps -C bench_helmholtz_hp_convergence \
  -o pid,etime,%cpu,rss,vsz,cmd --sort=-rss'
```

### Patch solver and thread tuning

The runners default to `HP_SOLVER=schur`. `DirectSchur` factors the sparse
patch Helmholtz block, applies it to all constraint and corrector right-hand
sides, then solves the small dense Schur complement. Use
`HP_SOLVER=saddle` only for a reference comparison or if a Schur case fails
the explicit residual/conditioning checks.

On the local WSL machine, the representative
`p=3,L_H=4,L_h=10,ell=3` case produced:

| solver | patch threads | wall time | peak RSS |
|---|---:|---:|---:|
| DirectSaddle | 1 | 98.16 s | 425 MiB |
| DirectSaddle | 8 | 34.70 s | 2.50 GiB |
| DirectSchur | 1 | 7.41 s | 168 MiB |
| DirectSchur | 8 | 1.97 s | 628 MiB |
| DirectSchur | 16 | 2.09 s | 1.07 GiB |

Eight-thread DirectSchur is 17.6x faster than eight-thread DirectSaddle and
49.8x faster than the original serial path for this case. Its exact, fine,
and LOD errors agree with DirectSaddle at printed precision. The CSV also
records `schur_residual`, `schur_rcond`, and `direct_fallbacks`; accepted runs
require small residuals and zero fallbacks.

The server has a different memory hierarchy, so pilot 8, 16, and 32 threads
on one representative case before the full matrix. Keep the fastest setting
that does not swap; more threads are not automatically faster.

## 7. Preregistered h=12 matrix

Run this before deeper fine spaces:

```bash
MODE=full JOBS=16 PATCH_THREADS=8 HP_SOLVER=schur \
  ./scripts/run_helmholtz_hp_convergence_server.sh
```

It performs:

- fine hp-FEM calibration through `L_h=12`;
- `ell=1,...,5` calibration at `L_H=4,L_h=10`;
- fixed-master scan `L_h=12`, `L_H=2,4,6,8`;
- coupled scan with `L_h-L_H=6`.

Monitor memory in another terminal:

```bash
watch -n 5 'free -h; ps -C bench_helmholtz_hp_convergence \
  -o pid,etime,%cpu,rss,vsz,cmd --sort=-rss'
```

## 8. Deep h=12,14 comparison

After the full matrix succeeds:

```bash
MODE=deep \
DEEP_FINE_LEVELS="12 14" \
DEEP_H_LEVELS=6,8 \
DEGREES=1,2,3 \
ELL_BY_P=3,3,3 \
JOBS=16 PATCH_THREADS=8 \
  ./scripts/run_helmholtz_hp_convergence_server.sh
```

The default deep scan intentionally uses L_H=6,8. Each degree is written
to its own CSV and .done marker. Combining very coarse
`L_H=2` with `L_h=14`, `p=3`, and a whole-domain patch creates much larger
direct sparse factorizations and is not the first deep-h experiment.

For a staged run, start with P1:

```bash
MODE=deep DEGREES=1 DEEP_FINE_LEVELS="12 14" \
  ./scripts/run_helmholtz_hp_convergence_server.sh
```

Then run P2 and P3 separately into distinct directories:

```bash
MODE=deep DEGREES=2 RESULT_DIR="$PWD/results/helmholtz_hp_p2" \
  ./scripts/run_helmholtz_hp_convergence_server.sh

MODE=deep DEGREES=3 RESULT_DIR="$PWD/results/helmholtz_hp_p3" \
  ./scripts/run_helmholtz_hp_convergence_server.sh
```

Stop escalation if the previous level swaps, approaches the machine memory
limit, or has not completed cleanly. Do not infer a completed case from a
partial `.csv.tmp`.

## 9. Manufactured-solution comparison

The convergence benchmark already uses

```text
u(x,y) = phi(x) phi(y) exp(i k x),
phi(t) = 16 t^2 (1-t)^2.
```

Both `u` and its normal derivative vanish on the square boundary, so it
exactly satisfies the homogeneous impedance condition. The dedicated runner
makes all three comparisons explicit:

```text
exact_*    = error between the LOD solution and the manufactured solution
fine_*     = error between the fine Pp FEM solution and the manufactured solution
lod_fine_* = error between the LOD solution and the fine Pp FEM solution
```

Start with:

```bash
chmod +x scripts/run_helmholtz_hp_manufactured_server.sh
MODE=smoke PATCH_THREADS=8 \
  ./scripts/run_helmholtz_hp_manufactured_server.sh
```

Then compare fixed master fine spaces `L_h=10,12,14`:

```bash
MODE=fixed \
MASTER_FINE_LEVELS="10 12 14" \
H_LEVELS=4,6,8 \
DEGREES=1,2,3 PATCH_THREADS=8 \
  ./scripts/run_helmholtz_hp_manufactured_server.sh
```

Run the coupled check separately:

```bash
MODE=coupled \
COUPLED_H_LEVELS=2,4,6 \
GAP=6 \
DEGREES=1,2,3 PATCH_THREADS=8 \
  ./scripts/run_helmholtz_hp_manufactured_server.sh
```

Use `MODE=all` to run both. Results are stored under
`results/helmholtz_hp_manufactured/`; `manufactured_solution.txt` records
the exact solution and the meaning of every error column. Existing `.done`
cases are skipped, so the runner can resume safely.

## 10. Result checks

For each CSV:

```bash
column -s, -t results/helmholtz_hp/H_deep_h12_p1.csv | less -S
tail -n 25 results/helmholtz_hp/H_deep_h12.time
```

Accept a row only when:

- `petrov_residual`, `corrector_residual`, `constraint_residual`, and
  `schur_residual` are below `1e-9`;
- `schur_rcond` is positive and not anomalously small compared with nearby
  cases, and `direct_fallbacks` is zero;
- errors are finite;
- the command exit status is zero and the matching `.done` exists;
- `/usr/bin/time -v` reports no swapping or termination.

The final H-order interpretation must also apply the fine-space and
localization-floor gates in
`HELMHOLTZ_HP_CORRECTOR_CONVERGENCE_PLAN.md`.

## 11. Return results

From the workstation:

```bash
scp -r user@server:~/code/LOD2d-CPP/results/helmholtz_hp \
  /home/qcxubuntu/learning/LOD2d-C++/results/
```

If separate P2/P3 directories were used, return those directories as well.
Keep the CSV, `.time`, and `.done` files together.
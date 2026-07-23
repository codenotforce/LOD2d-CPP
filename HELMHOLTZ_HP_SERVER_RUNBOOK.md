# Helmholtz hp LOD Server Runbook

This runbook builds and validates the continuous `P1/P2/P3` fine-space
Helmholtz LOD experiment, then runs fixed-master-fine scans on a Linux
server. The current correctness reference is the direct saddle patch solver.

## 1. Hardware and storage

The `deep` cases are intended for the AMD EPYC server with 377 GiB RAM.
Use a local Linux filesystem, not an NFS-mounted build directory. Keep at
least 20 GiB free for build products and result logs.

The implementation currently solves patch problems sequentially. `JOBS`
controls compilation only; increasing it does not increase corrector
parallelism.

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
corrector reuse.

## 6. Script permissions and smoke test

```bash
chmod +x scripts/run_helmholtz_hp_convergence_server.sh
MODE=smoke JOBS=16 \
  ./scripts/run_helmholtz_hp_convergence_server.sh
```

Expected files:

```text
results/helmholtz_hp/smoke_H4_h8.csv
results/helmholtz_hp/smoke_H4_h8.time
results/helmholtz_hp/smoke_H4_h8.done
```

The `.done` file is created only after a successful benchmark. The CSV is
first written to `.tmp` and renamed after success, so interrupted runs do
not look complete.

## 7. Preregistered h=12 matrix

Run this before deeper fine spaces:

```bash
MODE=full JOBS=16 \
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
JOBS=16 \
  ./scripts/run_helmholtz_hp_convergence_server.sh
```

The default deep scan intentionally uses `L_H=6,8`. Combining very coarse
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
MODE=smoke \
  ./scripts/run_helmholtz_hp_manufactured_server.sh
```

Then compare fixed master fine spaces `L_h=10,12,14`:

```bash
MODE=fixed \
MASTER_FINE_LEVELS="10 12 14" \
H_LEVELS=4,6,8 \
DEGREES=1,2,3 \
  ./scripts/run_helmholtz_hp_manufactured_server.sh
```

Run the coupled check separately:

```bash
MODE=coupled \
COUPLED_H_LEVELS=2,4,6 \
GAP=6 \
DEGREES=1,2,3 \
  ./scripts/run_helmholtz_hp_manufactured_server.sh
```

Use `MODE=all` to run both. Results are stored under
`results/helmholtz_hp_manufactured/`; `manufactured_solution.txt` records
the exact solution and the meaning of every error column. Existing `.done`
cases are skipped, so the runner can resume safely.

## 10. Result checks

For each CSV:

```bash
column -s, -t results/helmholtz_hp/H_deep_h12.csv | less -S
tail -n 25 results/helmholtz_hp/H_deep_h12.time
```

Accept a row only when:

- `petrov_residual`, `corrector_residual`, and `constraint_residual` are
  below `1e-9`;
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
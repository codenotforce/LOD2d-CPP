# Local absolute-energy convergence validation

This directory is a small, versioned schema and correctness sample for
`bench_helmholtz_H_convergence`. It is not the large-wave-number paper run.

The WSL command was:

```bash
OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PROC_BIND=spread OMP_PLACES=cores \
build-h-convergence-local/benchmarks/bench_helmholtz_H_convergence \
  --k=8 --h=11 --H-levels=4,5,6,7,8 --ell=3 \
  --solver=schur --threads=8 \
  --export-dir=results/helmholtz_H_convergence_local.tmp/snapshots \
  --export-fields --fine-reference \
  --summary-out=results/helmholtz_H_convergence_local.tmp/summary.csv \
  --check --format=csv
```

`summary.csv` contains five convergence rows. The `snapshots/` directory
contains every coarse mesh, the fixed fine mesh once, and fine-node fields at
every coarse level so that the visualization project can develop and test its
reader without waiting for the EPYC run.

The registered server parameters and interpretation rules are in
`HELMHOLTZ_H_CONVERGENCE_SERVER_RUNBOOK.md`.

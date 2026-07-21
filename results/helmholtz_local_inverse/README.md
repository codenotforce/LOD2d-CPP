# Helmholtz Local Inverse Experiment Results

This directory contains the 2026-07-21 WSL runs described in section 28 of
`DEVELOPMENT.md`.

Primary fixed-fine-space files:

- `summary_{fraction,single-chain,boundary-chain}.csv`: one row per iteration,
  basis, and denominator;
- `elements_*.csv`: per-coarse-element inverse constants and diagnostics;
- `mesh_*.csv`: stable element IDs, parent IDs, levels, and coordinates;
- `time_*.txt`, `run_*.log`, `metadata.txt`: commands and resources.

Auxiliary scans:

- `hscan_L11.csv` through `hscan_L14.csv`: fixed final coarse grid, varying
  master fine level;
- `time_hscan_L14.txt`: successful level-14 resource record;
- `time_hscan_L15.txt`: signal-9/OOM record at 11,645,504 KiB peak RSS;
- `threshold_1e-10.csv`, `threshold_1e-14.csv`: positive local-mass rank
  sensitivity (`1e-12` is the primary summary);
- `ell_2.csv`, `ell_4.csv`: oversampling sensitivity (`ell=3` is primary);
- `k_1.csv`, `k_4.csv`: wave-number sensitivity (`k=2` is primary);
- `fraction_0.10.csv`, `fraction_0.50.csv`: marking-fraction sensitivity
  (`0.25` is primary).

Feedback-driven refinement:

- `summary_feedback_argmax-{element,patch}_n{0,1}.csv`: each next mesh marks
  the current trial element- or patch-denominator maximizer; `n=0` marks only
  that element and `n=1` also marks its one-layer vertex neighbors;
- matching `elements_`, `mesh_`, `time_`, and `run_` files retain the local
  spectra, hierarchy trajectory, and resource diagnostics.

Oversampling-matched denominator and deep-`h` scans:

- `hmatched_patch3_L{11,12,13}.csv`: the same level-`3:5` coarse mesh with
  `ell=3` and denominator `omega_T^3`, giving `q_max=0.125,0.0884,0.0625`;
- `summary_feedback_argmax-patch_n0_matched_ell3_Lh13.csv`: successful
  feedback driven by the current `patch3` maximizer;
- files ending in `matched_ell3_Lh14` record the signal-9 endpoint at
  `11,668,680 KiB` peak RSS and are intentionally incomplete.

Every successful `--check` run verifies exact coarse-to-fine nesting through
nodal, element, DG, and quasi-interpolation prolongation identities. The empty
`hscan_L15.csv` is intentionally retained beside its time log to document the
failed high-memory endpoint; it is not a completed data point.

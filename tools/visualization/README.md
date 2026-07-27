# Helmholtz paper plots

This directory implements the first CSV-only visualization stage from
`VISUALIZATION_PLAN.md`. It is a post-processing layer: it does not invoke the
C++ solver, add model state, retain correctors, or alter benchmark timings.

## Figures

`plot_helmholtz_results.py` creates:

- `helmholtz_pollution_weighted_energy.{png,pdf}` from the manufactured
  fixed-resolution wave-number scans in:
  - `results/helmholtz_pollution_gap6/`;
  - `results/helmholtz_pollution/`;
  - `results/helmholtz_pollution_server/`;
- `manufactured_global_nvb_convergence.{png,pdf}` from
  `results/helmholtz_manufactured/validation.csv`;
- `helmholtz_plot_metrics.csv`, containing fitted rates and the numerical
  changes annotated in the figures.

The pollution figure uses fixed `kH=1`. Its primary panel uses `kh=1/8`; its
strict-reference panel uses `kh=1/16`. The coarse P1 FEM exact error is
independent of the auxiliary fine space. The fine-P1 exact error is included
to distinguish coarse-method pollution from a deteriorating reference.

The supported conclusion is limited to the checked-in manufactured scans:
standard coarse P1 FEM shows strong weighted-energy pollution through `k=64`,
while the tested two-sided LOD does not. This is numerical evidence, not a
theorem or a claim for arbitrary wave numbers.

## Reproduce in WSL

The archived environment uses Python 3.10.12, Matplotlib 3.5.1, and NumPy
1.21.5. No package installation is needed when those system packages are
already present.

```bash
cd /home/qcxubuntu/learning/LOD2d-C++
python3 tools/visualization/plot_helmholtz_results.py
python3 -m unittest tools/visualization/test_plot_helmholtz_results.py
```

To generate SVG in addition to PNG/PDF:

```bash
python3 tools/visualization/plot_helmholtz_results.py \
  --formats=png,pdf,svg
```

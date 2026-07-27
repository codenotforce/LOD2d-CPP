# Helmholtz paper plots

This directory implements stages 1--3 from `VISUALIZATION_PLAN.md`. The normal
solver and existing benchmarks remain export-free. The separate
`bench_helmholtz_visualization` executable explicitly performs the additional
reference solve and streams non-owning mesh/field views to VTU.

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

The convergence plot uses increasing degrees of freedom on the horizontal
axis. Consequently the fitted exponents are negative and every convergent
series runs from upper left to lower right. The metrics CSV names these values
`*_dof_fitted_slope`.

`plot_helmholtz_snapshot.py` reads the explicit VTU/JSON case and creates:

- coarse/fine mesh comparison;
- exact, fine FEM, LOD, coarse prolongation, fine-scale correction, and error;
- real part, imaginary part, magnitude, and phase of the complex LOD field;
- a fixed centerline comparison that preserves the oscillatory detail.

`plot_helmholtz_H_convergence.py` validates and plots the completed EPYC-server
experiment in `results/helmholtz_H_convergence_server/all_results.csv`. The
current archived run uses `k=32`, fine level 19, coarse levels 8--13, and
`ell=4`. It produces the absolute `k`-weighted energy error versus coarse DOFs,
the successive measured `H`-orders, and a CSV containing fitted DOF slopes and
residual maxima.

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
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target bench_helmholtz_visualization -j
./build/benchmarks/bench_helmholtz_visualization \
  --k=4 --H=4 --h=8 --ell=3 \
  --output-dir=results/visualization/helmholtz_manufactured_k4_H4_h8_ell3
python3 tools/visualization/plot_helmholtz_results.py
python3 tools/visualization/plot_helmholtz_snapshot.py
python3 tools/visualization/plot_helmholtz_H_convergence.py
python3 -m unittest \
  tools/visualization/test_plot_helmholtz_results.py \
  tools/visualization/test_helmholtz_snapshot.py \
  tools/visualization/test_helmholtz_H_convergence.py
```

To generate SVG in addition to PNG/PDF:

```bash
python3 tools/visualization/plot_helmholtz_results.py \
  --formats=png,pdf,svg
```

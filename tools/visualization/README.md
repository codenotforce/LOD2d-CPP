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

`plot_helmholtz_H_convergence.py` validates and plots the latest completed
EPYC-server experiment in
`results/helmholtz_H_convergence_server/all_results.csv`. It infers the common
wave number, fine level, coarse-level interval, oversampling level, and solver
from the rows. It produces a four-panel figure containing the absolute
`k`-weighted energy and `L^2` errors versus coarse DOFs and both sets of
successive measured `H`-orders. The metrics CSV contains fitted DOF slopes,
P1/LOD ratios, fine-reference floors, and residual maxima. When the fine exact
errors are available, both error panels display the corresponding floor.

`plot_helmholtz_method_comparison.py` compares completed case-S PALOD, SLOD,
uniform P1 FEM, and adaptive P1 FEM trajectories using the manufactured exact
solution. It requires every input run to be `success/TrajectoryComplete`,
checks the common case, wave number, initial coarse level, and initial coarse
DoF, and truncates each comparator at its first evaluated point at or below the
terminal PALOD exact relative weighted-energy error. The figure marks the
three PALOD reference-epoch starts and includes $N^{-1/2}$ and $N^{-1/3}$
reference slopes. Each method legend reports a whole-displayed-range log-log
DoF exponent. Repeated PALOD DoFs retain the later epoch value for this fit, so
vertical reference-refresh improvements are not counted as DoF convergence.
The script writes the plotted observations, first-crossing summary, fitted
rates, and PALOD-versus-$N^{-1/2}$ comparison beside the PNG and PDF outputs;
none of these post-processing quantities feed MARK or STOP.

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
python3 tools/visualization/plot_helmholtz_method_comparison.py \
  --palod=results/.../S_PALOD_run \
  --slod=results/.../S_SLOD_run \
  --ufem=results/.../S_UFEM_run \
  --afem=results/.../S_AFEM_run \
  --output=results/.../S-k16-error-vs-DoF.png
python3 -m unittest \
  tools/visualization/test_plot_helmholtz_results.py \
  tools/visualization/test_helmholtz_snapshot.py \
  tools/visualization/test_helmholtz_H_convergence.py \
  tools/visualization/test_plot_helmholtz_method_comparison.py
```

To generate SVG in addition to PNG/PDF:

```bash
python3 tools/visualization/plot_helmholtz_results.py \
  --formats=png,pdf,svg
```

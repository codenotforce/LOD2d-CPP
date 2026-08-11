# E0 R1 / kappa=16 calibration

This directory freezes the practical parameters used after the E0 gate.  It
is an implementation-study calibration, not a performance result and not a
fully verified certificate.

The accepted fixed hierarchy is coarse level 6, reference level 8, with the
ambient shadow refined until `rho_amb <= 0.25`.  The initially attempted
level-4/level-6 hierarchy was rejected because the ambient kernel was not
Helmholtz coercive at `kappa=16`.

Reproduce from a WSL Release build with:

```bash
experiments/helmholtz_adaptive_paper/run_e0_calibration.sh \
  "$PWD/build-release" "$PWD/build-release/e0-R1-k16"
```

The runner checks:

- `Theta_loc` for `ell=1,2,3,4` and its decrease with oversampling;
- direct ideal/localized reference-corrector perturbation against the
  ambient-to-reference one-sided certificate;
- local kernel constraints, Riesz stationarity, energy identities, and
  conservative element allocation;
- the local effectivity distribution for smooth R1, including a bounded P90
  to median ratio;
- a real `PracticalAdaptiveDriver` trajectory that performs
  `ell=1 -> IncreaseGlobalEll -> ell=2 -> Complete` without changing `H` or
  requesting a reference refresh.

`C0_usr` covers the observed well-localized multiplier
`||u_ref-U||_kappa / eta_H`.  `C1_usr` covers only the observed excess that
correlates with `Theta_loc`; both use a 1.25 safety factor.  These are empirical
parameters for the frozen paper setup, not theorem constants.  Reference
errors are computed only as calibration/post-processing diagnostics and never
enter the marking or stopping decisions of production trajectories.

Files:

- `e0-calibration.csv`: raw four-row calibration evidence;
- `frozen-parameters.json`: accepted parameters and driver-smoke evidence.

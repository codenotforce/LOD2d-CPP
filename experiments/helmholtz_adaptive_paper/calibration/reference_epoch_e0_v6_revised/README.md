# Revised schema-v6 reference-epoch E0

This directory freezes the implementation-study calibration used by the
revised E1/E2 configurations.  It was generated on 2026-08-24 with
`bench_helmholtz_reference_epoch_e0 --check` against manuscript SHA-256
`71f59581ea3a5e4cd65659055915715b9c86c582793db7154a5f0b3b31843ca8`.

The six numerical artifacts cover the R1 quadrature cross-check, localization/direct
defect bracket, reference Riesz identities, RT2/P2 reconstruction, candidate
dual lower bound, and the practical frozen values
`theta_loc_usr=2.6311141730106273` and
`C_rel_usr=1.4706194473212615`.

`06-calibration-provenance.json` binds them to generator commit
`042cdef892ba89027aa1043acedb21d4d203c6e9`, binary SHA-256
`9cae113b7172bf52e145ddccd5c01180897062432c28fadcae9a27a255278ea9`, the
manuscript SHA-256 above, GCC 11.4.0, and eight OpenMP threads.  The second Git
commit adds only these generated calibration artifacts, so the recorded first
commit is the exact source used to build the generator.

These constants are empirical and the production runs remain labelled
`implementation-study`; the files are not an interval-arithmetic proof of the
continuous or discrete stability constants.

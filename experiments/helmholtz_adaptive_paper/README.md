# Helmholtz adaptive paper experiment protocol

This directory contains the versioned, frozen inputs for the certified adaptive
Helmholtz LOD paper experiments. `schema-v1.json` defines the per-run
configuration. The C++ registry in `helmholtz/experiments/paper_config.h` is the
executable counterpart and rejects unknown fields and schema versions.

Run IDs have the form `case_method_kN_rN_hash`. The final component is the
lowercase 16-digit FNV-1a 64-bit hash of the canonical JSON configuration. The
canonical configuration includes the Git commit, build hash, and repeat index,
so changing any provenance field creates a different immutable run ID.

`HLOD-proxy` is intentionally present only as a diagnostic method and must not
be included in the six-method paper comparator matrix.

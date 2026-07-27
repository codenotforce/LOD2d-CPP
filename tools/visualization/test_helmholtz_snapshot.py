"""Schema, algebra, and headless-render tests for Helmholtz VTU snapshots."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY_ROOT = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

import plot_helmholtz_snapshot as snapshot


class HelmholtzSnapshotTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.case_dir = (
            REPOSITORY_ROOT
            / "results"
            / "visualization"
            / "helmholtz_manufactured_k4_H4_h8_ell3"
        )

    def test_schema_and_complex_algebra(self) -> None:
        manifest, coarse, fine = snapshot.load_case(self.case_dir)
        self.assertEqual(manifest["problem"], "manufactured_polynomial_plane_wave")
        self.assertGreater(fine[0].shape[0], coarse[0].shape[0])
        fields = fine[2]
        lod = snapshot.complex_field(fields, "u_lod")
        coarse_values = snapshot.complex_field(fields, "u_coarse")
        fine_scale = snapshot.complex_field(fields, "u_fine_scale")
        reference = snapshot.complex_field(fields, "u_reference")
        error = snapshot.complex_field(fields, "error_lod")
        np.testing.assert_allclose(lod, coarse_values + fine_scale, rtol=2e-14, atol=2e-14)
        np.testing.assert_allclose(error, lod - reference, rtol=2e-14, atol=2e-14)
        self.assertGreater(np.max(np.abs(fine_scale)), 1.0e-8)

    def test_headless_render(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            outputs = snapshot.generate_snapshot_figures(
                self.case_dir, Path(directory), ("png",)
            )
            self.assertEqual(len(outputs), 4)
            for path in outputs:
                self.assertTrue(path.is_file())
                self.assertGreater(path.stat().st_size, 5000)


if __name__ == "__main__":
    unittest.main()

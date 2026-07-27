"""Validation and headless-render tests for the server H-convergence plot."""

from __future__ import annotations

import math
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY_ROOT = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

import plot_helmholtz_H_convergence as convergence


class HelmholtzHConvergenceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.input_csv = (
            REPOSITORY_ROOT
            / "results"
            / "helmholtz_H_convergence_server"
            / "all_results.csv"
        )

    def test_server_rows_and_negative_dof_slopes(self) -> None:
        rows = convergence.load_and_validate(self.input_csv)
        self.assertEqual(len(rows), 5)
        self.assertTrue(all(convergence.integer(row, "k") == 16 for row in rows))
        self.assertTrue(all(convergence.integer(row, "h_level") == 18 for row in rows))
        self.assertEqual(
            [convergence.integer(row, "H_level") for row in rows],
            [9, 10, 11, 12, 13],
        )
        dofs = [convergence.integer(row, "coarse_nodes") for row in rows]
        p1 = [convergence.number(row, "p1_energy_abs") for row in rows]
        lod = [convergence.number(row, "lod_energy_abs") for row in rows]
        self.assertLess(convergence.power_slope(dofs, p1), 0.0)
        self.assertLess(convergence.power_slope(dofs, lod), 0.0)
        self.assertTrue(all(lod_value < p1_value for p1_value, lod_value in zip(p1, lod)))
        self.assertGreater(p1[-1] / lod[-1], 6.0)
        fine_floor = convergence.optional_number(rows[-1], "fine_energy_abs")
        self.assertIsNotNone(fine_floor)
        self.assertLess(lod[-1] / fine_floor, 1.05)

    def test_headless_render(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            outputs, metrics = convergence.generate(
                self.input_csv, Path(directory), ("png",)
            )
            self.assertEqual(len(outputs), 2)
            self.assertGreaterEqual(len(metrics), 6)
            for path in outputs:
                self.assertTrue(path.is_file())
                minimum_size = 100 if path.suffix == ".csv" else 1000
                self.assertGreater(path.stat().st_size, minimum_size)


if __name__ == "__main__":
    unittest.main()

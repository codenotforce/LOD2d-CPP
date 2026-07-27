"""Numerical/data-provenance tests for the Helmholtz paper plots."""

from __future__ import annotations

import math
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY_ROOT = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

import plot_helmholtz_results as plots


class HelmholtzPlotDataTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.results_root = REPOSITORY_ROOT / "results"

    def test_pollution_scan_supports_claim(self) -> None:
        datasets = plots.load_pollution_data(self.results_root)
        for name in ("main", "strict"):
            rows = datasets[name]
            self.assertEqual([plots.integer(row, "k") for row in rows], [4, 8, 16, 32, 64])
            fem = [plots.number(row, "fem_exact_energy_rel") for row in rows]
            lod = [plots.number(row, "lod_exact_energy_rel") for row in rows]
            self.assertGreater(plots.percent_change(fem[-2], fem[-1]), 80.0)
            self.assertGreater(fem[-1] / lod[-1], 25.0)

        strict_lod = [
            plots.number(row, "lod_exact_energy_rel")
            for row in datasets["strict"]
        ]
        self.assertLess(
            plots.percent_change(strict_lod[-2], strict_lod[-1]), 0.0
        )

    def test_manufactured_rates(self) -> None:
        fem, lod = plots.load_manufactured_data(self.results_root)
        fine_dofs = [plots.integer(row, "fine_nodes") for row in fem]
        fine_energy = [plots.number(row, "fem_exact_energy") for row in fem]
        fine_l2 = [plots.number(row, "fem_exact_l2") for row in fem]
        self.assertTrue(
            math.isclose(
                plots.power_rate(fine_dofs, fine_energy), -0.5, abs_tol=0.08
            )
        )
        self.assertTrue(
            math.isclose(plots.power_rate(fine_dofs, fine_l2), -1.0, abs_tol=0.08)
        )
        for ell in (2, 3):
            coarse_dofs = [plots.integer(row, "coarse_nodes") for row in lod[ell]]
            energy = [plots.number(row, "lod_exact_energy") for row in lod[ell]]
            l2_error = [plots.number(row, "lod_exact_l2") for row in lod[ell]]
            self.assertLess(plots.power_rate(coarse_dofs, energy), -0.65)
            self.assertLess(plots.power_rate(coarse_dofs, l2_error), -1.2)

    def test_headless_figure_smoke(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            outputs, metrics = plots.generate_figures(
                self.results_root, Path(directory), ("png",)
            )
            self.assertEqual(len(outputs), 3)
            self.assertGreaterEqual(len(metrics), 12)
            for path in outputs:
                self.assertTrue(path.is_file())
                self.assertGreater(path.stat().st_size, 1000)


if __name__ == "__main__":
    unittest.main()

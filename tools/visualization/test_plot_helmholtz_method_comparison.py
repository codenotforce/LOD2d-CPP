"""Tests for the manufactured-solution four-method comparison plot."""

from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import plot_helmholtz_method_comparison as plots


class MethodComparisonPlotTest(unittest.TestCase):
    def _run(self, root: Path, method: str, errors: list[float]) -> Path:
        directory = root / method
        directory.mkdir()
        metadata = {
            "case": "S",
            "method": method,
            "status": "success",
            "driver_state": "TrajectoryComplete",
            "config": {
                "wavenumber": 16,
                "boundary_beta": 1,
                "initial_coarse_level": 5,
                "manuscript_sha256": "test-manuscript",
                "quadrature": {"base_triangle_order": 7},
            },
        }
        (directory / "run.json").write_text(json.dumps(metadata), encoding="utf-8")
        fields = [
            "iteration", "reference_epoch", "N_H", "DoF_H",
            "relative_exact_energy_error",
        ]
        with (directory / "iterations.csv").open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            for index, error in enumerate(errors):
                writer.writerow(
                    {
                        "iteration": index,
                        "reference_epoch": min(index, 2) if method == "PALOD" else 0,
                        "N_H": 113 + 100 * index,
                        "DoF_H": 104 + 100 * index,
                        "relative_exact_energy_error": error,
                    }
                )
        return directory

    def test_truncates_at_first_palod_target_crossing(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            runs = [
                plots.load_method_run(self._run(root, "PALOD", [0.8, 0.4, 0.2])),
                plots.load_method_run(self._run(root, "SLOD", [0.9, 0.3, 0.15, 0.1])),
                plots.load_method_run(self._run(root, "UFEM", [1.1, 0.7, 0.4, 0.3])),
                plots.load_method_run(self._run(root, "AFEM", [1.0, 0.5, 0.18, 0.12])),
            ]
            output = root / "comparison.png"
            summary = plots.plot_method_comparison(runs, output)
            self.assertTrue(output.is_file())
            self.assertEqual(summary["first_target_point"]["SLOD"]["DoF_H"], 304)
            self.assertEqual(
                summary["first_target_point"]["SLOD"]["relative_exact_energy_error"],
                0.15,
            )
            self.assertIsNone(summary["first_target_point"]["UFEM"])
            self.assertEqual(summary["first_target_point"]["AFEM"]["DoF_H"], 304)
            self.assertTrue(output.with_suffix(".csv").is_file())
            self.assertTrue(output.with_suffix(".json").is_file())
            self.assertTrue(output.with_suffix(".pdf").is_file())

    def test_rejects_noncompleted_run(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = self._run(Path(temporary), "AFEM", [1.0, 0.5])
            path = directory / "run.json"
            metadata = json.loads(path.read_text(encoding="utf-8"))
            metadata["status"] = "interrupted"
            path.write_text(json.dumps(metadata), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "not success/TrajectoryComplete"):
                plots.load_method_run(directory)


if __name__ == "__main__":
    unittest.main()

"""Tests for exact-error fixed-horizon prefix validation."""

from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import compare_helmholtz_trajectory_prefix as prefixes


class TrajectoryPrefixTest(unittest.TestCase):
    def _run(self, root: Path, name: str, errors: list[float]) -> Path:
        directory = root / name
        directory.mkdir()
        (directory / "run.json").write_text(
            json.dumps(
                {
                    "case": "S",
                    "method": "SLOD",
                    "status": "success",
                    "driver_state": "TrajectoryComplete",
                    "config": {
                        "wavenumber": 16,
                        "boundary_beta": 1,
                        "initial_coarse_level": 5,
                        "quadrature": {"base_triangle_order": 7},
                    },
                }
            ),
            encoding="utf-8",
        )
        fields = [
            "reference_epoch",
            "N_H",
            "DoF_H",
            "ell",
            *prefixes.EXACT_FIELDS,
        ]
        with (directory / "iterations.csv").open(
            "w", encoding="utf-8", newline=""
        ) as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            for index, error in enumerate(errors):
                writer.writerow(
                    {
                        "reference_epoch": index,
                        "N_H": 100 + index,
                        "DoF_H": 90 + index,
                        "ell": 2,
                        "exact_energy_error": error,
                        "exact_L2_error": 0.5 * error,
                        "relative_exact_energy_error": 0.25 * error,
                        "relative_exact_L2_error": 0.125 * error,
                    }
                )
        return directory

    def test_accepts_reproduced_extended_prefix(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            baseline = self._run(root, "baseline", [1.0, 0.5])
            extended = self._run(root, "extended", [1.0, 0.5, 0.25])
            result = prefixes.compare_prefix(baseline, extended)
            self.assertEqual(result["status"], "prefix_reproduced")
            self.assertEqual(result["baseline_points"], 2)
            self.assertEqual(result["extended_points"], 3)

    def test_rejects_changed_prefix(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            baseline = self._run(root, "baseline", [1.0, 0.5])
            extended = self._run(root, "extended", [1.0, 0.4, 0.25])
            with self.assertRaisesRegex(ValueError, "differs"):
                prefixes.compare_prefix(baseline, extended)


if __name__ == "__main__":
    unittest.main()

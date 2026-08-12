"""Tests for adaptive-paper exact/reference error epoch plots."""

from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import plot_helmholtz_adaptive_epochs as plots


class AdaptiveEpochPlotTest(unittest.TestCase):
    def _run(self, root: Path, epoch: int) -> Path:
        directory = root / f"epoch-{epoch}"
        directory.mkdir()
        metadata = {
            "case": "S",
            "method": "PALOD",
            "status": "success",
            "config": {"reference_epoch": epoch, "wavenumber": 16},
        }
        (directory / "run.json").write_text(json.dumps(metadata), encoding="utf-8")
        fields = [
            "reference_epoch",
            "N_H",
            "DoF_H",
            "reference_energy_error",
            "reference_L2_error",
            "exact_energy_error",
            "exact_L2_error",
            "relative_exact_energy_error",
            "relative_exact_L2_error",
        ]
        with (directory / "iterations.csv").open(
            "w", encoding="utf-8", newline=""
        ) as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            for nodes, dofs, scale in ((225, 217, 1.0), (253, 245, 0.7)):
                writer.writerow(
                    {
                        "reference_epoch": epoch,
                        "N_H": nodes,
                        "DoF_H": dofs,
                        "reference_energy_error": 0.5 * scale,
                        "reference_L2_error": 0.3 * scale,
                        "exact_energy_error": 1.4 * scale,
                        "exact_L2_error": 0.08 * scale,
                        "relative_exact_energy_error": 0.55 * scale,
                        "relative_exact_L2_error": 0.32 * scale,
                    }
                )
        return directory

    def test_load_and_render_two_epochs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            runs = plots.load_epoch_runs((self._run(root, 0), self._run(root, 1)))
            self.assertEqual([run.epoch for run in runs], [0, 1])
            output = root / "error-vs-dof.png"
            plots.plot_epoch_errors(runs, output)
            self.assertTrue(output.is_file())
            self.assertGreater(output.stat().st_size, 1000)
            self.assertTrue(output.with_suffix(".csv").is_file())

    def test_exact_error_is_required_for_manufactured_case(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = self._run(Path(temporary), 0)
            path = directory / "iterations.csv"
            text = path.read_text(encoding="utf-8")
            path.write_text(text.replace("0.55", "", 1), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "missing exact errors"):
                plots.load_epoch_run(directory)


if __name__ == "__main__":
    unittest.main()

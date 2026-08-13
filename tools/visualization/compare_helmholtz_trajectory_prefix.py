#!/usr/bin/env python3
"""Verify that an extended fixed-horizon run reproduces an earlier prefix."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, List


Row = Dict[str, str]
EXACT_FIELDS = (
    "exact_energy_error",
    "exact_L2_error",
    "relative_exact_energy_error",
    "relative_exact_L2_error",
)
DISCRETE_FIELDS = ("reference_epoch", "N_H", "DoF_H", "ell")


def load_run(directory: Path) -> tuple[dict, List[Row]]:
    metadata = json.loads((directory / "run.json").read_text(encoding="utf-8"))
    if metadata.get("status") != "success" or metadata.get(
        "driver_state"
    ) != "TrajectoryComplete":
        raise ValueError(f"run is not success/TrajectoryComplete: {directory}")
    with (directory / "iterations.csv").open(
        "r", encoding="utf-8", newline=""
    ) as stream:
        rows = [
            row
            for row in csv.DictReader(stream)
            if row.get("relative_exact_energy_error", "") != ""
        ]
    if not rows:
        raise ValueError(f"run has no exact manufactured-solution errors: {directory}")
    return metadata, rows


def compare_prefix(
    baseline_directory: Path,
    extended_directory: Path,
    relative_tolerance: float = 5.0e-12,
    absolute_tolerance: float = 1.0e-13,
) -> dict:
    baseline_metadata, baseline = load_run(baseline_directory)
    extended_metadata, extended = load_run(extended_directory)
    identity_fields = ("case", "method")
    for field in identity_fields:
        if baseline_metadata.get(field) != extended_metadata.get(field):
            raise ValueError(f"run identity field {field!r} differs")
    for field in ("wavenumber", "boundary_beta", "initial_coarse_level"):
        if baseline_metadata["config"].get(field) != extended_metadata["config"].get(
            field
        ):
            raise ValueError(f"configuration field {field!r} differs")
    if baseline_metadata["config"].get("quadrature") != extended_metadata[
        "config"
    ].get("quadrature"):
        raise ValueError("quadrature policies differ")
    if len(extended) < len(baseline):
        raise ValueError("extended trajectory is shorter than the baseline")

    maximum_relative_difference = {field: 0.0 for field in EXACT_FIELDS}
    for index, baseline_row in enumerate(baseline):
        extended_row = extended[index]
        for field in DISCRETE_FIELDS:
            if baseline_row.get(field) != extended_row.get(field):
                raise ValueError(
                    f"prefix row {index} differs in {field}: "
                    f"{baseline_row.get(field)!r} != {extended_row.get(field)!r}"
                )
        for field in EXACT_FIELDS:
            first = float(baseline_row[field])
            second = float(extended_row[field])
            if not math.isclose(
                first,
                second,
                rel_tol=relative_tolerance,
                abs_tol=absolute_tolerance,
            ):
                raise ValueError(
                    f"prefix row {index} differs in {field}: {first:.17g} != "
                    f"{second:.17g}"
                )
            maximum_relative_difference[field] = max(
                maximum_relative_difference[field],
                abs(first - second) / max(absolute_tolerance, abs(first)),
            )
    return {
        "status": "prefix_reproduced",
        "case": baseline_metadata["case"],
        "method": baseline_metadata["method"],
        "baseline_points": len(baseline),
        "extended_points": len(extended),
        "relative_tolerance": relative_tolerance,
        "absolute_tolerance": absolute_tolerance,
        "maximum_relative_difference": maximum_relative_difference,
        "baseline_run": str(baseline_directory.resolve()),
        "extended_run": str(extended_directory.resolve()),
    }


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--extended", type=Path, required=True)
    parser.add_argument("--relative-tolerance", type=float, default=5.0e-12)
    parser.add_argument("--absolute-tolerance", type=float, default=1.0e-13)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    result = compare_prefix(
        arguments.baseline,
        arguments.extended,
        arguments.relative_tolerance,
        arguments.absolute_tolerance,
    )
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.output is not None:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(encoded, encoding="utf-8")
    print(encoded, end="")


if __name__ == "__main__":
    main()

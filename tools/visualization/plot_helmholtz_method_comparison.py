#!/usr/bin/env python3
"""Plot manufactured-solution error versus DoF for the S comparison methods.

All error values are post-processing quantities.  Comparator trajectories are
truncated at their first evaluated point at or below the terminal PALOD error;
they are never used to alter MARK or STOP decisions in the benchmark.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, ScalarFormatter

from paper_style import apply_paper_style


Row = Dict[str, str]


@dataclass(frozen=True)
class MethodRun:
    directory: Path
    provenance: str
    case: str
    method: str
    wavenumber: float
    boundary_beta: float
    initial_coarse_level: int
    manuscript_sha256: str
    quadrature_fingerprint: str
    status: str
    driver_state: str
    rows: List[Row]


def _number(row: Mapping[str, str], field: str) -> float:
    try:
        value = float(row[field])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"invalid numeric field {field!r} in row {row}") from error
    if not math.isfinite(value):
        raise ValueError(f"non-finite numeric field {field!r}")
    return value


def _integer(row: Mapping[str, str], field: str) -> int:
    value = _number(row, field)
    rounded = round(value)
    if not math.isclose(value, rounded, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValueError(f"non-integral field {field!r}: {value}")
    return int(rounded)


def load_method_run(directory: Path, provenance: Optional[str] = None) -> MethodRun:
    with (directory / "run.json").open("r", encoding="utf-8") as stream:
        metadata = json.load(stream)
    if metadata.get("status") != "success" or metadata.get("driver_state") != "TrajectoryComplete":
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
        raise ValueError(f"run has no manufactured exact-error points: {directory}")
    for row in rows:
        if _integer(row, "DoF_H") <= 0 or _number(
            row, "relative_exact_energy_error"
        ) <= 0.0:
            raise ValueError(f"invalid exact-error observation in {directory}")
    config = metadata["config"]
    return MethodRun(
        directory=directory,
        provenance=provenance or str(directory.resolve()),
        case=str(metadata["case"]),
        method=str(metadata["method"]),
        wavenumber=float(config["wavenumber"]),
        boundary_beta=float(config["boundary_beta"]),
        initial_coarse_level=int(config["initial_coarse_level"]),
        manuscript_sha256=str(config["manuscript_sha256"]),
        quadrature_fingerprint=json.dumps(
            config["quadrature"], sort_keys=True, separators=(",", ":")
        ),
        status=str(metadata["status"]),
        driver_state=str(metadata["driver_state"]),
        rows=rows,
    )


def truncate_at_target(run: MethodRun, target: float) -> tuple[List[Row], bool]:
    selected: List[Row] = []
    for row in run.rows:
        selected.append(row)
        if _number(row, "relative_exact_energy_error") <= target:
            return selected, True
    return selected, False


def validate_comparison(runs: Sequence[MethodRun]) -> None:
    if [run.method for run in runs] != ["PALOD", "SLOD", "UFEM", "AFEM"]:
        raise ValueError("runs must be supplied as PALOD, SLOD, UFEM, AFEM")
    identity = (
        runs[0].case,
        runs[0].wavenumber,
        runs[0].boundary_beta,
        runs[0].initial_coarse_level,
        runs[0].manuscript_sha256,
        runs[0].quadrature_fingerprint,
    )
    for run in runs[1:]:
        candidate = (
            run.case,
            run.wavenumber,
            run.boundary_beta,
            run.initial_coarse_level,
            run.manuscript_sha256,
            run.quadrature_fingerprint,
        )
        if candidate != identity:
            raise ValueError(
                "all methods must share the case, PDE, manuscript, quadrature, "
                "and initial coarse level"
            )
    initial_dofs = {_integer(run.rows[0], "DoF_H") for run in runs}
    if len(initial_dofs) != 1:
        raise ValueError("all methods must share the initial coarse DoF")
    palod_epochs = sorted(
        {_integer(row, "reference_epoch") for row in runs[0].rows}
    )
    if palod_epochs != [0, 1, 2]:
        raise ValueError("PALOD comparison must contain reference epochs 0, 1, and 2")


def _write_plot_data(
    runs: Sequence[MethodRun], selected: Mapping[str, Sequence[Row]], target: float, path: Path
) -> None:
    fields = [
        "case", "method", "kappa", "reference_epoch", "N_H", "DoF_H",
        "relative_exact_energy_error", "palod_terminal_target", "source_run_directory",
    ]
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        by_method = {run.method: run for run in runs}
        for method in ("PALOD", "SLOD", "UFEM", "AFEM"):
            run = by_method[method]
            for row in selected[method]:
                writer.writerow(
                    {
                        "case": run.case,
                        "method": method,
                        "kappa": f"{run.wavenumber:.17g}",
                        "reference_epoch": row["reference_epoch"],
                        "N_H": row["N_H"],
                        "DoF_H": row["DoF_H"],
                        "relative_exact_energy_error": row["relative_exact_energy_error"],
                        "palod_terminal_target": f"{target:.17g}",
                        "source_run_directory": run.provenance,
                    }
                )


def plot_method_comparison(
    runs: Sequence[MethodRun], output: Path, title: Optional[str] = None
) -> dict[str, object]:
    validate_comparison(runs)
    palod = runs[0]
    target = _number(palod.rows[-1], "relative_exact_energy_error")
    selected: Dict[str, List[Row]] = {"PALOD": list(palod.rows)}
    reached: Dict[str, bool] = {"PALOD": True}
    for run in runs[1:]:
        selected[run.method], reached[run.method] = truncate_at_target(run, target)

    apply_paper_style()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure, ax = plt.subplots(figsize=(5.4, 3.65))
    styles = {
        "PALOD": dict(marker="o", linestyle="-", linewidth=1.5),
        "SLOD": dict(marker="s", linestyle="--", linewidth=1.3),
        "UFEM": dict(marker="^", linestyle="-.", linewidth=1.2),
        "AFEM": dict(marker="D", linestyle=":", linewidth=1.3),
    }
    labels = {
        "PALOD": "PALOD (three reference epochs)",
        "SLOD": r"SLOD ($\ell=2$, fixed fine/coarse level gap)",
        "UFEM": "standard P1 FEM",
        "AFEM": "adaptive P1 FEM",
    }
    for run in runs:
        rows = selected[run.method]
        dofs = [_integer(row, "DoF_H") for row in rows]
        errors = [_number(row, "relative_exact_energy_error") for row in rows]
        suffix = "" if reached[run.method] else " (target not reached)"
        ax.plot(dofs, errors, label=labels[run.method] + suffix, **styles[run.method])

    seen_epochs = set()
    for row in selected["PALOD"]:
        epoch = _integer(row, "reference_epoch")
        if epoch in seen_epochs:
            continue
        seen_epochs.add(epoch)
        x = _integer(row, "DoF_H")
        y = _number(row, "relative_exact_energy_error")
        ax.scatter([x], [y], marker="*", s=80, edgecolors="black", linewidths=0.4, zorder=8)
        ax.annotate(
            f"epoch {epoch}", xy=(x, y), xytext=(5, 8), textcoords="offset points",
            fontsize=7.3, arrowprops={"arrowstyle": "-", "linewidth": 0.5},
        )

    all_dofs = [_integer(row, "DoF_H") for rows in selected.values() for row in rows]
    slope_x0 = max(min(all_dofs), min(max(all_dofs) / 8.0, max(all_dofs)))
    slope_x1 = min(max(all_dofs), slope_x0 * 4.0)
    anchor = target * 1.65
    slope_y1 = anchor * math.sqrt(slope_x0 / slope_x1)
    ax.plot([slope_x0, slope_x1], [anchor, slope_y1], color="black", linewidth=1.0, label=r"optimal $N^{-1/2}$")
    ax.axhline(target, color="0.35", linewidth=0.8, linestyle="--")
    ax.annotate(
        f"PALOD target {target:.3g}", xy=(max(all_dofs), target),
        xytext=(-4, 4), textcoords="offset points", ha="right", fontsize=7.3,
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Unconstrained coarse degrees of freedom")
    ax.set_ylabel("Relative weighted energy error vs exact solution")
    ticks = sorted(set(all_dofs))
    if len(ticks) > 5:
        ticks = [ticks[0], ticks[len(ticks) // 4], ticks[len(ticks) // 2], ticks[3 * len(ticks) // 4], ticks[-1]]
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.grid(True, which="major")
    ax.grid(True, which="minor", alpha=0.12)
    ax.legend(loc="best", frameon=False)
    figure.suptitle(title or f"{palod.case}, $\\kappa={palod.wavenumber:g}$: method comparison")
    figure.tight_layout()
    figure.savefig(output)
    if output.suffix.lower() != ".pdf":
        figure.savefig(output.with_suffix(".pdf"))
    plt.close(figure)
    _write_plot_data(runs, selected, target, output.with_suffix(".csv"))
    first_target_point = {}
    for method, rows in selected.items():
        if not reached[method]:
            first_target_point[method] = None
            continue
        row = rows[-1]
        first_target_point[method] = {
            "reference_epoch": _integer(row, "reference_epoch"),
            "iteration": _integer(row, "iteration"),
            "N_H": _integer(row, "N_H"),
            "DoF_H": _integer(row, "DoF_H"),
            "relative_exact_energy_error": _number(
                row, "relative_exact_energy_error"
            ),
        }
    summary = {
        "palod_terminal_target": target,
        "target_reached": reached,
        "first_target_point": first_target_point,
    }
    output.with_suffix(".json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    for method in ("palod", "slod", "ufem", "afem"):
        parser.add_argument(f"--{method}", type=Path, required=True)
        parser.add_argument(
            f"--{method}-source-label",
            help="optional canonical provenance path when plotting a copied run",
        )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--title")
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    runs = [
        load_method_run(
            getattr(arguments, method), getattr(arguments, f"{method}_source_label")
        )
        for method in ("palod", "slod", "ufem", "afem")
    ]
    summary = plot_method_comparison(runs, arguments.output, arguments.title)
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()

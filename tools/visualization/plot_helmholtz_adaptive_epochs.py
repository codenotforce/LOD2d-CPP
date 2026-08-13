#!/usr/bin/env python3
"""Plot practical adaptive-paper errors against true coarse degrees of freedom.

The plotter is deliberately post-processing only.  It reads completed run
directories, never invokes a solver, and marks the first evaluated point of
every explicit reference epoch.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter, ScalarFormatter

from paper_style import apply_paper_style


Row = Dict[str, str]


@dataclass(frozen=True)
class EpochRun:
    directory: Path
    case: str
    method: str
    wavenumber: float
    epoch: int
    status: str
    rows: List[Row]


def _number(row: Mapping[str, str], field: str) -> float:
    try:
        value = float(row[field])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"invalid numeric field {field!r} in row {row}") from error
    if not math.isfinite(value):
        raise ValueError(f"non-finite field {field!r} in row {row}")
    return value


def _integer(row: Mapping[str, str], field: str) -> int:
    value = _number(row, field)
    rounded = round(value)
    if not math.isclose(value, rounded, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValueError(f"field {field!r} is not integral: {value}")
    return int(rounded)


def _optional_number(row: Mapping[str, str], field: str) -> Optional[float]:
    text = row.get(field, "")
    if text == "":
        return None
    return _number(row, field)


def load_epoch_runs_from_directory(directory: Path) -> List[EpochRun]:
    run_path = directory / "run.json"
    iterations_path = directory / "iterations.csv"
    if not run_path.is_file() or not iterations_path.is_file():
        raise FileNotFoundError(
            f"run directory must contain run.json and iterations.csv: {directory}"
        )
    with run_path.open("r", encoding="utf-8") as stream:
        metadata = json.load(stream)
    with iterations_path.open("r", encoding="utf-8", newline="") as stream:
        all_rows = list(csv.DictReader(stream))

    config = metadata.get("config", {})
    rows = [
        row
        for row in all_rows
        if row.get("reference_energy_error", "") != ""
    ]
    if not rows:
        raise ValueError(f"run contains no evaluated error points: {directory}")
    initial_epoch = int(config["reference_epoch"])
    epochs = sorted({_integer(row, "reference_epoch") for row in rows})
    if not epochs or epochs[0] != initial_epoch or epochs != list(
        range(initial_epoch, initial_epoch + len(epochs))
    ):
        raise ValueError(
            f"evaluated epochs are not contiguous from run.json in {directory}"
        )
    for row in rows:
        dofs = _integer(row, "DoF_H")
        nodes = _integer(row, "N_H")
        if dofs <= 0 or dofs > nodes:
            raise ValueError(f"invalid coarse DoF/node counts in {directory}")
        for field in ("reference_energy_error", "reference_L2_error"):
            if not (_number(row, field) > 0.0):
                raise ValueError(f"non-positive {field} in {directory}")

    case = str(metadata["case"])
    exact_case = case in {"R1", "S"}
    exact_fields = (
        "exact_energy_error",
        "exact_L2_error",
        "relative_exact_energy_error",
        "relative_exact_L2_error",
    )
    for row in rows:
        present = [_optional_number(row, field) is not None for field in exact_fields]
        if exact_case and not all(present):
            raise ValueError(
                f"manufactured case {case} is missing exact errors in {directory}"
            )
        if not exact_case and any(present):
            raise ValueError(
                f"non-manufactured case {case} unexpectedly has exact errors"
            )

    return [
        EpochRun(
            directory=directory,
            case=case,
            method=str(metadata["method"]),
            wavenumber=float(config["wavenumber"]),
            epoch=epoch,
            status=str(metadata["status"]),
            rows=[row for row in rows if _integer(row, "reference_epoch") == epoch],
        )
        for epoch in epochs
    ]


def load_epoch_run(directory: Path) -> EpochRun:
    """Load a legacy single-epoch run, rejecting a continuous trajectory."""
    runs = load_epoch_runs_from_directory(directory)
    if len(runs) != 1:
        raise ValueError(f"run contains multiple reference epochs: {directory}")
    return runs[0]


def load_epoch_runs(directories: Iterable[Path]) -> List[EpochRun]:
    runs = sorted(
        (
            run
            for path in directories
            for run in load_epoch_runs_from_directory(path)
        ),
        key=lambda run: run.epoch,
    )
    if not runs:
        raise ValueError("at least one run directory is required")
    identity = (runs[0].case, runs[0].method, runs[0].wavenumber)
    epochs = set()
    for run in runs:
        if (run.case, run.method, run.wavenumber) != identity:
            raise ValueError("all epoch runs must share case, method, and wavenumber")
        if run.epoch in epochs:
            raise ValueError(f"duplicate reference epoch {run.epoch}")
        epochs.add(run.epoch)
    return runs


def _write_plot_data(runs: Sequence[EpochRun], path: Path) -> None:
    fields = [
        "case",
        "method",
        "kappa",
        "reference_epoch",
        "N_H",
        "DoF_H",
        "reference_energy_error",
        "reference_L2_error",
        "exact_energy_error",
        "exact_L2_error",
        "relative_exact_energy_error",
        "relative_exact_L2_error",
        "source_run_directory",
    ]
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for run in runs:
            for row in run.rows:
                writer.writerow(
                    {
                        "case": run.case,
                        "method": run.method,
                        "kappa": f"{run.wavenumber:.17g}",
                        "reference_epoch": run.epoch,
                        "N_H": row["N_H"],
                        "DoF_H": row["DoF_H"],
                        "reference_energy_error": row["reference_energy_error"],
                        "reference_L2_error": row["reference_L2_error"],
                        "exact_energy_error": row.get("exact_energy_error", ""),
                        "exact_L2_error": row.get("exact_L2_error", ""),
                        "relative_exact_energy_error": row.get(
                            "relative_exact_energy_error", ""
                        ),
                        "relative_exact_L2_error": row.get(
                            "relative_exact_L2_error", ""
                        ),
                        "source_run_directory": str(run.directory.resolve()),
                    }
                )


def plot_epoch_errors(
    runs: Sequence[EpochRun], output: Path, title: Optional[str] = None
) -> Path:
    if not runs:
        raise ValueError("cannot plot an empty run sequence")
    apply_paper_style()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(1, 2, figsize=(7.25, 3.25), sharex=False)
    colors = plt.get_cmap("tab10")
    panels = (
        (
            "reference_energy_error",
            "relative_exact_energy_error",
            "Relative weighted energy error",
        ),
        ("reference_L2_error", "relative_exact_L2_error", r"Relative $L^2$ error"),
    )

    has_exact = any(
        row.get("relative_exact_energy_error", "") != ""
        for run in runs
        for row in run.rows
    )
    for panel_index, (reference_field, exact_field, ylabel) in enumerate(panels):
        ax = axes[panel_index]
        for run_index, run in enumerate(runs):
            color = colors(run_index % 10)
            dofs = [_integer(row, "DoF_H") for row in run.rows]
            reference = [_number(row, reference_field) for row in run.rows]
            ax.plot(dofs, reference, color=color, marker="o", linestyle="-")
            exact = [_optional_number(row, exact_field) for row in run.rows]
            if all(value is not None for value in exact):
                exact_values = [float(value) for value in exact if value is not None]
                ax.plot(
                    dofs,
                    exact_values,
                    color=color,
                    marker="^",
                    linestyle="--",
                )
                start_y = exact_values[0]
            else:
                start_y = reference[0]
            ax.scatter(
                [dofs[0]],
                [start_y],
                color=color,
                marker="*",
                s=75,
                edgecolors="black",
                linewidths=0.45,
                zorder=6,
            )
            ax.annotate(
                f"epoch {run.epoch} start",
                xy=(dofs[0], start_y),
                xytext=(5, 8 + 11 * (run_index % 3)),
                textcoords="offset points",
                color=color,
                fontsize=7.2,
                arrowprops={"arrowstyle": "-", "color": color, "linewidth": 0.55},
            )

        ax.set_xscale("log")
        ax.set_yscale("log")
        all_dofs = sorted(
            {_integer(row, "DoF_H") for run in runs for row in run.rows}
        )
        if len(all_dofs) > 3:
            tick_dofs = [all_dofs[0], all_dofs[len(all_dofs) // 2], all_dofs[-1]]
        else:
            tick_dofs = all_dofs
        ax.set_xticks(tick_dofs)
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.grid(True, which="major")
        ax.grid(True, which="minor", alpha=0.12)
        ax.set_xlabel(r"Unconstrained coarse degrees of freedom")
        ax.set_ylabel(ylabel)

    style_handles = [
        Line2D([], [], color="black", marker="o", linestyle="-", label="vs reference FEM")
    ]
    if has_exact:
        style_handles.append(
            Line2D(
                [], [], color="black", marker="^", linestyle="--",
                label="vs manufactured exact solution"
            )
        )
    epoch_handles = [
        Line2D(
            [], [], color=colors(index % 10), marker="*", linestyle="none",
            label=f"reference epoch {run.epoch}"
        )
        for index, run in enumerate(runs)
    ]
    axes[0].legend(handles=style_handles + epoch_handles, loc="best", frameon=False)
    identity = runs[0]
    figure.suptitle(
        title
        or f"{identity.case} {identity.method}, $\\kappa={identity.wavenumber:g}$: error vs DoF"
    )
    figure.text(
        0.5,
        0.005,
        "Exploratory trajectory; reference adequacy is not frozen for paper use.",
        ha="center",
        va="bottom",
        fontsize=7.2,
    )
    figure.tight_layout(rect=(0.0, 0.04, 1.0, 0.95))
    figure.savefig(output)
    plt.close(figure)
    _write_plot_data(runs, output.with_suffix(".csv"))
    return output


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-dir",
        action="append",
        type=Path,
        required=True,
        help="completed adaptive-paper run directory; repeat for multiple epochs",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--title")
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    runs = load_epoch_runs(arguments.run_dir)
    path = plot_epoch_errors(runs, arguments.output, arguments.title)
    print(path)


if __name__ == "__main__":
    main()

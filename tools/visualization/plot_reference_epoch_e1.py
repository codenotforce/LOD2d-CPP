#!/usr/bin/env python3
"""Plot E1 R1 energy errors for mixed paper-runner schema v4/v5 outputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from paper_style import apply_paper_style


@dataclass(frozen=True)
class Observation:
    dofs: int
    seconds: float
    exact: float
    reference: float | None
    epoch: int
    ell: int | None
    skipped_correctors: int
    skipped_work_units: int


@dataclass(frozen=True)
class Run:
    label: str
    directory: Path
    state: str
    observations: list[Observation]


def _positive(row: dict[str, str], field: str) -> float | None:
    text = row.get(field, "")
    if text in ("", "NA"):
        return None
    value = float(text)
    return value if math.isfinite(value) and value > 0.0 else None


def load_run(directory: Path, label: str) -> Run:
    metadata = json.loads((directory / "run.json").read_text(encoding="utf-8"))
    schema = int(metadata["schema_version"])
    observations: list[Observation] = []
    last_skipped_correctors = 0
    last_skipped_work_units = 0
    with (directory / "iterations.csv").open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            if row.get("action") in ("AcceptCorrector", "IncreaseGlobalEll"):
                last_skipped_correctors = int(row.get("skipped_correctors", "0") or 0)
                last_skipped_work_units = int(
                    row.get("skipped_corrector_work_units", "0") or 0
                )
            if schema == 5:
                exact = _positive(row, "relative_exact_energy")
                reference = _positive(row, "relative_reference_energy")
                dofs = _positive(row, "N_H")
                epoch_text = row.get("epoch", "0")
            else:
                exact = _positive(row, "relative_exact_energy_error")
                reference = _positive(row, "reference_energy_error")
                dofs = _positive(row, "DoF_H")
                epoch_text = row.get("reference_epoch", "0")
            seconds = _positive(row, "time_total_cumulative")
            if exact is None or dofs is None or seconds is None:
                continue
            observations.append(
                Observation(
                    int(round(dofs)),
                    seconds,
                    exact,
                    reference,
                    int(epoch_text or 0),
                    int(row["ell"]) if row.get("ell", "") not in ("", "NA") else None,
                    last_skipped_correctors,
                    last_skipped_work_units,
                )
            )
    if not observations:
        raise ValueError(f"no complete positive error observations in {directory}")
    return Run(
        label,
        directory,
        str(metadata.get("driver_state", metadata.get("status", "unknown"))),
        observations,
    )


def _plot_panel(ax, runs: Iterable[Run], x_field: str, y_field: str) -> None:
    styles = {
        "Reference-epoch PALOD": dict(
            color="#1f77b4", marker="o", linestyle="-", linewidth=1.45
        ),
        r"Fixed LOD ($\ell=3$)": dict(
            color="#ff7f0e", marker="s", linestyle="--", linewidth=1.25
        ),
        "Hybrid reference-epoch PALOD": dict(
            color="#1f77b4", marker="o", linestyle="-", linewidth=1.45
        ),
        "Standard reference-epoch PALOD": dict(
            color="#9467bd", marker="v", linestyle="--", linewidth=1.3
        ),
        "UFEM": dict(
            color="#2ca02c", marker="^", linestyle="-.", linewidth=1.2
        ),
        "AFEM": dict(
            color="#d62728", marker="D", linestyle=":", linewidth=1.25
        ),
    }
    for run in runs:
        points = [
            point
            for point in run.observations
            if getattr(point, y_field) is not None
        ]
        if not points:
            continue
        x = [getattr(point, x_field) for point in points]
        y = [getattr(point, y_field) for point in points]
        segments = [points]
        if "PALOD" in run.label and y_field == "reference":
            segments = []
            for epoch in sorted({point.epoch for point in points}):
                segments.append([point for point in points if point.epoch == epoch])
        for segment_index, segment in enumerate(segments):
            ax.loglog(
                [getattr(point, x_field) for point in segment],
                [getattr(point, y_field) for point in segment],
                label=run.label if segment_index == 0 else None,
                markersize=3.5,
                markerfacecolor="white",
                markeredgewidth=0.8,
                **styles[run.label],
            )
        if "PALOD" in run.label:
            previous_epoch = points[0].epoch
            previous_ell = points[0].ell
            for point in points[1:]:
                if point.epoch != previous_epoch:
                    ax.scatter(
                        [getattr(point, x_field)],
                        [getattr(point, y_field)],
                        marker="*",
                        s=48,
                        zorder=5,
                        label=None,
                    )
                    previous_epoch = point.epoch
                if point.ell != previous_ell:
                    ax.scatter(
                        [getattr(point, x_field)],
                        [getattr(point, y_field)],
                        marker="x",
                        s=30,
                        linewidths=1.1,
                        zorder=6,
                        color="black",
                        label=None,
                    )
                    if x_field == "dofs" and y_field == "exact":
                        ax.annotate(
                            rf"$\ell={point.ell}$",
                            (getattr(point, x_field), getattr(point, y_field)),
                            xytext=(4, 5),
                            textcoords="offset points",
                            fontsize=6.5,
                        )
                    previous_ell = point.ell
    ax.grid(True, which="major", linewidth=0.45, alpha=0.45)
    ax.grid(True, which="minor", linewidth=0.25, alpha=0.2)


def tail_dof_exponent(run: Run, field: str, count: int = 4) -> float | None:
    observations = run.observations
    if "PALOD" in run.label:
        final_epoch = max(point.epoch for point in observations)
        observations = [point for point in observations if point.epoch == final_epoch]
    by_dof = {
        point.dofs: getattr(point, field)
        for point in observations
        if getattr(point, field) is not None
    }
    points = sorted(by_dof.items())[-count:]
    if len(points) < 2:
        return None
    x = [math.log(dofs) for dofs, _ in points]
    y = [math.log(value) for _, value in points]
    x_mean = sum(x) / len(x)
    y_mean = sum(y) / len(y)
    denominator = sum((value - x_mean) ** 2 for value in x)
    if denominator <= 0.0:
        return None
    slope = sum(
        (x_value - x_mean) * (y_value - y_mean)
        for x_value, y_value in zip(x, y)
    ) / denominator
    return -slope


def palod_epoch_dof_exponents(run: Run, field: str) -> dict[str, float | None]:
    if "PALOD" not in run.label:
        return {}
    result: dict[str, float | None] = {}
    for epoch in sorted({point.epoch for point in run.observations}):
        points = [
            point
            for point in run.observations
            if point.epoch == epoch and getattr(point, field) is not None
        ]
        by_dof = {
            point.dofs: getattr(point, field)
            for point in points
        }
        ordered = sorted(by_dof.items())
        if len(ordered) < 2:
            result[str(epoch)] = None
            continue
        x = [math.log(dofs) for dofs, _ in ordered]
        y = [math.log(value) for _, value in ordered]
        x_mean = sum(x) / len(x)
        y_mean = sum(y) / len(y)
        denominator = sum((value - x_mean) ** 2 for value in x)
        result[str(epoch)] = (
            None
            if denominator <= 0.0
            else -sum(
                (x_value - x_mean) * (y_value - y_mean)
                for x_value, y_value in zip(x, y)
            ) / denominator
        )
    return result


def plot_hybrid_saved_work(ax, runs: Iterable[Run]) -> None:
    hybrid = next(run for run in runs if run.label.startswith("Hybrid"))
    ax.semilogx(
        [point.dofs for point in hybrid.observations],
        [point.skipped_work_units for point in hybrid.observations],
        color="#1f77b4", marker="o", markersize=3.5,
        markerfacecolor="white", linewidth=1.35,
    )
    ax.set_xlabel("DoF")
    ax.set_ylabel("Skipped corrector work units")
    ax.grid(True, which="major", linewidth=0.45, alpha=0.45)
    ax.grid(True, which="minor", linewidth=0.25, alpha=0.2)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", choices=("E1", "E2"), default="E1")
    parser.add_argument("--palod", type=Path)
    parser.add_argument("--hybrid", type=Path)
    parser.add_argument("--standard", type=Path)
    parser.add_argument("--fixed-lod", type=Path)
    parser.add_argument("--ufem", type=Path)
    parser.add_argument("--afem", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    if arguments.experiment == "E1":
        if None in (arguments.palod, arguments.fixed_lod, arguments.ufem, arguments.afem):
            parser.error("E1 requires --palod, --fixed-lod, --ufem, and --afem")
        runs = [
            load_run(arguments.palod, "Reference-epoch PALOD"),
            load_run(arguments.fixed_lod, r"Fixed LOD ($\ell=3$)"),
            load_run(arguments.ufem, "UFEM"),
            load_run(arguments.afem, "AFEM"),
        ]
    else:
        if None in (arguments.hybrid, arguments.standard, arguments.afem):
            parser.error("E2 requires --hybrid, --standard, and --afem")
        runs = [
            load_run(arguments.hybrid, "Hybrid reference-epoch PALOD"),
            load_run(arguments.standard, "Standard reference-epoch PALOD"),
            load_run(arguments.afem, "AFEM"),
        ]
        if arguments.fixed_lod is not None:
            runs.insert(2, load_run(arguments.fixed_lod, r"Fixed LOD ($\ell=3$)"))
    apply_paper_style()
    figure, axes = plt.subplots(2, 2, figsize=(8.2, 6.4), constrained_layout=True)
    panels = [
        (axes[0, 0], "dofs", "exact", "DoF", "Relative exact energy error"),
        (axes[0, 1], "seconds", "exact", "Cumulative wall time [s]", "Relative exact energy error"),
        (axes[1, 0], "dofs", "reference", "DoF", "Relative reference energy error"),
        (axes[1, 1], "seconds", "reference", "Cumulative wall time [s]", "Relative reference energy error"),
    ]
    for ax, x_field, y_field, xlabel, ylabel in panels:
        _plot_panel(ax, runs, x_field, y_field)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
    if arguments.experiment == "E2":
        axes[1, 1].clear()
        plot_hybrid_saved_work(axes[1, 1], runs)
    axes[0, 0].legend(loc="best", fontsize=7.0, framealpha=0.82, handlelength=2.2)
    figure.suptitle(
        r"E1 (R1, $\kappa=16$, initial $H$ level 2)"
        if arguments.experiment == "E1"
        else r"E2 (L-shaped mixed boundary, $\kappa=16$, initial $H$ level 3)"
    )
    figure.text(
        0.5,
        0.005,
        r"PALOD: $\star$ reference refresh; $\times$ oversampling change",
        ha="center",
        fontsize=7.0,
    )
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(arguments.output.with_suffix(".png"), dpi=240)
    figure.savefig(arguments.output.with_suffix(".pdf"))
    plt.close(figure)

    summary = {
        run.label: {
            "directory": str(run.directory),
            "driver_state": run.state,
            "points": len(run.observations),
            "terminal_dofs": run.observations[-1].dofs,
            "terminal_seconds": run.observations[-1].seconds,
            "terminal_relative_exact_energy": run.observations[-1].exact,
            "terminal_relative_reference_energy": next(
                (
                    point.reference
                    for point in reversed(run.observations)
                    if point.reference is not None
                ),
                None,
            ),
            "epochs": sorted({point.epoch for point in run.observations}),
            "tail4_exact_dof_exponent": tail_dof_exponent(run, "exact"),
            "tail4_reference_dof_exponent": tail_dof_exponent(run, "reference"),
            "epoch_exact_dof_exponents": palod_epoch_dof_exponents(
                run, "exact"
            ),
            "epoch_reference_dof_exponents": palod_epoch_dof_exponents(
                run, "reference"
            ),
        }
        for run in runs
    }
    arguments.output.with_suffix(".json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()

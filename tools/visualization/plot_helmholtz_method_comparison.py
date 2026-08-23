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
from matplotlib.ticker import FuncFormatter, NullFormatter

from paper_style import apply_paper_style


Row = Dict[str, str]


def _format_dof_tick(value: float, _position: int) -> str:
    if value >= 1.0e6:
        return f"{value / 1.0e6:g}M"
    if value >= 1.0e3:
        return f"{value / 1.0e3:g}k"
    return f"{value:g}"


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
    initial_ell: int
    rows: List[Row]
    ell_changes: List[Row]


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
    with (directory / "ell_history.csv").open(
        "r", encoding="utf-8", newline=""
    ) as stream:
        ell_changes = [
            row for row in csv.DictReader(stream)
            if row.get("action") == "IncreaseGlobalEll"
        ]
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
        initial_ell=int(config["ell0"]),
        rows=rows,
        ell_changes=ell_changes,
    )


def truncate_at_target(run: MethodRun, target: float) -> tuple[List[Row], bool]:
    selected: List[Row] = []
    for row in run.rows:
        selected.append(row)
        if _number(row, "relative_exact_energy_error") <= target:
            return selected, True
    return selected, False


def fit_dof_rate(
    rows: Sequence[Row], tail_points: Optional[int] = None
) -> dict[str, float | int]:
    """Fit error = C * DoF**(-p), keeping the last value at repeated DoF."""
    by_dof: Dict[int, float] = {}
    for row in rows:
        by_dof[_integer(row, "DoF_H")] = _number(
            row, "relative_exact_energy_error"
        )
    observations = sorted(by_dof.items())
    if tail_points is not None:
        if tail_points < 2:
            raise ValueError("a tail rate fit requires at least two points")
        observations = observations[-tail_points:]
    if len(observations) < 2:
        raise ValueError("at least two distinct DoF values are required for a rate fit")
    x = [math.log(dof) for dof, _ in observations]
    y = [math.log(error) for _, error in observations]
    x_mean = sum(x) / len(x)
    y_mean = sum(y) / len(y)
    denominator = sum((value - x_mean) ** 2 for value in x)
    if denominator <= 0.0:
        raise ValueError("degenerate DoF range in rate fit")
    slope = sum(
        (x_value - x_mean) * (y_value - y_mean)
        for x_value, y_value in zip(x, y)
    ) / denominator
    fitted = [y_mean + slope * (value - x_mean) for value in x]
    residual = sum((value - estimate) ** 2 for value, estimate in zip(y, fitted))
    total = sum((value - y_mean) ** 2 for value in y)
    r_squared = 1.0 - residual / total if total > 0.0 else 1.0
    return {
        "exponent": -slope,
        "r_squared": r_squared,
        "distinct_dof_points": len(observations),
    }


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
    if palod_epochs != list(range(len(palod_epochs))):
        raise ValueError(
            "PALOD reference epochs must be consecutive and start at zero"
        )


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
    runs: Sequence[MethodRun],
    output: Path,
    title: Optional[str] = None,
    truncate_to_palod_target: bool = True,
    legend_tail_points: Optional[int] = None,
) -> dict[str, object]:
    validate_comparison(runs)
    palod = runs[0]
    target = _number(palod.rows[-1], "relative_exact_energy_error")
    selected: Dict[str, List[Row]] = {"PALOD": list(palod.rows)}
    reached: Dict[str, bool] = {"PALOD": True}
    for run in runs[1:]:
        target_rows, reached[run.method] = truncate_at_target(run, target)
        selected[run.method] = (
            target_rows if truncate_to_palod_target else list(run.rows)
        )
    rates = {method: fit_dof_rate(rows) for method, rows in selected.items()}
    tail_rates = {
        method: fit_dof_rate(rows, tail_points=3)
        for method, rows in selected.items()
    }
    if legend_tail_points is not None and legend_tail_points < 2:
        raise ValueError("legend_tail_points must be at least two")
    legend_rates = (
        {
            method: fit_dof_rate(rows, tail_points=legend_tail_points)
            for method, rows in selected.items()
        }
        if legend_tail_points is not None
        else rates
    )

    apply_paper_style()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure, ax = plt.subplots(figsize=(5.6, 4.65))
    styles = {
        "PALOD": dict(marker="o", linestyle="-", linewidth=1.5),
        "SLOD": dict(marker="s", linestyle="--", linewidth=1.3),
        "UFEM": dict(marker="^", linestyle="-.", linewidth=1.2),
        "AFEM": dict(marker="D", linestyle=":", linewidth=1.3),
    }
    labels = {
        "PALOD": "PALOD",
        "SLOD": r"SLOD ($\ell=2$)",
        "UFEM": "standard P1 FEM",
        "AFEM": "adaptive P1 FEM",
    }
    method_colors: Dict[str, str] = {}
    for run in runs:
        rows = selected[run.method]
        dofs = [_integer(row, "DoF_H") for row in rows]
        errors = [_number(row, "relative_exact_energy_error") for row in rows]
        rate = float(legend_rates[run.method]["exponent"])
        rate_name = (
            rf"p_{{\rm tail,{legend_tail_points}}}"
            if legend_tail_points is not None
            else "p"
        )
        suffix = "" if reached[run.method] else "; target not reached"
        line, = ax.plot(
            dofs, errors,
            label=labels[run.method]
            + rf" (${rate_name}\approx{rate:.2f}$)"
            + suffix,
            **styles[run.method],
        )
        method_colors[run.method] = line.get_color()

    ell_change_steps = {
        _integer(change, "H_step") for change in palod.ell_changes
    }
    seen_epochs = set()
    for row in selected["PALOD"]:
        epoch = _integer(row, "reference_epoch")
        if epoch in seen_epochs:
            continue
        seen_epochs.add(epoch)
        if _integer(row, "H_step") in ell_change_steps:
            continue
        x = _integer(row, "DoF_H")
        y = _number(row, "relative_exact_energy_error")
        ax.scatter([x], [y], marker="*", s=80, edgecolors="black", linewidths=0.4, zorder=8)
        epoch_offset = (6, 8) if epoch % 2 else (6, -15)
        ax.annotate(
            f"epoch {epoch}", xy=(x, y), xytext=epoch_offset,
            textcoords="offset points", fontsize=6.5,
            arrowprops={"arrowstyle": "-", "linewidth": 0.5},
        )

    palod_by_step: Dict[int, Row] = {}
    for row in selected["PALOD"]:
        palod_by_step[_integer(row, "H_step")] = row
    previous_ell = palod.initial_ell
    ell_change_summary: List[dict[str, int]] = []
    for change in palod.ell_changes:
        step = _integer(change, "H_step")
        new_ell = _integer(change, "ell")
        if step not in palod_by_step:
            raise ValueError(f"PALOD ell change at H-step {step} has no error point")
        point = palod_by_step[step]
        x = _integer(point, "DoF_H")
        y = _number(point, "relative_exact_energy_error")
        label = rf"$H$-step {step}: $\ell$ {previous_ell}$\to${new_ell}"
        ax.scatter(
            [x], [y], marker="P", s=62,
            facecolors=method_colors["PALOD"], edgecolors="black",
            linewidths=0.5, zorder=9,
        )
        ax.annotate(
            label, xy=(x, y), xytext=(18, 17), textcoords="offset points",
            fontsize=6.5, arrowprops={"arrowstyle": "->", "linewidth": 0.65},
        )
        ell_change_summary.append(
            {"H_step": step, "old_ell": previous_ell, "new_ell": new_ell}
        )
        previous_ell = new_ell

    all_dofs = [_integer(row, "DoF_H") for rows in selected.values() for row in rows]
    reference_slopes: Dict[str, dict[str, float | int | str]] = {}

    def add_reference_slope(method: str, exponent: float, label: str) -> None:
        # Anchor each guide at the method's terminal plotted point.  The guide
        # spans at most a factor four in DoF, so it is a local visual comparator
        # rather than a fitted or claimed asymptotic rate.
        by_dof = {
            _integer(row, "DoF_H"): _number(row, "relative_exact_energy_error")
            for row in selected[method]
        }
        anchor_dof, anchor_error = sorted(by_dof.items())[-1]
        first_dof = min(by_dof)
        left_dof = max(first_dof, anchor_dof / 4.0)
        left_error = anchor_error * (anchor_dof / left_dof) ** exponent
        ax.plot(
            [left_dof, anchor_dof], [left_error, anchor_error],
            color=method_colors[method], linewidth=1.0,
            linestyle=(0, (4, 2, 1, 2)), alpha=0.72, label=label,
        )
        reference_slopes[method] = {
            "exponent": exponent,
            "anchor_dof": anchor_dof,
            "anchor_error": anchor_error,
            "left_dof": left_dof,
            "left_error": left_error,
            "interpretation": "local visual guide anchored at the terminal plotted point",
        }

    add_reference_slope("PALOD", 0.5, r"$N^{-1/2}$, PALOD anchor")
    add_reference_slope("AFEM", 0.5, r"$N^{-1/2}$, AFEM anchor")
    add_reference_slope("SLOD", 1.0 / 3.0, r"$N^{-1/3}$, SLOD anchor")
    add_reference_slope("UFEM", 1.0 / 3.0, r"$N^{-1/3}$, standard FEM anchor")
    ax.axhline(target, color="0.35", linewidth=0.8, linestyle="--")
    ax.annotate(
        f"PALOD target {target:.3g}", xy=(max(all_dofs), target),
        xytext=(-4, -12), textcoords="offset points", ha="right", fontsize=7.3,
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Unconstrained coarse degrees of freedom")
    ax.set_ylabel("Relative weighted energy error vs exact solution")
    ticks = sorted(set(all_dofs))
    if len(ticks) > 5:
        ticks = [ticks[0], ticks[len(ticks) // 4], ticks[len(ticks) // 2], ticks[3 * len(ticks) // 4], ticks[-1]]
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(FuncFormatter(_format_dof_tick))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.grid(True, which="major")
    ax.grid(True, which="minor", alpha=0.12)
    handles, legend_labels = ax.get_legend_handles_labels()
    figure.legend(
        handles,
        legend_labels,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.015),
        frameon=False,
        ncol=2,
        fontsize=6.1,
        columnspacing=0.9,
        handlelength=2.0,
        handletextpad=0.45,
    )
    figure.suptitle(title or f"{palod.case}, $\\kappa={palod.wavenumber:g}$: method comparison")
    figure.tight_layout(rect=(0.0, 0.19, 1.0, 1.0))
    figure.savefig(output)
    if output.suffix.lower() != ".pdf":
        figure.savefig(output.with_suffix(".pdf"))
    plt.close(figure)
    _write_plot_data(runs, selected, target, output.with_suffix(".csv"))
    first_target_point = {}
    for method, rows in selected.items():
        crossing = next(
            (
                row
                for row in rows
                if _number(row, "relative_exact_energy_error") <= target
            ),
            None,
        )
        if crossing is None:
            first_target_point[method] = None
            continue
        first_target_point[method] = {
            "reference_epoch": _integer(crossing, "reference_epoch"),
            "iteration": _integer(crossing, "iteration"),
            "N_H": _integer(crossing, "N_H"),
            "DoF_H": _integer(crossing, "DoF_H"),
            "relative_exact_energy_error": _number(
                crossing, "relative_exact_energy_error"
            ),
        }
    summary = {
        "palod_terminal_target": target,
        "target_reached": reached,
        "trajectory_display_policy": (
            "truncate_at_palod_target"
            if truncate_to_palod_target
            else "full_fixed_horizon"
        ),
        "first_target_point": first_target_point,
        "empirical_dof_rate": rates,
        "palod_ell_changes": ell_change_summary,
        "empirical_dof_rate_last_three": tail_rates,
        "legend_dof_rate": legend_rates,
        "legend_rate_policy": (
            f"last_{legend_tail_points}_distinct_dof_points"
            if legend_tail_points is not None
            else "whole_displayed_range"
        ),
        "reference_slope_guides": reference_slopes,
        "palod_vs_optimal_half": {
            "empirical_exponent": rates["PALOD"]["exponent"],
            "last_three_exponent": tail_rates["PALOD"]["exponent"],
            "optimal_exponent": 0.5,
            "exponent_difference": rates["PALOD"]["exponent"] - 0.5,
            "last_three_exponent_difference": (
                tail_rates["PALOD"]["exponent"] - 0.5
            ),
            "interpretation": (
                "finite-range coarse-DoF fit; corrector and reference work excluded"
            ),
        },
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
    parser.add_argument(
        "--full-trajectories",
        action="store_true",
        help=(
            "plot every fixed-horizon point instead of truncating comparator "
            "curves at their first PALOD-target crossing"
        ),
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    runs = [
        load_method_run(
            getattr(arguments, method), getattr(arguments, f"{method}_source_label")
        )
        for method in ("palod", "slod", "ufem", "afem")
    ]
    summary = plot_method_comparison(
        runs,
        arguments.output,
        arguments.title,
        truncate_to_palod_target=not arguments.full_trajectories,
        legend_tail_points=4 if arguments.full_trajectories else None,
    )
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Plot the server Helmholtz H-convergence experiment."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

from paper_style import COLORS, apply_paper_style


Row = Dict[str, str]
Metric = Tuple[str, float]


def read_rows(path: Path) -> List[Row]:
    if not path.is_file():
        raise FileNotFoundError(f"missing H-convergence CSV: {path}")
    with path.open("r", encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) < 2:
        raise ValueError("H-convergence CSV must contain at least two rows")
    return rows


def number(row: Mapping[str, str], field: str) -> float:
    try:
        value = float(row[field])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"invalid numeric field {field!r}") from error
    if not math.isfinite(value):
        raise ValueError(f"non-finite required field {field!r}")
    return value


def integer(row: Mapping[str, str], field: str) -> int:
    value = number(row, field)
    rounded = round(value)
    if not math.isclose(value, rounded, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValueError(f"field {field!r} is not integral")
    return int(rounded)


def power_slope(x: Sequence[float], y: Sequence[float]) -> float:
    if len(x) != len(y) or len(x) < 2:
        raise ValueError("power fit needs equal sequences with at least two values")
    log_x = [math.log(value) for value in x]
    log_y = [math.log(value) for value in y]
    mean_x = sum(log_x) / len(log_x)
    mean_y = sum(log_y) / len(log_y)
    denominator = sum((value - mean_x) ** 2 for value in log_x)
    if denominator == 0.0:
        raise ValueError("power fit has zero x variance")
    return sum(
        (x_value - mean_x) * (y_value - mean_y)
        for x_value, y_value in zip(log_x, log_y)
    ) / denominator


def fitted_values(
    x: Sequence[float], y: Sequence[float], slope: float
) -> List[float]:
    intercept = sum(
        math.log(y_value) - slope * math.log(x_value)
        for x_value, y_value in zip(x, y)
    ) / len(x)
    return [math.exp(intercept) * value**slope for value in x]


def adjacent_order(
    previous_error: float, error: float, previous_h: float, h: float
) -> float:
    return math.log(previous_error / error) / math.log(previous_h / h)


def load_and_validate(path: Path) -> List[Row]:
    rows = sorted(read_rows(path), key=lambda row: integer(row, "coarse_nodes"))
    expected_levels = list(range(8, 14))
    if [integer(row, "H_level") for row in rows] != expected_levels:
        raise ValueError(
            "expected completed H levels 8,9,10,11,12,13 for the h=19 run"
        )
    for row in rows:
        if integer(row, "k") != 32:
            raise ValueError("server convergence row is not k=32")
        if integer(row, "h_level") != 19:
            raise ValueError("server convergence row is not h=19")
        if integer(row, "ell") != 4:
            raise ValueError("server convergence row is not ell=4")
        if row.get("solver") != "schur" or row.get("mode") != "two-sided":
            raise ValueError("unexpected solver or Petrov mode")
        for field in ("p1_energy_abs", "lod_energy_abs", "H_max"):
            if not number(row, field) > 0.0:
                raise ValueError(f"{field} must be positive")
        if number(row, "nesting_coordinate_residual") >= 1.0e-13:
            raise ValueError("nested-mesh coordinate check failed")
        for field in (
            "p1_residual",
            "petrov_residual",
            "corrector_residual",
            "constraint_residual",
            "schur_residual",
        ):
            if number(row, field) >= 1.0e-8:
                raise ValueError(f"{field} exceeds the runbook tolerance")
        if number(row, "schur_rcond") <= 1.0e-14:
            raise ValueError("Schur reciprocal condition estimate is too small")

    for field in ("p1_energy_abs", "lod_energy_abs"):
        values = [number(row, field) for row in rows]
        if any(current >= previous for previous, current in zip(values, values[1:])):
            raise ValueError(f"{field} is not strictly decreasing")
    if any(
        number(row, "lod_energy_abs") >= number(row, "p1_energy_abs")
        for row in rows
    ):
        raise ValueError("LOD is not more accurate than coarse P1 on every row")

    for index in range(1, len(rows)):
        previous = rows[index - 1]
        row = rows[index]
        for error_field, rate_field in (
            ("p1_energy_abs", "p1_energy_rate"),
            ("lod_energy_abs", "lod_energy_rate"),
        ):
            recomputed = adjacent_order(
                number(previous, error_field),
                number(row, error_field),
                number(previous, "H_max"),
                number(row, "H_max"),
            )
            if not math.isclose(
                recomputed,
                number(row, rate_field),
                rel_tol=2.0e-12,
                abs_tol=2.0e-12,
            ):
                raise ValueError(f"stored {rate_field} disagrees with errors/H_max")
    return rows


def plot_convergence(
    rows: Sequence[Row], output_dir: Path, formats: Sequence[str]
) -> Tuple[List[Path], List[Metric]]:
    apply_paper_style()
    dofs = [integer(row, "coarse_nodes") for row in rows]
    p1_error = [number(row, "p1_energy_abs") for row in rows]
    lod_error = [number(row, "lod_energy_abs") for row in rows]
    p1_slope = power_slope(dofs, p1_error)
    lod_slope = power_slope(dofs, lod_error)
    p1_orders = [number(row, "p1_energy_rate") for row in rows[1:]]
    lod_orders = [number(row, "lod_energy_rate") for row in rows[1:]]
    finer_dofs = dofs[1:]

    fig, axes = plt.subplots(
        1, 2, figsize=(7.25, 3.25), gridspec_kw={"width_ratios": (1.45, 1.0)}
    )
    ax = axes[0]
    ax.plot(
        dofs,
        p1_error,
        color=COLORS["fem"],
        marker="s",
        label=rf"Coarse P1 FEM ($s={p1_slope:.2f}$)",
        zorder=4,
    )
    ax.plot(
        dofs,
        lod_error,
        color=COLORS["lod"],
        marker="o",
        label=rf"Two-sided LOD ($s={lod_slope:.2f}$)",
        zorder=4,
    )
    ax.plot(
        dofs,
        fitted_values(dofs, p1_error, p1_slope),
        color=COLORS["fem"],
        linestyle=":",
        linewidth=1.0,
        alpha=0.8,
    )
    ax.plot(
        dofs,
        fitted_values(dofs, lod_error, lod_slope),
        color=COLORS["lod"],
        linestyle=":",
        linewidth=1.0,
        alpha=0.8,
    )
    reference = [lod_error[0] * (value / dofs[0]) ** -0.5 for value in dofs]
    ax.plot(
        dofs,
        reference,
        color=COLORS["reference"],
        linestyle="--",
        linewidth=1.0,
        label=r"$O(N_H^{-1/2})$",
        zorder=1,
    )
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.grid(True, which="major")
    ax.grid(True, which="minor", alpha=0.12)
    ax.set_xlabel(r"Coarse-space degrees of freedom $N_H$")
    ax.set_ylabel(r"Absolute $k$-weighted energy error")
    ax.set_title("(a) Error versus coarse DOFs", loc="left")
    ax.legend(frameon=False, fontsize=6.8)
    ratio_first = p1_error[0] / lod_error[0]
    ratio_last = p1_error[-1] / lod_error[-1]
    ax.text(
        0.98,
        0.04,
        f"P1/LOD: {ratio_first:.1f}× → {ratio_last:.1f}×",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7.0,
    )

    ax = axes[1]
    ax.plot(
        finer_dofs,
        p1_orders,
        color=COLORS["fem"],
        marker="s",
        label="Coarse P1 FEM",
    )
    ax.plot(
        finer_dofs,
        lod_orders,
        color=COLORS["lod"],
        marker="o",
        label="Two-sided LOD",
    )
    ax.axhline(1.0, color=COLORS["reference"], linestyle="--", linewidth=1.0)
    ax.set_xscale("log", base=2)
    ax.set_xticks(finer_dofs[::2])
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.grid(True, which="major")
    ax.grid(True, which="minor", alpha=0.12)
    ax.set_xlabel(r"Finer endpoint DOFs $N_H$")
    ax.set_ylabel(r"Adjacent $H$-order $p_j$")
    ax.set_title("(b) Successive observed orders", loc="left")
    ax.legend(frameon=False, fontsize=6.8)

    fig.suptitle(
        r"Helmholtz global-NVB convergence "
        r"($k=32$, $h$-level 19, $\ell=4$)",
        y=1.01,
    )
    fig.tight_layout()

    output_dir.mkdir(parents=True, exist_ok=True)
    outputs: List[Path] = []
    for extension in formats:
        path = output_dir / f"helmholtz_H_convergence_server_energy.{extension}"
        fig.savefig(path)
        outputs.append(path)
    plt.close(fig)

    metrics = [
        ("p1_dof_fitted_slope_all_6", p1_slope),
        ("lod_dof_fitted_slope_all_6", lod_slope),
        ("p1_over_lod_at_H8", ratio_first),
        ("p1_over_lod_at_H13", ratio_last),
        (
            "max_petrov_residual",
            max(number(row, "petrov_residual") for row in rows),
        ),
        (
            "max_corrector_residual",
            max(number(row, "corrector_residual") for row in rows),
        ),
    ]
    metrics_path = output_dir / "helmholtz_H_convergence_server_metrics.csv"
    with metrics_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("metric", "value", "source"))
        for name, value in metrics:
            writer.writerow(
                (name, f"{value:.17g}", "results/helmholtz_H_convergence_server/all_results.csv")
            )
    outputs.append(metrics_path)
    return outputs, metrics


def generate(
    input_csv: Path, output_dir: Path, formats: Sequence[str]
) -> Tuple[List[Path], List[Metric]]:
    unknown = set(formats) - {"png", "pdf", "svg"}
    if unknown:
        raise ValueError(f"unsupported output formats: {sorted(unknown)}")
    return plot_convergence(load_and_validate(input_csv), output_dir, formats)


def parse_args() -> argparse.Namespace:
    repository_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=repository_root
        / "results"
        / "helmholtz_H_convergence_server"
        / "all_results.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=repository_root / "figures" / "paper",
    )
    parser.add_argument("--formats", default="png,pdf")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    formats = tuple(
        item.strip().lower() for item in args.formats.split(",") if item.strip()
    )
    outputs, metrics = generate(args.input, args.output_dir, formats)
    for path in outputs:
        print(path)
    for name, value in metrics:
        print(f"{name}={value:.8g}")


if __name__ == "__main__":
    main()

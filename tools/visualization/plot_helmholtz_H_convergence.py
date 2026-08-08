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


def optional_number(row: Mapping[str, str], field: str) -> float | None:
    try:
        value = float(row[field])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"invalid optional numeric field {field!r}") from error
    return value if math.isfinite(value) else None


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
    levels = [integer(row, "H_level") for row in rows]
    if levels != list(range(levels[0], levels[-1] + 1)):
        raise ValueError("completed H levels are not consecutive")
    expected = {
        "k": integer(rows[0], "k"),
        "h_level": integer(rows[0], "h_level"),
        "ell": integer(rows[0], "ell"),
        "solver": rows[0].get("solver"),
        "mode": rows[0].get("mode"),
    }
    for row in rows:
        if integer(row, "k") != expected["k"]:
            raise ValueError("mixed wave numbers in one convergence CSV")
        if integer(row, "h_level") != expected["h_level"]:
            raise ValueError("mixed fine levels in one convergence CSV")
        if integer(row, "ell") != expected["ell"]:
            raise ValueError("mixed oversampling levels in one convergence CSV")
        if row.get("solver") != expected["solver"] or row.get("mode") != expected["mode"]:
            raise ValueError("mixed solver or Petrov modes")
        for field in (
            "p1_energy_abs",
            "lod_energy_abs",
            "p1_l2_abs",
            "lod_l2_abs",
            "H_max",
        ):
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
        if expected["solver"] == "schur" and number(row, "schur_rcond") <= 1.0e-14:
            raise ValueError("Schur reciprocal condition estimate is too small")

    for field in ("p1_energy_abs", "lod_energy_abs", "p1_l2_abs", "lod_l2_abs"):
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
            ("p1_l2_abs", "p1_l2_rate"),
            ("lod_l2_abs", "lod_l2_rate"),
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
    finer_dofs = dofs[1:]
    k = integer(rows[0], "k")
    fine_level = integer(rows[0], "h_level")
    ell = integer(rows[0], "ell")
    first_level = integer(rows[0], "H_level")
    last_level = integer(rows[-1], "H_level")
    series = {
        "energy": {
            "p1": [number(row, "p1_energy_abs") for row in rows],
            "lod": [number(row, "lod_energy_abs") for row in rows],
            "p1_orders": [number(row, "p1_energy_rate") for row in rows[1:]],
            "lod_orders": [number(row, "lod_energy_rate") for row in rows[1:]],
            "fine_field": "fine_energy_abs",
            "ylabel": r"Absolute $k$-weighted energy error",
            "reference_power": -0.5,
            "reference_label": r"$O(N_H^{-1/2})$",
        },
        "l2": {
            "p1": [number(row, "p1_l2_abs") for row in rows],
            "lod": [number(row, "lod_l2_abs") for row in rows],
            "p1_orders": [number(row, "p1_l2_rate") for row in rows[1:]],
            "lod_orders": [number(row, "lod_l2_rate") for row in rows[1:]],
            "fine_field": "fine_l2_abs",
            "ylabel": r"Absolute $L^2$ error",
            "reference_power": -1.0,
            "reference_label": r"$O(N_H^{-1})$",
        },
    }
    for values in series.values():
        values["p1_slope"] = power_slope(dofs, values["p1"])
        values["lod_slope"] = power_slope(dofs, values["lod"])
        floors = [
            value
            for row in rows
            if (value := optional_number(row, values["fine_field"])) is not None
        ]
        values["fine_floor"] = floors[-1] if floors else None

    fig, axes = plt.subplots(2, 2, figsize=(7.25, 5.35))

    def plot_error_panel(ax: plt.Axes, values: dict, title: str) -> None:
        p1_values = values["p1"]
        lod_values = values["lod"]
        ax.plot(
            dofs,
            p1_values,
            color=COLORS["fem"],
            marker="s",
            label=rf"Coarse P1 FEM ($s={values['p1_slope']:.2f}$)",
            zorder=4,
        )
        ax.plot(
            dofs,
            lod_values,
            color=COLORS["lod"],
            marker="o",
            label=rf"Two-sided LOD ($s={values['lod_slope']:.2f}$)",
            zorder=4,
        )
        for field, color in (("p1", COLORS["fem"]), ("lod", COLORS["lod"])):
            ax.plot(
                dofs,
                fitted_values(dofs, values[field], values[f"{field}_slope"]),
                color=color,
                linestyle=":",
                linewidth=1.0,
                alpha=0.8,
            )
        reference = [
            lod_values[0] * (value / dofs[0]) ** values["reference_power"]
            for value in dofs
        ]
        ax.plot(
            dofs,
            reference,
            color=COLORS["reference"],
            linestyle="--",
            linewidth=1.0,
            label=values["reference_label"],
            zorder=1,
        )
        if values["fine_floor"] is not None:
            ax.axhline(
                values["fine_floor"],
                color=COLORS["fine"],
                linestyle="-.",
                linewidth=1.0,
                label=rf"Fine P1 floor ({values['fine_floor']:.3g})",
                zorder=2,
            )
        ax.set_xscale("log", base=2)
        ax.set_yscale("log")
        ax.grid(True, which="major")
        ax.grid(True, which="minor", alpha=0.12)
        ax.set_xlabel(r"Coarse-space DOFs $N_H$")
        ax.set_ylabel(values["ylabel"])
        ax.set_title(title, loc="left")
        ax.legend(frameon=False, fontsize=6.2)
        ax.text(
            0.98,
            0.04,
            f"P1/LOD: {p1_values[0] / lod_values[0]:.1f}×"
            f" → {p1_values[-1] / lod_values[-1]:.1f}×",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=6.5,
        )

    def plot_order_panel(ax: plt.Axes, values: dict, title: str) -> None:
        ax.plot(
            finer_dofs,
            values["p1_orders"],
            color=COLORS["fem"],
            marker="s",
            label="Coarse P1 FEM",
        )
        ax.plot(
            finer_dofs,
            values["lod_orders"],
            color=COLORS["lod"],
            marker="o",
            label="Two-sided LOD",
        )
        ax.axhline(
            -2.0 * values["reference_power"],
            color=COLORS["reference"],
            linestyle="--",
            linewidth=1.0,
        )
        ax.set_xscale("log", base=2)
        ax.set_xticks(finer_dofs[::2])
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.grid(True, which="major")
        ax.grid(True, which="minor", alpha=0.12)
        ax.set_xlabel(r"Finer endpoint DOFs $N_H$")
        ax.set_ylabel(r"Adjacent $H$-order $p_j$")
        ax.set_title(title, loc="left")
        ax.legend(frameon=False, fontsize=6.2)

    plot_error_panel(axes[0, 0], series["energy"], "(a) Weighted energy error")
    plot_error_panel(axes[0, 1], series["l2"], r"(b) $L^2$ error")
    plot_order_panel(axes[1, 0], series["energy"], "(c) Energy-error orders")
    plot_order_panel(axes[1, 1], series["l2"], r"(d) $L^2$-error orders")

    fig.suptitle(
        r"Helmholtz global-NVB convergence "
        rf"($k={k}$, $h$-level {fine_level}, $\ell={ell}$)",
        y=1.01,
    )
    fig.tight_layout()

    output_dir.mkdir(parents=True, exist_ok=True)
    outputs: List[Path] = []
    figure_stem = f"helmholtz_H_convergence_k{k}_server_energy_l2"
    for extension in formats:
        path = output_dir / f"{figure_stem}.{extension}"
        fig.savefig(path)
        outputs.append(path)
    plt.close(fig)

    metrics = [
        (
            f"p1_energy_dof_fitted_slope_all_{len(rows)}",
            series["energy"]["p1_slope"],
        ),
        (
            f"lod_energy_dof_fitted_slope_all_{len(rows)}",
            series["energy"]["lod_slope"],
        ),
        (f"p1_energy_over_lod_at_H{first_level}", series["energy"]["p1"][0] / series["energy"]["lod"][0]),
        (f"p1_energy_over_lod_at_H{last_level}", series["energy"]["p1"][-1] / series["energy"]["lod"][-1]),
        (
            f"p1_l2_dof_fitted_slope_all_{len(rows)}",
            series["l2"]["p1_slope"],
        ),
        (
            f"lod_l2_dof_fitted_slope_all_{len(rows)}",
            series["l2"]["lod_slope"],
        ),
        (f"p1_l2_over_lod_at_H{first_level}", series["l2"]["p1"][0] / series["l2"]["lod"][0]),
        (f"p1_l2_over_lod_at_H{last_level}", series["l2"]["p1"][-1] / series["l2"]["lod"][-1]),
        (
            "max_petrov_residual",
            max(number(row, "petrov_residual") for row in rows),
        ),
        (
            "max_corrector_residual",
            max(number(row, "corrector_residual") for row in rows),
        ),
    ]
    for name, values in series.items():
        if values["fine_floor"] is not None:
            metrics.extend(
                [
                    (f"fine_exact_{name}_floor", values["fine_floor"]),
                    (
                        f"lod_{name}_over_fine_floor_at_last",
                        values["lod"][-1] / values["fine_floor"],
                    ),
                ]
            )
    metrics_path = output_dir / f"helmholtz_H_convergence_k{k}_server_metrics.csv"
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

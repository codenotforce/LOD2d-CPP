#!/usr/bin/env python3
"""Generate paper figures from checked-in Helmholtz CSV results.

This script is intentionally a pure post-processing step. It does not invoke
the C++ solver, recompute reference solutions, or retain any LOD model data.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

from paper_style import COLORS, apply_paper_style


Row = Dict[str, str]
Metric = Tuple[str, str, str, float, str]


def read_csv(path: Path) -> List[Row]:
    if not path.is_file():
        raise FileNotFoundError(f"missing input CSV: {path}")
    with path.open("r", encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"input CSV contains no rows: {path}")
    return rows


def number(row: Mapping[str, str], field: str) -> float:
    try:
        value = float(row[field])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"invalid numeric field {field!r} in row {row}") from error
    if not math.isfinite(value):
        raise ValueError(f"non-finite field {field!r} in row {row}")
    return value


def integer(row: Mapping[str, str], field: str) -> int:
    value = number(row, field)
    rounded = round(value)
    if not math.isclose(value, rounded, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValueError(f"field {field!r} is not integral: {value}")
    return int(rounded)


def merge_by_wavenumber(row_groups: Iterable[Iterable[Row]]) -> List[Row]:
    merged: MutableMapping[int, Row] = {}
    comparison_fields = (
        "fem_exact_energy_rel",
        "lod_exact_energy_rel",
        "fine_exact_energy_rel",
    )
    for rows in row_groups:
        for row in rows:
            if row.get("source") != "manufactured":
                continue
            k = integer(row, "k")
            if k in merged:
                for field in comparison_fields:
                    if not math.isclose(
                        number(row, field),
                        number(merged[k], field),
                        rel_tol=5.0e-11,
                        abs_tol=5.0e-13,
                    ):
                        raise ValueError(
                            f"duplicate k={k} disagrees in {field}: "
                            f"{row[field]} versus {merged[k][field]}"
                        )
                continue
            merged[k] = row
    return [merged[k] for k in sorted(merged)]


def load_pollution_data(results_root: Path) -> Dict[str, List[Row]]:
    main_local = read_csv(
        results_root / "helmholtz_pollution_gap6" / "scan_kH1p0.csv"
    )
    main_server = read_csv(
        results_root
        / "helmholtz_pollution_server"
        / "summary_main_k64_gap6_t32.csv"
    )
    strict_local = read_csv(
        results_root / "helmholtz_pollution" / "scan_kH1p0.csv"
    )
    strict_server_32 = read_csv(
        results_root
        / "helmholtz_pollution_server"
        / "summary_strict_k32_gap8_t32.csv"
    )
    strict_server_64 = read_csv(
        results_root
        / "helmholtz_pollution_server"
        / "summary_strict_k64_gap8_t32.csv"
    )

    datasets = {
        "main": merge_by_wavenumber((main_local, main_server)),
        "strict": merge_by_wavenumber(
            (strict_local, strict_server_32, strict_server_64)
        ),
    }
    expected_k = [4, 8, 16, 32, 64]
    expected_kh = {"main": 0.125, "strict": 0.0625}
    for name, rows in datasets.items():
        actual_k = [integer(row, "k") for row in rows]
        if actual_k != expected_k:
            raise ValueError(f"{name} pollution scan has k={actual_k}, expected {expected_k}")
        for row in rows:
            if not math.isclose(number(row, "kH"), 1.0, abs_tol=1.0e-12):
                raise ValueError(f"{name} scan is not at fixed kH=1")
            if not math.isclose(
                number(row, "kh"), expected_kh[name], abs_tol=1.0e-12
            ):
                raise ValueError(
                    f"{name} scan has kh={row['kh']}, expected {expected_kh[name]}"
                )

    # The coarse FEM solution must not depend on the auxiliary fine space.
    for main, strict in zip(datasets["main"], datasets["strict"]):
        if not math.isclose(
            number(main, "fem_exact_energy_rel"),
            number(strict, "fem_exact_energy_rel"),
            rel_tol=5.0e-9,
            abs_tol=1.0e-10,
        ):
            raise ValueError(
                f"coarse FEM exact error changes with fine gap at k={main['k']}"
            )
    return datasets


def percent_change(first: float, last: float) -> float:
    return 100.0 * (last / first - 1.0)


def power_rate(x: Sequence[float], y: Sequence[float]) -> float:
    if len(x) != len(y) or len(x) < 2:
        raise ValueError("power-rate fit requires equal sequences with at least two points")
    log_x = [math.log(value) for value in x]
    log_y = [math.log(value) for value in y]
    mean_x = sum(log_x) / len(log_x)
    mean_y = sum(log_y) / len(log_y)
    denominator = sum((value - mean_x) ** 2 for value in log_x)
    if denominator == 0.0:
        raise ValueError("power-rate fit has zero x variance")
    return sum(
        (x_value - mean_x) * (y_value - mean_y)
        for x_value, y_value in zip(log_x, log_y)
    ) / denominator


def configure_log_wavenumber_axis(ax: plt.Axes, ks: Sequence[int]) -> None:
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xticks(ks)
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.grid(True, which="major")
    ax.grid(True, which="minor", alpha=0.12)
    ax.set_xlabel(r"Wave number $k$ (fixed $kH=1$)")
    ax.set_ylabel(r"Relative weighted energy error")
    ax.set_ylim(8.0e-3, 1.05)


def plot_pollution(
    datasets: Mapping[str, Sequence[Row]], output_dir: Path, formats: Sequence[str]
) -> Tuple[List[Path], List[Metric]]:
    apply_paper_style()
    fig, axes = plt.subplots(1, 2, figsize=(7.25, 3.15), sharey=True)
    descriptions = {
        "main": r"(a) Primary fine space: $kh=1/8$",
        "strict": r"(b) Strict reference: $kh=1/16$",
    }
    metrics: List[Metric] = []

    for ax, name in zip(axes, ("main", "strict")):
        rows = list(datasets[name])
        ks = [integer(row, "k") for row in rows]
        fem = [number(row, "fem_exact_energy_rel") for row in rows]
        lod = [number(row, "lod_exact_energy_rel") for row in rows]
        fine = [number(row, "fine_exact_energy_rel") for row in rows]

        ax.plot(
            ks,
            fem,
            color=COLORS["fem"],
            marker="s",
            label="Coarse P1 FEM",
            zorder=3,
        )
        ax.plot(
            ks,
            lod,
            color=COLORS["lod"],
            marker="o",
            label="Two-sided LOD",
            zorder=4,
        )
        ax.plot(
            ks,
            fine,
            color=COLORS["fine"],
            marker="^",
            linestyle="--",
            label="Fine P1 floor",
            zorder=2,
        )
        configure_log_wavenumber_axis(ax, ks)
        ax.set_title(descriptions[name], loc="left")

        fem_32_64 = percent_change(fem[-2], fem[-1])
        lod_32_64 = percent_change(lod[-2], lod[-1])
        ratio_64 = fem[-1] / lod[-1]
        sign = "+" if lod_32_64 >= 0.0 else ""
        ax.text(
            0.035,
            0.055,
            rf"$k=32\!\to\!64$: FEM {fem_32_64:+.1f}%, "
            rf"LOD {sign}{lod_32_64:.1f}%"
            "\n"
            rf"at $k=64$: FEM/LOD = {ratio_64:.1f}",
            transform=ax.transAxes,
            va="bottom",
            ha="left",
            fontsize=7.4,
            bbox={
                "facecolor": "white",
                "edgecolor": "none",
                "alpha": 0.86,
                "pad": 1.2,
            },
        )

        source = (
            "results/helmholtz_pollution_gap6 + pollution_server"
            if name == "main"
            else "results/helmholtz_pollution + pollution_server"
        )
        metrics.extend(
            [
                ("pollution", name, "fem_change_k32_to_k64_percent", fem_32_64, source),
                ("pollution", name, "lod_change_k32_to_k64_percent", lod_32_64, source),
                ("pollution", name, "fem_over_lod_at_k64", ratio_64, source),
            ]
        )

    axes[0].legend(loc="upper left", frameon=False)
    fig.suptitle(
        r"Manufactured Helmholtz pollution scan "
        r"($\ell=\lceil\log_2 k\rceil$)",
        y=1.01,
    )
    fig.tight_layout()

    outputs: List[Path] = []
    for extension in formats:
        path = output_dir / f"helmholtz_pollution_weighted_energy.{extension}"
        fig.savefig(path)
        outputs.append(path)
    plt.close(fig)
    return outputs, metrics


def load_manufactured_data(results_root: Path) -> Tuple[List[Row], Dict[int, List[Row]]]:
    rows = read_csv(results_root / "helmholtz_manufactured" / "validation.csv")
    fem = sorted(
        (row for row in rows if row.get("study") == "fem"),
        key=lambda row: integer(row, "fine_nodes"),
    )
    lod: Dict[int, List[Row]] = {}
    for ell in (2, 3):
        lod[ell] = sorted(
            (
                row
                for row in rows
                if row.get("study") == "lod" and integer(row, "ell") == ell
            ),
            key=lambda row: integer(row, "coarse_nodes"),
        )
    if len(fem) != 3 or any(len(lod[ell]) != 4 for ell in lod):
        raise ValueError(
            "manufactured validation must contain 3 FEM rows and 4 LOD rows per ell"
        )
    return fem, lod


def add_reference_power(
    ax: plt.Axes,
    x: Sequence[float],
    anchor_error: float,
    power: float,
    label: str,
) -> None:
    x_min = min(x)
    reference = [anchor_error * (value / x_min) ** power for value in x]
    ax.plot(
        x,
        reference,
        color=COLORS["reference"],
        linestyle=":",
        linewidth=1.0,
        label=label,
        zorder=1,
    )


def configure_dof_axis(ax: plt.Axes, label: str) -> None:
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.grid(True, which="major")
    ax.grid(True, which="minor", alpha=0.12)
    ax.set_xlabel(label)
    ax.set_ylabel("Exact error")


def plot_manufactured_convergence(
    fem_rows: Sequence[Row],
    lod_rows: Mapping[int, Sequence[Row]],
    output_dir: Path,
    formats: Sequence[str],
) -> Tuple[List[Path], List[Metric]]:
    apply_paper_style()
    fig, axes = plt.subplots(1, 2, figsize=(7.25, 3.2))
    metrics: List[Metric] = []

    fine_dofs = [integer(row, "fine_nodes") for row in fem_rows]
    fine_energy = [number(row, "fem_exact_energy") for row in fem_rows]
    fine_l2 = [number(row, "fem_exact_l2") for row in fem_rows]
    fine_energy_rate = power_rate(fine_dofs, fine_energy)
    fine_l2_rate = power_rate(fine_dofs, fine_l2)

    axes[0].plot(
        fine_dofs,
        fine_energy,
        color=COLORS["energy"],
        marker="o",
        label=rf"Energy ($p={fine_energy_rate:.2f}$)",
        zorder=3,
    )
    axes[0].plot(
        fine_dofs,
        fine_l2,
        color=COLORS["l2"],
        marker="s",
        linestyle="--",
        label=rf"$L^2$ ($p={fine_l2_rate:.2f}$)",
        zorder=3,
    )
    add_reference_power(
        axes[0], fine_dofs, fine_energy[0], -0.5, r"$O(N_h^{-1/2})$"
    )
    add_reference_power(
        axes[0], fine_dofs, fine_l2[0], -1.0, r"$O(N_h^{-1})$"
    )
    configure_dof_axis(axes[0], r"Fine-space degrees of freedom $N_h$")
    axes[0].set_title("(a) Global-NVB fine P1 FEM", loc="left")
    axes[0].legend(
        frameon=False,
        ncol=2,
        fontsize=6.3,
        handlelength=1.5,
        columnspacing=0.8,
        handletextpad=0.4,
    )

    metrics.extend(
        [
            (
                "manufactured_convergence",
                "fine_fem",
                "energy_dof_fitted_slope",
                fine_energy_rate,
                "results/helmholtz_manufactured/validation.csv",
            ),
            (
                "manufactured_convergence",
                "fine_fem",
                "l2_dof_fitted_slope",
                fine_l2_rate,
                "results/helmholtz_manufactured/validation.csv",
            ),
        ]
    )

    markers = {2: "o", 3: "^"}
    line_styles = {2: "--", 3: "-"}
    all_dofs: List[float] = []
    ell3_energy: List[float] = []
    ell3_l2: List[float] = []
    for ell in (2, 3):
        rows = list(lod_rows[ell])
        coarse_dofs = [integer(row, "coarse_nodes") for row in rows]
        energy = [number(row, "lod_exact_energy") for row in rows]
        l2_error = [number(row, "lod_exact_l2") for row in rows]
        energy_rate = power_rate(coarse_dofs, energy)
        l2_rate = power_rate(coarse_dofs, l2_error)
        all_dofs = coarse_dofs
        if ell == 3:
            ell3_energy = energy
            ell3_l2 = l2_error

        axes[1].plot(
            coarse_dofs,
            energy,
            color=COLORS["energy"],
            marker=markers[ell],
            linestyle=line_styles[ell],
            label=rf"Energy, $\ell={ell}$ ($p={energy_rate:.2f}$)",
            zorder=3,
        )
        axes[1].plot(
            coarse_dofs,
            l2_error,
            color=COLORS["l2"],
            marker=markers[ell],
            linestyle=line_styles[ell],
            label=rf"$L^2$, $\ell={ell}$ ($p={l2_rate:.2f}$)",
            zorder=3,
        )
        metrics.extend(
            [
                (
                    "manufactured_convergence",
                    f"lod_ell_{ell}",
                    "energy_dof_fitted_slope",
                    energy_rate,
                    "results/helmholtz_manufactured/validation.csv",
                ),
                (
                    "manufactured_convergence",
                    f"lod_ell_{ell}",
                    "l2_dof_fitted_slope",
                    l2_rate,
                    "results/helmholtz_manufactured/validation.csv",
                ),
            ]
        )

    add_reference_power(
        axes[1], all_dofs, ell3_energy[0], -0.5, r"$O(N_H^{-1/2})$"
    )
    add_reference_power(
        axes[1], all_dofs, ell3_l2[0], -1.0, r"$O(N_H^{-1})$"
    )
    configure_dof_axis(axes[1], r"Coarse LOD degrees of freedom $N_H$")
    axes[1].set_title(r"(b) Two-sided LOD, fixed fine level $h=10$", loc="left")
    axes[1].legend(
        frameon=False,
        ncol=2,
        fontsize=5.7,
        handlelength=1.35,
        columnspacing=0.65,
        handletextpad=0.35,
        labelspacing=0.3,
    )

    fig.suptitle(r"Manufactured global-NVB convergence ($k=4$)", y=1.01)
    fig.tight_layout()

    outputs: List[Path] = []
    for extension in formats:
        path = output_dir / f"manufactured_global_nvb_convergence.{extension}"
        fig.savefig(path)
        outputs.append(path)
    plt.close(fig)
    return outputs, metrics


def write_metrics(path: Path, metrics: Sequence[Metric]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("figure", "series", "metric", "value", "source"))
        for figure, series, metric, value, source in metrics:
            writer.writerow((figure, series, metric, f"{value:.17g}", source))


def generate_figures(
    results_root: Path, output_dir: Path, formats: Sequence[str]
) -> Tuple[List[Path], List[Metric]]:
    allowed_formats = {"png", "pdf", "svg"}
    unknown = set(formats) - allowed_formats
    if unknown:
        raise ValueError(f"unsupported output formats: {sorted(unknown)}")
    output_dir.mkdir(parents=True, exist_ok=True)

    pollution = load_pollution_data(results_root)
    fem, lod = load_manufactured_data(results_root)
    pollution_outputs, pollution_metrics = plot_pollution(
        pollution, output_dir, formats
    )
    convergence_outputs, convergence_metrics = plot_manufactured_convergence(
        fem, lod, output_dir, formats
    )
    metrics = pollution_metrics + convergence_metrics
    metrics_path = output_dir / "helmholtz_plot_metrics.csv"
    write_metrics(metrics_path, metrics)
    return pollution_outputs + convergence_outputs + [metrics_path], metrics


def parse_args() -> argparse.Namespace:
    repository_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(
        description="Plot Helmholtz pollution and manufactured convergence results."
    )
    parser.add_argument(
        "--results-root",
        type=Path,
        default=repository_root / "results",
        help="directory containing checked-in result folders",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=repository_root / "figures" / "paper",
        help="destination for paper figures",
    )
    parser.add_argument(
        "--formats",
        default="png,pdf",
        help="comma-separated subset of png,pdf,svg",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    formats = tuple(
        item.strip().lower() for item in args.formats.split(",") if item.strip()
    )
    outputs, metrics = generate_figures(args.results_root, args.output_dir, formats)
    for path in outputs:
        print(path)
    print("metrics:")
    for figure, series, metric, value, _ in metrics:
        print(f"  {figure}.{series}.{metric}={value:.6g}")


if __name__ == "__main__":
    main()

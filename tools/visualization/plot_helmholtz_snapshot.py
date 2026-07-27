#!/usr/bin/env python3
"""Validate and plot an explicitly exported Helmholtz VTU snapshot."""

from __future__ import annotations

import argparse
import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, Mapping, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

from paper_style import apply_paper_style


ArrayMap = Dict[str, np.ndarray]


def _array(element: ET.Element, dtype: type = float) -> np.ndarray:
    return np.fromstring(element.text or "", sep=" ", dtype=dtype)


def read_vtu(path: Path) -> Tuple[np.ndarray, np.ndarray, ArrayMap, ArrayMap]:
    if not path.is_file():
        raise FileNotFoundError(f"missing VTU file: {path}")
    piece = ET.parse(path).getroot().find("./UnstructuredGrid/Piece")
    if piece is None:
        raise ValueError(f"VTU has no UnstructuredGrid/Piece: {path}")
    point_count = int(piece.attrib["NumberOfPoints"])
    cell_count = int(piece.attrib["NumberOfCells"])

    point_array = piece.find("./Points/DataArray")
    if point_array is None:
        raise ValueError("VTU has no point coordinate array")
    points = _array(point_array).reshape((-1, 3))[:, :2]

    cell_arrays = {
        item.attrib.get("Name", ""): item for item in piece.findall("./Cells/DataArray")
    }
    connectivity = _array(cell_arrays["connectivity"], int)
    offsets = _array(cell_arrays["offsets"], int)
    cell_types = _array(cell_arrays["types"], int)
    if connectivity.size != 3 * cell_count or not np.array_equal(
        offsets, 3 * np.arange(1, cell_count + 1)
    ):
        raise ValueError("VTU connectivity/offsets do not describe P1 triangles")
    if not np.all(cell_types == 5):
        raise ValueError("VTU contains a non-triangle cell")
    triangles = connectivity.reshape((-1, 3))

    def fields(section: str, expected: int) -> ArrayMap:
        result: ArrayMap = {}
        for item in piece.findall(f"./{section}/DataArray"):
            name = item.attrib.get("Name", "")
            values = _array(item, int if item.attrib.get("type", "").startswith(("Int", "UInt")) else float)
            if not name or values.size != expected:
                raise ValueError(f"invalid {section} field {name!r}")
            if name in result:
                raise ValueError(f"duplicate {section} field {name!r}")
            result[name] = values
        return result

    if points.shape != (point_count, 2):
        raise ValueError("VTU point count disagrees with coordinate data")
    if triangles.shape != (cell_count, 3):
        raise ValueError("VTU cell count disagrees with connectivity data")
    if triangles.size and (triangles.min() < 0 or triangles.max() >= point_count):
        raise ValueError("VTU triangle has an invalid point index")
    if not np.isfinite(points).all():
        raise ValueError("VTU coordinates contain NaN/Inf")
    point_data = fields("PointData", point_count)
    cell_data = fields("CellData", cell_count)
    for values in (*point_data.values(), *cell_data.values()):
        if not np.isfinite(values).all():
            raise ValueError("VTU field contains NaN/Inf")
    return points, triangles, point_data, cell_data


def complex_field(fields: Mapping[str, np.ndarray], name: str) -> np.ndarray:
    try:
        real = fields[f"{name}_real"]
        imag = fields[f"{name}_imag"]
    except KeyError as error:
        raise ValueError(f"missing complex VTU field {name!r}") from error
    return real + 1j * imag


def load_case(case_dir: Path) -> Tuple[dict, tuple, tuple]:
    manifest_path = case_dir / "run.json"
    with manifest_path.open("r", encoding="utf-8") as stream:
        manifest = json.load(stream)
    if manifest.get("schema") != "lod2d.helmholtz.visualization.v1":
        raise ValueError("unsupported visualization manifest schema")
    files = manifest["files"]
    coarse = read_vtu(case_dir / files["coarse_mesh"])
    fine = read_vtu(case_dir / files["fine_solution"])
    expected = manifest["mesh"]
    if coarse[0].shape[0] != expected["coarse_nodes"]:
        raise ValueError("manifest coarse node count disagrees with VTU")
    if coarse[1].shape[0] != expected["coarse_elements"]:
        raise ValueError("manifest coarse element count disagrees with VTU")
    if fine[0].shape[0] != expected["fine_nodes"]:
        raise ValueError("manifest fine node count disagrees with VTU")
    if fine[1].shape[0] != expected["fine_elements"]:
        raise ValueError("manifest fine element count disagrees with VTU")
    return manifest, coarse, fine


def _save(fig: plt.Figure, output_dir: Path, stem: str, formats: Sequence[str]) -> list[Path]:
    paths = []
    for extension in formats:
        path = output_dir / f"{stem}.{extension}"
        fig.savefig(path, dpi=400 if extension == "png" else None)
        paths.append(path)
    plt.close(fig)
    return paths


def plot_meshes(
    coarse: tuple, fine: tuple, output_dir: Path, formats: Sequence[str]
) -> list[Path]:
    apply_paper_style()
    fig, axes = plt.subplots(1, 2, figsize=(7.25, 3.35))
    for ax, data, title, width in (
        (axes[0], coarse, "(a) Coarse LOD mesh", 0.55),
        (axes[1], fine, "(b) Fine reference mesh", 0.18),
    ):
        points, triangles = data[:2]
        ax.triplot(points[:, 0], points[:, 1], triangles, color="#23313f", linewidth=width)
        ax.set_aspect("equal")
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.0)
        ax.set_title(title, loc="left")
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(r"$y$")
        ax.text(
            0.02,
            0.02,
            f"{points.shape[0]} nodes, {triangles.shape[0]} cells",
            transform=ax.transAxes,
            fontsize=7.5,
            va="bottom",
        )
    fig.tight_layout()
    return _save(fig, output_dir, "helmholtz_static_meshes", formats)


def _field_panel(
    ax: plt.Axes,
    triangulation: mtri.Triangulation,
    values: np.ndarray,
    title: str,
    cmap: str,
    vmin: float | None = None,
    vmax: float | None = None,
) -> None:
    image = ax.tripcolor(
        triangulation, values, shading="gouraud", cmap=cmap, vmin=vmin, vmax=vmax
    )
    ax.set_aspect("equal")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_title(title, loc="left")
    ax.set_xticks((0.0, 0.5, 1.0))
    ax.set_yticks((0.0, 0.5, 1.0))
    ax.figure.colorbar(image, ax=ax, fraction=0.046, pad=0.03)


def plot_decomposition(
    manifest: dict, fine: tuple, output_dir: Path, formats: Sequence[str]
) -> list[Path]:
    apply_paper_style()
    points, triangles, fields, _ = fine
    tri = mtri.Triangulation(points[:, 0], points[:, 1], triangles)
    exact = complex_field(fields, "u_exact")
    reference = complex_field(fields, "u_reference")
    lod = complex_field(fields, "u_lod")
    coarse = complex_field(fields, "u_coarse")
    fine_scale = complex_field(fields, "u_fine_scale")
    error = complex_field(fields, "error_lod")
    if not np.allclose(lod, coarse + fine_scale, rtol=2e-14, atol=2e-14):
        raise ValueError("exported LOD decomposition is algebraically inconsistent")
    if not np.allclose(error, lod - reference, rtol=2e-14, atol=2e-14):
        raise ValueError("exported LOD-reference error is inconsistent")

    real_limit = max(
        np.max(np.abs(exact.real)),
        np.max(np.abs(reference.real)),
        np.max(np.abs(lod.real)),
        np.max(np.abs(coarse.real)),
    )
    correction_limit = np.max(np.abs(fine_scale.real))
    fig, axes = plt.subplots(2, 3, figsize=(7.25, 5.0))
    panels = (
        (exact.real, "(a) Exact, real part", "RdBu_r", -real_limit, real_limit),
        (reference.real, "(b) Fine P1, real part", "RdBu_r", -real_limit, real_limit),
        (lod.real, "(c) LOD, real part", "RdBu_r", -real_limit, real_limit),
        (coarse.real, "(d) Coarse prolongation", "RdBu_r", -real_limit, real_limit),
        (
            fine_scale.real,
            "(e) LOD fine-scale correction",
            "RdBu_r",
            -correction_limit,
            correction_limit,
        ),
        (np.abs(error), r"(f) $|u_{\rm LOD}-u_h|$", "magma", 0.0, None),
    )
    for ax, (values, title, cmap, vmin, vmax) in zip(axes.flat, panels):
        _field_panel(ax, tri, values, title, cmap, vmin, vmax)
    parameters = manifest["parameters"]
    fig.suptitle(
        rf"Manufactured Helmholtz snapshot "
        rf"($k={parameters['k']:g}$, $\ell={parameters['ell']}$)",
        y=1.0,
    )
    fig.tight_layout()
    return _save(fig, output_dir, "helmholtz_solution_decomposition", formats)


def plot_complex_and_section(
    fine: tuple, output_dir: Path, formats: Sequence[str]
) -> list[Path]:
    apply_paper_style()
    points, triangles, fields, _ = fine
    tri = mtri.Triangulation(points[:, 0], points[:, 1], triangles)
    lod = complex_field(fields, "u_lod")
    amplitude = np.abs(lod)
    phase = np.ma.masked_where(amplitude < 1.0e-8 * amplitude.max(), np.angle(lod))
    limit = max(np.max(np.abs(lod.real)), np.max(np.abs(lod.imag)))

    fig, axes = plt.subplots(2, 2, figsize=(6.5, 5.2))
    _field_panel(axes[0, 0], tri, lod.real, "(a) Real part", "RdBu_r", -limit, limit)
    _field_panel(axes[0, 1], tri, lod.imag, "(b) Imaginary part", "RdBu_r", -limit, limit)
    _field_panel(axes[1, 0], tri, amplitude, "(c) Magnitude", "viridis", 0.0, None)
    _field_panel(axes[1, 1], tri, phase, "(d) Phase (zero-amplitude masked)", "twilight", -math.pi, math.pi)
    fig.tight_layout()
    outputs = _save(fig, output_dir, "helmholtz_lod_complex_field", formats)

    y_values = np.unique(points[:, 1])
    y_line = y_values[np.argmin(np.abs(y_values - 0.5))]
    mask = np.isclose(points[:, 1], y_line, rtol=0.0, atol=2.0e-13)
    order = np.argsort(points[mask, 0])
    x = points[mask, 0][order]
    fig, ax = plt.subplots(figsize=(7.25, 2.85))
    styles = (
        ("u_exact", "Exact", "#111111", "-"),
        ("u_reference", "Fine P1", "#E69F00", "--"),
        ("u_lod", "LOD", "#0072B2", "-"),
        ("u_coarse", "Coarse prolongation", "#009E73", ":"),
    )
    for name, label, color, style in styles:
        values = complex_field(fields, name)[mask][order].real
        ax.plot(x, values, color=color, linestyle=style, label=label)
    ax.set_xlabel(r"$x$ at $y=0.5$")
    ax.set_ylabel("Real part")
    ax.grid(True)
    ax.legend(frameon=False, ncol=4)
    ax.set_title("Fixed cross-section resolves the oscillatory correction", loc="left")
    fig.tight_layout()
    outputs.extend(_save(fig, output_dir, "helmholtz_centerline", formats))
    return outputs


def generate_snapshot_figures(
    case_dir: Path, output_dir: Path, formats: Sequence[str]
) -> list[Path]:
    unknown = set(formats) - {"png", "pdf", "svg"}
    if unknown:
        raise ValueError(f"unsupported output formats: {sorted(unknown)}")
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest, coarse, fine = load_case(case_dir)
    outputs = plot_meshes(coarse, fine, output_dir, formats)
    outputs.extend(plot_decomposition(manifest, fine, output_dir, formats))
    outputs.extend(plot_complex_and_section(fine, output_dir, formats))
    return outputs


def parse_args() -> argparse.Namespace:
    repository_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case-dir",
        type=Path,
        default=repository_root
        / "results"
        / "visualization"
        / "helmholtz_manufactured_k4_H4_h8_ell3",
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
    formats = tuple(item.strip().lower() for item in args.formats.split(",") if item.strip())
    for path in generate_snapshot_figures(args.case_dir, args.output_dir, formats):
        print(path)


if __name__ == "__main__":
    main()

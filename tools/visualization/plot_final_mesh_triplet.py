#!/usr/bin/env python3
"""Render final coarse/reference/candidate meshes with a corner zoom."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from paper_style import apply_paper_style
from plot_reference_epoch_meshes import mesh_edges, read_vtu


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def draw(ax, points, triangles, limits, linewidth: float, title: str) -> None:
    vertices = points[triangles]
    if limits is not None:
        xmin, xmax, ymin, ymax = limits
        selected = (
            (vertices[:, :, 0].max(axis=1) >= xmin)
            & (vertices[:, :, 0].min(axis=1) <= xmax)
            & (vertices[:, :, 1].max(axis=1) >= ymin)
            & (vertices[:, :, 1].min(axis=1) <= ymax)
        )
        triangles = triangles[selected]
    ax.add_collection(LineCollection(
        mesh_edges(points, triangles), colors="#26364a",
        linewidths=linewidth, rasterized=True,
    ))
    if limits is None:
        lower, upper = points.min(axis=0), points.max(axis=0)
        ax.set_xlim(lower[0], upper[0])
        ax.set_ylim(lower[1], upper[1])
    else:
        ax.set_xlim(limits[0], limits[1])
        ax.set_ylim(limits[2], limits[3])
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title, fontsize=8.2)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coarse", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--title", default="Final PALOD meshes")
    parser.add_argument(
        "--corner-zoom", type=float, nargs=4,
        default=(-0.3, 0.3, -0.3, 0.3),
        metavar=("XMIN", "XMAX", "YMIN", "YMAX"),
    )
    arguments = parser.parse_args()

    meshes = []
    for label, path in (
        ("coarse", arguments.coarse),
        ("reference", arguments.reference),
        ("candidate", arguments.candidate),
    ):
        points, triangles = read_vtu(path)
        meshes.append((label, path, points, triangles))
    reference_equals_candidate = digest(arguments.reference) == digest(
        arguments.candidate
    )

    apply_paper_style()
    figure, axes = plt.subplots(
        2, 3, figsize=(10.2, 6.7), constrained_layout=True,
    )
    for column, (label, _, points, triangles) in enumerate(meshes):
        suffix = " (same as reference)" if (
            label == "candidate" and reference_equals_candidate
        ) else ""
        draw(
            axes[0, column], points, triangles, None,
            0.22 if label == "coarse" else 0.045,
            f"{label.capitalize()}{suffix}\n"
            f"{len(points):,} vertices, {len(triangles):,} cells",
        )
        draw(
            axes[1, column], points, triangles,
            tuple(arguments.corner_zoom),
            0.24 if label == "coarse" else 0.07,
            "Detail zoom",
        )
    title = arguments.title
    if reference_equals_candidate:
        title += " (reference and candidate coincide)"
    figure.suptitle(title, fontsize=10.5)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(arguments.output.with_suffix(".png"), dpi=300)
    figure.savefig(arguments.output.with_suffix(".pdf"), dpi=300)
    plt.close(figure)


if __name__ == "__main__":
    main()

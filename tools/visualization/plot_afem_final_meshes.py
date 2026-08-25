#!/usr/bin/env python3
"""Render the final E1 and E2 AFEM meshes with problem-specific zooms."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from paper_style import apply_paper_style
from plot_reference_epoch_meshes import mesh_edges, read_vtu


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
    ax.set_title(title, fontsize=8.6)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--e1", type=Path, required=True)
    parser.add_argument("--e2", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--e1-zoom", type=float, nargs=4,
        default=(0.5, 0.9, 0.3, 0.7),
    )
    parser.add_argument(
        "--e2-zoom", type=float, nargs=4,
        default=(-0.35, 0.35, -0.35, 0.35),
    )
    arguments = parser.parse_args()
    cases = []
    for label, path, zoom in (
        ("E1", arguments.e1, tuple(arguments.e1_zoom)),
        ("E2", arguments.e2, tuple(arguments.e2_zoom)),
    ):
        points, triangles = read_vtu(path)
        cases.append((label, points, triangles, zoom))

    apply_paper_style()
    figure, axes = plt.subplots(
        2, 2, figsize=(8.2, 7.15), constrained_layout=True,
    )
    for column, (label, points, triangles, zoom) in enumerate(cases):
        draw(
            axes[0, column], points, triangles, None, 0.095,
            f"{label} AFEM final mesh\n"
            f"{len(points):,} vertices, {len(triangles):,} cells",
        )
        draw(
            axes[1, column], points, triangles, zoom, 0.13,
            f"{label} refinement-region zoom",
        )
    figure.suptitle("Final AFEM meshes", fontsize=10.5)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(arguments.output.with_suffix(".png"), dpi=300)
    figure.savefig(arguments.output.with_suffix(".pdf"), dpi=300)
    plt.close(figure)


if __name__ == "__main__":
    main()

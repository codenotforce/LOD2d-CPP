#!/usr/bin/env python3
"""Render one ASCII triangular VTU mesh as a paper-ready PNG/PDF pair."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from paper_style import apply_paper_style
from plot_reference_epoch_meshes import draw_mesh, read_vtu


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--title", required=True)
    parser.add_argument("--mesh-label", default="")
    parser.add_argument("--linewidth", type=float, default=0.16)
    arguments = parser.parse_args()

    _, triangles = read_vtu(arguments.mesh)
    title = arguments.title + f"\n{triangles.shape[0]} cells"
    apply_paper_style()
    figure, ax = plt.subplots(figsize=(5.2, 5.2), constrained_layout=True)
    draw_mesh(ax, arguments.mesh, title, arguments.linewidth)
    if arguments.mesh_label:
        ax.set_ylabel(arguments.mesh_label, fontsize=9)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(arguments.output.with_suffix(".png"), dpi=300)
    figure.savefig(arguments.output.with_suffix(".pdf"))
    plt.close(figure)


if __name__ == "__main__":
    main()

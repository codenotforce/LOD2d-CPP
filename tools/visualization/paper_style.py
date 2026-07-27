"""Shared, headless Matplotlib style for LOD2d paper figures."""

from __future__ import annotations

import matplotlib as mpl


COLORS = {
    "fem": "#C44E52",
    "lod": "#4C72B0",
    "fine": "#6C6C6C",
    "energy": "#4C72B0",
    "l2": "#DD8452",
    "reference": "#8A8A8A",
}


def apply_paper_style() -> None:
    """Apply a compact style that works with Matplotlib 3.5 in WSL."""

    mpl.rcParams.update(
        {
            "figure.dpi": 120,
            "savefig.dpi": 600,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.04,
            "font.family": "serif",
            "font.size": 9.0,
            "axes.labelsize": 9.0,
            "axes.titlesize": 9.5,
            "legend.fontsize": 7.6,
            "xtick.labelsize": 8.0,
            "ytick.labelsize": 8.0,
            "lines.linewidth": 1.6,
            "lines.markersize": 4.5,
            "axes.linewidth": 0.8,
            "grid.linewidth": 0.45,
            "grid.alpha": 0.28,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.unicode_minus": False,
        }
    )

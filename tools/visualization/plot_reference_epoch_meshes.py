#!/usr/bin/env python3
"""Render real coarse/reference/candidate VTU checkpoints by epoch."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Circle

from paper_style import apply_paper_style


@dataclass(frozen=True)
class MeshEntry:
    epoch: int
    h_step: int
    iteration: int
    stage: str
    role: str
    filename: str
    cells: int
    dofs: int
    reference_mesh_version: int | None


def read_manifest(path: Path, epoch: int | None = None) -> list[MeshEntry]:
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    entries = [
        MeshEntry(
            epoch=int(row["epoch"]),
            h_step=int(row.get("H_step", row["iteration"])),
            iteration=int(row["iteration"]), stage=row["stage"],
            role=row["mesh_role"], filename=row["filename"],
            cells=int(row["N_cells"]), dofs=int(row["N_dofs"]),
            reference_mesh_version=(
                int(row["reference_mesh_version"])
                if row.get("reference_mesh_version") not in (None, "", "NA")
                else None
            ),
        )
        for row in rows if epoch is None or int(row["epoch"]) == epoch
    ]
    if not entries:
        raise ValueError(f"mesh manifest contains no entries for epoch={epoch}")
    return entries


def read_vtu(path: Path) -> tuple[np.ndarray, np.ndarray]:
    root = ET.parse(path).getroot()
    piece = root.find("./UnstructuredGrid/Piece")
    if piece is None:
        raise ValueError(f"missing UnstructuredGrid/Piece in {path}")
    points_array = piece.find("./Points/DataArray")
    connectivity_array = piece.find("./Cells/DataArray[@Name='connectivity']")
    offsets_array = piece.find("./Cells/DataArray[@Name='offsets']")
    types_array = piece.find("./Cells/DataArray[@Name='types']")
    if any(item is None for item in (points_array, connectivity_array, offsets_array, types_array)):
        raise ValueError(f"incomplete VTU arrays in {path}")
    if points_array.get("format") != "ascii":
        raise ValueError(f"only ASCII VTU is supported: {path}")
    points = np.fromstring(points_array.text or "", sep=" ", dtype=float).reshape(-1, 3)[:, :2]
    connectivity = np.fromstring(connectivity_array.text or "", sep=" ", dtype=int)
    offsets = np.fromstring(offsets_array.text or "", sep=" ", dtype=int)
    types = np.fromstring(types_array.text or "", sep=" ", dtype=int)
    if not np.all(types == 5) or not np.all(np.diff(np.r_[0, offsets]) == 3):
        raise ValueError(f"expected only VTK triangle cells in {path}")
    return points, connectivity.reshape(-1, 3)


def mesh_edges(points: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    edges = np.concatenate(
        (triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]]), axis=0
    )
    edges.sort(axis=1)
    return points[np.unique(edges, axis=0)]


def draw_mesh(ax, path: Path, title: str, linewidth: float) -> None:
    points, triangles = read_vtu(path)
    ax.add_collection(
        LineCollection(mesh_edges(points, triangles), colors="#26364a", linewidths=linewidth)
    )
    lower, upper = points.min(axis=0), points.max(axis=0)
    padding = 0.015 * max(*(upper - lower), 1.0)
    ax.set_xlim(lower[0] - padding, upper[0] + padding)
    ax.set_ylim(lower[1] - padding, upper[1] + padding)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(title, fontsize=8.5)
    ax.set_xticks([])
    ax.set_yticks([])


def select_checkpoints(
    entries: list[MeshEntry], count: int, all_checkpoints: bool,
) -> list[tuple[MeshEntry, MeshEntry, MeshEntry]]:
    checkpoint_entries = [
        entry for entry in entries
        if entry.stage in {"epoch_start", "committed", "pre_switch"}
    ]
    by_key = {(entry.iteration, entry.role): entry for entry in checkpoint_entries}
    iterations = sorted(
        iteration for iteration, role in by_key
        if role == "coarse"
        and (iteration, "reference") in by_key
        and (iteration, "candidate") in by_key
    )
    if not iterations:
        raise ValueError("no matching coarse/reference/candidate checkpoints")
    if all_checkpoints:
        selected = iterations
    else:
        indices = np.linspace(0, len(iterations) - 1, min(count, len(iterations)))
        selected = sorted({iterations[int(round(index))] for index in indices})
    return [
        (
            by_key[(iteration, "coarse")],
            by_key[(iteration, "reference")],
            by_key[(iteration, "candidate")],
        )
        for iteration in selected
    ]


def save_reference(
    run_dir: Path, output_dir: Path, entry: MeshEntry, experiment: str
) -> None:
    figure, ax = plt.subplots(figsize=(4.0, 4.0), constrained_layout=True)
    draw_mesh(
        ax, run_dir / entry.filename,
        rf"Fixed reference mesh $\mathcal{{T}}_h^0$ — {entry.cells} cells", 0.18,
    )
    stem = f"{experiment}_epoch{entry.epoch}_reference_mesh"
    figure.savefig(output_dir / f"{stem}.png", dpi=260)
    figure.savefig(output_dir / f"{stem}.pdf")
    plt.close(figure)


def save_evolution(
    run_dir: Path, output_dir: Path,
    triplets: list[tuple[MeshEntry, MeshEntry, MeshEntry]], role_index: int,
    experiment: str,
) -> None:
    role = ("coarse", "reference", "candidate")[role_index]
    figure, axes = plt.subplots(
        1, len(triplets),
        figsize=(3.15 * len(triplets), 3.15), constrained_layout=True
    )
    for ax, triplet in zip(np.atleast_1d(axes), triplets):
        entry = triplet[role_index]
        draw_mesh(
            ax, run_dir / entry.filename,
            rf"$H$-step {entry.h_step}, $i={entry.iteration}$" + "\n"
            + f"{entry.cells} cells",
            0.42 if role == "coarse" else 0.12,
        )
    figure.suptitle(
        f"{experiment} epoch {triplets[0][0].epoch}: {role} mesh evolution",
        fontsize=10,
    )
    stem = f"{experiment}_epoch{triplets[0][0].epoch}_{role}_evolution"
    figure.savefig(output_dir / f"{stem}.png", dpi=260)
    figure.savefig(output_dir / f"{stem}.pdf")
    plt.close(figure)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def save_triplet_pages(
    run_dir: Path,
    output_dir: Path,
    triplets: list[tuple[MeshEntry, MeshEntry, MeshEntry]],
    experiment: str,
    page_columns: int,
) -> None:
    if page_columns < 1:
        raise ValueError("--page-columns must be positive")
    reference_entries = [triplet[1] for triplet in triplets]
    reference_versions = {
        entry.reference_mesh_version for entry in reference_entries
        if entry.reference_mesh_version is not None
    }
    reference_hashes = {
        file_sha256(run_dir / entry.filename) for entry in reference_entries
    }
    if len(reference_versions) > 1 or len(reference_hashes) != 1:
        raise ValueError(
            "reference mesh changed inside the selected epoch: "
            f"versions={sorted(reference_versions)}, hashes={sorted(reference_hashes)}"
        )

    epoch = triplets[0][0].epoch
    stem = f"{experiment}_epoch{epoch}_coarse_reference_candidate_by_H_step"
    audit = {
        "experiment": experiment,
        "epoch": epoch,
        "H_steps": [triplet[0].h_step for triplet in triplets],
        "checkpoint_stages": [triplet[0].stage for triplet in triplets],
        "reference_mesh_versions": sorted(reference_versions),
        "reference_sha256": sorted(reference_hashes),
        "reference_unchanged_within_epoch": True,
    }
    (output_dir / f"{stem}.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="utf-8"
    )

    roles = ("Coarse $\\mathcal{T}_H$", "Reference $\\mathcal{T}_h$",
             "Candidate $\\mathcal{T}_c$")
    widths = (0.42, 0.10, 0.10)
    with PdfPages(output_dir / f"{stem}.pdf") as pdf:
        for page, start in enumerate(range(0, len(triplets), page_columns), 1):
            chunk = triplets[start:start + page_columns]
            figure, axes = plt.subplots(
                3, len(chunk), figsize=(3.0 * len(chunk), 8.4),
                constrained_layout=True, squeeze=False,
            )
            for column, triplet in enumerate(chunk):
                for row, entry in enumerate(triplet):
                    stage_label = entry.stage.replace("_", "-")
                    title = (
                        rf"$H$-step {entry.h_step}, $i={entry.iteration}$"
                        + "\n" + f"{stage_label}; {entry.cells} cells"
                        if row == 0 else f"{entry.cells} cells"
                    )
                    draw_mesh(
                        axes[row, column], run_dir / entry.filename,
                        title, widths[row],
                    )
                    if column == 0:
                        axes[row, column].set_ylabel(roles[row], fontsize=8.5)
            figure.suptitle(
                f"{experiment} epoch {epoch}: matched mesh snapshots "
                f"(page {page})",
                fontsize=10,
            )
            pdf.savefig(figure)
            figure.savefig(output_dir / f"{stem}_p{page:02d}.png", dpi=240)
            plt.close(figure)


def save_multi_epoch_triplet_pages(
    run_dir: Path,
    output_dir: Path,
    epoch_triplets: list[tuple[int, list[tuple[MeshEntry, MeshEntry, MeshEntry]]]],
    experiment: str,
    page_columns: int,
) -> None:
    """Render matched three-mesh checkpoints from several epochs in one figure."""
    if page_columns < 1:
        raise ValueError("--page-columns must be positive")
    if not epoch_triplets:
        raise ValueError("no epoch checkpoints selected")

    audit_epochs = []
    flattened: list[tuple[MeshEntry, MeshEntry, MeshEntry]] = []
    final_candidates: dict[int, MeshEntry] = {}
    initial_references: dict[int, MeshEntry] = {}
    for epoch, triplets in epoch_triplets:
        if not triplets:
            raise ValueError(f"epoch {epoch} has no matched checkpoints")
        reference_entries = [triplet[1] for triplet in triplets]
        reference_versions = {
            entry.reference_mesh_version for entry in reference_entries
            if entry.reference_mesh_version is not None
        }
        reference_hashes = {
            file_sha256(run_dir / entry.filename) for entry in reference_entries
        }
        if len(reference_versions) > 1 or len(reference_hashes) != 1:
            raise ValueError(
                f"reference mesh changed inside epoch {epoch}: "
                f"versions={sorted(reference_versions)}, "
                f"hashes={sorted(reference_hashes)}"
            )
        initial_references[epoch] = triplets[0][1]
        final_candidates[epoch] = triplets[-1][2]
        audit_epochs.append({
            "epoch": epoch,
            "H_steps": [triplet[0].h_step for triplet in triplets],
            "iterations": [triplet[0].iteration for triplet in triplets],
            "checkpoint_stages": [triplet[0].stage for triplet in triplets],
            "reference_mesh_versions": sorted(reference_versions),
            "reference_sha256": sorted(reference_hashes),
            "reference_unchanged_within_epoch": True,
        })
        flattened.extend(triplets)

    epochs = [epoch for epoch, _ in epoch_triplets]
    promotions = []
    for previous, current in zip(epochs, epochs[1:]):
        candidate_hash = file_sha256(run_dir / final_candidates[previous].filename)
        reference_hash = file_sha256(run_dir / initial_references[current].filename)
        promotions.append({
            "from_epoch": previous,
            "to_epoch": current,
            "candidate_sha256": candidate_hash,
            "next_reference_sha256": reference_hash,
            "candidate_promoted_to_next_reference": candidate_hash == reference_hash,
        })

    epoch_token = "-".join(str(epoch) for epoch in epochs)
    stem = (
        f"{experiment}_epochs{epoch_token}_"
        "coarse_reference_candidate_by_H_step"
    )
    audit = {
        "experiment": experiment,
        "epochs": audit_epochs,
        "epoch_boundary_promotions": promotions,
    }
    (output_dir / f"{stem}.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="utf-8"
    )

    roles = ("Coarse $\\mathcal{T}_H$", "Reference $\\mathcal{T}_h$",
             "Candidate $\\mathcal{T}_c$")
    widths = (0.42, 0.10, 0.10)
    with PdfPages(output_dir / f"{stem}.pdf") as pdf:
        for page, start in enumerate(range(0, len(flattened), page_columns), 1):
            chunk = flattened[start:start + page_columns]
            figure, axes = plt.subplots(
                3, len(chunk), figsize=(3.0 * len(chunk), 8.4),
                constrained_layout=True, squeeze=False,
            )
            for column, triplet in enumerate(chunk):
                for row, entry in enumerate(triplet):
                    stage_label = entry.stage.replace("_", "-")
                    title = (
                        rf"Epoch {entry.epoch}, $H$-step {entry.h_step}, "
                        rf"$i={entry.iteration}$"
                        + "\n" + f"{stage_label}; {entry.cells} cells"
                        if row == 0 else f"{entry.cells} cells"
                    )
                    draw_mesh(
                        axes[row, column], run_dir / entry.filename,
                        title, widths[row],
                    )
                    if column == 0:
                        axes[row, column].set_ylabel(roles[row], fontsize=8.5)
            epoch_range = (
                str(epochs[0]) if len(epochs) == 1
                else f"{epochs[0]}--{epochs[-1]}"
            )
            figure.suptitle(
                f"{experiment} epochs {epoch_range}: matched mesh snapshots "
                f"(page {page})",
                fontsize=10,
            )
            pdf.savefig(figure)
            figure.savefig(output_dir / f"{stem}_p{page:02d}.png", dpi=240)
            plt.close(figure)


def save_unique_epoch_meshes(
    run_dir: Path,
    output_dir: Path,
    epoch_triplets: list[tuple[int, list[tuple[MeshEntry, MeshEntry, MeshEntry]]]],
    experiment: str,
) -> None:
    """Render each distinct mesh once per epoch while retaining all role labels."""
    role_names = ("Coarse $\\mathcal{T}_H$", "Reference $\\mathcal{T}_h$",
                  "Candidate $\\mathcal{T}_c$")
    role_keys = ("coarse", "reference", "candidate")
    widths = (0.42, 0.10, 0.10)
    audit = {"experiment": experiment, "epochs": []}

    pdf_path = output_dir / (
        f"{experiment}_epochs" + "-".join(str(epoch) for epoch, _ in epoch_triplets)
        + "_unique_coarse_reference_candidate.pdf"
    )
    with PdfPages(pdf_path) as pdf:
        for epoch, triplets in epoch_triplets:
            role_entries: list[list[MeshEntry]] = [[], [], []]
            role_hashes: list[list[str]] = [[], [], []]
            suppressed: list[dict[str, object]] = []

            for role_index in range(3):
                seen: dict[str, MeshEntry] = {}
                for triplet in triplets:
                    entry = triplet[role_index]
                    digest = file_sha256(run_dir / entry.filename)
                    if digest in seen:
                        suppressed.append({
                            "role": role_keys[role_index],
                            "filename": entry.filename,
                            "same_as": seen[digest].filename,
                            "reason": "identical within role",
                        })
                        continue
                    seen[digest] = entry
                    role_entries[role_index].append(entry)
                    role_hashes[role_index].append(digest)

            # At epoch start the candidate commonly equals the fixed reference.
            # Show that geometry only in the reference row and record the alias.
            reference_by_hash = {
                digest: entry
                for digest, entry in zip(role_hashes[1], role_entries[1])
            }
            filtered_candidates: list[MeshEntry] = []
            filtered_candidate_hashes: list[str] = []
            for digest, entry in zip(role_hashes[2], role_entries[2]):
                if digest in reference_by_hash:
                    suppressed.append({
                        "role": "candidate",
                        "filename": entry.filename,
                        "same_as": reference_by_hash[digest].filename,
                        "reason": "candidate geometry identical to reference",
                    })
                else:
                    filtered_candidates.append(entry)
                    filtered_candidate_hashes.append(digest)
            role_entries[2] = filtered_candidates
            role_hashes[2] = filtered_candidate_hashes

            columns = max(max(len(entries), 1) for entries in role_entries)
            figure, axes = plt.subplots(
                3, columns, figsize=(3.15 * columns, 8.4),
                constrained_layout=True, squeeze=False,
            )
            for role_index, entries in enumerate(role_entries):
                for column in range(columns):
                    ax = axes[role_index, column]
                    if column >= len(entries):
                        ax.axis("off")
                        note = None
                        if role_index == 1:
                            note = "Reference unchanged within epoch"
                        elif role_index == 2:
                            note = (
                                "Epoch-start candidate = reference\n"
                                "(identical mesh not repeated)"
                            )
                        if note is not None and column == len(entries):
                            ax.text(
                                0.5, 0.5, note,
                                ha="center", va="center", fontsize=8.5,
                                transform=ax.transAxes,
                            )
                        continue
                    entry = entries[column]
                    stage_label = entry.stage.replace("_", "-")
                    draw_mesh(
                        ax, run_dir / entry.filename,
                        rf"$H$-step {entry.h_step}, $i={entry.iteration}$"
                        + "\n" + f"{stage_label}; {entry.cells} cells",
                        widths[role_index],
                    )
                axes[role_index, 0].set_ylabel(role_names[role_index], fontsize=8.5)

            figure.suptitle(
                f"{experiment} epoch {epoch}: distinct mesh states",
                fontsize=10,
            )
            stem = (
                f"{experiment}_epoch{epoch}_unique_"
                "coarse_reference_candidate"
            )
            figure.savefig(output_dir / f"{stem}.png", dpi=260)
            figure.savefig(output_dir / f"{stem}.pdf")
            pdf.savefig(figure)
            plt.close(figure)
            audit["epochs"].append({
                "epoch": epoch,
                "distinct_by_role": {
                    role: [
                        {
                            "filename": entry.filename,
                            "stage": entry.stage,
                            "H_step": entry.h_step,
                            "iteration": entry.iteration,
                            "cells": entry.cells,
                            "sha256": digest,
                        }
                        for entry, digest in zip(entries, hashes)
                    ]
                    for role, entries, hashes in zip(
                        role_keys, role_entries, role_hashes
                    )
                },
                "suppressed_identical_snapshots": suppressed,
            })

    json_path = output_dir / (
        f"{experiment}_epochs" + "-".join(str(epoch) for epoch, _ in epoch_triplets)
        + "_unique_coarse_reference_candidate.json"
    )
    json_path.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")


def hybrid_physical_radius(run_dir: Path) -> float | None:
    run_json = run_dir / "run.json"
    if not run_json.exists():
        return None
    value = json.loads(run_json.read_text(encoding="utf-8")).get(
        "config", {}
    ).get("hybrid_minimum_physical_radius")
    if value is None or float(value) <= 0.0:
        return None
    return float(value)


def save_hybrid_final(run_dir: Path, output_dir: Path, entry: MeshEntry) -> None:
    path = run_dir / entry.filename
    points, triangles = read_vtu(path)
    root = ET.parse(path).getroot()
    field = root.find(
        "./UnstructuredGrid/Piece/CellData/DataArray[@Name='hybrid_region']"
    )
    if field is None:
        raise ValueError(f"hybrid_region cell field is missing in {path}")
    regions = np.fromstring(field.text or "", sep=" ", dtype=int)
    if regions.size != triangles.shape[0]:
        raise ValueError(f"hybrid_region has the wrong size in {path}")
    colors = np.asarray(["#f7f7f7", "#9ecae1", "#fdae6b"])[regions]
    physical_radius = hybrid_physical_radius(run_dir)
    figure, axes = plt.subplots(1, 2, figsize=(7.0, 3.35), constrained_layout=True)
    for ax, zoom in zip(axes, (False, True)):
        ax.add_collection(PolyCollection(
            points[triangles], facecolors=colors, edgecolors="#26364a",
            linewidths=0.12,
        ))
        ax.set_aspect("equal", adjustable="box")
        if physical_radius is not None:
            ax.add_patch(Circle(
                (0.0, 0.0), physical_radius, fill=False,
                edgecolor="#7b3294", linewidth=1.0, linestyle="--", zorder=5,
            ))
        if zoom:
            ax.set_xlim(-0.35, 0.35)
            ax.set_ylim(-0.35, 0.35)
            ax.set_title("Reentrant-corner zoom", fontsize=8.5)
        else:
            lower, upper = points.min(axis=0), points.max(axis=0)
            ax.set_xlim(lower[0] - 0.02, upper[0] + 0.02)
            ax.set_ylim(lower[1] - 0.02, upper[1] + 0.02)
            ax.set_title("Final solved hybrid mesh", fontsize=8.5)
        ax.set_xticks([])
        ax.set_yticks([])
    radius_text = (
        rf", $B_{{R_*}}$ dashed ($R_*={physical_radius:g}$)"
        if physical_radius is not None else ""
    )
    figure.suptitle(
        r"E2 matching regions: regular (white), $\Omega_F$ (blue), "
        r"$\Omega_S$ (orange)" + radius_text,
        fontsize=9.5,
    )
    figure.savefig(output_dir / "E2_final_hybrid_matching_mesh.png", dpi=260)
    figure.savefig(output_dir / "E2_final_hybrid_matching_mesh.pdf")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--checkpoints", type=int, default=4)
    parser.add_argument(
        "--all-checkpoints", action="store_true",
        help="render every saved H-step rather than an evenly spaced subset",
    )
    parser.add_argument("--page-columns", type=int, default=4)
    epoch_group = parser.add_mutually_exclusive_group()
    epoch_group.add_argument("--epoch", type=int, default=None)
    epoch_group.add_argument(
        "--epochs",
        help="comma-separated epochs rendered together, for example 0,1",
    )
    parser.add_argument("--experiment-label", default="E1")
    parser.add_argument(
        "--deduplicate-identical", action="store_true",
        help=(
            "write one figure per epoch and suppress repeated mesh geometries; "
            "candidate snapshots identical to the fixed reference are recorded "
            "in the JSON audit but not drawn twice"
        ),
    )
    arguments = parser.parse_args()
    if not arguments.all_checkpoints and arguments.checkpoints < 2:
        raise ValueError("--checkpoints must be at least two")
    manifest = arguments.run_dir / "mesh_manifest.csv"
    if arguments.epochs:
        epochs = [int(value.strip()) for value in arguments.epochs.split(",")]
        if not epochs or len(set(epochs)) != len(epochs):
            raise ValueError("--epochs must contain distinct epoch numbers")
    else:
        epochs = [0 if arguments.epoch is None else arguments.epoch]

    epoch_triplets = []
    for epoch in epochs:
        epoch_entries = read_manifest(manifest, epoch)
        epoch_triplets.append((
            epoch,
            select_checkpoints(
                epoch_entries, arguments.checkpoints, arguments.all_checkpoints
            ),
        ))

    entries = read_manifest(manifest, epochs[0])
    references = [entry for entry in entries if entry.role == "reference"]
    if not references:
        raise ValueError(f"epoch {epochs[0]} reference mesh is missing")
    triplets = epoch_triplets[0][1]
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    apply_paper_style()
    if arguments.deduplicate_identical:
        save_unique_epoch_meshes(
            arguments.run_dir, arguments.output_dir, epoch_triplets,
            arguments.experiment_label,
        )
    elif len(epoch_triplets) == 1:
        save_reference(
            arguments.run_dir, arguments.output_dir, references[0],
            arguments.experiment_label,
        )
        save_evolution(
            arguments.run_dir, arguments.output_dir, triplets, 0,
            arguments.experiment_label,
        )
        save_evolution(
            arguments.run_dir, arguments.output_dir, triplets, 1,
            arguments.experiment_label,
        )
        save_evolution(
            arguments.run_dir, arguments.output_dir, triplets, 2,
            arguments.experiment_label,
        )
        save_triplet_pages(
            arguments.run_dir, arguments.output_dir, triplets,
            arguments.experiment_label, arguments.page_columns,
        )
    else:
        save_multi_epoch_triplet_pages(
            arguments.run_dir, arguments.output_dir, epoch_triplets,
            arguments.experiment_label, arguments.page_columns,
        )
    if arguments.experiment_label == "E2":
        hybrid = [entry for entry in read_manifest(manifest) if entry.role == "hybrid_regions"]
        if hybrid:
            save_hybrid_final(arguments.run_dir, arguments.output_dir, hybrid[-1])


if __name__ == "__main__":
    main()

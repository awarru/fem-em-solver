"""Phantom-region field sampling, metrics, and export helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Any
import csv
import json

import numpy as np
from mpi4py import MPI


def _tagged_cells(cell_tags, tag: int) -> np.ndarray:
    cells = cell_tags.indices[cell_tags.values == int(tag)]
    return np.asarray(cells, dtype=np.int32)


def _cell_centroids(mesh, cells: np.ndarray) -> np.ndarray:
    if cells.size == 0:
        return np.zeros((0, 3), dtype=np.float64)

    gdim = mesh.geometry.dim
    x = mesh.geometry.x
    dofmap = mesh.geometry.dofmap

    centroids = np.zeros((cells.size, 3), dtype=np.float64)
    for i, cell in enumerate(cells):
        coords = x[dofmap[int(cell)], :gdim]
        centroids[i, :gdim] = np.mean(coords, axis=0)
    return centroids


def _evaluate_on_cells(field, points: np.ndarray, cells: np.ndarray) -> np.ndarray:
    value_shape = field.function_space.element.value_shape
    value_size = int(np.prod(value_shape)) if value_shape else 1

    if points.shape[0] == 0:
        return np.zeros((0, value_size), dtype=np.float64)

    values = field.eval(points, cells)
    values = np.asarray(values, dtype=np.float64)
    if values.ndim == 1:
        values = values.reshape(-1, 1)
    return values


def compute_tagged_vector_magnitude_stats(field, cell_tags, tag: int, *, comm: MPI.Intracomm | None = None) -> dict[str, float]:
    """Compute global min/max/mean |field| over centroids of cells with a given tag."""
    mesh = field.function_space.mesh
    comm = comm or mesh.comm

    cells = _tagged_cells(cell_tags, tag)
    points = _cell_centroids(mesh, cells)
    values = _evaluate_on_cells(field, points, cells)

    if values.shape[0] > 0:
        mags = np.linalg.norm(values, axis=1)
        local_count = int(mags.size)
        local_sum = float(np.sum(mags))
        local_min = float(np.min(mags))
        local_max = float(np.max(mags))
    else:
        local_count = 0
        local_sum = 0.0
        local_min = float("inf")
        local_max = float("-inf")

    global_count = comm.allreduce(local_count, op=MPI.SUM)
    global_sum = comm.allreduce(local_sum, op=MPI.SUM)
    global_min = comm.allreduce(local_min, op=MPI.MIN)
    global_max = comm.allreduce(local_max, op=MPI.MAX)

    if global_count == 0:
        raise ValueError(f"No cells found for tag={tag}; cannot compute tagged field metrics")

    return {
        "count": int(global_count),
        "min": float(global_min),
        "max": float(global_max),
        "mean": float(global_sum / global_count),
    }


def export_tagged_field_samples_csv(
    field,
    cell_tags,
    tag: int,
    output_path: str | Path,
    *,
    comm: MPI.Intracomm | None = None,
) -> Path | None:
    """Export centroid samples for a tagged region to CSV on rank 0.

    CSV columns: x,y,z,fx,fy,fz,mag
    """
    mesh = field.function_space.mesh
    comm = comm or mesh.comm

    cells = _tagged_cells(cell_tags, tag)
    points = _cell_centroids(mesh, cells)
    values = _evaluate_on_cells(field, points, cells)

    gathered_points = comm.gather(points, root=0)
    gathered_values = comm.gather(values, root=0)

    if comm.rank != 0:
        return None

    all_points = np.vstack(gathered_points) if gathered_points else np.zeros((0, 3), dtype=np.float64)
    all_values = np.vstack(gathered_values) if gathered_values else np.zeros((0, 3), dtype=np.float64)

    if all_points.shape[0] == 0:
        raise ValueError(f"No tagged samples available for tag={tag}; cannot export CSV")

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["x", "y", "z", "fx", "fy", "fz", "mag"])
        for point, vec in zip(all_points, all_values):
            mag = float(np.linalg.norm(vec))
            vec_pad = np.zeros(3, dtype=np.float64)
            vec_pad[: min(3, vec.size)] = vec[: min(3, vec.size)]
            writer.writerow([
                float(point[0]),
                float(point[1]),
                float(point[2]),
                float(vec_pad[0]),
                float(vec_pad[1]),
                float(vec_pad[2]),
                mag,
            ])

    return output_path


def compute_phantom_eb_metrics_and_export(
    e_field,
    b_field,
    cell_tags,
    *,
    phantom_tag: int = 3,
    output_dir: str | Path = "paraview_output",
    basename: str = "phantom_fields",
    comm: MPI.Intracomm | None = None,
) -> dict[str, Any]:
    """Compute phantom |E|/|B| stats and export phantom-only centroid samples."""
    mesh = e_field.function_space.mesh
    comm = comm or mesh.comm

    e_stats = compute_tagged_vector_magnitude_stats(e_field, cell_tags, phantom_tag, comm=comm)
    b_stats = compute_tagged_vector_magnitude_stats(b_field, cell_tags, phantom_tag, comm=comm)

    output_dir = Path(output_dir)
    e_csv = export_tagged_field_samples_csv(
        e_field,
        cell_tags,
        phantom_tag,
        output_dir / f"{basename}_phantom_E_samples.csv",
        comm=comm,
    )
    b_csv = export_tagged_field_samples_csv(
        b_field,
        cell_tags,
        phantom_tag,
        output_dir / f"{basename}_phantom_B_samples.csv",
        comm=comm,
    )

    summary = {
        "phantom_tag": int(phantom_tag),
        "E_magnitude": e_stats,
        "B_magnitude": b_stats,
        "exports": {
            "E_csv": str(e_csv) if e_csv is not None else None,
            "B_csv": str(b_csv) if b_csv is not None else None,
        },
    }

    if comm.rank == 0:
        summary_path = output_dir / f"{basename}_phantom_metrics.json"
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
        summary["exports"]["summary_json"] = str(summary_path)
    else:
        summary["exports"]["summary_json"] = None

    return summary

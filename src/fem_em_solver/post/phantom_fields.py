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


def _interior_tagged_cells(mesh, cell_tags, tag: int) -> np.ndarray:
    """Return tagged cells that are at least one facet away from non-tagged cells."""
    cells = _tagged_cells(cell_tags, tag)
    if cells.size == 0:
        return cells

    topology = mesh.topology
    tdim = topology.dim
    fdim = tdim - 1
    topology.create_connectivity(tdim, fdim)
    topology.create_connectivity(fdim, tdim)
    c_to_f = topology.connectivity(tdim, fdim)
    f_to_c = topology.connectivity(fdim, tdim)

    tag_lookup = {int(cell): int(value) for cell, value in zip(cell_tags.indices, cell_tags.values)}

    interior_cells: list[int] = []
    for cell in cells:
        c = int(cell)
        is_boundary_adjacent = False
        for facet in c_to_f.links(c):
            adjacent_cells = f_to_c.links(int(facet))
            for nbr in adjacent_cells:
                nbr_cell = int(nbr)
                if nbr_cell == c:
                    continue
                if tag_lookup.get(nbr_cell, None) != int(tag):
                    is_boundary_adjacent = True
                    break
            if is_boundary_adjacent:
                break

        if not is_boundary_adjacent:
            interior_cells.append(c)

    return np.asarray(interior_cells, dtype=np.int32)


def _evaluate_on_cells(field, points: np.ndarray, cells: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    value_shape = field.function_space.element.value_shape
    value_size = int(np.prod(value_shape)) if value_shape else 1

    if points.shape[0] == 0:
        return (
            np.zeros((0, value_size), dtype=np.float64),
            np.zeros((0, 3), dtype=np.float64),
            np.zeros((0,), dtype=np.int32),
            0,
        )

    try:
        values = field.eval(points, cells)
        values = np.asarray(values, dtype=np.float64)
        if values.ndim == 1:
            values = values.reshape(-1, 1)
        return values, np.asarray(points, dtype=np.float64), np.asarray(cells, dtype=np.int32), 0
    except Exception:
        # Fallback: evaluate point-by-point and skip invalid cell-point pairs.
        kept_values: list[np.ndarray] = []
        kept_points: list[np.ndarray] = []
        kept_cells: list[int] = []
        invalid_count = 0

        for point, cell in zip(points, cells):
            try:
                one = field.eval(np.asarray([point], dtype=np.float64), np.asarray([cell], dtype=np.int32))
                one = np.asarray(one, dtype=np.float64).reshape(-1)
                kept_values.append(one)
                kept_points.append(np.asarray(point, dtype=np.float64))
                kept_cells.append(int(cell))
            except Exception:
                invalid_count += 1

        if not kept_values:
            return (
                np.zeros((0, value_size), dtype=np.float64),
                np.zeros((0, 3), dtype=np.float64),
                np.zeros((0,), dtype=np.int32),
                invalid_count,
            )

        return (
            np.vstack(kept_values),
            np.vstack(kept_points),
            np.asarray(kept_cells, dtype=np.int32),
            invalid_count,
        )


def _sampling_cells_with_interface_guardrails(mesh, cell_tags, tag: int, *, prefer_interior: bool = True) -> tuple[np.ndarray, int]:
    tagged_cells = _tagged_cells(cell_tags, tag)
    if tagged_cells.size == 0:
        return tagged_cells, 0

    if not prefer_interior:
        return tagged_cells, 0

    interior_cells = _interior_tagged_cells(mesh, cell_tags, tag)
    dropped = int(tagged_cells.size - interior_cells.size)

    if interior_cells.size > 0:
        return interior_cells, dropped

    # Guardrail fallback: do not fail if every tagged cell touches an interface.
    return tagged_cells, 0


def compute_tagged_vector_magnitude_stats(
    field,
    cell_tags,
    tag: int,
    *,
    comm: MPI.Intracomm | None = None,
    prefer_interior_samples: bool = True,
) -> dict[str, float]:
    """Compute global min/max/mean |field| over robust tagged centroid samples."""
    mesh = field.function_space.mesh
    comm = comm or mesh.comm

    requested_cells = _tagged_cells(cell_tags, tag)
    sampling_cells, boundary_dropped_local = _sampling_cells_with_interface_guardrails(
        mesh,
        cell_tags,
        tag,
        prefer_interior=prefer_interior_samples,
    )
    points = _cell_centroids(mesh, sampling_cells)
    values, _valid_points, valid_cells, invalid_samples_local = _evaluate_on_cells(field, points, sampling_cells)

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

    global_requested = comm.allreduce(int(requested_cells.size), op=MPI.SUM)
    global_sampling = comm.allreduce(int(sampling_cells.size), op=MPI.SUM)
    global_boundary_dropped = comm.allreduce(int(boundary_dropped_local), op=MPI.SUM)
    global_invalid = comm.allreduce(int(invalid_samples_local), op=MPI.SUM)
    global_valid_cells = comm.allreduce(int(valid_cells.size), op=MPI.SUM)

    if global_count == 0:
        raise ValueError(
            f"No valid centroid samples found for tag={tag}; requested={global_requested}, "
            f"sampling={global_sampling}, invalid={global_invalid}"
        )

    return {
        "count": int(global_count),
        "min": float(global_min),
        "max": float(global_max),
        "mean": float(global_sum / global_count),
        "requested_cells": int(global_requested),
        "sampling_cells": int(global_sampling),
        "valid_sample_cells": int(global_valid_cells),
        "boundary_adjacent_cells_dropped": int(global_boundary_dropped),
        "invalid_samples_dropped": int(global_invalid),
    }


def export_tagged_field_samples_csv(
    field,
    cell_tags,
    tag: int,
    output_path: str | Path,
    *,
    comm: MPI.Intracomm | None = None,
    prefer_interior_samples: bool = True,
) -> Path | None:
    """Export centroid samples for a tagged region to CSV on rank 0.

    CSV columns: x,y,z,fx,fy,fz,mag
    """
    mesh = field.function_space.mesh
    comm = comm or mesh.comm

    sampling_cells, _ = _sampling_cells_with_interface_guardrails(
        mesh,
        cell_tags,
        tag,
        prefer_interior=prefer_interior_samples,
    )
    points = _cell_centroids(mesh, sampling_cells)
    values, valid_points, valid_cells, _ = _evaluate_on_cells(field, points, sampling_cells)

    gathered_points = comm.gather(valid_points, root=0)
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
    prefer_interior_samples: bool = True,
) -> dict[str, Any]:
    """Compute phantom |E|/|B| stats and export phantom-only centroid samples."""
    mesh = e_field.function_space.mesh
    comm = comm or mesh.comm

    e_stats = compute_tagged_vector_magnitude_stats(
        e_field,
        cell_tags,
        phantom_tag,
        comm=comm,
        prefer_interior_samples=prefer_interior_samples,
    )
    b_stats = compute_tagged_vector_magnitude_stats(
        b_field,
        cell_tags,
        phantom_tag,
        comm=comm,
        prefer_interior_samples=prefer_interior_samples,
    )

    output_dir = Path(output_dir)
    e_csv = export_tagged_field_samples_csv(
        e_field,
        cell_tags,
        phantom_tag,
        output_dir / f"{basename}_phantom_E_samples.csv",
        comm=comm,
        prefer_interior_samples=prefer_interior_samples,
    )
    b_csv = export_tagged_field_samples_csv(
        b_field,
        cell_tags,
        phantom_tag,
        output_dir / f"{basename}_phantom_B_samples.csv",
        comm=comm,
        prefer_interior_samples=prefer_interior_samples,
    )

    summary = {
        "phantom_tag": int(phantom_tag),
        "sampling": {
            "prefer_interior_samples": bool(prefer_interior_samples),
            "requested_cells": int(e_stats["requested_cells"]),
            "sampling_cells": int(e_stats["sampling_cells"]),
            "valid_sample_cells": int(e_stats["valid_sample_cells"]),
            "boundary_adjacent_cells_dropped": int(e_stats["boundary_adjacent_cells_dropped"]),
            "invalid_samples_dropped_e": int(e_stats["invalid_samples_dropped"]),
            "invalid_samples_dropped_b": int(b_stats["invalid_samples_dropped"]),
        },
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

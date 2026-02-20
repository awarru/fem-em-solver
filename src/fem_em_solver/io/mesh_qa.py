"""Mesh QA helpers for quick tag-integrity diagnostics."""

from __future__ import annotations

from typing import Dict, Mapping, Optional

import numpy as np
from mpi4py import MPI


def cell_tag_counts(cell_tags, comm: Optional[MPI.Intracomm] = None) -> Dict[int, int]:
    """Return global cell counts per cell-tag value.

    Parameters
    ----------
    cell_tags : dolfinx.mesh.MeshTags
        Cell tag data returned from mesh generation.
    comm : MPI.Intracomm, optional
        MPI communicator. Defaults to ``cell_tags.topology.comm`` when available.
    """
    if comm is None:
        comm = cell_tags.topology.comm

    local_values = np.asarray(cell_tags.values, dtype=np.int64)

    local_max = np.array([-1], dtype=np.int64)
    if local_values.size:
        local_max[0] = int(np.max(local_values))

    global_max = np.array([-1], dtype=np.int64)
    comm.Allreduce(local_max, global_max, op=MPI.MAX)

    if global_max[0] < 0:
        return {}

    local_hist = np.zeros(int(global_max[0]) + 1, dtype=np.int64)
    if local_values.size:
        local_hist = np.bincount(local_values, minlength=local_hist.size).astype(np.int64)

    global_hist = np.zeros_like(local_hist)
    comm.Allreduce(local_hist, global_hist, op=MPI.SUM)

    return {int(tag): int(count) for tag, count in enumerate(global_hist) if count > 0}


def format_cell_tag_summary(
    counts: Mapping[int, int],
    tag_names: Optional[Mapping[int, str]] = None,
) -> str:
    """Build a compact deterministic summary string for tag counts."""
    if not counts:
        return "(no tagged cells)"

    parts = []
    for tag in sorted(counts):
        name = tag_names.get(tag, f"tag_{tag}") if tag_names else f"tag_{tag}"
        parts.append(f"{name}={counts[tag]}")

    return ", ".join(parts)


def print_cell_tag_summary(
    cell_tags,
    tag_names: Optional[Mapping[int, str]] = None,
    comm: Optional[MPI.Intracomm] = None,
    prefix: str = "[mesh] ",
) -> Dict[int, int]:
    """Print global cell counts by tag on rank 0 and return the counts."""
    if comm is None:
        comm = cell_tags.topology.comm

    counts = cell_tag_counts(cell_tags, comm=comm)

    if comm.rank == 0:
        print(f"{prefix}cell-tag counts: {format_cell_tag_summary(counts, tag_names=tag_names)}")

    return counts

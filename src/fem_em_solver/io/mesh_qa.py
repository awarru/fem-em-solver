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


def format_expected_tag_counts(
    counts: Mapping[int, int],
    required_tags: Mapping[int, str],
    minimum_count: int = 1,
) -> str:
    """Format expected-vs-actual counts for required tags.

    Parameters
    ----------
    counts : Mapping[int, int]
        Observed global cell counts by tag.
    required_tags : Mapping[int, str]
        Required tags and display names.
    minimum_count : int, default=1
        Minimum required count for each required tag.
    """
    if not required_tags:
        return "(no required tags)"

    rendered = []
    for tag in sorted(required_tags):
        name = required_tags[tag]
        actual = int(counts.get(tag, 0))
        status = "OK" if actual >= minimum_count else "MISSING"
        rendered.append(
            f"{name}(tag={tag}): expected>={minimum_count}, actual={actual} [{status}]"
        )

    return "; ".join(rendered)


def print_required_tag_failure_summary(
    counts: Mapping[int, int],
    required_tags: Mapping[int, str],
    comm: Optional[MPI.Intracomm] = None,
    prefix: str = "[mesh-qa] ",
) -> None:
    """Print compact diagnostics for required-tag failures on rank 0.

    Includes expected-vs-actual required tag counts and a compact summary of
    all observed tags to speed up debugging for missing/empty tag failures.
    """
    if comm is not None and comm.rank != 0:
        return

    print(
        f"{prefix}required-tag expected vs actual: "
        f"{format_expected_tag_counts(counts, required_tags)}"
    )
    print(
        f"{prefix}observed-tag summary: "
        f"{format_cell_tag_summary(counts, tag_names=required_tags)}"
    )

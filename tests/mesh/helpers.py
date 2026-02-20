"""Test helpers for mesh tag integrity checks."""

from __future__ import annotations

from typing import Mapping

import numpy as np
from mpi4py import MPI

from fem_em_solver.io.mesh_qa import cell_tag_counts


REQUIRED_COIL_PHANTOM_TAGS = {
    1: "coil_1",
    2: "coil_2",
    3: "phantom",
    4: "air",
}


def assert_required_tags_nonempty(cell_tags, required_tags: Mapping[int, str], comm: MPI.Intracomm) -> None:
    """Assert that each required tag exists and has at least one cell."""
    counts = cell_tag_counts(cell_tags, comm=comm)

    missing = []
    for tag, name in required_tags.items():
        if counts.get(tag, 0) <= 0:
            missing.append(f"{name} (tag={tag})")

    assert not missing, f"Required mesh tags missing/empty: {', '.join(missing)}"


def compute_tag_cell_centroid(mesh, cell_tags, tag: int, comm: MPI.Intracomm) -> np.ndarray:
    """Compute global centroid of cell centroids for a given tag."""
    tagged_cells = cell_tags.indices[cell_tags.values == tag]

    local_sum = np.zeros(3, dtype=np.float64)
    local_count = np.array([0], dtype=np.int64)

    if tagged_cells.size:
        dofmap = mesh.geometry.dofmap
        coords = mesh.geometry.x
        centroids = np.array([coords[dofmap[cell]].mean(axis=0) for cell in tagged_cells])
        local_sum = centroids.sum(axis=0)
        local_count[0] = centroids.shape[0]

    global_sum = np.zeros(3, dtype=np.float64)
    global_count = np.array([0], dtype=np.int64)

    comm.Allreduce(local_sum, global_sum, op=MPI.SUM)
    comm.Allreduce(local_count, global_count, op=MPI.SUM)

    if global_count[0] == 0:
        raise AssertionError(f"Tag {tag} has no cells; centroid is undefined")

    return global_sum / float(global_count[0])


def assert_tags_distinct_by_centroid(
    mesh,
    cell_tags,
    tag_a: int,
    tag_b: int,
    comm: MPI.Intracomm,
    min_distance: float = 1.0e-4,
) -> None:
    """Assert that two tagged regions are spatially distinct."""
    centroid_a = compute_tag_cell_centroid(mesh, cell_tags, tag_a, comm=comm)
    centroid_b = compute_tag_cell_centroid(mesh, cell_tags, tag_b, comm=comm)

    distance = np.linalg.norm(centroid_a - centroid_b)
    assert distance > min_distance, (
        f"Tagged regions should be distinct but centroids are too close: "
        f"tag {tag_a} vs {tag_b}, distance={distance:.3e}"
    )

"""Utilities for robust field evaluation on distributed meshes."""

from __future__ import annotations

import numpy as np
from mpi4py import MPI
from dolfinx import geometry


def evaluate_vector_field_parallel(function, points: np.ndarray, comm: MPI.Intracomm | None = None):
    """Evaluate a vector-valued dolfinx function at points across MPI ranks.

    Parameters
    ----------
    function:
        dolfinx function to evaluate.
    points:
        Array with shape (N, 3).
    comm:
        MPI communicator. Defaults to ``function.function_space.mesh.comm``.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        values:
            Array of shape (N, value_size) with point evaluations.
            Invalid/outside points are filled with zeros.
        valid_mask:
            Boolean mask with shape (N,) marking which points were evaluated.
    """
    mesh = function.function_space.mesh
    comm = comm or mesh.comm

    points = np.asarray(points, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError(f"points must have shape (N, 3), got {points.shape}")

    n_points = points.shape[0]
    value_shape = function.function_space.element.value_shape
    value_size = int(np.prod(value_shape)) if len(value_shape) > 0 else 1

    bb_tree = geometry.bb_tree(mesh, mesh.topology.dim)
    candidate_cells = geometry.compute_collisions_points(bb_tree, points)
    colliding_cells = geometry.compute_colliding_cells(mesh, candidate_cells, points)

    local_point_idx: list[int] = []
    local_cells: list[int] = []
    for i in range(n_points):
        links = colliding_cells.links(i)
        if len(links) > 0:
            local_point_idx.append(i)
            local_cells.append(links[0])

    local_point_idx_arr = np.asarray(local_point_idx, dtype=np.int32)

    if local_point_idx_arr.size > 0:
        local_points = points[local_point_idx_arr]
        local_cells_arr = np.asarray(local_cells, dtype=np.int32)
        local_values = function.eval(local_points, local_cells_arr)
    else:
        local_values = np.zeros((0, value_size), dtype=np.float64)

    gathered_values = comm.gather(local_values, root=0)
    gathered_indices = comm.gather(local_point_idx_arr, root=0)

    if comm.rank == 0:
        values = np.zeros((n_points, value_size), dtype=np.float64)
        valid_mask = np.zeros(n_points, dtype=bool)
        for rank_values, rank_indices in zip(gathered_values, gathered_indices):
            if rank_indices.size > 0:
                values[rank_indices] = rank_values
                valid_mask[rank_indices] = True
    else:
        values = None
        valid_mask = None

    values = comm.bcast(values, root=0)
    valid_mask = comm.bcast(valid_mask, root=0)
    return values, valid_mask

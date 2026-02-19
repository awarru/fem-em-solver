"""Mesh test for two-torus Helmholtz prototype geometry."""

import numpy as np
from mpi4py import MPI

from fem_em_solver.io.mesh import MeshGenerator


def test_two_torus_mesh_generates_with_two_wire_volumes():
    """Generate two-torus mesh and verify wire/domain tags exist."""
    mesh, cell_tags, _ = MeshGenerator.two_torus_domain(
        separation=0.05,
        major_radius=0.02,
        minor_radius=0.005,
        resolution=0.01,
        comm=MPI.COMM_WORLD,
    )

    n_cells = mesh.topology.index_map(3).size_global
    assert n_cells > 0, "Mesh must contain cells"

    unique_tags = np.unique(cell_tags.values)
    assert 1 in unique_tags, "Missing wire_1 volume tag"
    assert 2 in unique_tags, "Missing wire_2 volume tag"
    assert 3 in unique_tags, "Missing domain volume tag"

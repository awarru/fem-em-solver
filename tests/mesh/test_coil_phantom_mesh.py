"""Mesh test for coarse coil + phantom + air geometry."""

import numpy as np
from mpi4py import MPI

from fem_em_solver.io.mesh import MeshGenerator


def test_coil_phantom_mesh_generates_required_tags():
    """Generate coil+phantom mesh and verify required region tags exist."""
    mesh, cell_tags, _ = MeshGenerator.coil_phantom_domain(
        coil_major_radius=0.08,
        coil_minor_radius=0.01,
        coil_separation=0.08,
        phantom_radius=0.04,
        phantom_height=0.10,
        air_padding=0.04,
        resolution=0.015,
        comm=MPI.COMM_WORLD,
    )

    n_cells = mesh.topology.index_map(3).size_global
    assert n_cells > 0, "Mesh must contain cells"
    assert n_cells < 50000, f"Mesh must stay under 50k cells, got {n_cells}"

    unique_tags = np.unique(cell_tags.values)
    assert 1 in unique_tags, "Missing coil_1 volume tag"
    assert 2 in unique_tags, "Missing coil_2 volume tag"
    assert 3 in unique_tags, "Missing phantom volume tag"
    assert 4 in unique_tags, "Missing air volume tag"

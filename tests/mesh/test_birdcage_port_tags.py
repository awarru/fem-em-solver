"""Chunk E2: coarse birdcage-like geometry fixture with explicit port tags."""

from __future__ import annotations

import numpy as np
from mpi4py import MPI

from fem_em_solver.io.mesh import MeshGenerator
from fem_em_solver.io.mesh_qa import print_cell_tag_summary


BIRDCAGE_CORE_TAGS = {
    1: "conductor",
    2: "air",
    3: "phantom",
}


def test_birdcage_like_mesh_has_core_and_port_tags():
    """Ensure coarse birdcage fixture provides all required core + port tags."""
    comm = MPI.COMM_WORLD

    n_legs = 4
    mesh, cell_tags, _ = MeshGenerator.birdcage_port_domain(
        n_legs=n_legs,
        ring_radius=0.07,
        leg_radius=0.006,
        leg_height=0.14,
        ring_minor_radius=0.004,
        phantom_radius=0.03,
        phantom_height=0.08,
        port_box_size=(0.010, 0.008, 0.010),
        air_padding=0.03,
        resolution=0.015,
        comm=comm,
    )

    n_cells = mesh.topology.index_map(3).size_global
    assert n_cells > 0, "Mesh must contain cells"
    assert n_cells < 50000, f"Mesh must stay under 50k cells, got {n_cells}"

    unique_tags = set(np.unique(cell_tags.values).tolist())

    expected_port_tags = [100 + i for i in range(1, n_legs + 1)]

    missing_core = [name for tag, name in BIRDCAGE_CORE_TAGS.items() if tag not in unique_tags]
    assert not missing_core, f"Missing core birdcage tags: {', '.join(missing_core)}"

    missing_ports = [f"P{i}" for i, tag in enumerate(expected_port_tags, start=1) if tag not in unique_tags]
    assert not missing_ports, f"Missing port tags for {', '.join(missing_ports)}"

    print_cell_tag_summary(
        cell_tags,
        tag_names={
            **BIRDCAGE_CORE_TAGS,
            **{tag: f"port_P{i}" for i, tag in enumerate(expected_port_tags, start=1)},
        },
        comm=comm,
        prefix="[birdcage-mesh] ",
    )

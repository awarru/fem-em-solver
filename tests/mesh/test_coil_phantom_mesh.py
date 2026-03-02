"""Mesh tests for coarse coil + phantom + air geometry presets."""

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver.io.mesh import MeshGenerator
from tests.mesh.helpers import REQUIRED_COIL_PHANTOM_TAGS, assert_required_tags_nonempty, compute_tag_cell_centroid


def _generate_default_mesh(**kwargs):
    return MeshGenerator.coil_phantom_domain(
        coil_major_radius=0.08,
        coil_minor_radius=0.01,
        coil_separation=0.08,
        phantom_radius=0.04,
        phantom_height=0.10,
        air_padding=0.04,
        resolution=0.015,
        comm=MPI.COMM_WORLD,
        **kwargs,
    )


def test_coil_phantom_mesh_generates_required_tags_centered_preset():
    """Centered preset should generate required tags, including non-empty phantom."""
    mesh, cell_tags, _ = _generate_default_mesh(phantom_placement_preset="centered")

    n_cells = mesh.topology.index_map(3).size_global
    assert n_cells > 0, "Mesh must contain cells"
    assert n_cells < 50000, f"Mesh must stay under 50k cells, got {n_cells}"

    assert_required_tags_nonempty(cell_tags, REQUIRED_COIL_PHANTOM_TAGS, comm=MPI.COMM_WORLD)

    phantom_centroid = compute_tag_cell_centroid(mesh, cell_tags, tag=3, comm=MPI.COMM_WORLD)
    assert abs(phantom_centroid[0]) < 1.0e-2
    assert abs(phantom_centroid[1]) < 1.0e-2


def test_coil_phantom_mesh_off_center_preset_moves_phantom_without_overlap():
    """Off-center preset should keep phantom non-empty and displaced in +x direction."""
    mesh, cell_tags, _ = _generate_default_mesh(phantom_placement_preset="off_center")

    assert_required_tags_nonempty(cell_tags, REQUIRED_COIL_PHANTOM_TAGS, comm=MPI.COMM_WORLD)

    phantom_centroid = compute_tag_cell_centroid(mesh, cell_tags, tag=3, comm=MPI.COMM_WORLD)
    assert phantom_centroid[0] > 5.0e-3, (
        "Off-center preset should produce positive x displacement for phantom centroid"
    )


def test_coil_phantom_mesh_rejects_overlapping_off_center_placement():
    """Explicit offsets that overlap the coil envelope should be rejected."""
    with pytest.raises(ValueError, match="Phantom overlaps coil conductor envelope"):
        _generate_default_mesh(
            phantom_placement_preset="off_center",
            phantom_offset_xy=(0.04, 0.0),
        )

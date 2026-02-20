"""Mesh QA checks for required and distinct coil+phantom tags."""

from mpi4py import MPI

from fem_em_solver.io.mesh import MeshGenerator
from fem_em_solver.io.mesh_qa import print_cell_tag_summary

from tests.mesh.helpers import (
    REQUIRED_COIL_PHANTOM_TAGS,
    assert_required_tags_nonempty,
    assert_tags_distinct_by_centroid,
)


def test_coil_phantom_mesh_tag_integrity():
    """Ensure coarse coil+phantom mesh has non-empty, distinct required tags."""
    comm = MPI.COMM_WORLD

    mesh, cell_tags, _ = MeshGenerator.coil_phantom_domain(
        coil_major_radius=0.08,
        coil_minor_radius=0.01,
        coil_separation=0.08,
        phantom_radius=0.04,
        phantom_height=0.10,
        air_padding=0.04,
        resolution=0.015,
        comm=comm,
    )

    print_cell_tag_summary(cell_tags, tag_names=REQUIRED_COIL_PHANTOM_TAGS, comm=comm, prefix="[mesh-qa] ")

    assert_required_tags_nonempty(cell_tags, REQUIRED_COIL_PHANTOM_TAGS, comm=comm)

    # Distinctness checks requested by roadmap chunk B2.
    assert_tags_distinct_by_centroid(mesh, cell_tags, tag_a=1, tag_b=3, comm=comm)
    assert_tags_distinct_by_centroid(mesh, cell_tags, tag_a=2, tag_b=3, comm=comm)

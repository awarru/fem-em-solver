"""Chunk E2/B2: birdcage-like geometry fixture with explicit robust port tags."""

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
    """Ensure parametric birdcage fixture provides all required core + port tags."""
    comm = MPI.COMM_WORLD

    leg_count = 4
    diagnostics = MeshGenerator.birdcage_port_layout_diagnostics(
        leg_count=leg_count,
        ring_radius=0.07,
        leg_width=0.012,
        ring_minor_radius=0.004,
        phantom_radius=0.03,
        port_box_size=(0.010, 0.008, 0.010),
    )

    if comm.rank == 0:
        print(
            "[birdcage-port] area/separation diagnostics: "
            f"port_face_area={diagnostics['port_face_area_m2']:.6e} m^2 "
            f"(required>={diagnostics['min_port_face_area_m2']:.6e}), "
            f"min_center_separation={diagnostics['min_port_center_separation_m']:.6e} m "
            f"(required>={diagnostics['required_port_center_separation_m']:.6e}), "
            f"conductor_clearance={diagnostics['conductor_radial_clearance_m']:.6e} m, "
            f"phantom_clearance={diagnostics['phantom_radial_clearance_m']:.6e} m"
        )

    mesh, cell_tags, _ = MeshGenerator.birdcage_port_domain(
        leg_count=leg_count,
        ring_radius=0.07,
        leg_width=0.012,
        leg_spacing=0.11,
        coil_length=0.14,
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

    expected_port_tags = [100 + i for i in range(1, leg_count + 1)]

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


def test_birdcage_port_layout_rejects_too_small_or_overlapping_port_regions():
    """Reject implausible port geometry before meshing."""
    # Too-small face area check.
    try:
        MeshGenerator.birdcage_port_layout_diagnostics(
            leg_count=4,
            ring_radius=0.07,
            leg_width=0.012,
            ring_minor_radius=0.004,
            phantom_radius=0.03,
            port_box_size=(0.003, 0.008, 0.003),
            min_port_face_area=2.5e-5,
        )
        raise AssertionError("Expected ValueError for too-small port face area")
    except ValueError as exc:
        assert "Port face area too small" in str(exc)

    # Overlap/separation check via unrealistically dense port packing.
    try:
        MeshGenerator.birdcage_port_layout_diagnostics(
            leg_count=32,
            ring_radius=0.03,
            leg_width=0.010,
            ring_minor_radius=0.004,
            phantom_radius=0.02,
            port_box_size=(0.010, 0.010, 0.010),
            min_port_center_separation=0.015,
            port_clearance=0.0,
        )
        raise AssertionError("Expected ValueError for insufficient port center separation")
    except ValueError as exc:
        assert "Port center separation too small" in str(exc)

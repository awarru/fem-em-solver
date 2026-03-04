"""Unit tests for coil+phantom geometry sanity report helpers."""

import numpy as np
from mpi4py import MPI

from fem_em_solver.io.mesh import MeshGenerator


class _DummyCellTags:
    """Minimal test double exposing a ``values`` array like dolfinx MeshTags."""

    def __init__(self, values):
        self.values = np.asarray(values, dtype=np.int32)


def test_coil_phantom_geometry_sanity_report_includes_expected_sections():
    comm = MPI.COMM_WORLD
    cell_tags = _DummyCellTags([1] * 6 + [2] * 6 + [3] * 30 + [4] * 958)

    report = MeshGenerator.coil_phantom_geometry_sanity_report(
        cell_tags=cell_tags,
        coil_major_radius=0.08,
        coil_minor_radius=0.01,
        coil_separation=0.08,
        phantom_radius=0.04,
        phantom_height=0.10,
        air_padding=0.04,
        comm=comm,
    )

    assert report["ok"] is True
    assert report["warnings"] == []

    assert set(report["required_tag_counts"].keys()) == {"coil_1", "coil_2", "phantom", "air"}
    assert set(report["expected_volume_ratios"].keys()) == {"coil_1", "coil_2", "phantom", "air"}
    assert set(report["observed_cell_ratios"].keys()) == {"coil_1", "coil_2", "phantom", "air"}


def test_coil_phantom_geometry_sanity_report_warns_for_missing_required_tag():
    comm = MPI.COMM_WORLD
    cell_tags = _DummyCellTags([1] * 10 + [2] * 10 + [4] * 980)

    report = MeshGenerator.coil_phantom_geometry_sanity_report(
        cell_tags=cell_tags,
        coil_major_radius=0.08,
        coil_minor_radius=0.01,
        coil_separation=0.08,
        phantom_radius=0.04,
        phantom_height=0.10,
        air_padding=0.04,
        comm=comm,
    )

    assert report["ok"] is False
    joined = "\n".join(report["warnings"])
    assert "missing required tags" in joined
    assert "phantom" in joined

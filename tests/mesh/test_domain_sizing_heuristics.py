"""Domain sizing heuristic checks for coil+phantom geometry."""

import pytest

from fem_em_solver.io.mesh import MeshGenerator


def test_coil_phantom_domain_sizing_defaults_are_not_undersized():
    diagnostics = MeshGenerator.coil_phantom_domain_sizing_diagnostics(
        coil_major_radius=0.08,
        coil_minor_radius=0.01,
        coil_separation=0.08,
        phantom_radius=0.04,
        phantom_height=0.10,
        air_padding=0.04,
        phantom_offset_xy=(0.0, 0.0),
    )

    assert diagnostics["is_domain_undersized"] is False
    assert diagnostics["effective_air_padding_m"] == pytest.approx(0.04)
    assert diagnostics["provided_air_padding_m"] >= diagnostics["recommended_min_air_padding_m"]


def test_coil_phantom_domain_sizing_detects_small_padding_and_recommends_floor():
    diagnostics = MeshGenerator.coil_phantom_domain_sizing_diagnostics(
        coil_major_radius=0.08,
        coil_minor_radius=0.01,
        coil_separation=0.08,
        phantom_radius=0.04,
        phantom_height=0.10,
        air_padding=0.005,
        phantom_offset_xy=(0.0, 0.0),
    )

    assert diagnostics["is_domain_undersized"] is True
    assert diagnostics["effective_air_padding_m"] == pytest.approx(
        diagnostics["recommended_min_air_padding_m"]
    )
    assert diagnostics["effective_air_padding_m"] > diagnostics["provided_air_padding_m"]


def test_coil_phantom_domain_sizing_accounts_for_off_center_phantom_extent():
    centered = MeshGenerator.coil_phantom_domain_sizing_diagnostics(
        coil_major_radius=0.08,
        coil_minor_radius=0.01,
        coil_separation=0.08,
        phantom_radius=0.04,
        phantom_height=0.10,
        air_padding=0.04,
        phantom_offset_xy=(0.0, 0.0),
    )
    shifted = MeshGenerator.coil_phantom_domain_sizing_diagnostics(
        coil_major_radius=0.08,
        coil_minor_radius=0.01,
        coil_separation=0.08,
        phantom_radius=0.04,
        phantom_height=0.10,
        air_padding=0.04,
        phantom_offset_xy=(0.03, 0.0),
    )

    assert shifted["radial_extent_without_padding_m"] > centered["radial_extent_without_padding_m"]
    assert shifted["recommended_domain_half_width_m"] > centered["recommended_domain_half_width_m"]


def test_coil_phantom_domain_sizing_rejects_negative_air_padding():
    with pytest.raises(ValueError, match="air_padding must be >= 0"):
        MeshGenerator.coil_phantom_domain_sizing_diagnostics(
            coil_major_radius=0.08,
            coil_minor_radius=0.01,
            coil_separation=0.08,
            phantom_radius=0.04,
            phantom_height=0.10,
            air_padding=-1.0e-3,
            phantom_offset_xy=(0.0, 0.0),
        )

"""Unit tests for coil+phantom region mesh resolution policy helpers."""

import pytest

from fem_em_solver.io.mesh import MeshGenerator


def test_region_resolution_policy_defaults_to_global_resolution():
    policy = MeshGenerator.coil_phantom_region_resolution_policy(resolution=0.015)

    assert policy["coil_resolution_m"] == pytest.approx(0.015)
    assert policy["phantom_resolution_m"] == pytest.approx(0.015)
    assert policy["air_resolution_m"] == pytest.approx(0.015)
    assert policy["min_resolution_m"] == pytest.approx(0.015)
    assert policy["max_resolution_m"] == pytest.approx(0.015)


def test_region_resolution_policy_prefers_explicit_kwargs_over_mapping():
    policy = MeshGenerator.coil_phantom_region_resolution_policy(
        resolution=0.02,
        coil_resolution=0.01,
        region_resolutions={"coil": 0.03, "phantom": 0.015, "air": 0.025},
    )

    assert policy["coil_resolution_m"] == pytest.approx(0.01)
    assert policy["phantom_resolution_m"] == pytest.approx(0.015)
    assert policy["air_resolution_m"] == pytest.approx(0.025)


def test_region_resolution_policy_rejects_unknown_mapping_keys():
    with pytest.raises(ValueError, match="unsupported keys"):
        MeshGenerator.coil_phantom_region_resolution_policy(
            resolution=0.015,
            region_resolutions={"coil": 0.01, "vacuum": 0.02},
        )

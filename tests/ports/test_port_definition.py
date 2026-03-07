"""Tests for lumped port data model + tag-validation contract (chunk E1)."""

from __future__ import annotations

import pytest

from fem_em_solver.ports import (
    DEFAULT_REFERENCE_IMPEDANCE_OHM,
    PORT_BOUNDARY_DIMENSION_FALLBACK,
    PORT_BOUNDARY_DIMENSION_PRIMARY,
    PORT_CALIBRATION_CHECKLIST_DOC,
    PORT_NEGATIVE_SUFFIX,
    PORT_POSITIVE_SUFFIX,
    PORT_TAG_NAME_TEMPLATE,
    PortDefinition,
    required_port_tags,
    run_port_calibration_checks,
    validate_required_port_tags_exist,
)


def test_port_definition_defaults_and_tag_name_contract():
    port = PortDefinition(
        port_id="P1",
        positive_tag=101,
        negative_tag=102,
        orientation="leg_1_to_leg_2",
    )
    port.validate()

    assert port.z0_ohm == DEFAULT_REFERENCE_IMPEDANCE_OHM
    assert PORT_BOUNDARY_DIMENSION_PRIMARY == "facet"
    assert PORT_BOUNDARY_DIMENSION_FALLBACK == "edge"
    assert port.positive_tag_name == PORT_TAG_NAME_TEMPLATE.format(
        port_id="P1", suffix=PORT_POSITIVE_SUFFIX
    )
    assert port.negative_tag_name == PORT_TAG_NAME_TEMPLATE.format(
        port_id="P1", suffix=PORT_NEGATIVE_SUFFIX
    )


def test_required_port_tags_collection_and_validation_success():
    ports = [
        PortDefinition(port_id="P1", positive_tag=11, negative_tag=12, orientation="cw"),
        PortDefinition(port_id="P2", positive_tag=13, negative_tag=14, orientation="cw"),
    ]

    assert required_port_tags(ports) == {11, 12, 13, 14}
    validate_required_port_tags_exist(ports, available_tags=[0, 3, 11, 12, 13, 14])


def test_validate_required_port_tags_exist_raises_for_missing_tags():
    ports = [PortDefinition(port_id="P3", positive_tag=31, negative_tag=32, orientation="ccw")]

    with pytest.raises(ValueError, match=r"missing required port tags: \[32\]"):
        validate_required_port_tags_exist(ports, available_tags=[31])


def test_port_definition_rejects_invalid_values():
    with pytest.raises(ValueError, match="non-empty"):
        PortDefinition(port_id="", positive_tag=1, negative_tag=2, orientation="x").validate()

    with pytest.raises(ValueError, match="distinct"):
        PortDefinition(port_id="P", positive_tag=1, negative_tag=1, orientation="x").validate()

    with pytest.raises(ValueError, match="positive"):
        PortDefinition(
            port_id="P", positive_tag=1, negative_tag=2, orientation="x", z0_ohm=0.0
        ).validate()


def test_run_port_calibration_checks_accepts_consistent_order_orientation_and_area():
    ports = [
        PortDefinition(port_id="P1", positive_tag=11, negative_tag=12, orientation="leg_1_to_leg_2"),
        PortDefinition(port_id="P2", positive_tag=13, negative_tag=14, orientation="leg_2_to_leg_3"),
    ]

    checks = run_port_calibration_checks(
        ports,
        expected_port_order=["P1", "P2"],
        port_face_areas={"P1": 1.00e-4, "P2": 1.10e-4},
        max_area_ratio=1.25,
    )

    assert checks.port_order == ("P1", "P2")
    assert checks.orientations == ("leg_1_to_leg_2", "leg_2_to_leg_3")
    assert checks.face_areas == {"P1": 1.00e-4, "P2": 1.10e-4}
    assert checks.max_face_area_ratio == pytest.approx(1.1)


def test_run_port_calibration_checks_rejects_order_mismatch_with_checklist_reference():
    ports = [
        PortDefinition(port_id="P1", positive_tag=11, negative_tag=12, orientation="leg_1_to_leg_2"),
        PortDefinition(port_id="P2", positive_tag=13, negative_tag=14, orientation="leg_2_to_leg_3"),
    ]

    with pytest.raises(ValueError, match="port ordering mismatch") as exc_info:
        run_port_calibration_checks(
            ports,
            expected_port_order=["P2", "P1"],
            port_face_areas={"P1": 1.0e-4, "P2": 1.0e-4},
        )

    assert PORT_CALIBRATION_CHECKLIST_DOC in str(exc_info.value)


def test_run_port_calibration_checks_rejects_missing_orientation_metadata():
    ports = [
        PortDefinition(port_id="P1", positive_tag=11, negative_tag=12, orientation="leg_1_to_leg_2"),
        PortDefinition(port_id="P2", positive_tag=13, negative_tag=14, orientation="   "),
    ]

    with pytest.raises(ValueError, match="missing orientation metadata") as exc_info:
        run_port_calibration_checks(
            ports,
            expected_port_order=["P1", "P2"],
            port_face_areas={"P1": 1.0e-4, "P2": 1.0e-4},
        )

    assert "P2" in str(exc_info.value)
    assert PORT_CALIBRATION_CHECKLIST_DOC in str(exc_info.value)


def test_run_port_calibration_checks_rejects_inconsistent_face_area_ratio():
    ports = [
        PortDefinition(port_id="P1", positive_tag=11, negative_tag=12, orientation="leg_1_to_leg_2"),
        PortDefinition(port_id="P2", positive_tag=13, negative_tag=14, orientation="leg_2_to_leg_3"),
    ]

    with pytest.raises(ValueError, match="area ratio exceeds threshold") as exc_info:
        run_port_calibration_checks(
            ports,
            expected_port_order=["P1", "P2"],
            port_face_areas={"P1": 1.0e-4, "P2": 3.0e-4},
            max_area_ratio=1.5,
        )

    assert PORT_CALIBRATION_CHECKLIST_DOC in str(exc_info.value)

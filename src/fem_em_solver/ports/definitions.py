"""Lumped-port definitions and calibration checklist validation helpers."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence


DEFAULT_REFERENCE_IMPEDANCE_OHM = 50.0
PORT_POSITIVE_SUFFIX = "positive"
PORT_NEGATIVE_SUFFIX = "negative"
PORT_SUFFIXES = (PORT_POSITIVE_SUFFIX, PORT_NEGATIVE_SUFFIX)

# Mesh tagging contract for lumped-port boundaries.
#
# Port terminal tags are expected on mesh facets (preferred) or edges (fallback),
# using explicit positive/negative polarity tags per port:
#   - positive terminal name: "port_<port_id>_positive"
#   - negative terminal name: "port_<port_id>_negative"
PORT_BOUNDARY_DIMENSION_PRIMARY = "facet"
PORT_BOUNDARY_DIMENSION_FALLBACK = "edge"
PORT_TAG_NAME_TEMPLATE = "port_{port_id}_{suffix}"

# Human calibration checklist document that these checks operationalize.
PORT_CALIBRATION_CHECKLIST_DOC = "docs/ports/human_port_calibration_checklist.md"


@dataclass(frozen=True)
class PortDefinition:
    """Explicit data model for one lumped port.

    Parameters
    ----------
    port_id
        Stable logical port identifier (e.g. ``P1``).
    positive_tag
        Integer mesh tag for the positive terminal boundary region.
    negative_tag
        Integer mesh tag for the negative terminal boundary region.
    orientation
        Free-form feed direction/orientation metadata (e.g. ``leg_1_to_leg_2``).
    z0_ohm
        Reference impedance in Ohms. Defaults to 50 Ω.
    """

    port_id: str
    positive_tag: int
    negative_tag: int
    orientation: str
    z0_ohm: float = DEFAULT_REFERENCE_IMPEDANCE_OHM

    def validate(self) -> None:
        """Raise ``ValueError`` when the port definition is invalid."""
        if not self.port_id or not self.port_id.strip():
            raise ValueError("port_id must be a non-empty string")
        if self.positive_tag < 0 or self.negative_tag < 0:
            raise ValueError("port tags must be non-negative integers")
        if self.positive_tag == self.negative_tag:
            raise ValueError("positive_tag and negative_tag must be distinct")
        if self.z0_ohm <= 0.0:
            raise ValueError("z0_ohm must be positive")

    @property
    def positive_tag_name(self) -> str:
        return PORT_TAG_NAME_TEMPLATE.format(port_id=self.port_id, suffix=PORT_POSITIVE_SUFFIX)

    @property
    def negative_tag_name(self) -> str:
        return PORT_TAG_NAME_TEMPLATE.format(port_id=self.port_id, suffix=PORT_NEGATIVE_SUFFIX)


def required_port_tags(ports: Sequence[PortDefinition]) -> set[int]:
    """Return the set of all required positive/negative terminal tags."""
    tags: set[int] = set()
    for port in ports:
        port.validate()
        tags.add(int(port.positive_tag))
        tags.add(int(port.negative_tag))
    return tags


def validate_required_port_tags_exist(
    ports: Sequence[PortDefinition],
    available_tags: Iterable[int],
) -> None:
    """Validate that all required port tags exist in ``available_tags``.

    Parameters
    ----------
    ports
        Port definitions to validate.
    available_tags
        Iterable of integer tag values present in the mesh tag data.
    """
    expected = required_port_tags(ports)
    observed = {int(tag) for tag in available_tags}
    missing = sorted(expected.difference(observed))
    if missing:
        raise ValueError(f"missing required port tags: {missing}")


def validate_port_ordering(
    ports: Sequence[PortDefinition],
    *,
    expected_port_order: Sequence[str],
) -> tuple[str, ...]:
    """Assert deterministic port ordering against an expected sequence.

    This maps to checklist section: "Z0 and normalization choices" where
    port ordering must stay stable and documented.
    """
    if not expected_port_order:
        raise ValueError("expected_port_order must be non-empty")

    observed = tuple(port.port_id for port in ports)
    expected = tuple(expected_port_order)
    if observed != expected:
        raise ValueError(
            "port ordering mismatch: "
            f"observed={list(observed)} expected={list(expected)}; "
            f"see {PORT_CALIBRATION_CHECKLIST_DOC}"
        )
    return observed


def validate_port_orientation_metadata(ports: Sequence[PortDefinition]) -> tuple[str, ...]:
    """Assert that every port carries non-empty orientation metadata.

    This maps to checklist section: "Port placement realism".
    Returns the normalized orientation tuple in port order.
    """
    orientations: list[str] = []
    missing: list[str] = []
    for port in ports:
        port.validate()
        normalized = port.orientation.strip() if isinstance(port.orientation, str) else ""
        if not normalized:
            missing.append(port.port_id)
        orientations.append(normalized)

    if missing:
        raise ValueError(
            "missing orientation metadata for ports: "
            f"{missing}; see {PORT_CALIBRATION_CHECKLIST_DOC}"
        )
    return tuple(orientations)


def validate_port_face_area_consistency(
    ports: Sequence[PortDefinition],
    port_face_areas: Mapping[str, float],
    *,
    max_area_ratio: float = 1.25,
) -> dict[str, float]:
    """Assert that per-port face areas are finite and roughly comparable.

    Parameters
    ----------
    ports
        Port definitions in solver order.
    port_face_areas
        Mapping ``port_id -> area`` from geometry/mesh QA routines.
    max_area_ratio
        Maximum allowed ratio ``max(area)/min(area)`` for consistency checks.

    Returns
    -------
    dict[str, float]
        Normalized area mapping in the same logical order as ``ports``.
    """
    if max_area_ratio < 1.0:
        raise ValueError("max_area_ratio must be >= 1.0")

    normalized: dict[str, float] = {}
    for port in ports:
        if port.port_id not in port_face_areas:
            raise ValueError(
                f"missing face area for port '{port.port_id}'; "
                f"see {PORT_CALIBRATION_CHECKLIST_DOC}"
            )

        area = float(port_face_areas[port.port_id])
        if area <= 0.0:
            raise ValueError(
                f"port '{port.port_id}' face area must be positive, got {area:.6e}; "
                f"see {PORT_CALIBRATION_CHECKLIST_DOC}"
            )
        normalized[port.port_id] = area

    values = tuple(normalized.values())
    area_ratio = max(values) / min(values)
    if area_ratio > max_area_ratio:
        raise ValueError(
            "port face area ratio exceeds threshold: "
            f"ratio={area_ratio:.6f}, max_area_ratio={max_area_ratio:.6f}; "
            f"see {PORT_CALIBRATION_CHECKLIST_DOC}"
        )

    return normalized


@dataclass(frozen=True)
class PortCalibrationChecks:
    """Summary of executable checks derived from the human checklist."""

    port_order: tuple[str, ...]
    orientations: tuple[str, ...]
    face_areas: dict[str, float]
    max_face_area_ratio: float


def run_port_calibration_checks(
    ports: Sequence[PortDefinition],
    *,
    expected_port_order: Sequence[str],
    port_face_areas: Mapping[str, float],
    max_area_ratio: float = 1.25,
) -> PortCalibrationChecks:
    """Run checklist-derived executable assertions for lumped ports."""
    port_order = validate_port_ordering(ports, expected_port_order=expected_port_order)
    orientations = validate_port_orientation_metadata(ports)
    face_areas = validate_port_face_area_consistency(
        ports,
        port_face_areas,
        max_area_ratio=max_area_ratio,
    )

    ratio = max(face_areas.values()) / min(face_areas.values())
    return PortCalibrationChecks(
        port_order=port_order,
        orientations=orientations,
        face_areas=face_areas,
        max_face_area_ratio=ratio,
    )

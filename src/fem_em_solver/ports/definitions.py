"""Lumped-port definitions and mesh-tag validation helpers (chunk E1)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence


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

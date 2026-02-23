"""Lumped-port data model and tag-contract helpers."""

from .definitions import (
    DEFAULT_REFERENCE_IMPEDANCE_OHM,
    PORT_BOUNDARY_DIMENSION_FALLBACK,
    PORT_BOUNDARY_DIMENSION_PRIMARY,
    PORT_NEGATIVE_SUFFIX,
    PORT_POSITIVE_SUFFIX,
    PORT_SUFFIXES,
    PORT_TAG_NAME_TEMPLATE,
    PortDefinition,
    required_port_tags,
    validate_required_port_tags_exist,
)
from .excitation import (
    PortVoltageCurrentEstimate,
    SinglePortExcitationResult,
    run_single_port_excitation_case,
)
from .sparameters import SParameterSweepResult, run_n_port_sparameter_sweep

__all__ = [
    "DEFAULT_REFERENCE_IMPEDANCE_OHM",
    "PORT_BOUNDARY_DIMENSION_FALLBACK",
    "PORT_BOUNDARY_DIMENSION_PRIMARY",
    "PORT_NEGATIVE_SUFFIX",
    "PORT_POSITIVE_SUFFIX",
    "PORT_SUFFIXES",
    "PORT_TAG_NAME_TEMPLATE",
    "PortDefinition",
    "required_port_tags",
    "validate_required_port_tags_exist",
    "PortVoltageCurrentEstimate",
    "SinglePortExcitationResult",
    "run_single_port_excitation_case",
    "SParameterSweepResult",
    "run_n_port_sparameter_sweep",
]

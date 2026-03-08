"""Lumped-port data model and tag-contract helpers."""

from .definitions import (
    DEFAULT_REFERENCE_IMPEDANCE_OHM,
    PORT_BOUNDARY_DIMENSION_FALLBACK,
    PORT_BOUNDARY_DIMENSION_PRIMARY,
    PORT_CALIBRATION_CHECKLIST_DOC,
    PORT_NEGATIVE_SUFFIX,
    PORT_POSITIVE_SUFFIX,
    PORT_SUFFIXES,
    PORT_TAG_NAME_TEMPLATE,
    PortCalibrationChecks,
    PortDefinition,
    required_port_tags,
    run_port_calibration_checks,
    validate_port_face_area_consistency,
    validate_port_ordering,
    validate_port_orientation_metadata,
    validate_required_port_tags_exist,
)
from .excitation import (
    PortSolveContext,
    PortVoltageCurrentEstimate,
    SinglePortExcitationResult,
    run_single_port_excitation_case,
)
from .sparameters import (
    SMatrixSanityReport,
    SParameterSweepResult,
    run_n_port_sparameter_sweep,
    summarize_sparameter_sanity,
)
from .touchstone import export_touchstone, load_touchstone

__all__ = [
    "DEFAULT_REFERENCE_IMPEDANCE_OHM",
    "PORT_BOUNDARY_DIMENSION_FALLBACK",
    "PORT_BOUNDARY_DIMENSION_PRIMARY",
    "PORT_CALIBRATION_CHECKLIST_DOC",
    "PORT_NEGATIVE_SUFFIX",
    "PORT_POSITIVE_SUFFIX",
    "PORT_SUFFIXES",
    "PORT_TAG_NAME_TEMPLATE",
    "PortCalibrationChecks",
    "PortDefinition",
    "required_port_tags",
    "run_port_calibration_checks",
    "validate_port_face_area_consistency",
    "validate_port_ordering",
    "validate_port_orientation_metadata",
    "validate_required_port_tags_exist",
    "PortSolveContext",
    "PortVoltageCurrentEstimate",
    "SinglePortExcitationResult",
    "run_single_port_excitation_case",
    "SMatrixSanityReport",
    "SParameterSweepResult",
    "run_n_port_sparameter_sweep",
    "summarize_sparameter_sanity",
    "export_touchstone",
    "load_touchstone",
]

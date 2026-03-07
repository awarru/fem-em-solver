"""FEM Electromagnetics Solver for MRI Coil Simulation.

This package provides finite element solvers for electromagnetic simulations
using FEniCSX/DolfinX, with a focus on MRI coil design and analysis.
"""

__version__ = "0.1.0"
__author__ = "Awarru"

# Core imports
from .core import (
    HomogeneousMaterial,
    LinearSolveDiagnostics,
    MagnetostaticProblem,
    MagnetostaticSolver,
    classify_residual_trend,
    TimeHarmonicBoundaryCondition,
    TimeHarmonicFields,
    TimeHarmonicProblem,
    TimeHarmonicSolver,
    build_material_fields,
    normalize_boundary_condition,
)
from .materials import GelledSalinePhantomMaterial
from .ports import (
    PortDefinition,
    PortVoltageCurrentEstimate,
    SParameterSweepResult,
    SinglePortExcitationResult,
    export_touchstone,
    load_touchstone,
    run_n_port_sparameter_sweep,
    run_single_port_excitation_case,
    validate_required_port_tags_exist,
)

__all__ = [
    "LinearSolveDiagnostics",
    "MagnetostaticProblem",
    "MagnetostaticSolver",
    "classify_residual_trend",
    "HomogeneousMaterial",
    "TimeHarmonicBoundaryCondition",
    "TimeHarmonicFields",
    "TimeHarmonicProblem",
    "TimeHarmonicSolver",
    "build_material_fields",
    "normalize_boundary_condition",
    "GelledSalinePhantomMaterial",
    "PortDefinition",
    "PortVoltageCurrentEstimate",
    "SParameterSweepResult",
    "SinglePortExcitationResult",
    "export_touchstone",
    "load_touchstone",
    "run_n_port_sparameter_sweep",
    "run_single_port_excitation_case",
    "validate_required_port_tags_exist",
]

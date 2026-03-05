"""Core FEM formulations and solvers."""

from .solvers import MagnetostaticProblem, MagnetostaticSolver
from .time_harmonic import (
    HomogeneousMaterial,
    TimeHarmonicBoundaryCondition,
    TimeHarmonicFields,
    TimeHarmonicProblem,
    TimeHarmonicSolver,
    build_material_fields,
    normalize_boundary_condition,
)

__all__ = [
    "MagnetostaticProblem",
    "MagnetostaticSolver",
    "HomogeneousMaterial",
    "TimeHarmonicBoundaryCondition",
    "TimeHarmonicFields",
    "TimeHarmonicProblem",
    "TimeHarmonicSolver",
    "build_material_fields",
    "normalize_boundary_condition",
]

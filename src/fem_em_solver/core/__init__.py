"""Core FEM formulations and solvers."""

from .solvers import MagnetostaticProblem, MagnetostaticSolver
from .time_harmonic import (
    HomogeneousMaterial,
    TimeHarmonicFields,
    TimeHarmonicProblem,
    TimeHarmonicSolver,
    build_material_fields,
)

__all__ = [
    "MagnetostaticProblem",
    "MagnetostaticSolver",
    "HomogeneousMaterial",
    "TimeHarmonicFields",
    "TimeHarmonicProblem",
    "TimeHarmonicSolver",
    "build_material_fields",
]

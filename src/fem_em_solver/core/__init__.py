"""Core FEM formulations and solvers."""

from .solvers import (
    LinearSolveDiagnostics,
    MagnetostaticProblem,
    MagnetostaticSolver,
    classify_residual_trend,
)
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
]

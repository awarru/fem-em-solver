"""FEM Electromagnetics Solver for MRI Coil Simulation.

This package provides finite element solvers for electromagnetic simulations
using FEniCSX/DolfinX, with a focus on MRI coil design and analysis.
"""

__version__ = "0.1.0"
__author__ = "Awarru"

# Core imports
from .core import (
    HomogeneousMaterial,
    MagnetostaticProblem,
    MagnetostaticSolver,
    TimeHarmonicFields,
    TimeHarmonicProblem,
    TimeHarmonicSolver,
)

__all__ = [
    "MagnetostaticProblem",
    "MagnetostaticSolver",
    "HomogeneousMaterial",
    "TimeHarmonicFields",
    "TimeHarmonicProblem",
    "TimeHarmonicSolver",
]

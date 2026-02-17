"""FEM Electromagnetics Solver for MRI Coil Simulation.

This package provides finite element solvers for electromagnetic simulations
using FEniCSX/DolfinX, with a focus on MRI coil design and analysis.
"""

__version__ = "0.1.0"
__author__ = "Awarru"

# Core imports
from .core.solvers import MagnetostaticSolver, TimeHarmonicSolver
from .core.function_spaces import create_function_spaces

# Coil imports
from .coils.loop import CircularLoop

# Material imports
from .materials.phantoms import GelledSalinePhantom

__all__ = [
    "MagnetostaticSolver",
    "TimeHarmonicSolver", 
    "create_function_spaces",
    "CircularLoop",
    "GelledSalinePhantom",
]

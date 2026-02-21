"""Minimal time-harmonic solver scaffold for frequency-domain EM workflows."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional, Sequence

import numpy as np
import dolfinx
from dolfinx import fem

from ..utils.constants import EPSILON_0, MU_0
from .solvers import MagnetostaticProblem, MagnetostaticSolver


@dataclass(frozen=True)
class HomogeneousMaterial:
    """Homogeneous linear material properties for MVP frequency-domain solves."""

    sigma: float
    epsilon_r: float
    mu_r: float = 1.0


@dataclass(frozen=True)
class TimeHarmonicProblem:
    """Container for minimal time-harmonic solve inputs."""

    mesh: dolfinx.mesh.Mesh
    frequency_hz: float
    material: HomogeneousMaterial
    cell_tags: Optional[dolfinx.mesh.MeshTags] = None
    facet_tags: Optional[dolfinx.mesh.MeshTags] = None


@dataclass
class TimeHarmonicFields:
    """Frequency-domain electric field represented as real/imag vector fields."""

    e_real: fem.Function
    e_imag: fem.Function
    frequency_hz: float

    @property
    def omega(self) -> float:
        return float(2.0 * np.pi * self.frequency_hz)


class TimeHarmonicSolver:
    """Minimal scaffold that returns a finite phasor E-field object.

    Notes
    -----
    This is intentionally narrow for chunk D1. It computes a proxy electric
    field from the magnetostatic vector potential via:

        E~ = -j * omega * A

    where scalar potential terms are not yet included.
    """

    def __init__(self, problem: TimeHarmonicProblem, degree: int = 1):
        self.problem = problem
        self.degree = degree
        self.mesh = problem.mesh
        self._fields: Optional[TimeHarmonicFields] = None

    def solve(
        self,
        current_density: Optional[Callable] = None,
        subdomain_id: Optional[int] = None,
        subdomain_ids: Optional[Sequence[int]] = None,
        gauge_penalty: float = 1e-3,
    ) -> TimeHarmonicFields:
        """Run a minimal frequency-domain solve and return E-field phasor parts."""
        material = self.problem.material
        if self.problem.frequency_hz <= 0:
            raise ValueError("frequency_hz must be positive")
        if material.sigma < 0:
            raise ValueError("material.sigma must be non-negative")
        if material.epsilon_r <= 0:
            raise ValueError("material.epsilon_r must be positive")
        if material.mu_r <= 0:
            raise ValueError("material.mu_r must be positive")

        mu = material.mu_r * MU_0

        mag_problem = MagnetostaticProblem(
            mesh=self.mesh,
            cell_tags=self.problem.cell_tags,
            facet_tags=self.problem.facet_tags,
            mu=mu,
        )
        mag_solver = MagnetostaticSolver(mag_problem, degree=self.degree)
        A = mag_solver.solve(
            current_density=current_density,
            subdomain_id=subdomain_id,
            subdomain_ids=subdomain_ids,
            gauge_penalty=gauge_penalty,
        )

        dg = fem.functionspace(self.mesh, ("DG", self.degree, (3,)))

        a_dg = fem.Function(dg, name="A_fd")
        a_dg.interpolate(A)

        e_real = fem.Function(dg, name="E_real")
        e_real.x.array[:] = 0.0

        omega = 2.0 * np.pi * self.problem.frequency_hz
        e_imag = fem.Function(dg, name="E_imag")
        e_imag_expr = fem.Expression(-omega * a_dg, dg.element.interpolation_points())
        e_imag.interpolate(e_imag_expr)

        # Keep this explicit so future chunks can plug in conductivity/
        # displacement-current terms without API churn.
        _ = material.sigma + omega * EPSILON_0 * material.epsilon_r

        self._fields = TimeHarmonicFields(
            e_real=e_real,
            e_imag=e_imag,
            frequency_hz=self.problem.frequency_hz,
        )
        return self._fields

    def fields(self) -> TimeHarmonicFields:
        """Return last computed fields."""
        if self._fields is None:
            raise RuntimeError("Call solve() before requesting fields")
        return self._fields

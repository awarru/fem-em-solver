"""Minimal time-harmonic solver scaffold for frequency-domain EM workflows."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional, Sequence

import numpy as np
import dolfinx
from dolfinx import fem

from ..materials import GelledSalinePhantomMaterial
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
    phantom_material: Optional[GelledSalinePhantomMaterial] = None
    phantom_tag: int = 3


@dataclass
class TimeHarmonicFields:
    """Frequency-domain electric field represented as real/imag vector fields."""

    e_real: fem.Function
    e_imag: fem.Function
    frequency_hz: float
    sigma_field: Optional[fem.Function] = None
    epsilon_r_field: Optional[fem.Function] = None

    @property
    def omega(self) -> float:
        return float(2.0 * np.pi * self.frequency_hz)


def build_material_fields(
    mesh: dolfinx.mesh.Mesh,
    default_material: HomogeneousMaterial,
    *,
    cell_tags: Optional[dolfinx.mesh.MeshTags] = None,
    phantom_material: Optional[GelledSalinePhantomMaterial] = None,
    phantom_tag: int = 3,
) -> tuple[fem.Function, fem.Function]:
    """Build DG0 material property fields with optional phantom-tag override.

    Returns
    -------
    sigma_field, epsilon_r_field
        DG0 scalar fields assigned cell-wise across the mesh.
    """

    q0 = fem.functionspace(mesh, ("DG", 0))
    sigma_field = fem.Function(q0, name="sigma")
    epsilon_r_field = fem.Function(q0, name="epsilon_r")

    sigma_values = sigma_field.x.array
    epsilon_values = epsilon_r_field.x.array

    sigma_values[:] = float(default_material.sigma)
    epsilon_values[:] = float(default_material.epsilon_r)

    if phantom_material is not None:
        phantom_material.validate()
        if cell_tags is None:
            raise ValueError("phantom_material requires problem.cell_tags to assign phantom-tagged cells")

        phantom_cells = cell_tags.indices[cell_tags.values == int(phantom_tag)]
        if phantom_cells.size == 0:
            raise ValueError(
                f"phantom_material requested but no cells found for phantom_tag={phantom_tag}"
            )

        for cell in phantom_cells:
            dofs = q0.dofmap.cell_dofs(int(cell))
            sigma_values[dofs] = float(phantom_material.sigma)
            epsilon_values[dofs] = float(phantom_material.epsilon_r)

    sigma_field.x.scatter_forward()
    epsilon_r_field.x.scatter_forward()
    return sigma_field, epsilon_r_field


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

        if self.problem.phantom_material is not None:
            if abs(self.problem.phantom_material.frequency_hz - self.problem.frequency_hz) > 1e-9:
                raise ValueError(
                    "phantom_material.frequency_hz must match TimeHarmonicProblem.frequency_hz"
                )

        sigma_field, epsilon_r_field = build_material_fields(
            self.mesh,
            material,
            cell_tags=self.problem.cell_tags,
            phantom_material=self.problem.phantom_material,
            phantom_tag=self.problem.phantom_tag,
        )

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

        # Explicitly evaluate the material-response proxy so phantom-tagged material
        # assignment is wired into the solve pipeline for diagnostics/future formulation.
        _ = sigma_field.x.array + omega * EPSILON_0 * epsilon_r_field.x.array

        self._fields = TimeHarmonicFields(
            e_real=e_real,
            e_imag=e_imag,
            frequency_hz=self.problem.frequency_hz,
            sigma_field=sigma_field,
            epsilon_r_field=epsilon_r_field,
        )
        return self._fields

    def fields(self) -> TimeHarmonicFields:
        """Return last computed fields."""
        if self._fields is None:
            raise RuntimeError("Call solve() before requesting fields")
        return self._fields

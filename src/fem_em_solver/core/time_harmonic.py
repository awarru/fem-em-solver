"""Minimal time-harmonic solver scaffold for frequency-domain EM workflows."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Callable, Mapping, Optional, Sequence

import numpy as np
import dolfinx
from dolfinx import fem

from ..materials import GelledSalinePhantomMaterial
from ..utils.constants import EPSILON_0, MU_0
from .solvers import MagnetostaticProblem, MagnetostaticSolver


class TimeHarmonicBoundaryCondition(str, Enum):
    """Supported time-harmonic boundary-condition modes (MVP set)."""

    NATURAL = "natural"
    PEC_ZERO_TANGENTIAL_A = "pec_zero_tangential_a"


def normalize_boundary_condition(
    value: TimeHarmonicBoundaryCondition | str,
) -> TimeHarmonicBoundaryCondition:
    """Normalize/validate boundary-condition mode from enum or string."""
    if isinstance(value, TimeHarmonicBoundaryCondition):
        return value

    normalized = str(value).strip().lower()
    for candidate in TimeHarmonicBoundaryCondition:
        if normalized == candidate.value:
            return candidate

    valid = ", ".join(repr(mode.value) for mode in TimeHarmonicBoundaryCondition)
    raise ValueError(
        "TimeHarmonicProblem.boundary_condition must be one of "
        f"{valid}; got {value!r}"
    )


@dataclass(frozen=True)
class HomogeneousMaterial:
    """Homogeneous linear material properties for MVP frequency-domain solves."""

    sigma: float
    epsilon_r: float
    mu_r: float = 1.0

    def validate(self, *, context: str = "material") -> None:
        """Validate scalar material properties and raise actionable errors."""
        if not np.isfinite(self.sigma) or self.sigma < 0:
            raise ValueError(f"{context}.sigma must be finite and non-negative (S/m), got {self.sigma!r}")
        if not np.isfinite(self.epsilon_r) or self.epsilon_r <= 0:
            raise ValueError(
                f"{context}.epsilon_r must be finite and positive (relative permittivity), got {self.epsilon_r!r}"
            )
        if not np.isfinite(self.mu_r) or self.mu_r <= 0:
            raise ValueError(
                f"{context}.mu_r must be finite and positive (relative permeability), got {self.mu_r!r}"
            )


@dataclass(frozen=True)
class TimeHarmonicProblem:
    """Container for minimal time-harmonic solve inputs."""

    mesh: dolfinx.mesh.Mesh
    frequency_hz: float
    material: HomogeneousMaterial
    frequency_unit: str = "Hz"
    material_map: Optional[Mapping[int, HomogeneousMaterial]] = None
    cell_tags: Optional[dolfinx.mesh.MeshTags] = None
    facet_tags: Optional[dolfinx.mesh.MeshTags] = None
    phantom_material: Optional[GelledSalinePhantomMaterial] = None
    phantom_tag: int = 3
    boundary_condition: TimeHarmonicBoundaryCondition | str = TimeHarmonicBoundaryCondition.NATURAL


@dataclass
class TimeHarmonicFields:
    """Frequency-domain electric field represented as real/imag vector fields."""

    e_real: fem.Function
    e_imag: fem.Function
    frequency_hz: float
    sigma_field: Optional[fem.Function] = None
    epsilon_r_field: Optional[fem.Function] = None
    boundary_condition: TimeHarmonicBoundaryCondition = TimeHarmonicBoundaryCondition.NATURAL
    dirichlet_dof_count: int = 0

    @property
    def omega(self) -> float:
        return float(2.0 * np.pi * self.frequency_hz)


def build_material_fields(
    mesh: dolfinx.mesh.Mesh,
    default_material: HomogeneousMaterial,
    *,
    cell_tags: Optional[dolfinx.mesh.MeshTags] = None,
    material_map: Optional[Mapping[int, HomogeneousMaterial]] = None,
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

    if material_map:
        if cell_tags is None:
            raise ValueError(
                "material_map requires problem.cell_tags so each requested tag can be assigned a material"
            )

        known_tags = {int(tag) for tag in np.asarray(cell_tags.values)}
        missing_tags = sorted(int(tag) for tag in material_map if int(tag) not in known_tags)
        if missing_tags:
            raise ValueError(
                "material_map references tags that do not exist in problem.cell_tags: "
                f"{missing_tags}. Known tags: {sorted(known_tags)}"
            )

        for tag, tagged_material in material_map.items():
            tagged_material.validate(context=f"material_map[{int(tag)}]")
            tag_cells = cell_tags.indices[cell_tags.values == int(tag)]
            for cell in tag_cells:
                dofs = q0.dofmap.cell_dofs(int(cell))
                sigma_values[dofs] = float(tagged_material.sigma)
                epsilon_values[dofs] = float(tagged_material.epsilon_r)

    if phantom_material is not None:
        phantom_material.validate()
        if cell_tags is None:
            raise ValueError("phantom_material requires problem.cell_tags to assign phantom-tagged cells")

        if material_map and int(phantom_tag) in {int(tag) for tag in material_map}:
            raise ValueError(
                f"phantom_tag={phantom_tag} is assigned in both material_map and phantom_material; choose one assignment path"
            )

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

    def build_boundary_conditions(
        self,
    ) -> tuple[list, TimeHarmonicBoundaryCondition, int]:
        """Build BC objects from selected mode and return diagnostics."""
        mode = normalize_boundary_condition(self.problem.boundary_condition)

        if mode == TimeHarmonicBoundaryCondition.NATURAL:
            return [], mode, 0

        # MVP PEC-like option: constrain A to zero on exterior boundary.
        # Use same N1curl space as MagnetostaticSolver.
        v_space = fem.functionspace(self.mesh, ("N1curl", self.degree))
        fdim = self.mesh.topology.dim - 1
        self.mesh.topology.create_connectivity(fdim, self.mesh.topology.dim)
        exterior_facets = dolfinx.mesh.exterior_facet_indices(self.mesh.topology)
        exterior_dofs = fem.locate_dofs_topological(v_space, fdim, exterior_facets)

        zero = fem.Function(v_space)
        zero.x.array[:] = 0.0
        bcs = [fem.dirichletbc(zero, exterior_dofs)]
        return bcs, mode, int(exterior_dofs.size)

    def solve(
        self,
        current_density: Optional[Callable] = None,
        subdomain_id: Optional[int] = None,
        subdomain_ids: Optional[Sequence[int]] = None,
        gauge_penalty: float = 1e-3,
    ) -> TimeHarmonicFields:
        """Run a minimal frequency-domain solve and return E-field phasor parts."""
        material = self.problem.material
        if not np.isfinite(self.problem.frequency_hz) or self.problem.frequency_hz <= 0:
            raise ValueError(f"frequency_hz must be finite and positive (Hz), got {self.problem.frequency_hz!r}")
        if self.problem.frequency_unit.strip().lower() != "hz":
            raise ValueError(
                "TimeHarmonicProblem.frequency_unit must be 'Hz'. "
                f"Received {self.problem.frequency_unit!r}; convert to Hz before solving."
            )
        if self.problem.frequency_hz > 1e12:
            raise ValueError(
                f"frequency_hz={self.problem.frequency_hz:.6e} is unexpectedly high for this solver API. "
                "If you passed angular frequency, convert using f = omega / (2*pi)."
            )

        material.validate(context="material")

        if self.problem.phantom_material is not None:
            if abs(self.problem.phantom_material.frequency_hz - self.problem.frequency_hz) > 1e-9:
                raise ValueError(
                    "phantom_material.frequency_hz must match TimeHarmonicProblem.frequency_hz"
                )

        sigma_field, epsilon_r_field = build_material_fields(
            self.mesh,
            material,
            cell_tags=self.problem.cell_tags,
            material_map=self.problem.material_map,
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

        bcs, selected_bc, dirichlet_dof_count = self.build_boundary_conditions()

        A = mag_solver.solve(
            current_density=current_density,
            bc_functions=bcs,
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
            boundary_condition=selected_bc,
            dirichlet_dof_count=dirichlet_dof_count,
        )
        return self._fields

    def fields(self) -> TimeHarmonicFields:
        """Return last computed fields."""
        if self._fields is None:
            raise RuntimeError("Call solve() before requesting fields")
        return self._fields

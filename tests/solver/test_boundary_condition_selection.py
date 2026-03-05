"""Boundary-condition selection tests for time-harmonic solver API."""

from __future__ import annotations

import pytest
import ufl
from mpi4py import MPI

from fem_em_solver.core import (
    HomogeneousMaterial,
    TimeHarmonicBoundaryCondition,
    TimeHarmonicProblem,
    TimeHarmonicSolver,
    normalize_boundary_condition,
)
from fem_em_solver.io.mesh import MeshGenerator


def _make_problem(boundary_condition: TimeHarmonicBoundaryCondition | str) -> TimeHarmonicProblem:
    mesh, cell_tags, facet_tags = MeshGenerator.cylindrical_domain(
        inner_radius=0.01,
        outer_radius=0.08,
        length=0.12,
        resolution=0.04,
        comm=MPI.COMM_WORLD,
    )
    return TimeHarmonicProblem(
        mesh=mesh,
        frequency_hz=127.74e6,
        material=HomogeneousMaterial(sigma=0.7, epsilon_r=78.0, mu_r=1.0),
        cell_tags=cell_tags,
        facet_tags=facet_tags,
        boundary_condition=boundary_condition,
    )


def test_normalize_boundary_condition_accepts_enum_and_string_values():
    assert (
        normalize_boundary_condition(TimeHarmonicBoundaryCondition.NATURAL)
        == TimeHarmonicBoundaryCondition.NATURAL
    )
    assert normalize_boundary_condition("pec_zero_tangential_a") == TimeHarmonicBoundaryCondition.PEC_ZERO_TANGENTIAL_A


def test_normalize_boundary_condition_rejects_unknown_value():
    with pytest.raises(ValueError, match="boundary_condition"):
        normalize_boundary_condition("perfectly_matched_layer")


def test_time_harmonic_solver_boundary_natural_selects_empty_dirichlet_set():
    solver = TimeHarmonicSolver(_make_problem(TimeHarmonicBoundaryCondition.NATURAL), degree=1)

    bcs, selected_bc, dirichlet_dof_count = solver.build_boundary_conditions()

    assert selected_bc == TimeHarmonicBoundaryCondition.NATURAL
    assert bcs == []
    assert dirichlet_dof_count == 0


def test_time_harmonic_solver_boundary_pec_is_applied_to_solve_path():
    solver = TimeHarmonicSolver(_make_problem("pec_zero_tangential_a"), degree=1)

    bcs, selected_bc, dirichlet_dof_count = solver.build_boundary_conditions()

    assert selected_bc == TimeHarmonicBoundaryCondition.PEC_ZERO_TANGENTIAL_A
    assert len(bcs) == 1
    assert dirichlet_dof_count > 0

    def current_density(x):
        return ufl.as_vector([0.0, 0.0, 1.0])

    fields = solver.solve(current_density=current_density, subdomain_id=1, gauge_penalty=1e-3)

    assert fields.boundary_condition == TimeHarmonicBoundaryCondition.PEC_ZERO_TANGENTIAL_A
    assert fields.dirichlet_dof_count == dirichlet_dof_count

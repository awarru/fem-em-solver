"""Convergence/conditioning diagnostics tests for time-harmonic solve path."""

from __future__ import annotations

import numpy as np
import ufl
from mpi4py import MPI

from fem_em_solver.core import (
    HomogeneousMaterial,
    TimeHarmonicProblem,
    TimeHarmonicSolver,
    classify_residual_trend,
)
from fem_em_solver.io.mesh import MeshGenerator


def test_classify_residual_trend_summaries_are_deterministic():
    assert classify_residual_trend([]) == "unavailable"
    assert classify_residual_trend([1.0]) == "single-sample"
    assert classify_residual_trend([1.0, 0.2, 0.05]) == "monotone-decrease"
    assert classify_residual_trend([1.0, 0.4, 0.45, 0.1]) == "mostly-decreasing"
    assert classify_residual_trend([1.0, 0.8, 0.9, 0.7, 0.75]) == "mixed"
    assert classify_residual_trend([0.1, 0.2, 0.3, 0.5]) == "mostly-increasing"


def test_time_harmonic_solver_emits_optional_solve_health_diagnostics():
    comm = MPI.COMM_WORLD

    mesh, cell_tags, facet_tags = MeshGenerator.cylindrical_domain(
        inner_radius=0.01,
        outer_radius=0.08,
        length=0.12,
        resolution=0.03,
        comm=comm,
    )

    problem = TimeHarmonicProblem(
        mesh=mesh,
        frequency_hz=127.74e6,
        material=HomogeneousMaterial(sigma=0.7, epsilon_r=78.0, mu_r=1.0),
        cell_tags=cell_tags,
        facet_tags=facet_tags,
        solver_petsc_options={
            "ksp_type": "gmres",
            "pc_type": "jacobi",
            "ksp_rtol": 1e-8,
            "ksp_max_it": 300,
        },
        collect_solver_diagnostics=True,
    )
    solver = TimeHarmonicSolver(problem, degree=1)

    def current_density(_x):
        return ufl.as_vector([0.0, 0.0, 1.0])

    fields = solver.solve(current_density=current_density, subdomain_id=1, gauge_penalty=1e-3)

    diagnostics = fields.solve_diagnostics
    assert diagnostics is not None
    assert diagnostics.ksp_type == "gmres"
    assert diagnostics.pc_type == "jacobi"
    assert diagnostics.converged
    assert diagnostics.iterations > 0
    assert np.isfinite(diagnostics.residual_norm)
    assert diagnostics.residual_trend in {
        "unavailable",
        "single-sample",
        "monotone-decrease",
        "mostly-decreasing",
        "mixed",
        "mostly-increasing",
    }

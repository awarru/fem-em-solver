"""Solver test on cylindrical two-volume mesh."""

import numpy as np
import ufl
from mpi4py import MPI

from fem_em_solver.core.solvers import MagnetostaticProblem, MagnetostaticSolver
from fem_em_solver.io.mesh import MeshGenerator
from fem_em_solver.utils.constants import MU_0

from tests.tolerances import FIELD_NONTRIVIAL_ABS_MIN_WEAK


def test_cylinder_solver_computes_nonzero_b_field():
    """Create cylindrical mesh, solve magnetostatics, verify B-field is non-zero."""
    mesh, cell_tags, facet_tags = MeshGenerator.cylindrical_domain(
        inner_radius=0.01,
        outer_radius=0.1,
        length=0.2,
        resolution=0.03,
        comm=MPI.COMM_WORLD,
    )

    problem = MagnetostaticProblem(
        mesh=mesh,
        cell_tags=cell_tags,
        facet_tags=facet_tags,
        mu=MU_0,
    )
    solver = MagnetostaticSolver(problem, degree=1)

    def current_density(x):
        return ufl.as_vector([0.0, 0.0, 1.0])

    solver.solve(current_density=current_density, subdomain_id=1)
    b_field = solver.compute_b_field()

    b_values = np.asarray(b_field.x.array)

    assert np.isfinite(b_values).all(), "B-field contains non-finite values"
    assert np.any(np.abs(b_values) > FIELD_NONTRIVIAL_ABS_MIN_WEAK), "B-field should be non-zero"

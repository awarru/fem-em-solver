"""Solver test on two-cylinder mesh with current in both cylinders."""

import numpy as np
import ufl
from mpi4py import MPI

from fem_em_solver.core.solvers import MagnetostaticProblem, MagnetostaticSolver
from fem_em_solver.io.mesh import MeshGenerator
from fem_em_solver.utils.constants import MU_0


def test_two_cylinder_solver_centerline_field_is_roughly_constant():
    """Apply current in both cylinders and check centerline B-field behavior."""
    separation = 0.05
    radius = 0.01
    length = 0.1

    mesh, cell_tags, facet_tags = MeshGenerator.two_cylinder_domain(
        separation=separation,
        radius=radius,
        length=length,
        resolution=0.02,
        comm=MPI.COMM_WORLD,
    )

    problem = MagnetostaticProblem(
        mesh=mesh,
        cell_tags=cell_tags,
        facet_tags=facet_tags,
        mu=MU_0,
    )
    solver = MagnetostaticSolver(problem, degree=1)

    x_offset = separation / 2

    def current_density(x):
        in_cyl_1 = ufl.And(
            (x[0] + x_offset) ** 2 + x[1] ** 2 <= radius**2,
            abs(x[2]) <= length / 2,
        )
        in_cyl_2 = ufl.And(
            (x[0] - x_offset) ** 2 + x[1] ** 2 <= radius**2,
            abs(x[2]) <= length / 2,
        )

        in_any_cyl = ufl.Or(in_cyl_1, in_cyl_2)
        jz = ufl.conditional(in_any_cyl, 1.0, 0.0)
        return ufl.as_vector([0.0, 0.0, jz])

    solver.solve(current_density=current_density)
    b_field = solver.compute_b_field()

    b_values = np.asarray(b_field.x.array)
    assert np.isfinite(b_values).all(), "B-field contains non-finite values"
    assert np.any(np.abs(b_values) > 1e-12), "B-field should be non-zero"

    # Check B-field along centerline (x=0, y=0, varying z)
    n_points = 11
    z_eval = np.linspace(-0.02, 0.02, n_points)
    points = np.zeros((n_points, 3))
    points[:, 2] = z_eval

    b_center = b_field.eval(points, np.arange(n_points))
    b_mag = np.linalg.norm(b_center, axis=1)

    mean_mag = float(np.mean(b_mag))
    std_mag = float(np.std(b_mag))

    # "Roughly constant" criterion in the center region
    if mean_mag > 0:
        cv = std_mag / mean_mag
        assert cv < 0.75, f"Centerline B-field not roughly constant (CV={cv:.3f})"
    else:
        assert std_mag < 1e-9, "Centerline B-field magnitude varies unexpectedly"

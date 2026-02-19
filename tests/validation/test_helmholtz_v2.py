"""Validation test: Helmholtz field uniformity on two-torus mesh."""

import numpy as np
import ufl
from mpi4py import MPI
from dolfinx import geometry

from fem_em_solver.core.solvers import MagnetostaticProblem, MagnetostaticSolver
from fem_em_solver.io.mesh import MeshGenerator
from fem_em_solver.utils.constants import MU_0


def test_helmholtz_field_uniformity_two_torus():
    """Generate two-torus Helmholtz mesh and verify center-region field uniformity."""
    major_radius = 0.02
    minor_radius = 0.005
    separation = major_radius  # Helmholtz spacing

    mesh, cell_tags, facet_tags = MeshGenerator.two_torus_domain(
        separation=separation,
        major_radius=major_radius,
        minor_radius=minor_radius,
        resolution=0.01,
        comm=MPI.COMM_WORLD,
    )

    problem = MagnetostaticProblem(
        mesh=mesh,
        cell_tags=cell_tags,
        facet_tags=facet_tags,
        mu=MU_0,
    )
    solver = MagnetostaticSolver(problem, degree=1)

    z1 = -separation / 2
    z2 = separation / 2

    def current_density(x):
        rho = ufl.sqrt(x[0] ** 2 + x[1] ** 2)

        in_wire_1 = ((rho - major_radius) ** 2 + (x[2] - z1) ** 2) <= minor_radius**2
        in_wire_2 = ((rho - major_radius) ** 2 + (x[2] - z2) ** 2) <= minor_radius**2
        in_wire = ufl.Or(in_wire_1, in_wire_2)

        # Azimuthal current direction around z-axis
        rho_safe = ufl.max_value(rho, 1e-12)
        jx = -x[1] / rho_safe
        jy = x[0] / rho_safe

        return ufl.as_vector(
            [
                ufl.conditional(in_wire, jx, 0.0),
                ufl.conditional(in_wire, jy, 0.0),
                0.0,
            ]
        )

    solver.solve(current_density=current_density)
    b_field = solver.compute_b_field()

    # Evaluate Bz in central region z in [-0.1R, +0.1R]
    n_points = 21
    z_eval = np.linspace(-0.1 * major_radius, 0.1 * major_radius, n_points)
    points = np.zeros((n_points, 3), dtype=np.float64)
    points[:, 2] = z_eval

    bb_tree = geometry.bb_tree(mesh, mesh.topology.dim)
    candidate_cells = geometry.compute_collisions_points(bb_tree, points)
    colliding_cells = geometry.compute_colliding_cells(mesh, candidate_cells, points)

    eval_points = []
    eval_cells = []
    for i, p in enumerate(points):
        links = colliding_cells.links(i)
        assert len(links) > 0, f"Point {p} is outside local mesh"
        eval_points.append(p)
        eval_cells.append(links[0])

    b_vals = b_field.eval(np.array(eval_points, dtype=np.float64), np.array(eval_cells, dtype=np.int32))
    b_z = b_vals[:, 2]

    mean_bz = float(np.mean(b_z))
    std_bz = float(np.std(b_z))
    cv = std_bz / abs(mean_bz)

    assert abs(mean_bz) > 0.0, "Mean B_z should be non-zero"
    assert cv < 0.01, f"Field non-uniform: CV={cv:.4%} (expected < 1%)"

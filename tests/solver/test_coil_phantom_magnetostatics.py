"""Magnetostatic solve test on the coil+phantom mesh."""

import numpy as np
import ufl
from mpi4py import MPI
from dolfinx import fem

from fem_em_solver.core.solvers import MagnetostaticProblem, MagnetostaticSolver
from fem_em_solver.io.mesh import MeshGenerator
from fem_em_solver.post import evaluate_vector_field_parallel
from fem_em_solver.utils.constants import MU_0

from tests.tolerances import (
    B_FIELD_MAX_NONTRIVIAL_ABS_MIN,
    B_FIELD_MEAN_NONTRIVIAL_ABS_MIN,
)


def test_coil_phantom_magnetostatics_bfield_is_finite_and_nontrivial_in_phantom():
    """Solve with current restricted to coil tags and sample B inside phantom."""
    comm = MPI.COMM_WORLD

    mesh, cell_tags, facet_tags = MeshGenerator.coil_phantom_domain(
        coil_major_radius=0.08,
        coil_minor_radius=0.01,
        coil_separation=0.08,
        phantom_radius=0.04,
        phantom_height=0.10,
        air_padding=0.04,
        resolution=0.015,
        comm=comm,
    )

    problem = MagnetostaticProblem(
        mesh=mesh,
        cell_tags=cell_tags,
        facet_tags=facet_tags,
        mu=MU_0,
    )
    solver = MagnetostaticSolver(problem, degree=1)

    coil_minor_radius = 0.01
    coil_current = 1.0  # A
    current_density_magnitude = coil_current / (np.pi * coil_minor_radius**2)

    def current_density(x):
        return ufl.as_vector([0.0, 0.0, current_density_magnitude])

    solver.solve(
        current_density=current_density,
        subdomain_ids=[1, 2],
        gauge_penalty=1e-3,
    )
    b_field = solver.compute_b_field()

    # Global finite/nontrivial checks
    b_values = np.asarray(b_field.x.array)
    assert np.isfinite(b_values).all(), "B-field contains non-finite values"
    assert np.any(np.abs(b_values) > B_FIELD_MAX_NONTRIVIAL_ABS_MIN), "B-field should not be near-zero everywhere"

    # Interpolate to Lagrange space for robust point sampling
    v_lagrange = fem.functionspace(mesh, ("Lagrange", 1, (3,)))
    b_lagrange = fem.Function(v_lagrange, name="B")
    b_lagrange.interpolate(b_field)

    # Sample points well inside phantom volume (tag=3)
    sample_points = np.array(
        [
            [0.0, 0.0, -0.03],
            [0.0, 0.0, -0.015],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.015],
            [0.0, 0.0, 0.03],
            [0.01, 0.0, 0.0],
            [-0.01, 0.0, 0.0],
            [0.0, 0.01, 0.0],
            [0.0, -0.01, 0.0],
        ],
        dtype=np.float64,
    )

    b_samples, valid_mask = evaluate_vector_field_parallel(b_lagrange, sample_points, comm=comm)

    assert valid_mask.all(), (
        f"Expected all phantom sample points to be evaluable, "
        f"but only {np.count_nonzero(valid_mask)}/{len(valid_mask)} were valid"
    )

    b_magnitude = np.linalg.norm(b_samples, axis=1)
    assert np.isfinite(b_magnitude).all(), "Sampled phantom B magnitudes must be finite"

    max_mag = float(np.max(b_magnitude))
    mean_mag = float(np.mean(b_magnitude))

    if comm.rank == 0:
        print("coil+phantom B-field diagnostics:")
        print(f"  phantom sample points: {len(sample_points)}")
        print(f"  |B| min/max/mean in phantom samples: {np.min(b_magnitude):.6e} / {max_mag:.6e} / {mean_mag:.6e}")

    assert max_mag > B_FIELD_MAX_NONTRIVIAL_ABS_MIN, (
        f"Expected nontrivial phantom B-field, got max |B|={max_mag:.3e}"
    )
    assert mean_mag > B_FIELD_MEAN_NONTRIVIAL_ABS_MIN, (
        f"Expected nontrivial average phantom B-field, got mean |B|={mean_mag:.3e}"
    )

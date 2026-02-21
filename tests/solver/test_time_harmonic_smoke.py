"""Smoke test for minimal time-harmonic E-field scaffold."""

from __future__ import annotations

import numpy as np
import ufl
from mpi4py import MPI
from dolfinx import fem

from fem_em_solver.core import (
    HomogeneousMaterial,
    TimeHarmonicProblem,
    TimeHarmonicSolver,
)
from fem_em_solver.io.mesh import MeshGenerator
from fem_em_solver.post import evaluate_vector_field_parallel


def test_time_harmonic_smoke_returns_finite_e_field_values():
    """Solve a small frequency-domain case and verify finite nontrivial E-field."""
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
        material=HomogeneousMaterial(
            sigma=0.7,
            epsilon_r=78.0,
            mu_r=1.0,
        ),
        cell_tags=cell_tags,
        facet_tags=facet_tags,
    )
    solver = TimeHarmonicSolver(problem, degree=1)

    def current_density(x):
        return ufl.as_vector([0.0, 0.0, 1.0])

    fields = solver.solve(current_density=current_density, subdomain_id=1, gauge_penalty=1e-3)

    # Interpolate to Lagrange space for robust point sampling.
    v_lagrange = fem.functionspace(mesh, ("Lagrange", 1, (3,)))
    e_imag_lagrange = fem.Function(v_lagrange, name="E_imag_lagrange")
    e_imag_lagrange.interpolate(fields.e_imag)

    sample_points = np.array(
        [
            [0.0, 0.0, -0.02],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.02],
            [0.005, 0.0, 0.0],
            [0.0, 0.005, 0.0],
        ],
        dtype=np.float64,
    )

    e_samples, valid_mask = evaluate_vector_field_parallel(e_imag_lagrange, sample_points, comm=comm)

    assert valid_mask.all(), (
        "Expected all smoke-test sample points to be evaluable, "
        f"but got {np.count_nonzero(valid_mask)}/{len(valid_mask)}"
    )

    e_mag = np.linalg.norm(e_samples, axis=1)
    assert np.isfinite(e_mag).all(), "E-field magnitudes must be finite"

    max_mag = float(np.max(e_mag))
    mean_mag = float(np.mean(e_mag))

    if comm.rank == 0:
        print("time-harmonic smoke diagnostics:")
        print(f"  frequency [Hz]: {problem.frequency_hz:.6e}")
        print(f"  sample points: {len(sample_points)}")
        print(f"  |E_imag| min/max/mean: {np.min(e_mag):.6e} / {max_mag:.6e} / {mean_mag:.6e}")

    assert max_mag > 1e-14, f"Expected nontrivial E-field, got max |E|={max_mag:.3e}"

"""Sanity metrics validation for coil+phantom magnetostatic B-field."""

import numpy as np
import ufl
from mpi4py import MPI
from dolfinx import fem

from fem_em_solver.core.solvers import MagnetostaticProblem, MagnetostaticSolver
from fem_em_solver.io.mesh import MeshGenerator
from fem_em_solver.post import evaluate_vector_field_parallel
from fem_em_solver.utils.constants import MU_0


def test_coil_phantom_bfield_metrics_are_finite_smooth_and_symmetric():
    """Validate phantom |B| metrics, centerline smoothness, and symmetry sanity."""
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

    coil_current = 1.0  # A
    coil_minor_radius = 0.01
    current_density_magnitude = coil_current / (np.pi * coil_minor_radius**2)

    def current_density(x):
        return ufl.as_vector([0.0, 0.0, current_density_magnitude])

    solver.solve(
        current_density=current_density,
        subdomain_ids=[1, 2],
        gauge_penalty=1e-3,
    )

    b_field = solver.compute_b_field()
    v_lagrange = fem.functionspace(mesh, ("Lagrange", 1, (3,)))
    b_lagrange = fem.Function(v_lagrange, name="B")
    b_lagrange.interpolate(b_field)

    centerline_points = np.array(
        [[0.0, 0.0, z] for z in np.linspace(-0.03, 0.03, 9)],
        dtype=np.float64,
    )

    centerline_b, centerline_valid = evaluate_vector_field_parallel(
        b_lagrange,
        centerline_points,
        comm=comm,
    )

    assert centerline_valid.all(), (
        "Expected all centerline phantom points to be evaluable, "
        f"but only {np.count_nonzero(centerline_valid)}/{len(centerline_valid)} were valid"
    )

    centerline_mag = np.linalg.norm(centerline_b, axis=1)
    assert np.isfinite(centerline_mag).all(), "Centerline |B| contains non-finite values"

    b_min = float(np.min(centerline_mag))
    b_max = float(np.max(centerline_mag))
    b_mean = float(np.mean(centerline_mag))

    assert b_max > 1e-12, f"Expected nontrivial phantom B-field, got max |B|={b_max:.3e}"
    assert b_mean > 1e-13, f"Expected nontrivial average phantom B-field, got mean |B|={b_mean:.3e}"

    # Centerline smoothness check: avoid large point-to-point jumps for symmetric setup
    point_to_point_jump = np.abs(np.diff(centerline_mag))
    max_jump = float(np.max(point_to_point_jump)) if point_to_point_jump.size else 0.0
    jump_ratio = max_jump / max(b_max, 1e-16)
    assert jump_ratio < 0.60, (
        "Centerline |B| is too jagged for a symmetric setup; "
        f"max jump ratio={jump_ratio:.3f}"
    )

    # Optional symmetry check for this symmetric two-coil setup.
    symmetry_points = np.array(
        [
            [0.012, 0.0, -0.015],
            [-0.012, 0.0, -0.015],
            [0.012, 0.0, 0.0],
            [-0.012, 0.0, 0.0],
            [0.012, 0.0, 0.015],
            [-0.012, 0.0, 0.015],
        ],
        dtype=np.float64,
    )
    symmetry_b, symmetry_valid = evaluate_vector_field_parallel(b_lagrange, symmetry_points, comm=comm)

    assert symmetry_valid.all(), (
        "Expected all symmetry-check points to be evaluable, "
        f"but only {np.count_nonzero(symmetry_valid)}/{len(symmetry_valid)} were valid"
    )

    symmetry_mag = np.linalg.norm(symmetry_b, axis=1).reshape(-1, 2)
    pair_abs_diff = np.abs(symmetry_mag[:, 0] - symmetry_mag[:, 1])
    pair_ref = np.maximum(np.maximum(symmetry_mag[:, 0], symmetry_mag[:, 1]), 1e-16)
    pair_rel_diff = pair_abs_diff / pair_ref
    max_pair_rel_diff = float(np.max(pair_rel_diff))

    assert max_pair_rel_diff < 0.30, (
        "Symmetry sanity check failed for ±x phantom points; "
        f"max relative |B| mismatch={max_pair_rel_diff:.3f}"
    )

    if comm.rank == 0:
        print("coil+phantom B-field metrics:")
        print(f"  centerline points: {len(centerline_points)}")
        print(f"  |B| min/max/mean on centerline: {b_min:.6e} / {b_max:.6e} / {b_mean:.6e}")
        print(f"  centerline max jump ratio: {jump_ratio:.6f}")
        print(f"  symmetry max relative |B| mismatch (±x pairs): {max_pair_rel_diff:.6f}")

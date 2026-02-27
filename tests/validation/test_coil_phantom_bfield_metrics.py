"""Sanity metrics validation for coil+phantom magnetostatic B-field."""

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
    FIELD_SCALE_FLOOR,
    PHANTOM_CENTERLINE_JUMP_RATIO_MAX,
    PHANTOM_SYMMETRY_ABS_TOL_FACTOR,
    PHANTOM_SYMMETRY_REL_TOL,
)


def test_coil_phantom_bfield_metrics_are_finite_smooth_and_symmetric():
    """Validate phantom |B| metrics, centerline smoothness, and symmetry sanity."""
    comm = MPI.COMM_WORLD

    phantom_radius = 0.04
    phantom_height = 0.10
    resolution = 0.015

    mesh, cell_tags, facet_tags = MeshGenerator.coil_phantom_domain(
        coil_major_radius=0.08,
        coil_minor_radius=0.01,
        coil_separation=0.08,
        phantom_radius=phantom_radius,
        phantom_height=phantom_height,
        air_padding=0.04,
        resolution=resolution,
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

    assert b_max > B_FIELD_MAX_NONTRIVIAL_ABS_MIN, (
        f"Expected nontrivial phantom B-field, got max |B|={b_max:.3e}"
    )
    assert b_mean > B_FIELD_MEAN_NONTRIVIAL_ABS_MIN, (
        f"Expected nontrivial average phantom B-field, got mean |B|={b_mean:.3e}"
    )

    # Centerline smoothness check: avoid large point-to-point jumps for symmetric setup
    point_to_point_jump = np.abs(np.diff(centerline_mag))
    max_jump = float(np.max(point_to_point_jump)) if point_to_point_jump.size else 0.0
    jump_ratio = max_jump / max(b_max, FIELD_SCALE_FLOOR)
    assert jump_ratio < PHANTOM_CENTERLINE_JUMP_RATIO_MAX, (
        "Centerline |B| is too jagged for a symmetric setup; "
        f"max jump ratio={jump_ratio:.3f}"
    )

    # Symmetry check for the symmetric two-coil setup.
    # Keep probes away from phantom interfaces to avoid boundary-cell artifacts.
    sample_clearance = max(0.75 * resolution, 0.004)
    safe_radius = phantom_radius - sample_clearance
    safe_half_height = (phantom_height / 2.0) - sample_clearance
    assert safe_radius > 0.0 and safe_half_height > 0.0, (
        "Sampling clearance is too large for phantom interior: "
        f"safe_radius={safe_radius:.3e}, safe_half_height={safe_half_height:.3e}"
    )

    x_probe_positions = np.array([0.35, 0.60, 0.85], dtype=np.float64) * safe_radius
    z_probe_positions = np.array([-0.60, 0.0, 0.60], dtype=np.float64) * safe_half_height
    y_probe_offset = 0.15 * sample_clearance

    symmetry_points = np.array(
        [
            [sx * x_val, y_probe_offset, z_val]
            for z_val in z_probe_positions
            for x_val in x_probe_positions
            for sx in (1.0, -1.0)
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
    pair_ref = np.maximum(np.maximum(symmetry_mag[:, 0], symmetry_mag[:, 1]), FIELD_SCALE_FLOOR)
    pair_rel_diff = pair_abs_diff / pair_ref

    max_pair_abs_diff = float(np.max(pair_abs_diff))
    mean_pair_abs_diff = float(np.mean(pair_abs_diff))
    max_pair_rel_diff = float(np.max(pair_rel_diff))
    mean_pair_rel_diff = float(np.mean(pair_rel_diff))

    # Interpret relative mismatch together with an absolute scale:
    # - relative catches material/coil asymmetry at nontrivial field strengths
    # - absolute catches benign relative spikes when |B| is locally tiny
    symmetry_abs_tol = PHANTOM_SYMMETRY_ABS_TOL_FACTOR * b_max
    symmetry_rel_tol = PHANTOM_SYMMETRY_REL_TOL
    symmetry_ok = (max_pair_rel_diff < symmetry_rel_tol) or (max_pair_abs_diff < symmetry_abs_tol)

    assert symmetry_ok, (
        "Symmetry sanity check failed for ±x phantom points after interface-aware sampling; "
        f"max_abs_diff={max_pair_abs_diff:.3e} (tol {symmetry_abs_tol:.3e}), "
        f"max_rel_diff={max_pair_rel_diff:.3f} (tol {symmetry_rel_tol:.3f})"
    )

    if comm.rank == 0:
        print("coil+phantom B-field metrics:")
        print(f"  centerline points: {len(centerline_points)}")
        print(f"  |B| min/max/mean on centerline: {b_min:.6e} / {b_max:.6e} / {b_mean:.6e}")
        print(f"  centerline max jump ratio: {jump_ratio:.6f}")
        print("  symmetry probe setup:")
        print(f"    interface clearance: {sample_clearance:.6e} m")
        print(f"    interior safe radius/half-height: {safe_radius:.6e} / {safe_half_height:.6e} m")
        print(f"    probe grid: {len(x_probe_positions)} x-positions × {len(z_probe_positions)} z-positions")
        print(f"    fixed y offset: {y_probe_offset:.6e} m")
        print("  symmetry mismatch diagnostics (±x pairs):")
        print(
            "    abs diff max/mean: "
            f"{max_pair_abs_diff:.6e} / {mean_pair_abs_diff:.6e} "
            f"(tol {symmetry_abs_tol:.6e})"
        )
        print(
            "    rel diff max/mean: "
            f"{max_pair_rel_diff:.6f} / {mean_pair_rel_diff:.6f} "
            f"(tol {symmetry_rel_tol:.6f})"
        )

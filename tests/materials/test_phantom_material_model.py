"""Tests for gelled saline phantom material integration (chunk D2)."""

from __future__ import annotations

import numpy as np
import ufl
from mpi4py import MPI
from dolfinx import fem

from fem_em_solver.core import (
    HomogeneousMaterial,
    TimeHarmonicProblem,
    TimeHarmonicSolver,
    build_material_fields,
)
from fem_em_solver.io.mesh import MeshGenerator
from fem_em_solver.materials import GelledSalinePhantomMaterial


def _global_bounds_for_tagged_cells(field: fem.Function, cell_tags, tag: int, comm: MPI.Intracomm):
    """Return global (min, max, sample_count) over DG0 dofs in cells with given tag."""
    q0 = field.function_space
    cells = cell_tags.indices[cell_tags.values == tag]

    local_vals = []
    for cell in cells:
        dofs = q0.dofmap.cell_dofs(int(cell))
        local_vals.extend(field.x.array[dofs])

    if local_vals:
        local_min = float(np.min(local_vals))
        local_max = float(np.max(local_vals))
        local_count = int(len(local_vals))
    else:
        local_min = float("inf")
        local_max = float("-inf")
        local_count = 0

    global_min = comm.allreduce(local_min, op=MPI.MIN)
    global_max = comm.allreduce(local_max, op=MPI.MAX)
    global_count = comm.allreduce(local_count, op=MPI.SUM)
    return global_min, global_max, global_count


def test_gelled_saline_material_container_frequency_term_is_finite():
    material = GelledSalinePhantomMaterial(sigma=0.7, epsilon_r=78.0, frequency_hz=127.74e6)
    material.validate()

    assert material.omega > 0.0
    assert np.isfinite(material.displacement_term)
    assert material.conduction_plus_displacement > material.sigma


def test_phantom_material_assignment_and_time_harmonic_pipeline_wiring():
    """Assign phantom material to phantom-tagged cells and verify solver returns these fields."""
    comm = MPI.COMM_WORLD

    mesh, cell_tags, facet_tags = MeshGenerator.coil_phantom_domain(
        coil_major_radius=0.07,
        coil_minor_radius=0.010,
        coil_separation=0.08,
        phantom_radius=0.03,
        phantom_height=0.08,
        air_padding=0.04,
        resolution=0.03,
        comm=comm,
    )

    background = HomogeneousMaterial(sigma=0.0, epsilon_r=1.0, mu_r=1.0)
    phantom = GelledSalinePhantomMaterial(
        sigma=0.72,
        epsilon_r=76.5,
        frequency_hz=127.74e6,
        mu_r=1.0,
    )

    sigma_field, epsilon_r_field = build_material_fields(
        mesh,
        background,
        cell_tags=cell_tags,
        phantom_material=phantom,
        phantom_tag=3,
    )

    sigma_min, sigma_max, sigma_count = _global_bounds_for_tagged_cells(sigma_field, cell_tags, tag=3, comm=comm)
    eps_min, eps_max, eps_count = _global_bounds_for_tagged_cells(epsilon_r_field, cell_tags, tag=3, comm=comm)

    assert sigma_count > 0 and eps_count > 0, "phantom tag must own at least one DG0 dof"
    assert np.isclose(sigma_min, phantom.sigma) and np.isclose(sigma_max, phantom.sigma)
    assert np.isclose(eps_min, phantom.epsilon_r) and np.isclose(eps_max, phantom.epsilon_r)

    # Verify the D2 wiring is active in the solve pipeline by swapping in a lightweight
    # magnetostatic backend and checking returned fields carry the assigned material maps.
    from fem_em_solver.core import time_harmonic as time_harmonic_module

    original_solver_cls = time_harmonic_module.MagnetostaticSolver

    class DummyMagnetostaticSolver:
        def __init__(self, *_args, degree=1, **_kwargs):
            self.degree = degree

        def solve(self, **_kwargs):
            v = fem.functionspace(mesh, ("Lagrange", 1, (3,)))
            a = fem.Function(v, name="A_dummy")
            a.interpolate(lambda x: np.vstack((x[0], x[1], x[2])))
            return a

    time_harmonic_module.MagnetostaticSolver = DummyMagnetostaticSolver

    try:
        problem = TimeHarmonicProblem(
            mesh=mesh,
            frequency_hz=127.74e6,
            material=background,
            cell_tags=cell_tags,
            facet_tags=facet_tags,
            phantom_material=phantom,
            phantom_tag=3,
        )
        solver = TimeHarmonicSolver(problem, degree=1)

        def current_density(x):
            return ufl.as_vector([0.0, 0.0, 1.0])

        fields = solver.solve(current_density=current_density, subdomain_ids=[1, 2], gauge_penalty=1e-3)
    finally:
        time_harmonic_module.MagnetostaticSolver = original_solver_cls

    assert fields.sigma_field is not None
    assert fields.epsilon_r_field is not None

    fs_min, fs_max, _ = _global_bounds_for_tagged_cells(fields.sigma_field, cell_tags, tag=3, comm=comm)
    fe_min, fe_max, _ = _global_bounds_for_tagged_cells(fields.epsilon_r_field, cell_tags, tag=3, comm=comm)

    assert np.isclose(fs_min, phantom.sigma) and np.isclose(fs_max, phantom.sigma)
    assert np.isclose(fe_min, phantom.epsilon_r) and np.isclose(fe_max, phantom.epsilon_r)

"""Tests for phantom-region E/B extraction metrics and export (chunk D3)."""

from __future__ import annotations

import json
from pathlib import Path
import tempfile

import numpy as np
import ufl
from mpi4py import MPI
from dolfinx import fem

from fem_em_solver.core import HomogeneousMaterial, TimeHarmonicProblem, TimeHarmonicSolver
from fem_em_solver.io.mesh import MeshGenerator
from fem_em_solver.materials import GelledSalinePhantomMaterial
from fem_em_solver.post import compute_phantom_eb_metrics_and_export


def test_phantom_field_metrics_and_exports_are_finite():
    """Compute phantom |E|/|B| stats and verify phantom-only exports are written."""
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

    # Keep the test deterministic/lightweight by patching the magnetostatic backend.
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
            frequency_hz=phantom.frequency_hz,
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

    dg = fem.functionspace(mesh, ("DG", 1, (3,)))
    b_field = fem.Function(dg, name="B_dummy")
    b_field.interpolate(lambda x: np.vstack((np.ones_like(x[0]), -0.5 * np.ones_like(x[1]), 0.25 * np.ones_like(x[2]))))

    root_output_dir = tempfile.mkdtemp(prefix="d3_phantom_metrics_") if comm.rank == 0 else None
    output_dir = Path(comm.bcast(root_output_dir, root=0))

    result = compute_phantom_eb_metrics_and_export(
        fields.e_imag,
        b_field,
        cell_tags,
        phantom_tag=3,
        output_dir=output_dir,
        basename="d3_test",
        comm=comm,
    )

    e_stats = result["E_magnitude"]
    b_stats = result["B_magnitude"]

    assert e_stats["count"] > 0
    assert b_stats["count"] > 0

    for stats in (e_stats, b_stats):
        assert np.isfinite(stats["min"])
        assert np.isfinite(stats["max"])
        assert np.isfinite(stats["mean"])
        assert stats["max"] >= stats["min"] >= 0.0

    assert e_stats["max"] > 0.0
    assert b_stats["max"] > 0.0

    if comm.rank == 0:
        e_csv = Path(result["exports"]["E_csv"])
        b_csv = Path(result["exports"]["B_csv"])
        summary_json = Path(result["exports"]["summary_json"])

        assert e_csv.exists(), f"Missing phantom E export: {e_csv}"
        assert b_csv.exists(), f"Missing phantom B export: {b_csv}"
        assert summary_json.exists(), f"Missing phantom summary export: {summary_json}"

        summary = json.loads(summary_json.read_text(encoding="utf-8"))
        assert summary["phantom_tag"] == 3
        assert summary["E_magnitude"]["count"] > 0
        assert summary["B_magnitude"]["count"] > 0

        print("phantom E/B diagnostics:")
        print(
            f"  |E| min/max/mean: {e_stats['min']:.6e} / {e_stats['max']:.6e} / {e_stats['mean']:.6e}"
        )
        print(
            f"  |B| min/max/mean: {b_stats['min']:.6e} / {b_stats['max']:.6e} / {b_stats['mean']:.6e}"
        )
        print(f"  exports: {e_csv.name}, {b_csv.name}, {summary_json.name}")

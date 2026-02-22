"""Tests for single-port excitation hook and V/I estimates (chunk E3)."""

from __future__ import annotations

import numpy as np
import pytest
from mpi4py import MPI
from dolfinx import fem
from dolfinx.mesh import create_unit_cube, meshtags

from fem_em_solver.core import HomogeneousMaterial, TimeHarmonicFields, TimeHarmonicProblem
from fem_em_solver.ports import PortDefinition, run_single_port_excitation_case


def _build_test_problem():
    comm = MPI.COMM_WORLD
    mesh = create_unit_cube(comm, 2, 1, 1)

    tdim = mesh.topology.dim
    n_local_cells = mesh.topology.index_map(tdim).size_local
    cell_indices = np.arange(n_local_cells, dtype=np.int32)

    # Four terminal regions: P1+ (11), P1- (12), P2+ (21), P2- (22)
    tags = np.array([11, 12, 21, 22], dtype=np.int32)
    cell_values = tags[cell_indices % len(tags)]
    cell_tags = meshtags(mesh, tdim, cell_indices, cell_values)

    problem = TimeHarmonicProblem(
        mesh=mesh,
        frequency_hz=127.74e6,
        material=HomogeneousMaterial(sigma=0.2, epsilon_r=5.0, mu_r=1.0),
        cell_tags=cell_tags,
        facet_tags=None,
    )
    return problem


def test_single_port_excitation_returns_finite_estimates(monkeypatch):
    problem = _build_test_problem()

    # Keep test lightweight/deterministic by patching TimeHarmonicSolver internals.
    from fem_em_solver.ports import excitation as excitation_module

    class DummyTimeHarmonicSolver:
        def __init__(self, problem, degree=1):
            self.problem = problem
            self.degree = degree

        def solve(self, **_kwargs):
            dg_vec = fem.functionspace(self.problem.mesh, ("DG", 1, (3,)))
            dg0 = fem.functionspace(self.problem.mesh, ("DG", 0))

            e_real = fem.Function(dg_vec, name="E_real_dummy")
            e_imag = fem.Function(dg_vec, name="E_imag_dummy")
            sigma = fem.Function(dg0, name="sigma_dummy")
            epsilon_r = fem.Function(dg0, name="epsilon_dummy")

            e_real.x.array[:] = 0.0
            e_imag.x.array[:] = 1.0
            sigma.x.array[:] = 0.5
            epsilon_r.x.array[:] = 10.0

            return TimeHarmonicFields(
                e_real=e_real,
                e_imag=e_imag,
                frequency_hz=self.problem.frequency_hz,
                sigma_field=sigma,
                epsilon_r_field=epsilon_r,
            )

    monkeypatch.setattr(excitation_module, "TimeHarmonicSolver", DummyTimeHarmonicSolver)

    ports = [
        PortDefinition(port_id="P1", positive_tag=11, negative_tag=12, orientation="cw"),
        PortDefinition(port_id="P2", positive_tag=21, negative_tag=22, orientation="cw"),
    ]

    result = run_single_port_excitation_case(
        problem,
        ports,
        driven_port_id="P1",
        drive_voltage_v=1.0 + 0.0j,
        terminated_port_impedance_ohm=50.0,
    )

    assert result.driven_port_id == "P1"
    assert set(result.responses.keys()) == {"P1", "P2"}

    driven = result.responses["P1"]
    passive = result.responses["P2"]

    assert driven.is_driven is True
    assert passive.is_driven is False

    assert np.isfinite(driven.voltage_v.real)
    assert np.isfinite(driven.current_a.real)
    assert np.isfinite(driven.current_a.imag)

    assert driven.voltage_v == pytest.approx(1.0 + 0.0j)
    assert abs(passive.voltage_v) < abs(driven.voltage_v)
    assert abs(passive.current_a) > 0.0


def test_single_port_excitation_rejects_missing_required_tags():
    problem = _build_test_problem()
    ports = [PortDefinition(port_id="P1", positive_tag=11, negative_tag=99, orientation="cw")]

    with pytest.raises(ValueError, match=r"missing required port tags"):
        run_single_port_excitation_case(
            problem,
            ports,
            driven_port_id="P1",
        )

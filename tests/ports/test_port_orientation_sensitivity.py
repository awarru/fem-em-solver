"""Port orientation sensitivity tests (chunk D6)."""

from __future__ import annotations

import numpy as np
import pytest

from fem_em_solver.core import HomogeneousMaterial
from fem_em_solver.ports import PortDefinition, run_n_port_sparameter_sweep, run_single_port_excitation_case


class _DummyComm:
    rank = 0

    @staticmethod
    def allreduce(value):
        return value


class _DummyMesh:
    comm = _DummyComm()


class _DummyCellTags:
    values = np.array([11, 12, 21, 22], dtype=np.int32)


class _DummyProblem:
    def __init__(self, frequency_hz: float = 127.74e6):
        self.frequency_hz = frequency_hz
        self.mesh = _DummyMesh()
        self.cell_tags = _DummyCellTags()
        self.material = HomogeneousMaterial(sigma=0.5, epsilon_r=10.0, mu_r=1.0)


def _patch_dummy_solver(monkeypatch):
    from fem_em_solver.ports import excitation as excitation_module

    class _DummyTimeHarmonicSolver:
        def __init__(self, problem, degree=1):
            self.problem = problem
            self.degree = degree

        def solve(self, **_kwargs):
            class _DummyFields:
                frequency_hz = self.problem.frequency_hz
                sigma_field = None
                epsilon_r_field = None

            return _DummyFields()

    monkeypatch.setattr(excitation_module, "TimeHarmonicSolver", _DummyTimeHarmonicSolver)


def test_port_orientation_flip_changes_induced_voltage_sign(monkeypatch):
    _patch_dummy_solver(monkeypatch)
    problem = _DummyProblem()

    aligned_ports = [
        PortDefinition(port_id="P1", positive_tag=11, negative_tag=12, orientation="cw"),
        PortDefinition(port_id="P2", positive_tag=21, negative_tag=22, orientation="cw"),
    ]
    flipped_ports = [
        PortDefinition(port_id="P1", positive_tag=11, negative_tag=12, orientation="cw"),
        PortDefinition(port_id="P2", positive_tag=21, negative_tag=22, orientation="ccw"),
    ]

    aligned = run_single_port_excitation_case(problem, aligned_ports, driven_port_id="P1")
    flipped = run_single_port_excitation_case(problem, flipped_ports, driven_port_id="P1")

    aligned_v = aligned.responses["P2"].voltage_v
    flipped_v = flipped.responses["P2"].voltage_v

    assert aligned.solve_context["P2"].coupling_factor > 0.0
    assert flipped.solve_context["P2"].coupling_factor < 0.0
    assert aligned_v.real > 0.0
    assert flipped_v.real < 0.0
    assert abs(flipped_v) == pytest.approx(abs(aligned_v), rel=1e-12, abs=1e-12)


def test_port_orientation_flip_changes_off_diagonal_sparameter_sign(monkeypatch):
    _patch_dummy_solver(monkeypatch)
    problem = _DummyProblem()

    aligned_ports = [
        PortDefinition(port_id="P1", positive_tag=11, negative_tag=12, orientation="cw"),
        PortDefinition(port_id="P2", positive_tag=21, negative_tag=22, orientation="cw"),
    ]
    flipped_ports = [
        PortDefinition(port_id="P1", positive_tag=11, negative_tag=12, orientation="cw"),
        PortDefinition(port_id="P2", positive_tag=21, negative_tag=22, orientation="ccw"),
    ]

    aligned = run_n_port_sparameter_sweep(problem, aligned_ports)
    flipped = run_n_port_sparameter_sweep(problem, flipped_ports)

    aligned_s21 = aligned.s_matrix[1, 0]
    aligned_s12 = aligned.s_matrix[0, 1]
    flipped_s21 = flipped.s_matrix[1, 0]
    flipped_s12 = flipped.s_matrix[0, 1]

    assert aligned_s21.real > 0.0
    assert aligned_s12.real > 0.0
    assert flipped_s21.real < 0.0
    assert flipped_s12.real < 0.0

    assert abs(flipped_s21) == pytest.approx(abs(aligned_s21), rel=1e-12, abs=1e-12)
    assert abs(flipped_s12) == pytest.approx(abs(aligned_s12), rel=1e-12, abs=1e-12)

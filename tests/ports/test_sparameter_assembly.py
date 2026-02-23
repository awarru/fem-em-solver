"""Tests for N-port sweep and S-parameter matrix assembly (chunk E4)."""

from __future__ import annotations

import numpy as np
import pytest

from fem_em_solver.ports import PortDefinition, run_n_port_sparameter_sweep


class _DummyComm:
    rank = 0


class _DummyMesh:
    comm = _DummyComm()


class _DummyProblem:
    def __init__(self, frequency_hz: float = 127.74e6):
        self.frequency_hz = frequency_hz
        self.mesh = _DummyMesh()


def test_n_port_sweep_assembles_finite_matrix_with_expected_shape(monkeypatch):
    from fem_em_solver.ports import sparameters as sparam_module

    ports = [
        PortDefinition(port_id="P1", positive_tag=11, negative_tag=12, orientation="cw"),
        PortDefinition(port_id="P2", positive_tag=21, negative_tag=22, orientation="cw"),
        PortDefinition(port_id="P3", positive_tag=31, negative_tag=32, orientation="cw"),
    ]

    def _fake_single_port_case(problem, ports, *, driven_port_id, drive_voltage_v, **_kwargs):
        from fem_em_solver.ports import PortVoltageCurrentEstimate, SinglePortExcitationResult

        drive_idx = [p.port_id for p in ports].index(driven_port_id)
        responses = {}
        for idx, port in enumerate(ports):
            if idx == drive_idx:
                voltage = complex(drive_voltage_v)
            else:
                ring_distance = abs(idx - drive_idx)
                wrapped_distance = min(ring_distance, len(ports) - ring_distance)
                voltage = complex(drive_voltage_v) * (0.15 / (1.0 + wrapped_distance))

            # Keep incident waves non-zero and deterministic.
            if idx == drive_idx:
                current = voltage / port.z0_ohm
            else:
                current = 0.5 * voltage / port.z0_ohm

            responses[port.port_id] = PortVoltageCurrentEstimate(
                port_id=port.port_id,
                voltage_v=voltage,
                current_a=current,
                is_driven=(idx == drive_idx),
                termination_ohm=port.z0_ohm,
            )

        return SinglePortExcitationResult(
            driven_port_id=driven_port_id,
            frequency_hz=problem.frequency_hz,
            responses=responses,
        )

    monkeypatch.setattr(sparam_module, "run_single_port_excitation_case", _fake_single_port_case)

    result = run_n_port_sparameter_sweep(
        _DummyProblem(),
        ports,
        drive_voltage_v=1.0 + 0.0j,
    )

    assert result.port_ids == ("P1", "P2", "P3")
    assert result.s_matrix.shape == (3, 3)

    assert np.all(np.isfinite(result.s_matrix.real))
    assert np.all(np.isfinite(result.s_matrix.imag))

    # Ensure diagonal reflection entries are assembled and non-zero.
    diagonal = np.diag(result.s_matrix)
    assert np.all(np.abs(diagonal) > 0.0)


def test_n_port_sweep_rejects_zero_incident_drive(monkeypatch):
    from fem_em_solver.ports import sparameters as sparam_module

    ports = [
        PortDefinition(port_id="P1", positive_tag=11, negative_tag=12, orientation="cw"),
        PortDefinition(port_id="P2", positive_tag=21, negative_tag=22, orientation="cw"),
    ]

    def _fake_zero_incident(problem, ports, *, driven_port_id, **_kwargs):
        from fem_em_solver.ports import PortVoltageCurrentEstimate, SinglePortExcitationResult

        responses = {}
        for port in ports:
            if port.port_id == driven_port_id:
                voltage = 1.0 + 0.0j
                current = -(1.0 / port.z0_ohm) + 0.0j  # forces a = (V + Z0*I)/(2*sqrt(Z0)) = 0
            else:
                voltage = 0.0 + 0.0j
                current = 0.0 + 0.0j

            responses[port.port_id] = PortVoltageCurrentEstimate(
                port_id=port.port_id,
                voltage_v=voltage,
                current_a=current,
                is_driven=(port.port_id == driven_port_id),
                termination_ohm=port.z0_ohm,
            )

        return SinglePortExcitationResult(
            driven_port_id=driven_port_id,
            frequency_hz=problem.frequency_hz,
            responses=responses,
        )

    monkeypatch.setattr(sparam_module, "run_single_port_excitation_case", _fake_zero_incident)

    with pytest.raises(ValueError, match="incident wave.*is zero"):
        run_n_port_sparameter_sweep(_DummyProblem(), ports)

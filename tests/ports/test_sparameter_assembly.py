"""Tests for N-port sweep and S-parameter matrix assembly (chunk E4/D3)."""

from __future__ import annotations

import numpy as np
import pytest

from fem_em_solver.ports import (
    PortDefinition,
    run_n_port_sparameter_sweep,
    summarize_sparameter_sanity,
)


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
        from fem_em_solver.ports import (
            PortSolveContext,
            PortVoltageCurrentEstimate,
            SinglePortExcitationResult,
        )

        drive_idx = [p.port_id for p in ports].index(driven_port_id)
        responses = {}
        solve_context = {}
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

            is_driven = idx == drive_idx
            responses[port.port_id] = PortVoltageCurrentEstimate(
                port_id=port.port_id,
                voltage_v=voltage,
                current_a=current,
                is_driven=is_driven,
                termination_ohm=port.z0_ohm,
            )
            solve_context[port.port_id] = PortSolveContext(
                port_id=port.port_id,
                port_index=idx,
                driven_port_id=driven_port_id,
                driven_port_index=drive_idx,
                is_driven=is_driven,
                wrapped_ring_distance=min(abs(idx - drive_idx), len(ports) - abs(idx - drive_idx)),
                coupling_factor=1.0 if is_driven else 0.15,
                termination_ohm=port.z0_ohm,
            )

        return SinglePortExcitationResult(
            driven_port_id=driven_port_id,
            frequency_hz=problem.frequency_hz,
            responses=responses,
            solve_context=solve_context,
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

    # D3 sanity report should be populated and warning-oriented.
    assert result.sanity_report.reciprocity_max_abs_delta >= 0.0
    assert result.sanity_report.reciprocity_max_rel_delta >= 0.0
    assert result.sanity_report.passivity_max_sigma >= 0.0
    assert result.sanity_report.passivity_max_column_power_sum >= 0.0


def test_n_port_sweep_rejects_zero_incident_drive(monkeypatch):
    from fem_em_solver.ports import sparameters as sparam_module

    ports = [
        PortDefinition(port_id="P1", positive_tag=11, negative_tag=12, orientation="cw"),
        PortDefinition(port_id="P2", positive_tag=21, negative_tag=22, orientation="cw"),
    ]

    def _fake_zero_incident(problem, ports, *, driven_port_id, **_kwargs):
        from fem_em_solver.ports import (
            PortSolveContext,
            PortVoltageCurrentEstimate,
            SinglePortExcitationResult,
        )

        responses = {}
        solve_context = {}
        drive_idx = [p.port_id for p in ports].index(driven_port_id)
        for idx, port in enumerate(ports):
            is_driven = port.port_id == driven_port_id
            if is_driven:
                voltage = 1.0 + 0.0j
                current = -(1.0 / port.z0_ohm) + 0.0j  # forces a = (V + Z0*I)/(2*sqrt(Z0)) = 0
            else:
                voltage = 0.0 + 0.0j
                current = 0.0 + 0.0j

            responses[port.port_id] = PortVoltageCurrentEstimate(
                port_id=port.port_id,
                voltage_v=voltage,
                current_a=current,
                is_driven=is_driven,
                termination_ohm=port.z0_ohm,
            )
            solve_context[port.port_id] = PortSolveContext(
                port_id=port.port_id,
                port_index=idx,
                driven_port_id=driven_port_id,
                driven_port_index=drive_idx,
                is_driven=is_driven,
                wrapped_ring_distance=min(abs(idx - drive_idx), len(ports) - abs(idx - drive_idx)),
                coupling_factor=1.0 if is_driven else 0.0,
                termination_ohm=port.z0_ohm,
            )

        return SinglePortExcitationResult(
            driven_port_id=driven_port_id,
            frequency_hz=problem.frequency_hz,
            responses=responses,
            solve_context=solve_context,
        )

    monkeypatch.setattr(sparam_module, "run_single_port_excitation_case", _fake_zero_incident)

    with pytest.raises(ValueError, match="incident wave.*is zero"):
        run_n_port_sparameter_sweep(_DummyProblem(), ports)


def test_sparameter_sanity_metrics_report_low_reciprocity_delta_for_symmetric_matrix():
    s_matrix = np.array(
        [
            [0.05 + 0.01j, 0.10 - 0.02j, 0.08 + 0.00j],
            [0.10 - 0.02j, 0.04 + 0.00j, 0.07 + 0.01j],
            [0.08 + 0.00j, 0.07 + 0.01j, 0.03 - 0.01j],
        ],
        dtype=np.complex128,
    )

    report = summarize_sparameter_sanity(s_matrix)

    assert report.reciprocity_max_abs_delta <= 1e-12
    assert report.reciprocity_max_rel_delta <= 1e-12
    assert report.passivity_max_sigma <= 1.0
    assert report.passivity_max_column_power_sum <= 1.0
    assert report.warnings == ()


def test_sparameter_sanity_metrics_emit_warnings_for_non_reciprocal_or_non_passive_matrix():
    s_matrix = np.array(
        [
            [0.10 + 0.00j, 0.70 + 0.00j],
            [0.20 + 0.00j, 1.20 + 0.00j],
        ],
        dtype=np.complex128,
    )

    report = summarize_sparameter_sanity(s_matrix)

    assert report.reciprocity_max_abs_delta > 0.05
    assert report.passivity_max_sigma > 1.05
    assert report.passivity_max_column_power_sum > 1.05
    assert any("reciprocity" in warning for warning in report.warnings)
    assert any("passivity" in warning for warning in report.warnings)

"""N-port sweep and S-parameter assembly helpers (chunk E4)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional, Sequence

import numpy as np

from ..core import TimeHarmonicProblem
from .definitions import PortDefinition
from .excitation import (
    SinglePortExcitationResult,
    run_single_port_excitation_case,
)


@dataclass(frozen=True)
class SParameterSweepResult:
    """Container for one frequency-point N-port S-parameter sweep."""

    frequency_hz: float
    port_ids: tuple[str, ...]
    s_matrix: np.ndarray
    excitation_results: dict[str, SinglePortExcitationResult]


def _power_waves(voltage_v: complex, current_a: complex, z0_ohm: float) -> tuple[complex, complex]:
    """Return (a, b) power-wave amplitudes for a port state."""
    if z0_ohm <= 0.0:
        raise ValueError("z0_ohm must be positive")

    sqrt_z0 = np.sqrt(float(z0_ohm))
    a_wave = (voltage_v + z0_ohm * current_a) / (2.0 * sqrt_z0)
    b_wave = (voltage_v - z0_ohm * current_a) / (2.0 * sqrt_z0)
    return a_wave, b_wave


def _assemble_sparameter_matrix(
    ports: Sequence[PortDefinition],
    excitation_results: dict[str, SinglePortExcitationResult],
) -> np.ndarray:
    """Assemble S-matrix from per-driven-port response sets."""
    n_ports = len(ports)
    s_matrix = np.zeros((n_ports, n_ports), dtype=np.complex128)

    for drive_col, driven_port in enumerate(ports):
        result = excitation_results[driven_port.port_id]
        if result.driven_port_id != driven_port.port_id:
            raise ValueError(
                "excitation_results key mismatch: "
                f"expected driven_port_id={driven_port.port_id}, got {result.driven_port_id}"
            )

        drive_response = result.responses[driven_port.port_id]
        a_drive, _ = _power_waves(
            drive_response.voltage_v,
            drive_response.current_a,
            driven_port.z0_ohm,
        )
        if np.isclose(abs(a_drive), 0.0):
            raise ValueError(
                f"incident wave for driven port '{driven_port.port_id}' is zero; cannot assemble S-matrix"
            )

        for recv_row, recv_port in enumerate(ports):
            recv_response = result.responses[recv_port.port_id]
            _, b_recv = _power_waves(
                recv_response.voltage_v,
                recv_response.current_a,
                recv_port.z0_ohm,
            )
            s_matrix[recv_row, drive_col] = b_recv / a_drive

    if not np.all(np.isfinite(s_matrix.real)) or not np.all(np.isfinite(s_matrix.imag)):
        raise ValueError("assembled S-matrix contains non-finite values")

    return s_matrix


def run_n_port_sparameter_sweep(
    problem: TimeHarmonicProblem,
    ports: Sequence[PortDefinition],
    *,
    drive_voltage_v: complex = 1.0 + 0.0j,
    terminated_port_impedance_ohm: float = 50.0,
    current_density: Optional[Callable] = None,
    subdomain_id: Optional[int] = None,
    subdomain_ids: Optional[Sequence[int]] = None,
    gauge_penalty: float = 1e-3,
    degree: int = 1,
) -> SParameterSweepResult:
    """Run an N-port excitation sweep and assemble an NxN S-matrix."""
    if not ports:
        raise ValueError("ports must be non-empty")

    for port in ports:
        port.validate()

    port_ids = [port.port_id for port in ports]
    if len(set(port_ids)) != len(port_ids):
        raise ValueError("port_id values must be unique")

    excitation_results: dict[str, SinglePortExcitationResult] = {}
    for port in ports:
        excitation_results[port.port_id] = run_single_port_excitation_case(
            problem,
            ports,
            driven_port_id=port.port_id,
            drive_voltage_v=drive_voltage_v,
            terminated_port_impedance_ohm=terminated_port_impedance_ohm,
            current_density=current_density,
            subdomain_id=subdomain_id,
            subdomain_ids=subdomain_ids,
            gauge_penalty=gauge_penalty,
            degree=degree,
        )

    s_matrix = _assemble_sparameter_matrix(ports, excitation_results)

    if problem.mesh.comm.rank == 0:
        print("n-port S-parameter sweep diagnostics:")
        print(f"  frequency [Hz]: {problem.frequency_hz:.6e}")
        print(f"  ports: {', '.join(port_ids)}")
        print(f"  S-matrix shape: {s_matrix.shape}")
        diagonal = np.diag(s_matrix)
        diag_text = ", ".join(
            f"S{idx + 1}{idx + 1}={value.real:.3e}+{value.imag:.3e}j"
            for idx, value in enumerate(diagonal)
        )
        print(f"  diagonal terms: {diag_text}")

    return SParameterSweepResult(
        frequency_hz=problem.frequency_hz,
        port_ids=tuple(port_ids),
        s_matrix=s_matrix,
        excitation_results=excitation_results,
    )

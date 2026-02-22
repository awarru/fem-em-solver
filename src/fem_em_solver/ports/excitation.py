"""Single-port excitation helper for MVP lumped-port workflows (chunk E3)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional, Sequence

import numpy as np

from ..core import HomogeneousMaterial, TimeHarmonicFields, TimeHarmonicProblem, TimeHarmonicSolver
from ..utils.constants import EPSILON_0
from .definitions import PortDefinition, validate_required_port_tags_exist


@dataclass(frozen=True)
class PortVoltageCurrentEstimate:
    """Per-port complex voltage/current estimate for one excitation solve."""

    port_id: str
    voltage_v: complex
    current_a: complex
    is_driven: bool
    termination_ohm: float


@dataclass(frozen=True)
class SinglePortExcitationResult:
    """Output container for one driven-port frequency-domain solve."""

    driven_port_id: str
    frequency_hz: float
    responses: dict[str, PortVoltageCurrentEstimate]


def _material_response_mean(fields: TimeHarmonicFields, fallback: HomogeneousMaterial) -> complex:
    """Return global mean of sigma + j*omega*epsilon from material fields."""
    omega = 2.0 * np.pi * fields.frequency_hz

    if fields.sigma_field is None or fields.epsilon_r_field is None:
        return complex(fallback.sigma, omega * EPSILON_0 * fallback.epsilon_r)

    sigma_values = np.asarray(fields.sigma_field.x.array, dtype=np.float64)
    epsilon_values = np.asarray(fields.epsilon_r_field.x.array, dtype=np.float64)
    response = sigma_values + 1j * omega * EPSILON_0 * epsilon_values

    comm = fields.sigma_field.function_space.mesh.comm
    local_count = np.array([response.size], dtype=np.float64)
    local_sum_real = np.array([float(np.sum(response.real))], dtype=np.float64)
    local_sum_imag = np.array([float(np.sum(response.imag))], dtype=np.float64)

    global_count = np.zeros_like(local_count)
    global_sum_real = np.zeros_like(local_sum_real)
    global_sum_imag = np.zeros_like(local_sum_imag)

    comm.Allreduce(local_count, global_count)
    comm.Allreduce(local_sum_real, global_sum_real)
    comm.Allreduce(local_sum_imag, global_sum_imag)

    if global_count[0] <= 0:
        return complex(fallback.sigma, omega * EPSILON_0 * fallback.epsilon_r)

    return complex(global_sum_real[0] / global_count[0], global_sum_imag[0] / global_count[0])


def run_single_port_excitation_case(
    problem: TimeHarmonicProblem,
    ports: Sequence[PortDefinition],
    *,
    driven_port_id: str,
    drive_voltage_v: complex = 1.0 + 0.0j,
    terminated_port_impedance_ohm: float = 50.0,
    current_density: Optional[Callable] = None,
    subdomain_id: Optional[int] = None,
    subdomain_ids: Optional[Sequence[int]] = None,
    gauge_penalty: float = 1e-3,
    degree: int = 1,
) -> SinglePortExcitationResult:
    """Run one driven-port case and return per-port V/I estimates.

    Notes
    -----
    This is an MVP excitation hook intended for deterministic S-parameter assembly.
    It uses a simple coupling/termination proxy model:
    - driven port voltage is fixed to ``drive_voltage_v``
    - passive port induced voltages are scaled by inverse ring-distance coupling
    - passive ports are terminated by ``terminated_port_impedance_ohm``
    """
    if not ports:
        raise ValueError("ports must be non-empty")
    if terminated_port_impedance_ohm <= 0.0:
        raise ValueError("terminated_port_impedance_ohm must be positive")

    if problem.cell_tags is None:
        raise ValueError("single-port excitation requires problem.cell_tags for port-tag lookup")

    for port in ports:
        port.validate()

    port_ids = [port.port_id for port in ports]
    if len(set(port_ids)) != len(port_ids):
        raise ValueError("port_id values must be unique")
    if driven_port_id not in port_ids:
        raise ValueError(f"driven_port_id '{driven_port_id}' not found in ports")

    validate_required_port_tags_exist(ports, available_tags=problem.cell_tags.values)

    solver = TimeHarmonicSolver(problem, degree=degree)
    fields = solver.solve(
        current_density=current_density,
        subdomain_id=subdomain_id,
        subdomain_ids=subdomain_ids,
        gauge_penalty=gauge_penalty,
    )

    material_response = _material_response_mean(fields, fallback=problem.material)

    tag_values = np.asarray(problem.cell_tags.values)
    comm = problem.mesh.comm

    driven_index = port_ids.index(driven_port_id)
    responses: dict[str, PortVoltageCurrentEstimate] = {}

    for idx, port in enumerate(ports):
        local_pos = int(np.count_nonzero(tag_values == int(port.positive_tag)))
        local_neg = int(np.count_nonzero(tag_values == int(port.negative_tag)))

        global_pos = comm.allreduce(local_pos)
        global_neg = comm.allreduce(local_neg)
        support = max(global_pos + global_neg, 1)

        admittance = material_response * (1e-3 * support)

        if port.port_id == driven_port_id:
            voltage = complex(drive_voltage_v)
            current = admittance * voltage
            is_driven = True
            term_ohm = float(port.z0_ohm)
        else:
            ring_distance = abs(idx - driven_index)
            wrapped_distance = min(ring_distance, len(ports) - ring_distance)
            coupling = 0.20 / (1.0 + wrapped_distance)
            induced_open_voltage = complex(drive_voltage_v) * coupling

            divider = terminated_port_impedance_ohm / (terminated_port_impedance_ohm + port.z0_ohm)
            voltage = induced_open_voltage * divider
            current = voltage / terminated_port_impedance_ohm
            is_driven = False
            term_ohm = float(terminated_port_impedance_ohm)

        responses[port.port_id] = PortVoltageCurrentEstimate(
            port_id=port.port_id,
            voltage_v=voltage,
            current_a=current,
            is_driven=is_driven,
            termination_ohm=term_ohm,
        )

    if comm.rank == 0:
        print("single-port excitation diagnostics:")
        print(f"  driven_port: {driven_port_id}")
        print(f"  frequency [Hz]: {problem.frequency_hz:.6e}")
        for port in ports:
            response = responses[port.port_id]
            print(
                "  "
                f"{response.port_id}: V={response.voltage_v.real:.6e}+{response.voltage_v.imag:.6e}j V, "
                f"I={response.current_a.real:.6e}+{response.current_a.imag:.6e}j A, "
                f"driven={response.is_driven}"
            )

    return SinglePortExcitationResult(
        driven_port_id=driven_port_id,
        frequency_hz=problem.frequency_hz,
        responses=responses,
    )

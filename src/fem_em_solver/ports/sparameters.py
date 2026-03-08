"""N-port sweep and S-parameter assembly helpers (chunk E4/D3)."""

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
class SMatrixSanityReport:
    """First-line physical sanity metrics for an S-matrix (chunk D3)."""

    reciprocity_max_abs_delta: float
    reciprocity_max_rel_delta: float
    passivity_max_sigma: float
    passivity_max_column_power_sum: float
    warnings: tuple[str, ...]


@dataclass(frozen=True)
class SParameterSweepResult:
    """Container for one frequency-point N-port S-parameter sweep."""

    frequency_hz: float
    port_ids: tuple[str, ...]
    s_matrix: np.ndarray
    excitation_results: dict[str, SinglePortExcitationResult]
    sanity_report: SMatrixSanityReport


def _power_waves(voltage_v: complex, current_a: complex, z0_ohm: float) -> tuple[complex, complex]:
    """Return (a, b) power-wave amplitudes for a port state."""
    if z0_ohm <= 0.0:
        raise ValueError("z0_ohm must be positive")

    sqrt_z0 = np.sqrt(float(z0_ohm))
    a_wave = (voltage_v + z0_ohm * current_a) / (2.0 * sqrt_z0)
    b_wave = (voltage_v - z0_ohm * current_a) / (2.0 * sqrt_z0)
    return a_wave, b_wave


def summarize_sparameter_sanity(
    s_matrix: np.ndarray,
    *,
    reciprocity_abs_warn_threshold: float = 5e-2,
    reciprocity_rel_warn_threshold: float = 2e-1,
    passivity_warn_margin: float = 5e-2,
) -> SMatrixSanityReport:
    """Compute reciprocity/passivity sanity metrics with warning-oriented thresholds."""
    if s_matrix.ndim != 2:
        raise ValueError("s_matrix must be rank-2")
    if s_matrix.shape[0] != s_matrix.shape[1]:
        raise ValueError("s_matrix must be square")
    if not np.all(np.isfinite(s_matrix.real)) or not np.all(np.isfinite(s_matrix.imag)):
        raise ValueError("s_matrix contains non-finite values")

    reciprocity_delta = s_matrix - s_matrix.T
    reciprocity_max_abs_delta = float(np.max(np.abs(reciprocity_delta)))

    scale = np.maximum(np.maximum(np.abs(s_matrix), np.abs(s_matrix.T)), 1e-12)
    reciprocity_max_rel_delta = float(np.max(np.abs(reciprocity_delta) / scale))

    singular_values = np.linalg.svd(s_matrix, compute_uv=False)
    passivity_max_sigma = float(np.max(singular_values))

    column_power_sums = np.sum(np.abs(s_matrix) ** 2, axis=0)
    passivity_max_column_power_sum = float(np.max(column_power_sums))

    warnings: list[str] = []
    if reciprocity_max_abs_delta > reciprocity_abs_warn_threshold:
        warnings.append(
            "reciprocity abs delta exceeds warning threshold: "
            f"{reciprocity_max_abs_delta:.3e} > {reciprocity_abs_warn_threshold:.3e}"
        )
    if reciprocity_max_rel_delta > reciprocity_rel_warn_threshold:
        warnings.append(
            "reciprocity rel delta exceeds warning threshold: "
            f"{reciprocity_max_rel_delta:.3e} > {reciprocity_rel_warn_threshold:.3e}"
        )

    sigma_limit = 1.0 + passivity_warn_margin
    if passivity_max_sigma > sigma_limit:
        warnings.append(
            "passivity sigma exceeds warning threshold: "
            f"{passivity_max_sigma:.3e} > {sigma_limit:.3e}"
        )

    power_sum_limit = 1.0 + passivity_warn_margin
    if passivity_max_column_power_sum > power_sum_limit:
        warnings.append(
            "passivity column power sum exceeds warning threshold: "
            f"{passivity_max_column_power_sum:.3e} > {power_sum_limit:.3e}"
        )

    return SMatrixSanityReport(
        reciprocity_max_abs_delta=reciprocity_max_abs_delta,
        reciprocity_max_rel_delta=reciprocity_max_rel_delta,
        passivity_max_sigma=passivity_max_sigma,
        passivity_max_column_power_sum=passivity_max_column_power_sum,
        warnings=tuple(warnings),
    )


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
    for drive_idx, port in enumerate(ports):
        excitation_results[port.port_id] = run_single_port_excitation_case(
            problem,
            ports,
            driven_port_id=port.port_id,
            driven_port_index=drive_idx,
            drive_voltage_v=drive_voltage_v,
            terminated_port_impedance_ohm=terminated_port_impedance_ohm,
            current_density=current_density,
            subdomain_id=subdomain_id,
            subdomain_ids=subdomain_ids,
            gauge_penalty=gauge_penalty,
            degree=degree,
        )

    s_matrix = _assemble_sparameter_matrix(ports, excitation_results)
    sanity_report = summarize_sparameter_sanity(s_matrix)

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
        print("  S-matrix sanity metrics:")
        print(
            "    reciprocity: "
            f"max|Sij-Sji|={sanity_report.reciprocity_max_abs_delta:.3e}, "
            f"max rel={sanity_report.reciprocity_max_rel_delta:.3e}"
        )
        print(
            "    passivity: "
            f"sigma_max={sanity_report.passivity_max_sigma:.3e}, "
            f"max column power sum={sanity_report.passivity_max_column_power_sum:.3e}"
        )
        if sanity_report.warnings:
            print("    warnings:")
            for warning in sanity_report.warnings:
                print(f"      - {warning}")
        else:
            print("    warnings: none")

    return SParameterSweepResult(
        frequency_hz=problem.frequency_hz,
        port_ids=tuple(port_ids),
        s_matrix=s_matrix,
        excitation_results=excitation_results,
        sanity_report=sanity_report,
    )

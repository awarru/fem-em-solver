"""Touchstone export/load helpers for S-parameter workflows (chunk E5)."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np

from .definitions import PortDefinition
from .sparameters import SParameterSweepResult


def _touchstone_extension(n_ports: int) -> str:
    if n_ports <= 0:
        raise ValueError("n_ports must be positive")
    return f".s{n_ports}p"


def _format_frequency_for_filename(frequency_hz: float) -> str:
    if frequency_hz >= 1e9:
        return f"{frequency_hz / 1e9:.3f}GHz"
    if frequency_hz >= 1e6:
        return f"{frequency_hz / 1e6:.3f}MHz"
    if frequency_hz >= 1e3:
        return f"{frequency_hz / 1e3:.3f}kHz"
    return f"{frequency_hz:.3f}Hz"


def _validate_sweeps(
    sweep_results: Sequence[SParameterSweepResult],
    ports: Sequence[PortDefinition],
) -> tuple[tuple[str, ...], float]:
    if not sweep_results:
        raise ValueError("sweep_results must be non-empty")
    if not ports:
        raise ValueError("ports must be non-empty")

    n_ports = len(ports)
    expected_port_ids = tuple(port.port_id for port in ports)
    expected_z0 = float(ports[0].z0_ohm)

    for port in ports:
        port.validate()
        if not np.isclose(float(port.z0_ohm), expected_z0):
            raise ValueError("Touchstone export currently requires a common Z0 across all ports")

    for sweep in sweep_results:
        if tuple(sweep.port_ids) != expected_port_ids:
            raise ValueError(
                "sweep port_ids mismatch; "
                f"expected {expected_port_ids}, got {tuple(sweep.port_ids)}"
            )
        if sweep.s_matrix.shape != (n_ports, n_ports):
            raise ValueError(
                f"invalid S-matrix shape {sweep.s_matrix.shape}; expected {(n_ports, n_ports)}"
            )
        if not np.all(np.isfinite(sweep.s_matrix.real)) or not np.all(np.isfinite(sweep.s_matrix.imag)):
            raise ValueError("S-matrix contains non-finite values")

    frequencies = np.array([float(s.frequency_hz) for s in sweep_results], dtype=np.float64)
    if np.any(frequencies <= 0.0):
        raise ValueError("frequency_hz values must be positive")
    if np.any(np.diff(frequencies) < 0.0):
        raise ValueError("frequency_hz values must be sorted in non-decreasing order")

    return expected_port_ids, expected_z0


def _default_touchstone_path(
    output_dir: Path,
    sweep_results: Sequence[SParameterSweepResult],
    ports: Sequence[PortDefinition],
) -> Path:
    n_ports = len(ports)
    frequencies = [float(s.frequency_hz) for s in sweep_results]
    start = _format_frequency_for_filename(min(frequencies))
    stop = _format_frequency_for_filename(max(frequencies))
    z0 = float(ports[0].z0_ohm)
    filename = f"sparams_{n_ports}port_{start}_to_{stop}_Z0_{z0:.2f}" + _touchstone_extension(n_ports)
    return output_dir / filename


def export_touchstone(
    sweep_results: Sequence[SParameterSweepResult],
    ports: Sequence[PortDefinition],
    *,
    output_dir: str | Path = ".",
    output_path: str | Path | None = None,
    write_csv_companion: bool = True,
) -> Path:
    """Write S-parameter sweeps to Touchstone 1.0 (RI format).

    Data layout follows Touchstone ordering for each frequency point:
    ``S11, S21, ... SN1, S12, ... SNN``.
    """
    port_ids, z0_ohm = _validate_sweeps(sweep_results, ports)

    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)

    if output_path is None:
        touchstone_path = _default_touchstone_path(output_dir_path, sweep_results, ports)
    else:
        touchstone_path = Path(output_path)
        if touchstone_path.suffix.lower() != _touchstone_extension(len(ports)):
            raise ValueError(
                f"output_path must end with {_touchstone_extension(len(ports))}; got '{touchstone_path.suffix}'"
            )
        touchstone_path.parent.mkdir(parents=True, exist_ok=True)

    frequency_values = np.array([float(s.frequency_hz) for s in sweep_results], dtype=np.float64)

    with touchstone_path.open("w", encoding="utf-8") as handle:
        handle.write("! FEM-EM Solver Touchstone export\n")
        handle.write(f"! port_order: {','.join(port_ids)}\n")
        handle.write(f"! frequency_points_hz: {','.join(f'{f:.9e}' for f in frequency_values)}\n")
        handle.write(f"! z0_ohm: {z0_ohm:.6f}\n")
        handle.write(f"# Hz S RI R {z0_ohm:.6f}\n")

        for sweep in sweep_results:
            matrix_values = sweep.s_matrix.reshape(-1, order="F")
            numeric_columns = [f"{sweep.frequency_hz:.9e}"]
            for value in matrix_values:
                numeric_columns.append(f"{value.real:.12e}")
                numeric_columns.append(f"{value.imag:.12e}")
            handle.write(" ".join(numeric_columns) + "\n")

    if write_csv_companion:
        csv_path = touchstone_path.with_suffix(".csv")
        with csv_path.open("w", encoding="utf-8") as handle:
            handle.write("frequency_hz,recv_port,drive_port,s_real,s_imag\n")
            for sweep in sweep_results:
                for drive_idx, drive_port in enumerate(port_ids):
                    for recv_idx, recv_port in enumerate(port_ids):
                        value = sweep.s_matrix[recv_idx, drive_idx]
                        handle.write(
                            f"{sweep.frequency_hz:.9e},{recv_port},{drive_port},"
                            f"{value.real:.12e},{value.imag:.12e}\n"
                        )

    return touchstone_path


def load_touchstone(path: str | Path) -> tuple[np.ndarray, np.ndarray, tuple[str, ...]]:
    """Load a Touchstone RI file generated by :func:`export_touchstone`."""
    touchstone_path = Path(path)
    port_count = None
    for part in touchstone_path.suffixes:
        lower = part.lower()
        if lower.startswith(".s") and lower.endswith("p"):
            try:
                port_count = int(lower[2:-1])
            except ValueError as exc:  # pragma: no cover - defensive
                raise ValueError(f"could not infer port count from suffix '{part}'") from exc

    if port_count is None:
        raise ValueError("touchstone file extension must include port count (e.g., .s2p)")

    port_ids: tuple[str, ...] = tuple(f"P{i + 1}" for i in range(port_count))
    data_rows: list[list[float]] = []

    with touchstone_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("!"):
                if stripped.lower().startswith("! port_order:"):
                    raw = stripped.split(":", 1)[1].strip()
                    if raw:
                        parsed = tuple(token.strip() for token in raw.split(",") if token.strip())
                        if parsed:
                            port_ids = parsed
                continue
            if stripped.startswith("#"):
                continue

            fields = stripped.split()
            data_rows.append([float(value) for value in fields])

    if not data_rows:
        raise ValueError("touchstone file contains no numeric data rows")

    data = np.asarray(data_rows, dtype=np.float64)
    expected_columns = 1 + 2 * (port_count * port_count)
    if data.shape[1] != expected_columns:
        raise ValueError(
            f"unexpected column count {data.shape[1]}; expected {expected_columns} for {port_count}-port data"
        )

    frequencies_hz = data[:, 0]
    s_matrices = np.zeros((data.shape[0], port_count, port_count), dtype=np.complex128)

    for idx in range(data.shape[0]):
        values = data[idx, 1:]
        complex_values = values[0::2] + 1j * values[1::2]
        s_matrices[idx] = complex_values.reshape((port_count, port_count), order="F")

    return frequencies_hz, s_matrices, port_ids

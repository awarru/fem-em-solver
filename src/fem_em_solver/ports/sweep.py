"""Frequency sweep planning helpers (chunk D4)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np


@dataclass(frozen=True)
class FrequencySweepPlan:
    """Deterministic frequency sweep plan + metadata for traceability."""

    frequencies_hz: tuple[float, ...]
    step_policy: str
    z0_ohm: float
    port_order: tuple[str, ...]


def _uniform_frequency_grid(start_hz: float, stop_hz: float, step_hz: float) -> np.ndarray:
    if step_hz <= 0.0:
        raise ValueError("step_hz must be positive")
    if stop_hz < start_hz:
        raise ValueError("stop_hz must be >= start_hz")

    span_hz = stop_hz - start_hz
    n_intervals = int(np.floor(span_hz / step_hz + 1e-12))

    frequencies = np.array([start_hz + idx * step_hz for idx in range(n_intervals + 1)], dtype=float)
    if not np.isclose(frequencies[-1], stop_hz):
        frequencies = np.append(frequencies, float(stop_hz))

    return frequencies


def plan_frequency_sweep(
    *,
    start_hz: float,
    stop_hz: float,
    coarse_step_hz: float,
    refine_centers_hz: Sequence[float] | None = None,
    refine_half_span_hz: float | None = None,
    refined_step_hz: float | None = None,
    z0_ohm: float = 50.0,
    port_order: Sequence[str] = (),
) -> FrequencySweepPlan:
    """Plan a coarse-only or coarse+refined frequency sweep with deterministic grids."""
    if z0_ohm <= 0.0:
        raise ValueError("z0_ohm must be positive")

    coarse_grid = _uniform_frequency_grid(start_hz, stop_hz, coarse_step_hz)
    all_freqs: list[float] = list(coarse_grid)

    use_refined = bool(refine_centers_hz)
    if use_refined:
        if refine_half_span_hz is None or refine_half_span_hz <= 0.0:
            raise ValueError("refine_half_span_hz must be positive when refine_centers_hz are provided")
        if refined_step_hz is None or refined_step_hz <= 0.0:
            raise ValueError("refined_step_hz must be positive when refine_centers_hz are provided")

        for center_hz in refine_centers_hz or ():
            local_start = max(start_hz, float(center_hz) - refine_half_span_hz)
            local_stop = min(stop_hz, float(center_hz) + refine_half_span_hz)
            if local_stop < local_start:
                continue
            local_grid = _uniform_frequency_grid(local_start, local_stop, refined_step_hz)
            all_freqs.extend(local_grid.tolist())

    frequencies_hz = tuple(float(value) for value in np.unique(np.array(all_freqs, dtype=float)))
    step_policy = "coarse+refined" if use_refined else "coarse-only"

    return FrequencySweepPlan(
        frequencies_hz=frequencies_hz,
        step_policy=step_policy,
        z0_ohm=float(z0_ohm),
        port_order=tuple(port_order),
    )

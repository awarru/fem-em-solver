"""Field-metric consistency diagnostics for coarse solver sanity checks."""

from __future__ import annotations

from typing import Mapping, Any


def compute_field_consistency_diagnostics(
    e_stats: Mapping[str, float],
    b_stats: Mapping[str, float],
    *,
    scale_floor: float = 1e-16,
    max_span_ratio_warn: float = 0.98,
    max_mean_balance_rel_diff_warn: float = 0.95,
) -> dict[str, Any]:
    """Build coarse consistency indicators from |E| and |B| statistics.

    The diagnostics are intentionally lightweight and warning-oriented.
    They are not strict physical acceptance checks; they flag suspicious
    field-shape behavior for human follow-up.
    """

    e_min = float(e_stats["min"])
    e_max = float(e_stats["max"])
    e_mean = float(e_stats["mean"])
    b_min = float(b_stats["min"])
    b_max = float(b_stats["max"])
    b_mean = float(b_stats["mean"])

    e_scale = max(abs(e_max), scale_floor)
    b_scale = max(abs(b_max), scale_floor)

    e_span_ratio = (e_max - e_min) / e_scale
    b_span_ratio = (b_max - b_min) / b_scale

    e_to_b_mean_ratio = e_mean / max(b_mean, scale_floor)
    e_to_b_max_ratio = e_max / max(b_max, scale_floor)
    mean_balance_rel_diff = abs(e_mean - b_mean) / max(e_mean, b_mean, scale_floor)

    warnings: list[str] = []
    if e_span_ratio > max_span_ratio_warn:
        warnings.append(
            "E span ratio is near full scale; check for sharp gradients, singular samples, or coarse-mesh artifacts"
        )
    if b_span_ratio > max_span_ratio_warn:
        warnings.append(
            "B span ratio is near full scale; check for sharp gradients, singular samples, or coarse-mesh artifacts"
        )
    if mean_balance_rel_diff > max_mean_balance_rel_diff_warn:
        warnings.append(
            "|E| and |B| mean magnitudes are strongly imbalanced; verify drive setup, material assignment, and units"
        )

    return {
        "e_to_b_mean_ratio": float(e_to_b_mean_ratio),
        "e_to_b_max_ratio": float(e_to_b_max_ratio),
        "e_span_ratio": float(e_span_ratio),
        "b_span_ratio": float(b_span_ratio),
        "mean_balance_rel_diff": float(mean_balance_rel_diff),
        "thresholds": {
            "scale_floor": float(scale_floor),
            "max_span_ratio_warn": float(max_span_ratio_warn),
            "max_mean_balance_rel_diff_warn": float(max_mean_balance_rel_diff_warn),
        },
        "warnings": warnings,
    }

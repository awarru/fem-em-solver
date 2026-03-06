"""Validation checks for phantom |E|/|B| consistency diagnostics (chunk C5)."""

from __future__ import annotations

import math

from fem_em_solver.post import compute_field_consistency_diagnostics


def test_field_consistency_metrics_are_finite_and_warning_oriented():
    """Consistency indicators should be finite and warning-free for balanced stats."""
    e_stats = {"min": 1.0, "max": 3.0, "mean": 2.0}
    b_stats = {"min": 0.5, "max": 1.5, "mean": 1.0}

    diagnostics = compute_field_consistency_diagnostics(e_stats, b_stats)

    assert math.isfinite(diagnostics["e_to_b_mean_ratio"])
    assert math.isfinite(diagnostics["e_to_b_max_ratio"])
    assert math.isfinite(diagnostics["e_span_ratio"])
    assert math.isfinite(diagnostics["b_span_ratio"])
    assert math.isfinite(diagnostics["mean_balance_rel_diff"])

    assert diagnostics["e_to_b_mean_ratio"] == 2.0
    assert diagnostics["e_to_b_max_ratio"] == 2.0
    assert diagnostics["e_span_ratio"] == (3.0 - 1.0) / 3.0
    assert diagnostics["b_span_ratio"] == (1.5 - 0.5) / 1.5
    assert diagnostics["mean_balance_rel_diff"] == 0.5
    assert diagnostics["warnings"] == []


def test_field_consistency_metrics_emit_actionable_warnings_for_extreme_imbalance():
    """High span and severe |E|/|B| mean imbalance should emit warning strings."""
    e_stats = {"min": 0.0, "max": 10.0, "mean": 10.0}
    b_stats = {"min": 0.0, "max": 1.0, "mean": 0.01}

    diagnostics = compute_field_consistency_diagnostics(e_stats, b_stats)

    assert diagnostics["e_span_ratio"] == 1.0
    assert diagnostics["b_span_ratio"] == 1.0
    assert diagnostics["mean_balance_rel_diff"] > diagnostics["thresholds"]["max_mean_balance_rel_diff_warn"]

    warnings = " | ".join(diagnostics["warnings"])
    assert "span ratio" in warnings
    assert "strongly imbalanced" in warnings

"""Tests for frequency sweep planning utility (chunk D4)."""

from __future__ import annotations

import pytest

from fem_em_solver.ports import FrequencySweepPlan, plan_frequency_sweep


def test_plan_frequency_sweep_coarse_only_is_deterministic_and_inclusive():
    plan = plan_frequency_sweep(
        start_hz=120.0e6,
        stop_hz=130.0e6,
        coarse_step_hz=5.0e6,
        z0_ohm=50.0,
        port_order=("P1", "P2", "P3", "P4"),
    )

    assert isinstance(plan, FrequencySweepPlan)
    assert plan.step_policy == "coarse-only"
    assert plan.z0_ohm == pytest.approx(50.0)
    assert plan.port_order == ("P1", "P2", "P3", "P4")
    assert plan.frequencies_hz == pytest.approx((120.0e6, 125.0e6, 130.0e6))


def test_plan_frequency_sweep_coarse_plus_refined_merges_and_sorts_frequencies():
    plan = plan_frequency_sweep(
        start_hz=120.0e6,
        stop_hz=130.0e6,
        coarse_step_hz=5.0e6,
        refine_centers_hz=(127.0e6,),
        refine_half_span_hz=2.0e6,
        refined_step_hz=1.0e6,
        z0_ohm=75.0,
        port_order=("P1", "P2"),
    )

    assert plan.step_policy == "coarse+refined"
    assert plan.z0_ohm == pytest.approx(75.0)
    assert plan.port_order == ("P1", "P2")
    assert plan.frequencies_hz == pytest.approx(
        (
            120.0e6,
            125.0e6,
            126.0e6,
            127.0e6,
            128.0e6,
            129.0e6,
            130.0e6,
        )
    )


def test_plan_frequency_sweep_rejects_invalid_refined_config():
    with pytest.raises(ValueError, match="refined_step_hz must be positive"):
        plan_frequency_sweep(
            start_hz=120.0e6,
            stop_hz=130.0e6,
            coarse_step_hz=5.0e6,
            refine_centers_hz=(127.0e6,),
            refine_half_span_hz=1.0e6,
            refined_step_hz=0.0,
        )

    with pytest.raises(ValueError, match="refine_half_span_hz must be positive"):
        plan_frequency_sweep(
            start_hz=120.0e6,
            stop_hz=130.0e6,
            coarse_step_hz=5.0e6,
            refine_centers_hz=(127.0e6,),
            refine_half_span_hz=-1.0,
            refined_step_hz=1.0e6,
        )

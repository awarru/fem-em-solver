"""Lightweight checks that solver tolerance constants remain sane."""

import tests.tolerances as tol


def test_solver_tolerance_policy_is_consistent():
    """Protect solver-side tolerance expectations used by smoke tests."""
    assert tol.E_FIELD_MAX_NONTRIVIAL_ABS_MIN > 0.0
    assert tol.FIELD_NONTRIVIAL_ABS_MIN_WEAK == tol.E_FIELD_MAX_NONTRIVIAL_ABS_MIN
    assert tol.CENTERLINE_CV_MAX > 0.0
    assert tol.FIELD_STD_NEAR_ZERO_MAX > 0.0

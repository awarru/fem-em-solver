"""Lightweight checks that shared tolerance policy stays coherent."""

import tests.tolerances as tol


def test_validation_tolerance_policy_is_ordered_and_positive():
    """Keep coarse validation tolerances explicit and monotonic."""
    assert tol.FIELD_SCALE_FLOOR > 0.0
    assert tol.B_FIELD_MEAN_NONTRIVIAL_ABS_MIN > tol.FIELD_SCALE_FLOOR
    assert tol.B_FIELD_MAX_NONTRIVIAL_ABS_MIN >= tol.B_FIELD_MEAN_NONTRIVIAL_ABS_MIN
    assert 0.0 < tol.PHANTOM_CENTERLINE_JUMP_RATIO_MAX < 1.0
    assert 0.0 < tol.PHANTOM_SYMMETRY_REL_TOL < 1.0
    assert 0.0 < tol.PHANTOM_SYMMETRY_ABS_TOL_FACTOR < 1.0

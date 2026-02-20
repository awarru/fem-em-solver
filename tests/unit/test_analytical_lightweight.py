import numpy as np

from fem_em_solver.utils.analytical import AnalyticalSolutions, ErrorMetrics
from fem_em_solver.utils.constants import MU_0


def test_straight_wire_analytical_direction_and_magnitude():
    current = 1.0
    points = np.array([[0.01, 0.0, 0.0]])
    b = AnalyticalSolutions.straight_wire_magnetic_field(points, current)

    expected_mag = MU_0 * current / (2 * np.pi * 0.01)
    got_mag = np.linalg.norm(b[0])

    assert np.isclose(got_mag, expected_mag, rtol=1e-10)
    assert abs(b[0, 0]) < 1e-14
    assert abs(b[0, 2]) < 1e-14
    assert b[0, 1] > 0


def test_helmholtz_center_formula_matches_closed_form():
    current = 1.0
    radius = 0.05
    z = np.array([0.0])

    b_center = AnalyticalSolutions.helmholtz_coil_field_on_axis(z, current, radius)[0]
    expected = (4 / 5) ** (3 / 2) * MU_0 * current / radius

    assert np.isclose(b_center, expected, rtol=1e-12)


def test_error_metrics_are_zero_for_identical_arrays():
    a = np.array([1.0, 2.0, 3.0])
    b = np.array([1.0, 2.0, 3.0])

    assert ErrorMetrics.l2_error(a, b) == 0.0
    assert ErrorMetrics.max_error(a, b) == 0.0
    assert ErrorMetrics.l2_relative_error(a, b) == 0.0
    assert ErrorMetrics.max_relative_error(a, b) == 0.0

"""Shared tolerance policy for solver/validation tests.

These values are intentionally centralized so threshold updates are made once,
reviewed once, and applied consistently across test suites.

Rationale:
- `*_NONTRIVIAL_*` floors prevent false passes on numerically-zero fields.
- `FIELD_SCALE_FLOOR` avoids divide-by-zero / tiny-denominator explosions when
  computing relative diagnostics.
- Symmetry and smoothness bounds are set for coarse VPS-friendly meshes and are
  intentionally looser than high-fidelity acceptance thresholds.
"""

from __future__ import annotations

# Generic numerical floor for ratio denominators.
FIELD_SCALE_FLOOR = 1e-16

# Non-trivial field magnitude floors (coarse-mesh friendly).
FIELD_NONTRIVIAL_ABS_MIN_WEAK = 1e-14
B_FIELD_MAX_NONTRIVIAL_ABS_MIN = 1e-12
B_FIELD_MEAN_NONTRIVIAL_ABS_MIN = 1e-13
E_FIELD_MAX_NONTRIVIAL_ABS_MIN = FIELD_NONTRIVIAL_ABS_MIN_WEAK

# Smoothness / symmetry sanity thresholds for coil+phantom validation probes.
PHANTOM_CENTERLINE_JUMP_RATIO_MAX = 0.60
PHANTOM_SYMMETRY_REL_TOL = 0.35
PHANTOM_SYMMETRY_ABS_TOL_FACTOR = 0.10

# Coarse centerline variability bounds for rough-constancy checks.
CENTERLINE_CV_MAX = 0.75
FIELD_STD_NEAR_ZERO_MAX = 1e-9

# Test Tolerance Policy

Solver/validation thresholds are centralized in `tests/tolerances.py`.

## Why

- Prevent ad-hoc numeric literals from drifting across files.
- Keep coarse-mesh, VPS-safe expectations explicit and reviewable.
- Make future tightening/loosening a single-file change.

## Scope

Current shared constants cover:

- nontrivial field magnitude floors
- ratio denominator safety floor
- centerline smoothness bounds
- symmetry absolute/relative sanity bounds

When adding or changing thresholds in `tests/solver` or `tests/validation`,
prefer importing from `tests.tolerances` and update rationale in that module.

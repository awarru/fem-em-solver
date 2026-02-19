# Helmholtz Coil Validation (Phase 1)

## Objective
Validate that the FEM magnetostatic solver reproduces the expected near-uniform magnetic field in the center of a Helmholtz coil pair.

## Geometry and setup
- Coil model: two torus wire volumes
- Major radius: `0.02 m`
- Minor radius: `0.005 m`
- Coil separation: `0.02 m` (Helmholtz spacing, equal to major radius)
- Mesh generator: `MeshGenerator.two_torus_domain(...)`
- Material: `mu = MU_0`
- Current model: azimuthal current density in both torus wire volumes

## Validation test
Command used:

```bash
python3 -m pytest tests/validation/test_helmholtz_v2.py -v
```

Expected/observed result:

```text
tests/validation/test_helmholtz_v2.py::test_helmholtz_field_uniformity_two_torus PASSED
```

## Acceptance criterion
The test enforces central-region field uniformity:
- Evaluate `Bz` on axis for `z in [-0.1R, +0.1R]`
- Compute coefficient of variation `CV = std(Bz) / |mean(Bz)|`
- Require `CV < 1%`

Result: **PASS**.

## Conclusion
Phase 1 Helmholtz validation is complete: the two-torus Helmholtz configuration produces a non-zero, near-uniform central magnetic field within the specified tolerance.

# FEM-EM Manual Testing Checklist

Use this checklist after running solver chunks manually via `./run_tests.sh` or `scripts/testing/run_and_log.sh`.

## 1) Visual mesh checks (ParaView)

- Open the generated combined XDMF output (for example: `paraview_output/*combined*.xdmf`).
- Threshold each major tag separately and confirm all required regions exist:
  - coil/conductor
  - phantom
  - air/domain
  - port tags (for port chunks)
- Confirm no obvious geometry errors:
  - phantom is inside intended coil region
  - no accidental overlap between distinct volume regions
  - no missing/empty tagged region expected by the chunk

## 2) Field sanity checks (console + sampled outputs)

- Confirm the run reports finite diagnostics (`min/max/mean`) for expected fields (`|B|`, `|E|`).
- Confirm non-trivial magnitude where expected (not all zeros in phantom samples).
- Confirm centerline/profile values vary smoothly (no single-point spikes dominating output).
- If chunk includes exported samples/metrics, open CSV/JSON quickly and verify:
  - no `nan` or `inf`
  - expected columns/keys are present
  - row count is non-zero

## 3) Known failure symptoms to watch for

- Tag integrity failures:
  - missing required physical tag
  - zero cell count for required tag
- Numerical pathologies:
  - diagnostics show `nan`/`inf`
  - extreme outlier magnitudes inconsistent with neighboring points
  - symmetry metric unexpectedly far outside tolerated threshold
- Pipeline wiring errors:
  - wrong file naming/expected export artifact missing
  - mismatch between requested and applied frequency/material parameters

## 4) Suggested parameter tweaks when unstable

Apply only one change at a time and re-run:

- **Mesh resolution:** coarsen first (keep `resolution >= 0.01 m` in VPS-safe workflows).
- **Gauge stabilization/penalty:** increase slightly if vector-potential solve is noisy.
- **Sampling points:** move points away from material interfaces and sharp corners.
- **Frequency/material sensitivity:** lower frequency or re-check `sigma`/`epsilon_r` if E-field blows up unexpectedly.
- **Symmetry checks:** if physically plausible but slightly outside threshold, inspect mesh asymmetry before relaxing limits.

## 5) Logging discipline

- Always run through logging wrapper so outputs are captured in `docs/testing/logs/*.log` and indexed in `docs/testing/test-results.md`.
- When reporting back, include:
  - chunk ID
  - pass/fail
  - log filename
  - one-line reason for any failure

# Human Port Calibration Checklist

Use this checklist before trusting birdcage lumped-port S-parameter results for tuning decisions.

## 1) Port placement realism

- [ ] Confirm each port tag is physically located between the intended leg pair.
- [ ] Confirm port orientation (`positive_tag` -> `negative_tag`) matches intended feed direction.
- [ ] Confirm port regions are not overlapping conductor bulk cells or phantom cells unexpectedly.
- [ ] In ParaView, verify all port faces are on comparable geometric areas (no accidental tiny/huge ports).

## 2) Termination assumptions

- [ ] Document the current termination model used for non-driven ports (ideal/simple termination in solver).
- [ ] Confirm this assumption is acceptable for the comparison target (bench fixture, VNA setup, or reference model).
- [ ] If comparing against measured data, note any baluns/cables/matching networks not represented in the FEM model.

## 3) Z0 and normalization choices

- [ ] Confirm reference impedance (`Z0`) value used during S-parameter normalization (default is usually 50 Ω).
- [ ] Confirm all ports use the same `Z0` unless you intentionally model mixed impedances.
- [ ] Recompute or rerun with alternate `Z0` (for example 50 Ω vs 75 Ω) when interpreting mismatch-sensitive trends.
- [ ] Record port ordering used in exported Touchstone metadata and keep it consistent across runs.

## 4) Comparison against known or bench measurements

- [ ] Compare simulated `S11` resonance region to known coil behavior or prior validated runs.
- [ ] Compare key coupling terms (for example adjacent-leg `S21`) for expected order of magnitude.
- [ ] Check if measured asymmetries can be explained by known hardware asymmetry (cable routing, connectors, loading).
- [ ] Log frequency grid, mesh resolution, and drive setup for each comparison run.

## Quick triage for suspicious S-parameter results

If results look nonphysical (for example `|Sii| > 1.05` over broad bands, unstable phase jumps, or unrealistic symmetry breaks), run this checklist:

1. Verify port tags exist and map to the intended physical regions.
2. Verify driven-port index and passive-port terminations for each sweep case.
3. Verify `Z0` normalization and port order in exported `.sNp` headers.
4. Re-run with a slightly finer mesh and compare trend stability (not just absolute values).
5. Cross-check one or two points against an independent expectation (analytic estimate, prior run, or bench datapoint).

## Minimum run notes to keep with each calibration attempt

- Geometry/mesh identifier and date
- Port ordering and orientation notes
- Termination model summary
- `Z0` used for normalization
- Frequency range and point count
- Any known caveats (fixtures, omitted components, solver limitations)

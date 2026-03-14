# MRI Birdcage + Phantom v1 Acceptance Checklist

This checklist defines minimum evidence required to declare **v1 achieved** for the coil + phantom workflow.

Status guidance:
- Mark item **PASS** only when linked logs/artifacts are present and criteria are met.
- Mark **BLOCKED** when a required prerequisite test cannot run (for example, environment unavailable).
- Keep this file aligned with `ROADMAP.md` and `docs/testing/pending-tests.md`.

---

## 1) Geometry validity

### Required chunk evidence
- B1 — Parametric birdcage geometry generator
- B2 — Port-face geometry robustness checks
- B3 — Phantom placement presets
- B4 — Air-box and boundary sizing heuristics
- B5 — Region-specific mesh resolution policy
- B6 — Geometry sanity report utility

### Acceptance criteria
- Human logs show required mesh tags are present with non-zero counts (`conductor`, `air`, `phantom`, required `port_*` tags).
- Human logs show no overlap/placement failures for centered and off-center phantom presets.
- Human logs include domain sizing diagnostics and no unhandled undersized-domain failures.
- Human logs include geometry sanity report output with expected sections and actionable warnings.

### Required artifacts/logs
- `docs/testing/test-results.md` rows for chunks B1..B6 with successful exit.
- Corresponding `docs/testing/logs/*_B1.log` .. `*_B6.log` entries.

---

## 2) Field plausibility and solve-health sanity

### Required chunk evidence
- C1 — Time-harmonic API hardening
- C2 — Phantom material model expansion
- C3 — Boundary-condition option set
- C4 — Interface-aware field extraction reliability
- C5 — Energy/consistency diagnostics
- C6 — Convergence/conditioning diagnostics
- E3 — Quick-look phantom metrics report

### Acceptance criteria
- Logs report finite, non-trivial phantom field metrics (no NaN/Inf in reported diagnostics).
- Material presets and frequency hooks execute without invalid-parameter failures.
- Boundary-condition selection is explicit and reflected in diagnostics.
- Consistency metrics and solve-health summaries are emitted and numerically finite.
- Quick-look report emits status (`OK` or `WARN`) with finite `|E|` and `|B|` summaries.

### Required artifacts/logs
- `docs/testing/test-results.md` rows for C1..C6 and E3 with successful exit.
- Corresponding `docs/testing/logs/*_C1.log` .. `*_C6.log` and `*_E3.log` files.
- Quick-look artifacts from example run (`mri_coil_phantom_quicklook.json`, `mri_coil_phantom_quicklook.md`).

---

## 3) S-parameter sanity and port-model credibility

### Required chunk evidence
- D1 — Calibration checklist to executable checks bridge
- D2 — Multi-port drive/termination consistency checks
- D3 — S-matrix reciprocity/passivity sanity metrics
- D4 — Frequency sweep orchestration utility
- D5 — Touchstone metadata completeness + parser cross-check
- D6 — Port-orientation sensitivity tests

### Acceptance criteria
- Port ordering/orientation/area checks are enforced and diagnostics reference calibration checklist.
- Logs show deterministic driven/passive port setup diagnostics.
- Reciprocity/passivity sanity metrics are emitted; suspicious cases produce warnings (not silent failure).
- Sweep planner produces deterministic, ordered frequency grids.
- Touchstone export includes metadata (`port_order`, `z0_ohm`, `frequency_points_hz`) and loader cross-check catches mismatches.
- Orientation sensitivity logs confirm expected sign-convention behavior when a port is flipped.

### Required artifacts/logs
- `docs/testing/test-results.md` rows for D1..D6 (successful exits for acceptance; blocked chunks must be resolved before signoff).
- Corresponding `docs/testing/logs/*_D1.log` .. `*_D6.log` files.
- At least one reviewed Touchstone output (`.sNp`) with metadata lines preserved.

---

## 4) Reproducibility and output bundle completeness

### Required chunk evidence
- E1 — Harden MRI example CLI/config
- E2 — Reproducible output bundle manifest
- E4 — Scenario presets (debug/dev/benchmark-lite)

### Acceptance criteria
- Example CLI help exposes expected control flags and exits cleanly.
- Manifest file is produced and contains at least: commit hash, key parameters, and artifact list.
- Preset selection is explicit in logs and effective parameter choices are printed.
- A debug preset run completes and produces expected output references.

### Required artifacts/logs
- `docs/testing/test-results.md` rows for E1, E2, and E4 with successful exit.
- Corresponding `docs/testing/logs/*_E1.log`, `*_E2.log`, `*_E4.log` files.
- Manifest artifact (`mri_coil_phantom_manifest.json`).

---

## 5) Signoff gate (v1 achieved)

Mark **v1 achieved** only when all statements below are true:

- [ ] No required chunk in sections 1-4 is in `🚫 BLOCKED` state.
- [ ] Every required chunk has a human-run logged result with successful exit in `docs/testing/test-results.md`.
- [ ] Required artifact files are present and referenced by logs.
- [ ] Any WARN-level diagnostics are reviewed and accepted as non-blocking with notes recorded in `docs/testing/pending-tests.md`.

## v1 signoff record

- Decision: `PENDING`
- Date (UTC): `TBD`
- Reviewer: `TBD`
- Notes: `TBD`

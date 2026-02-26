# ROADMAP.md — FEM-EM Solver (MRI Birdcage + Gelled Saline Phantom)

## Mission
Build a **reliable MRI-oriented EM solver workflow** that can:
1. Generate a realistic birdcage-like coil + gelled saline phantom geometry
2. Solve for magnetic and electric fields in/around the phantom
3. Produce robust diagnostics and ParaView-friendly outputs
4. Generate credible lumped-port S-parameters for downstream tuning workflows

This roadmap is organized into **~30 manageable chunks** designed for gradual progress via regular intervals.

---

## Ground Rules (Keep These)

### 1) Quality bar
Use `examples/magnetostatics/01_straight_wire.py` and existing validated tests as style references:
- clear diagnostics
- explicit assumptions
- deterministic checks where practical
- reproducible commands

### 2) Definition of done for a chunk (VPS-safe)
A chunk is done when:
- code/docs are committed
- manual test instruction is documented via `scripts/testing/run_and_log.sh`
- chunk status is set appropriately (`🧪` until human verification if heavy)

### 3) Cron-safe execution policy
During regular autonomous intervals:
- implement code + docs
- avoid heavy FEM/mesh solves unless explicitly lightweight
- update `docs/testing/pending-tests.md`
- provide exact human test command + expected pass signal

### 4) Resource constraints for human test runs
- Prefer `mpiexec -n 2`
- Keep memory under ~4GB for routine checks
- Keep problem sizes coarse for frequent iteration

### 5) Status legend
- ⬜ Not started
- 🟡 In progress
- 🧪 AWAITING-HUMAN-TEST
- ✅ Complete
- 🚫 Blocked

---

## Track A — Baseline Stabilization & Unblocks

### A1 — Resolve C2 symmetry metric strategy (sampling vs tolerance)
**Status:** ⬜
**Goal:** Unblock current C2 failure by separating physical asymmetry from numerical artifacts.

**Agent tasks:**
- Review `tests/validation/test_coil_phantom_bfield_metrics.py` sampling points and symmetry metric.
- Add interface-distance-aware sampling offsets to avoid boundary-induced spikes.
- Split metric output into absolute + relative components with clearer diagnostics.

**Manual test command:**
```bash
scripts/testing/run_and_log.sh A1 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/validation/test_coil_phantom_bfield_metrics.py -v'"
```

**Expected pass signal:**
- Symmetry test either passes with updated robust criterion, or fails with clear actionable diagnostics.

---

### A2 — Deterministic test tolerance policy
**Status:** ⬜
**Goal:** Reduce flaky fails by standardizing tolerances across solver/validation tests.

**Agent tasks:**
- Add shared tolerance constants/utilities.
- Refactor tests to consume shared tolerances.
- Document tolerance rationale in test comments/docs.

**Manual test command:**
```bash
scripts/testing/run_and_log.sh A2 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src python3 -m pytest tests/validation tests/solver -v -k tolerance'"
```

**Expected pass signal:**
- No ad-hoc tolerance literals in updated tests; tolerance behavior is consistent and documented.

---

### A3 — Lightweight smoke matrix for cron-safe confidence
**Status:** ⬜
**Goal:** Add a small test subset to detect regressions early without heavy runs.

**Agent tasks:**
- Create a curated smoke target (fast tests only).
- Add script alias (e.g., `run_tests.sh --smoke`).
- Ensure no heavy mesh/solve tasks included.

**Manual test command:**
```bash
scripts/testing/run_and_log.sh A3 "docker compose exec fem-em-solver bash -lc 'cd /workspace && ./run_tests.sh --smoke'"
```

**Expected pass signal:**
- Smoke command completes quickly and covers core module import + lightweight logic.

---

### A4 — Mesh-tag QA diagnostic hardening
**Status:** ⬜
**Goal:** Improve failure readability for missing/empty tags.

**Agent tasks:**
- Enrich `mesh_qa` outputs with expected vs actual counts and tag names.
- Add helper to print compact tag summary on failure.

**Manual test command:**
```bash
scripts/testing/run_and_log.sh A4 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/mesh/test_mesh_tag_integrity.py -v'"
```

**Expected pass signal:**
- Failures (or forced negative tests) show explicit tag-level diagnostics.

---

### A5 — Testing status dashboard section
**Status:** ⬜
**Goal:** Make pending/blocked/completed status obvious at a glance.

**Agent tasks:**
- Add/refresh summary table in `docs/testing/pending-tests.md`.
- Include chunk id, status, commit, last known log link.

**Manual test command:**
```bash
scripts/testing/run_and_log.sh A5 "docker compose exec fem-em-solver bash -lc 'cd /workspace && test -f docs/testing/pending-tests.md && echo OK'"
```

**Expected pass signal:**
- `OK` and dashboard section present with clear status columns.

---

## Track B — Geometry Realism for Birdcage + Phantom

### B1 — Parametric birdcage geometry generator
**Status:** ⬜
**Goal:** Introduce configurable birdcage geometry dimensions suitable for phased studies.

**Agent tasks:**
- Add parameters for ring radius, leg count, leg width, spacing, coil length.
- Keep defaults lightweight and meshing-safe.

**Manual test command:**
```bash
scripts/testing/run_and_log.sh B1 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/mesh/test_birdcage_port_tags.py -v'"
```

**Expected pass signal:**
- Tag integrity preserved under default param set.

---

### B2 — Port-face geometry robustness checks
**Status:** ⬜
**Goal:** Ensure port regions are valid, non-overlapping, and physically plausible.

**Agent tasks:**
- Add checks for minimum port area and separation.
- Add overlap detection between port and bulk conductor/phantom tags.

**Manual test command:**
```bash
scripts/testing/run_and_log.sh B2 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/mesh/test_birdcage_port_tags.py -v -k port'"
```

**Expected pass signal:**
- Port validation tests pass with explicit area/separation outputs.

---

### B3 — Phantom placement presets (centered/off-center)
**Status:** ⬜
**Goal:** Support realistic loaded scenarios beyond perfect symmetry.

**Agent tasks:**
- Add phantom placement presets.
- Validate no overlap and non-empty phantom tag in each preset.

**Manual test command:**
```bash
scripts/testing/run_and_log.sh B3 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/mesh/test_coil_phantom_mesh.py -v'"
```

**Expected pass signal:**
- Presets generate valid meshes and tag counts.

---

### B4 — Air-box and boundary sizing heuristics
**Status:** ⬜
**Goal:** Avoid boundary artifacts from too-tight domains.

**Agent tasks:**
- Add geometry heuristics for minimum air padding relative to coil size.
- Emit warning when user-provided domain is too small.

**Manual test command:**
```bash
scripts/testing/run_and_log.sh B4 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src python3 -m pytest tests/mesh -v -k domain'"
```

**Expected pass signal:**
- Domain sizing checks detect undersized setups predictably.

---

### B5 — Region-specific mesh resolution policy
**Status:** ⬜
**Goal:** Balance speed and fidelity by region (coil/phantom/air).

**Agent tasks:**
- Add per-region mesh size controls.
- Keep defaults cron-safe.

**Manual test command:**
```bash
scripts/testing/run_and_log.sh B5 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/mesh/test_mesh_tag_integrity.py -v'"
```

**Expected pass signal:**
- Mesh generation remains stable with region-specific sizing.

---

### B6 — Geometry sanity report utility
**Status:** ⬜
**Goal:** Generate one compact report of geometry/tag quality.

**Agent tasks:**
- Add utility/report with volume ratios + required tag counts + warnings.
- Wire into example/debug output.

**Manual test command:**
```bash
scripts/testing/run_and_log.sh B6 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src python3 -m pytest tests/mesh -v -k sanity'"
```

**Expected pass signal:**
- Report output appears with clear pass/warn sections.

---

## Track C — Frequency-Domain EM Solve Quality

### C1 — Time-harmonic API hardening
**Status:** ⬜
**Goal:** Clarify and enforce assumptions in the frequency-domain solver API.

**Agent tasks:**
- Add explicit parameter validation (frequency, units, material maps).
- Improve error messages for missing tag/material assignments.

**Manual test command:**
```bash
scripts/testing/run_and_log.sh C1 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/solver/test_time_harmonic_smoke.py -v'"
```

**Expected pass signal:**
- Smoke tests pass with improved validation coverage.

---

### C2 — Phantom material model expansion
**Status:** ⬜
**Goal:** Support practical gelled saline parameter presets for MRI-like frequencies.

**Agent tasks:**
- Add preset sets (low/mid/high conductivity variants).
- Add frequency-dependent helper hooks (MVP-level).

**Manual test command:**
```bash
scripts/testing/run_and_log.sh C2 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/materials/test_phantom_material_model.py -v'"
```

**Expected pass signal:**
- Material tests pass with finite derived terms across preset range.

---

### C3 — Boundary-condition option set
**Status:** ⬜
**Goal:** Compare and document BC choices impacting E/B stability.

**Agent tasks:**
- Add minimal BC option enum/config path.
- Add tests asserting BC selection is applied as intended.

**Manual test command:**
```bash
scripts/testing/run_and_log.sh C3 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src python3 -m pytest tests/solver -v -k boundary'"
```

**Expected pass signal:**
- BC variant tests pass; diagnostics include selected BC.

---

### C4 — Interface-aware field extraction reliability
**Status:** ⬜
**Goal:** Make E/B sampling robust near material boundaries.

**Agent tasks:**
- Add distance-from-interface sampling guardrails.
- Add fallback behavior for invalid sample cells.

**Manual test command:**
```bash
scripts/testing/run_and_log.sh C4 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/post/test_phantom_field_metrics.py -v'"
```

**Expected pass signal:**
- Metrics/export tests pass with reduced sampling pathologies.

---

### C5 — Energy/consistency diagnostics
**Status:** ⬜
**Goal:** Add solver-level consistency checks beyond finite/non-zero.

**Agent tasks:**
- Add diagnostic outputs for field norms and basic consistency indicators.
- Add threshold checks with clear warnings.

**Manual test command:**
```bash
scripts/testing/run_and_log.sh C5 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src python3 -m pytest tests/validation -v -k metrics'"
```

**Expected pass signal:**
- New diagnostics appear in logs and pass expected thresholds.

---

### C6 — Convergence/conditioning diagnostics
**Status:** ⬜
**Goal:** Surface conditioning issues early for human tuning.

**Agent tasks:**
- Add optional iterative diagnostics (residual trend summaries).
- Log solve health indicators in examples/tests.

**Manual test command:**
```bash
scripts/testing/run_and_log.sh C6 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src python3 -m pytest tests/solver -v -k convergence'"
```

**Expected pass signal:**
- Logs include stable residual/conditioning summaries where available.

---

## Track D — Port Model & S-Parameter Credibility

### D1 — Calibration checklist to executable checks bridge
**Status:** ⬜
**Goal:** Convert manual calibration guidance into partial automated assertions.

**Agent tasks:**
- Add helper checks for port ordering, area consistency, and orientation metadata presence.
- Link checks to checklist docs.

**Manual test command:**
```bash
scripts/testing/run_and_log.sh D1 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src python3 -m pytest tests/ports/test_port_definition.py -v'"
```

**Expected pass signal:**
- Port definition tests cover checklist-derived invariants.

---

### D2 — Multi-port drive/termination consistency checks
**Status:** ⬜
**Goal:** Ensure driven/passive configuration is deterministic and traceable.

**Agent tasks:**
- Add explicit validation on driven index and passive terminations.
- Emit consistent per-port solve context in diagnostics.

**Manual test command:**
```bash
scripts/testing/run_and_log.sh D2 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src python3 -m pytest tests/solver/test_single_port_excitation.py -v'"
```

**Expected pass signal:**
- Tests pass and logs clearly indicate drive/termination setup.

---

### D3 — S-matrix reciprocity/passivity sanity metrics
**Status:** ⬜
**Goal:** Add first-line physical sanity checks for computed S-matrices.

**Agent tasks:**
- Add reciprocity delta metrics (`Sij` vs `Sji`).
- Add passivity heuristics / warnings for suspicious norms.

**Manual test command:**
```bash
scripts/testing/run_and_log.sh D3 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src python3 -m pytest tests/ports/test_sparameter_assembly.py -v'"
```

**Expected pass signal:**
- Diagnostics include reciprocity/passivity summary values.

---

### D4 — Frequency sweep orchestration utility
**Status:** ⬜
**Goal:** Support practical coarse-to-refined frequency sweeps.

**Agent tasks:**
- Add simple sweep planner utility.
- Store sweep metadata (freqs, step policy, Z0, port order).

**Manual test command:**
```bash
scripts/testing/run_and_log.sh D4 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src python3 -m pytest tests/ports -v -k sweep'"
```

**Expected pass signal:**
- Sweep utility tests pass with deterministic frequency grids.

---

### D5 — Touchstone metadata completeness + parser cross-check
**Status:** ⬜
**Goal:** Make exported files robust for external tooling pipelines.

**Agent tasks:**
- Ensure header contains port ordering, Z0, frequency points, generation timestamp.
- Add stronger roundtrip parser assertions.

**Manual test command:**
```bash
scripts/testing/run_and_log.sh D5 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src python3 -m pytest tests/io/test_touchstone_export.py -v'"
```

**Expected pass signal:**
- Roundtrip loader confirms metadata integrity + value consistency.

---

### D6 — Port-orientation sensitivity tests
**Status:** ⬜
**Goal:** Quantify effect of orientation mistakes and catch sign convention regressions.

**Agent tasks:**
- Add tests comparing normal vs flipped orientation cases.
- Document expected qualitative effects in docs.

**Manual test command:**
```bash
scripts/testing/run_and_log.sh D6 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src python3 -m pytest tests/ports -v -k orientation'"
```

**Expected pass signal:**
- Tests confirm orientation changes produce expected directional/sign effects.

---

## Track E — End-to-End MRI Workflow & Outputs

### E1 — Harden MRI example CLI/config
**Status:** ⬜
**Goal:** Make `examples/mri/01_coil_phantom_fields.py` configurable and reproducible.

**Agent tasks:**
- Add basic CLI args (frequency, resolution preset, output dir).
- Preserve sensible defaults.

**Manual test command:**
```bash
scripts/testing/run_and_log.sh E1 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 examples/mri/01_coil_phantom_fields.py --help'"
```

**Expected pass signal:**
- CLI help shows options without errors.

---

### E2 — Reproducible output bundle manifest
**Status:** ⬜
**Goal:** Package outputs + metadata for easier review/comparison.

**Agent tasks:**
- Emit manifest JSON with commit hash, parameters, file list.
- Include key metrics references.

**Manual test command:**
```bash
scripts/testing/run_and_log.sh E2 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 examples/mri/01_coil_phantom_fields.py'"
```

**Expected pass signal:**
- Output directory includes manifest file with parameter + artifact entries.

---

### E3 — Quick-look phantom metrics report
**Status:** ⬜
**Goal:** Provide immediate human-readable summary after each run.

**Agent tasks:**
- Add compact report generation (console + optional markdown/json).
- Include |E| and |B| summaries and warning flags.

**Manual test command:**
```bash
scripts/testing/run_and_log.sh E3 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 examples/mri/01_coil_phantom_fields.py'"
```

**Expected pass signal:**
- Log includes clear quick-look section with finite metrics.

---

### E4 — Scenario presets (debug/dev/benchmark-lite)
**Status:** ⬜
**Goal:** Standardize run profiles for iterative work.

**Agent tasks:**
- Add named presets controlling mesh/sweep/detail levels.
- Ensure debug preset remains lightweight.

**Manual test command:**
```bash
scripts/testing/run_and_log.sh E4 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 examples/mri/01_coil_phantom_fields.py --preset debug'"
```

**Expected pass signal:**
- Preset runs complete; logs identify selected preset and effective parameters.

---

## Track F — Human-in-the-Loop Testing Operations

### F1 — Expand run-and-log metadata
**Status:** ⬜
**Goal:** Improve traceability of manual verification runs.

**Agent tasks:**
- Extend logging index with chunk, commit, elapsed time, key env metadata.
- Keep backward compatibility with existing logs.

**Manual test command:**
```bash
scripts/testing/run_and_log.sh F1 "docker compose exec fem-em-solver bash -lc 'cd /workspace && ./run_tests.sh --list'"
```

**Expected pass signal:**
- `test-results.md` entries include enhanced metadata fields.

---

### F2 — Guided pending-test queue helper
**Status:** ⬜
**Goal:** Surface the next most useful human test to run.

**Agent tasks:**
- Add helper that reads `pending-tests` and outputs prioritized next commands.
- Prioritize unblockers and dependencies.

**Manual test command:**
```bash
scripts/testing/run_and_log.sh F2 "docker compose exec fem-em-solver bash -lc 'cd /workspace && ./run_tests.sh --list'"
```

**Expected pass signal:**
- Output includes a clear “recommended next test” order.

---

### F3 — Define v1 milestone acceptance checklist
**Status:** ⬜
**Goal:** Establish explicit criteria for “MRI birdcage + phantom v1 achieved”.

**Agent tasks:**
- Add milestone checklist doc (geometry validity, field plausibility, S-parameter sanity, reproducibility).
- Link required logs/artifacts for signoff.

**Manual test command:**
```bash
scripts/testing/run_and_log.sh F3 "docker compose exec fem-em-solver bash -lc 'cd /workspace && test -f docs/testing/v1_acceptance_checklist.md && echo OK'"
```

**Expected pass signal:**
- `OK` and checklist includes measurable acceptance criteria.

---

## Immediate execution guidance
1. Start with **A1** (unblock C2) before adding new physics complexity.
2. Then progress **A2 → A3 → B1/B2** for stable foundations.
3. Move through Tracks C and D in parallel only when baseline tests are stable.
4. Keep E/F tracks active to maintain usable outputs and human verification flow.

---

## Current priority
**Next chunk: A1 — Resolve C2 symmetry metric strategy.**

Reason: Existing known blocker should be resolved before further broad expansion.
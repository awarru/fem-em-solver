# ROADMAP.md — FEM-EM Solver (Agent-Executable Chunks)

## Mission
Build a **reliable MRI-oriented EM solver workflow** that can:
1. Mesh a coil + phantom configuration (gelled saline phantom inside an MRI coil)
2. Solve for magnetic and electric fields in the phantom region
3. Export results that are easy to inspect and validate
4. Generate lumped-port S-parameters for downstream circuit tuning workflows

This roadmap is optimized for an automated coding agent working in ~30-minute chunks.

---

## Ground Rules (Read First)

### 1) Follow the straight-wire example standard
Use `examples/magnetostatics/01_straight_wire.py` as the quality template:
- Clear diagnostics and printouts
- Robust point-evaluation workflow (cell lookup before `eval`)
- Practical ParaView outputs
- Explicit numerical sanity checks
- Reproducible test commands

### 2) Definition of Done for every chunk (VPS-safe mode)
A chunk is complete when:
- Code changes are committed
- A **manual test note** is added with exact command + expected signal
- Agent does **not** run heavy FEM tests during cron cycles
- Chunk is marked `🧪 AWAITING-HUMAN-TEST` until you run and report results

### 3) When blocked
If blocked >30 min:
- Commit partial work to a branch-safe state (or stash in notes)
- Mark chunk as BLOCKED with reason
- Move to next chunk

### 4) Keep chunks practical
Small, testable, and likely to succeed in one agent run.

### 5) MPI and resource constraints
When you run tests manually, use constrained resources:
- **Use `mpiexec -n 2`** (2 cores)
- **Memory limit:** 4GB RAM maximum
- **Problem sizes:** Keep meshes coarse (resolution ≥ 0.01m, cell counts < 50k)

### 6) Cron execution policy (critical)
During cron cycles, agent should:
- Implement code + docs only
- **Not run heavy FEM/mesh/solver tests**
- Write test instructions to `docs/testing/pending-tests.md`
- Provide a logging command using `scripts/testing/run_and_log.sh`
- Mark chunk `🧪 AWAITING-HUMAN-TEST` with command and expected signal
- Wait for your reported results before marking fully complete

---

## Status Legend
- ⬜ Not started
- 🟡 In progress / partial
- 🧪 AWAITING-HUMAN-TEST
- ✅ Complete
- 🚫 Blocked (needs human)

---

## Phase A — Stabilize Core Infrastructure (before new physics)

### ✅ A0 — Baseline straight-wire quality reference
Keep straight-wire behavior as golden standard for diagnostics/evaluation/export.

### ✅ A1 — Add shared field-evaluation utility (2026-02-20, c01afde)
**Goal:** Remove duplicated fragile point-evaluation logic from examples/tests.

**Agent tasks:**
- Add utility module (e.g. `post/evaluation.py`) that:
  - locates containing cells for points
  - evaluates vector field robustly in parallel
  - returns values + mask for invalid points
- Refactor straight-wire example to use this utility (no behavior change expected).

**Human test command (run manually):**
```bash
cd ~/Development/fem-em-solver/docker
docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 examples/magnetostatics/01_straight_wire.py'
```

**Expected signal:**
- Example runs end-to-end
- Same style diagnostics printed
- Plot/output files still generated

**Human verification (YOU):**
- Quickly open output in ParaView and confirm visualization still works.

### ✅ A2 — Add shared ParaView export helper (2026-02-20, ff6bc8c)
**Goal:** Standardize exports (mesh tags + fields together) across examples.

**Agent tasks:**
- Implement/reinforce helper for combined XDMF output with tags + selected fields.
- Use helper in straight-wire and one coil example.

**Human test command (run manually):**
```bash
docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 examples/magnetostatics/01_straight_wire.py && ls -1 paraview_output'
```

**Expected signal:**
- Combined output file exists and can be opened
- No regression in existing exports

**Human verification (YOU):**
- Confirm you can threshold phantom/air/coil tags in ParaView.

---

## Phase B — Coil + Phantom Geometry Pipeline

### ✅ B1 — Add phantom-ready mesh generator (coil + cylindrical phantom + air box) (2026-02-20, 9de1bf8)
**Goal:** Create a robust geometry function for a first MRI-like test model.

**Constraints:**
- Mesh must generate with ≤ 4GB RAM
- Target: < 50k cells total (use resolution ≥ 0.01m)
- Must work with mpiexec -n 2

**Agent tasks:**
- Add mesh generator function (new method) for:
  - two-loop (or two-torus) coil conductors
  - cylindrical phantom volume centered in coil region
  - surrounding air/domain volume
- Ensure clean physical tags for at least:
  - coil_1, coil_2, phantom, air
- Use coarse resolution (0.01-0.02m) to stay within memory limits

**Human test command (run manually):**
```bash
docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/mesh/test_coil_phantom_mesh.py -v'
```

**Expected signal:**
- Test confirms mesh generation and required tags exist

**Human verification (YOU):**
- Inspect mesh visually once; verify phantom is actually inside coil region.

### 🧪 B2 — Harden mesh QA checks
**Goal:** Catch bad geometry early.

**Agent tasks:**
- Add test helpers that assert:
  - nonzero cell count per required tag
  - phantom volume is not empty
  - coil and phantom volumes are distinct
- Add quick mesh summary print utility (counts by tag)

**Human test command (run manually):**
```bash
docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/mesh/test_mesh_tag_integrity.py -v'
```

**Expected signal:**
- Deterministic pass/fail with clear failure messages

---

## Phase C — Magnetostatics on Coil + Phantom

### ⬜ C1 — Solve B-field on coil+phantom model
**Goal:** Reliable magnetostatic solve with source current in coil subdomains only.

**Agent tasks:**
- Build a test that:
  - assigns current to coil tags only
  - solves for A and computes B
  - checks finite + nontrivial B in phantom sample points
- Keep gauge stabilization behavior explicit and configurable

**Human test command (run manually):**
```bash
docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/solver/test_coil_phantom_magnetostatics.py -v'
```

**Expected signal:**
- Test passes
- B-field in phantom is finite and not near-zero everywhere

**Human verification (YOU):**
- Spot-check one centerline/profile plot looks physically plausible.

### ⬜ C2 — Add sanity validation metrics
**Goal:** Prevent silently wrong fields.

**Agent tasks:**
- Add simple metrics in tests/examples:
  - min/max/mean |B| in phantom
  - centerline smoothness check
  - optional symmetry check for symmetric setups

**Human test command (run manually):**
```bash
docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/validation/test_coil_phantom_bfield_metrics.py -v'
```

**Expected signal:**
- Threshold-based checks pass with printed metric values

---

## Phase D — Time-Harmonic E-field in Phantom (first practical version)

### ⬜ D1 — Introduce minimal frequency-domain solve scaffold
**Goal:** Add a narrow, testable path to compute E-field in phantom.

**Agent tasks:**
- Add a minimal time-harmonic module/API path (keep scope tight)
- Accept frequency + material properties
- Return an E-field object usable for sampling/export

**Human test command (run manually):**
```bash
docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/solver/test_time_harmonic_smoke.py -v'
```

**Expected signal:**
- Smoke test runs and returns finite E-field values

**Human verification (YOU):**
- If numerics look unstable, tune formulation choices manually.

### ⬜ D2 — Add gelled saline phantom material model (MVP)
**Goal:** Support phantom electrical properties needed for E-field estimates.

**Agent tasks:**
- Add material container for phantom with at least:
  - conductivity `sigma`
  - relative permittivity `epsilon_r`
  - frequency parameter
- Wire model into time-harmonic path for phantom-tagged cells

**Human test command (run manually):**
```bash
docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/materials/test_phantom_material_model.py -v'
```

**Expected signal:**
- Material assignment tests pass and are used in solve pipeline

### ⬜ D3 — E and B field extraction inside phantom
**Goal:** Compute both E and B metrics specifically in phantom region.

**Agent tasks:**
- Add post-processing helpers to sample/export inside phantom tag
- Compute summary stats for |E| and |B| in phantom

**Human test command (run manually):**
```bash
docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/post/test_phantom_field_metrics.py -v'
```

**Expected signal:**
- Tests confirm finite E/B values in phantom and export files created

**Human verification (YOU):**
- Validate magnitude order-of-growth vs expectation for your coil current/frequency.

---

## Phase E — Lumped Ports + S-Parameter Pipeline (Birdcage-oriented)

### ⬜ E1 — Define lumped port data model and tagging contract
**Goal:** Introduce explicit, testable representation of lumped ports between birdcage legs.

**Agent tasks:**
- Add port schema/data class (e.g. `PortDefinition`) with at least:
  - `port_id`, `positive_tag`, `negative_tag`
  - feed direction / orientation metadata
  - optional reference impedance (`Z0`, default 50Ω)
- Define mesh tagging contract for port faces/edges (documented constants).
- Add validation helper that checks required port tags exist.

**Human test command (run manually):**
```bash
docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/ports/test_port_definition.py -v'
```

**Expected signal:**
- Port definition and tag validation tests pass.

### ⬜ E2 — Add minimal birdcage-like test geometry with port tags
**Goal:** Create a lightweight geometry fixture that supports multiple lumped ports without blowing memory.

**Agent tasks:**
- Add a coarse "birdcage-like" geometry fixture (small cell count, not full production detail):
  - ring + simplified legs
  - explicit port regions between adjacent legs
- Ensure tags exist for: conductor, air, phantom (optional coarse), and each port.

**Constraints:**
- Keep mesh under 50k cells and runnable under 4GB RAM.

**Human test command (run manually):**
```bash
docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/mesh/test_birdcage_port_tags.py -v'
```

**Expected signal:**
- Test verifies all expected port tags and core regions exist.

**Human verification (YOU):**
- Visual check that each port sits between intended leg pair.

### ⬜ E3 — Implement port excitation hook (single-port solve)
**Goal:** Excite one lumped port at a time in frequency-domain solve.

**Agent tasks:**
- Add API to run a single-port excitation case:
  - one driven port, others terminated (initial simple termination model)
- Return per-port voltage/current estimates needed for S-parameter assembly.
- Keep implementation simple and explicit; prioritize deterministic behavior over sophistication.

**Human test command (run manually):**
```bash
docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/solver/test_single_port_excitation.py -v'
```

**Expected signal:**
- Test passes with finite voltage/current results for driven + passive ports.

### ⬜ E4 — Build N-port sweep and S-parameter assembly
**Goal:** Generate an S-matrix by sweeping driven port index.

**Agent tasks:**
- Add routine: for N ports, run N excitations and assemble NxN S-matrix.
- Include basic checks:
  - matrix shape correct
  - finite complex entries
  - diagonal reflection terms present
- Keep frequency list short/coarse for resource limits.

**Human test command (run manually):**
```bash
docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/ports/test_sparameter_assembly.py -v'
```

**Expected signal:**
- S-matrix tests pass and produce deterministic dimensions/outputs.

**Human verification (YOU):**
- Sanity-check S11/S21 trends against expectations for your geometry.

### ⬜ E5 — Export S-parameters for external circuit tuning workflow
**Goal:** Make solver output directly usable in separate tuning/circuit tools.

**Agent tasks:**
- Add export to at least Touchstone `.sNp` (and optional CSV companion).
- Include metadata in filename/header (frequency points, Z0, port ordering).
- Add tiny loader/roundtrip test to ensure exported file is parseable.

**Human test command (run manually):**
```bash
docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/io/test_touchstone_export.py -v'
```

**Expected signal:**
- `.sNp` file written and reloaded successfully in test.

### ⬜ E6 — Add "human calibration" checklist for port model assumptions
**Goal:** Explicitly call out assumptions agent cannot reliably tune alone.

**Agent tasks:**
- Add doc section/file describing what human should calibrate:
  - port placement realism
  - termination assumptions
  - Z0 choices and normalization
  - comparison vs known/bench measurements
- Add quick checklist for interpreting suspicious S-parameter results.

**Human test command (run manually):**
```bash
docker compose exec fem-em-solver bash -lc 'cd /workspace && test -f docs/ports/human_port_calibration_checklist.md && echo OK'
```

**Expected signal:**
- Checklist file exists and is actionable.

---

## Phase F — End-to-End Example + Documentation

### ⬜ F1 — New example: MRI coil with gelled saline phantom
**Goal:** One runnable script demonstrating end-to-end workflow.

**Agent tasks:**
- Add example script (e.g. `examples/mri/01_coil_phantom_fields.py`) that:
  - builds mesh
  - applies materials
  - solves fields
  - prints phantom metrics
  - exports ParaView files

**Human test command (run manually):**
```bash
docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 examples/mri/01_coil_phantom_fields.py'
```

**Expected signal:**
- Script finishes and writes output files
- Console shows phantom E/B summaries

### ⬜ F2 — Add “human test checklist” doc
**Goal:** Make it obvious what you should manually validate.

**Agent tasks:**
- Create `docs/testing_manual_checklist.md` with:
  - visual mesh checks
  - field sanity checks
  - known failure symptoms
  - suggested parameter tweaks when unstable

**Human test command (run manually):**
```bash
docker compose exec fem-em-solver bash -lc 'cd /workspace && test -f docs/testing_manual_checklist.md && echo OK'
```

**Expected signal:**
- Checklist file exists and is readable

---

## What You Should Personally Test/Fix (Human-in-the-loop)

These are intentionally marked as human checkpoints because agent debugging is limited:

1. **Mesh realism checks**
   - Coil physically surrounds phantom as intended
   - No accidental overlap/void artifacts
2. **Field plausibility checks**
   - Phantom |B| and |E| magnitudes in expected range
   - No obvious nonphysical spikes from bad sampling/BCs
3. **Stability tuning**
   - Gauge penalty, solver options, mesh resolution tradeoffs
4. **Physics assumptions**
   - Material values and frequency choice match your MRI use case

---

## Agent Execution Template (for cron job prompt)
For each run:
1. Find first ⬜ chunk (or 🧪 AWAITING-HUMAN-TEST chunk if you provided results)
2. Implement only that chunk
3. Do **not** run heavy FEM tests in cron mode
4. Append/update `docs/testing/pending-tests.md` with:
   - chunk id
   - manual test command (wrapped with `scripts/testing/run_and_log.sh`)
   - expected signal
   - files changed / commit hash
5. Human runs tests; results auto-log to:
   - `docs/testing/test-results.md` (index)
   - `docs/testing/logs/*.log` (full output)
6. Mark chunk as `🧪 AWAITING-HUMAN-TEST`
7. After you report test results, agent updates chunk to ✅ or 🚫 BLOCKED

---

## Immediate Next Chunk
**A2 — Add shared ParaView export helper**

Reason: with evaluation logic centralized, standardizing combined tag+field export is the next highest-leverage infrastructure step for reliable visual validation.

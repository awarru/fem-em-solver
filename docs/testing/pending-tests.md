# Pending Manual Tests

This file is updated by cron-driven coding runs in VPS-safe mode.

For each pending chunk, agent should append using this exact structure:
- Chunk: `<ID> — <title>`
- Status: `🧪 AWAITING-HUMAN-TEST`
- Commit: `<full-hash>`
- Files changed: bullet list
- Manual test command: **Exact command using** `scripts/testing/run_and_log.sh`
- Expected pass signal: bullet list of concrete `PASSED` lines or diagnostics
- Notes/blockers: concise line or `none`

Example:
```bash
# scripts/testing/run_and_log.sh A1 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 examples/magnetostatics/01_straight_wire.py'"
```

Logs are automatically written to:
- `docs/testing/test-results.md` (summary index)
- `docs/testing/logs/*.log` (full output)

Single-command entrypoint from repo root:
```bash
./run_tests.sh
```

Optional helpers:
```bash
./run_tests.sh --list
./run_tests.sh --chunk E3
```


## Testing Status Dashboard

| Chunk | Status | Commit | Last known log |
|---|---|---|---|
| C1 | ✅ COMPLETE | `09eb248f6e5ee161234d8a799692c75a63262efb` | `docs/testing/logs/20260226T164233Z_C1.log` |
| C2 | 🚫 BLOCKED | `7ac10f166e283ff7b6f15e20323b6402a4a49d65` | `docs/testing/logs/20260226T164235Z_C2.log` |
| A1 | 🧪 AWAITING-HUMAN-TEST | `79e5cb22abfb2ed757cd30937d6a4d97e5363b29` | none |
| A2 | 🧪 AWAITING-HUMAN-TEST | `529cc557998f51e48025a7fef4323cc54c259a2d` | none |
| A3 | 🧪 AWAITING-HUMAN-TEST | `9a61957e79936c9588d15805cfec10509afb76f3` | none |
| A4 | 🧪 AWAITING-HUMAN-TEST | `7c9b2c49cceb5f1035da23503e567ca242f6f821` | none |
| A5 | 🧪 AWAITING-HUMAN-TEST | `82d2b6e4a0bf00bc521d3751d727e79dbdb6208e` | none |

---

## C1 — Solve B-field on coil+phantom model (✅ COMPLETE)

- Commit: `09eb248f6e5ee161234d8a799692c75a63262efb`
- Files changed:
  - `ROADMAP.md`
  - `src/fem_em_solver/core/solvers.py`
  - `tests/solver/test_coil_phantom_magnetostatics.py`
- Manual test command:
  ```bash
  scripts/testing/run_and_log.sh C1 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/solver/test_coil_phantom_magnetostatics.py -v'"
  ```
- Expected pass signal:
  - Pytest reports `tests/solver/test_coil_phantom_magnetostatics.py::test_coil_phantom_magnetostatics_bfield_is_finite_and_nontrivial_in_phantom PASSED`
  - Test output prints `coil+phantom B-field diagnostics` with finite non-zero phantom `|B|` min/max/mean
- Human test log:
  - `docs/testing/logs/20260223T022337Z_C1.log` (exit 0)

---

## C2 — Add sanity validation metrics (🚫 BLOCKED)

- Commit: `7ac10f166e283ff7b6f15e20323b6402a4a49d65`
- Files changed:
  - `ROADMAP.md`
  - `docs/testing/pending-tests.md`
  - `tests/validation/test_coil_phantom_bfield_metrics.py`
- Manual test command:
  ```bash
  scripts/testing/run_and_log.sh C2 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/validation/test_coil_phantom_bfield_metrics.py -v'"
  ```
- Expected pass signal:
  - Pytest reports `tests/validation/test_coil_phantom_bfield_metrics.py::test_coil_phantom_bfield_metrics_are_finite_smooth_and_symmetric PASSED`
  - Test output prints `coil+phantom B-field metrics` with finite `|B|` min/max/mean and bounded smoothness/symmetry diagnostics
- Human test log:
  - `docs/testing/logs/20260223T022341Z_C2.log` (exit 1)
  - Failure: `Symmetry sanity check failed for ±x phantom points; max relative |B| mismatch=0.322` (limit `< 0.30`)

---

## D1 — Introduce minimal frequency-domain solve scaffold (🧪 AWAITING-HUMAN-TEST)

- Commit: `1b2186e0d87e1db87503ce273193fd94635fcde3`
- Files changed:
  - `ROADMAP.md`
  - `src/fem_em_solver/__init__.py`
  - `src/fem_em_solver/core/__init__.py`
  - `src/fem_em_solver/core/time_harmonic.py`
  - `tests/solver/test_time_harmonic_smoke.py`
- Manual test command:
  ```bash
  scripts/testing/run_and_log.sh D1 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/solver/test_time_harmonic_smoke.py -v'"
  ```
- Expected pass signal:
  - Pytest reports `tests/solver/test_time_harmonic_smoke.py::test_time_harmonic_smoke_returns_finite_e_field_values PASSED`
  - Test output prints `time-harmonic smoke diagnostics` with finite non-zero `|E_imag|` min/max/mean

## D2 — Add gelled saline phantom material model (MVP) (🧪 AWAITING-HUMAN-TEST)

- Commit: `c03f461e73a70c7fc8dd83291c4cf6531bd5b1c6`
- Files changed:
  - `ROADMAP.md`
  - `docs/testing/pending-tests.md`
  - `src/fem_em_solver/__init__.py`
  - `src/fem_em_solver/core/__init__.py`
  - `src/fem_em_solver/core/time_harmonic.py`
  - `src/fem_em_solver/materials/__init__.py`
  - `src/fem_em_solver/materials/phantom.py`
  - `tests/materials/test_phantom_material_model.py`
- Manual test command:
  ```bash
  scripts/testing/run_and_log.sh D2 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/materials/test_phantom_material_model.py -v'"
  ```
- Expected pass signal:
  - Pytest reports `tests/materials/test_phantom_material_model.py::test_gelled_saline_material_container_frequency_term_is_finite PASSED`
  - Pytest reports `tests/materials/test_phantom_material_model.py::test_phantom_material_assignment_and_time_harmonic_pipeline_wiring PASSED`
  - Test output has no `ValueError` related to phantom tag assignment or frequency mismatch

## D3 — E and B field extraction inside phantom (🧪 AWAITING-HUMAN-TEST)

- Commit: `6cb701ca509fdc69d63feecb7d300c220476d4b9`
- Files changed:
  - `ROADMAP.md`
  - `src/fem_em_solver/post/__init__.py`
  - `src/fem_em_solver/post/phantom_fields.py`
  - `tests/post/test_phantom_field_metrics.py`
- Manual test command:
  ```bash
  scripts/testing/run_and_log.sh D3 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/post/test_phantom_field_metrics.py -v'"
  ```
- Expected pass signal:
  - Pytest reports `tests/post/test_phantom_field_metrics.py::test_phantom_field_metrics_and_exports_are_finite PASSED`
  - Test output prints `phantom E/B diagnostics` with finite non-zero phantom `|E|` and `|B|` min/max/mean
  - Log shows exported files `d3_test_phantom_E_samples.csv`, `d3_test_phantom_B_samples.csv`, and `d3_test_phantom_metrics.json`

## E1 — Define lumped port data model and tagging contract (🧪 AWAITING-HUMAN-TEST)

- Commit: `c4234cb73b889e52a8f76f9ee66f8a93d9dc7756`
- Files changed:
  - `ROADMAP.md`
  - `src/fem_em_solver/__init__.py`
  - `src/fem_em_solver/ports/__init__.py`
  - `src/fem_em_solver/ports/definitions.py`
  - `tests/ports/test_port_definition.py`
- Manual test command:
  ```bash
  scripts/testing/run_and_log.sh E1 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/ports/test_port_definition.py -v'"
  ```
- Expected pass signal:
  - Pytest reports `tests/ports/test_port_definition.py` tests all `PASSED`
  - Log has no `ValueError` except inside expected `pytest.raises(...)` checks for invalid inputs

## E2 — Add minimal birdcage-like test geometry with port tags (🧪 AWAITING-HUMAN-TEST)

- Commit: `fd85b2ba40d84eede3dcce8cfb46c3a1feac1879`
- Files changed:
  - `ROADMAP.md`
  - `src/fem_em_solver/io/mesh.py`
  - `tests/mesh/test_birdcage_port_tags.py`
- Manual test command:
  ```bash
  scripts/testing/run_and_log.sh E2 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/mesh/test_birdcage_port_tags.py -v'"
  ```
- Expected pass signal:
  - Pytest reports `tests/mesh/test_birdcage_port_tags.py::test_birdcage_like_mesh_has_core_and_port_tags PASSED`
  - Test output prints `[birdcage-mesh]` tag summary including non-zero counts for `conductor`, `air`, `phantom`, and `port_P1`..`port_P4`

## E3 — Implement port excitation hook (single-port solve) (🧪 AWAITING-HUMAN-TEST)

- Commit: `1d6fdc1cec0fc742779601ba1d8df9d1caad365a`
- Files changed:
  - `ROADMAP.md`
  - `src/fem_em_solver/__init__.py`
  - `src/fem_em_solver/ports/__init__.py`
  - `src/fem_em_solver/ports/excitation.py`
  - `tests/solver/test_single_port_excitation.py`
- Manual test command:
  ```bash
  scripts/testing/run_and_log.sh E3 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/solver/test_single_port_excitation.py -v'"
  ```
- Expected pass signal:
  - Pytest reports `tests/solver/test_single_port_excitation.py::test_single_port_excitation_returns_finite_estimates PASSED`
  - Pytest reports `tests/solver/test_single_port_excitation.py::test_single_port_excitation_rejects_missing_required_tags PASSED`
  - Test output prints `single-port excitation diagnostics` with finite voltage/current estimates for driven and passive ports

## E4 — Build N-port sweep and S-parameter assembly (🧪 AWAITING-HUMAN-TEST)

- Commit: `593ca473f1cc53f47dce1bece8ea76cb17cc23b8`
- Files changed:
  - `ROADMAP.md`
  - `src/fem_em_solver/__init__.py`
  - `src/fem_em_solver/ports/__init__.py`
  - `src/fem_em_solver/ports/sparameters.py`
  - `tests/ports/test_sparameter_assembly.py`
- Manual test command:
  ```bash
  scripts/testing/run_and_log.sh E4 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/ports/test_sparameter_assembly.py -v'"
  ```
- Expected pass signal:
  - Pytest reports `tests/ports/test_sparameter_assembly.py::test_n_port_sweep_assembles_finite_matrix_with_expected_shape PASSED`
  - Pytest reports `tests/ports/test_sparameter_assembly.py::test_n_port_sweep_rejects_zero_incident_drive PASSED`
  - Test output prints `n-port S-parameter sweep diagnostics` with S-matrix shape `(3, 3)` and diagonal `S11/S22/S33` terms

---

- Chunk: E5 — Export S-parameters for external circuit tuning workflow
- Status: 🧪 AWAITING-HUMAN-TEST
- Commit: 531d82a7b852eee6af0f0340adf5a1df4c6a7f9a
- Files changed:
  - ROADMAP.md
  - src/fem_em_solver/__init__.py
  - src/fem_em_solver/ports/__init__.py
  - src/fem_em_solver/ports/touchstone.py
  - tests/io/test_touchstone_export.py
- Manual test command: scripts/testing/run_and_log.sh E5 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/io/test_touchstone_export.py -v'"
- Expected pass signal:
  - tests/io/test_touchstone_export.py::test_touchstone_export_and_roundtrip_loader PASSED
  - Output artifacts include one `.s2p` file and matching `.csv` companion in pytest tmp path
  - Touchstone header contains `! port_order: P1,P2`, `! frequency_points_hz: ...`, and `! z0_ohm: 50.000000`
- Notes/blockers: none

- Chunk: E6 — Add "human calibration" checklist for port model assumptions
- Status: 🧪 AWAITING-HUMAN-TEST
- Commit: b99dd6b8769347171bd99c9abc2c178241d6b192
- Files changed:
  - ROADMAP.md
  - docs/ports/human_port_calibration_checklist.md
- Manual test command: scripts/testing/run_and_log.sh E6 "docker compose exec fem-em-solver bash -lc 'cd /workspace && test -f docs/ports/human_port_calibration_checklist.md && echo OK'"
- Expected pass signal:
  - Log contains `OK`
  - Command exits with status `0`
- Notes/blockers: none

- Chunk: F1 — New example: MRI coil with gelled saline phantom
- Status: 🧪 AWAITING-HUMAN-TEST
- Commit: 800829bf61aab91f834ba55a5dd7940140264139
- Files changed:
  - ROADMAP.md
  - examples/mri/01_coil_phantom_fields.py
- Manual test command: scripts/testing/run_and_log.sh F1 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 examples/mri/01_coil_phantom_fields.py'"
- Expected pass signal:
  - Output contains `Example: MRI coil + gelled saline phantom fields`
  - Output contains `Phantom diagnostics:` plus finite `|E| min/max/mean` and `|B| min/max/mean` lines
  - Output lists `mri_coil_phantom_fields_combined.xdmf` and ends with `Example completed`
- Notes/blockers: none

- Chunk: F2 — Add “human test checklist” doc
- Status: 🧪 AWAITING-HUMAN-TEST
- Commit: a5a7350c779deade5ffeac1ab446baaf48e1bece
- Files changed:
  - ROADMAP.md
  - docs/testing/pending-tests.md
  - docs/testing_manual_checklist.md
- Manual test command: scripts/testing/run_and_log.sh F2 "docker compose exec fem-em-solver bash -lc 'cd /workspace && test -f docs/testing_manual_checklist.md && echo OK'"
- Expected pass signal:
  - Log contains `OK`
  - Command exits with status `0`
- Notes/blockers: awaiting human-run log in docs/testing/test-results.md (no new F2 result yet)

- Chunk: A1 — Resolve C2 symmetry metric strategy (sampling vs tolerance)
- Status: 🧪 AWAITING-HUMAN-TEST
- Commit: 79e5cb22abfb2ed757cd30937d6a4d97e5363b29
- Files changed:
  - ROADMAP.md
  - tests/validation/test_coil_phantom_bfield_metrics.py
  - docs/testing/pending-tests.md
- Manual test command: scripts/testing/run_and_log.sh A1 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/validation/test_coil_phantom_bfield_metrics.py -v'"
- Expected pass signal:
  - tests/validation/test_coil_phantom_bfield_metrics.py::test_coil_phantom_bfield_metrics_are_finite_smooth_and_symmetric PASSED
  - Output includes `symmetry probe setup:` with interface clearance and interior safe radius/half-height diagnostics
  - Output includes `symmetry mismatch diagnostics (±x pairs):` with both absolute and relative max/mean values and tolerances
- Notes/blockers: none

- Chunk: A2 — Deterministic test tolerance policy
- Status: 🧪 AWAITING-HUMAN-TEST
- Commit: 529cc557998f51e48025a7fef4323cc54c259a2d
- Files changed:
  - ROADMAP.md
  - docs/testing/tolerance-policy.md
  - tests/tolerances.py
  - tests/solver/test_coil_phantom_magnetostatics.py
  - tests/solver/test_cylinder.py
  - tests/solver/test_time_harmonic_smoke.py
  - tests/solver/test_tolerance_policy.py
  - tests/solver/test_two_cylinder.py
  - tests/validation/test_coil_phantom_bfield_metrics.py
  - tests/validation/test_tolerance_policy.py
- Manual test command: scripts/testing/run_and_log.sh A2 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src python3 -m pytest tests/validation tests/solver -v -k tolerance'"
- Expected pass signal:
  - tests/validation/test_tolerance_policy.py::test_validation_tolerance_policy_is_ordered_and_positive PASSED
  - tests/solver/test_tolerance_policy.py::test_solver_tolerance_policy_is_consistent PASSED
  - No failures mentioning undefined tolerance constants; updated solver/validation tests import thresholds from tests/tolerances.py
- Notes/blockers: none

- Chunk: A3 — Lightweight smoke matrix for cron-safe confidence
- Status: 🧪 AWAITING-HUMAN-TEST
- Commit: 9a61957e79936c9588d15805cfec10509afb76f3
- Files changed:
  - ROADMAP.md
  - run_tests.sh
  - scripts/run_tests.sh
  - scripts/testing/run_pending_tests.sh
- Manual test command: scripts/testing/run_and_log.sh A3 "docker compose exec fem-em-solver bash -lc 'cd /workspace && ./run_tests.sh --smoke'"
- Expected pass signal:
  - tests/unit/test_analytical_lightweight.py::test_straight_wire_analytical_direction_and_magnitude PASSED
  - tests/solver/test_tolerance_policy.py::test_solver_tolerance_policy_is_consistent PASSED
  - tests/validation/test_tolerance_policy.py::test_validation_tolerance_policy_is_ordered_and_positive PASSED
  - Pytest summary reports all selected smoke tests passed with no heavy mesh/solve commands executed
- Notes/blockers: none

- Chunk: A4 — Mesh-tag QA diagnostic hardening
- Status: 🧪 AWAITING-HUMAN-TEST
- Commit: 7c9b2c49cceb5f1035da23503e567ca242f6f821
- Files changed:
  - ROADMAP.md
  - src/fem_em_solver/io/__init__.py
  - src/fem_em_solver/io/mesh_qa.py
  - tests/io/test_mesh_qa_diagnostics.py
  - tests/mesh/helpers.py
- Manual test command: scripts/testing/run_and_log.sh A4 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/mesh/test_mesh_tag_integrity.py -v'"
- Expected pass signal:
  - tests/mesh/test_mesh_tag_integrity.py::test_coil_phantom_mesh_tag_integrity PASSED
  - On forced/missing-tag failures, log prints `[mesh-qa] required-tag expected vs actual:` with per-tag expected>=1 and actual counts
  - On forced/missing-tag failures, log prints `[mesh-qa] observed-tag summary:` including named required tags and unnamed tags as `tag_<id>`
- Notes/blockers: none


- Chunk: A5 — Testing status dashboard section
- Status: 🧪 AWAITING-HUMAN-TEST
- Commit: 82d2b6e4a0bf00bc521d3751d727e79dbdb6208e
- Files changed:
  - ROADMAP.md
  - docs/testing/pending-tests.md
- Manual test command: scripts/testing/run_and_log.sh A5 "docker compose exec fem-em-solver bash -lc 'cd /workspace && test -f docs/testing/pending-tests.md && echo OK'"
- Expected pass signal:
  - Log contains `OK`
  - Command exits with status `0`
  - `docs/testing/pending-tests.md` includes a `Testing Status Dashboard` table with columns: Chunk, Status, Commit, Last known log
- Notes/blockers: none

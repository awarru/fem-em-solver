# Pending Manual Tests

This file is updated by cron-driven coding runs in VPS-safe mode.

For each pending chunk, agent should append:
- Chunk ID / title
- Commit hash
- Files changed
- **Exact command using** `scripts/testing/run_and_log.sh`
- Expected pass signal
- Notes/blockers

Example:
```bash
scripts/testing/run_and_log.sh A1 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 examples/magnetostatics/01_straight_wire.py'"
```

Logs are automatically written to:
- `docs/testing/test-results.md` (summary index)
- `docs/testing/logs/*.log` (full output)

---

## B2 — Harden mesh QA checks (🧪 AWAITING-HUMAN-TEST)

- Commit: `3a600ec8036fc74db5382aa7dd5f822d2a2c5bfa`
- Files changed:
  - `ROADMAP.md`
  - `src/fem_em_solver/io/__init__.py`
  - `src/fem_em_solver/io/mesh_qa.py`
  - `tests/mesh/helpers.py`
  - `tests/mesh/test_mesh_tag_integrity.py`
- Manual test command:
  ```bash
  scripts/testing/run_and_log.sh B2 "docker compose exec fem-em-solver bash -lc 'cd /workspace && PYTHONPATH=/workspace/src mpiexec -n 2 python3 -m pytest tests/mesh/test_mesh_tag_integrity.py -v'"
  ```
- Expected pass signal:
  - Pytest reports `tests/mesh/test_mesh_tag_integrity.py::test_coil_phantom_mesh_tag_integrity PASSED`
  - Log contains a mesh QA summary line with non-zero counts for `coil_1`, `coil_2`, `phantom`, and `air`

---

## C1 — Solve B-field on coil+phantom model (🧪 AWAITING-HUMAN-TEST)

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

---

## C2 — Add sanity validation metrics (🧪 AWAITING-HUMAN-TEST)

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

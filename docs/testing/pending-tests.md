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

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

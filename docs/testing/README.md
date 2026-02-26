# Testing on a Fresh Clone (Other Computer)

Use this when you clone `fem-em-solver` on another machine and want to run pending manual tests reliably.

## 1) Clone and enter repo

```bash
git clone <repo-url>
cd fem-em-solver
```

## 2) One-command bootstrap + run (recommended)

```bash
./bootstrap_pending_tests.sh
```

This single command:
- starts/builds the Docker service (`fem-em-solver`)
- runs preflight checks
- executes all pending tests and logs results in-repo

## 3) Common options

List pending commands:
```bash
./bootstrap_pending_tests.sh --list
```

Run a single chunk:
```bash
./bootstrap_pending_tests.sh --chunk E3
```

Dry-run (no execution):
```bash
./bootstrap_pending_tests.sh --dry-run
```

(You can still use `./run_tests.sh` directly if environment is already prepared.)

---

## Output locations

- Full logs: `docs/testing/logs/*.log`
- Summary index: `docs/testing/test-results.md`
- Pending queue source: `docs/testing/pending-tests.md`

---

## Notes

- `run_tests.sh` routes to `scripts/testing/run_pending_tests.sh`.
- Per-test execution is wrapped by `scripts/testing/run_and_log.sh`.
- If your compose file is in the standard path (`docker/docker-compose.yml`), no extra setup is needed.
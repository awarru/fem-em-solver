"""Test configuration and fixtures.

Policy:
- tests/unit are lightweight and always run in CI
- all non-unit tests are marked `slow` automatically
  (they are intended for manual or larger machines)
"""

from __future__ import annotations

from pathlib import Path
import pytest


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line("markers", "slow: compute-intensive test")


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    for item in items:
        p = Path(str(item.fspath))
        # Any test outside tests/unit is considered heavy by default.
        if "tests" in p.parts and "unit" not in p.parts:
            item.add_marker(pytest.mark.slow)

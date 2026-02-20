"""Runtime path compatibility for container test commands.

The roadmap test commands set ``PYTHONPATH=/workspace/src`` which drops the
DolfinX path configured by the container environment. Add the known DolfinX
location back when available so imports remain stable.
"""

from __future__ import annotations

import os
import sys

_EXTRA_PATHS = [
    "/usr/local/dolfinx-real/lib/python3.10/dist-packages",
    "/usr/local/lib",
]

for _path in _EXTRA_PATHS:
    if os.path.isdir(_path) and _path not in sys.path:
        sys.path.append(_path)

"""Make the flat root modules and the scripts/ runners importable from tests.

The project ships flat top-level modules (``combicone``, ``reachability``, ...)
rather than a package, and ``scripts/`` is a plain directory, not a package. Both
directories therefore have to be on ``sys.path`` for the suite to import
``combicone`` and ``certificate_dossier`` alike. ``reproduce.sh`` exports
PYTHONPATH for the repo root, but ``pytest`` invoked directly does not, so this
is the single place that contract lives.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

for entry in (ROOT, ROOT / "scripts"):
    if str(entry) not in sys.path:
        sys.path.insert(0, str(entry))

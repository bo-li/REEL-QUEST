"""Pytest configuration for the repository root.

Puts the repo root on sys.path so tests can `import app`. app.py is the
Streamlit entry point and is deliberately not packaged (see pyproject.toml),
so it is not importable from an installed rel_quest alone.

This replaces the per-file sys.path.insert calls that used to head each test
module. Scripts under tests/manual/ keep their own, since they are run
standalone rather than through pytest, as does tests/_detail_probe_app.py,
which Streamlit's AppTest launches by file path.
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

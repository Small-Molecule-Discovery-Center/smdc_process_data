"""Shared pytest setup: make the repo's modules importable from tests/.

Run all tests from the repo root with:  pytest
"""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
for path in (REPO, REPO / "database_interactions"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

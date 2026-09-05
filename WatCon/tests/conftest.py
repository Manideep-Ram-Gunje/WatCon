"""Shared fixtures for the ConSurf parser test suite.

Fixture files live in ``WatCon/data/consurf/fixtures`` and are the *canonical*
copies of the real ConSurf outputs we hold.  There are four distinct datasets;
``1BRS_A_150`` additionally has a CRLF twin with byte-identical content, kept
deliberately so line-ending handling stays covered.

Full ConSurf result bundles (``.tar.gz``) stay under ``data/consurf/exploratory``
and are used only by the cross-check tests, which need the annotated PDBs.
"""

from __future__ import annotations

import tarfile
from pathlib import Path

import pytest

#: Repository root -- the directory containing the importable ``WatCon`` package.
REPO_ROOT = Path(__file__).resolve().parents[2]

DATA_ROOT = Path(__file__).resolve().parents[1] / "data" / "consurf"
FIXTURES = DATA_ROOT / "fixtures"
EXPLORATORY = DATA_ROOT / "exploratory"

# -- the four distinct real datasets ---------------------------------------
BRS_150 = FIXTURES / "1BRS_A_150.grades.txt"
BRS_150_CRLF = FIXTURES / "1BRS_A_150.crlf.grades.txt"
P00648_150 = FIXTURES / "P00648_150.grades.txt"
P00648_50 = FIXTURES / "P00648_50.grades.txt"
O7W_45 = FIXTURES / "7O7W_A_45.grades.txt"

#: Every distinct dataset, one entry per real ConSurf run we hold.
REAL_FIXTURES = [BRS_150, P00648_150, P00648_50, O7W_45]

#: Including the CRLF twin, for tests that should hold on every file on disk.
ALL_FIXTURES = REAL_FIXTURES + [BRS_150_CRLF]

#: Result bundles that contain an annotated PDB usable as a grade oracle.
#: (archive filename, basename of the members inside it)
BUNDLES = [
    ("1788241750_ConSurf.tar.gz", "7O7W_A"),
    ("1788241883_ConSurf.tar.gz", "P00648"),
    ("1788241897_ConSurf.tar.gz", "1BRS_A"),
]


def bundle_member(archive: str, member: str) -> bytes:
    """Read one member out of a ConSurf result bundle without extracting it."""
    with tarfile.open(EXPLORATORY / archive) as tar:
        handle = tar.extractfile(member)
        if handle is None:
            raise FileNotFoundError(f"{member} not found in {archive}")
        return handle.read()


@pytest.fixture(scope="session")
def fixtures_dir() -> Path:
    return FIXTURES


@pytest.fixture(scope="session")
def exploratory_dir() -> Path:
    return EXPLORATORY


# ---------------------------------------------------------------------------
# Optional-dependency handling
# ---------------------------------------------------------------------------
#
# Core WatCon needs MDAnalysis (and modeller); the consurf package deliberately
# needs neither.  Where MDAnalysis is absent, the core test modules fail at
# IMPORT time, and a collection error aborts the whole run -- including the
# consurf tests, which are perfectly runnable.
#
# Only modules that genuinely import MDAnalysis are skipped, and the skip is
# ANNOUNCED in the pytest header.  Silently dropping a module would hide real
# failures: test_inputs.py, for instance, needs no MDAnalysis at all and is
# currently failing for an unrelated reason, so it must keep running.
#
# No core WatCon module is modified.  When MDAnalysis is present, every module
# is collected exactly as before.

_REQUIRES_MDANALYSIS = [
    "test_general.py",
    "test_static.py",
]


def _has(module: str) -> bool:
    from importlib.util import find_spec

    try:
        return find_spec(module) is not None
    except (ImportError, ValueError):
        return False


_MDANALYSIS_AVAILABLE = _has("MDAnalysis")

collect_ignore = [] if _MDANALYSIS_AVAILABLE else list(_REQUIRES_MDANALYSIS)


def pytest_report_header(config):
    """Say out loud which modules are not being collected, and why."""
    if _MDANALYSIS_AVAILABLE:
        return "WatCon: MDAnalysis present -- collecting all test modules"
    return (
        "WatCon: MDAnalysis NOT installed -- skipping "
        + ", ".join(_REQUIRES_MDANALYSIS)
        + " (core network tests are therefore UNVERIFIED here)"
    )

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

# -- the protein tyrosine phosphatase family --------------------------------
#: One ConSurf run per member, chain A: {PDB id: (UniProt, gene)}.  Each has a
#: grades file and a CA-only extract of ConSurf's own annotated PDB (the grade
#: is written on every atom of a residue, so the CA line carries it; the
#: extract reproduces the full file's cross-check exactly, 1461/1461).
PTP_FAMILY = {
    "1AAX": ("P18031", "PTPN1"),
    "4GRZ": ("P29350", "PTPN6"),
    "1ZC0": ("P35236", "PTPN7"),
    "5HDE": ("Q05209", "PTPN12"),
    "3BRH": ("Q9Y2R2", "PTPN22"),
}


def ptp_grades(pdb_id: str) -> Path:
    return FIXTURES / f"{pdb_id}_A.grades.txt"


def ptp_annotated_ca(pdb_id: str) -> Path:
    return FIXTURES / f"{pdb_id}_A.consurf_ca.pdb"


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

# ---------------------------------------------------------------------------
# Data that a distribution deliberately does not carry
# ---------------------------------------------------------------------------
#
# The suite ships inside the wheel, so `pytest --pyargs WatCon.tests` is a thing
# a user can type.  Two sets of inputs are not in the wheel, on purpose:
#
#   * `data/consurf/exploratory/` -- 5 MB of raw ConSurf result bundles, kept in
#     the repository for provenance and excluded by MANIFEST.in.  The CI wheel
#     check asserts they are absent.
#   * `tests/inputs/` and `tests/water_dir/` -- raw structures used by the
#     older core tests.
#
# Run from an installed copy, the modules needing them produced 89 collection
# errors and 41 failures, which reads as a broken package rather than as
# "these tests need the repository".  They are now not collected there, and the
# header says so -- the same treatment MDAnalysis already gets above.

_EXPLORATORY_AVAILABLE = EXPLORATORY.is_dir() and any(EXPLORATORY.glob("*.tar.gz"))
_HERE = Path(__file__).resolve().parent
_INPUTS_AVAILABLE = (_HERE / "inputs").is_dir()
_WATER_DIR_AVAILABLE = (_HERE / "water_dir").is_dir()

#: Need the raw ConSurf bundles for their annotated-PDB oracle.
_REQUIRES_EXPLORATORY = [
    "test_consurf_crosscheck.py",
    "test_core_mapping.py",
    "test_end_to_end.py",
    "test_evolutionary.py",
    "test_evolutionary_clusters.py",
    "test_identity_check.py",
    "test_residue_index.py",
    "test_superpose.py",
]

#: Need the raw structures under tests/inputs/.
_REQUIRES_INPUTS = [
    "test_core_mapping.py",
    "test_evolutionary.py",
    "test_residue_index.py",
]

#: Need tests/water_dir/, a structure directory the older core tests analyse.
_REQUIRES_WATER_DIR = [
    "test_directed_networks.py",
    "test_general.py",
]


def _uncollectable():
    missing = []
    if not _MDANALYSIS_AVAILABLE:
        missing += _REQUIRES_MDANALYSIS
    if not _EXPLORATORY_AVAILABLE:
        missing += _REQUIRES_EXPLORATORY
    if not _INPUTS_AVAILABLE:
        missing += _REQUIRES_INPUTS
    if not _WATER_DIR_AVAILABLE:
        missing += _REQUIRES_WATER_DIR
    return sorted(set(missing))


collect_ignore = _uncollectable()

# `pytest_report_header` below is only called for a conftest at the rootdir, so
# running `pytest --pyargs WatCon.tests` from elsewhere -- which is exactly the
# installed case these skips exist for -- would drop modules silently. A warning
# surfaces in the summary wherever the suite is run from, and silence is the one
# outcome this file already says it does not want.
if collect_ignore:
    import warnings as _warnings

    _warnings.warn(
        "WatCon: not collecting %d test module(s) -- %s. %s"
        % (len(collect_ignore), ", ".join(collect_ignore),
           "Run the suite from a repository checkout to cover them."),
        UserWarning,
        stacklevel=2,
    )


def pytest_report_header(config):
    """Say out loud which modules are not being collected, and why."""
    lines = []
    if _MDANALYSIS_AVAILABLE:
        lines.append("WatCon: MDAnalysis present -- core network tests collected")
    else:
        lines.append(
            "WatCon: MDAnalysis NOT installed -- skipping "
            + ", ".join(_REQUIRES_MDANALYSIS)
            + " (core network tests are therefore UNVERIFIED here)"
        )

    if _EXPLORATORY_AVAILABLE and _INPUTS_AVAILABLE and _WATER_DIR_AVAILABLE:
        lines.append("WatCon: repository test data present -- full suite collected")
    else:
        absent = []
        if not _EXPLORATORY_AVAILABLE:
            absent.append("data/consurf/exploratory")
        if not _INPUTS_AVAILABLE:
            absent.append("tests/inputs")
        if not _WATER_DIR_AVAILABLE:
            absent.append("tests/water_dir")
        lines.append(
            "WatCon: running from an installed copy -- %s not distributed, so "
            "%d module(s) are not collected. Run the suite from a repository "
            "checkout to cover them."
            % (", ".join(absent),
               len(set(_REQUIRES_EXPLORATORY + _REQUIRES_INPUTS + _REQUIRES_WATER_DIR)))
        )
    return lines

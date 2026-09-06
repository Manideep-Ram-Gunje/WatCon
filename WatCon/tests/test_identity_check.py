"""Does the ConSurf file actually describe THIS structure?

The join key is ``(chain, resid, icode)``.  Nothing in it encodes which amino
acid the score belongs to, so a ConSurf run whose numbering differs from the
structure's will attach a score to every residue, report 100% coverage, and be
wrong everywhere.  That failure is invisible downstream -- the waters get
plausible numbers and the report looks fine.

These tests pin the only place it is detectable: comparing the residue ConSurf
recorded against the residue the structure actually contains.

The threshold is a rate rather than a boolean because two different things
produce mismatches.  A point mutant genuinely differs from the sequence ConSurf
aligned, and that is correct data.  A numbering offset differs almost
everywhere.  The measured separation on real barnase is wide -- 100% correct,
99.1% for one mutation, 2.8% for a one-residue shift -- and these tests assert
that both ends stay on the right side of it.
"""

from __future__ import annotations

import dataclasses

import pytest

from WatCon.evolutionary import (
    DEFAULT_IDENTITY_THRESHOLD,
    ConservationError,
    CoverageReport,
    enforce_identity,
    load_conservation,
)
from .conftest import FIXTURES

# ===========================================================================
# Fixtures: real barnase, real ConSurf run
# ===========================================================================

@pytest.fixture(scope="module")
def barnase_map():
    return load_conservation(
        "1BRS_A_150.grades.txt", FIXTURES, strict=False
    )


@pytest.fixture(scope="module")
def barnase_residues():
    """Chain A of 1BRS, read from the structure ConSurf itself annotated."""
    import os
    import tempfile

    from WatCon.residue_index import residues_from_pdb_file

    from .conftest import bundle_member

    path = os.path.join(tempfile.mkdtemp(), "1BRS.pdb")
    with open(path, "wb") as handle:
        handle.write(
            bundle_member(
                "1788241897_ConSurf.tar.gz",
                "1BRS_A_ATOMS_section_With_ConSurf.pdb",
            )
        )
    return residues_from_pdb_file(path, chain="A")


# ===========================================================================
# The three regimes, on real data
# ===========================================================================

def test_correct_structure_agrees_completely(barnase_map, barnase_residues):
    report = barnase_map.coverage(barnase_residues)
    assert report.identity_checked == 108
    assert report.identity_mismatched == 0
    assert report.identity_rate == pytest.approx(1.0)


def test_a_numbering_offset_is_caught(barnase_map, barnase_residues):
    """The failure this whole check exists for."""
    shifted = [
        dataclasses.replace(r, resid=r.resid + 1) for r in barnase_residues
    ]
    report = barnase_map.coverage(shifted)

    # Coverage alone would not reveal this -- almost everything still matches.
    assert report.matched > 100
    # Identity does.
    assert report.identity_rate < 0.10
    assert report.identity_rate < DEFAULT_IDENTITY_THRESHOLD


def test_a_point_mutant_still_passes(barnase_map, barnase_residues):
    """One real mutation must not be mistaken for a broken mapping."""
    index = next(
        i for i, r in enumerate(barnase_residues) if r.resname != "ALA"
    )
    mutated = list(barnase_residues)
    mutated[index] = dataclasses.replace(
        mutated[index], resname="ALA", one_letter="A"
    )

    report = barnase_map.coverage(mutated)
    assert report.identity_mismatched == 1
    assert report.identity_rate > DEFAULT_IDENTITY_THRESHOLD
    # And it says which residue, so the user can check it against the paper.
    (key, observed, expected) = report.mismatches[0]
    assert observed == "ALA"
    assert expected == barnase_residues[index].resname


def test_the_threshold_separates_the_two_regimes(barnase_map, barnase_residues):
    """A single mutation and a numbering shift must land either side of it."""
    index = next(
        i for i, r in enumerate(barnase_residues) if r.resname != "ALA"
    )
    mutated = list(barnase_residues)
    mutated[index] = dataclasses.replace(
        mutated[index], resname="ALA", one_letter="A"
    )
    shifted = [
        dataclasses.replace(r, resid=r.resid + 1) for r in barnase_residues
    ]

    mutant_rate = barnase_map.coverage(mutated).identity_rate
    shifted_rate = barnase_map.coverage(shifted).identity_rate

    assert shifted_rate < DEFAULT_IDENTITY_THRESHOLD <= mutant_rate


# ===========================================================================
# Degrading, rather than lying, when no comparison is possible
# ===========================================================================

def test_objects_without_residue_names_are_not_checked(barnase_map, barnase_residues):
    """Absent names must read as 'not checked', never as 'all mismatched'."""

    class Bare:
        def __init__(self, residue):
            self.chain = residue.chain
            self.resid = residue.resid
            self.icode = residue.icode

    report = barnase_map.coverage([Bare(r) for r in barnase_residues])
    assert report.matched == 108          # coverage still works
    assert report.identity_checked == 0   # identity simply was not assessed
    assert report.identity_rate is None
    assert "not checked" in report.describe_identity()


def test_one_letter_is_used_when_there_is_no_three_letter_name(barnase_map):
    """A ConSurf run with no model attached still supports a check."""
    record = next(
        r for r in barnase_map.result.records if r.pdb_residue is not None
    )
    key = record.pdb_residue.key

    class OneLetterOnly:
        def __init__(self, key, letter):
            self.chain, self.resid, self.icode = key
            self.one_letter = letter

    good = barnase_map.coverage([OneLetterOnly(key, record.sequence_residue)])
    assert good.identity_matched == 1

    wrong = "W" if record.sequence_residue != "W" else "G"
    bad = barnase_map.coverage([OneLetterOnly(key, wrong)])
    assert bad.identity_mismatched == 1


def test_each_residue_is_counted_once_not_once_per_atom(barnase_map):
    """Handed atoms, the rate must measure residues -- not atom counts."""
    record = next(
        r for r in barnase_map.result.records if r.pdb_residue is not None
    )
    chain, resid, icode = record.pdb_residue.key

    class Atom:
        def __init__(self):
            self.chain, self.resid, self.icode = chain, resid, icode
            self.resname = record.pdb_residue.name3

    report = barnase_map.coverage([Atom() for _ in range(9)])
    assert report.matched == 9            # nine atoms were covered
    assert report.identity_checked == 1   # of one residue


# ===========================================================================
# Enforcement
# ===========================================================================

def make_report(matched, mismatched):
    report = CoverageReport()
    report.identity_matched = matched
    report.identity_mismatched = mismatched
    report.mismatches = [
        (("A", i, None), "ALA", "VAL") for i in range(mismatched)
    ]
    return report


def test_strict_refuses_a_mismatched_file():
    with pytest.raises(ConservationError, match="does not describe"):
        enforce_identity(make_report(3, 97), label="1ABC.pdb", strict=True)


def test_the_error_names_the_structure_and_shows_examples():
    with pytest.raises(ConservationError) as excinfo:
        enforce_identity(make_report(3, 97), label="1ABC.pdb", strict=True)
    message = str(excinfo.value)
    assert "1ABC.pdb" in message
    assert "structure has ALA, ConSurf has VAL" in message
    assert "3 of 100" in message


def test_tolerant_mode_warns_instead_of_raising():
    with pytest.warns(RuntimeWarning, match="does not describe"):
        rate = enforce_identity(make_report(3, 97), strict=False)
    assert rate == pytest.approx(0.03)


def test_a_passing_rate_is_silent(recwarn):
    rate = enforce_identity(make_report(99, 1), strict=True)
    assert rate == pytest.approx(0.99)
    assert len(recwarn) == 0


def test_unchecked_identity_never_blocks(recwarn):
    """No residue names is not evidence of a mismatch."""
    assert enforce_identity(CoverageReport(), strict=True) is None
    assert len(recwarn) == 0


def test_threshold_is_configurable():
    report = make_report(90, 10)          # 0.90
    assert enforce_identity(report, strict=True, threshold=0.85) == pytest.approx(0.90)
    with pytest.raises(ConservationError):
        enforce_identity(report, strict=True, threshold=0.95)


def test_the_real_barnase_file_passes_enforcement(barnase_map, barnase_residues):
    report = barnase_map.coverage(barnase_residues)
    assert enforce_identity(report, label="1BRS", strict=True) == pytest.approx(1.0)

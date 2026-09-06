"""Pooling conservation from SEPARATE ConSurf runs, onto shared MSA columns.

The single-protein path needs none of this: one ConSurf run describes every
structure of that sequence.  A protein *family* is different -- each member has
its own run, and residue 40 of one protein is not residue 40 of another.  The
correspondence has to go through the alignment.

The trap these tests guard is not the join, which is simple.  It is the
interpretation.  ConSurf scores are z-normalised **within each run**, so pooled
values are only relatively comparable, and grades are per-run percentile bins
that must never be averaged.  The real cross-run degradation is measured here on
two genuine ConSurf runs of the same barnase sequence at different MSA depths --
the easiest possible case, and still not perfect agreement.
"""

from __future__ import annotations

import pytest

from WatCon.evolutionary import (
    CROSS_RUN_AGREEMENT,
    ColumnConservation,
    ConservationError,
    ConservationMap,
    conservation_by_msa_column,
    family_summary,
    load_conservation,
)
from WatCon.residue_index import ResidueIndex, StructureResidue

from .conftest import FIXTURES


# ===========================================================================
# Helpers
# ===========================================================================

def index_for(conservation_map, chain="A"):
    """A ResidueIndex over the residues a ConSurf file covers.

    The MSA column is the residue's ordinal, i.e. a self-alignment.  That is the
    correct alignment when pooling runs of the *same* sequence, which is what
    the real-data tests below do.
    """
    residues = []
    for record in conservation_map.result.records:
        if record.pdb_residue is None:
            continue
        residues.append(
            StructureResidue(
                chain=record.pdb_residue.chain,
                resid=record.pdb_residue.seq,
                icode=record.pdb_residue.icode,
                resname=record.pdb_residue.name3,
                one_letter=record.sequence_residue,
            )
        )
    return ResidueIndex(residues, list(range(len(residues))))


def fake_map(entries):
    """A ConservationMap from {(chain, resid, icode): (score, grade)}."""
    from WatCon.evolutionary import ResidueConservation

    conservation_map = ConservationMap.__new__(ConservationMap)
    conservation_map._by_residue = {
        key: ResidueConservation(
            score=score,
            grade=grade,
            low_confidence=False,
            msa_present=50,
            msa_total=50,
            confidence_lower=None,
            confidence_upper=None,
            buried_exposed=None,
            functional_structural=None,
            consurf_position=i + 1,
            amino_acid="A",
            residue_name="ALA",
            source="synthetic",
        )
        for i, (key, (score, grade)) in enumerate(entries.items())
    }
    conservation_map.result = None
    conservation_map.source = "synthetic"
    return conservation_map


def synthetic_member(label, scores):
    """(label, map, index) for residues 1..N of chain A with the given scores."""
    entries = {("A", i + 1, None): s for i, s in enumerate(scores)}
    residues = [
        StructureResidue("A", i + 1, None, "ALA", "A") for i in range(len(scores))
    ]
    return label, fake_map(entries), ResidueIndex(residues, list(range(len(scores))))


# ===========================================================================
# The join
# ===========================================================================

def test_members_meet_at_the_alignment_column():
    a = synthetic_member("a", [(-1.5, 9), (0.2, 5), (1.1, 2)])
    b = synthetic_member("b", [(-1.2, 9), (0.4, 4), (0.9, 3)])

    columns = conservation_by_msa_column([a, b])
    assert sorted(columns) == [0, 1, 2]
    assert columns[0].n_members == 2
    assert columns[0].scores == {"a": -1.5, "b": -1.2}
    assert columns[0].grades == {"a": 9, "b": 9}


def test_a_column_only_one_member_covers_is_kept_and_labelled():
    """Partial coverage is information, not something to drop silently."""
    a = synthetic_member("a", [(-1.5, 9), (0.2, 5)])
    b = synthetic_member("b", [(-1.2, 9)])

    columns = conservation_by_msa_column([a, b])
    assert columns[0].n_members == 2
    assert columns[1].n_members == 1
    assert set(columns[1].scores) == {"a"}


def test_pooling_without_an_alignment_is_refused():
    """Residue numbers are not comparable between different proteins."""
    label, conservation_map, index = synthetic_member("a", [(-1.0, 9)])
    no_msa = ResidueIndex(index.residues, None)

    with pytest.raises(ConservationError, match="no MSA mapping"):
        conservation_by_msa_column([(label, conservation_map, no_msa)])


def test_grades_are_never_averaged():
    """A grade is an ordinal bin; its mean is not a grade."""
    a = synthetic_member("a", [(-2.0, 9)])
    b = synthetic_member("b", [(0.5, 3)])

    column = conservation_by_msa_column([a, b])[0]
    assert column.max_grade == 9
    assert not hasattr(column, "mean_grade")


def test_min_score_follows_the_sign_convention():
    """More negative = more conserved, so the minimum is the most conserved."""
    a = synthetic_member("a", [(-2.0, 9)])
    b = synthetic_member("b", [(0.5, 3)])

    column = conservation_by_msa_column([a, b])[0]
    assert column.min_score == -2.0
    assert column.mean_score == pytest.approx(-0.75)


def test_spread_exposes_disagreement_between_runs():
    agree = conservation_by_msa_column(
        [synthetic_member("a", [(-1.5, 9)]), synthetic_member("b", [(-1.4, 9)])]
    )[0]
    disagree = conservation_by_msa_column(
        [synthetic_member("a", [(-2.0, 9)]), synthetic_member("b", [(1.0, 2)])]
    )[0]

    assert agree.spread == pytest.approx(0.1)
    assert disagree.spread == pytest.approx(3.0)
    assert disagree.spread > CROSS_RUN_AGREEMENT["grade_changed_fraction"]


def test_spread_is_undefined_for_a_single_member():
    assert conservation_by_msa_column([synthetic_member("a", [(-1.0, 9)])])[0].spread is None


def test_unanimous_conserved_uses_each_runs_own_verdict():
    """The one family claim that does not compare values across runs."""
    both = conservation_by_msa_column(
        [synthetic_member("a", [(-2.0, 9)]), synthetic_member("b", [(-1.1, 8)])]
    )[0]
    one = conservation_by_msa_column(
        [synthetic_member("a", [(-2.0, 9)]), synthetic_member("b", [(0.1, 5)])]
    )[0]

    assert both.unanimous_conserved
    assert not one.unanimous_conserved


# ===========================================================================
# Summary
# ===========================================================================

def test_family_summary_counts_what_it_says():
    columns = conservation_by_msa_column(
        [
            synthetic_member("a", [(-2.0, 9), (-1.5, 8), (0.9, 2)]),
            synthetic_member("b", [(-1.8, 9), (0.3, 5), (1.0, 2)]),
        ]
    )
    summary = family_summary(columns)
    assert summary["n_columns"] == 3
    assert summary["n_all_members"] == 3
    assert summary["n_unanimous_conserved"] == 1     # only the first column
    assert summary["median_spread"] is not None


def test_family_summary_of_nothing_is_not_an_error():
    summary = family_summary({})
    assert summary["n_columns"] == 0
    assert summary["median_spread"] is None


# ===========================================================================
# Real data: two genuine ConSurf runs of the SAME sequence
# ===========================================================================

@pytest.fixture(scope="module")
def barnase_two_runs():
    """P00648 at MSA depth 150 and 50 -- same protein, two independent runs."""
    deep = load_conservation("P00648_150.grades.txt", FIXTURES, strict=False)
    shallow = load_conservation("P00648_50.grades.txt", FIXTURES, strict=False)
    return [
        ("depth150", deep, index_for(deep)),
        ("depth50", shallow, index_for(shallow)),
    ]


def test_two_runs_of_one_sequence_pool_completely(barnase_two_runs):
    columns = conservation_by_msa_column(barnase_two_runs)
    assert len(columns) > 100
    assert all(c.n_members == 2 for c in columns.values())


def test_the_same_sequence_still_does_not_agree_perfectly(barnase_two_runs):
    """The honest floor on any family-level claim.

    Two runs of the *identical* sequence disagree; different proteins will do
    worse.  If this ever came out as perfect agreement, the two files would not
    be independent runs and the measured noise floor would be wrong.
    """
    columns = conservation_by_msa_column(barnase_two_runs)

    changed = sum(
        1 for c in columns.values() if len(set(c.grades.values())) > 1
    )
    fraction = changed / len(columns)

    assert 0.0 < fraction < 1.0
    # Matches the 43% recorded in CROSS_RUN_AGREEMENT, within a loose band.
    assert 0.25 < fraction < 0.60


def test_spread_across_real_runs_is_measurable(barnase_two_runs):
    columns = conservation_by_msa_column(barnase_two_runs)
    summary = family_summary(columns)

    assert summary["n_columns"] == len(columns)
    assert summary["median_spread"] > 0.0
    assert summary["n_unanimous_conserved"] > 0


def test_unanimity_is_stricter_than_either_run_alone(barnase_two_runs):
    """Requiring both runs to agree must discard some single-run calls."""
    columns = conservation_by_msa_column(barnase_two_runs)

    unanimous = sum(1 for c in columns.values() if c.unanimous_conserved)
    either = sum(
        1 for c in columns.values() if any(g >= 8 for g in c.grades.values())
    )
    assert 0 < unanimous < either

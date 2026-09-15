"""The five real PTP ConSurf runs: sound, poolable, and biologically sane.

These are the inputs to the first real multi-protein family analysis. Before
anything is pooled across them, three things have to hold, and each is checked
against something independent of our own parser's say-so:

* **Each run parses to what ConSurf itself wrote.** Grades are compared with the
  B-factors ConSurf wrote into its own annotated PDB -- an artefact produced
  separately from the grades table. 1461 of 1461 residues agree.
* **The runs are methodologically poolable.** Same dialect, same method, same
  homologue depth. Pooling a Bayesian run with a maximum-likelihood one, or a
  150-sequence alignment with a 50-sequence one, would mix scales that
  ``CROSS_RUN_AGREEMENT`` shows are already noisy for a single protein.
* **The catalytic machinery is graded as conserved.** This is a positive control
  from biochemistry, not from our code: PTPs share a P-loop H-C-x5-R, a WPD-loop
  general acid and a Q-loop glutamine. A pipeline that failed to grade those as
  conserved would be wrong regardless of what its tests said.

Three of the ConSurf *queries* are substrate-trapping mutants -- 1AAX C215S,
4GRZ C453S, 3BRH C227S with D195A, all confirmed by the PDB's SEQADV records.
ConSurf grades the alignment **column**, so those grades still describe the
family; the tests assert both facts so neither is forgotten downstream.
"""

from __future__ import annotations

import pytest

from WatCon.consurf import check_against_annotated_pdb, parse_consurf

from .conftest import PTP_FAMILY, ptp_annotated_ca, ptp_grades

IDS = sorted(PTP_FAMILY)


@pytest.fixture(scope="module")
def runs():
    return {pdb_id: parse_consurf(ptp_grades(pdb_id), strict=True) for pdb_id in IDS}


def _record(run, resid):
    return run.by_pdb_residue()[("A", resid, None)]


# ===========================================================================
# Each run is what ConSurf wrote
# ===========================================================================

#: (records, mapped to the structure, low-confidence) per run.
SHAPE = {
    "1AAX": (321, 297, 7),
    "4GRZ": (288, 282, 4),
    "1ZC0": (309, 286, 19),
    "5HDE": (307, 300, 9),
    "3BRH": (310, 296, 12),
}


@pytest.mark.parametrize("pdb_id", IDS)
def test_parses_strictly_with_nothing_malformed(runs, pdb_id):
    run = runs[pdb_id]
    records, mapped, low = SHAPE[pdb_id]
    assert len(run.records) == records
    assert len(run.mapped_records()) == mapped
    assert sum(r.low_confidence for r in run.records) == low
    assert run.malformed_lines == []


@pytest.mark.parametrize("pdb_id", IDS)
def test_every_grade_agrees_with_consurfs_own_annotated_pdb(runs, pdb_id):
    report = check_against_annotated_pdb(runs[pdb_id], ptp_annotated_ca(pdb_id))
    assert report.mismatches == []
    assert report.missing_from_pdb == []
    assert report.blank_in_pdb == []
    assert report.checked == report.agreements == SHAPE[pdb_id][1]


def test_1461_residues_agree_in_total(runs):
    total = sum(
        check_against_annotated_pdb(runs[p], ptp_annotated_ca(p)).agreements
        for p in IDS
    )
    assert total == 1461


def test_only_1zc0_reports_a_numbering_discontinuity(runs):
    """1ZC0 skips 178 -> 183 in its own numbering; ConSurf says so."""
    codes = {p: {w.code for w in runs[p].warnings} for p in IDS}
    assert "pdb_number_discontinuity" in codes["1ZC0"]
    assert all(not codes[p] for p in IDS if p != "1ZC0")


@pytest.mark.parametrize("pdb_id", IDS)
def test_fixture_line_endings_are_lf(pdb_id):
    """All five bundles were LF as downloaded; a CRLF fixture is a local artefact."""
    assert b"\r\n" not in ptp_grades(pdb_id).read_bytes()


# ===========================================================================
# Poolable
# ===========================================================================

def test_all_runs_share_dialect_method_and_depth(runs):
    provenance = {
        (r.provenance.dialect.value, r.provenance.method.value, r.provenance.msa_total)
        for r in runs.values()
    }
    assert provenance == {("webserver", "bayesian", 150)}


# ===========================================================================
# The catalytic machinery -- a positive control from biochemistry
# ===========================================================================

#: First residue (His) of the P-loop H-C/S-x5-R, and its grades in order.
#: The sixth position is the one variable site in the motif.
P_LOOP = {
    "1AAX": (214, [9, 9, 9, 9, 9, 8, 9, 9]),
    "4GRZ": (452, [9, 9, 9, 9, 9, 8, 9, 9]),
    "1ZC0": (269, [9, 9, 9, 9, 9, 8, 9, 9]),
    "5HDE": (230, [9, 9, 9, 9, 9, 7, 9, 9]),
    "3BRH": (226, [9, 9, 9, 9, 9, 7, 9, 9]),
}


@pytest.mark.parametrize("pdb_id", IDS)
def test_p_loop_is_conserved(runs, pdb_id):
    start, expected = P_LOOP[pdb_id]
    grades = [_record(runs[pdb_id], start + i).grade for i in range(8)]
    assert grades == expected
    assert _record(runs[pdb_id], start).sequence_residue == "H"
    assert _record(runs[pdb_id], start + 7).sequence_residue == "R"


#: The WPD-loop general-acid column: residue number, and the amino acid in the
#: ConSurf query (3BRH is D195A).
WPD_ACID = {
    "1AAX": (181, "D"),
    "4GRZ": (419, "D"),
    "1ZC0": (236, "D"),
    "5HDE": (199, "D"),
    "3BRH": (195, "A"),
}


@pytest.mark.parametrize("pdb_id", IDS)
def test_wpd_general_acid_column_is_grade_9_and_aspartate_in_homologues(runs, pdb_id):
    resid, query = WPD_ACID[pdb_id]
    record = _record(runs[pdb_id], resid)
    assert record.grade == 9
    assert record.sequence_residue == query
    assert record.residue_variety_raw.startswith("D")


#: Q-loop glutamine, which positions the catalytic water.
Q_LOOP_GLN = {"1AAX": 262, "4GRZ": 500, "1ZC0": 314, "5HDE": 278, "3BRH": 274}


@pytest.mark.parametrize("pdb_id", IDS)
def test_q_loop_glutamine_is_grade_9(runs, pdb_id):
    record = _record(runs[pdb_id], Q_LOOP_GLN[pdb_id])
    assert record.sequence_residue == "Q"
    assert record.grade == 9
    assert record.residue_variety_raw.startswith("Q")


# ===========================================================================
# Substrate traps -- recorded, so nothing downstream forgets them
# ===========================================================================

#: Catalytic nucleophile: residue number, query amino acid, residue name
#: ConSurf read from its PDB.
NUCLEOPHILE = {
    "1AAX": (215, "S", "SER"),
    "4GRZ": (453, "S", "SER"),
    "1ZC0": (270, "C", "CYS"),
    "5HDE": (231, "C", "CSP"),     # phosphocysteine intermediate, not a mutation
    "3BRH": (227, "S", "SER"),
}


@pytest.mark.parametrize("pdb_id", IDS)
def test_nucleophile_identity_is_what_the_pdb_entry_deposited(runs, pdb_id):
    resid, query, name3 = NUCLEOPHILE[pdb_id]
    record = _record(runs[pdb_id], resid)
    assert record.sequence_residue == query
    assert record.pdb_residue.name3 == name3


@pytest.mark.parametrize("pdb_id", IDS)
def test_a_mutant_query_still_grades_the_family_column(runs, pdb_id):
    """ConSurf grades the column: a C->S query still reads Cys-conserved."""
    resid = NUCLEOPHILE[pdb_id][0]
    record = _record(runs[pdb_id], resid)
    assert record.grade == 9
    assert record.residue_variety_raw.startswith("C")


def test_5hde_phosphocysteine_is_graded_by_consurf(runs):
    """ConSurf scores CSP231. If WatCon later loses it, that is WatCon's bug."""
    record = _record(runs["5HDE"], 231)
    assert record.pdb_residue.name3 == "CSP"
    assert record.grade == 9

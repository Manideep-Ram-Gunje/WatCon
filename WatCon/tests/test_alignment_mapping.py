"""Placing real PTP structures on the authors' family alignment.

The fixtures are real: CA atoms of chain A from the RCSB entries (native
numbering) and the alignment published with WatCon (Zenodo
10.5281/zenodo.15213225, CC-BY-4.0), all twenty-four rows of it. The tests here
exercise the five proteins that have a structure pair, because the defects worth
catching only show up when two structures of one protein disagree.

That alignment has two defects, kept deliberately because the mapper must catch
them rather than reproduce them:

* the **5HDE** row omits CSP231, the catalytic phosphocysteine;
* the **3O4U** row slides E124, D125 and E241 across disordered loops, putting
  E241 in the column every other PTP fills with the WPD general-acid aspartate.

Anchor columns are never hard-coded: they are read off PTP1B's own residues
(D181, Q262, C215) through the mapping, so a change in column convention cannot
make these tests pass by accident.
"""

from __future__ import annotations

import os

import pytest

pytest.importorskip("Bio", reason="sequence mapping needs Biopython")

from WatCon.alignment import (
    AlignmentError,
    consensus_columns,
    map_structure_to_row,
    read_alignment,
    residue_index_for,
)
from WatCon.residue_index import ResidueIndex, residues_from_pdb_file

PACKAGE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FAMILY = os.path.join(PACKAGE, "data", "examples", "ptp_family")
ALIGNMENT = os.path.join(FAMILY, "ptp_family_alignment.pir")

MEMBERS = {
    "PTPN1": ("2F71", "8U1E"),
    "PTPN6": ("4GRZ", "4HJP"),
    "PTPN7": ("1ZC0", "3O4U"),
    "PTPN12": ("5HDE", "5J8R"),
    "PTPN22": ("3BRH", "3OLR"),
}

pytestmark = pytest.mark.skipif(not os.path.isfile(ALIGNMENT),
                                reason="PTP family fixtures absent")


@pytest.fixture(scope="module")
def rows():
    return read_alignment(ALIGNMENT)


def _row_for(rows, pdb_id):
    return next((name, seq) for name, seq in rows.items() if name.startswith(pdb_id))


@pytest.fixture(scope="module")
def residues():
    return {pdb_id: residues_from_pdb_file(os.path.join(FAMILY, "%s_A_ca.pdb" % pdb_id), chain="A")
            for pair in MEMBERS.values() for pdb_id in pair}


@pytest.fixture(scope="module")
def mappings(rows, residues):
    out = {}
    for pair in MEMBERS.values():
        for pdb_id in pair:
            name, row = _row_for(rows, pdb_id)
            out[pdb_id] = map_structure_to_row(residues[pdb_id], row, label=pdb_id, row_name=name)
    return out


@pytest.fixture(scope="module")
def consensus(mappings):
    return {protein: consensus_columns(protein, [mappings[a], mappings[b]])
            for protein, (a, b) in MEMBERS.items()}


# ===========================================================================
# Reading
# ===========================================================================

def test_reads_every_row_at_one_length(rows):
    """One row per structure in the family: fifteen proteins, twenty-four rows."""
    assert len(rows) == 24
    assert {len(s) for s in rows.values()} == {347}


def test_a_ragged_file_is_not_an_alignment(tmp_path):
    bad = tmp_path / "bad.fa"
    bad.write_text(">a\nACDE\n>b\nAC\n")
    with pytest.raises(AlignmentError, match="not an alignment"):
        read_alignment(str(bad))


# ===========================================================================
# Each structure onto its own row
# ===========================================================================

@pytest.mark.parametrize("pdb_id", [p for pair in MEMBERS.values() for p in pair])
def test_each_structure_matches_its_own_row_exactly(mappings, pdb_id):
    assert mappings[pdb_id].identity == 1.0
    assert mappings[pdb_id].mismatches == []


def test_only_5hde_csp231_has_no_letter_in_its_row(mappings):
    unmapped = {p: m.unmapped for p, m in mappings.items() if m.unmapped}
    assert unmapped == {"5HDE": [(231, None)]}


def test_the_wrong_row_is_refused(rows, residues):
    """1ZC0's residues against PTP1B's row: a different protein, not a mapping."""
    name, row = _row_for(rows, "2F71")
    with pytest.raises(AlignmentError, match="does not match"):
        map_structure_to_row(residues["1ZC0"], row, label="1ZC0", row_name=name)


# ===========================================================================
# Consensus within a protein -- where the defects are caught
# ===========================================================================

def test_rows_agree_everywhere_for_four_proteins(consensus):
    for protein in ("PTPN1", "PTPN6", "PTPN12", "PTPN22"):
        assert consensus[protein].conflicts == {}, protein


def test_the_3o4u_slides_are_caught_as_conflicts(consensus):
    assert set(consensus["PTPN7"].conflicts) == {(124, None), (125, None), (241, None)}
    for position in ((124, None), (125, None), (241, None)):
        assert consensus["PTPN7"].column(*position) is None


def test_5hde_csp231_takes_its_column_from_5j8r(consensus):
    ptpn12 = consensus["PTPN12"]
    assert ptpn12.column(231) is not None
    assert ptpn12.single_source[(231, None)] == "5J8R"


def test_no_column_is_shared_by_two_residues_of_one_protein(consensus):
    for protein, c in consensus.items():
        values = list(c.columns.values())
        assert len(values) == len(set(values)), protein


# ===========================================================================
# Anchors: catalytic residues of different proteins meet in one column
# ===========================================================================

def _letter_at(residues, consensus_for, column):
    by_column = {consensus_for.column(r.resid, r.icode): r for r in residues}
    residue = by_column.get(column)
    return None if residue is None else residue.one_letter


def test_wpd_general_acid_column(residues, consensus):
    column = consensus["PTPN1"].column(181)
    assert _letter_at(residues["2F71"], consensus["PTPN1"], column) == "D"
    expected = {"8U1E": "D", "4GRZ": "D", "4HJP": "D", "1ZC0": "D",
                "3O4U": None,      # disordered loop; its slid E241 was excluded
                "5HDE": "D", "5J8R": "D",
                "3BRH": "A",       # engineered D195A
                "3OLR": "D"}
    for protein, pair in MEMBERS.items():
        for pdb_id in pair:
            if pdb_id in expected:
                assert _letter_at(residues[pdb_id], consensus[protein], column) == expected[pdb_id], pdb_id


def test_q_loop_glutamine_column(residues, consensus):
    column = consensus["PTPN1"].column(262)
    for protein, pair in MEMBERS.items():
        for pdb_id in pair:
            assert _letter_at(residues[pdb_id], consensus[protein], column) == "Q", pdb_id


def test_nucleophile_column_holds_cys_or_its_engineered_serine(residues, consensus):
    column = consensus["PTPN1"].column(215)
    expected = {"2F71": "C", "8U1E": "C", "4GRZ": "S", "4HJP": "C", "1ZC0": "C",
                "3O4U": "C", "5HDE": "C",  # CSP231, via 5J8R's row
                "5J8R": "C", "3BRH": "S", "3OLR": "S"}
    for protein, pair in MEMBERS.items():
        for pdb_id in pair:
            assert _letter_at(residues[pdb_id], consensus[protein], column) == expected[pdb_id], pdb_id


# ===========================================================================
# Into the existing family join
# ===========================================================================

def test_residue_index_carries_none_for_unplaced_residues(residues, consensus):
    index = residue_index_for(residues["3O4U"], consensus["PTPN7"])
    assert index.has_msa
    assert index.msa_column("A", 241) is None
    assert index.msa_column("A", 234) is not None


def test_a_plain_residue_index_still_accepts_none():
    from WatCon.residue_index import StructureResidue

    residues = [StructureResidue("A", 1, None, "GLY", "G"), StructureResidue("A", 2, None, "ALA", "A")]
    index = ResidueIndex.build(residues, [5, None])
    assert index.msa_column("A", 1) == 5
    assert index.msa_column("A", 2) is None

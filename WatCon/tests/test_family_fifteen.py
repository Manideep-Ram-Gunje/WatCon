"""The family at fifteen: fifteen ConSurf runs, twenty-four structures.

The five-protein family in :mod:`test_family` establishes that the machinery is
correct. This module asks the different question fifteen proteins make possible:
does the family statement *survive* being asked of three times as many
independently normalised runs?

It is a real risk. ``unanimous_conserved`` means every run graded the column 8 or
9, so each protein added is another chance for a column to fall out. If the
five-protein number were an artefact of a small, closely related sample, it would
collapse here.

The members are the fifteen human PTPs the WatCon authors assembled (Zenodo
10.5281/zenodo.15213225, CC-BY-4.0), each with its own ConSurf run over 150
homologues, all chain A, all Bayesian. Two of them are worth naming:

* **PTPRN2** (2QEP) is a pseudophosphatase. Its P-loop reads ``CSDGAGR`` -- the
  nucleophile is there but the family-invariant positions are not -- which is why
  IA-2beta has no catalytic activity. It is kept, and the tests check that
  keeping it changes nothing about the conservation statement.
* **PTPN9** (6KZQ) carries a C->A rather than the usual C->S trap, so the
  nucleophile column is not simply "C or S".
"""

from __future__ import annotations

import os

import pytest

pytest.importorskip("Bio", reason="family mapping needs Biopython")

from WatCon.family import FamilyProtein, FamilyStructure, build_family_conservation

from .conftest import ptp_grades

PACKAGE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FAMILY_DIR = os.path.join(PACKAGE, "data", "examples", "ptp_family")
ALIGNMENT = os.path.join(FAMILY_DIR, "ptp_family_alignment.pir")

pytestmark = pytest.mark.skipif(not os.path.isfile(ALIGNMENT), reason="PTP family fixtures absent")

#: protein -> (ConSurf run, reference structure, structures)
LAYOUT = {
    "PTPN1": ("1AAX", "2F71", ("2F71", "8U1E")),
    "PTPN3": ("4S0G", "4S0G", ("4S0G", "2B49")),
    "PTPN4": ("2I75", "2I75", ("2I75",)),
    "PTPN5": ("8SLS", "8SLS", ("8SLS",)),
    "PTPN6": ("4GRZ", "4GRZ", ("4GRZ", "4HJP")),
    "PTPN7": ("1ZC0", "1ZC0", ("1ZC0", "3O4U")),
    "PTPN9": ("6KZQ", "6KZQ", ("6KZQ", "4GE6")),
    "PTPN11": ("3ZM1", "3ZM1", ("3ZM1",)),
    "PTPN12": ("5HDE", "5HDE", ("5HDE", "5J8R")),
    "PTPN13": ("1WCH", "1WCH", ("1WCH",)),
    "PTPN14": ("6IWD", "6IWD", ("6IWD",)),
    "PTPN18": ("4GFU", "4GFU", ("4GFU", "2OC3")),
    "PTPN21": ("8GVV", "8GVV", ("8GVV", "8GWH")),
    "PTPN22": ("3BRH", "3BRH", ("3BRH", "3OLR")),
    "PTPRN2": ("2QEP", "2QEP", ("2QEP",)),
}

#: The nucleophile each entry actually models, from its own SEQADV record. The
#: point of listing them is that the column below is conserved *despite* this.
NUCLEOPHILE = {
    "2F71": "C", "8U1E": "C", "4S0G": "S", "2B49": "S", "2I75": "C", "8SLS": "C",
    "4GRZ": "S", "4HJP": "C", "1ZC0": "C", "3O4U": "C", "6KZQ": "A", "4GE6": "A",
    "3ZM1": "C", "5HDE": "C", "5J8R": "C", "1WCH": "C", "6IWD": "C", "4GFU": "S",
    "2OC3": "S", "8GVV": "S", "8GWH": "S", "3BRH": "S", "3OLR": "S", "2QEP": "C",
}


def _proteins(layout=LAYOUT):
    return [
        FamilyProtein(
            name=name,
            consurf_path=str(ptp_grades(run)),
            structures=[FamilyStructure(pdb_id, os.path.join(FAMILY_DIR, "%s_A_ca.pdb" % pdb_id))
                        for pdb_id in structures],
            reference=reference,
        )
        for name, (run, reference, structures) in layout.items()
    ]


@pytest.fixture(scope="module")
def family():
    return build_family_conservation(_proteins(), ALIGNMENT)


@pytest.fixture(scope="module")
def without_pseudophosphatase():
    layout = {k: v for k, v in LAYOUT.items() if k != "PTPRN2"}
    return build_family_conservation(_proteins(layout), ALIGNMENT)


# ===========================================================================
# Everything is placed
# ===========================================================================

def test_all_fifteen_proteins_and_twenty_four_structures(family):
    assert len(family.proteins) == 15
    assert len(family.structures) == 24


def test_every_structure_matches_its_own_alignment_row_exactly(family):
    """Row identity 1.000 for all twenty-four, with one residue unplaced.

    That residue is 5HDE's CSP231, which the published alignment row omits; it
    is recovered from 5J8R by the consensus, so it is not lost.
    """
    unplaced = {r.pdb_id: r.unplaced for r in family.structures if r.unplaced}
    for report in family.structures:
        assert report.row_identity == 1.0, report.pdb_id
    assert unplaced == {"5HDE": [(231, None)]}


def test_every_structure_is_the_protein_its_consurf_run_describes(family):
    for report in family.structures:
        assert report.identity is not None and report.identity >= 0.99, report.pdb_id


def test_the_only_disagreements_are_declared_mutations(family):
    """Each named difference is an engineered substitution or a modified residue.

    Listed explicitly rather than counted, because a silent change here is
    exactly the failure the identity check exists to catch.
    """
    found = {r.pdb_id: sorted("%s%d %s>%s" % (k[0], k[1], a, b)
                              for k, a, b in r.identity_mismatches)
             for r in family.structures if r.identity_mismatches}
    assert found == {
        "2F71": ["A215 CYS>SER"],            # run is 1AAX, the C215S trap
        "8U1E": ["A215 CYS>SER"],
        "2B49": ["A811 ASP>ALA", "A842 CYS>SER"],
        "4HJP": ["A453 CYS>SER"],            # run is 4GRZ, the C453S trap
        "3O4U": ["A72 ASP>SER"],
        "4GE6": ["A470 ASP>ALA", "A515 CYS>ALA"],
        "5J8R": ["A231 CYS>CSP", "A61 ARG>LYS"],   # run is 5HDE, phospho-Cys
        "2OC3": ["A229 CYS>SER"],            # run is 4GFU, the C229S trap
        "3OLR": ["A195 ASP>ALA"],            # run is 3BRH, which also has D195A
    }


def test_residues_two_structures_disagree_about_are_excluded(family):
    """Four residues in the whole family, all in disordered loops.

    PTPN7's three are the published alignment's known slide (E124, D125, E241).
    PTPN18's single one appears only now that 2OC3 joins 4GFU. Excluding them is
    the right outcome: a residue whose column two structures of one protein
    disagree about cannot be pooled.
    """
    assert family.conflicts() == {
        "PTPN7": [(124, None), (125, None), (241, None)],
        "PTPN18": [(44, None)],
    }


# ===========================================================================
# The statement survives the larger family
# ===========================================================================

def test_the_pooled_counts(family):
    summary = family.summary()
    assert summary["n_proteins"] == 15
    assert len(family.columns) == 344
    assert summary["n_columns_all_proteins"] == 229
    assert summary["n_unanimous_all_proteins"] == 58


def test_unanimity_is_not_an_artefact_of_a_small_family(family):
    """The five-protein family gave 64 of 248; fifteen gives 58 of 229.

    Tripling the number of independently normalised runs moves the *fraction*
    of universally covered columns that every run calls conserved by less than a
    percentage point. That is the whole reason for pooling verdicts rather than
    separately normalised scores.
    """
    summary = family.summary()
    fraction = summary["n_unanimous_all_proteins"] / summary["n_columns_all_proteins"]
    assert 0.24 < fraction < 0.27
    assert abs(fraction - 64 / 248) < 0.01


def test_dropping_the_pseudophosphatase_changes_nothing(family, without_pseudophosphatase):
    """PTPRN2 is catalytically dead, and it costs the family no conserved column."""
    assert without_pseudophosphatase.summary()["n_unanimous_all_proteins"] == 58
    assert len(without_pseudophosphatase.columns) == len(family.columns) == 344


# ===========================================================================
# Biochemistry the code was never told
# ===========================================================================

@pytest.mark.parametrize("protein,resid", [
    ("PTPN1", 215),    # the nucleophile, by PTP1B numbering
    ("PTPN1", 221),    # the P-loop arginine
    ("PTPN1", 262),    # the Q-loop glutamine
])
def test_catalytic_columns_are_unanimously_conserved_across_fifteen(family, protein, resid):
    column = family.at(protein, resid)
    assert column is not None
    assert column.n_members == 15
    assert column.unanimous_conserved, sorted(column.grades.items())


def test_the_nucleophile_column_holds_every_engineered_substitution(family):
    """Cys in the wild types, Ser in the traps, Ala in PTPN9's -- one column.

    The column is located through PTP1B's C215 and read out per structure, so
    this fails if the placement drifts, not merely if a grade changes.
    """
    column = family.column_for("PTPN1", 215)
    seen = {}
    for protein, residues in family.reference_residues.items():
        for residue in residues:
            if family.column_for(protein, residue.resid, residue.icode) == column:
                seen[protein] = residue.one_letter
    assert len(seen) == 15, sorted(seen)
    for protein, letter in seen.items():
        reference = LAYOUT[protein][1]
        assert letter == NUCLEOPHILE[reference], protein


def test_six_of_the_seven_p_loop_columns_are_graded_nine_by_all_fifteen(family):
    """The strongest single statement the family makes.

    Not "conserved on average" and not "conserved in our reference" -- fifteen
    separately normalised runs, each over its own 150 homologues, every one of
    them putting the strictest available grade on the same six columns.
    """
    graded_nine = [resid for resid in range(215, 222)
                   if min(family.at("PTPN1", resid).grades.values()) == 9]
    assert graded_nine == [215, 216, 217, 218, 220, 221]


def test_the_variable_p_loop_position_is_the_one_that_is_not_unanimous(family):
    """Position 219 is I, V or C across the family, and it is the only gap."""
    not_unanimous = [resid for resid in range(215, 222)
                     if not family.at("PTPN1", resid).unanimous_conserved]
    assert not_unanimous == [219]


def test_conservation_measures_constraint_not_identity(family):
    """PTPRN2 reads CSDGAGR where the family reads CSAGIGR, and is still graded 9.

    The pseudophosphatase is the case that separates the two ideas. At PTP1B's
    A217 the family has alanine and PTPRN2 has aspartate -- a different residue,
    in a protein that cannot catalyse -- yet every run including PTPRN2's own
    grades that column 9, because it is invariant among PTPRN2 orthologues too.
    A tool that reported identity would call this position divergent. ConSurf
    calls it constrained, and it is right to.
    """
    letters = {}
    for protein, residues in family.reference_residues.items():
        for residue in residues:
            column = family.column_for(protein, residue.resid, residue.icode)
            if column is not None:
                letters.setdefault(column, {})[protein] = residue.one_letter

    column = family.column_for("PTPN1", 217)
    assert letters[column]["PTPRN2"] == "D"
    assert {v for k, v in letters[column].items() if k != "PTPRN2"} == {"A"}
    assert family.columns[column].grades["PTPRN2"] == 9
    assert family.columns[column].unanimous_conserved

    # And where PTPRN2 does diverge in grade, it is the variable position.
    variable = family.at("PTPN1", 219)
    assert variable.grades["PTPRN2"] == min(variable.grades.values()) == 7

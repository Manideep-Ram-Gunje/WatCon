"""The PTP family, end to end: five ConSurf runs, ten structures, one alignment.

Real data throughout: the five ConSurf grades files, CA-only chain-A extracts of
the ten RCSB structures (native numbering), and the ten corresponding rows of the
alignment published with WatCon (Zenodo 10.5281/zenodo.15213225, CC-BY-4.0).

What must hold, and why each is a meaningful check rather than a restatement:

* **Every structure is the protein its ConSurf run describes.** The only residue
  disagreements are the ones the PDB entries declare as engineered (SEQADV) or
  modified (MODRES) -- nothing else.
* **The alignment's two defects are handled, not reproduced.** 3O4U's slid
  residues are excluded; 5HDE's missing CSP231 is placed via 5J8R.
* **Catalytic machinery is conserved in every protein independently.** The
  P-loop His, the Cys nucleophile, the P-loop Arg, the WPD general acid and the
  Q-loop Gln are graded 9 by all five runs. This is biochemistry, not our code.
* **The runs agree although they are independent.** Their homologue sets are
  nearly disjoint (Jaccard <= 0.01), and their scores still correlate above 0.9.
"""

from __future__ import annotations

import os

import pytest

pytest.importorskip("Bio", reason="family mapping needs Biopython")

from WatCon.evolutionary import ConservationError
from WatCon.family import FamilyProtein, FamilyStructure, build_family_conservation

from .conftest import ptp_grades

PACKAGE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FAMILY_DIR = os.path.join(PACKAGE, "data", "examples", "ptp_family")
ALIGNMENT = os.path.join(FAMILY_DIR, "ptp_family_alignment.pir")

pytestmark = pytest.mark.skipif(not os.path.isfile(ALIGNMENT), reason="PTP family fixtures absent")

#: protein -> (ConSurf run, reference structure, structures)
LAYOUT = {
    "PTPN1": ("1AAX", "2F71", ("2F71", "8U1E")),
    "PTPN6": ("4GRZ", "4GRZ", ("4GRZ", "4HJP")),
    "PTPN7": ("1ZC0", "1ZC0", ("1ZC0", "3O4U")),
    "PTPN12": ("5HDE", "5HDE", ("5HDE", "5J8R")),
    "PTPN22": ("3BRH", "3BRH", ("3BRH", "3OLR")),
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


# ===========================================================================
# Every structure is the protein its ConSurf run describes
# ===========================================================================

def test_all_ten_structures_pass_the_identity_check(family):
    assert len(family.structures) == 10
    assert all(report.identity >= 0.99 for report in family.structures)


#: Residue disagreements with the protein's ConSurf query, per structure:
#: {pdb id: {resid: (structure residue, ConSurf residue)}}. Every one is declared
#: by a PDB entry -- SEQADV engineered mutations, or MODRES for 5HDE's CSP231.
EXPECTED_MISMATCHES = {
    "2F71": {215: ("CYS", "SER")},          # 1AAX query is C215S
    "8U1E": {215: ("CYS", "SER")},
    "4GRZ": {},
    "4HJP": {453: ("CYS", "SER")},          # 4GRZ query is C453S
    "1ZC0": {},
    "3O4U": {72: ("ASP", "SER")},           # 3O4U S72D (UniProt 93)
    "5HDE": {},
    "5J8R": {61: ("ARG", "LYS"), 231: ("CYS", "CSP")},  # 5J8R K61R; 5HDE's CSP231
    "3BRH": {},
    "3OLR": {195: ("ASP", "ALA")},          # 3BRH query is D195A
}


def test_the_only_disagreements_are_declared_by_the_pdb_entries(family):
    observed = {
        report.pdb_id: {key[1]: (structure, consurf) for key, structure, consurf in report.identity_mismatches}
        for report in family.structures
    }
    assert observed == EXPECTED_MISMATCHES


def test_every_structure_matches_its_own_alignment_row(family):
    assert all(report.row_identity == 1.0 for report in family.structures)


def test_the_wrong_consurf_run_is_refused():
    """PTPN7's structures against PTPN6's ConSurf run: different protein, stop."""
    wrong = FamilyProtein(
        name="PTPN7",
        consurf_path=str(ptp_grades("4GRZ")),
        structures=[FamilyStructure("1ZC0", os.path.join(FAMILY_DIR, "1ZC0_A_ca.pdb"))],
    )
    with pytest.raises(ConservationError, match="does not describe"):
        build_family_conservation([wrong], ALIGNMENT)


# ===========================================================================
# The alignment's defects are handled, not reproduced
# ===========================================================================

def test_only_the_3o4u_slides_are_excluded(family):
    assert family.conflicts() == {"PTPN7": [(124, None), (125, None), (241, None)]}


def test_5hde_csp231_is_pooled_through_5j8rs_row(family):
    column = family.at("PTPN12", 231)
    assert column is not None
    assert column.grades["PTPN12"] == 9


# ===========================================================================
# Catalytic machinery -- conserved in every protein independently
# ===========================================================================

#: Placed through PTP1B's own residues, so no column number is hard-coded.
CATALYTIC = {
    "P-loop His": 214,
    "nucleophile Cys": 215,
    "P-loop Arg": 221,
    "WPD general acid": 181,
    "Q-loop Gln": 262,
}


@pytest.mark.parametrize("label", sorted(CATALYTIC))
def test_catalytic_column_is_grade_9_in_all_five_runs(family, label):
    column = family.at("PTPN1", CATALYTIC[label])
    assert column is not None, label
    assert column.n_members == 5
    assert set(column.grades.values()) == {9}
    assert column.unanimous_conserved


def test_a_substrate_loop_residue_is_not_unanimous(family):
    """Tyr46 is conserved, but PTPN12's run grades it 7: unanimity is a real bar."""
    column = family.at("PTPN1", 46)
    assert column.grades["PTPN12"] == 7
    assert not column.unanimous_conserved


# ===========================================================================
# Independent runs, and they agree
# ===========================================================================

def test_independent_runs_agree_across_proteins(family):
    correlations = family.pairwise_spearman()
    assert len(correlations) == 10
    assert all(n >= 250 for n, _ in correlations.values())
    assert all(rho > 0.90 for _, rho in correlations.values())


def test_summary_reports_columns_every_protein_covers(family):
    """Pinned. Identical whether built from these CA extracts or from prepared
    full structures, which is itself a check that the extracts are faithful."""
    summary = family.summary()
    assert summary["n_proteins"] == 5
    assert summary["n_columns"] == 337
    assert summary["n_columns_all_proteins"] == 248
    assert summary["n_unanimous_all_proteins"] == 64
    assert summary["median_spread_all_proteins"] == pytest.approx(0.365, abs=0.001)


# ===========================================================================
# A member whose ConSurf run matches nothing
#
# Found by walking the manual-test guide: giving one protein another's grades
# file -- a plausible slip in a members file -- produced a complete,
# plausible-looking family result with that protein contributing no
# conservation at all, and said nothing. Coverage of zero passes
# `enforce_identity` vacuously, because with no matched residues there is no
# identity to check.
#
# The static and dynamic paths warn about this. Here it is an error: a member
# is listed precisely so that its conservation is pooled.
# ===========================================================================

def test_a_member_whose_run_matches_nothing_is_refused():
    """PTPN1's structures with PTPN13's run: no residue numbers in common."""
    layout = dict(LAYOUT)
    layout["PTPN1"] = ("1WCH",) + LAYOUT["PTPN1"][1:]

    with pytest.raises(ConservationError, match="matched none of its residues"):
        build_family_conservation(_proteins(layout), ALIGNMENT)


def test_that_refusal_says_the_chains_are_not_the_problem():
    """Both are chain A, so suggesting a chain map would misdirect."""
    layout = dict(LAYOUT)
    layout["PTPN1"] = ("1WCH",) + LAYOUT["PTPN1"][1:]

    with pytest.raises(ConservationError) as raised:
        build_family_conservation(_proteins(layout), ALIGNMENT)

    message = str(raised.value)
    assert "Those agree, so the cause is not chain labelling" in message
    assert "members file pairs each protein with its own ConSurf run" in message
    assert "consurf_chain_map" not in message


def test_tolerant_mode_warns_instead_of_stopping(capsys):
    """--tolerant is documented as warn-instead-of-stop; this is no exception."""
    layout = dict(LAYOUT)
    layout["PTPN1"] = ("1WCH",) + LAYOUT["PTPN1"][1:]

    family = build_family_conservation(_proteins(layout), ALIGNMENT, strict=False)
    assert "matched none of its residues" in capsys.readouterr().out
    assert len(family.proteins) == 5

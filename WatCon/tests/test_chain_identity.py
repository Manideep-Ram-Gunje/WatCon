"""Two ways the ConSurf join could attach nothing and not say so.

The join key is ``(chain, resid, icode)``. This package already refuses a run
whose residue *numbering* disagrees with its ConSurf file, because attaching
scores to the wrong residues silently is the failure mode that matters. These
tests cover the other half of the same key -- the chain -- and the residue name
the check compares, both of which could fail quietly instead of loudly.

Found by running the dynamic path on the Zenodo MD system:

* **The chain was read from the wrong attribute.** ``chainID`` is an *atom*-level
  attribute in MDAnalysis, so reading it off a residue always gave ``None`` and
  the code fell through to ``segid``. For a file straight from the PDB that is
  right by luck, because MDAnalysis fills segid from the chain column. For a PDB
  *written* by MDAnalysis the segid is the invented ``SYST`` while the real chain
  sits in ``chainID``, so every lookup missed and every residue came back
  unscored. 127 of 127 active-site atoms lost their grade, with no warning.

* **Protonation states read as sequence differences.** A force field writes HID,
  HIE or HIP where a crystal structure writes HIS. Those are the same amino
  acid, but the identity check compared three-letter names, so an MD system
  looked like it disagreed with its own sequence.
"""

from __future__ import annotations

import os

import pytest

pytest.importorskip("MDAnalysis", reason="chain identity is read through MDAnalysis")

import MDAnalysis as mda

from WatCon.consurf import parse_consurf
from WatCon.evolutionary import ConservationMap, enforce_identity
from WatCon.residue_index import (
    PROTONATION_VARIANTS,
    _universe_chain_id,
    residues_from_pdb_file,
    standard_residue_name,
)

PACKAGE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ENSEMBLE = os.path.join(PACKAGE, "data", "examples", "ptp1b_ensemble")
MD_SITE = os.path.join(ENSEMBLE, "md_active_site_h.pdb")
FAMILY = os.path.join(PACKAGE, "data", "examples", "ptp_family")
GRADES = os.path.join(PACKAGE, "data", "consurf", "fixtures", "1AAX_A.grades.txt")

#: The MD system is PTP1B renumbered from 1, measured: +1 gives identity 1.000
#: over 295 residues and every other offset is below 0.09. Declared here rather
#: than applied anywhere silently.
MD_OFFSET = 1

pytestmark = pytest.mark.skipif(not os.path.isfile(MD_SITE),
                                reason="MD active-site fixture absent")


# ===========================================================================
# The chain comes from the chain, whatever wrote the file
# ===========================================================================

def _column_22(path):
    """The chain exactly as the PDB text says it, which is what the join uses."""
    for line in open(path, encoding="utf-8", errors="replace"):
        if line.startswith(("ATOM", "HETATM")):
            return line[21].strip()
    return ""


@pytest.mark.parametrize("fixture", [
    os.path.join(FAMILY, "2F71_A_ca.pdb"),      # straight from the PDB: chain A
    MD_SITE,                                    # written by MDAnalysis: chain X
])
def test_both_readers_agree_about_the_chain(fixture):
    """The invariant the join depends on.

    WatCon reads structures two ways -- a text parser for residue identity and
    MDAnalysis for the networks -- and the ConSurf join is keyed on the chain
    both of them produce. If they disagree, every lookup misses and nothing
    says so. That is exactly what happened: MDAnalysis-written files resolved
    to the invented segid 'SYST' while the text parser read the real chain.
    """
    universe = mda.Universe(fixture)
    residue = universe.select_atoms("protein").residues[0]
    from_text = residues_from_pdb_file(fixture)[0].chain

    assert _universe_chain_id(residue) == from_text == _column_22(fixture)


def test_the_invented_segid_is_never_taken_for_a_chain():
    """The MD fixture is the case: segid SYST, real chain X."""
    universe = mda.Universe(MD_SITE)
    assert set(universe.atoms.segids) == {"SYST"}
    resolved = _universe_chain_id(universe.select_atoms("protein").residues[0])
    assert resolved == "X"
    assert resolved not in ("SYST", "SYSTEM")


def test_a_blank_chain_column_reports_no_chain(tmp_path):
    """Nothing is better than a guess: a blank column reads as no chain."""
    stripped = str(tmp_path / "stripped.pdb")
    with open(stripped, "w", encoding="utf-8", newline="\n") as handle:
        for line in open(MD_SITE, encoding="utf-8"):
            if line.startswith(("ATOM", "HETATM")) and len(line) > 22:
                line = line[:21] + " " + line[22:]
            handle.write(line)
    universe = mda.Universe(stripped)
    assert _universe_chain_id(universe.select_atoms("protein").residues[0]) == ""


# ===========================================================================
# A protonation state is not a sequence difference
# ===========================================================================

def test_the_variants_map_to_their_parent():
    assert standard_residue_name("HID") == "HIS"
    assert standard_residue_name("HIE") == "HIS"
    assert standard_residue_name("HIP") == "HIS"
    assert standard_residue_name("CYM") == "CYS"
    assert standard_residue_name("GLH") == "GLU"
    assert standard_residue_name("hid") == "HIS"


def test_a_chemical_modification_is_left_alone():
    """CSP really is not CYS: 5HDE's phosphocysteine must stay reportable."""
    assert standard_residue_name("CSP") == "CSP"
    assert "CSP" not in PROTONATION_VARIANTS
    assert standard_residue_name("SER") == "SER"


def test_the_md_system_agrees_with_its_own_sequence(tmp_path):
    """The whole point, end to end.

    Renumbered by the declared +1 offset, the MD active site should disagree
    with the 1AAX run at exactly one residue -- 1AAX is the C215S trap and the
    MD system models a real cysteine. Before the fixes this reported two
    differences, the extra one being HID against HIS.
    """
    renumbered = str(tmp_path / "md_site_renumbered.pdb")
    with open(renumbered, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("REMARK   Renumbered by %+d and chain set to A, both declared,\n"
                     "REMARK   to match the 1AAX ConSurf run.\n" % MD_OFFSET)
        for line in open(MD_SITE, encoding="utf-8"):
            if line.startswith(("ATOM", "HETATM")):
                line = (line[:21] + "A"
                        + "%4d" % (int(line[22:26]) + MD_OFFSET) + line[26:])
            handle.write(line)

    conservation = ConservationMap.build(parse_consurf(GRADES, strict=True))
    residues = residues_from_pdb_file(renumbered, chain="A")
    coverage = conservation.coverage(residues)

    named = ["%s%d %s>%s" % (k[0], k[1], a, b) for k, a, b in coverage.mismatches]
    assert named == ["A215 CYS>SER"]
    assert enforce_identity(coverage, strict=False) > 0.97

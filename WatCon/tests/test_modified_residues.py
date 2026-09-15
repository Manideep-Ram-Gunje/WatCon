"""In-chain modified residues are protein, and their waters count.

5HDE's catalytic nucleophile is **CSP231**, the phosphocysteine intermediate of
the PTP reaction, deposited as a HETATM. ConSurf grades it 9. WatCon lost it
three separate ways:

* the residue walk read ATOM records only, so 5HDE had 299 of its 300 residues;
* MDAnalysis's ``protein`` keyword does not include CSP (0 of its atoms), so the
  network builders never made it a node;
* no front end passed a ``custom_selection`` that could have added it back.

So the phosphate that hydrogen-bonds three waters (563, 598, 601, at 2.76-3.11 A)
contributed nothing -- no conservation, no contacts -- and nothing said so.

The fixture is real: 5HDE trimmed to 12 A around CSP231, native numbering, so
the real 5HDE ConSurf run attaches to it directly.
"""

from __future__ import annotations

import os

import pytest

from WatCon.residue_index import (
    AMINO_ACID_3TO1,
    MODIFIED_RESIDUES,
    protein_selection,
    residues_from_pdb_file,
)

from .conftest import ptp_grades

PACKAGE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SITE = os.path.join(PACKAGE, "data", "examples", "ptp_family", "5HDE_A_csp_site.pdb")

needs_site = pytest.mark.skipif(not os.path.isfile(SITE), reason="5HDE fixture absent")


# ===========================================================================
# The residue walk
# ===========================================================================

@needs_site
def test_the_walk_reads_the_phosphocysteine():
    residues = {r.key: r for r in residues_from_pdb_file(SITE)}
    csp = residues[("A", 231, None)]
    assert csp.resname == "CSP"
    assert csp.one_letter == "C"


@needs_site
def test_alternate_conformers_still_collapse_to_one_residue():
    """CSP231 is deposited with altlocs A and B; it is one residue."""
    keys = [r.key for r in residues_from_pdb_file(SITE)]
    assert keys.count(("A", 231, None)) == 1


@needs_site
def test_waters_are_still_not_residues():
    assert all(r.resname != "HOH" for r in residues_from_pdb_file(SITE))


def test_a_free_standard_name_ligand_stays_a_ligand(tmp_path):
    """Only curated in-chain modifications are read from HETATM by default."""
    pdb = tmp_path / "ligand.pdb"
    pdb.write_text(
        "ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  CA  TYR A 900      10.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    3  CA  MSE A   2       3.800   0.000   0.000  1.00  0.00           C\n"
        "END\n"
    )
    names = [r.resname for r in residues_from_pdb_file(str(pdb))]
    assert names == ["GLY", "MSE"]
    with_all = [r.resname for r in residues_from_pdb_file(str(pdb), include_hetatm=True)]
    assert with_all == ["GLY", "TYR", "MSE"]


def test_modified_residues_map_to_their_parent():
    assert MODIFIED_RESIDUES["CSP"] == "C"
    assert MODIFIED_RESIDUES["PTR"] == "Y"
    assert MODIFIED_RESIDUES["TPO"] == "T"
    assert all(AMINO_ACID_3TO1[name] == parent for name, parent in MODIFIED_RESIDUES.items())


def test_the_chromophore_is_still_not_a_modified_residue():
    """7O7W's fused PIA must stay excluded -- it is not a single residue."""
    assert "PIA" not in MODIFIED_RESIDUES
    assert "PIA" not in AMINO_ACID_3TO1


# ===========================================================================
# The selection the network builders use
# ===========================================================================

def test_protein_selection_names_every_modified_residue():
    clause = protein_selection()
    assert clause.startswith("(protein or resname ")
    for name in MODIFIED_RESIDUES:
        assert name in clause


@needs_site
def test_mdanalysis_protein_alone_misses_csp_and_the_helper_does_not():
    """Documents the bug, and the fix, against MDAnalysis itself."""
    mda = pytest.importorskip("MDAnalysis")
    universe = mda.Universe(SITE)
    assert len(universe.select_atoms("protein and resname CSP")) == 0
    polar = universe.select_atoms(
        "%s and resname CSP and (name N* or name O* or name P* or name S*)"
        % protein_selection())
    assert {"P", "O1P", "O2P", "O3P", "SG"} <= set(polar.names)


@needs_site
def test_no_water_or_ion_is_pulled_in_as_protein():
    mda = pytest.importorskip("MDAnalysis")
    universe = mda.Universe(SITE)
    assert set(universe.select_atoms(protein_selection()).resnames).isdisjoint({"HOH", "PO4"})


# ===========================================================================
# End to end: conservation and water contacts reach CSP231
# ===========================================================================

@pytest.fixture(scope="module")
def network():
    pytest.importorskip("MDAnalysis")
    if not os.path.isfile(SITE):
        pytest.skip("5HDE fixture absent")
    from WatCon.consurf import parse_consurf
    from WatCon.evolutionary import ConservationMap
    from WatCon.generate_static_networks import extract_objects

    conservation = ConservationMap.build(parse_consurf(ptp_grades("5HDE"), strict=True))
    return extract_objects(
        SITE, "water-protein", None,
        active_region_reference=None, active_region_COM=False,
        active_region_radius=8.0, water_name=None, msa_indexing=None,
        max_connection_distance=3.3, conservation_map=conservation,
    ), conservation


def test_csp231_atoms_are_network_nodes_carrying_grade_9(network):
    net, _ = network
    csp = [a for a in net.protein_atoms if a.resname == "CSP"]
    assert csp, "CSP231 contributed no atoms to the network"
    assert all(int(a.resid) == 231 for a in csp)
    assert all(a.evolutionary is not None and a.evolutionary.grade == 9 for a in csp)


def test_csp231_is_in_the_water_contacts(network):
    from WatCon.evolutionary import water_residue_contacts

    net, _ = network
    contacted = {key for residues in water_residue_contacts(net).values() for key in residues}
    assert ("A", 231, None) in contacted


def test_identity_check_passes_with_csp_compared_as_csp(network):
    from WatCon.evolutionary import enforce_identity

    net, conservation = network
    coverage = conservation.coverage(net.protein_atoms)
    assert enforce_identity(coverage, label="5HDE site") == 1.0
    assert coverage.identity_mismatched == 0

"""Alternate conformers: drop a position only when another is more populated.

MDAnalysis keeps every alternate position as its own atom, so a residue modelled
in two conformers gave WatCon two copies of each hydrogen-bonding atom and two
edges for one water contact -- about 70% of water-protein edges in the PanDDA
ensemble models of PTP1B.

The rule: keep every conformer tied for the highest mean occupancy, drop only
those strictly lower. Ties are the norm -- 92% of alternate-conformer residues in
the PTP1B set, 95% in barnase -- so a "first-listed wins" tie-break would decide
almost everything arbitrarily. On barnase it removed a grade-9 site's only lining
residue, Asp54, whose water touched conformer B at 0.50.

The real-fixture test does not hard-code the answer: it re-derives it from the
PDB text independently and compares.
"""

from __future__ import annotations

import collections
import os

import pytest

mda = pytest.importorskip("MDAnalysis")

from WatCon.conformers import (
    conformer_choices,
    has_alternate_conformers,
    highest_occupancy_atoms,
)

PACKAGE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SITE = os.path.join(PACKAGE, "data", "examples", "ptp_family", "5HDE_A_csp_site.pdb")


def _atom(record, serial, name, altloc, resname, chain, resid, x, y, z, occ, element):
    field = name if len(name) == 4 else " %-3s" % name
    return ("%-6s%5d %4s%1s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s"
            % (record, serial, field, altloc, resname, chain, resid, x, y, z, occ, 10.0, element))


def _write(tmp_path, lines, name="model.pdb"):
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\nEND\n")
    return str(path)


def _network(path):
    from WatCon.generate_static_networks import extract_objects

    return extract_objects(path, "water-protein", None, active_region_reference=None,
                           active_region_COM=False, active_region_radius=8.0,
                           water_name=None, msa_indexing=None, max_connection_distance=3.3)


# ===========================================================================
# The choice
# ===========================================================================

def test_the_more_populated_conformer_wins_even_when_listed_second(tmp_path):
    path = _write(tmp_path, [
        _atom("ATOM", 1, "CA", "", "SER", "A", 10, 0, 0, 0, 1.00, "C"),
        _atom("ATOM", 2, "OG", "A", "SER", "A", 10, 1, 0, 0, 0.30, "O"),
        _atom("ATOM", 3, "OG", "B", "SER", "A", 10, 0, 1, 0, 0.70, "O"),
    ])
    atoms = mda.Universe(path).atoms
    assert conformer_choices(atoms) == {("A", 10, "", "SER"): ("B",)}
    kept = highest_occupancy_atoms(atoms)
    assert sorted(kept.names) == ["CA", "OG"]
    assert list(kept.select_atoms("name OG").altLocs) == ["B"]


def test_tied_conformers_are_all_kept(tmp_path):
    """0.50/0.50 is no evidence for either; neither is discarded."""
    path = _write(tmp_path, [
        _atom("ATOM", 1, "OG", "A", "SER", "A", 10, 1, 0, 0, 0.50, "O"),
        _atom("ATOM", 2, "OG", "B", "SER", "A", 10, 0, 1, 0, 0.50, "O"),
    ])
    atoms = mda.Universe(path).atoms
    assert conformer_choices(atoms) == {("A", 10, "", "SER"): ("A", "B")}
    assert sorted(highest_occupancy_atoms(atoms).altLocs) == ["A", "B"]


def test_only_the_strictly_less_populated_position_is_dropped(tmp_path):
    """A and B tied at 0.40, C at 0.20: C goes, A and B stay."""
    path = _write(tmp_path, [
        _atom("ATOM", 1, "OG", "A", "SER", "A", 10, 1, 0, 0, 0.40, "O"),
        _atom("ATOM", 2, "OG", "B", "SER", "A", 10, 0, 1, 0, 0.40, "O"),
        _atom("ATOM", 3, "OG", "C", "SER", "A", 10, 0, 0, 1, 0.20, "O"),
    ])
    assert sorted(highest_occupancy_atoms(mda.Universe(path).atoms).altLocs) == ["A", "B"]


def test_residues_are_decided_independently(tmp_path):
    path = _write(tmp_path, [
        _atom("ATOM", 1, "OG", "A", "SER", "A", 10, 0, 0, 0, 0.80, "O"),
        _atom("ATOM", 2, "OG", "B", "SER", "A", 10, 0, 1, 0, 0.20, "O"),
        _atom("ATOM", 3, "OG1", "A", "THR", "A", 11, 5, 0, 0, 0.40, "O"),
        _atom("ATOM", 4, "OG1", "B", "THR", "A", 11, 5, 1, 0, 0.60, "O"),
    ])
    choices = conformer_choices(mda.Universe(path).atoms)
    assert choices == {("A", 10, "", "SER"): ("A",), ("A", 11, "", "THR"): ("B",)}


def test_occupancy_is_compared_per_conformer_not_per_atom(tmp_path):
    """B's mean (0.55) beats A's (0.45) although A has the single highest atom."""
    path = _write(tmp_path, [
        _atom("ATOM", 1, "CB", "A", "SER", "A", 10, 0, 0, 0, 0.60, "C"),
        _atom("ATOM", 2, "OG", "A", "SER", "A", 10, 0, 0, 1, 0.30, "O"),
        _atom("ATOM", 3, "CB", "B", "SER", "A", 10, 1, 0, 0, 0.55, "C"),
        _atom("ATOM", 4, "OG", "B", "SER", "A", 10, 1, 0, 1, 0.55, "O"),
    ])
    assert conformer_choices(mda.Universe(path).atoms) == {("A", 10, "", "SER"): ("B",)}


def test_a_structure_without_altlocs_is_returned_untouched(tmp_path):
    path = _write(tmp_path, [
        _atom("ATOM", 1, "CA", "", "GLY", "A", 1, 0, 0, 0, 1.00, "C"),
        _atom("ATOM", 2, "CA", "", "GLY", "A", 2, 3.8, 0, 0, 1.00, "C"),
    ])
    atoms = mda.Universe(path).atoms
    assert not has_alternate_conformers(atoms)
    assert highest_occupancy_atoms(atoms) is atoms


# ===========================================================================
# Real data, against an independent reading of the file
# ===========================================================================

def _oracle(path):
    """Kept labels per residue (all tied for highest mean occupancy), from PDB columns."""
    groups = collections.defaultdict(collections.OrderedDict)
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")) and line[16] != " ":
            key = (line[21], int(line[22:26]), line[26].strip(), line[17:20].strip())
            groups[key].setdefault(line[16], []).append(float(line[54:60]))
    out = {}
    for key, by in groups.items():
        mean = {c: sum(v) / len(v) for c, v in by.items()}
        top = max(mean.values())
        out[key] = tuple(c for c in by if top - mean[c] <= 1e-6)
    return out


@pytest.mark.skipif(not os.path.isfile(SITE), reason="5HDE fixture absent")
def test_real_5hde_choices_match_an_independent_reading():
    assert conformer_choices(mda.Universe(SITE).atoms) == _oracle(SITE)


@pytest.mark.skipif(not os.path.isfile(SITE), reason="5HDE fixture absent")
def test_real_5hde_has_a_residue_where_first_listed_is_the_minority():
    """The case that rules out 'first conformer wins' is present in real data."""
    firsts = {}
    for line in open(SITE):
        if line.startswith(("ATOM", "HETATM")) and line[16] != " ":
            key = (line[21], int(line[22:26]), line[26].strip(), line[17:20].strip())
            firsts.setdefault(key, line[16])
    assert any(firsts[key] not in labels for key, labels in _oracle(SITE).items())


@pytest.mark.skipif(not os.path.isfile(SITE), reason="5HDE fixture absent")
def test_real_5hde_csp231_keeps_its_phosphate():
    """CSP231: A at 0.70 carries the phosphate; B at 0.30 is dropped."""
    kept = highest_occupancy_atoms(mda.Universe(SITE).atoms).select_atoms("resname CSP")
    assert set(kept.altLocs) <= {"", "A"}
    assert {"P", "O1P", "O2P", "O3P", "SG"} <= set(kept.names)
    assert list(kept.names).count("SG") == 1


# ===========================================================================
# In the network builder
# ===========================================================================

def test_one_water_contact_is_one_edge_when_a_conformer_dominates(tmp_path):
    """A water within H-bond distance of both OG positions of one Ser, 0.35/0.65."""
    path = _write(tmp_path, [
        _atom("ATOM", 1, "N", "", "SER", "A", 10, -1.5, 0.0, 0.0, 1.00, "N"),
        _atom("ATOM", 2, "CA", "", "SER", "A", 10, 0.0, 0.0, 0.0, 1.00, "C"),
        _atom("ATOM", 3, "OG", "A", "SER", "A", 10, 1.0, 1.0, 0.0, 0.35, "O"),
        _atom("ATOM", 4, "OG", "B", "SER", "A", 10, 1.0, -1.0, 0.0, 0.65, "O"),
        _atom("HETATM", 5, "O", "", "HOH", "A", 500, 3.0, 0.0, 0.0, 1.00, "O"),
    ])
    net = _network(path)
    assert len([a for a in net.protein_atoms if a.name == "OG"]) == 1
    assert len([c for c in net.connections if c[3] == "WAT-PROT"]) == 1


def test_a_contact_made_only_by_a_tied_conformer_survives(tmp_path):
    """The barnase Asp54 case: the water reaches only conformer B, tied at 0.50."""
    from WatCon.evolutionary import water_residue_contacts

    path = _write(tmp_path, [
        _atom("ATOM", 1, "N", "", "ASP", "A", 54, 0.0, 0.0, 0.0, 1.00, "N"),
        _atom("ATOM", 2, "CA", "", "ASP", "A", 54, 1.5, 0.0, 0.0, 1.00, "C"),
        _atom("ATOM", 3, "OD1", "A", "ASP", "A", 54, 3.0, 5.0, 0.0, 0.50, "O"),
        _atom("ATOM", 4, "OD1", "B", "ASP", "A", 54, 3.0, -5.0, 0.0, 0.50, "O"),
        _atom("HETATM", 5, "O", "", "HOH", "A", 600, 3.0, -7.8, 0.0, 1.00, "O"),
    ])
    net = _network(path)
    contacted = {key for residues in water_residue_contacts(net).values() for key in residues}
    assert ("A", 54, None) in contacted


def test_water_position_follows_the_more_populated_oxygen(tmp_path):
    """Regression: the water loop read mol.atoms and brought filtered positions back."""
    path = _write(tmp_path, [
        # A polar protein atom far from the water: the builder needs at least one.
        _atom("ATOM", 1, "N", "", "GLY", "A", 1, 20.0, 20.0, 20.0, 1.00, "N"),
        _atom("HETATM", 2, "O", "A", "HOH", "A", 700, 0.0, 0.0, 0.0, 0.30, "O"),
        _atom("HETATM", 3, "O", "B", "HOH", "A", 700, 2.0, 0.0, 0.0, 0.70, "O"),
    ])
    net = _network(path)
    assert len(net.water_molecules) == 1
    assert tuple(round(float(c), 3) for c in net.water_molecules[0].O.coordinates) == (2.0, 0.0, 0.0)


def test_a_tied_water_uses_its_first_surviving_oxygen(tmp_path):
    """Documented limit: a water node holds one oxygen, so a tie keeps the first."""
    path = _write(tmp_path, [
        # A polar protein atom far from the water: the builder needs at least one.
        _atom("ATOM", 1, "N", "", "GLY", "A", 1, 20.0, 20.0, 20.0, 1.00, "N"),
        _atom("HETATM", 2, "O", "A", "HOH", "A", 700, 0.0, 0.0, 0.0, 0.50, "O"),
        _atom("HETATM", 3, "O", "B", "HOH", "A", 700, 2.0, 0.0, 0.0, 0.50, "O"),
    ])
    net = _network(path)
    assert len(net.water_molecules) == 1
    assert tuple(round(float(c), 3) for c in net.water_molecules[0].O.coordinates) == (0.0, 0.0, 0.0)

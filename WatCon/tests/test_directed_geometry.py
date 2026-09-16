"""The directed path: hydrogen-bond geometry, on real hydrogens.

WatCon can build either an oxygen-only network, where an edge means two polar
atoms are close, or a *directed* one, where an edge means a hydrogen points from
a donor at an acceptor within an angle criterion. The second is the more
chemically meaningful of the two and, until this module, was covered by no test
at all -- because crystal structures have no hydrogens, so nothing in the
fixtures could exercise it.

The fixture here is the active site of the Zenodo MD system (Amber ff14SB),
which is the only real hydrogen-bearing PTP1B available. Its numbering is the MD
system's own -- **resid 214 is PTP1B's Cys215** -- and that offset is stated
wherever the fixture is used rather than silently applied.
"""

from __future__ import annotations

import collections
import os

import pytest

pytest.importorskip("MDAnalysis", reason="the directed path needs MDAnalysis")

from WatCon import generate_dynamic_networks as gdn

PACKAGE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ENSEMBLE = os.path.join(PACKAGE, "data", "examples", "ptp1b_ensemble")
MD_SITE = "md_active_site_h.pdb"

pytestmark = pytest.mark.skipif(
    not os.path.isfile(os.path.join(ENSEMBLE, MD_SITE)),
    reason="MD active-site fixture absent")

CHEAP = {
    "density": "on", "connected_components": "on", "interaction_counts": "on",
    "per_residue_interactions": "on", "characteristic_path_length": "off",
    "graph_entropy": "on", "clustering_coefficient": "on", "shortest_path": "off",
}


def _run(directed, angle=None, distance=3.0):
    _metrics, networks, _coords = gdn.initialize_network(
        topology_file=MD_SITE, trajectory_file=MD_SITE,
        structure_directory=ENSEMBLE, network_type="water-protein",
        water_name="WAT", msa_indexing=False, return_network=True,
        num_workers=1, include_hydrogens=directed, angle_criteria=angle,
        max_distance=distance, analysis_conditions=CHEAP)
    return networks[0]


@pytest.fixture(scope="module")
def oxygen():
    return _run(False, distance=3.0)


@pytest.fixture(scope="module")
def directed():
    return _run(True, angle=150, distance=2.5)


# ===========================================================================
# The fixture is what it claims to be
# ===========================================================================

def test_the_system_has_real_hydrogens(oxygen):
    """Without these there is nothing for the directed path to read."""
    text = open(os.path.join(ENSEMBLE, MD_SITE), encoding="utf-8").read()
    hydrogens = sum(1 for line in text.splitlines()
                    if line.startswith(("ATOM", "HETATM"))
                    and line[12:16].strip().startswith("H"))
    assert hydrogens > 300


def test_waters_are_named_wat_and_carry_their_hydrogens(directed):
    assert len(directed.water_molecules) == 6
    assert all(water.H1 is not None and water.H2 is not None
               for water in directed.water_molecules)


# ===========================================================================
# The two paths differ in the way they should
# ===========================================================================

def test_the_oxygen_path_is_undirected_and_the_other_is_not(oxygen, directed):
    assert not oxygen.graph.is_directed()
    assert directed.graph.is_directed()


def test_both_paths_find_hydrogen_bonds(oxygen, directed):
    assert oxygen.graph.number_of_edges() > 0
    assert directed.graph.number_of_edges() > 0


def test_the_directed_path_sees_more_protein_atoms(oxygen, directed):
    """It has to: it keeps the hydrogens as well as the polar heavy atoms."""
    assert len(directed.protein_atoms) > len(oxygen.protein_atoms)


def test_every_connection_is_classified(oxygen, directed):
    for network in (oxygen, directed):
        kinds = collections.Counter(c[3] for c in network.connections)
        assert set(kinds) <= {"WAT-WAT", "WAT-PROT"}
        assert sum(kinds.values()) == len(network.connections)


# ===========================================================================
# The angle criterion actually does something
# ===========================================================================

def test_a_stricter_angle_cannot_add_edges():
    """Monotonic by definition: tightening the criterion only removes bonds."""
    loose = _run(True, angle=90, distance=2.5).graph.number_of_edges()
    strict = _run(True, angle=170, distance=2.5).graph.number_of_edges()
    assert strict <= loose


def test_a_shorter_cutoff_cannot_add_edges():
    far = _run(True, angle=150, distance=2.5).graph.number_of_edges()
    near = _run(True, angle=150, distance=1.8).graph.number_of_edges()
    assert near <= far

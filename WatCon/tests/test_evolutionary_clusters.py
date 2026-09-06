"""Cross-structure join: conserved water SITES to residue conservation.

This is the piece that makes the original scientific question answerable --
*across a family, are structurally conserved water sites lined by
evolutionarily conserved residues?*  Until now conservation was attached
per structure and never reached WatCon's cluster centres.

Synthetic geometry is used where exact membership must be asserted; real
structures (7O7W, 1BRS) are used for the end-to-end behaviour.

NOT tested here, and deliberately so: whether conserved water sites *are*
lined by conserved residues.  That is the research question, and answering it
needs a real protein family with ConSurf data for every member.  We hold four
ConSurf datasets covering two proteins.  Mechanical correctness is verified;
the biology is not.
"""

from __future__ import annotations

import csv

import pytest

from WatCon.consurf import parse_consurf
from WatCon.evolutionary import (
    REPORT_COLUMNS,
    ClusterConservation,
    ConservationMap,
    aggregate_site,
    conservation_of_clusters,
    conservation_summary,
    water_residue_contacts,
    write_conservation_report,
)

from .conftest import O7W_45, bundle_member

pytest.importorskip("MDAnalysis", reason="cluster join tests need MDAnalysis")


# ---------------------------------------------------------------------------
# Minimal fakes, for geometry we control exactly
# ---------------------------------------------------------------------------

class FakeAtom:
    def __init__(self, index, chain, resid, conservation, icode=None):
        self.index = index
        self.chain = chain
        self.resid = resid
        self.icode = icode
        self.evolutionary = conservation


class FakeOxygen:
    def __init__(self, index, xyz):
        self.index = index
        self.coordinates = xyz


class FakeWater:
    def __init__(self, index, xyz):
        self.O = FakeOxygen(index, xyz)
        self.evolutionary = None


class FakeNetwork:
    """Enough of WaterNetwork for the join: waters, atoms, connections."""

    def __init__(self, waters, atoms, connections):
        self.water_molecules = waters
        self.protein_atoms = atoms
        self.connections = connections


def conservation(score, grade, low_confidence=False):
    from WatCon.evolutionary import ResidueConservation

    return ResidueConservation(
        score=score, grade=grade, low_confidence=low_confidence,
        msa_present=50, msa_total=50,
        confidence_lower=None, confidence_upper=None,
        buried_exposed=None, functional_structural=None,
        consurf_position=1, source="test",
    )


def wat_prot(protein_index, water_index):
    """A WAT-PROT connection tuple in WatCon's undirected 6-element shape."""
    return (protein_index, water_index, "O", "WAT-PROT", "None", "side-chain")


# ===========================================================================
# The join
# ===========================================================================

def test_single_structure_single_site():
    water = FakeWater(100, (0.0, 0.0, 0.0))
    atom = FakeAtom(1, "A", 10, conservation(-1.2, 9))
    net = FakeNetwork([water], [atom], [wat_prot(1, 100)])

    clusters = conservation_of_clusters([net], [(0.0, 0.0, 0.0)], dist_cutoff=1.5)

    assert set(clusters) == {0}
    record = clusters[0]
    assert record.occupancy == 1
    assert record.n_structures_occupied == 1
    assert record.n_structures_total == 1
    assert record.min_score == pytest.approx(-1.2)
    assert record.max_grade == 9
    assert record.n_residues == 1
    assert record.residue_keys == (("A", 10, None),)


def test_distance_cutoff_is_respected():
    near = FakeWater(100, (0.0, 0.0, 1.0))
    far = FakeWater(101, (0.0, 0.0, 9.0))
    atoms = [FakeAtom(1, "A", 10, conservation(-1.0, 9)),
             FakeAtom(2, "A", 20, conservation(0.5, 3))]
    net = FakeNetwork([near, far], atoms, [wat_prot(1, 100), wat_prot(2, 101)])

    clusters = conservation_of_clusters([net], [(0.0, 0.0, 0.0)], dist_cutoff=1.5)

    assert clusters[0].occupancy == 1
    assert clusters[0].residue_keys == (("A", 10, None),)


def test_residue_reached_by_two_waters_counts_once():
    """De-duplication across waters within a site."""
    waters = [FakeWater(100, (0.0, 0.0, 0.0)), FakeWater(101, (0.5, 0.0, 0.0))]
    atom = FakeAtom(1, "A", 10, conservation(-1.0, 9))
    net = FakeNetwork(waters, [atom], [wat_prot(1, 100), wat_prot(1, 101)])

    record = conservation_of_clusters([net], [(0.0, 0.0, 0.0)], 1.5)[0]

    assert record.occupancy == 2          # two waters occupy the site
    assert record.n_residues == 1         # but one residue lines it
    assert record.residue_keys == (("A", 10, None),)


def test_residue_shared_across_structures_counts_once():
    """De-duplication ACROSS structures, not just within one."""
    nets = []
    for offset in (0.0, 0.3):
        water = FakeWater(100, (offset, 0.0, 0.0))
        atom = FakeAtom(1, "A", 10, conservation(-1.0, 9))
        nets.append(FakeNetwork([water], [atom], [wat_prot(1, 100)]))

    record = conservation_of_clusters(nets, [(0.0, 0.0, 0.0)], 1.5)[0]

    assert record.occupancy == 2
    assert record.n_structures_occupied == 2
    assert record.n_structures_total == 2
    assert record.n_residues == 1
    assert record.occupancy_fraction == pytest.approx(1.0)


def test_site_occupied_in_a_subset_of_structures():
    occupied = FakeNetwork(
        [FakeWater(100, (0.0, 0.0, 0.0))],
        [FakeAtom(1, "A", 10, conservation(-1.0, 9))],
        [wat_prot(1, 100)],
    )
    empty = FakeNetwork([FakeWater(100, (20.0, 0.0, 0.0))], [], [])
    also_occupied = FakeNetwork(
        [FakeWater(100, (0.2, 0.0, 0.0))],
        [FakeAtom(1, "A", 11, conservation(0.0, 5))],
        [wat_prot(1, 100)],
    )

    record = conservation_of_clusters(
        [occupied, empty, also_occupied], [(0.0, 0.0, 0.0)], 1.5
    )[0]

    assert record.n_structures_occupied == 2
    assert record.n_structures_total == 3
    assert record.occupancy_fraction == pytest.approx(2 / 3)
    assert record.n_residues == 2


def test_unscored_residues_are_counted_not_dropped():
    water = FakeWater(100, (0.0, 0.0, 0.0))
    atoms = [FakeAtom(1, "A", 10, conservation(-1.0, 9)),
             FakeAtom(2, "A", 20, None)]
    net = FakeNetwork([water], atoms, [wat_prot(1, 100), wat_prot(2, 100)])

    record = conservation_of_clusters([net], [(0.0, 0.0, 0.0)], 1.5)[0]

    assert record.n_residues == 1
    assert record.n_unscored == 1
    assert record.min_score == pytest.approx(-1.0)


def test_site_with_no_scored_residue_reports_none_not_zero():
    """None means 'no data'.  Zero would read as 'not conserved'."""
    water = FakeWater(100, (0.0, 0.0, 0.0))
    net = FakeNetwork([water], [FakeAtom(1, "A", 10, None)], [wat_prot(1, 100)])

    record = conservation_of_clusters([net], [(0.0, 0.0, 0.0)], 1.5)[0]

    assert record.min_score is None
    assert record.mean_score is None
    assert record.max_grade is None
    assert record.has_conservation is False
    assert record.n_unscored == 1
    assert record.occupancy == 1          # the site IS occupied, just unscored


def test_empty_site_still_reported():
    net = FakeNetwork([FakeWater(100, (50.0, 0.0, 0.0))], [], [])
    record = conservation_of_clusters([net], [(0.0, 0.0, 0.0)], 1.5)[0]

    assert record.occupancy == 0
    assert record.n_structures_occupied == 0
    assert record.occupancy_fraction == 0.0
    assert record.has_conservation is False


def test_low_confidence_counted_but_kept():
    water = FakeWater(100, (0.0, 0.0, 0.0))
    atoms = [FakeAtom(1, "A", 10, conservation(-1.0, 9, low_confidence=True)),
             FakeAtom(2, "A", 20, conservation(0.0, 5))]
    net = FakeNetwork([water], atoms, [wat_prot(1, 100), wat_prot(2, 100)])

    record = conservation_of_clusters([net], [(0.0, 0.0, 0.0)], 1.5)[0]

    assert record.n_residues == 2
    assert record.n_low_confidence == 1
    assert record.min_score == pytest.approx(-1.0)


def test_chain_is_part_of_the_residue_key():
    water = FakeWater(100, (0.0, 0.0, 0.0))
    atoms = [FakeAtom(1, "A", 10, conservation(-1.0, 9)),
             FakeAtom(2, "B", 10, conservation(0.5, 3))]
    net = FakeNetwork([water], atoms, [wat_prot(1, 100), wat_prot(2, 100)])

    record = conservation_of_clusters([net], [(0.0, 0.0, 0.0)], 1.5)[0]

    assert record.n_residues == 2
    assert set(record.residue_keys) == {("A", 10, None), ("B", 10, None)}


def test_residue_zero_participates():
    """Regression: residue number 0 is falsy and was once silently dropped."""
    water = FakeWater(100, (0.0, 0.0, 0.0))
    atom = FakeAtom(1, "A", 0, conservation(-0.8, 8))
    net = FakeNetwork([water], [atom], [wat_prot(1, 100)])

    record = conservation_of_clusters([net], [(0.0, 0.0, 0.0)], 1.5)[0]

    assert record.n_residues == 1
    assert record.residue_keys == (("A", 0, None),)


def test_negative_residue_numbers_participate():
    water = FakeWater(100, (0.0, 0.0, 0.0))
    atom = FakeAtom(1, "A", -5, conservation(-0.5, 7))
    net = FakeNetwork([water], [atom], [wat_prot(1, 100)])

    record = conservation_of_clusters([net], [(0.0, 0.0, 0.0)], 1.5)[0]
    assert record.residue_keys == (("A", -5, None),)


def test_centers_accepted_as_dict_or_array():
    """WatCon passes both shapes around; the join must take either."""
    water = FakeWater(100, (0.0, 0.0, 0.0))
    net = FakeNetwork([water], [FakeAtom(1, "A", 10, conservation(-1.0, 9))],
                      [wat_prot(1, 100)])

    as_list = conservation_of_clusters([net], [(0.0, 0.0, 0.0)], 1.5)
    as_dict = conservation_of_clusters([net], {7: (0.0, 0.0, 0.0)}, 1.5)

    assert set(as_list) == {0}
    assert set(as_dict) == {7}
    assert as_dict[7].min_score == as_list[0].min_score


# ===========================================================================
# Site-level aggregation
# ===========================================================================

def test_aggregate_site_over_all_waters():
    waters = [FakeWater(100, (0.0, 0.0, 0.0)), FakeWater(101, (5.0, 0.0, 0.0))]
    atoms = [FakeAtom(1, "A", 10, conservation(-1.0, 9)),
             FakeAtom(2, "A", 20, conservation(0.5, 3))]
    net = FakeNetwork(waters, atoms, [wat_prot(1, 100), wat_prot(2, 101)])

    result = aggregate_site(net)
    assert result.n_residues == 2
    assert result.min_score == pytest.approx(-1.0)
    assert result.max_grade == 9


def test_aggregate_site_over_a_subset():
    waters = [FakeWater(100, (0.0, 0.0, 0.0)), FakeWater(101, (5.0, 0.0, 0.0))]
    atoms = [FakeAtom(1, "A", 10, conservation(-1.0, 9)),
             FakeAtom(2, "A", 20, conservation(0.5, 3))]
    net = FakeNetwork(waters, atoms, [wat_prot(1, 100), wat_prot(2, 101)])

    result = aggregate_site(net, waters=[waters[0]])
    assert result.n_residues == 1
    assert result.max_grade == 9


def test_aggregate_site_with_nothing_scored():
    net = FakeNetwork([FakeWater(100, (0.0, 0.0, 0.0))],
                      [FakeAtom(1, "A", 10, None)], [wat_prot(1, 100)])
    assert aggregate_site(net) is None


def test_water_residue_contacts_handles_missing_connections():
    net = FakeNetwork([], [], None)
    assert water_residue_contacts(net) == {}


# ===========================================================================
# Report
# ===========================================================================

def test_report_header_and_rows(tmp_path):
    water = FakeWater(100, (0.0, 0.0, 0.0))
    net = FakeNetwork([water], [FakeAtom(1, "A", 10, conservation(-1.0, 9))],
                      [wat_prot(1, 100)])
    clusters = conservation_of_clusters([net], [(0.0, 0.0, 0.0), (50.0, 0.0, 0.0)], 1.5)

    path = tmp_path / "report.csv"
    n = write_conservation_report(clusters, path)

    assert n == 2
    rows = list(csv.DictReader(open(path)))
    assert list(rows[0]) == REPORT_COLUMNS
    assert rows[0]["evo_max_grade"] == "9"
    # The unoccupied site reports NA, not 0 -- absence of data must survive.
    assert rows[1]["evo_min_score"] == "NA"
    assert rows[1]["evo_max_grade"] == "NA"


def test_report_carries_both_measures_side_by_side():
    """Structural and evolutionary conservation must both be present, unblended."""
    structural = [c for c in REPORT_COLUMNS if not c.startswith("evo_")]
    evolutionary = [c for c in REPORT_COLUMNS if c.startswith("evo_")]

    assert "occupancy" in structural
    assert "evo_min_score" in evolutionary

    # No column blends the two measures into one number.
    blended = [c for c in REPORT_COLUMNS
               if any(w in c for w in ("combined", "composite", "weighted", "overall"))]
    assert blended == []

    # Every column belongs to exactly one side; the evo_ prefix is the boundary.
    assert set(structural) & set(evolutionary) == set()


def test_conservation_summary_is_sorted():
    net = FakeNetwork([FakeWater(100, (0.0, 0.0, 0.0))],
                      [FakeAtom(1, "A", 10, conservation(-1.0, 9))],
                      [wat_prot(1, 100)])
    clusters = conservation_of_clusters(
        [net], {5: (0.0, 0.0, 0.0), 2: (9.0, 0.0, 0.0)}, 1.5
    )
    rows = conservation_summary(clusters)
    assert [r["cluster_id"] for r in rows] == [2, 5]


# ===========================================================================
# Real structures, end to end
# ===========================================================================

@pytest.fixture(scope="module")
def o7w_network(tmp_path_factory):
    from WatCon.generate_static_networks import extract_objects

    path = tmp_path_factory.mktemp("clusters") / "7O7W.pdb"
    path.write_bytes(
        bundle_member(
            "1788241750_ConSurf.tar.gz", "7O7W_A_ATOMS_section_With_ConSurf.pdb"
        )
    )
    cmap = ConservationMap.build(parse_consurf(O7W_45))
    return extract_objects(
        str(path), "water-protein", None, None, False, 8.0, None, None,
        conservation_map=cmap,
    )


def test_real_structure_self_consistency(o7w_network):
    """Every real water used as its own centre must be found at that centre."""
    centers = [w.O.coordinates for w in o7w_network.water_molecules[:20]]
    clusters = conservation_of_clusters([o7w_network], centers, dist_cutoff=1.5)

    assert len(clusters) == 20
    assert all(c.n_structures_occupied == 1 for c in clusters.values())
    assert all(c.occupancy >= 1 for c in clusters.values())


def test_real_structure_residue_keys_are_plain_ints(o7w_network):
    """MDAnalysis yields numpy ints; they must not leak into the data model."""
    centers = [w.O.coordinates for w in o7w_network.water_molecules[:10]]
    clusters = conservation_of_clusters([o7w_network], centers, 1.5)

    keys = [k for c in clusters.values() for k in c.residue_keys]
    assert keys
    assert all(type(k[1]) is int for k in keys)


def test_real_structure_produces_conservation(o7w_network):
    centers = [w.O.coordinates for w in o7w_network.water_molecules[:20]]
    clusters = conservation_of_clusters([o7w_network], centers, 1.5)

    scored = [c for c in clusters.values() if c.has_conservation]
    assert scored, "7O7W is fully covered by ConSurf; some sites must score"
    for record in scored:
        assert 1 <= record.max_grade <= 9
        assert record.n_residues >= 1
        assert record.min_score <= record.mean_score + 1e-9


def test_far_away_centre_is_empty_on_a_real_structure(o7w_network):
    clusters = conservation_of_clusters([o7w_network], [(999.0, 999.0, 999.0)], 1.5)
    assert clusters[0].occupancy == 0
    assert clusters[0].has_conservation is False

"""ConSurf conservation attached to WatCon residues, waters and graphs.

Real-data tests use the ConSurf bundles already in the repository:

* **7O7W_A** -- expression tag numbered -5..0, a chromophore ConSurf omitted,
  41 fully-conserved positions, shallow (45-sequence) MSA.
* **1BRS_A** -- six chains, ConSurf run on chain A only.
* **P00648** x2 -- the same protein at MSA depth 150 and 50, for reproducibility.
* **1AKI** -- no ConSurf data at all; the "everything unscored" path.

The unit-level parts (``ConservationMap``, ``aggregate_water``) need no
MDAnalysis; the end-to-end parts do and are skipped without it.
"""

from __future__ import annotations

import io

import pytest

from WatCon.consurf import LookupStatus, parse_consurf
from WatCon.evolutionary import (
    ConservationError,
    ConservationMap,
    ResidueConservation,
    WaterConservation,
    aggregate_water,
    find_consurf_file,
    load_conservation,
)

from .conftest import BRS_150, FIXTURES, O7W_45, P00648_50, P00648_150, bundle_member


@pytest.fixture(scope="module")
def o7w_map():
    return ConservationMap.build(parse_consurf(O7W_45))


@pytest.fixture(scope="module")
def brs_map():
    return ConservationMap.build(parse_consurf(BRS_150))


# ===========================================================================
# ConservationMap
# ===========================================================================

def test_map_built_from_real_file(o7w_map):
    assert len(o7w_map) == 236
    assert o7w_map.chains() == ["A"]


def test_lookup_by_identity(o7w_map):
    """The expression tag: negative residue numbers must resolve."""
    tag = o7w_map.for_residue("A", -5)
    assert tag is not None
    assert tag.grade == 7
    assert tag.score == pytest.approx(-0.509)
    assert tag.low_confidence is True
    assert tag.msa_present == 9 and tag.msa_total == 45
    assert tag.consurf_position == 8


def test_lookup_zero_residue_number(o7w_map):
    """resid 0 is falsy -- it must still be found."""
    assert o7w_map.for_residue("A", 0) is not None
    assert ("A", 0, None) in o7w_map


def test_offset_numbering_lookup():
    """P00648 is UniProt-numbered: ConSurf POS 1 is residue 48."""
    cmap = ConservationMap.build(parse_consurf(P00648_150))
    residue = cmap.for_residue("A", 48)
    assert residue is not None
    assert residue.consurf_position == 1


def test_unscored_residue_returns_none_not_a_guess(o7w_map):
    """7O7W's chromophore PIA:68:A has no ConSurf record."""
    assert o7w_map.for_residue("A", 68) is None
    assert o7w_map.status("A", 68) is LookupStatus.NO_SCORE_IN_FILE


def test_out_of_range_is_distinct_from_unscored(o7w_map):
    assert o7w_map.status("A", 9999) is LookupStatus.NOT_IN_FILE
    assert o7w_map.status("A", 68) is LookupStatus.NO_SCORE_IN_FILE
    assert o7w_map.status("A", -5) is LookupStatus.SCORED


def test_chain_isolation(brs_map):
    """1BRS was run on chain A; chain B must not borrow chain A's scores."""
    assert brs_map.for_residue("A", 10) is not None
    assert brs_map.for_residue("B", 10) is None
    assert brs_map.status("B", 10) is LookupStatus.NOT_IN_FILE


def test_chain_map_translation():
    """A ConSurf run labelled 'A' can be attached to structure chain 'X'."""
    result = parse_consurf(BRS_150)
    translated = ConservationMap.build(result, chain_map={"A": "X"})
    assert translated.for_residue("X", 10) is not None
    assert translated.for_residue("A", 10) is None
    assert translated.chains() == ["X"]


def test_sign_convention_is_preserved(o7w_map):
    """More negative score = more conserved; grade 9 = most conserved."""
    values = [o7w_map.for_residue(*k) for k in o7w_map.keys()]
    most_conserved = min(values, key=lambda c: c.score)
    least_conserved = max(values, key=lambda c: c.score)
    assert most_conserved.grade == 9
    assert least_conserved.grade == 1
    assert most_conserved.is_conserved
    assert not least_conserved.is_conserved


def test_low_confidence_is_kept_not_filtered(o7w_map):
    flagged = [
        o7w_map.for_residue(*k)
        for k in o7w_map.keys()
        if o7w_map.for_residue(*k).low_confidence
    ]
    assert flagged, "7O7W has low-confidence positions"
    assert all(c.score is not None and c.grade in range(1, 10) for c in flagged)


def test_msa_fraction(o7w_map):
    tag = o7w_map.for_residue("A", -5)
    assert tag.msa_fraction == pytest.approx(9 / 45)


# ===========================================================================
# Coverage
# ===========================================================================

def test_coverage_against_the_real_structure(o7w_map, tmp_path):
    from WatCon.residue_index import residues_from_pdb_file

    path = tmp_path / "7O7W.pdb"
    path.write_bytes(
        bundle_member(
            "1788241750_ConSurf.tar.gz", "7O7W_A_ATOMS_section_With_ConSurf.pdb"
        )
    )
    report = o7w_map.coverage(residues_from_pdb_file(path))

    assert report.matched == 236
    assert report.fraction == pytest.approx(1.0)
    assert "A" in report.per_chain


def test_coverage_detects_a_chain_mismatch(brs_map, tmp_path):
    """Six-chain structure, ConSurf run on one chain -- coverage must show it."""
    from WatCon.residue_index import residues_from_pdb_file

    path = tmp_path / "1BRS.pdb"
    path.write_bytes(
        bundle_member(
            "1788241897_ConSurf.tar.gz", "1BRS_A_ATOMS_section_With_ConSurf.pdb"
        )
    )
    report = brs_map.coverage(residues_from_pdb_file(path))

    assert report.matched == 108
    assert report.fraction < 0.25
    assert len(report.per_chain) == 6
    assert "108" in report.describe()


# ===========================================================================
# Water aggregation
# ===========================================================================

def _conservation(score, grade, low_confidence=False):
    return ResidueConservation(
        score=score, grade=grade, low_confidence=low_confidence,
        msa_present=50, msa_total=50,
        confidence_lower=None, confidence_upper=None,
        buried_exposed=None, functional_structural=None,
        consurf_position=1, source="test",
    )


def test_aggregate_basic_statistics():
    result = aggregate_water([_conservation(-1.0, 9), _conservation(0.5, 4)])
    assert result.min_score == pytest.approx(-1.0)
    assert result.mean_score == pytest.approx(-0.25)
    assert result.max_grade == 9
    assert result.n_residues == 2
    assert result.n_unscored == 0


def test_aggregate_counts_unscored_without_dropping_them():
    result = aggregate_water([_conservation(-1.0, 9), None, None])
    assert result.n_residues == 1
    assert result.n_unscored == 2
    assert result.min_score == pytest.approx(-1.0)


def test_aggregate_returns_none_when_nothing_is_scored():
    """None means 'no data', and must not be read as zero conservation."""
    assert aggregate_water([None, None]) is None
    assert aggregate_water([]) is None


def test_aggregate_counts_low_confidence_but_keeps_it():
    result = aggregate_water(
        [_conservation(-1.0, 9, low_confidence=True), _conservation(0.0, 5)]
    )
    assert result.n_residues == 2
    assert result.n_low_confidence == 1
    assert result.min_score == pytest.approx(-1.0)


def test_water_conservation_to_dict_uses_evo_prefix():
    result = aggregate_water([_conservation(-1.0, 9)])
    assert set(result.to_dict()) == {
        "evo_min_score", "evo_mean_score", "evo_max_grade",
        "evo_n_residues", "evo_n_low_confidence", "evo_n_unscored",
    }


# ===========================================================================
# File discovery and loading
# ===========================================================================

def test_find_consurf_file_matches_by_name():
    found = find_consurf_file("7O7W_A.pdb", FIXTURES)
    assert found is not None and "7O7W" in found.name


def test_find_consurf_file_returns_none_when_absent():
    assert find_consurf_file("1AKI.pdb", FIXTURES) is None


def test_load_conservation_disabled_returns_none():
    assert load_conservation("anything.pdb", None) is None


def test_load_conservation_strict_raises_on_missing_file():
    with pytest.raises(ConservationError, match="No ConSurf grades file"):
        load_conservation("1AKI.pdb", FIXTURES, strict=True)


def test_load_conservation_tolerant_returns_none(capsys):
    assert load_conservation("1AKI.pdb", FIXTURES, strict=False) is None
    assert "Warning" in capsys.readouterr().out


def test_load_conservation_success():
    cmap = load_conservation("7O7W_A.pdb", FIXTURES)
    assert cmap is not None and len(cmap) == 236


# ===========================================================================
# Reproducibility across ConSurf runs
# ===========================================================================

def test_same_protein_two_runs_score_the_same_residues():
    deep = ConservationMap.build(parse_consurf(P00648_150))
    shallow = ConservationMap.build(parse_consurf(P00648_50))
    assert set(deep.keys()) == set(shallow.keys())


def test_score_rank_is_reproducible_across_runs():
    """Spearman >= 0.94 between MSA depth 150 and 50 (measured 0.955)."""
    deep = ConservationMap.build(parse_consurf(P00648_150))
    shallow = ConservationMap.build(parse_consurf(P00648_50))
    keys = sorted(set(deep.keys()) & set(shallow.keys()))

    a = [deep.for_residue(*k).score for k in keys]
    b = [shallow.for_residue(*k).score for k in keys]

    def rank(values):
        order = sorted(range(len(values)), key=lambda i: values[i])
        ranks = [0] * len(values)
        for position, index in enumerate(order):
            ranks[index] = position
        return ranks

    ra, rb = rank(a), rank(b)
    n = len(keys)
    d2 = sum((ra[i] - rb[i]) ** 2 for i in range(n))
    assert 1 - 6 * d2 / (n * (n * n - 1)) >= 0.94


# ===========================================================================
# End to end through the network builders
# ===========================================================================

pytest.importorskip("MDAnalysis", reason="end-to-end tests require MDAnalysis")


@pytest.fixture(scope="module")
def o7w_network(tmp_path_factory):
    from WatCon.generate_static_networks import extract_objects

    path = tmp_path_factory.mktemp("evo") / "7O7W.pdb"
    path.write_bytes(
        bundle_member(
            "1788241750_ConSurf.tar.gz", "7O7W_A_ATOMS_section_With_ConSurf.pdb"
        )
    )
    cmap = ConservationMap.build(parse_consurf(O7W_45))
    network = extract_objects(
        str(path), "water-protein", None, None, False, 8.0, None, None,
        conservation_map=cmap,
    )
    return cmap, network


def test_every_atom_matches_the_map(o7w_network):
    cmap, network = o7w_network
    assert network.protein_atoms
    for atom in network.protein_atoms:
        assert atom.evolutionary is cmap.for_residue(atom.chain, atom.resid, atom.icode)


def test_atoms_of_a_residue_share_one_object(o7w_network):
    """Per-atom storage must not duplicate the conservation record."""
    _, network = o7w_network
    tag = [a for a in network.protein_atoms if a.resid == -5]
    assert len(tag) > 1
    assert len({id(a.evolutionary) for a in tag}) == 1


def test_waters_receive_aggregates(o7w_network):
    _, network = o7w_network
    annotated = [w for w in network.water_molecules if w.evolutionary is not None]
    assert annotated
    assert any(w.evolutionary.n_residues > 1 for w in annotated), (
        "some water should bridge more than one residue"
    )


def test_water_aggregates_deduplicate_residues(o7w_network):
    """n_residues counts DISTINCT residues, never contacts."""
    _, network = o7w_network
    atoms_by_index = {a.index: a for a in network.protein_atoms}
    water_by_oxygen = {w.O.index: w for w in network.water_molecules}

    contacts = {}
    for connection in network.connections:
        if connection[3] != "WAT-PROT":
            continue
        first, second = connection[0], connection[1]
        if first in atoms_by_index and second in water_by_oxygen:
            atom, water_index = atoms_by_index[first], second
        elif second in atoms_by_index and first in water_by_oxygen:
            atom, water_index = atoms_by_index[second], first
        else:
            continue
        contacts.setdefault(water_index, set()).add(
            (atom.chain, atom.resid, atom.icode)
        )

    for water_index, residues in contacts.items():
        water = water_by_oxygen[water_index]
        if water.evolutionary is None:
            continue
        total = water.evolutionary.n_residues + water.evolutionary.n_unscored
        assert total == len(residues)


def test_graph_attributes_agree_with_objects(o7w_network):
    _, network = o7w_network
    by_index = {a.index: a for a in network.protein_atoms}
    water_by_oxygen = {w.O.index: w for w in network.water_molecules}

    for node, data in network.graph.nodes(data=True):
        if data.get("atom_category") == "PROTEIN":
            atom = by_index[node]
            expected = None if atom.evolutionary is None else atom.evolutionary.score
            assert data["evo_score"] == expected
        elif data.get("atom_category") == "WAT":
            water = water_by_oxygen[node]
            expected = None if water.evolutionary is None else water.evolutionary.min_score
            assert data["evo_min_score"] == expected


def test_no_consurf_leaves_everything_unscored(tmp_path):
    """1AKI has no ConSurf data: nothing must be fabricated."""
    from WatCon.generate_static_networks import extract_objects
    from pathlib import Path

    pdb = Path(__file__).resolve().parent / "inputs" / "1AKI.pdb"
    network = extract_objects(
        str(pdb), "water-protein", None, None, False, 8.0, None, None,
        conservation_map=None,
    )
    assert network.protein_atoms
    assert all(a.evolutionary is None for a in network.protein_atoms)
    assert all(w.evolutionary is None for w in network.water_molecules)

    for _, data in network.graph.nodes(data=True):
        if data.get("atom_category") == "PROTEIN":
            assert data["evo_score"] is None
        elif data.get("atom_category") == "WAT":
            assert data["evo_min_score"] is None
            assert data["evo_n_residues"] == 0


def test_structural_and_evolutionary_conservation_do_not_collide(o7w_network):
    """The two 'conservations' must never share a key."""
    _, network = o7w_network
    for _, data in network.graph.nodes(data=True):
        evo_keys = {k for k in data if k.startswith("evo_")}
        assert evo_keys
        assert not (evo_keys & {"conservation", "commonality", "cluster"})


# ===========================================================================
# Classification CSV
# ===========================================================================
#
# plot_interactions_from_angles and identify_clustered_angles read this CSV by
# COLUMN POSITION, so the evolutionary fields must be appended at the end and
# the header must stay in step with the rows.

CSV_HEADER = (
    "PDB ID,Resid,MSA_Resid,Index_1,Index_2,Protein_Atom,Classification,"
    "Protein_Coords,Water_Coords,Angle_1,Angle_2,Evo_Score,Evo_Grade,Evo_LowConf"
)


def test_csv_header_matches_the_writer(o7w_network):
    """Header field count must equal the row field count."""
    from WatCon.residue_analysis import classify_waters

    _, network = o7w_network
    rows = classify_waters(network, ref1_coords=None, ref2_coords=None)
    assert rows

    key, value = next(iter(rows.items()))
    row = f"7O7W,{key},{value[0]},{value[1]},{value[2]}"
    assert len(row.split(",")) == len(CSV_HEADER.split(","))


def test_csv_header_literal_is_present_in_both_builders():
    """Guard against the header and this test drifting apart."""
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    for module in ("generate_static_networks.py", "generate_dynamic_networks.py"):
        text = (root / module).read_text(encoding="utf-8")
        assert "Angle_1,Angle_2,Evo_Score,Evo_Grade,Evo_LowConf" in text


def test_csv_reports_na_when_unscored(tmp_path):
    from WatCon.generate_static_networks import extract_objects
    from WatCon.residue_analysis import classify_waters
    from pathlib import Path

    pdb = Path(__file__).resolve().parent / "inputs" / "1AKI.pdb"
    network = extract_objects(
        str(pdb), "water-protein", None, None, False, 8.0, None, None,
        conservation_map=None,
    )
    rows = classify_waters(network, ref1_coords=None, ref2_coords=None)
    assert rows
    assert all(value[2] == "NA,NA,NA" for value in rows.values())


def test_existing_csv_columns_did_not_move():
    """Evolutionary columns are APPENDED; nothing before them shifted."""
    original = (
        "PDB ID,Resid,MSA_Resid,Index_1,Index_2,Protein_Atom,Classification,"
        "Protein_Coords,Water_Coords,Angle_1,Angle_2"
    )
    assert CSV_HEADER.startswith(original)
    assert CSV_HEADER[len(original):] == ",Evo_Score,Evo_Grade,Evo_LowConf"

"""Preparing a folder of structures into one analysable set.

This is the step that used to be manual, and the step where a dataset silently
goes wrong: the wrong chain gets picked out of a complex, a mis-numbered entry
is admitted, or a structure vanishes without anyone noticing.

The tests use the bundled barnase example, which was chosen to contain those
cases rather than to be convenient -- 1BRN carries barnase as chain **L**, 1BRS
is a barnase-barstar complex, and 1BSA and 1RNB start at different residue
numbers from the reference.
"""

from __future__ import annotations

import os

import pytest

from WatCon.prepare import (
    PreparationError,
    detect_offset,
    extract_chain_with_waters,
    find_target_chain,
    identity_by_resid,
    prepare_directory,
)
from WatCon.residue_index import StructureResidue, residues_from_pdb_file

EXAMPLES = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "data", "examples", "barnase", "structures",
)

pytestmark = pytest.mark.skipif(
    not os.path.isdir(EXAMPLES), reason="bundled example data not present"
)


@pytest.fixture(scope="module")
def reference_map():
    residues = residues_from_pdb_file(os.path.join(EXAMPLES, "1A2P.pdb"))
    chain = max(
        sorted({r.chain for r in residues}),
        key=lambda c: len([r for r in residues if r.chain == c]),
    )
    return {r.resid: r.one_letter for r in residues if r.chain == chain}


def residues(*items):
    return [StructureResidue("A", resid, None, "ALA", letter)
            for resid, letter in items]


# ===========================================================================
# Identity, keyed on residue number
# ===========================================================================

def test_identity_is_keyed_on_residue_number():
    reference = {3: "V", 4: "I", 5: "N"}
    assert identity_by_resid(residues((3, "V"), (4, "I"), (5, "N")), reference) == 1.0
    # Same letters, different numbers -- must NOT match.
    assert identity_by_resid(residues((1, "V"), (2, "I"), (3, "N")), reference) < 0.5


def test_a_fragment_cannot_score_highly_by_agreeing_on_a_few_residues():
    """The denominator is the reference length, deliberately."""
    reference = {i: "A" for i in range(1, 101)}
    assert identity_by_resid(residues((1, "A"), (2, "A")), reference) == pytest.approx(0.02)


def test_empty_reference_is_zero_not_an_error():
    assert identity_by_resid(residues((1, "A")), {}) == 0.0


def test_offset_is_detected_but_never_applied():
    reference = {3: "V", 4: "I", 5: "N", 6: "T"}
    shifted = residues((1, "V"), (2, "I"), (3, "N"), (4, "T"))

    offset, score = detect_offset(shifted, reference)
    assert offset == 2
    assert score == pytest.approx(1.0)
    # The input is untouched: detection is a diagnosis, not a correction.
    assert [r.resid for r in shifted] == [1, 2, 3, 4]


# ===========================================================================
# Chain selection, on real complexes
# ===========================================================================

def test_barnase_is_found_as_chain_L_in_1BRN(reference_map):
    """The protein of interest is not reliably chain 'A'."""
    chain, score, note = find_target_chain(
        os.path.join(EXAMPLES, "1BRN.pdb"), reference_map
    )
    assert chain == "L"
    assert score > 0.95
    assert note == ""


def test_barstar_does_not_win_in_the_1BRS_complex(reference_map):
    chain, score, _ = find_target_chain(
        os.path.join(EXAMPLES, "1BRS.pdb"), reference_map
    )
    assert score > 0.95
    residues_here = residues_from_pdb_file(os.path.join(EXAMPLES, "1BRS.pdb"))
    chosen = [r for r in residues_here if r.chain == chain]
    assert len(chosen) > 100, "barstar is ~89 residues; barnase ~110"


def test_structures_numbered_from_a_different_start_still_match(reference_map):
    """1BSA starts at residue 4 and 1RNB at 2, against the reference's 3.

    Compared positionally these score 0.03-0.05 and would be thrown away.
    """
    for pdb_id in ("1BSA", "1RNB"):
        chain, score, _ = find_target_chain(
            os.path.join(EXAMPLES, pdb_id + ".pdb"), reference_map
        )
        assert chain is not None, pdb_id
        assert score > 0.95, pdb_id


def test_an_unrelated_structure_is_rejected_with_a_score(reference_map):
    chain, score, _ = find_target_chain(
        os.path.join(EXAMPLES, "1A2P.pdb"), {1: "W", 2: "W", 3: "W"}
    )
    assert chain is None
    assert score < 0.8


# ===========================================================================
# Extraction
# ===========================================================================

def test_waters_are_selected_by_proximity_not_chain_label(reference_map):
    """Crystallographic waters are often deposited under their own chain."""
    protein, waters = extract_chain_with_waters(
        os.path.join(EXAMPLES, "1A2P.pdb"), "A"
    )
    assert protein
    assert waters, "selecting waters by chain label would find none here"


def test_water_cutoff_is_respected():
    path = os.path.join(EXAMPLES, "1A2P.pdb")
    _, near = extract_chain_with_waters(path, "A", water_cutoff=3.0)
    _, far = extract_chain_with_waters(path, "A", water_cutoff=8.0)
    assert len(near) < len(far)


# ===========================================================================
# The whole directory
# ===========================================================================

@pytest.fixture(scope="module")
def prepared(tmp_path_factory):
    out = tmp_path_factory.mktemp("prepared")
    report = prepare_directory(
        EXAMPLES, str(out / "structures"), reference="1A2P", verbose=False
    )
    return report, out / "structures"


def test_every_example_structure_prepares(prepared):
    report, out_dir = prepared
    assert len(report.prepared) == 6
    assert not report.rejected
    assert len(os.listdir(out_dir)) == 6


def test_everything_lands_in_one_frame(prepared):
    """Superposition is the point; assert it actually happened."""
    report, _ = prepared
    for outcome in report.prepared:
        assert outcome.rmsd is not None
        assert outcome.rmsd < 1.5, outcome.pdb_id
        assert outcome.n_shared > 100


def test_waters_travel_with_their_protein(prepared):
    report, out_dir = prepared
    assert report.n_waters > 500
    for outcome in report.prepared:
        assert outcome.n_waters > 0, outcome.pdb_id


def test_chains_are_relabelled_so_one_consurf_file_serves_all(prepared):
    _, out_dir = prepared
    for name in os.listdir(out_dir):
        chains = set()
        for line in open(os.path.join(out_dir, name)):
            if line.startswith(("ATOM", "HETATM")) and len(line) > 21:
                chains.add(line[21])
        assert chains == {"A"}, name


def test_the_reference_superposes_onto_itself_exactly(prepared):
    report, _ = prepared
    reference = next(o for o in report.prepared if o.pdb_id == "1A2P")
    assert reference.rmsd == pytest.approx(0.0, abs=1e-6)


def test_output_directory_is_cleared_not_merged(tmp_path):
    """A structure rejected on this run must not survive from a previous one.

    This happened for real: a stale copy of a mis-numbered structure reached the
    analysis and was caught only by the ConSurf identity check.
    """
    out = tmp_path / "out"
    out.mkdir()
    stale = out / "STALE.pdb"
    stale.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000\n")

    prepare_directory(EXAMPLES, str(out), reference="1A2P", verbose=False)
    assert not stale.exists()


def test_reference_selection_is_deterministic(tmp_path):
    """1A2P has three identical chains; the choice must not vary by run.

    ``max`` over a set returns whichever equal-length chain iteration reaches
    first, and set order for strings depends on PYTHONHASHSEED.
    """
    chosen = set()
    for i in range(3):
        report = prepare_directory(
            EXAMPLES, str(tmp_path / ("run%d" % i)), reference="1A2P", verbose=False
        )
        chosen.add(report.reference_chain)
    assert len(chosen) == 1


# ===========================================================================
# Failing loudly
# ===========================================================================

def test_a_missing_directory_is_refused(tmp_path):
    with pytest.raises(PreparationError, match="no such directory"):
        prepare_directory(str(tmp_path / "nope"), str(tmp_path / "out"))


def test_an_empty_directory_is_refused(tmp_path):
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(PreparationError, match="no .pdb files"):
        prepare_directory(str(empty), str(tmp_path / "out"))


def test_an_unknown_reference_is_refused(tmp_path):
    with pytest.raises(PreparationError, match="not among the structures"):
        prepare_directory(EXAMPLES, str(tmp_path / "out"), reference="9XYZ")


def test_a_dataset_that_matches_nothing_fails_rather_than_returning_empty(tmp_path):
    """An empty output directory would look like success to the next step."""
    lonely = tmp_path / "lonely"
    lonely.mkdir()
    with open(lonely / "a.pdb", "w") as handle:
        for i in range(5):
            handle.write(
                "ATOM  %5d  CA  TRP X%4d    %8.3f   0.000   0.000  1.00  0.00           C\n"
                % (i + 1, 900 + i, i * 3.8)
            )
    with open(lonely / "b.pdb", "w") as handle:
        for i in range(5):
            handle.write(
                "ATOM  %5d  CA  GLY Y%4d    %8.3f   0.000   0.000  1.00  0.00           C\n"
                % (i + 1, 500 + i, i * 3.8)
            )
    with pytest.raises(PreparationError, match="no structure could be prepared"):
        prepare_directory(str(lonely), str(tmp_path / "out"), verbose=False)


def test_the_report_records_every_structure(prepared, tmp_path):
    report, _ = prepared
    path = tmp_path / "preparation.csv"
    assert report.write_csv(path) == 6

    import csv
    rows = list(csv.DictReader(open(path)))
    assert {r["pdb_id"] for r in rows} == {
        "1A2P", "1BRN", "1BRS", "1BSA", "1BSE", "1RNB"
    }
    assert all(r["status"] == "ok" for r in rows)

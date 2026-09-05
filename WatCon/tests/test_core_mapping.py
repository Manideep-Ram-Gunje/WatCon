"""Live residue -> MSA mapping inside the network builders.

These tests need MDAnalysis and exercise the real code path in
``generate_static_networks``.  They are skipped automatically where MDAnalysis
is absent (see ``conftest.pytest_report_header``).

MODELLER is *not* needed: ``sequence_processing.generate_msa_alignment`` only
invokes it when the alignment file does not already exist, so the fixtures under
``tests/inputs/msa`` let the whole MSA path run in pure Python.

What is being defended
----------------------

The old code computed ``msa_indices[atm.resid - 1]`` -- positional-list
arithmetic driven by an absolute residue number.  It is correct only for dense,
1-origin, single-chain numbering, and wrong (silently) otherwise.  These tests
pin the replacement: lookup by ``(chain, resid, icode)``.
"""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("MDAnalysis", reason="core mapping tests require MDAnalysis")

import MDAnalysis as mda  # noqa: E402

from WatCon.generate_static_networks import (  # noqa: E402
    _build_residue_index,
    extract_objects,
)
from WatCon.residue_index import (  # noqa: E402
    ResidueIndex,
    ResidueIndexError,
    atom_identity,
    build_fasta_sequence,
    residues_from_pdb_file,
    residues_from_universe,
)
from WatCon.sequence_processing import (  # noqa: E402
    generate_msa_alignment,
    pdb_to_fastas,
)

from .conftest import bundle_member  # noqa: E402

INPUTS = Path(__file__).resolve().parent / "inputs"
PDB_1AKI = INPUTS / "1AKI.pdb"
MSA_DIR = INPUTS / "msa"
FASTA_1AKI = MSA_DIR / "1AKI.fa"
ALIGNMENT = MSA_DIR / "alignment.txt"


@pytest.fixture(scope="module")
def msa_indices():
    """MSA columns for 1AKI, from the committed gapped alignment fixture."""
    return generate_msa_alignment(str(ALIGNMENT), None, str(FASTA_1AKI))


@pytest.fixture(scope="module")
def bundle_pdb(tmp_path_factory):
    def _get(archive, member):
        path = tmp_path_factory.mktemp("struct") / member
        path.write_bytes(bundle_member(archive, member))
        return path

    return _get


# ===========================================================================
# The MDAnalysis adapter, previously never executed
# ===========================================================================

def test_universe_walk_matches_file_walk_on_1aki():
    universe = mda.Universe(str(PDB_1AKI))
    assert [r.key for r in residues_from_universe(universe)] == [
        r.key for r in residues_from_pdb_file(PDB_1AKI)
    ]


@pytest.mark.parametrize(
    "archive,member",
    [
        ("1788241750_ConSurf.tar.gz", "7O7W_A_ATOMS_section_With_ConSurf.pdb"),
        ("1788241897_ConSurf.tar.gz", "1BRS_A_ATOMS_section_With_ConSurf.pdb"),
        ("1788241883_ConSurf.tar.gz", "P00648_ATOMS_section_With_ConSurf.pdb"),
    ],
    ids=["7O7W_negative_resids", "1BRS_six_chains", "P00648_model"],
)
def test_universe_walk_matches_file_walk_on_hard_structures(
    bundle_pdb, archive, member
):
    """Negative residue numbers, six chains, and a predicted model."""
    path = bundle_pdb(archive, member)
    universe = mda.Universe(str(path))
    assert [r.key for r in residues_from_universe(universe)] == [
        r.key for r in residues_from_pdb_file(path)
    ]


def test_atom_identity_reads_chain_and_icode():
    universe = mda.Universe(str(PDB_1AKI))
    atom = universe.select_atoms("protein")[0]
    chain, icode = atom_identity(atom)
    assert chain == "A"
    assert icode is None


def test_atom_identity_agrees_with_the_residue_walk():
    """Every protein atom must resolve to a residue the walk produced."""
    universe = mda.Universe(str(PDB_1AKI))
    keys = {r.key for r in residues_from_pdb_file(PDB_1AKI)}
    for atom in universe.select_atoms("protein").atoms:
        chain, icode = atom_identity(atom)
        assert (chain, int(atom.resid), icode) in keys


# ===========================================================================
# _build_residue_index -- the guard, in the live path
# ===========================================================================

def test_build_residue_index_from_real_alignment(msa_indices):
    index = _build_residue_index(str(PDB_1AKI), msa_indices)
    assert len(index) == 129
    assert index.has_msa
    assert index.msa_column("A", 1) == msa_indices[0]
    assert index.msa_column("A", 129) == msa_indices[-1]


def test_build_residue_index_rejects_a_mismatched_alignment():
    """A FASTA describing a different structure must stop the run."""
    with pytest.raises(ResidueIndexError, match="length mismatch"):
        _build_residue_index(str(PDB_1AKI), list(range(1, 50)))


def test_build_residue_index_returns_none_without_msa():
    assert _build_residue_index(str(PDB_1AKI), None) is None


def test_build_residue_index_resid_fallback():
    """The dynamic pipeline's no-MSA fallback maps each residue to its number."""
    index = _build_residue_index(str(PDB_1AKI), None, fallback_to_resid=True)
    assert index is not None
    assert index.msa_column("A", 1) == 1
    assert index.msa_column("A", 129) == 129


# ===========================================================================
# End to end through extract_objects
# ===========================================================================

def _network_with_msa(msa_indices):
    index = _build_residue_index(str(PDB_1AKI), msa_indices)
    return index, extract_objects(
        str(PDB_1AKI),
        network_type="water-protein",
        custom_selection=None,
        active_region_reference=None,
        active_region_COM=False,
        active_region_radius=8.0,
        water_name=None,
        msa_indexing=index,
    )


def test_protein_atoms_carry_identity(msa_indices):
    _, network = _network_with_msa(msa_indices)
    assert network.protein_atoms
    for atom in network.protein_atoms:
        assert atom.chain == "A"
        assert atom.icode is None


def test_every_atom_msa_resid_matches_the_index(msa_indices):
    """The property the old arithmetic could not guarantee."""
    index, network = _network_with_msa(msa_indices)
    for atom in network.protein_atoms:
        assert atom.msa_resid == index.msa_column(atom.chain, atom.resid, atom.icode)


def test_msa_columns_are_not_the_trivial_identity(msa_indices):
    """The alignment is gapped, so resid and column genuinely differ."""
    index, network = _network_with_msa(msa_indices)
    assert msa_indices != list(range(1, 130))
    shifted = [
        atom for atom in network.protein_atoms if atom.msa_resid != atom.resid
    ]
    assert shifted, "fixture alignment should shift at least some residues"


def test_graph_nodes_agree_with_atoms(msa_indices):
    """Step 3d replaced 7 duplicate lookups with a read of OtherAtom.msa_resid."""
    _, network = _network_with_msa(msa_indices)
    by_index = {atom.index: atom for atom in network.protein_atoms}

    checked = 0
    for node, data in network.graph.nodes(data=True):
        if data.get("atom_category") != "PROTEIN":
            continue
        assert data["MSA"] == by_index[node].msa_resid
        checked += 1
    assert checked > 0


def test_no_msa_leaves_columns_unset():
    """msa_indexing=None must yield None, never a fabricated column."""
    network = extract_objects(
        str(PDB_1AKI),
        network_type="water-protein",
        custom_selection=None,
        active_region_reference=None,
        active_region_COM=False,
        active_region_radius=8.0,
        water_name=None,
        msa_indexing=None,
    )
    assert network.protein_atoms
    assert all(atom.msa_resid is None for atom in network.protein_atoms)
    for _, data in network.graph.nodes(data=True):
        if data.get("atom_category") == "PROTEIN":
            assert data["MSA"] is None


# ===========================================================================
# The failure mode this whole phase exists to remove
# ===========================================================================

def test_negative_residue_numbers_map_correctly(bundle_pdb):
    """7O7W's expression tag: HIS:-5:A is the FIRST residue, column 1.

    The old ``msa_indices[atm.resid - 1]`` evaluated ``msa_indices[-6]``, which
    Python resolves by wrapping to the end of the list -- no exception, just a
    column from the opposite end of the protein.
    """
    path = bundle_pdb(
        "1788241750_ConSurf.tar.gz", "7O7W_A_ATOMS_section_With_ConSurf.pdb"
    )
    residues = residues_from_pdb_file(path)
    columns = list(range(1, len(residues) + 1))
    index = ResidueIndex.build(residues, columns, source=str(path))

    assert index.msa_column("A", -5) == 1
    assert index.msa_column("A", 0) == 6

    legacy = index.legacy_indexing_report()
    assert legacy.negative_wraparound == 6
    assert legacy.correct == 0


def test_chain_is_part_of_the_key(bundle_pdb):
    """1BRS has six chains; residue 10 exists in each and they are different."""
    path = bundle_pdb(
        "1788241897_ConSurf.tar.gz", "1BRS_A_ATOMS_section_With_ConSurf.pdb"
    )
    residues = residues_from_pdb_file(path)
    index = ResidueIndex.build(residues, list(range(1, len(residues) + 1)))

    assert len(index.chains()) == 6
    columns = {
        chain: index.msa_column(chain, 10)
        for chain in index.chains()
        if index.ordinal(chain, 10) is not None
    }
    assert len(set(columns.values())) == len(columns), (
        "each chain's residue 10 must map to a distinct column"
    )


# ===========================================================================
# FASTA generation now shares the residue walk
# ===========================================================================

def test_pdb_to_fastas_matches_the_shared_walk(tmp_path):
    pdb_to_fastas(str(PDB_1AKI), str(tmp_path), name="1AKI")
    written = (tmp_path / "1AKI.fa").read_text()
    expected = ">1AKI\n" + build_fasta_sequence(residues_from_pdb_file(PDB_1AKI))
    assert written == expected


def test_pdb_to_fastas_reproduces_the_committed_fixture(tmp_path):
    """The fixture the MSA tests rely on must stay in step with the code."""
    pdb_to_fastas(str(PDB_1AKI), str(tmp_path), name="1AKI")
    assert (tmp_path / "1AKI.fa").read_text() == FASTA_1AKI.read_text()


def test_pdb_to_fastas_can_restrict_to_one_chain(tmp_path, bundle_pdb):
    path = bundle_pdb(
        "1788241897_ConSurf.tar.gz", "1BRS_A_ATOMS_section_With_ConSurf.pdb"
    )
    pdb_to_fastas(str(path), str(tmp_path), name="all")
    pdb_to_fastas(str(path), str(tmp_path), name="chainA", chain="A")

    everything = (tmp_path / "all.fa").read_text().split("\n")[1]
    chain_a = (tmp_path / "chainA.fa").read_text().split("\n")[1]

    assert len(chain_a) < len(everything)
    assert len(chain_a) == len(residues_from_pdb_file(path, chain="A"))


def test_fasta_and_index_stay_in_step(tmp_path, msa_indices):
    """The guard only works if FASTA and index come from the same walk."""
    pdb_to_fastas(str(PDB_1AKI), str(tmp_path), name="1AKI")
    sequence = (tmp_path / "1AKI.fa").read_text().split("\n")[1]
    index = _build_residue_index(str(PDB_1AKI), msa_indices)
    assert len(sequence) == len(index) == len(msa_indices)

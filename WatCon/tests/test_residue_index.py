"""Tests for the residue-identity and MSA-mapping layer.

Real structures used here:

* ``WatCon/tests/inputs/1AKI.pdb`` -- lysozyme, chain A, residues 1..129, dense.
  The one case where the old ``msa_indices[resid - 1]`` arithmetic happens to be
  correct, which is why the bug survived so long.
* the ConSurf bundles under ``data/consurf/exploratory`` -- 7O7W carries an
  expression tag numbered -5..0, P00648 is UniProt-numbered from 48, and 1BRS
  has two chains.

The MDAnalysis adapter :func:`residues_from_universe` is NOT covered here; it
cannot run without MDAnalysis.  See ``docs/CONSURF_INTEGRATION.md``.
"""

from __future__ import annotations

import pytest

from WatCon.residue_index import (
    AMINO_ACID_3TO1,
    ResidueIndex,
    ResidueIndexError,
    StructureResidue,
    build_fasta_sequence,
    check_sequence_consistency,
    legacy_indexing_report,
    residues_from_pdb_file,
)

from .conftest import EXPLORATORY, REPO_ROOT, bundle_member

PDB_1AKI = None  # resolved by the fixture below


@pytest.fixture(scope="module")
def aki_path():
    from pathlib import Path

    return Path(__file__).resolve().parent / "inputs" / "1AKI.pdb"


@pytest.fixture(scope="module")
def aki(aki_path):
    return residues_from_pdb_file(aki_path)


def _bundle_pdb(tmp_path_factory, archive, member):
    path = tmp_path_factory.mktemp("struct") / member
    path.write_bytes(bundle_member(archive, member))
    return path


@pytest.fixture(scope="module")
def o7w_path(tmp_path_factory):
    return _bundle_pdb(
        tmp_path_factory,
        "1788241750_ConSurf.tar.gz",
        "7O7W_A_ATOMS_section_With_ConSurf.pdb",
    )


@pytest.fixture(scope="module")
def brs_path(tmp_path_factory):
    return _bundle_pdb(
        tmp_path_factory,
        "1788241897_ConSurf.tar.gz",
        "1BRS_A_ATOMS_section_With_ConSurf.pdb",
    )


@pytest.fixture(scope="module")
def p00648_path(tmp_path_factory):
    return _bundle_pdb(
        tmp_path_factory,
        "1788241883_ConSurf.tar.gz",
        "P00648_ATOMS_section_With_ConSurf.pdb",
    )


# ===========================================================================
# Reading residues
# ===========================================================================

def test_1aki_residue_walk(aki):
    assert len(aki) == 129
    assert aki[0].resid == 1 and aki[0].resname == "LYS"
    assert aki[-1].resid == 129 and aki[-1].resname == "LEU"
    assert {r.chain for r in aki} == {"A"}
    assert all(r.icode is None for r in aki)


def test_residue_order_is_file_order(aki):
    assert [r.resid for r in aki] == list(range(1, 130))


def test_negative_and_zero_residue_numbers_are_read(o7w_path):
    residues = residues_from_pdb_file(o7w_path)
    tag = [r for r in residues if r.resid <= 0]
    assert [r.resid for r in tag] == [-5, -4, -3, -2, -1, 0]
    assert [r.resname for r in tag] == ["HIS", "HIS", "HIS", "THR", "ASP", "PRO"]


def test_multiple_chains_are_kept_separate(brs_path):
    everything = residues_from_pdb_file(brs_path)
    chains = {r.chain for r in everything}
    assert len(chains) > 1, "1BRS is expected to be multi-chain"

    only_a = residues_from_pdb_file(brs_path, chain="A")
    assert {r.chain for r in only_a} == {"A"}
    assert len(only_a) < len(everything)


def test_chain_filter_changes_the_ordinals(brs_path):
    """The bug class this module removes: identity must include the chain."""
    everything = residues_from_pdb_file(brs_path)
    only_a = residues_from_pdb_file(brs_path, chain="A")

    all_index = ResidueIndex.build(everything, list(range(1, len(everything) + 1)))
    a_index = ResidueIndex.build(only_a, list(range(1, len(only_a) + 1)))

    # Same residue number, different chain -> different residues.
    assert ("A", 10, None) in all_index
    assert ("B", 10, None) in all_index
    assert ("B", 10, None) not in a_index


def test_waters_and_ligands_are_excluded(o7w_path):
    residues = residues_from_pdb_file(o7w_path)
    assert all(r.resname != "HOH" for r in residues)
    # The chromophore PIA is a HETATM and not in the residue table.
    assert all(r.resname != "PIA" for r in residues)


def test_hetatm_opt_in(o7w_path):
    without = residues_from_pdb_file(o7w_path)
    with_het = residues_from_pdb_file(o7w_path, include_hetatm=True)
    assert len(with_het) >= len(without)


def test_custom_residue_table_adds_a_name():
    """A modified residue name can be brought into the table."""
    from WatCon.residue_index import _resolve_table

    assert "PIA" not in _resolve_table(None)
    assert _resolve_table({"PIA": "X"})["PIA"] == "X"


def test_chromophore_cannot_be_recovered_by_the_ca_walk(o7w_path):
    """Documented limitation: the walk keys on an atom named exactly 'CA'.

    7O7W's chromophore PIA:68:A is a fused tripeptide whose alpha carbons are
    named CA1, CA2 and CA3, so neither include_hetatm nor a custom residue name
    can bring it in.  ConSurf omits it too, so the two agree -- but anything
    that needs such a residue must locate it by another route.
    """
    base = residues_from_pdb_file(o7w_path, include_hetatm=True)
    extended = residues_from_pdb_file(
        o7w_path, include_hetatm=True, custom_residues={"PIA": "X"}
    )
    assert len(extended) == len(base)
    assert all(r.resname != "PIA" for r in extended)


def test_sequence_uses_the_residue_table(aki):
    sequence = build_fasta_sequence(aki)
    assert len(sequence) == 129
    assert set(sequence) <= set(AMINO_ACID_3TO1.values())


# ===========================================================================
# Parity with the historical FASTA walk
# ===========================================================================

def _legacy_fasta(path):
    """Reproduce sequence_processing.pdb_to_fastas exactly, warts and all.

    Cannot be imported directly: sequence_processing does
    ``from modeller import *`` at module scope.
    """
    sequence = []
    for line in open(path, encoding="utf-8", errors="replace"):
        if ("ATOM" in line) and ("CA" in line) and ("ANISOU" not in line):
            resname = line[17:20]
            if resname in AMINO_ACID_3TO1:
                sequence.append(AMINO_ACID_3TO1[resname])
    return "".join(sequence)


def test_fasta_parity_on_well_formed_single_chain(aki_path, aki):
    """On a clean single-chain PDB the new walk reproduces the old sequence."""
    assert build_fasta_sequence(aki) == _legacy_fasta(aki_path)


def test_new_walk_differs_where_the_old_one_was_wrong(brs_path):
    """The legacy predicate is a substring test and merges every chain.

    ``'ATOM' in line`` also matches HETATM, and there is no chain filter, so on a
    multi-chain structure the legacy sequence is not the chain-A sequence.
    """
    chain_a = build_fasta_sequence(residues_from_pdb_file(brs_path, chain="A"))
    legacy = _legacy_fasta(brs_path)
    assert legacy != chain_a
    assert len(legacy) > len(chain_a)


# ===========================================================================
# Consistency checking
# ===========================================================================

def test_consistency_report_matches(aki, aki_path):
    report = check_sequence_consistency(aki, _legacy_fasta(aki_path))
    assert report.matches
    assert report.structure_length == report.fasta_length == 129
    assert "consistent" in report.describe()


def test_consistency_report_detects_length_drift(aki):
    report = check_sequence_consistency(aki, build_fasta_sequence(aki)[:-5])
    assert not report.matches
    assert report.structure_length == 129
    assert report.fasta_length == 124
    assert "length mismatch" in report.describe()


def test_consistency_report_detects_substitution(aki):
    sequence = list(build_fasta_sequence(aki))
    sequence[40] = "W" if sequence[40] != "W" else "G"
    report = check_sequence_consistency(aki, "".join(sequence))
    assert not report.matches
    assert report.first_mismatch == 40
    assert "differ from ordinal 40" in report.describe()


def test_consistency_tolerates_whitespace_and_case(aki):
    wrapped = "\n".join(
        build_fasta_sequence(aki)[i : i + 60]
        for i in range(0, 129, 60)
    ).lower()
    assert check_sequence_consistency(aki, wrapped).matches


# ===========================================================================
# ResidueIndex.build -- the guard
# ===========================================================================

def test_length_mismatch_is_refused(aki):
    with pytest.raises(ResidueIndexError, match="length mismatch"):
        ResidueIndex.build(aki, list(range(1, 50)))


def test_length_mismatch_message_names_the_source(aki):
    with pytest.raises(ResidueIndexError, match="my_structure"):
        ResidueIndex.build(aki, [1, 2, 3], source="my_structure.pdb")


def test_duplicate_identity_is_refused():
    duplicate = [
        StructureResidue("A", 1, None, "GLY", "G"),
        StructureResidue("A", 1, None, "ALA", "A"),
    ]
    with pytest.raises(ResidueIndexError, match="Duplicate residue identities"):
        ResidueIndex.build(duplicate, [1, 2])


def test_index_without_msa_is_allowed(aki):
    index = ResidueIndex.build(aki)
    assert not index.has_msa
    assert index.ordinal("A", 5) == 4
    assert index.msa_column("A", 5) is None
    with pytest.raises(ResidueIndexError, match="no MSA indices"):
        index.require_msa_column("A", 5)


# ===========================================================================
# Lookup
# ===========================================================================

def test_lookup_by_identity(aki):
    index = ResidueIndex.build(aki, [10 * (i + 1) for i in range(len(aki))])
    assert index.ordinal("A", 1) == 0
    assert index.msa_column("A", 1) == 10
    assert index.ordinal("A", 129) == 128
    assert index.msa_column("A", 129) == 1290


def test_lookup_of_unknown_residue_returns_none_not_a_guess(aki):
    index = ResidueIndex.build(aki, list(range(1, len(aki) + 1)))
    assert index.ordinal("A", 9999) is None
    assert index.msa_column("A", 9999) is None
    assert index.msa_column("Z", 1) is None
    with pytest.raises(ResidueIndexError, match="No MSA column"):
        index.require_msa_column("A", 9999)


def test_negative_residue_lookup_is_correct(o7w_path):
    """The case the old arithmetic silently got wrong."""
    residues = residues_from_pdb_file(o7w_path)
    index = ResidueIndex.build(residues, list(range(1, len(residues) + 1)))

    assert index.ordinal("A", -5) == 0
    assert index.msa_column("A", -5) == 1
    assert index.ordinal("A", 0) == 5
    assert index.msa_column("A", 0) == 6


def test_offset_numbering_lookup_is_correct(p00648_path):
    """P00648 is UniProt-numbered; the first modelled residue is 1, not 48."""
    residues = residues_from_pdb_file(p00648_path)
    index = ResidueIndex.build(residues, list(range(1, len(residues) + 1)))
    first = residues[0]
    assert index.ordinal(first.chain, first.resid) == 0
    assert index.msa_column(first.chain, first.resid) == 1


def test_insertion_code_is_part_of_identity():
    residues = [
        StructureResidue("A", 100, None, "HIS", "H"),
        StructureResidue("A", 100, "A", "HIS", "H"),
    ]
    index = ResidueIndex.build(residues, [7, 8])
    assert index.msa_column("A", 100) == 7
    assert index.msa_column("A", 100, "A") == 8
    assert index.ordinal("A", 100, "B") is None


def test_index_introspection(aki):
    index = ResidueIndex.build(aki, list(range(1, len(aki) + 1)))
    assert len(index) == 129
    assert index.chains() == ["A"]
    assert index.residue_at(0).resid == 1
    assert len(index.keys()) == 129
    assert index.sequence() == build_fasta_sequence(aki)
    assert ("A", 1, None) in index


def test_from_pdb_file_convenience(aki_path):
    index = ResidueIndex.from_pdb_file(aki_path, list(range(1, 130)))
    assert len(index) == 129
    assert index.msa_column("A", 1) == 1
    assert "1AKI" in (index.source or "")


# ===========================================================================
# Legacy-arithmetic diagnostic
# ===========================================================================

def test_legacy_indexing_is_safe_on_1aki(aki):
    """1AKI is dense and 1-origin -- the one shape the old code handled."""
    report = legacy_indexing_report(aki, list(range(1, len(aki) + 1)))
    assert report.is_safe
    assert report.correct == 129
    assert report.wrong == 0
    assert report.negative_wraparound == 0


def test_legacy_indexing_breaks_on_7o7w(o7w_path):
    """Quantifies the bug: every residue wrong, six silently wrapped."""
    residues = residues_from_pdb_file(o7w_path)
    report = legacy_indexing_report(residues, list(range(1, len(residues) + 1)))

    assert not report.is_safe
    assert report.correct == 0
    assert report.negative_wraparound == 6

    # The wrap-around returns a column from the far end of the sequence.
    keys = {key: legacy for key, _, legacy in report.examples}
    assert keys[("A", -5, None)] > len(residues) - 10


def test_legacy_indexing_is_safe_on_a_dense_model(p00648_path):
    """P00648's PDB is an AlphaFold model numbered 1..157, dense and 1-origin.

    Its notorious +47 offset is between ConSurf POS and the PDB residue number,
    not inside the structure -- so the legacy arithmetic is safe *here* while
    still being wrong for the ConSurf join.  The two coordinate systems must not
    be conflated; the POS side is the parser's problem and is tested there.
    """
    residues = residues_from_pdb_file(p00648_path)
    assert [r.resid for r in residues] == list(range(1, len(residues) + 1))

    report = legacy_indexing_report(residues, list(range(1, len(residues) + 1)))
    assert report.is_safe


def test_legacy_indexing_breaks_when_numbering_does_not_start_at_one():
    """A structure whose first residue is not 1 shifts every lookup."""
    residues = [
        StructureResidue("A", n, None, "GLY", "G") for n in range(3, 13)
    ]
    report = legacy_indexing_report(residues, list(range(1, 11)))
    assert not report.is_safe
    assert report.wrong + report.out_of_range == len(residues)


def test_legacy_report_via_index(aki):
    index = ResidueIndex.build(aki, list(range(1, len(aki) + 1)))
    assert index.legacy_indexing_report().is_safe
    assert "correct" in index.legacy_indexing_report().describe()


# ===========================================================================
# Module import hygiene
# ===========================================================================

def test_residue_index_does_not_require_mdanalysis():
    """Importing this module must never pull in MDAnalysis or modeller.

    Checked in a subprocess: once any other test module has imported
    MDAnalysis, an in-process ``sys.modules`` check only reflects collection
    order, not what this module actually imports.
    """
    import subprocess
    import sys

    probe = (
        "import sys; import WatCon.residue_index; "
        "assert 'MDAnalysis' not in sys.modules, 'residue_index pulled in MDAnalysis'; "
        "assert 'modeller' not in sys.modules, 'residue_index pulled in modeller'"
    )
    result = subprocess.run(
        [sys.executable, "-c", probe],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr


def test_sequence_processing_does_not_require_modeller():
    """MODELLER is licensed and often absent; it must be imported lazily.

    A module-level ``from modeller import *`` used to make the whole network
    module chain unimportable without it.
    """
    import subprocess
    import sys

    probe = (
        "import sys; import WatCon.generate_static_networks; "
        "assert 'modeller' not in sys.modules, 'modeller imported at module scope'"
    )
    result = subprocess.run(
        [sys.executable, "-c", probe],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr


# ===========================================================================
# Alternate conformers
# ===========================================================================

def test_residues_with_only_non_primary_altlocs_are_kept(o7w_path):
    """Regression: 11 residues of 7O7W have only altloc B and C conformers.

    An earlier version whitelisted altloc {"", "A", "1"} and dropped them
    silently, which desynchronised the walk from both the FASTA and ConSurf.
    """
    residues = residues_from_pdb_file(o7w_path)
    found = {r.resid for r in residues}
    for resid in (65, 69, 146, 147, 148, 149, 150, 151, 152, 204, 223):
        assert resid in found, f"residue {resid} (altloc B/C only) was dropped"


def test_alternate_conformers_yield_one_residue_each(o7w_path):
    """Multiple conformers of one residue must collapse to a single entry."""
    residues = residues_from_pdb_file(o7w_path)
    keys = [r.key for r in residues]
    assert len(keys) == len(set(keys))


# ===========================================================================
# Cross-layer agreement with ConSurf
# ===========================================================================
#
# This is the check that matters for the eventual join: the identity the
# structure walk produces must be the identity ConSurf reports.  It is what
# caught the altloc bug above.

@pytest.mark.parametrize(
    "archive,base,expect_extra_in_walk",
    [
        ("1788241750_ConSurf.tar.gz", "7O7W_A", 0),
        ("1788241897_ConSurf.tar.gz", "1BRS_A", 0),
        # P00648 is an AlphaFold model numbered 1..157; ConSurf covered 48..153,
        # so 51 modelled residues legitimately have no ConSurf record.
        ("1788241883_ConSurf.tar.gz", "P00648", 51),
    ],
)
def test_walk_agrees_with_consurf_identities(
    tmp_path_factory, archive, base, expect_extra_in_walk
):
    import io

    from WatCon.consurf import parse_consurf

    pdb_path = _bundle_pdb(
        tmp_path_factory, archive, f"{base}_ATOMS_section_With_ConSurf.pdb"
    )
    grades = bundle_member(archive, f"{base}_consurf_grades.txt").decode()
    result = parse_consurf(io.StringIO(grades), strict=True)

    walk = {r.key for r in residues_from_pdb_file(pdb_path, chain="A")}
    consurf = set(result.by_pdb_residue())

    # Every residue ConSurf mapped must exist in the structure walk.
    assert consurf - walk == set()
    assert len(walk - consurf) == expect_extra_in_walk


def test_consurf_lookup_through_the_index(tmp_path_factory):
    """End-to-end: structure identity -> index ordinal -> ConSurf record."""
    import io

    from WatCon.consurf import LookupStatus, parse_consurf

    pdb_path = _bundle_pdb(
        tmp_path_factory,
        "1788241750_ConSurf.tar.gz",
        "7O7W_A_ATOMS_section_With_ConSurf.pdb",
    )
    grades = bundle_member(
        "1788241750_ConSurf.tar.gz", "7O7W_A_consurf_grades.txt"
    ).decode()
    result = parse_consurf(io.StringIO(grades), strict=True)

    residues = residues_from_pdb_file(pdb_path, chain="A")
    index = ResidueIndex.build(residues, list(range(1, len(residues) + 1)))

    # The expression tag: correct ordinal, and a real ConSurf record.
    assert index.ordinal("A", -5) == 0
    lookup = result.lookup_structural("A", -5)
    assert lookup.status is LookupStatus.SCORED
    assert lookup.record.pdb_residue.name3 == "HIS"

    # Every residue in the walk resolves to a ConSurf grade.
    for residue in residues:
        hit = result.lookup_structural(*residue.key)
        assert hit.status is LookupStatus.SCORED
        assert hit.record.pdb_residue.name3 == residue.resname


# ===========================================================================
# Multi-chain reality
# ===========================================================================

def test_1brs_has_six_chains(brs_path):
    """1BRS is a barnase-barstar complex with six chains, not one."""
    residues = residues_from_pdb_file(brs_path)
    index = ResidueIndex.build(residues, list(range(1, len(residues) + 1)))
    assert index.chains() == ["A", "B", "C", "D", "E", "F"]
    assert len(residues) == 588


def test_legacy_indexing_is_catastrophic_on_multichain(brs_path):
    """With six chains sharing residue numbers, resid-1 indexing is meaningless."""
    residues = residues_from_pdb_file(brs_path)
    report = legacy_indexing_report(residues, list(range(1, len(residues) + 1)))
    assert report.correct == 0
    assert not report.is_safe

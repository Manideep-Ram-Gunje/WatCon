"""Independent validation against ConSurf's own annotated PDB output.

This is REAL-output validation, and the strongest evidence we have that the
residue-identity mapping is right.  ConSurf writes the conservation grade into
the B-factor column of ``*_ATOMS_section_With_ConSurf.pdb`` and leaves it blank
for residues it did not score.  That artifact is produced by ConSurf itself and
is entirely independent of the grades table, so agreement between the two is a
genuine cross-check rather than a restatement.

Deliberately NOT used as an oracle: ``msa_aa_variety_percentage.csv``.  Its
``pos`` and ``ConSurf grade`` columns are indexed on the grades POS, but its
amino-acid composition columns are shifted by one row from the point where
ConSurf omitted a residue.  See ``docs/CONSURF_INTEGRATION.md``.
"""

from __future__ import annotations

import io

import pytest

from WatCon.consurf import (
    LookupStatus,
    check_against_annotated_pdb,
    parse_consurf,
    read_annotated_pdb,
)

from .conftest import BUNDLES, bundle_member


@pytest.fixture(scope="module", params=BUNDLES, ids=[b[1] for b in BUNDLES])
def bundle(request, tmp_path_factory):
    """Parse a real ConSurf bundle and lay its annotated PDB out on disk."""
    archive, base = request.param

    grades = bundle_member(archive, f"{base}_consurf_grades.txt").decode()
    result = parse_consurf(io.StringIO(grades), strict=True)

    pdb_bytes = bundle_member(archive, f"{base}_ATOMS_section_With_ConSurf.pdb")
    pdb_path = tmp_path_factory.mktemp("consurf") / f"{base}.pdb"
    pdb_path.write_bytes(pdb_bytes)

    return base, result, pdb_path


def test_grades_match_annotated_pdb_bfactor(bundle):
    """Every mapped record's grade must equal the B-factor ConSurf wrote."""
    base, result, pdb_path = bundle
    report = check_against_annotated_pdb(result, pdb_path)

    assert report.mismatches == [], f"{base}: {report.mismatches[:5]}"
    assert report.missing_from_pdb == [], f"{base}: {report.missing_from_pdb[:5]}"
    assert report.checked == len(result.mapped_records())
    assert report.agreements == report.checked
    assert report.ok


def test_every_mapped_record_is_present_in_the_structure(bundle):
    base, result, pdb_path = bundle
    pdb_residues = read_annotated_pdb(pdb_path)
    for record in result.mapped_records():
        assert record.pdb_residue.key in pdb_residues


def test_residue_names_agree_with_the_structure(bundle):
    """The three-letter code in the grades table must match the PDB."""
    base, result, pdb_path = bundle

    names = {}
    for line in pdb_path.read_text(errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")) and len(line) > 26:
            try:
                resseq = int(line[22:26].strip())
            except ValueError:
                continue
            key = (line[21], resseq, line[26].strip() or None)
            names.setdefault(key, line[17:20].strip())

    for record in result.mapped_records():
        assert names[record.pdb_residue.key] == record.pdb_residue.name3


def test_unscored_structural_residues_are_reported(bundle):
    """Residues in the structure that ConSurf never scored must be visible."""
    base, result, pdb_path = bundle
    report = check_against_annotated_pdb(result, pdb_path)

    # Every bundle has some: waters, and either a modified residue or the part
    # of a model that fell outside the analysed region.
    assert report.unscored_residues

    for key in report.unscored_residues:
        chain, seq, icode = key
        assert result.lookup_structural(chain, seq, icode).status in (
            LookupStatus.NO_SCORE_IN_FILE,
            LookupStatus.NOT_IN_FILE,
        )


def test_7o7w_chromophore_is_present_but_unscored(tmp_path):
    """The rsEGFP2 chromophore PIA:68:A exists in the structure with no grade.

    This is the concrete case that motivates the three-valued lookup: POS is not
    a dense sequence index, and a structurally important residue can be absent
    from the grades table entirely.
    """
    grades = bundle_member(
        "1788241750_ConSurf.tar.gz", "7O7W_A_consurf_grades.txt"
    ).decode()
    result = parse_consurf(io.StringIO(grades), strict=True)

    pdb_path = tmp_path / "7O7W.pdb"
    pdb_path.write_bytes(
        bundle_member(
            "1788241750_ConSurf.tar.gz", "7O7W_A_ATOMS_section_With_ConSurf.pdb"
        )
    )

    pdb_residues = read_annotated_pdb(pdb_path)

    # Present in the structure, blank B-factor, no ConSurf record.
    assert ("A", 68, None) in pdb_residues
    assert pdb_residues[("A", 68, None)] is None
    assert ("A", 68, None) not in result.by_pdb_residue()

    lookup = result.lookup_structural("A", 68)
    assert lookup.status is LookupStatus.NO_SCORE_IN_FILE
    assert lookup.record is None

    # And it is exactly what the reported numbering gap points at.
    assert result.position_gaps == [(65, 69)]


def test_1brs_second_chain_is_outside_the_run(tmp_path):
    """1BRS was run on chain A; chain B exists in the structure but is unscored."""
    grades = bundle_member(
        "1788241897_ConSurf.tar.gz", "1BRS_A_consurf_grades.txt"
    ).decode()
    result = parse_consurf(io.StringIO(grades), strict=True)

    assert result.chains() == ["A"]
    assert result.lookup_structural("B", 10).status is LookupStatus.NOT_IN_FILE


def test_p00648_model_region_outside_the_run(tmp_path):
    """The AlphaFold model is numbered 1-157; ConSurf scored only 48-153."""
    grades = bundle_member(
        "1788241883_ConSurf.tar.gz", "P00648_consurf_grades.txt"
    ).decode()
    result = parse_consurf(io.StringIO(grades), strict=True)

    assert result.covered_range("A") == (48, 153)
    assert result.lookup_structural("A", 1).status is LookupStatus.NOT_IN_FILE
    assert result.lookup_structural("A", 48).status is LookupStatus.SCORED


def test_blank_bfactor_is_not_counted_as_a_mismatch(tmp_path):
    """A bundle with no grades written into the PDB must not fail the check."""
    grades = bundle_member(
        "1788241897_ConSurf.tar.gz", "1BRS_A_consurf_grades.txt"
    ).decode()
    result = parse_consurf(io.StringIO(grades), strict=True)

    stripped = tmp_path / "blank.pdb"
    source = bundle_member(
        "1788241897_ConSurf.tar.gz", "1BRS_A_ATOMS_section_With_ConSurf.pdb"
    ).decode(errors="replace")
    stripped.write_text(
        "\n".join(
            line[:60] + " " * 6 + line[66:] if len(line) >= 66 else line
            for line in source.splitlines()
        )
    )

    report = check_against_annotated_pdb(result, stripped)
    assert report.checked == 0
    assert report.mismatches == []
    assert len(report.blank_in_pdb) == len(result.mapped_records())


def test_mismatch_is_actually_detected(tmp_path):
    """Guard against the cross-check silently passing on everything."""
    grades = bundle_member(
        "1788241897_ConSurf.tar.gz", "1BRS_A_consurf_grades.txt"
    ).decode()
    result = parse_consurf(io.StringIO(grades), strict=True)

    corrupted = tmp_path / "corrupt.pdb"
    source = bundle_member(
        "1788241897_ConSurf.tar.gz", "1BRS_A_ATOMS_section_With_ConSurf.pdb"
    ).decode(errors="replace")

    lines = []
    for line in source.splitlines():
        if line.startswith("ATOM") and line[22:26].strip() == "3" and len(line) >= 66:
            line = line[:60] + "     1" + line[66:]
        lines.append(line)
    corrupted.write_text("\n".join(lines))

    report = check_against_annotated_pdb(result, corrupted)
    assert report.mismatches, "corrupting a grade must be detected"
    assert not report.ok

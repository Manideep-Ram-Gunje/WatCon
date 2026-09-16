"""``watcon family`` -- the command, on the repository's own real data.

The members file is tab-separated rather than ``name:dir:file`` because Windows
paths contain colons: a format that breaks on ``C:/structures`` is no format at
all. These tests cover that file's failure modes as much as the happy path,
since a user's first attempt is where a bad format shows up.

Members here are the **full-length** CA extracts of the ten structures. That is
deliberate, and it is what the command requires: placing a structure on an
alignment row needs the whole chain. A structure trimmed to a pocket is several
disconnected segments that no pairwise alignment can place, and the command says
so rather than guessing -- ``test_a_trimmed_structure_is_refused_with_a_reason``.

The water-site half of the command needs structures that still contain waters,
which the CA extracts do not; that path is covered on real data in
``test_family_sites.py``. Here it is checked that missing waters are reported
clearly rather than crashing.
"""

from __future__ import annotations

import os
import shutil

import pytest

pytest.importorskip("Bio", reason="family mapping needs Biopython")
pytest.importorskip("MDAnalysis", reason="building networks needs MDAnalysis")

from WatCon.cli import build_parser, main

PACKAGE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FAMILY_DIR = os.path.join(PACKAGE, "data", "examples", "ptp_family")
FIXTURES = os.path.join(PACKAGE, "data", "consurf", "fixtures")
ALIGNMENT = os.path.join(FAMILY_DIR, "ptp_family_alignment.pir")

#: protein -> (ConSurf run, reference, both structures)
MEMBERS = {
    "PTPN1": ("1AAX", "2F71", ("2F71", "8U1E")),
    "PTPN6": ("4GRZ", "4GRZ", ("4GRZ", "4HJP")),
    "PTPN7": ("1ZC0", "1ZC0", ("1ZC0", "3O4U")),
    "PTPN12": ("5HDE", "5HDE", ("5HDE", "5J8R")),
    "PTPN22": ("3BRH", "3BRH", ("3BRH", "3OLR")),
}

pytestmark = pytest.mark.skipif(not os.path.isfile(ALIGNMENT),
                                reason="PTP family fixtures absent")


def _members_file(tmp_path, source=lambda pdb_id: "%s_A_ca.pdb" % pdb_id):
    """Lay the fixtures out as a user would, and write the members file."""
    lines = ["# protein\tstructures\tconsurf\treference"]
    for protein, (run, reference, structures) in MEMBERS.items():
        directory = tmp_path / protein
        directory.mkdir()
        for pdb_id in structures:
            # named by PDB id, which is how the alignment row is found
            shutil.copyfile(os.path.join(FAMILY_DIR, source(pdb_id)),
                            str(directory / ("%s.pdb" % pdb_id)))
        lines.append("\t".join([protein, str(directory),
                                os.path.join(FIXTURES, "%s_A.grades.txt" % run), reference]))
    path = tmp_path / "members.tsv"
    path.write_text("\n".join(lines) + "\n")
    return str(path)


# ===========================================================================
# The interface
# ===========================================================================

def test_family_is_an_advertised_subcommand():
    parser = build_parser()
    commands = set()
    for action in parser._actions:
        if getattr(action, "choices", None):
            commands.update(action.choices)
    assert "family" in commands


def test_members_and_alignment_are_required():
    with pytest.raises(SystemExit):
        build_parser().parse_args(["family"])


# ===========================================================================
# The members file, where a first attempt goes wrong
# ===========================================================================

def test_a_missing_members_file_is_reported(tmp_path, capsys):
    code = main(["family", "--members", str(tmp_path / "nope.tsv"),
                 "--alignment", ALIGNMENT])
    assert code == 1
    assert "nope.tsv" in capsys.readouterr().err


def test_a_short_line_names_the_line_and_what_was_expected(tmp_path, capsys):
    path = tmp_path / "members.tsv"
    path.write_text("PTPN1\tonly-two-fields\n")
    code = main(["family", "--members", str(path), "--alignment", ALIGNMENT])
    assert code == 1
    error = capsys.readouterr().err
    assert "line 1" in error and "three tab-separated fields" in error


def test_a_reference_that_is_not_there_is_refused(tmp_path, capsys):
    members = _members_file(tmp_path)
    # read fully first: open(..., "w") truncates before the read would happen
    text = open(members).read().replace("\t2F71\n", "\t9XYZ\n")
    open(members, "w").write(text)
    code = main(["family", "--members", members, "--alignment", ALIGNMENT])
    assert code == 1
    assert "9XYZ" in capsys.readouterr().err


def test_a_bad_state_label_is_refused(tmp_path, capsys):
    code = main(["family", "--members", _members_file(tmp_path), "--alignment", ALIGNMENT,
                 "--state", "3OLR-open"])
    assert code == 1
    assert "PDBID=LABEL" in capsys.readouterr().err


def test_a_trimmed_structure_is_refused_with_a_reason(tmp_path, capsys):
    """A pocket is disconnected segments; no alignment can place it on a row."""
    trimmed = {"2F71": "2F71_A_site.pdb", "8U1E": "8U1E_A_ca.pdb"}
    members = _members_file(tmp_path, source=lambda p: trimmed.get(p, "%s_A_ca.pdb" % p))
    code = main(["family", "--members", members, "--alignment", ALIGNMENT, "--no-sites"])
    assert code == 1
    error = capsys.readouterr().err
    assert "supply the whole chain" in error


# ===========================================================================
# End to end
# ===========================================================================

def test_conservation_run_reports_the_audit(tmp_path, capsys):
    code = main(["family", "--members", _members_file(tmp_path), "--alignment", ALIGNMENT,
                 "--out-dir", str(tmp_path / "out"), "--no-sites"])
    assert code == 0
    out = capsys.readouterr().out
    for protein in MEMBERS:
        assert protein in out
    # the engineered residues are named, not hidden
    assert "differs from the ConSurf query at A215" in out     # 1AAX query is C215S
    assert "differs from the ConSurf query at A195" in out     # 3BRH query is D195A
    # the alignment's own defect is reported
    assert "124, 125, 241" in out
    assert "248 covered by all 5 proteins" in out
    assert "64 of those unanimously conserved" in out


def test_structures_without_waters_are_reported_not_crashed(tmp_path, capsys):
    """The CA extracts have no waters. Say so; do not fail inside the clustering."""
    code = main(["family", "--members", _members_file(tmp_path), "--alignment", ALIGNMENT,
                 "--out-dir", str(tmp_path / "out")])
    assert code == 1
    error = capsys.readouterr().err
    assert "water" in error.lower()

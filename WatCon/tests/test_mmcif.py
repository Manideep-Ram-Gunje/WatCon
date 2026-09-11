"""mmCIF input, and fetching structures by id.

WatCon's preparation step parses ATOM/HETATM text directly, so it only ever
understood PDB format::

    ValueError: 'CIF' isn't a valid topology format

That stopped being an acceptable limitation. RCSB no longer issues PDB files for
large or recent entries -- fetching the PTP1B benchmark set, **31 of 287
selected entries had no PDB file at all**. 7GSA returns 404 for ``.pdb`` and
560 kB for ``.cif``. MDAnalysis 2.10 does not read mmCIF either, so the
conversion has to happen before any reader sees the file.

Network tests are opt-in (``WATCON_NETWORK_TESTS=1``): CI must not fail because
RCSB is briefly down. Everything that can be tested offline, is.
"""

from __future__ import annotations

import gzip
import os
import shutil

import pytest

from WatCon.structure_io import (
    ConversionError,
    convert_to_pdb,
    is_readable_structure,
    needs_conversion,
    pdb_name_for,
)

PACKAGE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CIF = os.path.join(PACKAGE, "data", "examples", "mmcif", "1A2P_A.cif")
PDB_EXAMPLES = os.path.join(PACKAGE, "data", "examples", "barnase", "structures")

needs_fixture = pytest.mark.skipif(not os.path.isfile(CIF),
                                   reason="mmCIF fixture not present")
needs_network = pytest.mark.skipif(
    os.environ.get("WATCON_NETWORK_TESTS") != "1",
    reason="set WATCON_NETWORK_TESTS=1 to run tests that contact RCSB")


def _atom_summary(path):
    """(n_atoms, chains, n_waters) straight from the PDB text."""
    chains, waters, count = set(), 0, 0
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")):
            count += 1
            chains.add(line[21])
            if line[17:20].strip() == "HOH":
                waters += 1
    return count, chains, waters


# ===========================================================================
# Recognising formats -- no gemmi needed
# ===========================================================================

@pytest.mark.parametrize("name,expected", [
    ("1A2P.pdb", False), ("1A2P.ent", False), ("1A2P.pdb.gz", False),
    ("7GSA.cif", True), ("7GSA.mmcif", True), ("7GSA.cif.gz", True),
    ("7GSA.CIF", True), ("notes.txt", False),
])
def test_needs_conversion(name, expected):
    assert needs_conversion(name) is expected


@pytest.mark.parametrize("name", ["1A2P.pdb", "7GSA.cif", "7GSA.cif.gz",
                                  "1A2P.pdb.gz", "1A2P.ent"])
def test_readable_structures_are_recognised(name):
    assert is_readable_structure(name)


def test_non_structures_are_not_readable():
    assert not is_readable_structure("README.txt")
    assert not is_readable_structure("results.csv")


@pytest.mark.parametrize("name,expected", [
    ("7GSA.cif", "7GSA.pdb"),
    ("7GSA.cif.gz", "7GSA.pdb"),
    ("1A2P.pdb", "1A2P.pdb"),
    # Regression: the network builder used to truncate at the first dot, so
    # 1A2P.run1.pdb became "1A2P" and its ConSurf file was never found.
    ("1A2P.run1.cif", "1A2P.run1.pdb"),
])
def test_converted_names_keep_their_dots(name, expected):
    assert pdb_name_for(name) == expected


# ===========================================================================
# Converting
# ===========================================================================

@needs_fixture
def test_an_mmcif_becomes_readable_pdb(tmp_path):
    out = convert_to_pdb(CIF, str(tmp_path / "out.pdb"))
    count, chains, waters = _atom_summary(out)
    assert count > 500, "conversion produced almost no atoms"
    assert chains == {"A"}
    assert waters > 0, "waters were lost in conversion -- they are the point"


@needs_fixture
def test_the_waters_survive(tmp_path):
    """WatCon exists to analyse waters; a converter that drops them is useless."""
    out = convert_to_pdb(CIF, str(tmp_path / "out.pdb"))
    _, _, waters = _atom_summary(out)
    assert waters >= 100


@needs_fixture
def test_a_gzipped_mmcif_converts(tmp_path):
    gz = tmp_path / "1A2P_A.cif.gz"
    with open(CIF, "rb") as source, gzip.open(str(gz), "wb") as target:
        shutil.copyfileobj(source, target)
    out = convert_to_pdb(str(gz), str(tmp_path / "out.pdb"))
    assert _atom_summary(out)[0] > 500


def test_a_plain_pdb_is_copied_not_round_tripped(tmp_path):
    """Reformatting a file WatCon already reads could only introduce a change."""
    source = os.path.join(PDB_EXAMPLES, "1A2P.pdb")
    out = convert_to_pdb(source, str(tmp_path / "same.pdb"))
    assert open(out, "rb").read() == open(source, "rb").read()


def test_a_gzipped_pdb_is_decompressed(tmp_path):
    source = os.path.join(PDB_EXAMPLES, "1A2P.pdb")
    gz = tmp_path / "1A2P.pdb.gz"
    with open(source, "rb") as s, gzip.open(str(gz), "wb") as t:
        shutil.copyfileobj(s, t)
    out = convert_to_pdb(str(gz), str(tmp_path / "out.pdb"))
    assert open(out, "rb").read() == open(source, "rb").read()


def test_rubbish_named_cif_fails_naming_the_file(tmp_path):
    bad = tmp_path / "broken.cif"
    bad.write_text("this is not an mmCIF file at all\n")
    with pytest.raises(ConversionError) as excinfo:
        convert_to_pdb(str(bad), str(tmp_path / "out.pdb"))
    assert "broken.cif" in str(excinfo.value)


# ===========================================================================
# Through prepare, which is the point
# ===========================================================================

@needs_fixture
def test_prepare_accepts_a_folder_containing_mmcif(tmp_path):
    """The end-to-end claim: an mmCIF input prepares like any other."""
    pytest.importorskip("numpy")
    from WatCon.prepare import prepare_directory

    raw = tmp_path / "raw"
    raw.mkdir()
    shutil.copyfile(CIF, str(raw / "1A2P_A.cif"))
    shutil.copyfile(os.path.join(PDB_EXAMPLES, "1BRS.pdb"), str(raw / "1BRS.pdb"))

    report = prepare_directory(str(raw), str(tmp_path / "prepared"),
                               reference="1A2P_A", verbose=False)
    assert len(report.prepared) == 2
    assert report.n_waters > 0
    # Output is PDB whatever went in.
    written = sorted(os.listdir(str(tmp_path / "prepared")))
    assert all(f.endswith(".pdb") for f in written), written


@needs_fixture
def test_the_temporary_conversion_directory_is_cleaned_up(tmp_path):
    """A run that leaves converted copies behind would fill a disk quietly."""
    import tempfile

    from WatCon.prepare import prepare_directory

    raw = tmp_path / "raw"
    raw.mkdir()
    shutil.copyfile(CIF, str(raw / "1A2P_A.cif"))

    before = {d for d in os.listdir(tempfile.gettempdir())
              if d.startswith("watcon_convert_")}
    prepare_directory(str(raw), str(tmp_path / "prepared"), verbose=False)
    after = {d for d in os.listdir(tempfile.gettempdir())
             if d.startswith("watcon_convert_")}
    assert after == before


# ===========================================================================
# Fetching
# ===========================================================================

def test_a_bad_id_is_rejected_without_contacting_rcsb(tmp_path):
    from WatCon.fetch import FetchError, fetch_structure

    with pytest.raises(FetchError) as excinfo:
        fetch_structure("not-an-id", str(tmp_path))
    assert "PDB id" in str(excinfo.value)


def test_an_already_downloaded_file_is_reused(tmp_path):
    """Re-running must be cheap, and must work offline."""
    from WatCon.fetch import fetch_structure

    existing = tmp_path / "1AAX.pdb"
    existing.write_text("ATOM\nEND\n")
    assert fetch_structure("1aax", str(tmp_path)) == str(existing)


def test_fetching_nothing_says_so(tmp_path):
    from WatCon.fetch import FetchError, fetch_structures

    with pytest.raises(FetchError):
        fetch_structures([], str(tmp_path))


@needs_network
def test_an_mmcif_only_entry_is_fetched_as_cif(tmp_path):
    """7GSA has no PDB file. This is the case that motivated the whole step."""
    from WatCon.fetch import fetch_structure

    path = fetch_structure("7GSA", str(tmp_path))
    assert path.endswith(".cif")


@needs_network
def test_one_bad_id_does_not_lose_the_good_ones(tmp_path):
    from WatCon.fetch import fetch_structures

    paths, failures = fetch_structures(["1AAX", "ZZZZ"], str(tmp_path),
                                       verbose=False)
    assert len(paths) == 1
    assert [pdb_id for pdb_id, _ in failures] == ["ZZZZ"]

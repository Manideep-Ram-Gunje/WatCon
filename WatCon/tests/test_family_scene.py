"""The family PyMOL session: it must load what it colours, and open.

Same three traps as the single-protein session, each of which shipped once:
a ``.pml`` with no ``load`` (an empty window, every command succeeding), a
semicolon inside a ``#`` comment (PyMOL splits on it and ran prose as Python),
and a group left enabled while its members were disabled.

These run without PyMOL, on the real five-protein fixtures. The headless PyMOL
check lives beside them and is skipped where PyMOL is absent, as in
``test_view.py``.
"""

from __future__ import annotations

import os
import shutil
import subprocess

import pytest

pytest.importorskip("Bio", reason="family mapping needs Biopython")
pytest.importorskip("MDAnalysis", reason="building networks needs MDAnalysis")

from WatCon.family import FamilyProtein, FamilyStructure, build_family_conservation
from WatCon.family_scene import write_family_session
from WatCon.family_sites import build_family_sites

from .conftest import ptp_grades

PACKAGE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FAMILY_DIR = os.path.join(PACKAGE, "data", "examples", "ptp_family")
ALIGNMENT = os.path.join(FAMILY_DIR, "ptp_family_alignment.pir")

SITES = {
    "PTPN1": ("1AAX", "2F71", "2F71_A_site.pdb"),
    "PTPN6": ("4GRZ", "4GRZ", "4GRZ_A_site.pdb"),
    "PTPN7": ("1ZC0", "1ZC0", "1ZC0_A_site.pdb"),
    "PTPN12": ("5HDE", "5HDE", "5HDE_A_csp_site.pdb"),
    "PTPN22": ("3BRH", "3BRH", "3BRH_A_site.pdb"),
}
PAIRS = {
    "PTPN1": ("2F71", "8U1E"), "PTPN6": ("4GRZ", "4HJP"), "PTPN7": ("1ZC0", "3O4U"),
    "PTPN12": ("5HDE", "5J8R"), "PTPN22": ("3BRH", "3OLR"),
}

pytestmark = pytest.mark.skipif(
    not all(os.path.isfile(os.path.join(FAMILY_DIR, f)) for _r, _c, f in SITES.values()),
    reason="PTP family site fixtures absent")

HAVE_PYMOL = shutil.which("pymol") is not None


@pytest.fixture(scope="module")
def session(tmp_path_factory):
    family = build_family_conservation(
        [FamilyProtein(name, str(ptp_grades(SITES[name][0])),
                       [FamilyStructure(p, os.path.join(FAMILY_DIR, "%s_A_ca.pdb" % p))
                        for p in pair], reference=SITES[name][1])
         for name, pair in PAIRS.items()],
        ALIGNMENT)
    out = tmp_path_factory.mktemp("family_session")
    sites = build_family_sites(
        [FamilyProtein(name, str(ptp_grades(run)),
                       [FamilyStructure(closed, os.path.join(FAMILY_DIR, fixture))],
                       reference=closed)
         for name, (run, closed, fixture) in SITES.items()],
        family, reference="2F71", out_dir=str(out / "superposed"), min_cluster_samples=2)
    return sites, write_family_session(sites, out_dir=str(out))


# ===========================================================================
# The script
# ===========================================================================

def test_it_loads_what_it_colours(session):
    _sites, pml = session
    loads = [l for l in open(pml).read().splitlines() if l.startswith("load ")]
    assert any(", sites" in l for l in loads)
    assert sum(1 for l in loads if l.endswith(("2F71", "4GRZ", "1ZC0", "5HDE", "3BRH"))) == 5


def test_every_loaded_path_exists(session):
    _sites, pml = session
    for line in open(pml):
        if line.startswith("load "):
            path = line[len("load "):].split(",")[0].strip()
            assert os.path.isfile(path), path


def test_no_line_contains_a_semicolon(session):
    """PyMOL splits on ';' even inside a '#' comment."""
    _sites, pml = session
    assert [l for l in open(pml).read().splitlines() if ";" in l] == []


def test_the_shared_sites_are_a_subset_of_all_sites(session):
    sites, pml = session
    directory = os.path.dirname(pml)

    def coordinates(name):
        return {line[30:54] for line in open(os.path.join(directory, name))
                if line.startswith("ATOM")}

    everything = coordinates("family_sites_all.pdb")
    shared = coordinates("family_sites_shared.pdb")
    assert shared
    assert shared <= everything
    assert len(everything) == len(sites.sites)


def test_columns_carry_occupancy_and_grade(session):
    _sites, pml = session
    rows = [l for l in open(os.path.join(os.path.dirname(pml), "family_sites_all.pdb"))
            if l.startswith("ATOM")]
    assert rows
    for line in rows:
        occupancy, bfactor = float(line[54:60]), float(line[60:66])
        assert 0.0 < occupancy <= 1.0
        assert 0.0 <= bfactor <= 9.0


def test_sphere_size_follows_how_many_proteins_hold_a_water(session):
    _sites, pml = session
    text = open(pml).read()
    assert "alter sites, vdw=" in text and "*q" in text


# ===========================================================================
# It opens
# ===========================================================================

@pytest.mark.skipif(not HAVE_PYMOL, reason="PyMOL is not installed")
def test_pymol_opens_it_and_shows_atoms(session):
    _sites, pml = session
    probe = ("print('OBJECTS', ','.join(cmd.get_names()));"
             "print('VISIBLE', cmd.count_atoms('visible'))")
    result = subprocess.run([shutil.which("pymol"), "-cq", str(pml), "-d", probe],
                            capture_output=True, text=True, timeout=600)
    output = result.stdout + result.stderr
    assert "Error" not in output and "SyntaxError" not in output, output[-800:]
    objects = next((l for l in output.splitlines() if l.startswith("OBJECTS ")), "")
    visible = next((l for l in output.splitlines() if l.startswith("VISIBLE ")), "")
    assert "sites_shared" in objects, output[-500:]
    assert int(visible.split()[1]) > 100, output[-500:]

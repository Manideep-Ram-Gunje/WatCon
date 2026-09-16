"""The plugin's Family tab.

Same split as the rest of the plugin: every decision is made by
:mod:`WatCon.family` and :mod:`WatCon.family_sites`, which are tested without
PyMOL; the tab only collects paths and draws results. So these tests check the
collecting and the drawing -- above all that the tab refuses input the analysis
could not honour, before starting a worker thread that would fail minutes later.

PyMOL and Qt are needed here, so the whole module is skipped where they are
absent, as CI is.
"""

from __future__ import annotations

import os

import pytest

pytest.importorskip("Bio", reason="family mapping needs Biopython")

try:
    import pymol                            # noqa: F401
    from pymol.Qt import QtWidgets          # noqa: F401
    HAVE_PYMOL_QT = True
except Exception:                           # noqa: BLE001
    HAVE_PYMOL_QT = False

pytestmark = pytest.mark.skipif(not HAVE_PYMOL_QT, reason="needs PyMOL with Qt")

PACKAGE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FAMILY_DIR = os.path.join(PACKAGE, "data", "examples", "ptp_family")
FIXTURES = os.path.join(PACKAGE, "data", "consurf", "fixtures")
ALIGNMENT = os.path.join(FAMILY_DIR, "ptp_family_alignment.pir")

MEMBERS = {
    "PTPN1": ("1AAX", "2F71", ("2F71", "8U1E")),
    "PTPN6": ("4GRZ", "4GRZ", ("4GRZ", "4HJP")),
    "PTPN7": ("1ZC0", "1ZC0", ("1ZC0", "3O4U")),
    "PTPN12": ("5HDE", "5HDE", ("5HDE", "5J8R")),
    "PTPN22": ("3BRH", "3BRH", ("3BRH", "3OLR")),
}
SITE_FIXTURES = {
    "PTPN1": ("1AAX", "2F71", "2F71_A_site.pdb"),
    "PTPN6": ("4GRZ", "4GRZ", "4GRZ_A_site.pdb"),
    "PTPN7": ("1ZC0", "1ZC0", "1ZC0_A_site.pdb"),
    "PTPN12": ("5HDE", "5HDE", "5HDE_A_csp_site.pdb"),
    "PTPN22": ("3BRH", "3BRH", "3BRH_A_site.pdb"),
}


@pytest.fixture(scope="module")
def app():
    from pymol.Qt import QtWidgets

    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


@pytest.fixture
def dialog(app):
    from WatCon.pymol_plugin.dialog import WatConDialog

    return WatConDialog()


def _fill_members(dialog, tmp_path, use_sites=False):
    """Lay structures out per protein and fill the members table, as a user would."""
    import shutil

    from pymol.Qt import QtWidgets

    for protein, entry in (SITE_FIXTURES if use_sites else MEMBERS).items():
        run, reference = entry[0], entry[1]
        directory = tmp_path / protein
        directory.mkdir(exist_ok=True)
        sources = [entry[2]] if use_sites else ["%s_A_ca.pdb" % p for p in entry[2]]
        ids = [reference] if use_sites else list(entry[2])
        for pdb_id, source in zip(ids, sources):
            shutil.copyfile(os.path.join(FAMILY_DIR, source), str(directory / ("%s.pdb" % pdb_id)))
        row = dialog.members_table.rowCount()
        dialog.members_table.insertRow(row)
        for column, value in enumerate((protein, str(directory),
                                        os.path.join(FIXTURES, "%s_A.grades.txt" % run),
                                        reference)):
            dialog.members_table.setItem(row, column, QtWidgets.QTableWidgetItem(value))
    dialog.alignment_edit.setText(ALIGNMENT)


# ===========================================================================
# The tab exists and does not disturb the original
# ===========================================================================

def test_the_window_has_both_tabs(dialog):
    assert dialog.tabs.count() == 2
    assert dialog.tabs.tabText(0) == "One protein"
    assert dialog.tabs.tabText(1) == "Family"


def test_the_run_button_says_what_it_will_run(dialog):
    dialog.tabs.setCurrentIndex(1)
    assert "family" in dialog.run_button.text().lower()
    dialog.tabs.setCurrentIndex(0)
    assert dialog.run_button.text() == "Run"


def test_the_single_protein_tab_still_works(dialog, tmp_path):
    """The original path must be untouched by the restructuring."""
    with pytest.raises(ValueError, match="folder of structures"):
        dialog._options()


# ===========================================================================
# Refusing input the analysis could not honour
# ===========================================================================

def test_one_member_is_not_a_family(dialog, tmp_path):
    _fill_members(dialog, tmp_path)
    while dialog.members_table.rowCount() > 1:
        dialog.members_table.removeRow(1)
    with pytest.raises(ValueError, match="at least two proteins"):
        dialog._family_options()


def test_an_alignment_is_required(dialog, tmp_path):
    _fill_members(dialog, tmp_path)
    dialog.alignment_edit.setText("")
    with pytest.raises(ValueError, match="alignment"):
        dialog._family_options()


def test_a_missing_folder_is_named(dialog, tmp_path):
    from pymol.Qt import QtWidgets

    _fill_members(dialog, tmp_path)
    dialog.members_table.setItem(0, 1, QtWidgets.QTableWidgetItem(str(tmp_path / "gone")))
    with pytest.raises(ValueError, match="no such folder"):
        dialog._family_options()


def test_a_reference_not_among_the_structures_is_refused(dialog, tmp_path):
    from pymol.Qt import QtWidgets

    _fill_members(dialog, tmp_path)
    dialog.members_table.setItem(0, 3, QtWidgets.QTableWidgetItem("9XYZ"))
    with pytest.raises(ValueError, match="9XYZ"):
        dialog._family_options()


def test_a_malformed_state_label_is_refused(dialog, tmp_path):
    _fill_members(dialog, tmp_path)
    dialog.states_edit.setText("3OLR-open")
    with pytest.raises(ValueError, match="3OLR=open"):
        dialog._family_options()


def test_good_input_becomes_five_members(dialog, tmp_path):
    _fill_members(dialog, tmp_path)
    dialog.states_edit.setText("8U1E=open 2F71=closed")
    options = dialog._family_options()
    assert [p.name for p in options["proteins"]] == list(MEMBERS)
    assert all(len(p.structures) == 2 for p in options["proteins"])
    assert options["states"] == {"8U1E": "open", "2F71": "closed"}
    assert options["reference"] == "2F71"


# ===========================================================================
# Drawing the result
# ===========================================================================

@pytest.fixture(scope="module")
def built(tmp_path_factory):
    """A real family result, built the way the worker builds it."""
    from WatCon.family import FamilyProtein, FamilyStructure, build_family_conservation
    from WatCon.family_scene import write_family_session
    from WatCon.family_sites import build_family_sites

    family = build_family_conservation(
        [FamilyProtein(name, os.path.join(FIXTURES, "%s_A.grades.txt" % run),
                       [FamilyStructure(p, os.path.join(FAMILY_DIR, "%s_A_ca.pdb" % p))
                        for p in structures], reference=reference)
         for name, (run, reference, structures) in MEMBERS.items()],
        ALIGNMENT)
    out = tmp_path_factory.mktemp("plugin_family")
    sites = build_family_sites(
        [FamilyProtein(name, os.path.join(FIXTURES, "%s_A.grades.txt" % run),
                       [FamilyStructure(reference, os.path.join(FAMILY_DIR, fixture))],
                       reference=reference)
         for name, (run, reference, fixture) in SITE_FIXTURES.items()],
        family, reference="2F71", out_dir=str(out / "superposed"), min_cluster_samples=2)
    return family, sites, write_family_session(sites, out_dir=str(out))


def test_the_table_fills_and_filters(dialog, built):
    family, sites, _session = built
    dialog._sites = sites
    dialog.shared_only.setChecked(True)
    dialog._fill_family_table()
    shared = dialog.family_table.rowCount()
    dialog.shared_only.setChecked(False)
    dialog._fill_family_table()
    assert 0 < shared < dialog.family_table.rowCount() == len(sites.sites)
    assert "in every protein" in dialog.family_summary.text()


def test_clicking_a_row_before_running_explains_itself(dialog, built):
    from pymol import cmd

    _family, sites, _session = built
    cmd.delete("all")
    dialog._sites = sites
    dialog._fill_family_table()
    dialog.family_table.selectRow(0)
    assert "not loaded" in dialog.status.text()


def test_completion_draws_the_session_and_reports_the_audit(dialog, built):
    from pymol import cmd

    family, sites, session = built
    dialog._on_family_completed(family, sites, session)

    objects = cmd.get_names("objects")
    assert "sites_shared" in objects and "2F71" in objects
    assert cmd.count_atoms("visible") > 100

    audit = dialog.audit.toPlainText()
    assert "identity" in audit
    assert "differs at A215" in audit          # the 1AAX query's C215S, named
    assert "PTPN7" in audit and "124" in audit  # the excluded 3O4U slides

    dialog.family_table.clearSelection()
    dialog.family_table.selectRow(0)
    assert "held by" in dialog.status.text()
    cmd.delete("all")

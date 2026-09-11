"""The PyMOL plugin.

CI has no PyMOL, and neither does most of the world, so the tests are in two
groups:

* those that run **anywhere**, including one that imports the plugin package
  with ``pymol`` deliberately blocked -- because a plugin that cannot be
  imported on a machine without PyMOL cannot be checked by the test suite at
  all, and a broken one would reach users;
* those that need a real PyMOL and a Qt application, skipped otherwise.

The second group is where the interesting bug was found: ``build_scene``
computed everything correctly but never wrote its ``.pml``, so the plugin ran
``@<path>`` on a file that did not exist and failed *after* the analysis had
succeeded. ``watcon view`` never noticed because it wrote the file itself.
"""

from __future__ import annotations

import builtins
import importlib
import os
import sys

import pytest

PACKAGE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLES = os.path.join(PACKAGE, "data", "examples", "barnase", "structures")
GRADES = os.path.join(PACKAGE, "data", "consurf", "fixtures", "1BRS_A_150.grades.txt")

try:
    import pymol                            # noqa: F401
    from pymol.Qt import QtWidgets          # noqa: F401
    HAVE_PYMOL_QT = True
except Exception:                           # noqa: BLE001
    HAVE_PYMOL_QT = False

needs_pymol = pytest.mark.skipif(not HAVE_PYMOL_QT,
                                 reason="needs PyMOL with Qt")
needs_examples = pytest.mark.skipif(not os.path.isdir(EXAMPLES),
                                    reason="bundled example data not present")


# ===========================================================================
# Must be importable with no PyMOL at all
# ===========================================================================

def test_the_plugin_package_imports_without_pymol(monkeypatch):
    """The whole point of keeping every heavy import inside a function.

    PyMOL also imports every module in its startup directory at launch, so a
    plugin that drags in scikit-learn at import time slows every session
    whether or not anyone opens it.
    """
    real_import = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name == "pymol" or name.startswith("pymol."):
            raise ImportError("pymol is blocked for this test")
        return real_import(name, *args, **kwargs)

    for module in [m for m in sys.modules if m.startswith("WatCon.pymol_plugin")]:
        del sys.modules[module]

    monkeypatch.setattr(builtins, "__import__", blocked)
    module = importlib.import_module("WatCon.pymol_plugin")
    assert callable(module.__init_plugin__)
    assert callable(module.run_plugin_gui)
    assert module.MENU_LABEL


def test_nothing_heavy_is_imported_at_module_scope():
    import inspect

    from WatCon import pymol_plugin

    source = inspect.getsource(pymol_plugin)
    top_level = [line for line in source.splitlines()
                 if line.startswith(("import ", "from ")) and "__future__" not in line]
    assert top_level == [], top_level


def test_the_dialog_module_is_not_imported_until_asked():
    """Opening the dialog is what pulls in Qt, MDAnalysis and scikit-learn."""
    import inspect

    from WatCon import pymol_plugin

    source = inspect.getsource(pymol_plugin.run_plugin_gui)
    assert "from .dialog import" in source


# ===========================================================================
# The installer
# ===========================================================================

def test_the_shim_is_valid_python_and_imports_the_plugin():
    import ast

    from WatCon.cli import PLUGIN_SHIM

    ast.parse(PLUGIN_SHIM)
    assert "from WatCon.pymol_plugin import __init_plugin__" in PLUGIN_SHIM


def test_the_shim_is_a_shim_not_a_copy():
    """A copied plugin goes stale the moment WatCon is upgraded."""
    from WatCon.cli import PLUGIN_SHIM

    assert len(PLUGIN_SHIM.splitlines()) < 15
    assert "QtWidgets" not in PLUGIN_SHIM


@needs_pymol
def test_install_and_uninstall_round_trip(tmp_path, monkeypatch, capsys):
    from WatCon import cli

    monkeypatch.setattr(cli, "plugin_directory",
                        lambda: (str(tmp_path), "a test directory"))
    target = tmp_path / cli.PLUGIN_FILENAME

    parser = cli.build_parser()
    assert cli.cmd_plugin(parser.parse_args(["plugin", "--install"])) == 0
    assert target.is_file()

    # Installing twice must not silently overwrite, nor fail.
    assert cli.cmd_plugin(parser.parse_args(["plugin", "--install"])) == 0
    assert "Already installed" in capsys.readouterr().out

    assert cli.cmd_plugin(parser.parse_args(["plugin", "--uninstall"])) == 0
    assert not target.exists()


@needs_pymol
def test_the_chosen_directory_is_one_pymol_scans():
    from pymol.plugins import get_startup_path

    from WatCon.cli import plugin_directory

    directory, explanation = plugin_directory()
    assert explanation
    assert directory in get_startup_path() or os.path.isdir(directory)


@needs_pymol
def test_pymol_would_discover_an_installed_plugin(tmp_path):
    """`findPlugins` is what PyMOL uses at startup, so ask it directly."""
    from pymol.plugins import findPlugins

    from WatCon.cli import PLUGIN_FILENAME, PLUGIN_SHIM

    (tmp_path / PLUGIN_FILENAME).write_text(PLUGIN_SHIM)
    found = findPlugins([str(tmp_path)])
    assert os.path.splitext(PLUGIN_FILENAME)[0] in found


# ===========================================================================
# The dialog
# ===========================================================================

@pytest.fixture(scope="module")
def app():
    from pymol.Qt import QtWidgets

    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


@pytest.fixture(scope="module")
def prepared(tmp_path_factory):
    from WatCon.prepare import prepare_directory

    root = tmp_path_factory.mktemp("plugin")
    prepare_directory(EXAMPLES, str(root / "prepared"), reference="1A2P",
                      verbose=False)
    return str(root / "prepared")


@needs_pymol
@needs_examples
def test_the_dialog_constructs(app):
    from WatCon.pymol_plugin.dialog import WatConDialog

    dialog = WatConDialog()
    assert dialog.windowTitle() == "WatCon + ConSurf"


@needs_pymol
@needs_examples
def test_the_reference_and_chain_dropdowns_fill_themselves(app, prepared):
    """The user asked to pick the structure and the chain, not type them."""
    from WatCon.pymol_plugin.dialog import WatConDialog

    dialog = WatConDialog()
    dialog.structures_edit.setText(prepared)
    dialog._refresh_structures()

    references = [dialog.reference_combo.itemText(i)
                  for i in range(dialog.reference_combo.count())]
    assert references == ["1A2P", "1BRN", "1BRS", "1BSA", "1BSE", "1RNB"]
    chains = [dialog.chain_combo.itemText(i)
              for i in range(dialog.chain_combo.count())]
    assert chains == ["A"], "prepared structures are relabelled to one chain"


@needs_pymol
@needs_examples
def test_missing_input_is_refused_with_a_reason(app, prepared):
    from WatCon.pymol_plugin.dialog import WatConDialog

    dialog = WatConDialog()
    with pytest.raises(ValueError, match="folder of structures"):
        dialog._options()

    dialog.structures_edit.setText(prepared)
    with pytest.raises(ValueError, match="ConSurf"):
        dialog._options()


@needs_pymol
@needs_examples
def test_the_dialog_does_not_offer_to_contact_consurf(app, prepared):
    """WatCon never talks to the ConSurf server, and must not imply it does."""
    from WatCon.pymol_plugin.dialog import WatConDialog

    dialog = WatConDialog()
    dialog.structures_edit.setText(prepared)
    try:
        dialog._options()
    except ValueError as error:
        assert "does not contact the ConSurf server" in str(error)


@needs_pymol
@needs_examples
def test_the_table_fills_from_a_scene_and_can_be_filtered(app, prepared, tmp_path):
    from WatCon.pymol_plugin.dialog import WatConDialog
    from WatCon.scene import build_scene

    scene = build_scene(prepared, GRADES, out_dir=str(tmp_path / "view"),
                        reference="1A2P", verbose=False)
    dialog = WatConDialog()
    dialog._scene = scene

    dialog.conserved_only.setChecked(True)
    dialog._fill_table()
    assert dialog.table.rowCount() == len(scene.conserved_sites)

    dialog.conserved_only.setChecked(False)
    dialog._fill_table()
    assert dialog.table.rowCount() == len(scene.sites)
    assert "clusters" in dialog.summary.text()


@needs_pymol
@needs_examples
def test_clicking_a_row_before_running_explains_itself(app, prepared, tmp_path):
    """The table outlives `delete all`, so this is reachable, not theoretical."""
    from pymol import cmd

    from WatCon.pymol_plugin.dialog import WatConDialog
    from WatCon.scene import build_scene

    cmd.delete("all")
    scene = build_scene(prepared, GRADES, out_dir=str(tmp_path / "view"),
                        reference="1A2P", verbose=False)
    dialog = WatConDialog()
    dialog._scene = scene
    dialog._fill_table()
    dialog.table.selectRow(0)
    assert "not loaded" in dialog.status.text()


@needs_pymol
@needs_examples
def test_running_the_scene_draws_it_and_focusing_works(app, prepared, tmp_path):
    """The end-to-end claim, through the plugin's own code path."""
    from pymol import cmd

    from WatCon.pymol_plugin.dialog import WatConDialog
    from WatCon.scene import build_scene

    scene = build_scene(prepared, GRADES, out_dir=str(tmp_path / "view"),
                        reference="1A2P", verbose=False)
    dialog = WatConDialog()
    dialog._on_completed(scene)

    objects = cmd.get_names("objects")
    assert "protein" in objects and "sites" in objects
    assert cmd.count_atoms("visible") > 500

    dialog.table.clearSelection()
    dialog.table.selectRow(0)
    assert cmd.count_atoms("watcon_site") == 1
    assert "Lined by" in dialog.status.text()

    cmd.delete("all")


# ===========================================================================
# The bug the plugin found
# ===========================================================================

@needs_examples
def test_build_scene_writes_its_pml(tmp_path):
    """Regression: `pml_path` promised a file that build_scene never wrote.

    `watcon view` never noticed because it wrote the file itself. The plugin
    ran `@<path>` and raised FileNotFoundError after the analysis had already
    succeeded -- the worst moment to fail.
    """
    pytest.importorskip("MDAnalysis")
    from WatCon.prepare import prepare_directory
    from WatCon.scene import build_scene

    prepare_directory(EXAMPLES, str(tmp_path / "prepared"), reference="1A2P",
                      verbose=False)
    scene = build_scene(str(tmp_path / "prepared"), GRADES,
                        out_dir=str(tmp_path / "view"), reference="1A2P",
                        verbose=False)
    for path in (scene.pml_path, scene.sites_all_path,
                 scene.sites_conserved_path, scene.report_path):
        assert os.path.isfile(path), path

"""WatCon + ConSurf as a PyMOL plugin.

Everything the tool does was reachable only from a terminal: change the
reference structure or a cutoff and you re-ran a command and re-opened a file.
PyMOL is already open in front of the user, so the controls belong there too.

Installing
----------
``watcon plugin --install`` writes a three-line shim into PyMOL's startup
directory that imports this package. Nothing is copied, so upgrading WatCon
upgrades the plugin.

It works because PyMOL and WatCon share an interpreter in a normal pip install
-- the plugin simply ``import``\\ s WatCon. Where they do not share one, install
WatCon into the Python that PyMOL uses.

Why the imports are where they are
----------------------------------
Nothing at module scope imports PyMOL, Qt, MDAnalysis or anything else heavy.
Two reasons, both load-bearing:

* **CI has no PyMOL.** ``import WatCon.pymol_plugin`` must succeed anyway, or
  the package's own test suite cannot check that this file is well formed --
  which is how a broken plugin would reach users.
* PyMOL imports every module in its startup directory at launch. A plugin that
  drags in scikit-learn at import time slows every PyMOL session, whether or not
  anyone opens it.
"""

from __future__ import annotations

__all__ = ["__init_plugin__", "run_plugin_gui"]

#: Shown in PyMOL's Plugin menu.
MENU_LABEL = "WatCon + ConSurf"

_dialog = None


def __init_plugin__(app=None):
    """Called by PyMOL at startup. Adds the menu entry and returns."""
    from pymol.plugins import addmenuitemqt

    addmenuitemqt(MENU_LABEL, run_plugin_gui)


def run_plugin_gui():
    """Open the dialog, reusing it if it is already open.

    Kept as a module-level reference: a QDialog with no surviving reference is
    garbage-collected and vanishes from the screen, which looks like a crash.
    """
    global _dialog

    from .dialog import WatConDialog

    if _dialog is None:
        _dialog = WatConDialog()
    _dialog.show()
    _dialog.raise_()
    _dialog.activateWindow()
    return _dialog

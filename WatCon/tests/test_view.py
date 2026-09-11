"""The session must actually open and show something.

`pymol_project_evolutionary` shipped writing a .pml with no ``load`` command:
~200 colour commands against whatever happened to be open. Run on its own it
produced an empty window --

    objects loaded: []
    atoms visible: 0

-- and every command succeeded, so nothing reported a problem. The data was
right the whole time; only the viewing was broken, and no test looked.

These tests check the two things that were never checked: that the script loads
what it colours, and that opening it leaves atoms on screen.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile

import pytest

pytest.importorskip("MDAnalysis", reason="building a session needs MDAnalysis")

from WatCon.view import HIGHLY_CONSERVED, build_session

PACKAGE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLES = os.path.join(PACKAGE, "data", "examples", "barnase", "structures")
GRADES = os.path.join(PACKAGE, "data", "consurf", "fixtures", "1BRS_A_150.grades.txt")

needs_examples = pytest.mark.skipif(
    not os.path.isdir(EXAMPLES), reason="bundled example data not present")

HAVE_PYMOL = shutil.which("pymol") is not None


@pytest.fixture(scope="module")
def session(tmp_path_factory):
    """Build one session from the bundled example, as `watcon view` does."""
    from WatCon.prepare import prepare_directory

    root = tmp_path_factory.mktemp("view")
    prepare_directory(EXAMPLES, str(root / "prepared"), reference="1A2P",
                      verbose=False)
    pml = build_session(str(root / "prepared"), GRADES,
                        out_dir=str(root / "session"), verbose=False)
    return root / "session", pml


# ===========================================================================
# The script is self-contained
# ===========================================================================

@needs_examples
def test_the_script_loads_what_it_colours(session):
    """Regression: the original .pml had no `load` at all."""
    _, pml = session
    text = open(pml).read()
    loads = [line for line in text.splitlines() if line.startswith("load ")]
    assert loads, "the script colours things it never loads"
    assert any("protein" in line for line in loads)
    assert any("sites" in line for line in loads)


@needs_examples
def test_every_loaded_path_exists(session):
    """A load line pointing at a missing file fails silently in PyMOL."""
    _, pml = session
    for line in open(pml):
        if line.startswith("load "):
            path = line[len("load "):].split(",")[0].strip()
            assert os.path.isfile(path), "load target missing: %s" % path


@needs_examples
def test_the_outputs_are_written(session):
    out, _ = session
    for name in ("watcon_view.pml", "sites_all.pdb", "sites_conserved.pdb",
                 "conservation.csv"):
        assert (out / name).is_file(), name


@needs_examples
def test_grades_reach_the_bfactor_column(session):
    """The grade is what colours the spheres, so it must be in the file."""
    out, _ = session
    grades = set()
    for line in open(out / "sites_all.pdb"):
        if line.startswith("ATOM"):
            grades.add(round(float(line[60:66])))
    assert grades - {0}, "no site carries a ConSurf grade"
    assert max(grades) <= 9 and min(grades) >= 0


@needs_examples
def test_the_conserved_subset_is_a_subset(session):
    """sites_conserved must be drawn from sites, at grade >= the cutoff."""
    out, _ = session

    def coords(path):
        return {(line[30:38], line[38:46], line[46:54])
                for line in open(path) if line.startswith("ATOM")}

    def grades(path):
        return [round(float(line[60:66]))
                for line in open(path) if line.startswith("ATOM")]

    assert coords(out / "sites_conserved.pdb") <= coords(out / "sites_all.pdb")
    assert all(g >= HIGHLY_CONSERVED for g in grades(out / "sites_conserved.pdb"))


# ===========================================================================
# It opens
# ===========================================================================

@needs_examples
@pytest.mark.skipif(not HAVE_PYMOL, reason="PyMOL is not installed")
def test_pymol_opens_it_and_shows_atoms(session):
    """The check whose absence let a file that displays nothing ship.

    Runs PyMOL headlessly on the generated script and asserts that objects
    loaded and that atoms are actually visible.
    """
    _, pml = session
    # No %-formatting here: PyMOL's -d consumes '%' for its own substitution,
    # so "print('X=%s' % y)" reaches Python as "print('X= s'   y)" and raises a
    # SyntaxError that looks like a broken session rather than a broken probe.
    probe = (
        "print('OBJECTS', ','.join(cmd.get_names()));"
        "print('VISIBLE', cmd.count_atoms('visible'))"
    )
    result = subprocess.run(
        [shutil.which("pymol"), "-cq", str(pml), "-d", probe],
        capture_output=True, text=True, timeout=600,
    )
    output = result.stdout + result.stderr
    objects = next((l for l in output.splitlines()
                    if l.startswith("OBJECTS ")), "")
    visible = next((l for l in output.splitlines()
                    if l.startswith("VISIBLE ")), "")

    assert "protein" in objects, output[-800:]
    assert "sites" in objects, output[-800:]

    count = int(visible.split()[1])
    assert count > 0, "the session opens but shows nothing: %s" % output[-800:]
    assert count > 500, "expected the protein to be visible, saw %d atoms" % count


@needs_examples
@pytest.mark.skipif(not HAVE_PYMOL, reason="PyMOL is not installed")
def test_a_script_without_a_structure_is_marked_as_inert(tmp_path):
    """Colour-only output is still allowed -- but must say what it is.

    `pymol_project_evolutionary` without a structure produces a file that
    colours an already-open session. That is a legitimate use; silently
    producing a file that shows nothing is not.
    """
    from WatCon.visualize_structures import pymol_project_evolutionary

    class Atom:
        def __init__(self):
            self.chain, self.resid, self.icode = "A", 1, None
            self.evolutionary = None

    class Net:
        protein_atoms = [Atom()]

    path = pymol_project_evolutionary(Net(), filename="bare.pml",
                                      out_path=str(tmp_path))
    text = open(path).read()
    assert "load" not in text.split("bg white")[0].replace("Load your", "")
    assert "NOTE" in text and "show nothing" in text

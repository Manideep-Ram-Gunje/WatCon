"""The scene, checked without PyMOL running.

``watcon view`` writes a ``.pml`` and the PyMOL plugin drives a live session.
Both consume :func:`WatCon.scene.build_scene`, so this is where the picture is
actually pinned down -- and it is pinned here rather than in ``test_view.py``
because **CI has no PyMOL**. A scene bug that only a running PyMOL could catch
would reach users.

The semicolon test below is not hypothetical. PyMOL treats ``;`` as a command
separator *even inside a ``#`` comment*, so this line::

    # WatCon chose the residues; PyMOL only decides which atom pair each

executed everything after the semicolon as Python and raised ``SyntaxError:
invalid syntax`` in the middle of an otherwise working session.
"""

from __future__ import annotations

import os

import pytest

pytest.importorskip("MDAnalysis", reason="building a scene needs MDAnalysis")

from WatCon.scene import (
    CONTACT_CUTOFF,
    HIGHLY_CONSERVED,
    NO_DATA_COLOR,
    Scene,
    Site,
    build_scene,
)

PACKAGE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLES = os.path.join(PACKAGE, "data", "examples", "barnase", "structures")
GRADES = os.path.join(PACKAGE, "data", "consurf", "fixtures", "1BRS_A_150.grades.txt")

needs_examples = pytest.mark.skipif(
    not os.path.isdir(EXAMPLES), reason="bundled example data not present")


@pytest.fixture(scope="module")
def scene(tmp_path_factory):
    """One scene from the bundled example, as `watcon view` builds it."""
    from WatCon.prepare import prepare_directory

    root = tmp_path_factory.mktemp("scene")
    prepare_directory(EXAMPLES, str(root / "prepared"), reference="1A2P",
                      verbose=False)
    built = build_scene(str(root / "prepared"), GRADES,
                        out_dir=str(root / "session"), reference="1A2P",
                        verbose=False)
    built.write_pml()
    return built


# ===========================================================================
# Commands PyMOL can actually run
# ===========================================================================

@needs_examples
def test_no_command_contains_a_semicolon(scene):
    """Regression: PyMOL splits on ';' even inside a '#' comment.

    A semicolon in a prose comment ran the rest of the sentence as Python.
    """
    offenders = [line for line in scene.commands if ";" in line]
    assert not offenders, offenders


@needs_examples
def test_every_loaded_path_exists(scene):
    """A load line pointing at a missing file fails silently in PyMOL."""
    for line in scene.commands:
        if line.startswith("load "):
            path = line[len("load "):].split(",")[0].strip()
            assert os.path.isfile(path), "load target missing: %s" % path


@needs_examples
def test_the_scene_loads_what_it_colours(scene):
    """Regression: the original .pml had no `load` at all."""
    loads = [line for line in scene.commands if line.startswith("load ")]
    assert any("protein" in line for line in loads)
    assert any("sites" in line for line in loads)


@needs_examples
def test_the_scene_is_deterministic(scene, tmp_path):
    """Same input, same picture -- or a screenshot cannot be reproduced."""
    from WatCon.prepare import prepare_directory

    prepare_directory(EXAMPLES, str(tmp_path / "prepared"), reference="1A2P",
                      verbose=False)
    again = build_scene(str(tmp_path / "prepared"), GRADES,
                        out_dir=str(tmp_path / "session"), reference="1A2P",
                        verbose=False)
    assert [c for c in again.commands if not c.startswith(("load ", "#"))] == \
           [c for c in scene.commands if not c.startswith(("load ", "#"))]


@needs_examples
def test_writing_the_pml_round_trips(scene, tmp_path):
    path = scene.write_pml(str(tmp_path / "copy.pml"))
    assert open(path).read().splitlines() == scene.commands


# ===========================================================================
# What the picture claims
# ===========================================================================

@needs_examples
def test_the_bundled_example_reproduces_its_numbers(scene):
    """The figures quoted in the README and changelog, asserted."""
    assert scene.n_clusters == 193
    assert len(scene.sites) == 190
    assert scene.n_with_conservation == 165
    assert len(scene.conserved_sites) == 57


@needs_examples
def test_no_data_is_yellow_and_not_grey(scene):
    """grey70 was almost the near-white of an average grade.

    'We have no score for this residue' and 'this residue is averagely
    conserved' are opposite claims and must not look alike.
    """
    text = "\n".join(scene.commands)
    assert "consurf_nodata" in text
    assert "set_color consurf_nodata, [%.3f, %.3f, %.3f]" % NO_DATA_COLOR in text
    # Only actual commands -- the comment explaining the change names grey70.
    colours = [c for c in scene.commands if c.startswith("color ")]
    assert not [c for c in colours if "grey70" in c or "grey60" in c], colours
    assert any("consurf_nodata, protein" in c for c in colours)
    assert any("consurf_nodata, sites" in c for c in colours)


@needs_examples
def test_radius_tracks_occupancy(scene):
    """A site seen in every structure must not look like one seen once."""
    alters = [c for c in scene.commands if c.startswith("alter ") and "vdw=" in c]
    assert alters, "sphere radius is constant"
    assert any("*q" in c for c in alters), "radius does not use the occupancy column"


@needs_examples
def test_the_contacts_are_built_but_off(scene):
    """The network is the point of WatCon, but 57 sites' worth buries the view."""
    text = "\n".join(scene.commands)
    assert "distance site_contacts" in text
    assert "%.1f, mode=2" % CONTACT_CUTOFF in text
    assert "group WatCon_contacts" in text
    # The group carries its own enabled flag; disabling only the members leaves
    # `enable WatCon_contacts` appearing to do nothing.
    assert "disable WatCon_contacts" in text


@needs_examples
def test_only_the_intended_objects_are_shown_on_opening(scene):
    text = "\n".join(scene.commands)
    assert "show spheres, sites_conserved" in text
    assert "disable sites" in text


# ===========================================================================
# The site records the plugin's table is built from
# ===========================================================================

@needs_examples
def test_sites_know_which_residues_line_them(scene):
    """`ClusterConservation.residue_keys` held this and nothing surfaced it."""
    with_residues = [s for s in scene.sites if s.residues]
    assert len(with_residues) >= 100
    site = max(scene.sites, key=lambda s: len(s.residues))
    chain, resid, _icode, grade = site.residues[0]
    assert isinstance(resid, int)
    assert grade is None or 1 <= grade <= 9
    assert site.residue_label()


@needs_examples
def test_occupancy_fraction_is_a_fraction(scene):
    for site in scene.sites:
        assert 0.0 < site.occupancy_fraction <= 1.0


@needs_examples
def test_conserved_is_a_subset_at_the_threshold(scene):
    assert all(s.grade >= scene.highly_conserved for s in scene.conserved_sites)
    assert set(id(s) for s in scene.conserved_sites) <= set(id(s) for s in scene.sites)


@needs_examples
def test_a_zero_grade_means_no_data_not_low_conservation(scene):
    unscored = [s for s in scene.sites if s.grade == 0]
    assert unscored, "expected some sites with no ConSurf-scored lining residue"
    assert all(not s.has_conservation for s in unscored)
    assert all(s.grade < scene.highly_conserved for s in unscored)


@needs_examples
def test_the_threshold_is_a_parameter(scene, tmp_path):
    """The dialog exposes it, so it must not be baked in."""
    from WatCon.prepare import prepare_directory

    prepare_directory(EXAMPLES, str(tmp_path / "prepared"), reference="1A2P",
                      verbose=False)
    strict = build_scene(str(tmp_path / "prepared"), GRADES,
                         out_dir=str(tmp_path / "s9"), reference="1A2P",
                         highly_conserved=9, verbose=False)
    assert len(strict.conserved_sites) < len(scene.conserved_sites)
    assert all(s.grade == 9 for s in strict.conserved_sites)


# ===========================================================================
# The site PDB, whose columns PyMOL reads
# ===========================================================================

@needs_examples
def test_grade_and_occupancy_sit_in_the_right_columns(scene):
    """The previous writer put the B-factor one column late.

    It parsed only because the trailing digit fell off the end of the field --
    fine for 9.00, wrong the moment a value changes width.
    """
    rows = [line for line in open(scene.sites_all_path)
            if line.startswith("ATOM")]
    assert rows
    for line in rows:
        occupancy = float(line[54:60])
        bfactor = float(line[60:66])
        assert 0.0 < occupancy <= 1.0
        assert 0.0 <= bfactor <= 9.0
        assert float(bfactor).is_integer()


@needs_examples
def test_the_conserved_pdb_is_a_subset(scene):
    def coords(path):
        return {(line[30:38], line[38:46], line[46:54])
                for line in open(path) if line.startswith("ATOM")}

    assert coords(scene.sites_conserved_path) <= coords(scene.sites_all_path)


@needs_examples
def test_the_remarks_say_what_zero_means(scene):
    text = open(scene.sites_all_path).read()
    assert "NO DATA" in text
    assert "occupancy" in text


# ===========================================================================
# The module must stay importable where PyMOL is not
# ===========================================================================

def test_scene_does_not_import_pymol():
    """CI has no PyMOL. Neither does most of the world."""
    import inspect

    from WatCon import scene as scene_module

    source = inspect.getsource(scene_module)
    assert "import pymol" not in source
    assert "from pymol" not in source


def test_a_site_with_no_structures_has_zero_occupancy_fraction():
    """Guard against dividing by zero in the table."""
    site = Site(cluster_id=1, centre=(0.0, 0.0, 0.0), grade=0, occupancy=0,
                n_structures_occupied=0, n_structures_total=0)
    assert site.occupancy_fraction == 0.0
    assert not site.has_conservation


def test_the_default_threshold_is_the_documented_one():
    assert HIGHLY_CONSERVED == 8

"""What happens when the input directory is not perfectly clean.

Every case here was measured against the code as published, and every message
quoted in a docstring is the message a user actually got. None of them named the
file responsible or said what to do, and one of them blamed ConSurf for a
filename problem.

The tests assert on the *content* of the new messages, not just that something
was raised -- a clear error is the feature, so a vague one is a regression.
"""

from __future__ import annotations

import os

import pytest

from WatCon.residue_index import (
    NoStructuresFound,
    list_structure_files,
    require_structure_files,
)


@pytest.fixture
def folder(tmp_path):
    """A structure folder with the debris a real one accumulates."""
    for name in ("1A2P.pdb", "1BRS.pdb", "1A2P.run1.pdb"):
        (tmp_path / name).write_text("ATOM\nEND\n")
    (tmp_path / "README.txt").write_text("notes\n")
    (tmp_path / ".DS_Store").write_bytes(b"\x00")
    (tmp_path / "1A2P.pdb.swp").write_text("editor backup\n")
    (tmp_path / "results").mkdir()
    return tmp_path


# ===========================================================================
# Choosing what to read
# ===========================================================================

def test_a_readme_is_skipped_not_parsed(folder):
    """Was: ValueError: 'TXT' isn't a valid topology format."""
    structures, skipped = list_structure_files(str(folder))
    assert "README.txt" not in structures
    assert "README.txt" in skipped


def test_a_subdirectory_is_skipped(folder):
    """Was: ValueError: '' isn't a valid topology format."""
    structures, skipped = list_structure_files(str(folder))
    assert "results" not in structures
    assert "results" in skipped


def test_editor_backups_and_dotfiles_are_skipped(folder):
    structures, skipped = list_structure_files(str(folder))
    assert set(skipped) >= {".DS_Store", "1A2P.pdb.swp"}
    assert not any(s.endswith(".swp") for s in structures)


def test_the_structures_are_found(folder):
    structures, _ = list_structure_files(str(folder))
    assert structures == ["1A2P.pdb", "1A2P.run1.pdb", "1BRS.pdb"]


def test_the_order_is_deterministic(folder):
    assert list_structure_files(str(folder)) == list_structure_files(str(folder))


def test_gzipped_structures_are_recognised(tmp_path):
    (tmp_path / "8XYZ.cif.gz").write_bytes(b"\x1f\x8b")
    structures, _ = list_structure_files(str(tmp_path))
    assert structures == ["8XYZ.cif.gz"]


# ===========================================================================
# Saying so when there is nothing to read
# ===========================================================================

def test_an_empty_directory_names_the_problem(tmp_path):
    """Was: ValueError: not enough values to unpack (expected 2, got 0)."""
    with pytest.raises(NoStructuresFound) as excinfo:
        require_structure_files(str(tmp_path))
    message = str(excinfo.value)
    assert str(tmp_path) in message
    assert ".pdb" in message


def test_a_directory_of_non_structures_lists_them(tmp_path):
    (tmp_path / "notes.txt").write_text("x")
    (tmp_path / "table.csv").write_text("x")
    with pytest.raises(NoStructuresFound) as excinfo:
        require_structure_files(str(tmp_path))
    message = str(excinfo.value)
    assert "notes.txt" in message and "table.csv" in message


def test_a_missing_directory_says_so(tmp_path):
    with pytest.raises(NoStructuresFound) as excinfo:
        require_structure_files(str(tmp_path / "nope"))
    assert "No such directory" in str(excinfo.value)


def test_require_returns_what_list_returns(folder):
    assert require_structure_files(str(folder)) == list_structure_files(str(folder))


# ===========================================================================
# Filenames with dots
# ===========================================================================

def test_a_dotted_filename_keeps_its_name(folder):
    """Was: 1A2P.run1.pdb -> "1A2P", so the ConSurf lookup searched for the
    wrong name and failed blaming ConSurf."""
    structures, _ = list_structure_files(str(folder))
    stems = [os.path.splitext(f)[0] for f in structures]
    assert "1A2P.run1" in stems
    assert stems.count("1A2P") == 1


def test_the_builder_no_longer_truncates_at_the_first_dot():
    """Pin the call site itself, not only the helper."""
    import inspect

    from WatCon import generate_static_networks

    source = inspect.getsource(generate_static_networks.initialize_network)
    assert "f.split('.')[0] for f in os.listdir" not in source
    assert "require_structure_files" in source


def test_the_builder_no_longer_reads_everything_in_the_folder():
    import inspect

    from WatCon import generate_static_networks

    source = inspect.getsource(generate_static_networks.initialize_network)
    assert "'swp' not in f" not in source


# ===========================================================================
# Structures with no waters
# ===========================================================================

def test_zero_waters_names_water_not_sklearn():
    """Was: Found array with 0 sample(s) (shape=(0, 3)) while a minimum of 1
    is required by HDBSCAN."""
    from WatCon.find_conserved_networks import (
        NoWaterCoordinates,
        cluster_coordinates_only,
    )

    with pytest.raises(NoWaterCoordinates) as excinfo:
        cluster_coordinates_only([], source=["1A2P", "1BRS"])
    message = str(excinfo.value)
    assert "water" in message.lower()
    assert "1A2P" in message and "1BRS" in message


def test_zero_waters_without_a_source_still_explains():
    from WatCon.find_conserved_networks import (
        NoWaterCoordinates,
        cluster_coordinates_only,
    )

    with pytest.raises(NoWaterCoordinates) as excinfo:
        cluster_coordinates_only([])
    assert "HOH" in str(excinfo.value)


def test_a_ragged_coordinate_array_is_rejected():
    """The old bare `except` printed a warning and carried on with the
    un-reshaped array, so a (N,)-shaped array reached sklearn."""
    from WatCon.find_conserved_networks import (
        NoWaterCoordinates,
        cluster_coordinates_only,
    )

    with pytest.raises(NoWaterCoordinates) as excinfo:
        cluster_coordinates_only([1.0, 2.0, 3.0, 4.0])
    assert "3-dimensional" in str(excinfo.value)


def test_real_coordinates_still_cluster():
    """The guard must not have changed the working path."""
    import numpy as np

    from WatCon.find_conserved_networks import cluster_coordinates_only

    rng = np.random.default_rng(0)
    blob = np.concatenate([rng.normal(loc, 0.2, size=(40, 3))
                           for loc in (0.0, 10.0, 20.0)])
    labels, centers = cluster_coordinates_only(blob, min_samples=10)
    assert len(labels) == 120
    assert len(centers) == 3

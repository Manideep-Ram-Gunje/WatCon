"""Rigid-body superposition without MODELLER.

WatCon could only superpose structures through MODELLER, which is licensed and
often absent.  ``align_with_waters`` -- which already moves waters along with
their protein -- just needs rotation and translation matrices, so this supplies
them with numpy.

Two failures are worth more than the arithmetic here, and both have named tests:
a superposition onto a *mirror image* (an improper rotation the SVD will happily
return), and a correspondence keyed on the chain letter, which makes two copies
of the same protein look like they share no residues at all.
"""

from __future__ import annotations

import os
import tempfile

import numpy as np
import pytest

from WatCon.superpose import (
    SuperpositionError,
    ca_coordinates,
    kabsch,
    rmsd,
    superpose_structures,
)

from .conftest import bundle_member


# ===========================================================================
# Fixtures: the three barnase copies inside 1BRS
# ===========================================================================

@pytest.fixture(scope="module")
def barnase_copies():
    """Chains A, B and C of 1BRS, split into one file each.

    Three crystallographically independent copies of the same protein, in one
    coordinate frame but at different positions -- exactly what superposition is
    for, and real rather than synthetic.
    """
    raw = bundle_member(
        "1788241897_ConSurf.tar.gz", "1BRS_A_ATOMS_section_With_ConSurf.pdb"
    ).decode(errors="replace")

    directory = tempfile.mkdtemp()
    for chain in "ABC":
        path = os.path.join(directory, "copy_%s.pdb" % chain)
        with open(path, "w") as handle:
            for line in raw.splitlines():
                if (
                    line.startswith(("ATOM", "HETATM"))
                    and len(line) > 21
                    and line[21] == chain
                ):
                    handle.write(line + "\n")
            handle.write("END\n")
    return directory


# ===========================================================================
# kabsch
# ===========================================================================

def test_identity_gives_no_transform():
    points = np.random.default_rng(0).normal(size=(40, 3)) * 10
    rotation, translation = kabsch(points, points)
    assert np.allclose(rotation, np.eye(3))
    assert np.allclose(translation, 0.0, atol=1e-9)
    assert rmsd(points @ rotation.T + translation, points) < 1e-9


def test_a_known_transform_is_recovered_exactly():
    points = np.random.default_rng(1).normal(size=(50, 3)) * 8
    angle = 0.7
    expected_r = np.array(
        [
            [np.cos(angle), -np.sin(angle), 0.0],
            [np.sin(angle), np.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    expected_t = np.array([5.0, -3.0, 2.0])

    rotation, translation = kabsch(points, points @ expected_r.T + expected_t)
    assert np.allclose(rotation, expected_r)
    assert np.allclose(translation, expected_t)


def test_a_mirror_image_is_never_produced():
    """The reflection correction.

    Without it the SVD returns an improper rotation (det -1) for a mirrored
    target, superposing the protein onto its own mirror image -- not a rigid
    motion, and silently the wrong handedness.
    """
    points = np.random.default_rng(2).normal(size=(30, 3)) * 6
    mirrored = points.copy()
    mirrored[:, 2] *= -1

    rotation, _ = kabsch(points, mirrored)
    assert np.linalg.det(rotation) == pytest.approx(1.0)


def test_the_fit_is_optimal():
    """No other rigid transform may beat the one returned."""
    rng = np.random.default_rng(3)
    mobile = rng.normal(size=(40, 3)) * 5
    target = mobile @ np.eye(3) + np.array([1.0, 2.0, 3.0]) + rng.normal(
        size=(40, 3)
    ) * 0.2

    rotation, translation = kabsch(mobile, target)
    best = rmsd(mobile @ rotation.T + translation, target)

    for _ in range(20):
        perturbation = rng.normal(size=3) * 0.05
        q, _ = np.linalg.qr(np.eye(3) + rng.normal(size=(3, 3)) * 0.02)
        worse = rmsd(mobile @ (rotation @ q).T + translation + perturbation, target)
        assert best <= worse + 1e-9


def test_too_few_points_is_refused():
    with pytest.raises(SuperpositionError, match="at least 3"):
        kabsch(np.zeros((2, 3)), np.zeros((2, 3)))


def test_mismatched_shapes_are_refused():
    with pytest.raises(SuperpositionError, match="same shape"):
        kabsch(np.zeros((5, 3)), np.zeros((4, 3)))


# ===========================================================================
# Reading coordinates
# ===========================================================================

def test_ca_coordinates_reads_one_point_per_residue(barnase_copies):
    coords = ca_coordinates(os.path.join(barnase_copies, "copy_A.pdb"))
    assert len(coords) == 108
    assert all(v.shape == (3,) for v in coords.values())


def test_chain_filter_is_respected(barnase_copies):
    path = os.path.join(barnase_copies, "copy_B.pdb")
    assert ca_coordinates(path, chain="B")
    assert ca_coordinates(path, chain="A") == {}


# ===========================================================================
# Whole structures, on real data
# ===========================================================================

def test_copies_superpose_onto_the_reference(barnase_copies):
    result = superpose_structures(barnase_copies, reference="copy_A.pdb", chain=None)

    assert result.reference == "copy_A.pdb"
    assert len(result.rotations) == 2
    assert len(result.translations) == 2

    for name in ("copy_B.pdb", "copy_C.pdb"):
        value, count = result.quality[name]
        assert count == 108
        assert value < 1.0, "independent copies of one protein should fit closely"


def test_superposition_actually_moves_the_copies_together(barnase_copies):
    """The transforms must be applied, not merely computed."""
    result = superpose_structures(barnase_copies, reference="copy_A.pdb", chain=None)
    reference = {k[1:]: v for k, v in
                 ca_coordinates(os.path.join(barnase_copies, "copy_A.pdb")).items()}

    for i, name in enumerate(result.names[1:]):
        mobile = {k[1:]: v for k, v in
                  ca_coordinates(os.path.join(barnase_copies, name)).items()}
        shared = sorted(set(reference) & set(mobile))
        moving = np.array([mobile[k] for k in shared])
        fixed = np.array([reference[k] for k in shared])

        before = rmsd(moving, fixed)
        after = rmsd(moving @ result.rotations[i].T + result.translations[i], fixed)

        assert before > 10.0, "the copies start far apart"
        assert after < 1.0
        assert after < before / 10


def test_correspondence_ignores_the_chain_letter(barnase_copies):
    """Regression: the copies are chains A, B and C.

    Keying the correspondence on ``(chain, resid, icode)`` made them share zero
    residues and appear un-superposable.  A chain letter is a label the
    depositor chose, not part of a residue's identity across structures.
    """
    result = superpose_structures(barnase_copies, reference="copy_A.pdb", chain=None)

    # Each copy carries a different chain letter, yet all 108 residues correspond.
    chains = {
        key[0]
        for name in ("copy_A.pdb", "copy_B.pdb", "copy_C.pdb")
        for key in ca_coordinates(os.path.join(barnase_copies, name))
    }
    assert chains == {"A", "B", "C"}, "the copies really do differ by chain"

    for name in ("copy_B.pdb", "copy_C.pdb"):
        _, matched = result.quality[name]
        assert matched == 108


def test_ambiguous_chains_are_reported_not_guessed():
    """A file with two chains has two residues numbered 1; say so."""
    directory = tempfile.mkdtemp()
    lines = []
    for chain in "AB":
        for resid in range(1, 6):
            lines.append(
                "ATOM  %5d  CA  ALA %s%4d    %8.3f%8.3f%8.3f  1.00  0.00           C"
                % (len(lines) + 1, chain, resid, resid * 1.0, 0.0, 0.0)
            )
    for name in ("a.pdb", "b.pdb"):
        with open(os.path.join(directory, name), "w") as handle:
            handle.write("\n".join(lines) + "\nEND\n")

    with pytest.raises(SuperpositionError, match="more than one chain"):
        superpose_structures(directory, chain=None)


def test_unrelated_structures_are_refused(barnase_copies):
    """Too few shared residue numbers must fail loudly, not align on noise.

    A silently skipped structure would vanish from the water clustering without
    appearing anywhere in the results.
    """
    directory = tempfile.mkdtemp()
    import shutil

    shutil.copyfile(
        os.path.join(barnase_copies, "copy_A.pdb"), os.path.join(directory, "a.pdb")
    )
    with open(os.path.join(directory, "z.pdb"), "w") as handle:
        for i in range(5):
            handle.write(
                "ATOM  %5d  CA  ALA X%4d    %8.3f%8.3f%8.3f  1.00  0.00           C\n"
                % (i + 1, 9000 + i, i * 1.0, 0.0, 0.0)
            )
        handle.write("END\n")

    with pytest.raises(SuperpositionError, match="shares only"):
        superpose_structures(directory, reference="a.pdb", chain=None)


def test_max_rmsd_is_enforced(barnase_copies):
    with pytest.raises(SuperpositionError, match="above the"):
        superpose_structures(
            barnase_copies, reference="copy_A.pdb", chain=None, max_rmsd=0.01
        )


def test_a_missing_reference_is_refused(barnase_copies):
    with pytest.raises(SuperpositionError, match="not among the structures"):
        superpose_structures(barnase_copies, reference="nope.pdb")


# ===========================================================================
# Interoperating with WatCon's existing alignment path
# ===========================================================================

def test_as_dict_matches_the_shape_watcon_already_uses(barnase_copies):
    result = superpose_structures(barnase_copies, reference="copy_A.pdb", chain=None)
    transforms = result.as_dict()

    assert set(transforms) == {"Rot", "Trans"}
    assert len(transforms["Rot"]) == len(transforms["Trans"]) == 2
    assert transforms["Rot"][0].shape == (3, 3)
    assert transforms["Trans"][0].shape == (3,)


def test_a_reference_that_does_not_sort_first_is_refused(barnase_copies):
    """align_with_waters always treats the sorted-first file as the reference.

    A different reference would pair every structure with the wrong transform,
    and the result would still look like a superposition -- just a wrong one.
    """
    result = superpose_structures(barnase_copies, reference="copy_C.pdb", chain=None)
    # Computing it is fine; handing it to align_with_waters is not.
    assert len(result.rotations) == 2
    with pytest.raises(SuperpositionError, match="sorts first"):
        result.as_dict()


def test_describe_reports_rmsd_and_atom_counts(barnase_copies):
    text = superpose_structures(
        barnase_copies, reference="copy_A.pdb", chain=None
    ).describe()
    assert "reference: copy_A.pdb" in text
    assert "copy_B.pdb" in text
    assert "RMSD" in text
    assert "108 CA" in text

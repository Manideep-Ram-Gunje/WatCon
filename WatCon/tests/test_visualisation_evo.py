"""Projections of EVOLUTIONARY conservation.

WatCon already writes STRUCTURAL water conservation into the B-factor column
(``project_clusters``).  These outputs do the same for evolutionary
conservation, into separate files, so a user can load both and compare rather
than being handed one blended number.

The tests below therefore care about two things: that the numbers land in the
right column, and that the two kinds of conservation never overwrite each other.
"""

from __future__ import annotations

import os

import pytest

from WatCon.evolutionary import ResidueConservation, conservation_of_clusters
from WatCon.visualize_structures import (
    CONSURF_GRADE_COLORS,
    project_clusters,
    project_clusters_by_conservation,
    pymol_project_evolutionary,
)

from .test_evolutionary_clusters import (
    FakeAtom,
    FakeNetwork,
    FakeWater,
    conservation,
    wat_prot,
)


def sample_clusters():
    """Two sites: one scored (grade 9), one with no waters at all."""
    water = FakeWater(100, (1.0, 2.0, 3.0))
    net = FakeNetwork(
        [water], [FakeAtom(1, "A", 10, conservation(-1.25, 9))], [wat_prot(1, 100)]
    )
    centers = [(1.0, 2.0, 3.0), (50.0, 50.0, 50.0)]
    return net, centers, conservation_of_clusters([net], centers, 1.5)


def bfactors(path):
    values = []
    for line in open(path):
        if line.startswith("ATOM"):
            values.append(float(line[60:66]))
    return values


# ===========================================================================
# Cluster projection
# ===========================================================================

def test_grade_lands_in_the_bfactor_column(tmp_path):
    _, centers, clusters = sample_clusters()
    path = project_clusters_by_conservation(
        clusters, centers, filename_base="evo", out_dir=str(tmp_path), value="grade"
    )
    assert bfactors(path) == [9.0, 0.0]


def test_score_option_writes_the_continuous_measure(tmp_path):
    _, centers, clusters = sample_clusters()
    path = project_clusters_by_conservation(
        clusters, centers, filename_base="evo_s", out_dir=str(tmp_path), value="score"
    )
    values = bfactors(path)
    assert values[0] == pytest.approx(-1.25)
    assert values[1] == 0.0


def test_coordinates_are_preserved(tmp_path):
    _, centers, clusters = sample_clusters()
    path = project_clusters_by_conservation(
        clusters, centers, filename_base="evo", out_dir=str(tmp_path)
    )
    first = [line for line in open(path) if line.startswith("ATOM")][0]
    assert float(first[30:38]) == pytest.approx(1.0)
    assert float(first[38:46]) == pytest.approx(2.0)
    assert float(first[46:54]) == pytest.approx(3.0)


def test_no_data_is_documented_not_disguised(tmp_path):
    """B-factor 0 must be explained, so it is not read as 'not conserved'."""
    _, centers, clusters = sample_clusters()
    path = project_clusters_by_conservation(
        clusters, centers, filename_base="evo", out_dir=str(tmp_path)
    )
    text = open(path).read()
    assert "NO DATA" in text
    assert "sites with no conservation data: 1" in text


def test_sign_convention_is_stated_for_score(tmp_path):
    _, centers, clusters = sample_clusters()
    path = project_clusters_by_conservation(
        clusters, centers, filename_base="evo_s", out_dir=str(tmp_path), value="score"
    )
    assert "MORE NEGATIVE = MORE CONSERVED" in open(path).read()


def test_rejects_an_unknown_value(tmp_path):
    _, centers, clusters = sample_clusters()
    with pytest.raises(ValueError, match="grade.*score"):
        project_clusters_by_conservation(
            clusters, centers, out_dir=str(tmp_path), value="nonsense"
        )


def test_dict_centers_are_accepted(tmp_path):
    net, _, _ = sample_clusters()
    centers = {3: (1.0, 2.0, 3.0)}
    clusters = conservation_of_clusters([net], centers, 1.5)
    path = project_clusters_by_conservation(
        clusters, centers, filename_base="evo_d", out_dir=str(tmp_path)
    )
    assert bfactors(path) == [9.0]


# ===========================================================================
# Structural and evolutionary outputs must stay separate
# ===========================================================================

def test_structural_and_evolutionary_files_differ(tmp_path):
    """The two conservations go to different files and carry different numbers."""
    _, centers, clusters = sample_clusters()

    evo = project_clusters_by_conservation(
        clusters, centers, filename_base="separate", out_dir=str(tmp_path)
    )
    # project_clusters writes structural conservation, hardcoded to cluster_pdbs/
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        project_clusters(centers, filename_base="separate", b_factors=[0.25, 0.75])
        structural = os.path.join("cluster_pdbs", "separate.pdb")
        assert os.path.abspath(structural) != os.path.abspath(evo)
        assert bfactors(structural) == [0.25, 0.75]
    finally:
        os.chdir(cwd)

    assert bfactors(evo) == [9.0, 0.0]


# ===========================================================================
# Residue colouring
# ===========================================================================

def test_pml_uses_the_consurf_scale(tmp_path):
    net, _, _ = sample_clusters()
    path = pymol_project_evolutionary(
        net, filename="evo.pml", out_path=str(tmp_path)
    )
    text = open(path).read()
    assert "set_color consurf_9" in text
    assert "color consurf_9, resi 10 and chain A" in text


def test_pml_reports_unscored_residues_without_colouring_them(tmp_path):
    water = FakeWater(100, (0.0, 0.0, 0.0))
    net = FakeNetwork(
        [water],
        [FakeAtom(1, "A", 10, conservation(-1.0, 9)), FakeAtom(2, "A", 20, None)],
        [wat_prot(1, 100), wat_prot(2, 100)],
    )
    text = open(pymol_project_evolutionary(
        net, filename="u.pml", out_path=str(tmp_path))).read()

    assert "resi 10" in text
    assert "color consurf" not in text.split("resi 20")[0].split("\n")[-1]
    assert "no ConSurf score: 1" in text
    assert "residues coloured: 1" in text


def test_pml_grade_cutoff(tmp_path):
    water = FakeWater(100, (0.0, 0.0, 0.0))
    net = FakeNetwork(
        [water],
        [FakeAtom(1, "A", 10, conservation(-1.0, 9)),
         FakeAtom(2, "A", 20, conservation(0.9, 2))],
        [wat_prot(1, 100), wat_prot(2, 100)],
    )
    text = open(pymol_project_evolutionary(
        net, filename="c.pml", out_path=str(tmp_path), grade_cutoff=8)).read()

    assert "resi 10" in text
    assert "resi 20" not in text
    assert "residues coloured: 1" in text


def test_pml_counts_residues_not_atoms(tmp_path):
    """One residue reached through several atoms is coloured once."""
    shared = conservation(-1.0, 9)
    net = FakeNetwork(
        [FakeWater(100, (0.0, 0.0, 0.0))],
        [FakeAtom(1, "A", 10, shared), FakeAtom(2, "A", 10, shared),
         FakeAtom(3, "A", 10, shared)],
        [wat_prot(1, 100)],
    )
    text = open(pymol_project_evolutionary(
        net, filename="d.pml", out_path=str(tmp_path))).read()

    assert text.count("color consurf_9, resi 10") == 1
    assert "residues coloured: 1" in text


def test_consurf_scale_covers_all_grades():
    assert sorted(CONSURF_GRADE_COLORS) == list(range(1, 10))
    for rgb in CONSURF_GRADE_COLORS.values():
        assert len(rgb) == 3
        assert all(0.0 <= channel <= 1.0 for channel in rgb)

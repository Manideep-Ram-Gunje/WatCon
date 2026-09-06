"""The documented user workflow, run end to end on real data.

Every other test exercises a layer in isolation.  This one runs what a user
actually runs:

    structures on disk + ConSurf results on disk
        -> initialize_network(consurf_directory=...)
        -> cluster the water coordinates
        -> conservation_of_clusters
        -> report + projections

It exists because the layers can all be individually correct while the
assembled pipeline is broken -- wrong file-naming convention, a kwarg that never
reaches the builder, an output written to the wrong place.

Real structures: 1BRS (six chains, ConSurf run on chain A only) and 7O7W
(negative residue numbers, an omitted chromophore, fully ConSurf-covered).
"""

from __future__ import annotations

import csv
import os
import shutil

import pytest

pytest.importorskip("MDAnalysis", reason="end-to-end test needs MDAnalysis")

from .conftest import FIXTURES, bundle_member

PAIRS = [
    ("1788241897_ConSurf.tar.gz", "1BRS_A", "1BRS_A_150.grades.txt"),
    ("1788241750_ConSurf.tar.gz", "7O7W_A", "7O7W_A_45.grades.txt"),
]


@pytest.fixture(scope="module")
def workspace(tmp_path_factory):
    """A directory laid out the way the documentation tells a user to."""
    root = tmp_path_factory.mktemp("e2e")
    (root / "structures").mkdir()
    (root / "consurf").mkdir()

    for archive, base, fixture in PAIRS:
        (root / "structures" / (base + ".pdb")).write_bytes(
            bundle_member(archive, base + "_ATOMS_section_With_ConSurf.pdb")
        )
        # ConSurf files are matched to structures by the same naming convention
        # WatCon already uses for FASTA files.
        shutil.copyfile(
            FIXTURES / fixture,
            root / "consurf" / (base + "_consurf_grades.txt"),
        )
    return root


@pytest.fixture(scope="module")
def built(workspace):
    """Run initialize_network exactly as a user would."""
    from WatCon.generate_static_networks import initialize_network

    cwd = os.getcwd()
    os.chdir(workspace)
    try:
        metrics, networks, centers, names = initialize_network(
            "structures",
            network_type="water-protein",
            msa_indexing=False,
            classify_water=False,
            return_network=True,
            num_workers=1,
            max_distance=3.3,
            consurf_directory="consurf",
            consurf_strict=True,
        )
    finally:
        os.chdir(cwd)
    return metrics, networks, names


# ===========================================================================
# The pipeline runs
# ===========================================================================

def test_networks_are_built_for_every_structure(built):
    _, networks, names = built
    assert sorted(names) == ["1BRS_A", "7O7W_A"]
    assert len(networks) == 2
    assert all(net.protein_atoms for net in networks)


def test_conservation_reaches_the_atoms(built):
    """The kwarg must actually arrive at the builder, not be silently dropped."""
    _, networks, names = built
    o7w = networks[names.index("7O7W_A")]
    scored = [a for a in o7w.protein_atoms if a.evolutionary is not None]
    assert len(scored) == len(o7w.protein_atoms), "7O7W is fully ConSurf-covered"


def test_waters_are_annotated(built):
    _, networks, names = built
    o7w = networks[names.index("7O7W_A")]
    annotated = [w for w in o7w.water_molecules if w.evolutionary is not None]
    assert annotated
    assert all(w.evolutionary.n_residues >= 1 for w in annotated)


def test_coverage_metric_reports_the_chain_mismatch(built):
    """1BRS has six chains; ConSurf covered one.  Coverage must show it."""
    metrics, _, names = built
    brs = metrics[names.index("1BRS_A")]["evolutionary_coverage"]
    o7w = metrics[names.index("7O7W_A")]["evolutionary_coverage"]

    assert o7w["fraction"] == pytest.approx(1.0)
    assert brs["fraction"] < 0.5
    assert brs["not_in_file_atoms"] > 0     # the five uncovered chains
    assert brs["unscored_atoms"] == 0       # nothing inside chain A was missed


def test_strict_mode_rejects_a_missing_consurf_file(workspace, tmp_path):
    """A structure with no ConSurf file must fail loudly under strict."""
    from WatCon.evolutionary import ConservationError, load_conservation

    with pytest.raises(ConservationError, match="No ConSurf grades file"):
        load_conservation("NOSUCH.pdb", workspace / "consurf", strict=True)


def test_tolerant_mode_continues_without_conservation(workspace):
    from WatCon.evolutionary import load_conservation

    assert load_conservation("NOSUCH.pdb", workspace / "consurf", strict=False) is None


# ===========================================================================
# Clustering through to the report
# ===========================================================================

@pytest.fixture(scope="module")
def clustered(built):
    import numpy as np

    from WatCon.evolutionary import conservation_of_clusters
    from WatCon.find_conserved_networks import cluster_coordinates_only

    _, networks, names = built
    net = networks[names.index("7O7W_A")]
    coords = np.array([w.O.coordinates for w in net.water_molecules])
    _, centers = cluster_coordinates_only(
        coords, cluster="hdbscan", min_samples=2, eps=0.0
    )
    return net, centers, conservation_of_clusters([net], centers, dist_cutoff=1.5)


def test_clustering_produces_sites(clustered):
    _, centers, clusters = clustered
    assert len(centers) > 0
    assert len(clusters) == len(centers)


def test_occupancy_matches_the_geometry(clustered):
    """Sites reported empty must genuinely have no water within the cutoff.

    Cluster centres are the mean of their members, so a loose cluster's centroid
    can land in a gap.  That is correct, and this asserts the join agrees with
    the actual distances rather than dropping sites for some other reason.
    """
    import numpy as np

    net, centers, clusters = clustered
    coords = np.array([w.O.coordinates for w in net.water_molecules])

    for cluster_id, record in clusters.items():
        centre = np.array(centers[cluster_id])
        nearest = float(np.min(np.linalg.norm(coords - centre, axis=1)))
        if record.occupancy == 0:
            assert nearest > 1.5, (
                "site %d reported empty but a water is %.2f A away"
                % (cluster_id, nearest)
            )
        else:
            assert nearest <= 1.5


def test_report_is_written_and_well_formed(clustered, tmp_path):
    from WatCon.evolutionary import REPORT_COLUMNS, write_conservation_report

    _, _, clusters = clustered
    path = tmp_path / "conservation.csv"
    n = write_conservation_report(clusters, path)

    assert n == len(clusters)
    rows = list(csv.DictReader(open(path)))
    assert list(rows[0]) == REPORT_COLUMNS
    assert len(rows) == len(clusters)

    for row in rows:
        # Unoccupied sites must read NA, never 0 -- 0 would say "not conserved".
        if int(row["occupancy"]) == 0:
            assert row["evo_min_score"] == "NA"


def test_projection_bfactors_match_the_report(clustered, tmp_path):
    """The number in the file must be the number in the table."""
    from WatCon.visualize_structures import project_clusters_by_conservation

    _, centers, clusters = clustered
    path = project_clusters_by_conservation(
        clusters, centers, filename_base="e2e", out_dir=str(tmp_path)
    )

    written = [float(line[60:66]) for line in open(path) if line.startswith("ATOM")]
    expected = [
        float(clusters[cid].max_grade) if clusters[cid].has_conservation else 0.0
        for cid in sorted(clusters)
    ]
    assert written == expected


def test_pml_is_written(clustered, tmp_path):
    from WatCon.visualize_structures import pymol_project_evolutionary

    net, _, _ = clustered
    path = pymol_project_evolutionary(
        net, filename="e2e.pml", out_path=str(tmp_path)
    )
    text = open(path).read()
    assert "set_color consurf_" in text
    assert "residues coloured:" in text


# ===========================================================================
# Without ConSurf, nothing changes
# ===========================================================================

def test_pipeline_without_consurf_is_unaffected(workspace):
    """consurf_directory=None must leave the original behaviour intact."""
    from WatCon.generate_static_networks import initialize_network

    cwd = os.getcwd()
    os.chdir(workspace)
    try:
        metrics, networks, _, names = initialize_network(
            "structures",
            network_type="water-protein",
            msa_indexing=False,
            classify_water=False,
            return_network=True,
            num_workers=1,
            max_distance=3.3,
        )
    finally:
        os.chdir(cwd)

    assert len(networks) == 2
    for net in networks:
        assert net.protein_atoms
        assert all(a.evolutionary is None for a in net.protein_atoms)
        assert all(w.evolutionary is None for w in net.water_molecules)
    for m in metrics:
        assert "evolutionary_coverage" not in m

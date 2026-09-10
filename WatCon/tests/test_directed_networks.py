"""Directed water networks -- the path that had never executed.

``extract_objects`` passed ``max_neighbors=max_neigbhbors``, a name that was
never defined, so every call with ``directed=True`` raised ``NameError``. Static
directed networks therefore could not have run in any released version, and
nothing downstream of that typo has ever been exercised.

The typo is fixed. These tests keep the path executing, and cover the
hydrogen-bond geometry behind it -- including the branch whose cosine was being
normalised by an undefined vector.

Directed networks are built from explicit hydrogens (H -> O directionality), so
these use the water box with hydrogens rather than a crystal structure.
"""

from __future__ import annotations

import os
import warnings

import pytest

pytest.importorskip("MDAnalysis", reason="network building needs MDAnalysis")

from WatCon.generate_static_networks import initialize_network

WATER_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "water_dir")


def build(**kwargs):
    """Build a network from the bundled water box, directed unless told otherwise."""
    defaults = dict(
        network_type="water-water",
        msa_indexing=False,
        classify_water=False,
        return_network=True,
        num_workers=1,
        max_distance=3.3,
        water_name="HOH",
        include_hydrogens=True,      # this is what selects the directed path
    )
    defaults.update(kwargs)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _, networks, _, _ = initialize_network(WATER_DIR, **defaults)
    return networks[0]


# ===========================================================================
# The path runs at all
# ===========================================================================

def test_the_directed_path_executes():
    """Regression for the `max_neigbhbors` typo: this used to raise NameError."""
    network = build()
    assert network.graph.is_directed()


def test_undirected_is_still_undirected():
    network = build(include_hydrogens=False)
    assert not network.graph.is_directed()


def test_both_paths_find_the_same_waters():
    """Hydrogens change the edges, not which molecules exist."""
    assert (len(build().water_molecules)
            == len(build(include_hydrogens=False).water_molecules))


def test_directed_edges_are_real():
    network = build()
    assert network.graph.number_of_nodes() > 0
    assert network.graph.number_of_edges() > 0
    for u, v in network.graph.edges():
        assert u != v, "a water may not hydrogen-bond to itself"


# ===========================================================================
# Hydrogen-bond geometry
# ===========================================================================

def test_angle_criteria_can_only_remove_edges():
    """A geometric filter is a filter: it never invents a bond.

    This is the branch where the cosine was normalised by `water1`, a name
    undefined in that scope. Any edge surviving the filter must also exist
    without it.
    """
    unfiltered = build(angle_criteria=None)
    filtered = build(angle_criteria=120)

    assert filtered.graph.number_of_edges() <= unfiltered.graph.number_of_edges()
    assert set(filtered.graph.edges()) <= set(unfiltered.graph.edges())


def test_a_stricter_angle_is_never_more_permissive():
    """Monotonic in the threshold, and it must actually discriminate.

    The waters in this box form near-linear hydrogen bonds, so only strict
    thresholds cut anything -- 150 deg keeps all six, 179 deg keeps none. An
    earlier version of this test used 90-150, where nothing changes, and so
    passed while the filter was doing nothing at all.
    """
    counts = [
        build(angle_criteria=threshold).graph.number_of_edges()
        for threshold in (150, 160, 170, 179)
    ]
    assert counts == sorted(counts, reverse=True), counts
    assert counts[0] > counts[-1], "the filter must discriminate, not just not crash"


def test_an_impossible_angle_removes_everything():
    """Nothing can satisfy a 180-degree requirement in a real box."""
    assert build(angle_criteria=180).graph.number_of_edges() == 0


def test_angle_criteria_reaches_the_builder():
    """Regression: it was dropped entirely in the non-active-region path.

    `generate_directed_network` called `find_directed_connections` twice; the
    default branch hardcoded its arguments and never passed `angle_criteria`, so
    every threshold from 60 to 180 produced an identical graph.
    """
    assert (build(angle_criteria=179).graph.number_of_edges()
            < build(angle_criteria=None).graph.number_of_edges())


# ===========================================================================
# Parameters the hardcoded call site used to discard
# ===========================================================================

def test_max_distance_reaches_the_builder():
    """Regression: the cutoff was hardcoded to 2.5 A regardless of the request."""
    counts = [build(max_distance=d).graph.number_of_edges() for d in (2.0, 3.3, 8.0)]
    assert counts == sorted(counts), counts
    assert counts[0] < counts[-1], "a wider cutoff must admit more bonds"


def test_max_neighbors_reaches_the_builder():
    """Regression: hardcoded to 10, so the argument had no effect."""
    counts = [
        build(max_distance=8.0, max_neighbors=k).graph.number_of_edges()
        for k in (1, 2, 5, 10)
    ]
    assert counts == sorted(counts), counts
    assert counts[0] < counts[-1]


def test_a_single_neighbour_does_not_crash():
    """Regression: cKDTree.query(k=1) returns scalars, not length-1 rows.

    Every per-neighbour loop then raised "'numpy.float64' object is not
    iterable". It was unreachable while max_neighbors was hardcoded to 10.
    """
    build(max_distance=8.0, max_neighbors=1)
    build(max_distance=3.3, max_neighbors=1, include_hydrogens=False)


def test_metrics_are_produced_for_a_directed_graph():
    """Graph metrics must cope with directedness, not just undirected input."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        metrics, _, _, _ = initialize_network(
            WATER_DIR,
            network_type="water-water",
            msa_indexing=False,
            classify_water=False,
            return_network=True,
            num_workers=1,
            max_distance=3.3,
            water_name="HOH",
            include_hydrogens=True,
            analysis_conditions="all",
            analysis_selection="all",
        )
    assert metrics and isinstance(metrics[0], dict)
    density = metrics[0].get("density")
    assert density is None or density >= 0.0

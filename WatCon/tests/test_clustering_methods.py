"""Which clustering methods exist, and what happens when one does not.

`cluster_coordinates_only` dispatches on a string. Three methods are
implemented; a fourth, ``kmeans``, was advertised in the docstrings of *both*
`initialize_network` entry points and implemented nowhere. An unrecognised name
fell through every branch and then raised ``UnboundLocalError: cannot access
local variable 'clustering'`` -- which names neither the parameter nor the
accepted values, so a typo looked like an internal failure.

``kmeans`` was removed from the documentation rather than implemented: it needs
a *k* the API has no parameter for, and it cannot express the noise label (-1)
that every caller here relies on. Implementing it would have meant inventing a
specification, which is the same call made for
`residue_analysis.get_all_water_distances`.
"""

from __future__ import annotations

import numpy as np
import pytest

from WatCon.find_conserved_networks import (
    CLUSTERING_METHODS,
    cluster_coordinates_only,
)


@pytest.fixture(scope="module")
def coordinates():
    """Three well-separated blobs, so every method has something to find."""
    rng = np.random.RandomState(0)
    return np.vstack([rng.normal(centre, 0.3, size=(25, 3))
                      for centre in ([0, 0, 0], [10, 0, 0], [0, 10, 0])])


@pytest.mark.parametrize("method", CLUSTERING_METHODS)
def test_every_advertised_method_runs(method, coordinates):
    labels, centres = cluster_coordinates_only(coordinates, method, 5, None, 1)
    assert len(labels) == len(coordinates)
    assert isinstance(centres, dict)


@pytest.mark.parametrize("method", ["hdbscan", "optics"])
def test_the_density_methods_find_the_three_blobs(method, coordinates):
    """Not just "it ran" -- the answer has to be right."""
    _labels, centres = cluster_coordinates_only(coordinates, method, 5, None, 1)
    assert len(centres) == 3


def test_kmeans_is_refused_by_name():
    """It was advertised for a long time, so the message must be clear."""
    with pytest.raises(ValueError, match="Unknown clustering method 'kmeans'"):
        cluster_coordinates_only(np.zeros((10, 3)), "kmeans", 3, None, 1)


def test_a_typo_names_the_alternatives_instead_of_crashing():
    with pytest.raises(ValueError) as raised:
        cluster_coordinates_only(np.zeros((10, 3)), "hbdscan", 3, None, 1)
    message = str(raised.value)
    for method in CLUSTERING_METHODS:
        assert method in message


def test_kmeans_is_not_advertised_anywhere():
    """The docstrings promised it; that promise is what made it a defect."""
    import inspect

    from WatCon import generate_dynamic_networks, generate_static_networks

    for module in (generate_static_networks, generate_dynamic_networks):
        doc = inspect.getdoc(module.initialize_network) or ""
        assert "kmeans" not in doc, module.__name__


def test_eps_defaults_per_algorithm(coordinates):
    """None means 'this algorithm's own default', not a literal passed to sklearn."""
    _l, from_none = cluster_coordinates_only(coordinates, "hdbscan", 5, None, 1)
    _l, from_zero = cluster_coordinates_only(coordinates, "hdbscan", 5, 0.0, 1)
    assert len(from_none) == len(from_zero)

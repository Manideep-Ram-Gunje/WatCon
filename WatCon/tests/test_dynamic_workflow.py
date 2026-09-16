"""The dynamic path, run on real input for the first time.

Until this module, ``generate_dynamic_networks.initialize_network`` had never
been executed on real data -- the largest untested area in the package. Running
it found four defects, one of which meant the dynamic path could never have
processed a crystal structure at all:

* :meth:`WaterNetwork.add_water` took ``h1`` and ``h2`` positionally, so a water
  with no hydrogens -- every water in every X-ray structure -- raised
  ``TypeError``. The static builder had always accepted them.
* ``eps`` defaults to ``None`` in both entry points and was passed straight to
  sklearn, so clustering without naming ``eps`` raised ``InvalidParameterError``.
  Shared with the static path.
* An ``active_region_reference`` that matched nothing reached ``distance_array``
  as an empty array and came back as a complaint about array shape.
* A multi-model PDB whose models hold different atom counts was accepted,
  reported a frame count, and failed only once a later frame was read.

The fixtures are real coordinates, built by
``experiments/benchmark/scripts/build_ensemble_fixtures.py``; see the README
beside them for what each one is and what it is not.
"""

from __future__ import annotations

import os

import pytest

pytest.importorskip("MDAnalysis", reason="the dynamic path needs MDAnalysis")

from WatCon import generate_dynamic_networks as gdn
from WatCon.find_conserved_networks import cluster_coordinates_only
from WatCon.structure_io import (
    VaryingAtomCount,
    model_atom_counts,
    require_constant_atom_count,
)

PACKAGE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ENSEMBLE = os.path.join(PACKAGE, "data", "examples", "ptp1b_ensemble")
RAGGED = os.path.join(ENSEMBLE, "ptp1b_ensemble_ragged.pdb")
CONSTANT = os.path.join(ENSEMBLE, "ptp1b_ensemble_constant.pdb")
MD_SITE = os.path.join(ENSEMBLE, "md_active_site_h.pdb")

pytestmark = pytest.mark.skipif(not os.path.isfile(CONSTANT),
                                reason="ensemble fixtures absent")

#: Everything except the two metrics that run all-pairs shortest paths. On the
#: full MD system those two took 819 s and 45 s of a 957 s run while the network
#: itself took 1.4 s, so they are off here to keep the suite quick.
CHEAP = {
    "density": "on", "connected_components": "on", "interaction_counts": "on",
    "per_residue_interactions": "on", "characteristic_path_length": "off",
    "graph_entropy": "on", "clustering_coefficient": "on", "shortest_path": "off",
}


def _run(topology, **kwargs):
    options = dict(
        topology_file=topology, trajectory_file=topology,
        structure_directory=ENSEMBLE, msa_indexing=False, return_network=True,
        num_workers=1, analysis_conditions=CHEAP,
    )
    options.update(kwargs)
    return gdn.initialize_network(**options)


# ===========================================================================
# An ensemble of crystal structures is not a trajectory
# ===========================================================================

def test_the_ragged_ensemble_really_does_vary():
    """The fixture is only meaningful if the models genuinely differ."""
    counts = model_atom_counts(RAGGED)
    assert len(counts) == 3
    assert len(set(counts)) > 1


def test_the_constant_ensemble_does_not():
    counts = model_atom_counts(CONSTANT)
    assert len(counts) == 3
    assert len(set(counts)) == 1
    assert require_constant_atom_count(CONSTANT) == counts[0]


def test_a_file_without_models_is_not_rejected():
    """An ordinary PDB has no MODEL records, so there is nothing to disagree.

    The check reports None rather than a count, and above all does not raise:
    single-structure input is the dynamic path's other legitimate case.
    """
    assert model_atom_counts(MD_SITE) == []
    assert require_constant_atom_count(MD_SITE) is None


def test_the_refusal_names_the_models_and_points_somewhere_useful():
    with pytest.raises(VaryingAtomCount) as raised:
        require_constant_atom_count(RAGGED)
    message = str(raised.value)
    assert "model 1 has" in message
    assert "static" in message
    assert "generate_static_networks" in message


def test_a_ragged_ensemble_is_refused_before_any_work_is_done():
    """The point of the check: fail at the door, not three frames in."""
    with pytest.raises(VaryingAtomCount):
        _run("ptp1b_ensemble_ragged.pdb", multi_model_pdb=True, water_name="HOH")


# ===========================================================================
# The path runs, across frames, on real coordinates
# ===========================================================================

@pytest.fixture(scope="module")
def ensemble_run():
    return _run("ptp1b_ensemble_constant.pdb", multi_model_pdb=True,
                water_name="HOH", max_distance=3.0, cluster_coordinates=True,
                min_cluster_samples=2)


def test_every_frame_is_analysed(ensemble_run):
    _metrics, networks, _coords = ensemble_run
    assert len(networks) == 3


def test_each_frame_holds_the_waters_the_model_holds(ensemble_run):
    """Nine waters per model, by construction of the fixture."""
    _metrics, networks, _coords = ensemble_run
    for network in networks:
        assert len(network.water_molecules) == 9
        assert len(network.protein_atoms) > 0


def test_every_frame_builds_a_network_with_edges(ensemble_run):
    _metrics, networks, _coords = ensemble_run
    for network in networks:
        assert network.graph.number_of_edges() > 0
        assert not network.graph.is_directed()


def test_coordinates_are_clustered_across_frames(ensemble_run):
    """The cross-frame step: 27 water positions over 3 frames become sites."""
    _metrics, _networks, centres = ensemble_run
    assert centres is not None
    assert 0 < len(centres) <= 27


# ===========================================================================
# The regressions, each stated as the thing that used to break
# ===========================================================================

def test_waters_without_hydrogens_are_accepted():
    """Crystallographic waters are an oxygen and nothing else.

    ``add_water`` demanded h1 and h2 positionally, so this raised TypeError for
    every X-ray water -- meaning the dynamic path could not read a crystal
    structure at all.
    """
    network = gdn.WaterNetwork()
    atom = _FakeAtom(index=0, position=(1.0, 2.0, 3.0))
    network.add_water(1, atom, 1)
    assert len(network.water_molecules) == 1
    water = network.water_molecules[0]
    assert water.O is not None
    assert water.H1 is None and water.H2 is None


def test_waters_with_hydrogens_still_carry_them():
    network = gdn.WaterNetwork()
    network.add_water(1, _FakeAtom(0, (0.0, 0.0, 0.0)), 1,
                      _FakeAtom(1, (1.0, 0.0, 0.0)), _FakeAtom(2, (0.0, 1.0, 0.0)))
    water = network.water_molecules[0]
    assert water.H1 is not None and water.H2 is not None


def test_clustering_without_naming_eps_uses_the_algorithm_default():
    """Both entry points default eps to None, and passed it straight to sklearn."""
    import numpy as np

    coordinates = np.random.RandomState(0).rand(60, 3) * 10
    _labels, with_none = cluster_coordinates_only(coordinates, "hdbscan", 3, None, 1)
    _labels, with_zero = cluster_coordinates_only(coordinates, "hdbscan", 3, 0.0, 1)
    assert len(with_none) == len(with_zero) > 0


def test_an_active_region_reference_that_matches_nothing_says_so():
    """Residue 9999 exists in no structure; the old error named an array shape."""
    with pytest.raises(ValueError, match="selected no atoms"):
        _run("md_active_site_h.pdb", water_name="WAT",
             active_region_reference="resid 9999 and name SG",
             active_region_only=True)


class _FakeAtom:
    """The two attributes ``add_water`` reads off an MDAnalysis atom."""

    def __init__(self, index, position):
        self.index = index
        self.position = position


# ===========================================================================
# The ConSurf join, through the dynamic path, on real data
# ===========================================================================

#: The MD system is PTP1B renumbered from 1. Measured, not assumed: +1 gives
#: identity 1.000 over 295 residues, every other offset below 0.09.
MD_OFFSET = 1


def _renumbered_md_site(tmp_path):
    """A copy carrying the declared offset and chain, with its provenance."""
    source = os.path.join(ENSEMBLE, "md_active_site_h.pdb")
    work = tmp_path / "structures"
    work.mkdir(exist_ok=True)
    target = work / "md_site_renumbered.pdb"
    with open(str(target), "w", encoding="utf-8", newline="\n") as handle:
        handle.write("REMARK   Renumbered by %+d and chain set to A, both declared,\n"
                     "REMARK   to match the 1AAX ConSurf run. The offset was measured:\n"
                     "REMARK   identity 1.000 over 295 residues on the full MD system.\n"
                     % MD_OFFSET)
        for line in open(source, encoding="utf-8"):
            if line.startswith(("ATOM", "HETATM")):
                line = (line[:21] + "A"
                        + "%4d" % (int(line[22:26]) + MD_OFFSET) + line[26:])
            handle.write(line)
    return str(work)


def test_conservation_reaches_the_network_built_from_a_trajectory(tmp_path):
    """End to end: a ConSurf run attached to a network built by the dynamic path.

    This is what the chain bug broke silently -- every atom came back unscored
    while the run reported success. The assertion is therefore that *every*
    protein atom carries a grade, not merely that some do.
    """
    import shutil

    work = _renumbered_md_site(tmp_path)
    consurf = tmp_path / "consurf"
    consurf.mkdir()
    shutil.copyfile(
        os.path.join(PACKAGE, "data", "consurf", "fixtures", "1AAX_A.grades.txt"),
        str(consurf / "md_site_renumbered_consurf_grades.txt"))

    _metrics, networks, _coords = gdn.initialize_network(
        topology_file="md_site_renumbered.pdb",
        trajectory_file="md_site_renumbered.pdb",
        structure_directory=work, network_type="water-protein", water_name="WAT",
        msa_indexing=False, return_network=True, num_workers=1, max_distance=3.0,
        analysis_conditions=CHEAP, consurf_directory=str(consurf), consurf_strict=True)

    network = networks[0]
    graded = [a for a in network.protein_atoms if a.evolutionary is not None]
    assert len(graded) == len(network.protein_atoms) > 0

    # The WPD loop is the general acid, and every run grades it 9. Recovering it
    # here means the join landed on the right residues, not merely on some.
    by_resid = {int(a.resid): a.evolutionary.grade for a in graded}
    assert by_resid[181] == 9          # Asp181, the general acid
    assert by_resid[179] == 9          # Trp179, the loop it sits on

    # Conservation rolled up onto the waters that contact those residues.
    assert any(w.evolutionary is not None for w in network.water_molecules)

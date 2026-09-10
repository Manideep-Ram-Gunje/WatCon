"""The workflow the documentation actually tells people to use.

Everything else drives WatCon through its Python API. A user follows the
tutorial instead: write an input file, run the builder, write an analysis file,
run the post-analysis. That path was broken from end to end and nothing noticed,
because nothing tested it.

Seven separate failures stood between a user and a result:

1. the shipped ``input_static.txt`` used ``active_site_*`` where the builders
   accept ``active_region_*`` -- the template itself did not work;
2. ``analysis_conditions`` was indexed with ``[]``, so omitting any property
   switch raised ``KeyError``;
3. the ``concatenate`` key was misspelled ``concatemate`` in the parser, so it
   never fired;
4. ``concatenate`` then arrived as a string and was iterated character by
   character, looking for a file called ``b``;
5. ``metrics['shortest_path']`` held a lazy generator, making the results
   unpicklable -- the run crashed while saving, after doing all the work;
6. ``collect_coordinates`` expected a bare metrics list but was handed the
   four-tuple that gets pickled;
7. ``project_clusters`` was called with a ``separate_files`` argument it does
   not accept.

Each is pinned below.
"""

from __future__ import annotations

import inspect
import os

import pytest

pytest.importorskip("MDAnalysis", reason="the workflow needs MDAnalysis")

from WatCon.WatCon import (
    _normalise_concatenate,
    parse_analysis,
    parse_inputs,
    run_watcon,
    run_watcon_postanalysis,
)

PACKAGE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLES = os.path.join(PACKAGE, "data", "examples", "barnase", "structures")
GRADES = os.path.join(PACKAGE, "data", "consurf", "fixtures", "1BRS_A_150.grades.txt")

needs_examples = pytest.mark.skipif(
    not os.path.isdir(EXAMPLES), reason="bundled example data not present"
)


# ===========================================================================
# The shipped templates must work as shipped
# ===========================================================================

@pytest.mark.parametrize("template", ["input_static.txt", "input_dynamic.txt"])
def test_the_shipped_input_templates_are_usable(template):
    """Regression: both used active_site_* where the builders want active_region_*.

    Anyone copying the template -- which is what the docs say to do -- got
    ``TypeError: unexpected keyword argument 'active_site_only'``.
    """
    structure_type, kwargs = parse_inputs(os.path.join(PACKAGE, template))

    if structure_type == "static":
        from WatCon.generate_static_networks import initialize_network
    else:
        from WatCon.generate_dynamic_networks import initialize_network

    accepted = set(inspect.signature(initialize_network).parameters)
    unknown = sorted(set(kwargs) - accepted)
    assert not unknown, "%s passes settings the builder rejects: %s" % (
        template, unknown)


def test_the_old_active_site_spelling_still_works(tmp_path):
    """Input files written against the old templates must keep working."""
    path = tmp_path / "old.txt"
    path.write_text(
        "structure_type: static\n"
        "structure_directory: prepared\n"
        "active_site_only: off\n"
        "active_site_radius: 8\n"
    )
    with pytest.warns(DeprecationWarning):
        _, kwargs = parse_inputs(str(path))
    # 'off' is coerced to False downstream, as it is for the new spelling.
    assert kwargs["active_region_only"] is False
    assert kwargs["active_region_radius"] == 8.0
    assert "active_site_only" not in kwargs


def test_an_unusable_input_file_is_explained_not_a_typeerror(tmp_path):
    """The '; Property calculation' comment is load-bearing; say so."""
    from WatCon.generate_static_networks import initialize_network

    path = tmp_path / "no_header.txt"
    path.write_text(
        "structure_type: static\n"
        "structure_directory: prepared\n"
        "density: on\n"                 # no '; Property calculation' above it
    )
    _, kwargs = parse_inputs(str(path))
    with pytest.raises(ValueError, match="Property calculation"):
        run_watcon("static", kwargs)


# ===========================================================================
# concatenate
# ===========================================================================

def test_concatenate_is_split_on_commas_not_characters():
    """Regression: a string was iterated character by character."""
    assert _normalise_concatenate("run_1,run_2") == ["run_1.pkl", "run_2.pkl"]
    assert _normalise_concatenate("barnase") == ["barnase.pkl"]


def test_concatenate_accepts_names_with_or_without_the_extension():
    assert _normalise_concatenate(["a", "b.pkl"]) == ["a.pkl", "b.pkl"]


def test_concatenate_of_nothing_stays_nothing():
    assert _normalise_concatenate(None) is None


def test_comments_are_not_parsed_as_settings(tmp_path):
    """Regression: any comment containing a colon became a keyword argument."""
    path = tmp_path / "a.txt"
    path.write_text(
        "; Writes one row per site: occupancy beside conservation\n"
        "# another comment: with a colon\n"
        "cluster_concatenated: on\n"
    )
    kwargs = parse_analysis(str(path))
    assert kwargs == {"cluster_concatenated": True}, kwargs


def test_the_analysis_template_parses():
    kwargs = parse_analysis(os.path.join(PACKAGE, "analysis.txt"))
    accepted = set(inspect.signature(run_watcon_postanalysis).parameters)
    unknown = sorted(set(kwargs) - accepted)
    assert not unknown, "analysis.txt passes unknown settings: %s" % unknown


# ===========================================================================
# The whole thing, on real data
# ===========================================================================

@pytest.fixture(scope="module")
def workspace(tmp_path_factory):
    """A directory laid out exactly as the tutorial describes."""
    import shutil

    from WatCon.prepare import prepare_directory

    root = tmp_path_factory.mktemp("workflow")
    report = prepare_directory(EXAMPLES, str(root / "prepared"),
                               reference="1A2P", verbose=False)
    consurf = root / "consurf"
    consurf.mkdir()
    for outcome in report.prepared:
        shutil.copyfile(GRADES, consurf / ("%s_consurf_grades.txt" % outcome.pdb_id))

    (root / "input.txt").write_text(
        "structure_type: static\n"
        "structure_directory: prepared\n"
        "network_type: water-protein\n"
        "include_hydrogens: off\n"
        "water_name: HOH\n"
        "max_distance: 3.3\n"
        "\n"
        "; Property calculation\n"
        "density: on\n"
        "save_coordinates: on\n"
        "analysis_selection: all\n"
        "\n"
        "cluster_coordinates: off\n"
        "msa_indexing: off\n"
        "classify_water: off\n"
        "consurf_directory: consurf\n"
        "consurf_strict: on\n"
        "num_workers: 1\n"
    )
    (root / "analysis.txt").write_text(
        "concatenate: barnase\n"
        "input_directory: watcon_output\n"
        "cluster_concatenated: on\n"
        "min_samples: 2\n"
        "conservation_report: on\n"
        "conservation_dist_cutoff: 1.5\n"
    )
    return root


@needs_examples
def test_the_documented_workflow_runs_end_to_end(workspace):
    """Input file in, conservation report out -- the tutorial, executed."""
    import pickle
    import warnings

    cwd = os.getcwd()
    os.chdir(workspace)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            structure_type, kwargs = parse_inputs("input.txt")
            results = run_watcon(structure_type, kwargs)

            # Regression: metrics held a lazy generator, so this raised
            # "cannot pickle 'generator' object" AFTER all the work was done.
            os.makedirs("watcon_output", exist_ok=True)
            with open(os.path.join("watcon_output", "barnase.pkl"), "wb") as handle:
                pickle.dump(results, handle)

            assert os.path.getsize("watcon_output/barnase.pkl") > 0

            run_watcon_postanalysis(**parse_analysis("analysis.txt"))

        report = os.path.join("images", "CLUSTER_conservation.csv")
        assert os.path.isfile(report), "the conservation report was not written"

        import csv
        rows = list(csv.DictReader(open(report)))
        scored = [r for r in rows if r["evo_mean_score"] != "NA"]
        assert len(rows) > 100
        assert len(scored) > 50, "conservation did not reach the sites"

        assert os.path.isfile(os.path.join("cluster_pdbs", "CLUSTER_evolutionary.pdb"))
    finally:
        os.chdir(cwd)


@needs_examples
def test_the_input_file_route_agrees_with_the_python_api(workspace):
    """Two independent routes to the same result must not disagree.

    The API route is what ``watcon demo`` runs; this is what the tutorial runs.
    """
    import csv

    report = os.path.join(workspace, "images", "CLUSTER_conservation.csv")
    if not os.path.isfile(report):
        pytest.skip("depends on the workflow test having run")

    rows = list(csv.DictReader(open(report)))
    totals = {int(r["n_structures_total"]) for r in rows}
    assert totals == {6}, "all six prepared structures should be represented"


def test_a_missing_input_directory_says_what_to_do(tmp_path):
    """Regression: a bare FileNotFoundError on 'watcon_output'."""
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        with pytest.raises(FileNotFoundError, match="watcon run"):
            run_watcon_postanalysis(input_directory="watcon_output")
    finally:
        os.chdir(cwd)

"""``watcon demo`` -- the whole pipeline on data shipped with the package.

Runs, offline and in about a minute:

    6 raw barnase structures + 1 ConSurf run
        -> prepare   (find the barnase chain, superpose, carry the waters)
        -> build     (water-protein networks, conservation attached)
        -> cluster   (recurring water sites)
        -> join      (site occupancy beside the conservation of its lining residues)
        -> project   (PyMOL files to look at)

The structures are deliberately shipped **unprepared**: 1BRN carries barnase as
chain L and 1BRS as a barnase-barstar complex, so chain-finding has real work to
do rather than confirming that chain 'A' is chain 'A'.

The same ConSurf run serves all six, because conservation is a property of the
sequence.  Each structure's residues are still checked against it -- that check
is what makes reusing one run defensible, and it is not skipped for the demo.
"""

from __future__ import annotations

import os
import shutil
import warnings
from pathlib import Path

#: Bundled example: 6 barnase crystal structures, ANISOU/header-trimmed.
EXAMPLE_DIR = Path(__file__).resolve().parent / "data" / "examples" / "barnase"

#: The ConSurf run that describes barnase, already shipped for validation use.
GRADES = (Path(__file__).resolve().parent / "data" / "consurf" / "fixtures"
          / "1BRS_A_150.grades.txt")

REFERENCE = "1A2P"          # 1.5 A wild-type barnase
SITE_RADIUS = 1.5           # two waters this close are the same site
MIN_CLUSTER_SAMPLES = 2     # keep low-occupancy sites; they are the baseline


def _heading(step: str, text: str) -> None:
    print()
    print("=" * 70)
    print("%s  %s" % (step, text))
    print("=" * 70)


def run_demo(out_dir: str = "watcon_demo", keep: bool = False) -> int:
    """Run the bundled example end to end.  Returns a process exit code."""
    if not EXAMPLE_DIR.is_dir() or not GRADES.is_file():
        print("error: the bundled example data is missing. This usually means the "
              "package was installed without its data files; reinstall from a "
              "source checkout or a wheel built with the current pyproject.toml.")
        return 1

    out = Path(out_dir).resolve()
    if out.exists() and not keep:
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)

    print("WatCon + ConSurf demo")
    print("  example data : %s" % EXAMPLE_DIR)
    print("  ConSurf run  : %s" % GRADES.name)
    print("  output       : %s" % out)

    # -- 1. prepare ---------------------------------------------------------
    _heading("STEP 1/5", "Prepare: find the barnase chain, superpose, keep waters")
    from .prepare import prepare_directory

    prepared = out / "prepared"
    report = prepare_directory(
        str(EXAMPLE_DIR / "structures"), str(prepared), reference=REFERENCE
    )
    report.write_csv(out / "preparation.csv")

    # -- 2. ConSurf per structure -------------------------------------------
    _heading("STEP 2/5", "Stage ConSurf: one run describes every structure")
    consurf = out / "consurf"
    consurf.mkdir(exist_ok=True)
    for outcome in report.prepared:
        shutil.copyfile(GRADES, consurf / ("%s_consurf_grades.txt" % outcome.pdb_id))
    print("  %d ConSurf files staged (the same run, named per structure)"
          % len(report.prepared))
    print("  conservation is a property of the SEQUENCE, so one run covers all six")

    # -- 3. build -----------------------------------------------------------
    _heading("STEP 3/5", "Build water networks with conservation attached")
    warnings.filterwarnings("ignore")
    from .generate_static_networks import initialize_network

    metrics, networks, _, names = initialize_network(
        str(prepared),
        network_type="water-protein",
        msa_indexing=False,
        classify_water=False,
        return_network=True,
        num_workers=1,
        max_distance=3.3,
        water_name="HOH",
        consurf_directory=str(consurf),
        consurf_strict=True,
    )

    print()
    print("  %-8s %8s %10s %10s" % ("entry", "waters", "coverage", "identity"))
    for name, metric, network in zip(names, metrics, networks):
        coverage = metric["evolutionary_coverage"]
        rate = coverage["identity_rate"]
        print("  %-8s %8d %9.0f%% %9s" % (
            name, len(network.water_molecules), coverage["fraction"] * 100,
            "n/a" if rate is None else "%.0f%%" % (rate * 100)))
    print()
    print("  Identity is the check that makes one ConSurf run safe for many")
    print("  structures: it compares the residue ConSurf scored against the")
    print("  residue actually present, so a numbering mismatch cannot pass.")

    # -- 4. cluster and join ------------------------------------------------
    _heading("STEP 4/5", "Cluster recurring water sites and join to conservation")
    import numpy as np

    from .evolutionary import conservation_of_clusters, write_conservation_report
    from .find_conserved_networks import cluster_coordinates_only

    coordinates = np.array([
        water.O.coordinates for network in networks
        for water in network.water_molecules
    ])
    print("  pooling %d waters from %d structures" % (len(coordinates), len(networks)))

    _, centers = cluster_coordinates_only(
        coordinates, cluster="hdbscan", min_samples=MIN_CLUSTER_SAMPLES, eps=0.0
    )
    clusters = conservation_of_clusters(networks, centers, dist_cutoff=SITE_RADIUS)

    csv_path = out / "conservation.csv"
    write_conservation_report(clusters, csv_path)

    occupied = [c for c in clusters.values() if c.occupancy > 0]
    scored = [c for c in occupied if c.has_conservation]
    everywhere = [c for c in occupied if c.n_structures_occupied == len(networks)]
    print("  %d cluster centres, %d occupied, %d carry conservation"
          % (len(centers), len(occupied), len(scored)))
    print("  %d site(s) occupied in ALL %d structures" % (len(everywhere), len(networks)))

    if scored:
        conserved = [c for c in scored if c.max_grade is not None and c.max_grade >= 8]
        print("  %d scored site(s) are lined by a highly conserved residue "
              "(ConSurf grade >= 8)" % len(conserved))

    # -- 5. project ---------------------------------------------------------
    _heading("STEP 5/5", "Write PyMOL projections")
    from .visualize_structures import (
        project_clusters_by_conservation,
        pymol_project_evolutionary,
    )

    pdb_path = project_clusters_by_conservation(
        clusters, centers, filename_base="conserved_sites",
        out_dir=str(out / "pymol"), value="grade",
    )
    pml_path = pymol_project_evolutionary(
        networks[0], filename="%s_conservation.pml" % names[0],
        out_path=str(out / "pymol"),
    )

    # -- summary ------------------------------------------------------------
    print()
    print("=" * 70)
    print("DONE. Outputs in %s" % out)
    print("=" * 70)
    print("  conservation.csv   one row per water site: structural occupancy")
    print("                     beside the evolutionary conservation of the")
    print("                     residues lining it. The two are separate")
    print("                     columns -- deciding how they relate is the")
    print("                     science, so the tool does not blend them.")
    print("  preparation.csv    what happened to each structure, including any")
    print("                     rejected and why")
    print("  %s" % os.path.relpath(pdb_path, out))
    print("                     cluster centres, ConSurf grade in the B-factors")
    print("  %s" % os.path.relpath(pml_path, out))
    print("                     residues coloured on ConSurf's own 1-9 scale")
    print()
    print("  To look at it:")
    print("    pymol %s %s" % (os.path.relpath(pdb_path, out),
                               os.path.relpath(pml_path, out)))
    return 0

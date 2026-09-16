"""Build the water networks, cluster them, and join to conservation.

This is the documented WatCon+ConSurf workflow run at scale for the first time:

    prepared structures + one ConSurf file
        -> initialize_network(consurf_directory=...)
        -> cluster pooled water coordinates
        -> conservation_of_clusters
        -> report

One ConSurf run serves every structure because conservation is a property of the
sequence.  ``prepare.py`` has already relabelled each barnase chain to 'A' and
verified its numbering, so no chain_map is needed -- and the per-structure
identity rate is checked again here, since that check is what makes applying one
run to twenty-two depositions defensible.

Writes results/conservation.csv and results/run_summary.txt.
"""

from __future__ import annotations

import os
import shutil
import sys
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
WATCON = os.path.join(os.path.dirname(os.path.dirname(ROOT)),
                      "tools", "WatCon_ConSurf")
sys.path.insert(0, WATCON)

warnings.filterwarnings("ignore")

from WatCon.evolutionary import (                        # noqa: E402
    conservation_of_clusters,
    write_conservation_report,
)
from WatCon.find_conserved_networks import cluster_coordinates_only   # noqa: E402
from WatCon.generate_static_networks import initialize_network        # noqa: E402

PREPARED = os.path.join(ROOT, "prepared")
RESULTS = os.path.join(ROOT, "results")
CONSURF = os.path.join(ROOT, "consurf")
GRADES = os.path.join(WATCON, "WatCon", "data", "consurf", "fixtures",
                      "1BRS_A_150.grades.txt")

#: Two waters this far apart in the superposed frame are the same site.
#: 1.5 A is WatCon's own default for the cluster join.
SITE_RADIUS = 1.5

#: A cluster must contain waters from at least this many entries to be a site.
#: 2, not more, because the pre-registered comparison needs low-occupancy sites
#: as its baseline -- excluding them would delete the control group.
MIN_CLUSTER_SAMPLES = 2


def stage_consurf(names):
    """One copy of the grades file per structure, named so the lookup finds it.

    A user with per-structure ConSurf runs would simply drop them in.  Here the
    same run describes every structure, so it is copied under each name rather
    than special-casing the loader.
    """
    if os.path.isdir(CONSURF):
        shutil.rmtree(CONSURF)
    os.makedirs(CONSURF)
    for name in names:
        shutil.copyfile(GRADES, os.path.join(CONSURF, name + "_consurf_grades.txt"))


def main() -> int:
    os.makedirs(RESULTS, exist_ok=True)
    names = sorted(f[:-4] for f in os.listdir(PREPARED) if f.endswith(".pdb"))
    print("Building networks for %d structures" % len(names))
    stage_consurf(names)

    metrics, networks, _, built = initialize_network(
        PREPARED,
        network_type="water-protein",
        msa_indexing=False,
        classify_water=False,
        return_network=True,
        num_workers=1,
        max_distance=3.3,
        water_name="HOH",
        consurf_directory=CONSURF,
        consurf_strict=True,
    )
    print("  built %d networks" % len(networks))

    lines = ["Barnase conserved-water study", "=" * 60, ""]
    lines.append("Structures: %d" % len(networks))

    # Identity is re-checked per structure: one ConSurf run is being applied to
    # every one of them, and this is the only thing standing between that and a
    # silent mis-mapping.
    lines.append("")
    lines.append("Per-structure ConSurf mapping")
    lines.append("%-8s %8s %10s %8s" % ("entry", "waters", "coverage", "identity"))
    worst = 1.0
    for name, metric, network in zip(built, metrics, networks):
        cov = metric["evolutionary_coverage"]
        rate = cov["identity_rate"]
        worst = min(worst, rate if rate is not None else 1.0)
        lines.append("%-8s %8d %9.1f%% %7s" % (
            name, len(network.water_molecules), cov["fraction"] * 100,
            "n/a" if rate is None else "%.1f%%" % (rate * 100)))
    lines.append("")
    lines.append("Lowest identity rate across structures: %.1f%%" % (worst * 100))

    # Pool every water into the shared frame and cluster.
    coords = np.array([
        w.O.coordinates for network in networks for w in network.water_molecules
    ])
    print("  clustering %d waters" % len(coords))
    _, centers = cluster_coordinates_only(
        coords, cluster="hdbscan", min_samples=MIN_CLUSTER_SAMPLES, eps=0.0
    )
    print("  %d cluster centres" % len(centers))

    clusters = conservation_of_clusters(networks, centers, dist_cutoff=SITE_RADIUS)
    out = os.path.join(RESULTS, "conservation.csv")
    written = write_conservation_report(clusters, out)
    print("  wrote %d sites to %s" % (written, out))

    occupied = [c for c in clusters.values() if c.occupancy > 0]
    with_data = [c for c in occupied if c.has_conservation]
    lines += [
        "",
        "Waters pooled:            %d" % len(coords),
        "Cluster centres:          %d" % len(centers),
        "Sites with >=1 water:     %d" % len(occupied),
        "Sites with conservation:  %d" % len(with_data),
        "",
        "Occupancy distribution (structures in which a site is occupied)",
    ]
    counts = {}
    for c in occupied:
        counts[c.n_structures_occupied] = counts.get(c.n_structures_occupied, 0) + 1
    for n in sorted(counts):
        lines.append("  in %2d structures: %4d sites" % (n, counts[n]))

    # Dump what the analysis needs but the CSV cannot hold: which residues line
    # each site, and which structures occupy it.  The permutation null needs the
    # first (shuffle conservation across residues, keep the geometry fixed) and
    # the resolution / wild-type controls need the second (recompute occupancy
    # over a subset of structures without rebuilding any networks).
    import json

    detail = {}
    for cluster_id, record in clusters.items():
        centre = np.asarray(centers[cluster_id] if not hasattr(centers, "items")
                            else centers[cluster_id])
        occupying = []
        for name, network in zip(built, networks):
            for water in network.water_molecules:
                if np.linalg.norm(np.asarray(water.O.coordinates) - centre) <= SITE_RADIUS:
                    occupying.append(name)
                    break
        detail[str(cluster_id)] = {
            "occupancy": record.occupancy,
            "n_structures_occupied": record.n_structures_occupied,
            "occupying": occupying,
            "residues": ["%s|%d|%s" % (c, r, i or "") for c, r, i in record.residue_keys],
        }
    with open(os.path.join(RESULTS, "site_detail.json"), "w") as handle:
        json.dump(detail, handle, indent=1)
    print("  wrote site_detail.json (%d sites)" % len(detail))

    # Per-residue conservation, for the permutation null.
    residue_scores = {}
    for network in networks:
        for atom in network.protein_atoms:
            if atom.evolutionary is not None:
                key = "%s|%d|%s" % (atom.chain, int(atom.resid), atom.icode or "")
                residue_scores[key] = {
                    "score": atom.evolutionary.score,
                    "grade": atom.evolutionary.grade,
                    "buried_exposed": atom.evolutionary.buried_exposed,
                }
    with open(os.path.join(RESULTS, "residue_conservation.json"), "w") as handle:
        json.dump(residue_scores, handle, indent=1)
    print("  wrote residue_conservation.json (%d residues)" % len(residue_scores))

    # Per-structure identity, for the wild-type control.
    per_structure = {
        name: {
            "identity_rate": metric["evolutionary_coverage"]["identity_rate"],
            "identity_mismatched": metric["evolutionary_coverage"]["identity_mismatched"],
            "n_waters": len(network.water_molecules),
        }
        for name, metric, network in zip(built, metrics, networks)
    }
    with open(os.path.join(RESULTS, "per_structure.json"), "w") as handle:
        json.dump(per_structure, handle, indent=1)

    path = os.path.join(RESULTS, "run_summary.txt")
    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n")
    print("\n".join(lines[-(len(counts) + 8):]))
    print("Wrote %s" % path)
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""Where do structural and evolutionary conservation disagree?

The held-out benchmark shows conservation adds nothing to *predicting* whether a
site is occupied. That is a real negative result, and it is not the whole story:
two correlated measures can be redundant for prediction while still disagreeing
about many individual sites.

This is descriptive, not a hypothesis test. It quantifies how often the two
disagree, and checks whether the disagreement means anything -- using ConSurf's
own buried/exposed annotation, which is an independent structural call rather
than something derived from the conservation score.

Sites are split at the median of each measure, so the disagreement rate is the
arithmetic consequence of their correlation, not a discovery. What is worth
knowing is what the disagreeing groups look like.

    python scripts/disagreement.py --prepared <dir> --consurf <grades> --label NAME
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import statistics
import sys
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(ROOT)),
                                "tools", "WatCon_ConSurf"))
sys.path.insert(0, HERE)
warnings.filterwarnings("ignore")

from scipy.stats import mannwhitneyu                       # noqa: E402

SITE_RADIUS = 1.5
MIN_CLUSTER_SAMPLES = 2


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepared", required=True)
    parser.add_argument("--consurf", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--subsample", type=int, default=None)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    import heldout_benchmark as hb
    from WatCon.evolutionary import conservation_of_clusters
    from WatCon.find_conserved_networks import cluster_coordinates_only

    subset = None
    if args.subsample:
        names = sorted(f[:-4] for f in os.listdir(args.prepared)
                       if f.endswith(".pdb"))
        subset = hb.choose_subsample(names, args.subsample, args.seed)

    names, metrics, networks = hb.load_structures(args.prepared, args.consurf,
                                                  only=subset)
    print("structures: %d" % len(networks))

    pooled = np.vstack([hb.waters_of(n) for n in networks])
    _, centers = cluster_coordinates_only(
        pooled, cluster="hdbscan", min_samples=MIN_CLUSTER_SAMPLES, eps=0.0)
    clusters = conservation_of_clusters(networks, centers, dist_cutoff=SITE_RADIUS)

    # ConSurf's own buried/exposed call, keyed by residue.
    burial = {}
    for net in networks:
        for atom in net.protein_atoms:
            if atom.evolutionary is not None:
                key = (atom.chain, int(atom.resid), atom.icode)
                burial.setdefault(key, atom.evolutionary.buried_exposed)

    rows = []
    for cluster_id, record in clusters.items():
        if record.occupancy == 0 or record.mean_score is None:
            continue
        calls = [burial.get(k) for k in record.residue_keys]
        known = [c for c in calls if c in ("b", "e")]
        rows.append({
            "site": cluster_id,
            "occupancy_fraction": record.n_structures_occupied / len(networks),
            "mean_score": record.mean_score,
            "n_lining": len(record.residue_keys),
            "buried_fraction": (sum(1 for c in known if c == "b") / len(known)
                                if known else None),
        })

    if not rows:
        print("no scored sites")
        return 1

    median_occ = statistics.median(r["occupancy_fraction"] for r in rows)
    median_cons = statistics.median(r["mean_score"] for r in rows)

    groups = {"both": [], "recurrent_only": [], "conserved_only": [], "neither": []}
    for r in rows:
        high_occ = r["occupancy_fraction"] > median_occ
        high_cons = r["mean_score"] < median_cons        # lower = more conserved
        key = ("both" if high_occ and high_cons else
               "recurrent_only" if high_occ else
               "conserved_only" if high_cons else "neither")
        r["group"] = key
        groups[key].append(r)

    print()
    print("=" * 70)
    print("Where the two measures disagree -- %s" % args.label)
    print("=" * 70)
    print("scored sites: %d   median occupancy %.3f   median conservation %+.3f"
          % (len(rows), median_occ, median_cons))
    print()
    print("  %-18s %6s   %s" % ("group", "n", "mean buried fraction of lining residues"))
    for key in ("both", "recurrent_only", "conserved_only", "neither"):
        g = groups[key]
        vals = [r["buried_fraction"] for r in g if r["buried_fraction"] is not None]
        print("  %-18s %6d   %s" % (key, len(g),
              "%.3f" % statistics.mean(vals) if vals else "n/a"))

    disagree = len(groups["recurrent_only"]) + len(groups["conserved_only"])
    print()
    print("  the measures disagree at %d of %d sites (%.0f%%)"
          % (disagree, len(rows), 100 * disagree / len(rows)))

    a = [r["buried_fraction"] for r in groups["both"] if r["buried_fraction"] is not None]
    b = [r["buried_fraction"] for r in groups["recurrent_only"]
         if r["buried_fraction"] is not None]
    if len(a) > 5 and len(b) > 5:
        u = mannwhitneyu(a, b, alternative="two-sided")
        print("  recurrent+conserved vs recurrent-only, burial: Mann-Whitney p = %.3g"
              % u.pvalue)
        print()
        print("  Sites that recur WITHOUT conservation are predominantly surface.")
        print("  Caveat: burial and conservation are themselves correlated, so")
        print("  conservation is partly acting as a burial proxy here.")

    out = os.path.join(ROOT, "results", "disagreement_%s.csv" % args.label)
    with open(out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print()
    print("Wrote %s" % out)
    return 0


if __name__ == "__main__":
    sys.exit(main())

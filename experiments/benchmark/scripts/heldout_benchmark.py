"""Phase B: does conservation predict anything occupancy does not already?

Implements ../PREREGISTRATION.md exactly. Nothing here chooses a model, a
threshold or a direction -- all of that was fixed in writing and committed
before this file existed.

    Leave one structure out. Build water sites from the other N-1. Predict
    whether the held-out structure has a water at each site.

    M1 baseline (what WatCon gives you today) : occupancy fraction, n lining residues
    M2 ours                                    : M1 + mean ConSurf score

    Supported only if mean dAUC >= 0.02 AND paired Wilcoxon p < 0.05.

Usage:
    python scripts/heldout_benchmark.py --prepared <dir> --consurf <grades file> \
                                        --label barnase [--resolutions file.json]
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
WATCON = os.path.join(os.path.dirname(os.path.dirname(ROOT)),
                      "tools", "WatCon_ConSurf")
sys.path.insert(0, WATCON)
warnings.filterwarnings("ignore")

from scipy.stats import wilcoxon                                   # noqa: E402
from sklearn.linear_model import LogisticRegression                # noqa: E402
from sklearn.metrics import roc_auc_score                          # noqa: E402
from sklearn.preprocessing import StandardScaler                   # noqa: E402

# Pre-registered constants -- do not tune.
SITE_RADIUS = 1.5
MIN_CLUSTER_SAMPLES = 2
MIN_EFFECT = 0.02
ALPHA = 0.05
MIN_FOLDS = 10
MIN_SITES_PER_FOLD = 200
N_PERMUTATIONS = 200
SEED = 0


def load_structures(prepared, consurf, strict=True, only=None):
    """Build one network per prepared structure, with conservation attached.

    ``only`` restricts which structures are read at all. Building networks is the
    dominant cost, so loading 253 of them to use 50 wastes most of the run.
    """
    import shutil
    import tempfile

    from WatCon.generate_static_networks import initialize_network

    names = sorted(f[:-4] for f in os.listdir(prepared) if f.endswith(".pdb"))
    if only is not None:
        names = [n for n in names if n in only]

    staging = tempfile.mkdtemp()
    structures = tempfile.mkdtemp()
    for name in names:
        shutil.copyfile(os.path.join(prepared, name + ".pdb"),
                        os.path.join(structures, name + ".pdb"))
        shutil.copyfile(consurf, os.path.join(staging, name + "_consurf_grades.txt"))

    metrics, networks, _, built = initialize_network(
        structures, network_type="water-protein", msa_indexing=False,
        classify_water=False, return_network=True, num_workers=1,
        max_distance=3.3, water_name="HOH",
        consurf_directory=staging, consurf_strict=strict,
    )
    return built, metrics, networks


def waters_of(network):
    return np.array([w.O.coordinates for w in network.water_molecules], dtype=float)


def build_fold(train_networks, held_network):
    """Sites from the training structures; labels from the held-out one.

    Returns (occupancy_fraction, n_lining_residues, mean_conservation, label)
    for every site that has at least one lining residue carrying conservation.
    """
    from WatCon.evolutionary import conservation_of_clusters
    from WatCon.find_conserved_networks import cluster_coordinates_only

    pooled = np.vstack([waters_of(n) for n in train_networks])
    _, centers = cluster_coordinates_only(
        pooled, cluster="hdbscan", min_samples=MIN_CLUSTER_SAMPLES, eps=0.0
    )
    clusters = conservation_of_clusters(train_networks, centers,
                                        dist_cutoff=SITE_RADIUS)

    held = waters_of(held_network)
    ordered = (centers.items() if hasattr(centers, "items") else enumerate(centers))

    rows, residue_sets = [], []
    for cluster_id, centre in ordered:
        record = clusters.get(cluster_id)
        if record is None or not record.residue_keys:
            continue
        # Features come only from the training structures.
        occupancy = record.n_structures_occupied / len(train_networks)
        n_lining = len(record.residue_keys)
        if record.mean_score is None:
            continue
        # Label: does the held-out structure have a water here?
        distance = np.min(np.linalg.norm(held - np.asarray(centre), axis=1))
        rows.append((occupancy, n_lining, record.mean_score,
                     1 if distance <= SITE_RADIUS else 0))
        residue_sets.append(record.residue_keys)
    return np.array(rows, dtype=float), residue_sets


def auc_out_of_fold(train_x, train_y, test_x, test_y):
    """AUC on a fold, from a model fitted on the OTHER folds.

    Fitting and scoring on the same rows makes a third feature look useful
    whether or not it carries information -- more parameters always fit better
    in sample. A first version of this script did exactly that, reported
    dAUC = +0.0017 at p = 0.0008, and was only caught because the
    shuffled-conservation null centred on +0.0016 rather than zero. The extra
    feature was buying fit, not prediction.

    Training on the other folds makes the comparison a genuine prediction, and
    drives the null back to zero where it belongs.
    """
    if len(np.unique(test_y)) < 2 or len(np.unique(train_y)) < 2:
        return None
    scaler = StandardScaler().fit(train_x)
    model = LogisticRegression(max_iter=2000)
    model.fit(scaler.transform(train_x), train_y)
    scores = model.predict_proba(scaler.transform(test_x))[:, 1]
    return roc_auc_score(test_y, scores)


def choose_subsample(names, size, seed):
    """A reproducible random subset of structures.

    Leave-one-out over 253 PTP1B structures means 253 folds, each re-clustering
    ~60 000 waters -- hours of compute for folds that differ by one structure out
    of 253 and are therefore nearly identical to each other.

    The tempting shortcut is to cluster once and reuse the centres across folds.
    That is wrong: the held-out structure would help define the sites it is then
    scored against, leaking in exactly the direction that inflates the result.

    Subsampling keeps every fold honest and is stated as a deviation. The
    safeguard is running it twice with different seeds -- if the verdict is not
    stable across independent subsamples, it is not a result.
    """
    rng = np.random.default_rng(seed)
    chosen = rng.choice(len(names), size=min(size, len(names)), replace=False)
    return {names[i] for i in sorted(chosen)}


def run(label, prepared, consurf, subset=None):
    print("=" * 74)
    print("Held-out benchmark: %s" % label)
    print("=" * 74)

    names, metrics, networks = load_structures(prepared, consurf, only=subset)
    print("structures: %d" % len(names))

    # Pre-registered exclusion: identity below 95% means the conservation
    # mapping is untrustworthy. Applied BEFORE any modelling.
    keep = []
    for i, (name, metric) in enumerate(zip(names, metrics)):
        rate = metric.get("evolutionary_coverage", {}).get("identity_rate")
        if rate is not None and rate < 0.95:
            print("  excluded %s: identity %.1f%%" % (name, rate * 100))
        else:
            keep.append(i)
    names = [names[i] for i in keep]
    networks = [networks[i] for i in keep]
    print("structures after the identity filter: %d" % len(names))

    if len(names) < MIN_FOLDS:
        print("INCONCLUSIVE: %d folds, floor is %d" % (len(names), MIN_FOLDS))
        return None

    # Build every fold first, so each can be scored by a model trained on the
    # others -- see auc_out_of_fold for why in-sample scoring is not enough.
    built_folds = []
    for i, name in enumerate(names):
        train = [n for j, n in enumerate(networks) if j != i]
        data, residue_sets = build_fold(train, networks[i])
        if len(data) < MIN_SITES_PER_FOLD or len(np.unique(data[:, 3])) < 2:
            continue
        built_folds.append({"structure": name, "data": data,
                            "residues": residue_sets})

    per_fold = []
    for k, fold in enumerate(built_folds):
        others = np.vstack([f["data"] for j, f in enumerate(built_folds) if j != k])
        data = fold["data"]
        cols1, cols2 = [0, 1], [0, 1, 2]

        m1 = auc_out_of_fold(others[:, cols1], others[:, 3], data[:, cols1], data[:, 3])
        m2 = auc_out_of_fold(others[:, cols2], others[:, 3], data[:, cols2], data[:, 3])
        cons_only = auc_out_of_fold(others[:, [2]], others[:, 3],
                                    data[:, [2]], data[:, 3])
        if m1 is None or m2 is None:
            continue
        per_fold.append({"structure": fold["structure"], "n_sites": len(data),
                         "positives": int(data[:, 3].sum()),
                         "m1": m1, "m2": m2, "cons_only": cons_only,
                         "data": data, "residues": fold["residues"]})
        print("  %-10s sites=%4d pos=%4d  M1=%.4f  M2=%.4f  d=%+.4f"
              % (fold["structure"], len(data), int(data[:, 3].sum()), m1, m2, m2 - m1))

    if len(per_fold) < MIN_FOLDS:
        print("INCONCLUSIVE: %d usable folds, floor is %d" % (len(per_fold), MIN_FOLDS))
        return None

    m1 = np.array([f["m1"] for f in per_fold])
    m2 = np.array([f["m2"] for f in per_fold])
    delta = m2 - m1

    print()
    print("folds                 %d" % len(per_fold))
    print("M1 (occupancy + size) %.4f  (mean AUC)" % m1.mean())
    print("M2 (+ conservation)   %.4f" % m2.mean())
    print("conservation alone    %.4f"
          % np.mean([f["cons_only"] for f in per_fold if f["cons_only"] is not None]))
    print("mean dAUC             %+.4f" % delta.mean())

    stat = wilcoxon(m2, m1) if np.any(delta != 0) else None
    p = stat.pvalue if stat is not None else 1.0
    print("paired Wilcoxon p     %.3g" % p)

    supported = delta.mean() >= MIN_EFFECT and p < ALPHA
    if supported:
        verdict = "CONSERVATION ADDS MEASURABLE VALUE"
    elif p < ALPHA and 0 < delta.mean() < MIN_EFFECT:
        verdict = ("detectable but negligible (dAUC %.4f < %.2f)"
                   % (delta.mean(), MIN_EFFECT))
    elif delta.mean() <= 0:
        verdict = "conservation adds nothing"
    else:
        verdict = "not significant at p < %.2f" % ALPHA
    print("-> %s" % verdict)

    return {"label": label, "folds": per_fold, "m1": m1, "m2": m2,
            "delta": delta, "p": p, "verdict": verdict}


def permutation_null(result, rng_seed=SEED):
    """Shuffle conservation across residues; the null must centre on zero."""
    print()
    print("CONTROL -- shuffled-conservation null (%d permutations)" % N_PERMUTATIONS)
    rng = np.random.default_rng(rng_seed)
    folds = result["folds"]

    # One conservation value per distinct residue, pooled across folds.
    residues = {}
    for f in folds:
        for keys, value in zip(f["residues"], f["data"][:, 2]):
            for k in keys:
                residues.setdefault(k, value)
    keys = list(residues)
    values = np.array([residues[k] for k in keys])

    null = []
    for _ in range(N_PERMUTATIONS):
        permuted = dict(zip(keys, rng.permutation(values)))
        shuffled = {}
        for j, f in enumerate(folds):
            shuffled[j] = np.column_stack([
                f["data"][:, 0], f["data"][:, 1],
                np.array([np.mean([permuted[k] for k in ks]) for ks in f["residues"]]),
            ])
        deltas = []
        for j, f in enumerate(folds):
            others_x = np.vstack([shuffled[i] for i in shuffled if i != j])
            others_y = np.concatenate([folds[i]["data"][:, 3]
                                       for i in shuffled if i != j])
            m2 = auc_out_of_fold(others_x, others_y, shuffled[j], f["data"][:, 3])
            if m2 is not None:
                deltas.append(m2 - f["m1"])
        if deltas:
            null.append(np.mean(deltas))
    null = np.array(null)
    observed = result["delta"].mean()
    tail = float(np.mean(np.abs(null) >= abs(observed)))
    print("  null dAUC mean %+.4f, sd %.4f  (an unbiased null is ~0)"
          % (null.mean(), null.std()))
    print("  observed %+.4f, empirical p = %.3f" % (observed, tail))
    print("  -> observed effect is %s the null"
          % ("OUTSIDE" if tail < 0.05 else "INSIDE"))
    return {"null_mean": float(null.mean()), "null_sd": float(null.std()),
            "empirical_p": tail}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepared", required=True)
    parser.add_argument("--consurf", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--resolutions", default=None)
    parser.add_argument("--max-resolution", type=float, default=None)
    parser.add_argument("--subsample", type=int, default=None,
                        help="use a random subset of this many structures")
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()

    subset = None
    if args.resolutions and args.max_resolution:
        table = json.load(open(args.resolutions))
        subset = {k for k, v in table.items()
                  if v is not None and v <= args.max_resolution}

    if args.subsample:
        available = sorted(f[:-4] for f in os.listdir(args.prepared)
                           if f.endswith(".pdb"))
        if subset:
            available = [n for n in available if n in subset]
        subset = choose_subsample(available, args.subsample, args.seed)
        print("subsampled %d of %d structures (seed %d)"
              % (len(subset), len(available), args.seed))

    result = run(args.label, args.prepared, args.consurf, subset)
    if result is None:
        return 1
    control = permutation_null(result)

    out = os.path.join(ROOT, "results", "benchmark_%s.json" % args.label)
    with open(out, "w") as handle:
        json.dump({
            "label": args.label,
            "n_folds": len(result["folds"]),
            "m1_mean_auc": float(result["m1"].mean()),
            "m2_mean_auc": float(result["m2"].mean()),
            "mean_delta_auc": float(result["delta"].mean()),
            "wilcoxon_p": float(result["p"]),
            "verdict": result["verdict"],
            "null": control,
            "per_fold": [{k: v for k, v in f.items()
                          if k not in ("data", "residues")} for f in result["folds"]],
        }, handle, indent=1)
    print()
    print("Wrote %s" % out)
    return 0


if __name__ == "__main__":
    sys.exit(main())

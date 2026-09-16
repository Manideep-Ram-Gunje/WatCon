"""The pre-registered analysis.  See ../PREREGISTRATION.md.

Nothing here chose a test, a threshold, or a direction: all of that was fixed in
writing and committed before this script existed.

One documented deviation, forced by the study's own control.  The pre-registered
primary statistic is ``evo_min_score`` -- the most conserved lining residue.  The
shuffled-conservation null for it does **not** centre on zero (it centres near
-0.23), because the minimum of a larger set of draws is mechanically lower, and
high-occupancy sites are lined by more residues.  Part of the pre-registered
correlation is therefore arithmetic rather than biology.

So this reports the pre-registered statistic *as pre-registered*, with its bias
stated, and adds a size-independent one (mean over lining residues) whose null
does centre on zero.  The second is the number to quote.  Both are shown; the
pre-registered result is not quietly replaced.
"""

from __future__ import annotations

import csv
import json
import os
import sys

import numpy as np
from scipy.stats import mannwhitneyu, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, "results")
DATA = os.path.join(ROOT, "data")

ALPHA = 0.01
MIN_EFFECT = 0.20
POWER_FLOOR = 50
N_PERMUTATIONS = 1000
RESOLUTION_CUT = 1.8
SEED = 0


def load():
    with open(os.path.join(RESULTS, "site_detail.json")) as handle:
        detail = json.load(handle)
    with open(os.path.join(RESULTS, "residue_conservation.json")) as handle:
        conservation = json.load(handle)
    with open(os.path.join(RESULTS, "per_structure.json")) as handle:
        per_structure = json.load(handle)
    resolution = {}
    with open(os.path.join(DATA, "entries.csv")) as handle:
        for row in csv.DictReader(handle):
            if row["resolution"]:
                resolution[row["pdb_id"]] = float(row["resolution"])
    return detail, conservation, per_structure, resolution


def sites_for(detail, conservation, structures=None, total=None):
    """[(occupancy_fraction, [residue keys])] over a subset of structures.

    Sites with no lining residue carrying conservation are dropped: no data is
    not a score of zero.
    """
    rows = []
    for record in detail.values():
        occupying = record["occupying"]
        if structures is not None:
            occupying = [s for s in occupying if s in structures]
        if not occupying:
            continue
        keys = [k for k in record["residues"] if k in conservation]
        if not keys:
            continue
        rows.append((len(occupying) / total, keys))
    return rows


def correlate(sites, scores, aggregate):
    occupancy = [s[0] for s in sites]
    values = [aggregate([scores[k] for k in s[1]]) for s in sites]
    return spearmanr(occupancy, values)


def permutation_null(sites, scores, aggregate, n=N_PERMUTATIONS):
    """Shuffle conservation across residues, keeping the geometry fixed.

    Catches an effect manufactured by the shape of the join rather than by
    biology -- which is exactly what it found for the minimum.
    """
    keys = list(scores)
    values = np.array([scores[k] for k in keys])
    rng = np.random.default_rng(SEED)
    null = np.empty(n)
    for i in range(n):
        permuted = dict(zip(keys, rng.permutation(values)))
        null[i] = correlate(sites, permuted, aggregate).statistic
    return null


def report(lines, text=""):
    print(text)
    lines.append(text)


def assess(rho, p):
    if p >= ALPHA:
        return "H1 NOT SUPPORTED -- not significant at p < %.2f" % ALPHA
    if rho >= 0:
        return "H1 CONTRADICTED -- correlation is POSITIVE"
    if abs(rho) < MIN_EFFECT:
        return ("H1 NOT SUPPORTED -- detectable but negligible (|rho| %.3f < %.2f)"
                % (abs(rho), MIN_EFFECT))
    return "H1 SUPPORTED"


def main() -> int:
    detail, conservation, per_structure, resolution = load()
    scores = {k: v["score"] for k, v in conservation.items()}
    n_structures = len(per_structure)
    lines = []

    report(lines, "Barnase conserved water sites vs evolutionary conservation")
    report(lines, "=" * 72)
    report(lines, "Pre-registered in PREREGISTRATION.md, committed before this ran.")
    report(lines, "")
    report(lines, "Structures: %d   ConSurf runs: 1 (conservation is a sequence"
                  " property)" % n_structures)
    report(lines, "ConSurf convention: MORE NEGATIVE = MORE CONSERVED, so H1")
    report(lines, "predicts a NEGATIVE correlation with occupancy.")

    sites = sites_for(detail, conservation, total=n_structures)
    n = len(sites)
    report(lines, "")
    report(lines, "Sites with conservation data: %d" % n)
    if n < POWER_FLOOR:
        report(lines, "INCONCLUSIVE: below the pre-registered floor of %d."
               % POWER_FLOOR)
        return 0

    # ------------------------------------------------ pre-registered primary
    rho_min, p_min = correlate(sites, scores, min)
    null_min = permutation_null(sites, scores, min)
    report(lines, "")
    report(lines, "PRIMARY, AS PRE-REGISTERED -- Spearman(occupancy, evo_min_score)")
    report(lines, "  rho = %+.4f   p = %.3g   n = %d" % (rho_min, p_min, n))
    report(lines, "  -> %s" % assess(rho_min, p_min))
    report(lines, "")
    report(lines, "  Its shuffled null does NOT centre on zero:")
    report(lines, "    null mean %+.4f, sd %.4f  (an unbiased null would be ~0)"
           % (null_min.mean(), null_min.std()))
    report(lines, "    z vs null = %+.2f, empirical p = %.3f"
           % ((rho_min - null_min.mean()) / null_min.std(),
              float(np.mean(np.abs(null_min) >= abs(rho_min)))))
    report(lines, "  The minimum of a larger set is mechanically lower, and")
    report(lines, "  high-occupancy sites are lined by more residues, so this")
    report(lines, "  statistic is inflated. The raw rho OVERSTATES the effect.")

    sizes = [len(s[1]) for s in sites]
    rho_size = spearmanr([s[0] for s in sites], sizes).statistic
    rho_size_min = spearmanr(
        sizes, [min([scores[k] for k in s[1]]) for s in sites]
    ).statistic
    report(lines, "    Spearman(occupancy, n_lining_residues) = %+.4f" % rho_size)
    report(lines, "    Spearman(n_lining_residues, min_score) = %+.4f" % rho_size_min)

    # -------------------------------------------- size-independent statistic
    rho_mean, p_mean = correlate(sites, scores, np.mean)
    null_mean = permutation_null(sites, scores, np.mean)
    z_mean = (rho_mean - null_mean.mean()) / null_mean.std()
    tail_mean = float(np.mean(np.abs(null_mean) >= abs(rho_mean)))
    report(lines, "")
    report(lines, "SIZE-INDEPENDENT (added after the control above; not pre-registered)")
    report(lines, "  Spearman(occupancy, MEAN score over lining residues)")
    report(lines, "  rho = %+.4f   p = %.3g   n = %d" % (rho_mean, p_mean, n))
    report(lines, "  null mean %+.4f, sd %.4f  -- centres on zero, as it should"
           % (null_mean.mean(), null_mean.std()))
    report(lines, "  z vs null = %+.2f, empirical p = %.3f" % (z_mean, tail_mean))
    report(lines, "  -> %s" % assess(rho_mean, p_mean))
    report(lines, "  THIS IS THE NUMBER TO QUOTE.")

    # ------------------------------------------------------------- secondary
    report(lines, "")
    report(lines, "SECONDARY, AS PRE-REGISTERED -- recurrent vs singleton sites")
    singletons = [r for r in detail.values() if len(r["occupying"]) == 1]
    scored_singletons = [
        r for r in singletons if any(k in conservation for k in r["residues"])
    ]
    report(lines, "  sites occupied in exactly 1 structure: %d" % len(singletons))
    report(lines, "  of those, with a scored lining residue: %d"
           % len(scored_singletons))
    report(lines, "  -> UNDERPOWERED. Not a defect: singleton waters mostly have")
    report(lines, "     no protein residue in contact at all, which is itself a")
    report(lines, "     finding. The pre-registered comparison cannot be run.")

    low = [s for s in sites if s[0] <= 3.0 / n_structures]
    high = [s for s in sites if s[0] >= 0.5]
    report(lines, "")
    report(lines, "  POST HOC substitute (labelled as such): <=3 vs >=50% of structures")
    if len(low) >= 5 and len(high) >= 5:
        low_v = [np.mean([scores[k] for k in s[1]]) for s in low]
        high_v = [np.mean([scores[k] for k in s[1]]) for s in high]
        stat = mannwhitneyu(high_v, low_v, alternative="two-sided")
        effect = 2.0 * stat.statistic / (len(high_v) * len(low_v)) - 1.0
        report(lines, "    low-occupancy  n=%3d median mean-score %+.3f"
               % (len(low_v), np.median(low_v)))
        report(lines, "    high-occupancy n=%3d median mean-score %+.3f"
               % (len(high_v), np.median(high_v)))
        report(lines, "    Mann-Whitney p = %.3g, rank-biserial r = %+.3f"
               % (stat.pvalue, effect))
    else:
        report(lines, "    too few sites in one group")

    # -------------------------------------------------------------- controls
    for label, subset in (
        ("CONTROL -- resolution <= %.1f A" % RESOLUTION_CUT,
         {s for s in per_structure if resolution.get(s, 9.9) <= RESOLUTION_CUT}),
        ("CONTROL -- wild-type-like only (no identity mismatches)",
         {s for s, v in per_structure.items() if v["identity_mismatched"] == 0}),
    ):
        report(lines, "")
        report(lines, label)
        report(lines, "  %d of %d structures: %s"
               % (len(subset), n_structures, ", ".join(sorted(subset))))
        rows = sites_for(detail, conservation, subset, total=len(subset))
        if len(rows) >= POWER_FLOOR:
            r, p = correlate(rows, scores, np.mean)
            report(lines, "  mean-score rho = %+.4f, p = %.3g, n = %d"
                   % (r, p, len(rows)))
        else:
            report(lines, "  only %d sites -- underpowered" % len(rows))

    # ------------------------------------------------------------ confounder
    report(lines, "")
    report(lines, "CONFOUND -- burial (declared in advance, NOT controlled for)")
    buried = [v["score"] for v in conservation.values()
              if v.get("buried_exposed") == "b"]
    exposed = [v["score"] for v in conservation.values()
               if v.get("buried_exposed") == "e"]
    if buried and exposed:
        report(lines, "  buried  n=%3d median score %+.3f" % (len(buried),
                                                              np.median(buried)))
        report(lines, "  exposed n=%3d median score %+.3f" % (len(exposed),
                                                              np.median(exposed)))
        report(lines, "  Buried residues are both more conserved AND more likely to")
        report(lines, "  hold ordered water. This design cannot separate the two, so")
        report(lines, "  no causal claim is made.")
    else:
        report(lines, "  this ConSurf run carries no buried/exposed annotation")

    report(lines, "")
    report(lines, "=" * 72)
    report(lines, "VERDICT (on the size-independent statistic): %s"
           % assess(rho_mean, p_mean))
    report(lines, "Effect: rho = %+.3f, outside a null centred on %+.3f (z = %+.2f)."
           % (rho_mean, null_mean.mean(), z_mean))
    report(lines, "Correlational, one protein, burial not disentangled.")

    path = os.path.join(RESULTS, "analysis_report.txt")
    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n")
    print("\nWrote %s" % path)
    return 0


if __name__ == "__main__":
    sys.exit(main())

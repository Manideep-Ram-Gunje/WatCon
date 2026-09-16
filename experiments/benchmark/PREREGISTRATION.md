# Pre-registration — does evolutionary conservation add predictive value?

**Written before the analysis script existed.** `scripts/heldout_benchmark.py` is
committed *after* this file; the git history is the evidence. The pKa experiment
(`../pka/`) and the barnase study (`../barnase_waters/`) were pre-registered the
same way, and the pKa one returned a negative result that was reported as such.

---

## The question

WatCon finds water sites that recur across structures. We added ConSurf
evolutionary conservation. The honest question is not "is our tool better" —
the water-network algorithm is unchanged and there is no shared metric — but:

> **Does conservation predict anything that occupancy does not already?**

If it does not, the tool still *reports* something WatCon cannot, but it adds no
predictive power, and we say so.

## Design

Leave-one-structure-out cross-validation. For each held-out structure *i*:

1. Build water sites by clustering the pooled waters of the **other N−1**
   structures (HDBSCAN, `min_samples=2`, as in the barnase study).
2. Compute every feature from those N−1 structures only.
3. **Label**: does structure *i* have a water within **1.5 Å** of the site centre?

The held-out structure contributes nothing to the sites or the features, so the
label is a genuine prediction.

## Models

Logistic regression — interpretable, and the point is the difference, not the
absolute ceiling.

| model | features |
|---|---|
| **M1 — baseline (what WatCon gives you today)** | training occupancy fraction; **number of lining residues** |
| **M2 — ours** | M1 **+** mean ConSurf score over the lining residues |

Reported alongside, for interpretation only: conservation alone.

**Why `n_lining_residues` is in the baseline.** The barnase study measured
Spearman(occupancy, lining-residue count) = **+0.555**, and that dependence
inflated its pre-registered statistic. Including the count in M1 means any gain
credited to conservation cannot be site size wearing a disguise.

**Mean, never minimum.** The minimum of a larger set is mechanically lower; that
is the same artefact.

## Primary statistic

ΔAUC = AUC(M2) − AUC(M1), computed **per fold** and paired across folds.

## Success criteria — fixed in advance

Conservation **adds measurable value** only if **both** hold:

| | threshold |
|---|---|
| effect | mean ΔAUC ≥ **0.02** |
| significance | paired Wilcoxon signed-rank **p < 0.05** |

- Significant but ΔAUC < 0.02 → **"detectable but negligible"**, not support.
- ΔAUC ≤ 0 → **conservation adds nothing**, reported plainly.

**Power floor.** Fewer than **10 folds** or fewer than **200 labelled sites per
fold** → reported **inconclusive**, not negative.

## Controls — all run regardless of outcome

1. **Shuffled-conservation null.** Permute conservation across residues, keep the
   geometry, recompute ΔAUC (200 permutations). **The null must centre on zero.**
   If it does not, the statistic is measuring geometry and must be replaced —
   exactly what this control caught in the barnase study.
2. **Two independent proteins.** PTP1B and barnase, analysed separately. A result
   on one protein only is reported as such.
3. **Resolution restriction.** Repeat at ≤1.8 Å.

## Datasets

| | |
|---|---|
| **PTP1B** | X-ray, UniProt P18031, ≤2.0 Å, prepared with `watcon prepare`; conservation from **one** ConSurf run on 1AAX chain A |
| **Barnase** | the 22 structures from `../barnase_waters/` |

One ConSurf run per protein suffices because conservation is a property of the
sequence.

## Confounds, declared in advance

- **Occupancy is already strongly predictive of occupancy.** The bar for
  conservation is incremental, and it is a high bar. This is deliberate.
- **Sites are not independent** — neighbouring sites share lining residues.
  Folds are by *structure*, not by site, which limits but does not remove this.
- **Burial** remains unresolved from the barnase study: buried residues are both
  more conserved and more likely to hold ordered water. No causal claim is made.
- **PTP1B structures are ligand-bound in many cases**, which displaces water. Not
  controlled for; noted as a limit.

## What would make this wrong

Stated now so it can be checked later: if the identity match rate for any
structure falls below 95%, its conservation mapping is untrustworthy and it must
be excluded **before** the primary analysis, not after seeing its effect.

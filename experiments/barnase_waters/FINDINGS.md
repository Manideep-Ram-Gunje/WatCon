# Findings — barnase conserved water sites vs evolutionary conservation

**Result: H1 supported, at a moderate effect size, with one important
correction to the pre-registered statistic that the study's own control forced.**

Pre-registered in [PREREGISTRATION.md](PREREGISTRATION.md), committed before
`analyse.py` was written. Reproduce with:

```
python scripts/fetch_structures.py   # 47 barnase X-ray entries -> 24 at <=2.0 A
python scripts/prepare.py            # find barnase chain, superpose, carry waters
python scripts/run_study.py          # networks + ConSurf, cluster, join
python scripts/analyse.py            # the pre-registered analysis
```

---

## The headline

> Water sites that recur across many barnase crystal structures are lined by
> residues that are more conserved across evolution.
> **Spearman ρ = −0.375**, p = 6.6×10⁻¹⁷, n = 462 sites, against a permutation
> null centred on +0.033 (z = −4.47).

Correlational. One protein. Burial is not disentangled. Details below.

## Dataset

| | |
|---|---|
| Barnase X-ray entries (UniProt P00648) | 47 |
| At ≤ 2.0 Å | 24 |
| Passed preparation | **22** |
| Waters pooled | 2 753 |
| Cluster sites | 594 |
| Sites with ≥1 water | 590 |
| **Sites with conservation data** | **462** (floor was 50) |
| ConSurf runs used | **1** |

One ConSurf run covers all 22 structures because conservation is a property of
the *sequence*. That is what made the study possible with no new ConSurf data.

Superposition RMSD ranged 0.00–1.06 Å (median ≈ 0.30 Å) over 107–108 Cα.
Lowest per-structure ConSurf identity rate: **98.1%**.

## The correction that matters

The pre-registered primary statistic was `evo_min_score` — the *most* conserved
lining residue. It gives **ρ = −0.501, p = 9.4×10⁻³¹**, which looks emphatic.

Its own shuffled-conservation control says otherwise. **The null does not centre
on zero — it centres on −0.226.** Shuffling conservation at random still
produces a correlation of −0.23.

The reason is arithmetic, not biology:

| | |
|---|---|
| Spearman(occupancy, number of lining residues) | **+0.555** |
| Spearman(number of lining residues, min score) | **−0.553** |

The minimum of a larger set of draws is mechanically lower, and high-occupancy
sites are lined by more residues. So a good part of −0.50 is the statistic
measuring its own set size.

Replacing the minimum with the **mean** over lining residues removes that
dependence, and its null does centre on zero (+0.033):

| statistic | observed ρ | null mean | z | verdict |
|---|---|---|---|---|
| min (pre-registered) | −0.501 | **−0.226** | −3.29 | inflated |
| **mean (size-independent)** | **−0.375** | **+0.033** | **−4.47** | **the number to quote** |

Both are reported. The pre-registered result is not quietly replaced — it is
shown, and shown to be biased.

## Controls

| control | ρ (mean score) | p | n |
|---|---|---|---|
| all 22 structures | −0.375 | 6.6×10⁻¹⁷ | 462 |
| resolution ≤ 1.8 Å (7 structures) | −0.332 | 8.7×10⁻¹³ | 441 |
| wild-type-like only (6 structures) | −0.338 | 8.2×10⁻¹³ | 426 |

The effect is not an artefact of poorly resolved water, and it is not driven by
the mutant entries.

## The pre-registered secondary test could not be run

It compared recurrent sites against **singleton** sites. Of 21 sites occupied in
exactly one structure, only **1** has any scored lining residue — most singleton
waters have no protein residue in contact at all.

That is a finding rather than a defect: transient waters sit away from the
protein surface. But it leaves the pre-registered comparison with n = 1, so it is
reported as underpowered. A post-hoc substitute (≤3 structures vs ≥50%),
**labelled as post-hoc**, gives median mean-score **+0.401 vs −0.072**,
p = 1.6×10⁻⁹, rank-biserial r = −0.44.

## Confound, declared in advance and still not resolved

Buried residues are both more conserved and more likely to hold ordered water:

| | n | median ConSurf score |
|---|---|---|
| buried | 27 | **−0.694** |
| exposed | 81 | **+0.180** |

This design cannot separate "conserved residues hold water" from "buried
positions are conserved *and* hold water". **No causal claim is made.**
Disentangling it needs a burial-matched comparison, which is the obvious next
experiment.

## Limits

- **One protein.** This is barnase. It is not a family-level result.
- **Correlational**, with a known unresolved confound (above).
- **Crystallographic water only** — ordered, cryo-temperature, model-dependent.
- Some entries contain several barnase copies sharing crystallisation
  conditions, so their waters are not fully independent.
- 2F56 and 2F5M were **excluded**: their residue numbering is shifted by +2, so
  every ConSurf score would land on the wrong residue. Detected automatically,
  twice (see below). Renumbering them would recover two structures and is an
  untried sensitivity check.

## What the tooling caught

Worth recording, because these are the failures the integration was built to
prevent, and they occurred in real data rather than in tests:

1. **2F56 / 2F5M numbering offset.** Both entries number the residue everyone
   else calls 3 as 1. A positional sequence comparison scored them a perfect
   1.00 and let them in; they superposed at **6.05 Å**. Comparing by residue
   number — the way the ConSurf join compares — rejected them with the
   diagnosis *"matches at 1.00 if renumbered by +2"*.
2. **The same two, caught independently.** A stale `2F56.pdb` survived into the
   study run from an earlier pass. `enforce_identity` stopped it:
   *"only 2 of 106 compared residues agree (1.9%) … structure has ASN, ConSurf
   has VAL"*. Two mechanisms, built for different reasons, agreed.
3. **Three good structures nearly lost.** 1BSA, 1B2S and 1RNB start at residues
   4, 1 and 2 against the reference's 3. Positionally they scored 0.03–0.05 and
   were rejected; by residue number they agree at 0.98–0.99. Fixing the
   comparison recovered them.

## Next

- **Burial-matched re-analysis** — the one control that could turn this from a
  correlation into something stronger.
- **A real family.** The scaffold (`conservation_by_msa_column`) is built and
  tested; it needs one ConSurf run per family member.
  **Correction (2026-09-16):** this note previously gave a "measured ceiling" of
  ρ ≈ 0.37 for two runs of the same barnase sequence from different starting
  structures. That figure came from pairing the runs by ConSurf position, and the
  two are numbered differently (1BRS from 3, P00648 from 48), so half the pairs
  were different amino acids. Paired by sequence they agree at ρ = 0.97 (depth
  150) and 0.94 (depth 50). There is no 0.37 ceiling. The study's own result,
  ρ = −0.375 between site occupancy and conservation, is a different quantity and
  is unaffected.
- **Renumber 2F56/2F5M** and confirm the result is unchanged with 24 structures.

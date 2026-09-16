# Pre-registration — barnase conserved water sites vs evolutionary conservation

**Written before any analysis was run.** Committed before `analyse.py` was
executed for the first time; the commit history is the evidence. The pKa
experiment in `../pka/` was pre-registered the same way and returned a negative
result, which was reported as a result. This one is bound by the same rule.

---

## Question

WatCon finds water sites that recur across many structures of a protein —
*structural* conservation. ConSurf scores how conserved each residue is across
evolution — *evolutionary* conservation. The whole point of joining them is one
question:

> **Are structurally conserved water sites lined by evolutionarily conserved
> residues?**

It has been assumed in the literature far more often than it has been measured.

## Data

- **Protein**: barnase (UniProt P00648), a small, exhaustively crystallised
  model protein.
- **Structures**: all X-ray entries for P00648 at **≤ 2.0 Å** — 24 of 47 —
  selected mechanically by `scripts/fetch_structures.py` and listed in
  `data/selected.txt`.
- **Conservation**: one ConSurf run, `1BRS_A_150.grades.txt` (150 sequences).
  Conservation is a property of the *sequence*, so one run describes all 24
  structures. This is what makes the study possible without new ConSurf data.
- **Alignment**: Kabsch superposition (`WatCon.superpose`) onto the
  highest-resolution barnase chain, waters carried along.

## Hypothesis

**H1.** Water sites occupied in many structures are lined by residues with more
negative ConSurf scores (more conserved) than sites occupied in few.

**H0.** Occupancy and lining-residue conservation are unrelated.

## Primary analysis

Spearman rank correlation between, across all cluster sites carrying
conservation data:

- `occupancy_fraction` — fraction of structures in which the site is occupied
  (structural conservation), and
- `evo_min_score` — the most conserved lining residue's ConSurf score
  (evolutionary conservation).

**Predicted sign: negative.** ConSurf scores are more negative when more
conserved, so if H1 holds, higher occupancy pairs with lower (more negative)
score.

## Secondary analysis

Mann-Whitney U comparing `evo_min_score` between:

- **recurrent sites** — occupied in ≥ 50% of structures, and
- **singleton sites** — occupied in exactly one structure.

Reported with rank-biserial correlation as the effect size.

## Success criteria — fixed in advance

H1 is **supported** only if **both** hold on the primary analysis:

| | threshold |
|---|---|
| significance | p < 0.01 |
| effect size | \|ρ\| ≥ 0.2, **negative** |

A result that is significant but with \|ρ\| < 0.2 will be reported as "detectable
but negligible", not as support. A positive ρ of any size **contradicts** H1 and
will be reported as such.

**Power floor.** If fewer than **50** cluster sites carry conservation data, the
study is underpowered and will be reported as **inconclusive** rather than
negative — absence of evidence, not evidence of absence.

## Controls — all three run regardless of the primary outcome

1. **Shuffled-conservation null.** Permute conservation values across residues
   and recompute, 1000 times. The observed ρ must fall outside the null
   distribution. This catches an effect manufactured by the geometry of the join
   rather than by biology.
2. **Resolution restriction.** Repeat at ≤ 1.8 Å only. Water placement is
   resolution-dependent, and a result that only exists at poorer resolution is a
   modelling artefact.
3. **Wild-type only.** Many of the 24 are point mutants. Repeat on entries whose
   identity match rate is 100%. A result driven by mutants is a result about
   mutation, not conservation.

## Known confounds, stated in advance

- **Buried residues are both more conserved and more likely to coordinate
  ordered water.** Burial could produce a correlation with no causal link
  between conservation and water. ConSurf's own buried/exposed annotation will
  be reported alongside, but this study **cannot** disentangle it, and no causal
  claim will be made either way.
- **One protein.** Whatever comes out describes barnase. It is not a
  family-level result and will not be presented as one.
- **Copies within one crystal.** Some entries contain several barnase chains
  sharing crystallisation conditions; their waters are not fully independent.
- **Occupancy is bounded by coverage.** A site cannot be occupied in more
  structures than actually resolve that region.

## What would make this wrong

Stated so it can be checked later: if the identity match rate is below 95% for
any structure, its conservation mapping is not trustworthy and it must be
excluded before the primary analysis — not after seeing its effect on the
result.

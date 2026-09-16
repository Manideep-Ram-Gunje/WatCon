# Phase A — residue→alignment mapping: published WatCon vs ours

**Result: the bug is real, severe, and silent — but it does not bite in the
workflow the authors used, which is why it survived publication.** It does bite
on one of their own three published datasets, and on any raw PDB file.

Reproduce:

```
python scripts/fetch_data.py
python scripts/compare_mapping.py <structure directory> --label NAME
```

Compared: published WatCon at `89e992d` (`main` of the fork, untouched) against
this work. Dataset: the authors' own, Zenodo
[10.5281/zenodo.15213225](https://doi.org/10.5281/zenodo.15213225), CC-BY-4.0,
from Brownless, Harrison-Rawn & Kamerlin, *JACS Au* 2025.

---

## The mechanism

Published WatCon assigns a residue's alignment column with

```python
MSA_index = MSA_indices[molecule.resid-1]      # lines 714, 798, 816
```

`MSA_indices` comes from `generate_msa_alignment` and is **positional** — element
*k* is the alignment column of the *k*-th residue. So the correct index for the
residue at ordinal *k* is `MSA_indices[k]`, and `MSA_indices[resid-1]` is right
**only when `resid == k+1`**: numbering starts at 1, no gaps, single chain.

Nothing checks this, and nothing documents it as a requirement.

Four outcomes are possible, and only the last is visible to a user:

| outcome | what happens |
|---|---|
| correct | `resid-1 == k` |
| wrong, in range | returns a real column for the wrong residue — **silent** |
| wrong, negative wrap | `resid < 1` indexes from the end of the list — **silent** |
| out of range | `IndexError` — the only visible failure |

## Results

| dataset | source | first resid | **published correct** | ours |
|---|---|---|---|---|
| PTP1B (`pdbs_all`) | authors', aligned | 1 | **99.9%** (62 929/62 996) | 0 wrong |
| PTPs combined | authors', aligned | 1 | **99.0%** (6 785/6 855) | 0 wrong |
| **TPIs combined** | **authors', as shipped** | −7 … 7 | **11.2%** (1 487/13 315) | 0 wrong |
| **Raw PTP1B 1AAX** | RCSB, unmodified | 2 | **0%** (0/297) | 0 wrong |
| **Raw barnase, 24** | RCSB, unmodified | 1–4 | **3.9%** (326/8 351) | 0 wrong |

## Why it survived publication

The authors' PTP1B structures **are renumbered by their preprocessing**. Compare
their `1AAX_aligned.pdb` with the raw RCSB entry:

```
their aligned 1AAX : 297 residues, first resid 1
raw RCSB 1AAX      : 297 residues, first resid 2
identical (resid, residue) sequence?  False
  first 5 theirs: [(1,'E'), (2,'M'), (3,'E'), (4,'K'), (5,'E')]
  first 5 raw   : [(2,'E'), (3,'M'), (4,'E'), (5,'K'), (6,'E')]
```

Same protein, same residues, shifted by one. Their pipeline produces exactly the
input the arithmetic requires, so within their workflow the code is correct and
the published PTP1B results are unaffected.

**Their TPI dataset is not renumbered**, and there the same code maps 11.2% of
residues correctly — 3 of 22 structures perfectly, one structure numbered from
−7. That dataset is part of the same publication.

## What this does and does not say

**Does:**
- Published WatCon is correct **only for input renumbered to start at 1**, and
  that requirement is neither documented nor enforced.
- A user pointing it at raw PDB files gets silently wrong alignment columns —
  0% on raw 1AAX, 3.9% on raw barnase.
- The failure is invisible: a real column is returned for the wrong residue.
- This work removes the requirement. `residue_index.ResidueIndex` pairs residues
  with columns by ordinal and **refuses to build** when the lengths disagree.
  0 residues mapped incorrectly across all 91 517 tested here.

**Does not:**
- Say the published PTP1B or PTP results are wrong. They are not — their
  preprocessing satisfies the precondition.
- Say anything about water-network quality. That algorithm is theirs and
  unchanged; this measures one indexing step.
- Establish that our approach is better *designed*. It fixes a specific defect.

## Correction (2026-09-16): how the PTP figure was derived

The **99.0%** row above is right, but it was derived the wrong way and its
companion column was empty of content. Both are corrected here.

**How it was measured before.** `compare_mapping.py` classified a residue as
correctly mapped when `resid - 1 == ordinal`, walking residues with a reader that
skipped HETATM records. In 5HDE that reader missed CSP230, the catalytic
phosphocysteine, so the walk had a hole at 230 and every later residue failed the
`resid - 1 == ordinal` test. The arithmetic was about our own reader's numbering,
not about the alignment.

**Measured properly.** Ground truth is now the column each residue actually
occupies, obtained by aligning each structure to its own alignment row by
sequence (`WatCon.alignment.map_structure_to_row`, identity 1.000 for every one
of the 24 rows), with modified residues read. Published WatCon's
`MSA_indices[resid-1]` is then compared against it:

| dataset | residues | correct | wrong | IndexError |
|---|---|---|---|---|
| PTPs combined | 6,856 | **6,785 (98.96%)** | 70 | 1 |

Every error is in **5HDE**, and the cause is in the authors' own alignment: its
5HDE row omits CSP230, so the column list is one short. Residues 231 onward are
assigned the previous residue's column, and the last residue raises `IndexError`.
The count is the same 6,785 as before, reached for the right reason.

**The "ours 0 wrong" column meant less than it looked.** `ours_is_correct` handed
`ResidueIndex` a stand-in column list, `range(n)`, and asked whether it returned
the ordinal it was given. That tests the identity-to-ordinal lookup, not agreement
with a real alignment. Against the real rows our mapping places 5HDE's 299 other
residues correctly and gives CSP231 **no** column rather than a wrong one, which
is the behaviour the anchor tests in `WatCon/tests/test_alignment_mapping.py`
check independently.

## A note on fairness

`compare_mapping.py` grades published WatCon against a ground truth derived the
way our code does it, which is not self-evidently neutral. Two guards:

1. The measurement is **arithmetic, not judgement** — whether `resid-1 == k`.
   The worked example printed for each dataset can be checked by hand.
2. The script **asserts our own mapping is correct** rather than assuming it, by
   querying `ResidueIndex` for every residue and comparing against the ordinal.
   It reports our error count in the same table.

Verified by hand on raw 1AAX: residues are numbered 2,3,4,… at ordinals 0,1,2,…
so `resid-1` is always one past the correct slot — 0/297, as reported.

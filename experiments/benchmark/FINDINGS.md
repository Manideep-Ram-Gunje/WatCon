# WatCon-ConSurf vs published WatCon — what actually holds up

**Bottom line.** The predictive claim fails on two independent proteins.
Conservation adds nothing to predicting where water sits. What survives is
narrower and real: a correctness fix, and a way to *classify* water sites that
occupancy alone conflates.

Reproduce:

```
python scripts/fetch_data.py
python scripts/compare_mapping.py <structures> --label NAME     # Phase A
python scripts/heldout_benchmark.py --prepared ... --consurf ... # Phase B
python scripts/disagreement.py --prepared ... --consurf ...      # descriptive
```

Compared against published WatCon at `89e992d`, on the authors' own dataset
(Zenodo [10.5281/zenodo.15213225](https://doi.org/10.5281/zenodo.15213225),
CC-BY-4.0, Brownless, Harrison-Rawn & Kamerlin, *JACS Au* 2025).

---

## 1. "Is our tool better?" is not one question

The water-network algorithm is **theirs and unchanged**. We added an orthogonal
data layer. There is no shared metric on which the two compete at the same task,
and any single "we win" number would be manufactured. Three separable questions
were asked instead, and they have different answers.

---

## 2. Correctness — a real defect, latent in the authors' own workflow

Published WatCon assigns alignment columns with `MSA_indices[molecule.resid-1]`
(lines 714, 798, 816). `MSA_indices` is **positional**, so this is right only
when `resid == ordinal+1`: numbering starts at 1, no gaps, single chain.

| dataset | source | first resid | **published** | ours |
|---|---|---|---|---|
| PTP1B | authors', aligned | 1 | **99.9%** (62 929/62 996) | 0 wrong |
| PTPs combined | authors', aligned | 1 | **99.0%** | 0 wrong |
| **TPIs combined** | **authors', as shipped** | −7 … 7 | **11.2%** (1 487/13 315) | 0 wrong |
| Raw 1AAX | RCSB, unmodified | 2 | **0%** (0/297) | 0 wrong |
| Raw barnase ×24 | RCSB, unmodified | 1–4 | **3.9%** | 0 wrong |

**Their published PTP1B results are not wrong.** Their preprocessing renumbers
structures — their `1AAX_aligned.pdb` starts at residue 1 where the raw RCSB
entry starts at 2 — so their pipeline produces exactly the input the arithmetic
needs. Their **TPI** dataset is not renumbered, and there the same code maps
11.2% of residues correctly. That dataset is part of the same publication.

The honest claim is narrow: *published WatCon is correct only for input
renumbered to start at 1; that requirement is neither documented nor enforced;
and the failure is silent — a real alignment column returned for the wrong
residue.*

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

## 3. Prediction — conservation adds nothing

Pre-registered in [PREREGISTRATION.md](PREREGISTRATION.md), committed before the
analysis script existed. Leave-one-structure-out: build sites from N−1
structures, predict whether the held-out one has a water there.

| | barnase | PTP1B seed 0 | PTP1B seed 1 |
|---|---|---|---|
| structures | 22 | 50 of 253 | 50 of 253 |
| M1 — occupancy + site size | 0.8676 | 0.9006 | 0.9030 |
| M2 — **+ conservation** | 0.8679 | 0.9007 | 0.9029 |
| **ΔAUC** | **+0.0004** | **+0.0001** | **−0.0002** |
| Wilcoxon p | 0.388 | 0.404 | 0.0029 |
| shuffled null | −0.0002 | −0.0001 | −0.0001 |
| observed vs null | **inside** | **inside** | **inside** |
| conservation **alone** | 0.6514 | 0.5694 | 0.5671 |

Pre-registered threshold: ΔAUC ≥ 0.02 **and** p < 0.05. Not met, not close, on
two independent proteins and two independent subsamples.

Seed 1 is instructive: p = 0.0029 on a **negative** delta. With 50 paired folds,
significance says almost nothing — the effect size and the null are what matter.

**Conservation is redundant, not uninformative.** Alone it scores AUC 0.65
(barnase) and 0.57 (PTP1B), clearly above chance. It carries real signal about
where water sits; it carries nothing that occupancy does not already carry.

### The control earned its place twice

1. **It caught a bug in this analysis.** The first implementation fitted and
   scored the regression on the same fold, so a third feature improved fit
   regardless of content. It reported ΔAUC = **+0.0017 at p = 0.0008** — which
   reads as a win. The shuffled null centred on **+0.0016**, essentially
   identical, which is what exposed it. Each fold is now scored by a model
   trained on the other folds, and the null moved to −0.0002.
2. Without it, this document would have claimed a positive result.

## 4. Interpretation — what the tool is actually for

Two correlated measures can be redundant for prediction while still disagreeing
about individual sites. Splitting barnase's 462 scored sites at the median of
each measure:

| group | barnase, n | buried fraction | PTP1B, n | buried fraction |
|---|---|---|---|---|
| recurrent **and** conserved | 136 | **0.257** | 615 | **0.219** |
| conserved only | 95 | 0.119 | 450 | 0.107 |
| **recurrent only** | 72 | **0.059** | 454 | **0.036** |
| neither | 159 | 0.010 | 621 | 0.029 |
| **disagreement** | | **36%** (167/462) | | **42%** (904/2140) |
| Mann-Whitney p | | 2.5×10⁻⁸ | | **1.3×10⁻³³** |

Both proteins agree: sites that recur *without* conservation are overwhelmingly
**surface** sites, by ConSurf's own buried/exposed call. Those are packing and
crystallisation waters. Occupancy alone cannot tell them apart from buried,
functionally-lined waters -- and this is the one thing the tool does that
published WatCon cannot.

**Two caveats, stated next to the result rather than below it:**

- Sites are split at the **median**, so the 36% disagreement rate is the
  arithmetic consequence of the ρ = −0.375 correlation, not a discovery. What is
  worth knowing is what the disagreeing groups look like, not how many there are.
- **Burial and conservation are correlated** (buried median −0.694 vs exposed
  +0.180). Conservation is partly acting as a **burial proxy**, and a user might
  get comparable separation from a burial calculation alone. This is the same
  confound declared at the start of the barnase study, still unresolved.

## 5. Engineering — a defect found only at scale

`conservation_of_clusters` compared every centre against every water of every
structure in a Python triple loop: 50 PTP1B structures give ~2 300 centres over
~11 000 waters, i.e. **27 million distance evaluations per join**. Measured
**40–60 s before, 1.55 s after** a KD-tree per structure — about **30×**.

Results are unchanged, not merely similar: the bundled demo returns the same
193/190/165 sites, the same 30 occupied in all six structures, the same 57 lined
by a grade ≥ 8 residue. 476 tests pass; CI green on 9 platform/version
combinations.

It surfaced only by running on 253 structures instead of the handful the tool had
been tested on.

---

## 6. The family at fifteen (2026-09-16)

Ten more ConSurf runs took the family from five proteins to fifteen (24
structures). The runs are in `data/consurf_raw/`, mapped to proteins in
`MANIFEST.md`; the analysis is reproducible from `ptp_family15/members15.tsv`.

**Verification first.** All fifteen runs parse strictly, share identical
settings (webserver, Bayesian, 150 homologues, chain A), and agree with
ConSurf's own annotated PDBs **4,289 / 4,289**. Each of the ten new runs
reproduces its row of the authors' published alignment exactly, residue for
residue. Across all 24 structures the identity check flags **only** declared
engineered substitutions — the C→S traps, PTPN9's C→A, the D→A general-acid
mutants, 5J8R's phosphocysteine — and no false alarms.

**The conservation statistic does not decay with more runs.** This was the real
question: `unanimous_conserved` counts columns every run independently grades
≥8, so each protein added is another chance for a column to drop out.

| | 5 proteins | 15 proteins |
|---|---|---|
| alignment columns | 337 | 344 |
| covered by all | 248 | 229 |
| unanimously conserved | 64 (25.8%) | 58 (25.3%) |

Tripling the number of independently normalised runs moved the fraction by half
a percentage point. The five-protein result was not an artefact of a small,
closely related sample.

**The site-level separation gets much sharper.**

| | 5 proteins | 15 proteins |
|---|---|---|
| occupied sites | 308 | 556 |
| lined by a unanimous column | 81 | 117 |
| mean proteins holding those | 3.43 | **7.34** |
| mean proteins holding the rest | 2.66 | **4.08** |
| Mann–Whitney one-sided p | 7.6 × 10⁻⁸ | **1.0 × 10⁻²⁴** |

**What does not survive is universality.** At five proteins, 22 sites were held
by all five, including the catalytic water. At fifteen, **exactly one** site is
held by all fifteen — and it is a buried structural water lined by PTP1B 56/57/67,
not a catalytic one. The catalytic sites are still found and still
conserved-lined (WPD Asp181 in 7–8 proteins, Q-loop Gln262 in 8, P-loop Arg221
in 6, all grade 9), but none is universal across a family spanning open and
closed WPD-loop states, ten catalytically dead entries and a pseudophosphatase.

The earlier "a water position held by all five, lined by the residues each
protein uses for catalysis" was true of five closely related classical PTPs and
**does not generalise**. It is corrected here rather than quietly dropped.

**PTPRN2 makes no difference.** 2QEP is a pseudophosphatase (`CSDGAGR` where the
family reads `CSAGIGR`). Run both ways: 58 unanimous columns with or without it;
p = 1.0e-24 against 4.5e-25. Excluding it restores three universal sites (4 of
14 rather than 1 of 15), which is the only thing it changes. It is kept, because
it is the clearest case in the set of conservation measuring *constraint* rather
than *identity*: at PTP1B's A217 the family has alanine and PTPRN2 aspartate,
and every run — including PTPRN2's own — grades that column 9.

The authors' dataset labels 2QEP "PTPN2". The entry is **PTPRN2** (Q92932); real
PTPN2 is P17706.

## 7. The dynamic path, run for the first time (2026-09-16)

Running `generate_dynamic_networks.initialize_network` on real input found four
defects, one of which meant it could never have read a crystal structure:
`add_water` required hydrogens positionally. Details in
`tools/WatCon_ConSurf/docs/CONSURF_CHANGELOG.md`, Phase 27.

Two results worth keeping here:

**Cost is not where it looks.** Profiling one frame of the Zenodo MD system over
957 s: `get_shortest_path` 819 s, `get_CPL` 45 s, **building the network 1.4 s**.
Both metrics run all-pairs shortest paths and are on by default. Off, the
identical network takes 1.0 s instead of 757 s.

**An ensemble of crystal structures is not a trajectory.** Each entry resolves a
different number of waters, so no constant atom count exists; MDAnalysis refuses
it, and the tool now refuses it up front with a pointer to the static path. No
MD trajectory exists in this project to substitute — Zenodo ships starting
structures only — so the directed path is validated on one real MD frame.

### The join could attach nothing and report success

Two defects in the identity key, both found by attaching the 1AAX run to the MD
system rather than by reading the code.

**The chain was read from the wrong attribute.** `chainID` is an *atom*-level
attribute in MDAnalysis, so reading it off a residue always returned `None` and
the code fell back to `segid`. That is right by luck for files straight from the
PDB, where MDAnalysis fills segid from the chain column — which is why every
result in this project survived it. For a PDB *written* by MDAnalysis the segid
is the invented `SYST` while the real chain sits in `chainID`, so lookups missed
entirely: **0 of 127** active-site atoms scored, run reported success.

**Protonation states read as sequence differences.** A force field writes `HID`,
`HIE`, `HIP` where a crystal structure writes `HIS`; the check compared
three-letter names. Chemical modifications are deliberately still reported, so
5HDE's phosphocysteine remains `CYS>CSP`.

After both fixes, on the MD active site renumbered by the declared +1 offset:
identity 0.9545 → **0.9773**, the only remaining difference being the real one
(1AAX is the C215S trap, the MD system models a cysteine), and **127 of 127**
atoms scored — recovering the WPD loop at grade 9, Trp179/Pro180/**Asp181**.

**No number in this document moves.** Every result here was measured on RCSB
files, where the old path happened to give the right chain. The bug mattered for
inputs this project had not previously used.

A third case now warns instead of passing silently: a ConSurf run whose chain
matches nothing scored 0 of 894 atoms and said nothing, because coverage is
reported rather than enforced. It still is — low coverage is legitimate when a
run covers one chain of several — but *zero* coverage now prints what the two
chains are and how to reconcile them.

## What to claim, and what not to

**Claim:**
- Published WatCon silently mis-maps residues unless input is renumbered from 1.
  Measured on their own data: 11.2% correct on their TPI set, 0% on raw 1AAX.
- This work removes that requirement — 0 residues mis-mapped across 91 517 tested.
- Conservation and occupancy disagree at about a third of sites, and the
  disagreement separates surface packing waters from buried ones.
- The tool runs at family scale: 253 structures, 60 750 waters; and across
  **fifteen** different proteins with one ConSurf run each.
- Sites lined by a unanimously conserved column are held by more proteins
  than the rest: 7.34 against 4.08 of 15, p = 1.0e-24.

**Do not claim:**
- That conservation improves prediction of water occupancy. **It does not** —
  ΔAUC ≈ 0 on two proteins.
- That our water networks are better. That algorithm is theirs and unchanged.
- That the published PTP1B or PTP results are wrong. They are not.
- That the site classification is causal. Burial is not disentangled.
- That a catalytic water site is common to the whole PTP family. At five
  proteins it looked that way; at fifteen, only one site is universal and it
  is structural.

## Limits

- **Two proteins**, both single-domain, both with a single ConSurf run each.
- PTP1B used **50 of 253** structures per benchmark run, for tractability;
  clustering once and reusing centres would have been faster but would leak the
  held-out structure into its own site definitions. Two seeds agree.
- **4.1% of PTP1B waters carry altloc flags** (63 of 253 structures); WatCon keeps
  the first instance, discarding real occupancy information.
- 4 of 253 structures superpose worse than 2.0 Å (median 1.00 Å) — the open/closed
  WPD-loop difference the authors themselves studied.
- Many PTP1B entries are **ligand-bound**, which displaces water. Not controlled.
- The family comparison shares the burial confound: no burial control was run.
- The **trajectory path has not been run on real molecular dynamics** — no
  trajectory exists in this project. Coverage is a three-frame crystal
  ensemble and one MD frame.

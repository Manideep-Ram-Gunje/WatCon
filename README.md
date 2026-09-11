# WatCon-ConSurf

[![CI](https://github.com/Manideep-Ram-Gunje/WatCon/actions/workflows/CI.yaml/badge.svg?branch=consurf-integration)](https://github.com/Manideep-Ram-Gunje/WatCon/actions/workflows/CI.yaml)

**Conserved water networks, joined to evolutionary conservation.**

WatCon finds water sites that recur across many structures of a protein —
*structural* conservation. [ConSurf](https://consurf.tau.ac.il/) scores how
conserved each residue is across evolution — *evolutionary* conservation. This
package joins the two, so you can ask a question that has been assumed far more
often than it has been measured:

> Are structurally conserved water sites lined by evolutionarily conserved
> residues?

On 22 barnase crystal structures, the answer is **yes, moderately**:
Spearman ρ = −0.375 (p = 6.6×10⁻¹⁷, n = 462 sites) against a permutation null
centred on +0.033. [Full result, method and caveats →](experiments/barnase_waters/FINDINGS.md)

---

## This is a derivative work

This package **extends** [WatCon](https://github.com/kamerlinlab/WatCon) by
Brownless, Harrison-Rawn & Kamerlin, published in *JACS Au* 2025. WatCon is
theirs; the water-network analysis, clustering, MSA machinery and PyMOL
projections are their work, and it is excellent. Please cite it:

> Brownless A-LR, Harrison-Rawn T, Kamerlin SCL. *WatCon: A Python Tool for
> Analysis of Conserved Water Networks Across Protein Families.* JACS Au 2025.
> DOI: [10.1021/jacsau.5c00447](https://pubs.acs.org/doi/10.1021/jacsau.5c00447)

**What is added here** is the ConSurf integration and the infrastructure it
needed: a ConSurf parser, residue-identity mapping, conservation attached to
residues and waters, the site↔conservation join, MODELLER-free superposition, a
dataset preparation step, and a command-line interface.
Licensed GPL-3.0, like the original. See [CHANGES](docs/CONSURF_CHANGELOG.md).

---

## Install

```bash
pip install "git+https://github.com/Manideep-Ram-Gunje/WatCon.git@consurf-integration"
```

> **The `@consurf-integration` part is required.** This fork's `main` is a clean
> mirror of upstream WatCon, which does *not* contain the ConSurf integration —
> an install URL without the branch gives you plain WatCon.

Python ≥3.10. Everything needed is installed automatically.

**MODELLER is optional.** Upstream WatCon required it for structural alignment;
it is licensed and conda-only. `WatCon.superpose` now does the same job in numpy
for structures of one protein, so you only need MODELLER for MSA-based alignment
across a protein *family*. See [installation docs](docs/installation.rst).

## Try it in one command

```bash
watcon demo
```

Runs the whole pipeline on six real barnase crystal structures bundled with the
package — offline, in about ten seconds:

```
STEP 1/5  Prepare: find the barnase chain, superpose, keep waters
  1A2P     chain A  identity 1.00   148 waters  RMSD 0.00 A over 108 CA
  1BRN     chain L  identity 1.00   120 waters  RMSD 0.45 A over 108 CA
  ...
STEP 4/5  Cluster recurring water sites and join to conservation
  193 cluster centres, 190 occupied, 165 carry conservation
  30 site(s) occupied in ALL 6 structures
  57 scored site(s) are lined by a highly conserved residue (ConSurf grade >= 8)
```

Note that 1BRN's barnase is chain **L**, not A — the example is deliberately
messy, because real depositions are.

## The commands

```bash
watcon prepare  --input-dir raw/ --out-dir prepared/ --reference 1A2P
watcon run      --input input.txt --analysis analysis.txt
watcon validate --consurf my_run_consurf_grades.txt
watcon view     --prepared prepared/ --consurf grades.txt
watcon demo
```

## See it

```bash
watcon view --prepared prepared/ --consurf my_run_consurf_grades.txt
pymol watcon_view/watcon_view.pml
```

The protein is coloured on **ConSurf's own 1-9 scale** (maroon conserved, cyan
variable) and the water sites lined by a highly conserved residue are marked.
Type `enable sites` in PyMOL to add every occupied site, coloured by grade.

The script is self-contained -- it loads what it colours, so opening it is all
you do.

`python WatCon/WatCon.py --input input.txt` still works exactly as before.

## Using your own data

1. **Get structures.** Any PDB files of the same protein.
2. **Get ConSurf results.** There is no public ConSurf API, so this step is
   manual — [step-by-step guide](docs/consurf_data.rst). Keep the
   `*_consurf_grades.txt` file and name it after the structure.
   Conservation is a property of the **sequence**, so *one* ConSurf run
   describes every structure of that protein.
3. **Check it before you rely on it:** `watcon validate --consurf FILE`
4. **Prepare:** `watcon prepare --input-dir raw/ --out-dir prepared/`
5. **Run:** set `consurf_directory:` in your input file and
   `conservation_report: on` in the analysis file, then `watcon run`.

You get a CSV with one row per conserved water site: its structural occupancy
beside the evolutionary conservation of the residues lining it — **two separate
columns**. Deciding how they relate is the science, so the tool does not blend
them into one number.

## What it checks for you

The join key is `(chain, resid, icode)`, which says nothing about *which amino
acid* a score belongs to. A ConSurf run numbered differently from your structure
would attach a score to every residue, report 100% coverage, and be wrong
everywhere — invisibly. So the identity of every matched residue is verified:

```
ConSurf data does not describe 2F56.pdb: only 2 of 106 compared residues
agree (1.9%, threshold 95%). ... A3: structure has ASN, ConSurf has VAL
```

That is a real message from the barnase study — 2F56 numbers the residue
everyone else calls 3 as 1. It is a *rate*, not a pass/fail, because a point
mutant genuinely differs at one position (99.1%) while a numbering offset
collapses to 2.8%.

## Limits

Stated plainly, because they matter more than the headline:

- The barnase result is **one protein** and **correlational**.
- **Burial is not disentangled**: buried residues are both more conserved
  (median −0.694 vs +0.180) and more likely to hold ordered water. No causal
  claim is made.
- Cross-run ConSurf agreement is the ceiling on any family-level claim: two runs
  of the *same* barnase sequence from different starting structures agree at only
  ρ ≈ 0.37. Effects smaller than that are inside ConSurf's own noise.
- The multi-protein family path (`conservation_by_msa_column`) is built and
  tested, but has never been run on a real family — the join is tested, the
  biology is not.

## Documentation

- [Getting ConSurf data](docs/consurf_data.rst) — the manual step, in detail
- [Conservation tutorial](docs/tutorials/consurf_conservation.rst)
- [Design and decisions](docs/CONSURF_INTEGRATION.md)
- [Change log](docs/CONSURF_CHANGELOG.md)
- Upstream: https://watcon.readthedocs.io/

## Tests

```bash
pip install -e ".[test]"
python -m pytest WatCon/tests -q      # 476 tests
```

## License

GPL-3.0-only, inherited from WatCon.
Copyright (c) 2025 Alfie-Louise Brownless, Shina Caroline Lynn Kamerlin, and
contributors to this extension.

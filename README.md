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

## Use it inside PyMOL

```bash
watcon plugin --install
```

Restart PyMOL, then **Plugin → WatCon + ConSurf**. Choose a folder of
structures (or type PDB ids and let it fetch them), pick the reference structure
and chain from the dropdowns, point it at your ConSurf file, and press **Run**.

The session is drawn in the viewport you are already looking at, and a sortable
table of every site appears beside it — site, ConSurf grade, how many structures
it is occupied in, and which residues line it. **Click a row** and the camera
flies to that site and shows its lining side chains.

The analysis runs on a background thread, so PyMOL stays usable while it works.

This needs PyMOL and WatCon in the same Python, which is what a normal
`pip install` of both gives you.

## The commands

```bash
watcon fetch    --ids 1AAX 7GSA --out-dir raw/
watcon prepare  --input-dir raw/ --out-dir prepared/ --reference 1AAX
watcon run      --input input.txt --analysis analysis.txt
watcon validate --consurf my_run_consurf_grades.txt
watcon view     --prepared prepared/ --consurf grades.txt
watcon family   --members members.tsv --alignment alignment.pir
watcon plugin   --install
watcon demo
```

PDB and **mmCIF** are both read, gzipped or not. That matters more than it
sounds: RCSB no longer issues PDB files for large or recent entries — 31 of the
287 PTP1B structures used to test this have no PDB file at all.

## See it without the plugin

```bash
watcon view --prepared prepared/ --consurf my_run_consurf_grades.txt
pymol watcon_view/watcon_view.pml
```

The protein is coloured on **ConSurf's own 1-9 scale** — maroon conserved, cyan
variable, and **yellow where ConSurf has no score**, which is not the same as
low conservation. Water sites lined by a highly conserved residue are marked,
and sphere size tracks how many structures hold a water there.

```
enable sites              every occupied site, coloured by grade
enable WatCon_contacts    the residues lining them, and their polar contacts
```

The script is self-contained — it loads what it colours, so opening it is all
you do. The plugin runs this identical file, so the two cannot disagree.

## Does it work on real data?

Yes, and here is the honest version. On **253 PTP1B crystal structures** with one
real ConSurf run, in about two minutes, it finds 295 recurring water sites, 92 of
them lined by a residue ConSurf grades 8 or 9. Among those:

* a site occupied in **213/253** structures, lined by the entire P-loop (Cys215
  the nucleophile, Ser216, Ala217, Gly218, Ile219, Gly220) **and Gln262** — the
  catalytic water position;
* two sites occupied in **192/253**, lined by **Gln262**, whose job is to
  position the catalytic water;
* sites occupied in **197** and **191/253** on the WPD loop at **Asp181**, the
  general acid, together with Arg221.

It found those with no knowledge of PTP1B's chemistry. But conservation does
**not** put them at the top: ranked by conservation then occupancy, the first
ten sites are buried structural waters, and the catalytic ones only accumulate
by depth 50 (19/50, against 4/50 by occupancy alone). Conservation re-ranks
toward the active site; it is not a shortcut to it.

### A family of proteins

The same machinery runs across *different* proteins, each with its own ConSurf
run, placed on a shared alignment:

```bash
watcon family --members members.tsv --alignment alignment.pir
```

On **fifteen** protein tyrosine phosphatases (twenty-four structures, one
ConSurf run each), it puts them in a single frame at 0.65–0.97 Å.

Across the 556 occupied sites, those lined by a residue **every** ConSurf run
calls conserved are held by far more proteins than the rest: **7.34 against
4.08** on average, p = 1.0 × 10⁻²⁴. The separation is much sharper than it was
at five proteins (3.43 vs 2.66, p = 7.6 × 10⁻⁸).

The conservation statement holds up too. Going from five runs to fifteen leaves
the fraction of universally covered alignment columns that *every* run calls
conserved essentially unchanged — 25.8% (64/248) to 25.3% (58/229) — although
unanimity across fifteen independently normalised runs is a far stricter bar.

What does **not** survive is universality. At five proteins, 22 water sites were
held by all five, including the catalytic one. At fifteen, exactly one site is
held by all fifteen, and it is a buried structural water. The catalytic waters
are still found and still conserved-lined — the WPD aspartate, the Q-loop
glutamine and the P-loop arginine, all grade 9 — but none is universal across a
family spanning open and closed states, ten catalytically dead entries and a
pseudophosphatase. Burial is still not disentangled, so part of the headline may
be a burial effect.

Site numbers are cluster labels from one run, not residue numbers, so they are
not quoted here. Residues with alternate conformations use the most populated
position, keeping every position tied for most populated. Full numbers in
[docs/CONSURF_CHANGELOG.md](docs/CONSURF_CHANGELOG.md).

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
- ConSurf scores are normalised within each run, so separate runs are only
  relatively comparable. Measured agreement between independent runs of barnase,
  paired by sequence: ρ = 0.955 (same sequence, two MSA depths) and 0.94–0.97
  (different starting structures). An earlier version of this README gave
  ρ ≈ 0.37 for the latter; that figure came from pairing the two runs by
  position across different numbering, so half the pairs were different amino
  acids. There is no 0.37 floor.
- The family result is **fifteen proteins**, all protein tyrosine phosphatases,
  and its conserved-versus-other comparison shares the burial confound above.
- One of the fifteen, PTPRN2, is a **pseudophosphatase**. It is kept because
  excluding it changes the conservation result not at all (58 unanimous columns
  either way); both versions are reported in the change log.
- Alternate conformers: positions tied for most populated are all kept, so
  WatCon's own edge counts stay inflated for structures modelling many equally
  occupied positions.
- The **trajectory path has not been run on real molecular dynamics**, because
  no trajectory was available: it is covered by a three-frame crystal ensemble
  and a single MD frame. An ensemble of crystal structures is *not* a
  trajectory — each entry resolves different waters — and the tool now refuses
  one with a message pointing at the static path.
- Two graph metrics, `shortest_path` and `characteristic_path_length`, run
  all-pairs shortest paths and are on by default. On a 4,851-water system they
  took 819 s and 45 s against **1.4 s** to build the network itself. Turn them
  off for anything large; the network is unchanged.

## Documentation

- [The barnase study](experiments/barnase_waters/) — the evidence for the
  headline result, with the scripts that reproduce it
- [The PTP1B and family benchmark](experiments/benchmark/) — findings, method
  and numbers for the 253-structure and fifteen-protein results
- [Getting ConSurf data](docs/consurf_data.rst) — the manual step, in detail
- [Conservation options](docs/faq/conservation_options.rst) — every setting the
  extension adds, and what each one refuses to guess
- [Conservation tutorial](docs/tutorials/consurf_conservation.rst)
- [Design and decisions](docs/CONSURF_INTEGRATION.md)
- [Change log](docs/CONSURF_CHANGELOG.md)
- Upstream: https://watcon.readthedocs.io/

## Tests

```bash
pip install -e ".[test]"
python -m pytest WatCon/tests -q      # 792 tests
```

Run it from a **repository checkout**. The suite also ships inside the wheel, so
`pytest --pyargs WatCon.tests` works on an installed copy, but ten modules need
inputs a distribution deliberately does not carry -- the raw ConSurf bundles and
the older core tests' structures. Those are not collected there, and a warning
says which.

## License

GPL-3.0-only, inherited from WatCon.
Copyright (c) 2025 Alfie-Louise Brownless, Shina Caroline Lynn Kamerlin, and
contributors to this extension.

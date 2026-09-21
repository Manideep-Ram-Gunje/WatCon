# Testing WatCon-ConSurf by hand

Every command below was run before this page was written, and every expected
output is copied from what actually appeared. If something differs on your
machine, that is a real difference worth looking at — not a typo here.

Three parts: a **10-minute smoke test** that proves the whole thing works, a
**full pass** over every command and every guard, and a **demo sequence** for
showing someone.

Work in a scratch directory, not the repository:

```bash
mkdir ~/watcon-testing && cd ~/watcon-testing
```

---

## Part 1 — Smoke test (~10 minutes)

### 1. The package is installed and knows its own version

```bash
watcon --version
```

```
watcon-consurf 0.9.0
```

A bare `1+unknown` means you are running an older build.

### 2. The whole pipeline, on data that ships with the package

```bash
watcon demo
```

Takes about ten seconds, offline. The parts to check:

```
STEP 1/5  Prepare: find the barnase chain, superpose, keep waters
  1A2P     chain A  identity 1.00   148 waters  RMSD 0.00 A over 108 CA
  1BRN     chain L  identity 1.00   120 waters  RMSD 0.45 A over 108 CA
  1BRS     chain A  identity 1.00   174 waters  RMSD 0.46 A over 108 CA
  1BSA     chain A  identity 0.98    81 waters  RMSD 0.16 A over 107 CA
  1BSE     chain A  identity 0.99    59 waters  RMSD 0.24 A over 108 CA
  1RNB     chain A  identity 0.98    89 waters  RMSD 0.44 A over 108 CA
...
STEP 4/5  Cluster recurring water sites and join to conservation
  193 cluster centres, 190 occupied, 165 carry conservation
  30 site(s) occupied in ALL 6 structures
  57 scored site(s) are lined by a highly conserved residue (ConSurf grade >= 8)
```

**Three things worth noticing.** 1BRN's barnase is chain **L**, not A — the
example is deliberately awkward because real depositions are. 1BSA and 1RNB
score identity 0.98 because they are point mutants, and are kept; a numbering
error would collapse to near zero instead. And `193 → 190 → 165` is not
attrition: 190 sites are occupied, and 165 of those have a residue near enough
to carry conservation.

### 3. Look at it in PyMOL

```bash
pymol watcon_demo/view/watcon_view.pml
```

The protein is coloured on **ConSurf's own 1–9 scale** — maroon conserved, cyan
variable, **yellow where ConSurf has no score**, which is not the same as low
conservation. Then, in the PyMOL prompt:

```
enable sites              every occupied site, sphere size = how many structures hold it
enable WatCon_contacts    the residues lining them, and their polar contacts
```

No display handy? Check it loads without one:

```bash
pymol -cq watcon_demo/view/watcon_view.pml -d "print(cmd.get_names('objects'))"
```

```
['protein', 'sites', 'sites_conserved', 'site_contacts', 'WatCon_contacts']
```

### 4. The plugin

```bash
watcon plugin --install
```

Restart PyMOL, then **Plugin → WatCon + ConSurf**. In the dialog:

* **Structures** → `watcon_demo/prepared`
* **ConSurf file** → `watcon_demo/consurf/1A2P_consurf_grades.txt`
* press **Run**

The analysis runs on a background thread, so PyMOL stays usable. When it
finishes you get the session from step 3 plus a sortable table of every site.

### 5. Click a row — this is the part worth seeing

The camera flies to that site and shows the side chains lining it.

**This is the answer to "you cannot compare 3D structures by eye".** You are not
meant to. The table tells you *which* site is worth looking at — grade,
occupancy, lining residues — and clicking takes you there. Sort by grade, then
by occupancy, and walk down the list.

### 6. A whole protein family

```bash
watcon family \
  --members  <path>/ptp_family15/members15.tsv \
  --alignment WatCon/data/examples/ptp_family/ptp_family_alignment.pir \
  --out-dir   family15
```

About 18 seconds for fifteen proteins and twenty-four structures:

```
Conservation: 344 columns, 229 covered by all 15 proteins, 58 of those unanimously conserved (every run grade >= 8)
...
Water sites: 612 clusters over 5913 waters, 556 occupied, 505 in two or more proteins, 1 in every protein, 117 lined by a unanimously conserved column
```

The per-structure table above that should show every structure superposing
between **0.65 and 0.97 Å**, all twenty-four.

Do not have the family data? Rebuild it from scratch — it fetches the
twenty-four entries itself:

```bash
python experiments/benchmark/scripts/build_family15.py
```

### 7. Fifteen proteins in one frame

```bash
pymol family15/watcon_family.pml
```

Twenty-four structures superposed, water sites coloured by how many proteins
hold them. `sites_shared` is the subset held by several.

**What to look for:** one single site is held by all fifteen — and it is a
buried structural water, not the catalytic one. That is a real result and a
deliberately humbling one; see the Limits section of the README.

---

## Part 2 — Full pass

### Every subcommand, in dependency order

**`fetch`** — including an entry the PDB no longer serves in PDB format:

```bash
watcon fetch --ids 1AAX 7GSA --out-dir raw
```

```
  1AAX   pdb (252 kB)
  7GSA   cif (547 kB)
Fetched 2/2 structures into .../raw
```

7GSA falls back to **mmCIF** because it has no PDB file at all — 31 of the 287
PTP1B entries are like this.

**`prepare`** — converts, finds the right chain, superposes, keeps waters:

```bash
watcon prepare --input-dir raw --out-dir prepared_ptp --reference 1AAX
```

```
Converted 1 file(s) to PDB format: 7GSA.cif
Preparing 2 structures (reference 1AAX.pdb chain A, 297 residues)
  1AAX     chain A  identity 1.00   233 waters  RMSD 0.00 A over 297 CA
  7GSA     chain A  identity 0.94   250 waters  RMSD 1.15 A over 283 CA
```

**`validate`** — check a ConSurf file *before* you rely on it:

```bash
watcon validate --consurf WatCon/data/consurf/fixtures/1AAX_A.grades.txt
```

```
  dialect              webserver
  method               bayesian
  msa_total            150
  records              321
  mapped               297
  unmapped             24
```

`records` counts rows in the grades table; `mapped` counts those ConSurf could
place on a structure. The difference is residues the crystal did not resolve.

With the annotated PDB, it cross-checks every grade against ConSurf's own:

```bash
watcon validate --consurf .../1AAX_A.grades.txt --pdb .../1AAX_A.consurf_ca.pdb
```

```
cross-check vs 1AAX_A.consurf_ca.pdb: checked=297 agree=297 mismatch=0 ...
```

**`view`** — the PyMOL session on its own:

```bash
watcon view --prepared watcon_demo/prepared \
            --consurf watcon_demo/consurf/1A2P_consurf_grades.txt \
            --out-dir view_test
```

```
  sites_all.pdb         190 occupied sites, grade in B-factor
  sites_conserved.pdb   57 sites lined by a grade>=8 residue
```

**`run`** — the input-file route, and `python WatCon/WatCon.py --input …`, which
still works exactly as upstream.

### Flags that change the science

| Flag | What it does |
|---|---|
| `family --no-sites` | pool conservation only, skip water sites — much faster |
| `family --tolerant` | warn instead of stopping when a structure disagrees |
| `family --state 3OLR=open` | label a structure, so per-state occupancy is reported |
| `view --site-radius` | how far from a site a residue counts as lining it |
| `--min-cluster-samples` | how many waters make a site; **changes what a site is** |
| `demo --keep` | reuse the output directory instead of clearing it |

### The trajectory and directed paths

Both use fixtures committed to the repository:

```python
from WatCon import generate_dynamic_networks as gdn
D = "WatCon/data/examples/ptp1b_ensemble"
CHEAP = {"density":"on","connected_components":"on","interaction_counts":"on",
         "per_residue_interactions":"on","characteristic_path_length":"off",
         "graph_entropy":"on","clustering_coefficient":"on","shortest_path":"off"}

_m, nets, centres = gdn.initialize_network(
    topology_file="ptp1b_ensemble_constant.pdb",
    trajectory_file="ptp1b_ensemble_constant.pdb",
    structure_directory=D, multi_model_pdb=True, water_name="HOH",
    msa_indexing=False, return_network=True, num_workers=1,
    cluster_coordinates=True, min_cluster_samples=2, analysis_conditions=CHEAP)
```

```
frames: 3
  frame 0: 9 waters, 14 edges
  frame 1: 9 waters,  9 edges
  frame 2: 9 waters, 12 edges
cluster centres across frames: 10
```

Swap in `md_active_site_h.pdb` with `include_hydrogens=True, angle_criteria=150,
max_distance=2.5` for the **directed** path, which needs real hydrogens:

```
oxygen-only 3.0 A          edges=11   undirected
directed 2.5 A angle 150   edges=13   directed
```

**`CHEAP` is not decoration.** `shortest_path` and
`characteristic_path_length` run all-pairs shortest paths and are on when
`analysis_conditions='all'`. On a 4,851-water system they took **819 s** and
**45 s** against **1.4 s** to build the network itself. Turn them off for
anything large; the network is identical.

### Six guards that must fire

A tool that only works when everything is right is not worth trusting. Each of
these is reproducible from repository fixtures, and each **must** fail.

**1. An unknown clustering method**

```python
from WatCon.find_conserved_networks import cluster_coordinates_only
cluster_coordinates_only(np.zeros((10, 3)), "kmeans", 3, None, 1)
```

```
ValueError: Unknown clustering method 'kmeans'. Choose one of: dbscan, hdbscan, optics.
```

**2. A crystal ensemble offered as a trajectory** — `ptp1b_ensemble_ragged.pdb`
with `multi_model_pdb=True`:

```
VaryingAtomCount: ... cannot be read as a trajectory: its models do not all
contain the same atoms (model 1 has 397, model 2 has 395; 3 distinct counts
across 3 models). ... Analyse it with the static path instead
```

**3. One protein paired with another protein's ConSurf run** — edit
`members15.tsv` so PTPN1 points at `1WCH_A.grades.txt`:

```
error: ConSurf data was supplied for 2F71 (PTPN1) but matched none of its
residues, so nothing is scored. The structure has chain(s) A and the ConSurf
run describes chain(s) A. Those agree, so the cause is not chain labelling:
either the residue numbering disagrees, or this run describes a different
protein. Check that the members file pairs each protein with its own ConSurf run.
```

Exit code 1. Without this guard the run completed and reported 114 unanimously
conserved columns instead of 117 — a wrong published number, silently.

**4. An active-site reference that matches nothing** —
`active_region_reference="resid 9999 and name SG"`:

```
ValueError: active_region_reference selected no atoms. Check the residue
numbering of this system -- a renumbered model will not answer to the numbering
of the structure the reference was written for.
```

**5. A ConSurf run whose chain matches nothing** — the MD fixture is chain `X`,
the run is chain `A`:

```
ConSurf data was supplied for md site but matched none of its residues ...
The structure has chain(s) X and the ConSurf run describes chain(s) A. Those
differ: pass consurf_chain_map to say which is which ...
```

Note the message differs from guard 3 — it suggests a chain map only when the
chains actually differ.

**6. A structure whose numbering disagrees with its run** — the identity check,
which is the whole point of the package:

```
ConSurf data does not describe 2F56.pdb: only 2 of 106 compared residues
agree (1.9%, threshold 95%). ... A3: structure has ASN, ConSurf has VAL
```

It is a *rate*, not pass/fail: a point mutant genuinely differs at one position
(99.1%) while a numbering offset collapses to 2.8%.

### Both environments

```bash
python -m pytest WatCon/tests -q            # from a checkout: 795 passed, 2 skipped
```

```bash
cd /somewhere/else
python -m pytest --pyargs WatCon.tests -q   # installed: 569 passed, 28 skipped
```

The installed run collects fewer modules and **says so**:

```
UserWarning: WatCon: not collecting 10 test module(s) -- test_consurf_crosscheck.py,
... Run the suite from a repository checkout to cover them.
```

Those ten need the raw ConSurf bundles and test structures a distribution
deliberately does not carry. Expected; silence would not be.

---

## Part 3 — Showing it to someone

Twelve minutes, in this order.

1. **`watcon demo`** — "the whole pipeline, six real crystal structures, ten
   seconds, offline. One of them has barnase as chain L, and it finds it."
2. **Open the PyMOL session** — "ConSurf's own colour scale. Yellow is *no
   data*, not low conservation — that distinction is ours, and most tools lose
   it."
3. **`enable sites`** — "each sphere is a water position several structures
   agree on. Size is how many."
4. **The plugin, and click a row** — "you cannot eyeball which of 190 sites
   matters. The table ranks them; clicking flies you there."
5. **`watcon family` on fifteen proteins** — "fifteen proteins, one ConSurf run
   each, one shared alignment. Eighteen seconds."
6. **The family session** — "twenty-four structures in one frame, under 1 Å."

Then the three things worth saying out loud:

- **The strongest result.** Sites lined by a residue *every* run calls conserved
  are held by 7.34 proteins on average against 4.08, p = 1.0 × 10⁻²⁴. Going from
  five runs to fifteen barely moved the unanimous fraction — 25.8% to 25.3% —
  so it was not a small-sample artefact.
- **The retraction.** At five proteins, 22 sites were held by all five including
  the catalytic water. At fifteen, exactly one is universal and it is a buried
  structural water. The earlier claim did not survive, and is corrected rather
  than dropped.
- **The limit, before anyone asks.** Burial is not disentangled. Buried residues
  are both more conserved and more likely to hold ordered water, so no causal
  claim is made. And conservation does **not** predict occupancy — ΔAUC ≈ 0 on
  two proteins. It re-ranks toward the active site; it is not a shortcut to it.

The honest framing: this is a tool that reports two measures side by side and
refuses to blend them, plus a set of checks that stop the run rather than
attach a score to the wrong residue.

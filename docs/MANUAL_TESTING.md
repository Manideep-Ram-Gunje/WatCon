# Testing WatCon-ConSurf by hand

Every command here was run before this page was written, and every expected
output is copied from what actually appeared. If something differs on your
machine, that is a real difference worth looking at — not a typo here.

**This page is self-contained.** You can paste any section into a chat with an
assistant that cannot see the repository and it will have enough to help you.

| Part | What | Time |
|---|---|---|
| [0](#part-0--from-a-machine-with-nothing-on-it) | From a machine with nothing on it | 15 min |
| [1](#part-1--smoke-test) | Smoke test — does it work at all | 10 min |
| [2](#part-2--pymol-and-the-interface) | PyMOL and the interface, in detail | 30 min |
| [3](#part-3--every-command) | Every command and every flag | 30 min |
| [4](#part-4--six-guards-that-must-fire) | Six guards that must fire | 15 min |
| [5](#part-5--when-something-goes-wrong) | Troubleshooting | — |
| [6](#part-6--showing-it-to-someone) | Demo script | 12 min |

Work in a scratch directory, never inside the repository:

```bash
mkdir watcon-testing
cd watcon-testing
```

---

## Part 0 — From a machine with nothing on it

The path a borrowed computer takes. Check each step before starting the next —
it is much easier to fix a problem where it happens.

### 0.1 Python 3.10 or newer

```bash
python --version
```

Expect `Python 3.10.x` or higher. If `python` is not found, try `python3`. On
Windows, install from [python.org](https://www.python.org/downloads/) and tick
**"Add Python to PATH"** during setup — almost every "command not found" problem
later traces back to that box.

### 0.2 Install the package

```bash
pip install "git+https://github.com/Manideep-Ram-Gunje/WatCon.git@consurf-integration"
```

This pulls in numpy, scipy, MDAnalysis, networkx, scikit-learn, biopython and
gemmi. Two to five minutes on a normal connection.

> **The `@consurf-integration` part is required.** This fork's `main` branch is
> a clean mirror of the original WatCon, which contains none of this work. Without
> the branch you get plain WatCon and `watcon demo` does not exist.

```bash
watcon --version
```

```
watcon-consurf 0.9.0+7.gbdce413
```

Installing from the branch puts you a few commits past the `v0.9.0` tag, so the
version carries a `+N.gHASH` suffix naming the commit. That is correct and
useful — it identifies exactly what you have. A bare `0.9.0` means you installed
from the tag. **`1+unknown` is the one to worry about**; see
[5.2](#52-watcon---version-prints-1unknown).

### 0.3 Prove it works

```bash
watcon demo
```

If this finishes, you have a working installation — whatever happens later.
Details in Part 1.

### 0.4 PyMOL — only for the interactive plugin

Everything else works without it. The analysis, the CSV outputs and the `.pml`
session files are all produced with no PyMOL installed; you only need PyMOL to
*open* those sessions.

```bash
conda install -c conda-forge pymol-open-source
```

No conda? On Windows the installer from [pymol.org](https://pymol.org/) also
works. On Linux, `apt install pymol` installs a PyMOL that often cannot see a
pip-installed WatCon — conda is the reliable route.

**Both of these must work in the same Python:**

```bash
python -c "import pymol; print(pymol.__file__)"
python -c "import WatCon; print(WatCon.__version__)"
```

If the first succeeds and the second fails, PyMOL has its own Python. See
[5.3](#53-pymol-opens-but-the-plugin-menu-is-missing).

### 0.5 Install the plugin

```bash
watcon plugin --install
```

```
Installed the PyMOL plugin into ...:
  .../pmg_tk/startup/watcon_consurf_plugin.py

Restart PyMOL, then:  Plugin  >  WatCon + ConSurf

The file is a three-line shim -- the plugin itself stays in the
installed package, so upgrading WatCon upgrades the plugin.
Remove it with:  watcon plugin --uninstall
```

Already installed? It says so and refuses rather than overwriting; add
`--force` if you mean to replace it.

Restart PyMOL. **Plugin → WatCon + ConSurf** should now exist.

### Checklist

| Step | Check | Needed for |
|---|---|---|
| Python ≥ 3.10 | `python --version` | everything |
| Package | `watcon --version` → `0.9.0...`, not `1+unknown` | everything |
| Demo | `watcon demo` completes | proves the install |
| PyMOL | `python -c "import pymol"` silent | viewing sessions |
| Plugin | menu entry exists | the interactive table |

---

## Part 1 — Smoke test

### 1.1 The bundled demo

```bash
watcon demo
```

Ten seconds, offline, on six real barnase crystal structures shipped inside the
package. The parts to check:

```
STEP 1/5  Prepare: find the barnase chain, superpose, keep waters
  1A2P     chain A  identity 1.00   148 waters  RMSD 0.00 A over 108 CA
  1BRN     chain L  identity 1.00   120 waters  RMSD 0.45 A over 108 CA
  1BRS     chain A  identity 1.00   174 waters  RMSD 0.46 A over 108 CA
  1BSA     chain A  identity 0.98    81 waters  RMSD 0.16 A over 107 CA
  1BSE     chain A  identity 0.99    59 waters  RMSD 0.24 A over 108 CA
  1RNB     chain A  identity 0.98    89 waters  RMSD 0.44 A over 108 CA
STEP 2/5  Stage ConSurf: one run describes every structure
STEP 3/5  Build water networks with conservation attached
STEP 4/5  Cluster recurring water sites and join to conservation
  193 cluster centres, 190 occupied, 165 carry conservation
  30 site(s) occupied in ALL 6 structures
  57 scored site(s) are lined by a highly conserved residue (ConSurf grade >= 8)
STEP 5/5  Write PyMOL projections
```

**Four things worth understanding, because they recur everywhere:**

**1BRN is chain L.** The demo set is deliberately awkward — in one entry the
barnase is not chain A, and 1BRS is half a complex. Finding the right chain is
part of the job.

**identity 0.98 is fine, 0.02 is not.** That column is how much of the structure
agrees with the ConSurf run, residue by residue. 1BSA and 1RNB are point mutants
and score 0.98. A numbering mismatch collapses to near zero, and the run stops.

**193 → 190 → 165 is not attrition.** 193 cluster centres; 190 have a water in
at least one structure; 165 of those have a protein residue near enough to carry
a conservation score.

**One ConSurf run covers all six.** Conservation is a property of the sequence.

### 1.2 What it wrote

```bash
ls watcon_demo
```

```
conservation.csv  consurf  preparation.csv  prepared  pymol  view
```

| File | What |
|---|---|
| `conservation.csv` | one row per site — occupancy and conservation as separate columns |
| `preparation.csv` | what happened to each structure, including any rejected and why |
| `prepared/` | the six structures in one frame, waters kept |
| `consurf/` | the grades file, staged once per structure |
| `view/watcon_view.pml` | the PyMOL session |

```bash
head -3 watcon_demo/conservation.csv
```

The columns are `cluster_id, occupancy, n_structures_occupied,
n_structures_total, occupancy_fraction, evo_min_score, evo_mean_score,
evo_max_grade, evo_n_residues, evo_n_low_confidence, evo_n_unscored`.

**Occupancy and conservation are separate columns and are never blended.**
Deciding how they relate is the science.

### 1.3 Open it

```bash
pymol watcon_demo/view/watcon_view.pml
```

Covered in detail in Part 2.

---

## Part 2 — PyMOL and the interface

### 2.1 What loads

```bash
pymol watcon_demo/view/watcon_view.pml
```

Five objects appear in the panel on the right:

| Object | What it is | On by default |
|---|---|---|
| `protein` | the reference structure, coloured by conservation | yes |
| `sites` | every occupied water site, as spheres | no |
| `sites_conserved` | only sites lined by a grade ≥ 8 residue, in red | yes |
| `site_contacts` | the residues lining those sites | no |
| `WatCon_contacts` | polar contacts between sites and lining residues | no |

Click a name to toggle it. `sites` is off at the start because 190 spheres at
once hides the protein.

### 2.2 The colour scale

The protein is coloured on **ConSurf's own 1–9 scale**:

- **maroon** — grade 9, most conserved
- **white/pale** — the middle grades
- **cyan** — grade 1, most variable
- **yellow** — **no ConSurf score for this residue**

**Yellow is not low conservation.** It means ConSurf did not score that
position — the residue is outside the analysed region, or is a modified residue
the run skipped. Reading yellow as "variable" is a real misreading, which is why
it is given its own colour rather than being blended into the scale.

### 2.3 Commands worth typing

At the PyMOL prompt (`PyMOL>`):

```
enable sites                 every occupied site; sphere size = how many structures hold it
enable WatCon_contacts       the lining residues and their polar contacts
enable site_contacts         just the lining residues
disable sites                hide them again

show labels, sites_conserved     put the site number on each red sphere
hide labels

zoom sites_conserved             frame the conserved sites
orient protein                   back to the whole protein

iterate sites_conserved, print(resi, b)    site number and its ConSurf grade
```

**Sphere size encodes occupancy** — a big sphere is a water position many
structures agree on. **B-factor carries the grade**, which is why `iterate`
prints it. A B-factor of **0 means no conservation data**, not grade zero.

### 2.4 No display? Check it headless

Useful over SSH, and as a scriptable test:

```bash
pymol -cq watcon_demo/view/watcon_view.pml -d "print(cmd.get_names('objects'))"
```

```
['protein', 'sites', 'sites_conserved', 'site_contacts', 'WatCon_contacts']
```

```bash
pymol -cq watcon_demo/view/watcon_view.pml -d "print(cmd.count_atoms('visible'))"
```

Expect around `1085`. If the object list is right, the session is fine —
`-c` means no GUI, `-q` quiet.

### 2.5 The plugin — One protein tab

**Plugin → WatCon + ConSurf.** Two tabs; this one first.

| Field | Put in | Notes |
|---|---|---|
| **Structures** | `watcon_demo/prepared` | a folder of already-prepared structures |
| **PDB ids** | *(leave empty)* | alternative to the above — type ids and it fetches them |
| **ConSurf** | `watcon_demo/consurf/1A2P_consurf_grades.txt` | one run describes every structure |
| **chain** | *(leave empty)* | only needed when the run and the structures label chains differently |
| **Reference** | *(leave empty)* | defaults sensibly; pick one to fix the frame |

Press **Run**. A progress bar moves; PyMOL stays usable throughout, because the
analysis runs on a background thread.

**If Run is greyed out or errors immediately**, the dialog is refusing input it
could not honour — the message says which field. That check is deliberate: it is
better than failing three minutes into an analysis.

### 2.6 The results table

Five columns:

| Column | What it means | Sort by it to find |
|---|---|---|
| **Site** | cluster id — a label, not a residue number | nothing; it is an identifier |
| **Grade** | highest ConSurf grade among lining residues (0 = no data) | the most evolutionarily constrained sites |
| **Structures** | how many structures have a water here | the most reproducible positions |
| **Waters** | how many water molecules went into the cluster | dense, well-defined sites |
| **Lining residues** | which residues are near it | a residue you care about |

Click a column header to sort; click again to reverse.

**The two orderings disagree, and that is the point.** Sort by Structures and
you get the most reproducible water positions — many of them surface packing
water. Sort by Grade and you get sites against conserved residues. Sites high in
*both* are the interesting ones.

### 2.7 Click a row

The camera flies to that site and shows the side chains lining it.

**This is the answer to "you cannot compare 3D structures by eye."** You are not
meant to. With 190 sites, no one can see which matters. The table ranks them and
clicking takes you there.

Try this: sort by **Grade** descending, click the top row, look at what lines it.
Then sort by **Structures** descending and do the same. In barnase the two lists
overlap only partly.

### 2.8 The plugin — Family tab

For several different proteins, each with its own ConSurf run.

The members table has one row per protein — **Protein**, **Structures folder**,
**ConSurf file**, **Reference** — plus three fields below:

| Field | Put in |
|---|---|
| **Alignment** | a PIR/FASTA alignment covering every structure |
| **Frame reference** | the structure everything is superposed onto |
| **States** | optional labels like `3OLR=open 2F71=closed` |

Its results table has six columns: **Site**, **Proteins**, **Waters**,
**States**, **Conserved columns**, **Residues by protein**.

**Proteins** is the one to sort by — it is how many different proteins hold a
water at that position, which is a much stronger statement than how many
structures do.

### 2.9 The family session

If you have the fifteen-protein data:

```bash
watcon family \
  --members <path>/ptp_family15/members15.tsv \
  --alignment <path>/WatCon/data/examples/ptp_family/ptp_family_alignment.pir \
  --out-dir family15
```

About 18 seconds:

```
Conservation: 344 columns, 229 covered by all 15 proteins, 58 of those unanimously conserved (every run grade >= 8)
Water sites: 612 clusters over 5913 waters, 556 occupied, 505 in two or more proteins, 1 in every protein, 117 lined by a unanimously conserved column
```

Every structure should superpose between **0.65 and 0.97 Å**.

Do not have the data? It fetches everything itself:

```bash
python experiments/benchmark/scripts/build_family15.py
```

Then:

```bash
pymol family15/watcon_family.pml
```

Twenty-four structures in one frame, plus `sites_all` and `sites_shared`.
Toggle individual proteins on and off in the object panel to see how the same
water position is held by different proteins.

**What to look for:** exactly one site is held by all fifteen — and it is a
buried structural water, not the catalytic one. That is a real result, and a
deliberately humbling one.

---

## Part 3 — Every command

### fetch

```bash
watcon fetch --ids 1AAX 7GSA --out-dir raw
```

```
  1AAX   pdb (252 kB)
  7GSA   cif (547 kB)
Fetched 2/2 structures into .../raw
```

7GSA arrives as **mmCIF** because it has no PDB file at all — 31 of the 287
PTP1B entries are like that. Needs internet; everything else here does not.

### prepare

```bash
watcon prepare --input-dir raw --out-dir prepared_ptp --reference 1AAX
```

```
Converted 1 file(s) to PDB format: 7GSA.cif
Preparing 2 structures (reference 1AAX.pdb chain A, 297 residues)
  1AAX     chain A  identity 1.00   233 waters  RMSD 0.00 A over 297 CA
  7GSA     chain A  identity 0.94   250 waters  RMSD 1.15 A over 283 CA
```

### validate

```bash
watcon validate --consurf <package>/WatCon/data/consurf/fixtures/1AAX_A.grades.txt
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

With `--pdb` it cross-checks every grade against ConSurf's own annotated PDB:

```
cross-check vs 1AAX_A.consurf_ca.pdb: checked=297 agree=297 mismatch=0 ...
```

### view

```bash
watcon view --prepared watcon_demo/prepared \
            --consurf watcon_demo/consurf/1A2P_consurf_grades.txt \
            --out-dir view_test
```

```
  sites_all.pdb         190 occupied sites, grade in B-factor
  sites_conserved.pdb   57 sites lined by a grade>=8 residue
```

### run, family, plugin, demo

`run` takes the original input-file interface — `watcon run --input input.txt`.
`python WatCon/WatCon.py --input input.txt` still works identically.
`family` is in [2.9](#29-the-family-session); `plugin` in
[0.5](#05-install-the-plugin); `demo` in [1.1](#11-the-bundled-demo).

### Flags that change the result

| Flag | Effect |
|---|---|
| `family --no-sites` | conservation only, skip water sites — much faster |
| `family --tolerant` | warn instead of stopping when a structure disagrees |
| `family --state 3OLR=open` | label a structure; per-state occupancy is then reported |
| `view --site-radius` | how far from a site a residue counts as lining it |
| `--min-cluster-samples` | how many waters make a site — **changes what a site is** |
| `demo --keep` | reuse the output directory instead of clearing it |

Try it on the demo and watch what happens:

```bash
watcon view --prepared watcon_demo/prepared             --consurf watcon_demo/consurf/1A2P_consurf_grades.txt             --out-dir v4 --min-cluster-samples 4
```

```
  sites_all.pdb          72 occupied sites, grade in B-factor      (was 190)
  sites_conserved.pdb    31 sites lined by a grade>=8 residue      (was 57)
```

**190 sites become 72.** That is correct, not a bug: `min_cluster_samples` is
how many waters must agree before a position counts as a site, so it *defines*
what a site is. Raising it demands more evidence per site. Any number you quote
from this tool is only meaningful alongside the setting that produced it.

### The trajectory and directed paths

Both use fixtures inside the package:

```python
from WatCon import generate_dynamic_networks as gdn
D = "<package>/WatCon/data/examples/ptp1b_ensemble"
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

Swap in `md_active_site_h.pdb` with `include_hydrogens=True,
angle_criteria=150, max_distance=2.5` for the **directed** path, which needs
real hydrogens:

```
oxygen-only 3.0 A          edges=11   undirected
directed 2.5 A angle 150   edges=13   directed
```

> **`CHEAP` is not decoration.** `shortest_path` and
> `characteristic_path_length` run all-pairs shortest paths and are on when
> `analysis_conditions='all'`. On a 4,851-water system they took **819 s** and
> **45 s** against **1.4 s** to build the network. Off, the network is identical.

---

## Part 4 — Six guards that must fire

A tool that only works when everything is right is not worth trusting. Each of
these **must** fail, and each is reproducible from files inside the package.

### 4.1 An unknown clustering method

```python
import numpy as np
from WatCon.find_conserved_networks import cluster_coordinates_only
cluster_coordinates_only(np.zeros((10, 3)), "kmeans", 3, None, 1)
```

```
ValueError: Unknown clustering method 'kmeans'. Choose one of: dbscan, hdbscan, optics.
```

### 4.2 A crystal ensemble offered as a trajectory

Run the Part 3 trajectory snippet against `ptp1b_ensemble_ragged.pdb`:

```
VaryingAtomCount: ... cannot be read as a trajectory: its models do not all
contain the same atoms (model 1 has 397, model 2 has 395; 3 distinct counts
across 3 models). ... Analyse it with the static path instead
```

Each crystal structure resolves a different number of waters, so no constant
atom count exists. The tool says so before starting rather than failing three
frames in.

### 4.3 A protein paired with another protein's ConSurf run

Copy `members15.tsv`, change PTPN1's grades file to `1WCH_A.grades.txt`
(PTPN13's), and run `watcon family` on it:

```
error: ConSurf data was supplied for 2F71 (PTPN1) but matched none of its
residues, so nothing is scored. The structure has chain(s) A and the ConSurf
run describes chain(s) A. Those agree, so the cause is not chain labelling:
either the residue numbering disagrees, or this run describes a different
protein. Check that the members file pairs each protein with its own ConSurf run.
```

Exit code 1. **Without this guard the run completed** and reported 114
unanimously conserved columns instead of 117 — a wrong published number, with no
warning. It was found by writing this page.

### 4.4 An active-site reference that matches nothing

Pass `active_region_reference="resid 9999 and name SG"`:

```
ValueError: active_region_reference selected no atoms. Check the residue
numbering of this system -- a renumbered model will not answer to the numbering
of the structure the reference was written for.
```

### 4.5 A ConSurf run whose chain matches nothing

The MD fixture is chain `X`; the 1AAX run is chain `A`:

```
ConSurf data was supplied for md site but matched none of its residues ...
The structure has chain(s) X and the ConSurf run describes chain(s) A. Those
differ: pass consurf_chain_map to say which is which ...
```

Note the difference from 4.3 — it suggests a chain map only when the chains
actually differ.

### 4.6 A structure whose numbering disagrees with its run

The check the whole package exists for:

```
ConSurf data does not describe 2F56.pdb: only 2 of 106 compared residues
agree (1.9%, threshold 95%). ... A3: structure has ASN, ConSurf has VAL
```

A **rate**, not pass/fail: a point mutant genuinely differs at one position
(99.1%) while a numbering offset collapses to 2.8%.

### The test suite

```bash
python -m pytest WatCon/tests -q          # from a checkout: 795 passed, 2 skipped
```

From an installed copy, in some other directory:

```bash
python -m pytest --pyargs WatCon.tests -q  # 569 passed, 28 skipped
```

```
UserWarning: WatCon: not collecting 10 test module(s) -- ...
Run the suite from a repository checkout to cover them.
```

Ten modules need raw data a distribution does not ship. Both numbers are
correct; the warning says which is which.

---

## Part 5 — When something goes wrong

### 5.1 `watcon: command not found`

The package installed but its script directory is not on PATH. Use it as a
module:

```bash
python -m WatCon.cli --version
```

If that works, PATH is the only problem. On Windows this is usually the
"Add Python to PATH" box at install time. `pip show -f watcon-consurf` shows
where things went.

### 5.2 `watcon --version` prints `1+unknown`

Built without git metadata — usually installed from a downloaded zip rather than
the git URL. The tool works; reinstall from the URL in [0.2](#02-install-the-package)
for a real version number.

### 5.3 PyMOL opens but the plugin menu is missing

Almost always PyMOL and WatCon in **different Pythons**. Inside PyMOL:

```
PyMOL> import WatCon; print(WatCon.__file__)
```

If that fails, install WatCon into PyMOL's Python:

```
PyMOL> import sys; print(sys.executable)
```

then, in a terminal, with that interpreter:

```bash
<that path> -m pip install "git+https://github.com/Manideep-Ram-Gunje/WatCon.git@consurf-integration"
<that path> -m watcon plugin --install
```

Restart PyMOL. If `import WatCon` works but the menu is still absent, run
`watcon plugin --install --force`.

### 5.4 The `.pml` opens but the protein is one flat colour

You loaded `sites_all.pdb` on its own instead of the `.pml`. Open
`watcon_view.pml` — the script loads what it colours.

### 5.5 "ConSurf data does not describe …"

Working as intended: your structure and your ConSurf run disagree about residue
numbering. The message names the first few disagreements. Either the run is for
a different protein, or the structure is renumbered. **Do not bypass this** — the
alternative is every conservation score landing on the wrong residue, invisibly.

### 5.6 The family run says a member "matched none of its residues"

That member's row in the members file points at the wrong grades file. See 4.3.

### 5.7 A run is taking many minutes

Almost certainly `shortest_path` and `characteristic_path_length`. Turn them off
in `analysis_conditions`; the network is unchanged. See the box in Part 3.

### 5.8 `watcon family` cannot find the family data

The fifteen-protein structures are not in the package — they are fetched:

```bash
python experiments/benchmark/scripts/build_family15.py
```

Needs internet and about a minute.

---

## Part 6 — Showing it to someone

Twelve minutes.

1. **`watcon demo`** — "the whole pipeline, six real crystal structures, ten
   seconds, offline. One of them has barnase as chain L, and it finds it."
2. **Open the session** — "ConSurf's own colour scale. Yellow means *no data*,
   not low conservation — most tools lose that distinction."
3. **`enable sites`** — "each sphere is a water position several structures
   agree on. Size is how many."
4. **The plugin, and click a row** — "you cannot eyeball which of 190 sites
   matters. The table ranks them; clicking flies you there."
5. **Sort by Grade, then by Structures** — "the two orderings disagree. That
   disagreement is the finding."
6. **`watcon family`** — "fifteen proteins, one ConSurf run each, one shared
   alignment. Eighteen seconds."
7. **The family session** — "twenty-four structures in one frame, all under 1 Å."

Then three things worth saying out loud:

**The strongest result.** Sites lined by a residue *every* run calls conserved
are held by 7.34 proteins on average against 4.08, p = 1.0 × 10⁻²⁴. Going from
five runs to fifteen barely moved the unanimous fraction — 25.8% to 25.3% — so
it was not a small-sample artefact.

**The retraction.** At five proteins, 22 sites were held by all five including
the catalytic water. At fifteen, exactly one is universal and it is a buried
structural water. The earlier claim did not survive, and is corrected rather
than dropped.

**The limit, before anyone asks.** Burial is not disentangled — buried residues
are both more conserved and more likely to hold ordered water, so no causal
claim is made. And conservation does **not** predict occupancy: ΔAUC ≈ 0 on two
proteins. It re-ranks toward the active site; it is not a shortcut to it.

The honest framing: a tool that reports two measures side by side, refuses to
blend them, and stops rather than attaching a score to the wrong residue.

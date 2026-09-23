# Code map — what lives where

The one page to open when you know *what* you want and not *where* it is.
34 modules, roughly 16,000 lines. This tells you which file, which function, and
what to read first.

Three documents pair with this one: [ARCHITECTURE.md](ARCHITECTURE.md) for how
the pieces fit together, [developer_guide.rst](developer_guide.rst) for how to
work on them, and [CONSURF_INTEGRATION.md](CONSURF_INTEGRATION.md) for why each
decision was made.

---

## The fifteen capabilities

Each entry: what it does, the file that owns it, and the function to read first.

### 1. Reading ConSurf grades files

| | |
|---|---|
| **Owns** | `WatCon/consurf/` — a subpackage of six modules |
| **Start at** | `consurf/parser.py` → `parse_consurf(source, strict=True)` |
| **Returns** | `ConSurfParseResult` — records, provenance, warnings |
| **Also** | `dialects.py` tells webserver from ConSurf-DB · `model.py` holds the dataclasses · `crosscheck.py` checks grades against ConSurf's own annotated PDB · `errors.py` the exception and warning types · `validate.py` backs `watcon validate` |

Depends on **nothing** — no MDAnalysis, no numpy. You can parse a grades file
with the scientific stack uninstalled, which is why these tests still run when
MDAnalysis is missing.

### 2. Residue identity — the join key

| | |
|---|---|
| **Owns** | `WatCon/residue_index.py` |
| **Start at** | `residues_from_pdb_file(path, chain=None)` and `atom_identity(atom)` |
| **Key idea** | every join is keyed on `(chain, resid, icode)`, never on array position |

Also here: `AMINO_ACID_3TO1` · `MODIFIED_RESIDUES` (CSP, PTR, MSE and two dozen
others) · `PROTONATION_VARIANTS` (Amber's HID/HIE/HIP, CHARMM's HSD/HSE/HSP) ·
`protein_selection()`, the MDAnalysis selection string that includes modified
residues, which a plain `protein` selection excludes.

**If two readers disagree about the chain, every lookup misses in silence.**
That happened — see Phase 28 in the changelog.

### 3. Attaching conservation to residues and waters

| | |
|---|---|
| **Owns** | `WatCon/evolutionary.py` — 1,263 lines, the largest module we wrote |
| **Start at** | `ConservationMap.build(parse_result)`, then `.for_residue(chain, resid, icode)` |
| **The guard** | `enforce_identity(coverage, strict=True)` — stops the run when the structure is not the protein the ConSurf run describes |
| **Also** | `describe_unmatched_coverage()` for the zero-coverage case · `conservation_of_clusters()` joins sites to conservation · `conservation_by_msa_column()` pools across a family |

### 4. Static water networks

| | |
|---|---|
| **Owns** | `WatCon/generate_static_networks.py` — **inherited**, roughly 2,000 lines |
| **Start at** | `initialize_network(structure_directory=..., ...)` |
| **Ours inside it** | the `protein_selection()` calls, the conformer filter, the conservation hand-off |

### 5. Dynamic and directed networks

| | |
|---|---|
| **Owns** | `WatCon/generate_dynamic_networks.py` — **inherited**, roughly 2,200 lines |
| **Start at** | `initialize_network(topology_file=..., trajectory_file=...)` |
| **Directed** | `include_hydrogens=True` with `angle_criteria` builds H→O edges instead of proximity edges |

**Two metrics dominate the cost.** `shortest_path` and
`characteristic_path_length` run all-pairs shortest paths: 819 s and 45 s
against 1.4 s to build the network itself, on a 4,851-water system. Switch them
off in `analysis_conditions` for anything large — the network is identical.

### 6. Clustering waters into recurring sites

| | |
|---|---|
| **Owns** | `WatCon/find_conserved_networks.py` — **inherited** |
| **Start at** | `cluster_coordinates_only(coords, cluster, min_samples, eps, n_jobs)` |
| **Methods** | `CLUSTERING_METHODS = ("hdbscan", "dbscan", "optics")` — anything else raises, naming the alternatives |

### 7. Superposition without MODELLER

| | |
|---|---|
| **Owns** | `WatCon/superpose.py` |
| **Start at** | `superpose_structures(mobile, reference)` → `SuperpositionResult` |
| **How** | Kabsch with reflection correction, over CA atoms paired by residue identity |

Upstream needed MODELLER for this — licensed, conda-only. For structures of one
protein this replaces it in numpy; MODELLER is still needed for MSA-based
alignment across a family.

### 8. Structure I/O and fetching

| | |
|---|---|
| **Owns** | `WatCon/structure_io.py`, `WatCon/fetch.py` |
| **Start at** | `convert_to_pdb(path, out)` · `fetch_structures(ids, out_dir)` |
| **Why** | RCSB no longer issues PDB files for large or recent entries — 31 of the 287 PTP1B entries have none. `fetch` tries PDB, falls back to mmCIF; `structure_io` converts with gemmi |
| **Also** | `require_constant_atom_count()` refuses a ragged multi-model PDB before the run starts, rather than three frames in |

### 9. Alternate conformers

| | |
|---|---|
| **Owns** | `WatCon/conformers.py` |
| **Start at** | `conformer_choices(atoms)` → the labels to keep |
| **Policy** | keep the most populated; where occupancies tie, keep every tied position |

### 10. Placing structures on an alignment

| | |
|---|---|
| **Owns** | `WatCon/alignment.py` |
| **Start at** | `map_structure_to_row(residues, row)` → `RowMapping` |
| **Then** | `consensus_columns(protein, mappings)` cross-checks structures of the same protein against each other |
| **Thresholds** | `MIN_COVERAGE = 0.80`, `DEFAULT_MIN_IDENTITY = 0.95` |

Mapping is **by sequence**, never by position — a published alignment row can
omit a residue or slide across a disordered loop, and both happen in the
authors' own PTP alignment.

### 11. Conservation pooled across a family

| | |
|---|---|
| **Owns** | `WatCon/family.py` |
| **Start at** | `build_family_conservation(proteins, alignment_path)` |
| **Input** | `read_members(path)` reads the tab-separated members file |
| **Output** | `FamilyConservation` — `.columns`, `.summary()`, `.conflicts()`, `.pairwise_spearman()` |

Pools **verdicts**, not scores. `unanimous_conserved` means every run
independently graded the column 8 or 9. ConSurf normalises within each run, so
averaging raw scores across runs would compare numbers that were never on a
common scale.

### 12. Water sites shared across a family

| | |
|---|---|
| **Owns** | `WatCon/family_sites.py` |
| **Start at** | `build_family_sites(proteins, family, reference, out_dir)` |
| **How** | iterative trimmed-core superposition through alignment columns — `TRIM_CUTOFF=2.0`, `MAX_CORE_RMSD=1.5`, `MIN_CORE_COLUMNS=50` — then cluster the pooled waters |

### 13. Scene — a result as data

| | |
|---|---|
| **Owns** | `WatCon/scene.py`, `WatCon/view.py`, `WatCon/family_scene.py` |
| **Start at** | `build_scene(...)` → a `Scene` of `Site` objects |
| **Why it matters** | the CLI and the plugin render the **same** `.pml`, so they cannot disagree about what a result looks like |

### 14. The PyMOL plugin

| | |
|---|---|
| **Owns** | `WatCon/pymol_plugin/` |
| **Start at** | `dialog.py` → `WatConDialog` · `__init__.py` → `__init_plugin__()` |
| **Threads** | `_Worker` and `_FamilyWorker` run the analysis off the UI thread, so PyMOL stays usable |
| **Install** | `watcon plugin --install` writes a shim into PyMOL's startup directory; the plugin itself stays in the package, so upgrading WatCon upgrades the plugin |

### 15. The command line

| | |
|---|---|
| **Owns** | `WatCon/cli.py`, `WatCon/demo.py` |
| **Start at** | `build_parser()`, then the `cmd_*` handler for the subcommand |
| **Subcommands** | `fetch` · `prepare` · `run` · `validate` · `view` · `family` · `plugin` · `demo` |

The legacy entry point `python WatCon/WatCon.py --input input.txt` still works
exactly as upstream.

---

## Reverse lookup — "I want to change X"

| If you want to… | Touch |
|---|---|
| change what counts as **highly conserved** | `scene.py` `HIGHLY_CONSERVED`; `family.py` uses grade ≥ 8 for `unanimous_conserved` |
| add a **CLI subcommand** | `cli.py` — write `cmd_yours`, register it in `build_parser`, set `func=`; then a dispatch test in `tests/test_cli.py` |
| support a new **ConSurf dialect** | `consurf/dialects.py` `detect_dialect()`, then a fixture under `data/consurf/fixtures/` |
| change **what counts as a site** | `--min-cluster-samples`, which reaches `find_conserved_networks.cluster_coordinates_only` |
| change **which residues line a site** | the site radius in `scene.py` and `family_sites.py` |
| add a **modified residue** | `residue_index.py` `MODIFIED_RESIDUES` — but a *protonation* variant belongs in `PROTONATION_VARIANTS` instead, because those are the same amino acid |
| change the **identity threshold** | `evolutionary.DEFAULT_IDENTITY_THRESHOLD`, `alignment.DEFAULT_MIN_IDENTITY` |
| change **what PyMOL draws** | `scene.py` `_commands_for()` — not the plugin, which only runs what scene emits |
| make a run **faster** | turn off `shortest_path` and `characteristic_path_length` in `analysis_conditions` |
| add a **family member** | the members TSV. No code change — every function here takes n members |

---

## Inherited from upstream, versus added here

This fork extends [WatCon](https://github.com/kamerlinlab/WatCon) by Brownless,
Harrison-Rawn and Kamerlin (*JACS Au* 2025). Knowing which is which tells you how
freely to edit: **inherited files carry the original authors' design and may go
back upstream one day; added files are ours.**

**Inherited**, edited only where noted:

| File | What we changed |
|---|---|
| `generate_static_networks.py` | `protein_selection()`, the conformer filter, the conservation hand-off, a guard on the directed path |
| `generate_dynamic_networks.py` | the same, plus the trajectory pre-flight check and `add_water` accepting hydrogen-free waters |
| `find_conserved_networks.py` | `NoWaterCoordinates`, `CLUSTERING_METHODS` validation, `eps=None` handling |
| `sequence_processing.py` | `seq_similarity` ported off the deprecated `Bio.pairwise2` |
| `residue_analysis.py` | `get_all_water_distances` made to raise rather than fail confusingly |
| `visualize_structures.py` | conservation-aware projections |
| `WatCon.py` | untouched — the legacy entry point |

**Added by this fork:** `consurf/` (six modules) · `evolutionary.py` ·
`residue_index.py` · `alignment.py` · `family.py` · `family_sites.py` ·
`family_scene.py` · `scene.py` · `view.py` · `superpose.py` · `prepare.py` ·
`structure_io.py` · `fetch.py` · `conformers.py` · `cli.py` · `demo.py` ·
`pymol_plugin/` · all 44 test modules.

---

## Where the data lives

| Path | What |
|---|---|
| `WatCon/data/examples/barnase/` | six real structures — the `watcon demo` dataset |
| `WatCon/data/examples/ptp_family/` | 24 CA extracts and the 24-row alignment |
| `WatCon/data/examples/ptp1b_ensemble/` | trajectory and hydrogen-bearing fixtures |
| `WatCon/data/consurf/fixtures/` | 15 ConSurf runs, with CA extracts of their annotated PDBs |
| `WatCon/data/consurf/exploratory/` | full result bundles, deliberately kept out of the wheel |
| `experiments/barnase_waters/`, `experiments/benchmark/` | the two studies: findings, method, numbers, scripts |

---

## Suggested reading order

**To understand the science** — about 30 minutes. The README, then
`experiments/barnase_waters/FINDINGS.md`, then the family section of
`experiments/benchmark/FINDINGS.md`.

**To understand the code** — about half a day:

1. [ARCHITECTURE.md](ARCHITECTURE.md) — the four invariants
2. `consurf/model.py` — the data shapes everything else moves around
3. `residue_index.py` — the join key, and why it is not an array index
4. `evolutionary.py`, class `ConservationMap` — where conservation meets structure
5. `scene.py` — how a result becomes something you can look at
6. `cli.py` — how a user reaches any of it

**To judge whether to trust it** — `tests/test_family_fifteen.py` and
`tests/test_chain_identity.py`. Both are written to explain the defect they
exist to prevent, so they read as documentation of the failure modes.

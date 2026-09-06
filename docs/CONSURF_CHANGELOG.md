# ConSurf Integration — Engineering Log

Append-only. Newest first. Six lines per entry.
Standing specification lives in [CONSURF_INTEGRATION.md](CONSURF_INTEGRATION.md).

Template:

```markdown
## YYYY-MM-DD — <summary>
- **Change:**   what was done
- **Files:**    paths touched
- **Reason:**   audit finding / spec rule
- **Tests:**    command -> result
- **Limits:**   known gaps left open
- **Decision:** architectural choice + rejected alternative
```

---

## 2026-09-06 — Phase 4: the consumer end -- cluster join, projections, report

- **Change:** Closed the gap that made the evolutionary layer a dead end. Before
  this, `find_conserved_networks.py` and `visualize_structures.py` had ZERO
  references to conservation -- it was attached per structure and never reached
  WatCon's cross-structure conserved water sites, was never drawn, never
  reported. Added `conservation_of_clusters` (the cross-structure join),
  `ClusterConservation`, `aggregate_site`, `write_conservation_report`,
  `project_clusters_by_conservation`, `pymol_project_evolutionary`, and a
  `conservation_report` post-analysis option. Extracted `water_residue_contacts`
  so the contact walk has one implementation instead of being duplicated in both
  builders.
- **Files:** `WatCon/evolutionary.py`, `WatCon/visualize_structures.py`,
  `WatCon/WatCon.py`, `WatCon/analysis.txt`,
  `WatCon/generate_{static,dynamic}_networks.py`;
  `WatCon/tests/test_evolutionary_clusters.py` (new),
  `WatCon/tests/test_visualisation_evo.py` (new); `docs/`.
  `find_conserved_networks.py` deliberately untouched -- the join consumes its
  output rather than modifying it, so structural conservation stays clean.
- **Reason:** The integration was unusable for its actual scientific purpose.
  Conservation existed on residues and waters but nothing joined it to the
  conserved water SITES that WatCon computes across a family, which is the level
  the original question is asked at.
- **Tests:** `python -m pytest .` from `WatCon/tests/` -> **355 passed, 1 failed**
  (up from 318; the failure is the pre-existing `test_inputs.py` missing
  `import sys`). New: 24 cluster-join, 13 projection. Verified on real
  structures: 7O7W self-consistency (every water used as its own centre is found
  there), residue keys are plain ints, far-away centres report empty.
- **Limits:** **The biology is not validated.** Whether conserved water sites are
  lined by conserved residues needs a real protein family with ConSurf for every
  member; we hold four datasets covering two proteins. Aggregation stays
  unweighted (distance weighting measured below ConSurf noise). No combined
  score. Static directed networks still blocked by the `max_neigbhbors` typo.
- **Decision:** Structural and evolutionary conservation are reported side by
  side and **never blended** -- separate CSV columns, separate projection files.
  Correlating them is the research question, and a tool that pre-computes one
  number would be answering it on the user's behalf. Rejected: a single
  "combined conservation" column, which was the obvious convenience and would
  have destroyed the ability to ask the question.
  Also fixed three bugs found on the way: `get_density` raised
  `UnboundLocalError` for any non-'all' selection (one line in each builder, not
  six methods as earlier notes said); `project_clusters` wrote the B-factor one
  column right of the PDB spec, so 0.25 read back as 0.2; and MDAnalysis numpy
  integers were leaking into residue keys.

## 2026-09-03 — Phase 3: ConSurf conservation attached to residues and waters

- **Change:** New `WatCon/evolutionary.py` joining parsed ConSurf results to
  WatCon residues on `(chain, resid, icode)`: `ResidueConservation`,
  `WaterConservation`, `ConservationMap` (with `coverage()`), `aggregate_water`,
  and supplied-file loading (`find_consurf_file`, `load_conservation`). Wired
  through both builders: `OtherAtom.evolutionary` (one shared immutable object
  per residue), `WaterMolecule.evolutionary`,
  `WaterNetwork.annotate_water_conservation()`, graph attributes `evo_score` /
  `evo_grade` / `evo_min_score` / `evo_n_residues`, three `initialize_network`
  kwargs, three `parse_inputs` keys, three appended CSV columns, and an
  `evolutionary_coverage` metric. Also fixed a falsy-zero bug in
  `get_per_residue_interactions` that silently dropped residue number 0.
- **Files:** `WatCon/evolutionary.py` (new), `WatCon/tests/test_evolutionary.py`
  (new); modified `WatCon/generate_static_networks.py`,
  `WatCon/generate_dynamic_networks.py`, `WatCon/residue_analysis.py`,
  `WatCon/WatCon.py`, `WatCon/input_static.txt`, `WatCon/input_dynamic.txt`;
  `docs/CONSURF_{INTEGRATION,CHANGELOG}.md`. `find_conserved_networks.py`
  deliberately **untouched** — structural conservation stays separate.
- **Reason:** Phase 3 of the ConSurf integration. Phases 1-2 produced a validated
  parser and a real residue identity; this is the join that makes evolutionary
  conservation usable alongside water-network structure.
- **Tests:** `python -m pytest .` from `WatCon/tests/` -> **318 passed, 1 failed**
  (the failure is the pre-existing `test_inputs.py` missing `import sys`).
  New: 37 evolutionary tests. Groups: REAL ConSurf 125, REAL identity+mapping 63,
  REAL evolutionary 37, synthetic 92, core regression passing. End-to-end on
  7O7W: 805/805 protein atoms scored across 236 residues, atoms of a residue
  share one object, 53/56 waters aggregated, 28 waters bridge >1 residue,
  coverage 236/236. CSV verified at 14 header fields == 14 row fields, `NA,NA,NA`
  when unscored.
- **Limits:** **No automatic ConSurf submission** — WatCon reads results you
  supply and never contacts the ConSurf server; the future hook is documented at
  `CONSURF_INTEGRATION.md` section 12, seam = `find_consurf_file`. No combined
  structural+evolutionary score (Phase 5). Aggregation is unweighted because the
  connection tuples carry no contact distance (Phase 4). Dynamic-pipeline
  conservation is wired but exercised only through the shared code path; no
  trajectory fixture exists here.
- **Decision:** Carry both score and grade with **score primary**. Measured on
  the same protein at two MSA depths, Spearman 0.955 (score) vs 0.947 (grade),
  and grade changed for 43% of positions — binning amplifies small score moves,
  so grade is the less stable representation. Rejected grade-only (loses
  resolution, less reproducible) and score-only (no human-readable bin, nothing
  to check against ConSurf's own B-factor output). Second decision: waters get
  several statistics rather than one number, because which summary is
  biologically right is still open and baking one in would force a re-run to
  change it. Third: low-confidence positions are flagged, never filtered —
  filtering moved cross-run agreement only 0.505 -> 0.495 and costs real data.

## 2026-09-03 — Phase 2: live residue -> MSA mapping wired in

- **Change:** Replaced `msa_indices[atm.resid - 1]` throughout the live path with
  identity lookup on `(chain, resid, icode)`. MODELLER import made lazy so the
  network modules import without it; `pdb_to_fastas` now delegates to
  `residue_index.residues_from_pdb_file` so the FASTA and the MSA index come
  from one walk; `OtherAtom` carries `chain`/`icode`; `ResidueIndex` built once
  per structure via `_build_residue_index` (raises on residue/MSA length
  mismatch); 7 duplicate graph-node lookups now read `molecule.msa_resid`;
  `convert_msa_to_individual` does exact column matching instead of
  `np.argmin` nearest-neighbour. Added `residue_index.atom_identity`.
- **Files:** `WatCon/sequence_processing.py`, `WatCon/generate_static_networks.py`,
  `WatCon/generate_dynamic_networks.py`, `WatCon/residue_index.py`;
  `WatCon/tests/test_core_mapping.py` (new), `WatCon/tests/inputs/msa/` (new
  fixtures), `WatCon/tests/{conftest,test_residue_index}.py`;
  `docs/CONSURF_{INTEGRATION,CHANGELOG}.md`.
- **Reason:** `msa_indices` is a POSITIONAL list indexed by absolute residue
  number. Measured on real structures: 1AKI 129/129 and P00648 157/157 correct
  by luck (dense, 1-origin, single-chain), but **7O7W 0/236 with 6 silent
  negative wrap-arounds** (`resid -5` -> `msa_indices[-6]`, which Python resolves
  by wrapping to the end of the list) and **1BRS 0/588** across its six chains.
- **Tests:** `python -m pytest WatCon/tests -q` from `WatCon/tests/` ->
  **281 passed, 1 failed**. The failure is pre-existing: `test_inputs.py:5` uses
  `sys.modules` without importing `sys` (tracked, commit `00d1a2c`). From the
  repo root it is 280/2 because `test_general.py` hardcodes the relative path
  `water_dir` and only passes when run from `WatCon/tests/` — also pre-existing
  (commit `e821a1b`). New: 21 core-mapping tests. Regression: `test_general.py`
  (end-to-end `initialize_network`, `msa_indexing=False`) passes.
  Installed: MDAnalysis 2.10.0, networkx 3.6.1, joblib 1.6.0, matplotlib 3.11.1,
  pandas 3.0.5, scikit-learn 1.9.0, biopython 1.88.
- **Limits:** `msa_with_modeller` and `perform_structure_alignment` remain
  **untested** — real MODELLER is licensed and not pip-installable (the PyPI
  `modeller` is an unrelated package). Everything else runs without it because
  `generate_msa_alignment` only calls MODELLER when the alignment file is
  absent. No ConSurf data is attached yet (Phase 3). Existing
  `msa_classification/*.csv` are invalid and must be regenerated.
- **Decision:** Kept the dynamic pipeline's "no MSA" fallback but made it honest:
  it used `u.residues.resids` (all residues, waters included) as if those were
  alignment columns. Now `fallback_to_resid=True` maps each residue to its own
  number, which is what the original intended. Rejected: silently returning None,
  which would have changed the meaning of `msa_resid` on that path. Separately,
  proved the old `pdb_to_fastas` crashes with `KeyError: 'RY '` on this repo's
  own `1AKI.pdb`, because its line filter accepted a `REMARK 500` record.

## 2026-09-03 — Residue-identity layer; parser follow-ups; conftest correction

- **Change:** Added `WatCon/residue_index.py` — explicit `(chain, resid, icode)`
  identity, a stdlib PDB residue walk, FASTA generation, a walk-vs-FASTA
  consistency check, and `ResidueIndex` mapping identity to MSA column. Its
  `build()` **raises** on a residue/MSA length mismatch or duplicate identity.
  Added `legacy_indexing_report` to quantify the old arithmetic per structure.
  Parser: bare-letter duplicate residue keys in RESIDUE VARIETY now warn and keep
  the first value instead of silently overwriting; added `POS < 1` and
  non-ascending-POS warnings. conftest: stopped ignoring `test_inputs.py` (it
  needs no MDAnalysis) and added a `pytest_report_header` line naming any module
  skipped for a missing dependency.
- **Files:** `WatCon/residue_index.py` (new); `WatCon/tests/test_residue_index.py`
  (new); `WatCon/consurf/{errors,parser}.py`; `WatCon/tests/conftest.py`;
  `docs/CONSURF_{INTEGRATION,CHANGELOG}.md`. **No core WatCon module modified.**
- **Reason:** Phase 2 identity layer. `msa_indices[atm.resid - 1]` is positional
  arithmetic on an absolute residue number; on 7O7W `resid = -5` becomes
  `msa_indices[-6]`, which Python resolves silently. The conftest entry for
  `test_inputs.py` was hiding a real, runnable, failing test.
- **Tests:** `python -m pytest WatCon/tests -q` -> **258 passed, 1 failed**.
  The failure is pre-existing and not ours: `test_inputs.py:5` uses `sys.modules`
  without importing `sys` (tracked, commit `00d1a2c`); left unfixed because core
  WatCon is out of scope this phase. Split: REAL-output consurf 125, REAL-output
  identity 41, synthetic 92. Legacy-arithmetic sweep: 1AKI 129/129 correct,
  P00648 157/157, 1BRS chain A 0/108, 1BRS all chains 0/588, **7O7W 0/225 with 6
  silent negative wrap-arounds**. Cross-layer identity agreement with ConSurf:
  7O7W 236/236, 1BRS_A 108/108, P00648 106 of 157 (51 outside ConSurf's range).
- **Limits:** The five live call sites are **not** rewired — that needs
  MDAnalysis, which is not installed here, so it would ship unverified.
  `residues_from_universe` is written but **unexecuted** for the same reason.
  `pdb_to_fastas` and `convert_msa_to_individual` are untouched. Modified
  residues whose alpha carbon is not named `CA` (7O7W's `PIA:68:A`) cannot be
  read by the walk; ConSurf omits them too, so the two stay consistent.
- **Decision:** Made the residue walk stdlib-only with MDAnalysis confined to one
  lazy adapter, so the entire identity layer is testable in this environment
  against real structures. Rejected: building directly on MDAnalysis, which would
  have made every assertion here unverifiable. Also: the walk keeps the **first**
  conformer of each residue whatever its altloc label — an earlier whitelist of
  `{"", "A", "1"}` silently dropped 11 residues of 7O7W that carry only `B`/`C`
  conformers, found by the cross-layer check against ConSurf.

## 2026-09-03 — Parser rebuilt for robustness; core WatCon untouched

- **Change:** Rewrote the ConSurf parser around a typed data model and two row
  grammars. Bare-letter RESIDUE VARIETY now parses as 100% (was silently
  dropped); confidence interval is optional so maximum-likelihood runs parse;
  ambiguity codes (`X`,`B`,`Z`,`J`,`O`) warn instead of aborting; ATOM
  identifiers support both dialects plus negative/zero numbers and insertion
  codes; grade-layer thresholds are captured; PDB-number discontinuities are
  reported; added `by_pdb_residue()` and three-valued `lookup_structural()`;
  added a B-factor cross-check against ConSurf's annotated PDB. Reorganised
  fixtures into `data/consurf/fixtures` (5 canonical files, 4 distinct datasets)
  and moved raw bundles to `data/consurf/exploratory`.
- **Files:** `WatCon/consurf/{errors,model,dialects,parser,crosscheck,validate,__init__}.py`;
  `WatCon/tests/{conftest,test_consurf_parser,test_consurf_dialects,test_consurf_adversarial,test_consurf_crosscheck}.py`;
  `WatCon/data/consurf/{fixtures,exploratory}/`; `docs/CONSURF_{INTEGRATION,CHANGELOG}.md`.
  `errors.py`, `model.py`, `dialects.py`, `crosscheck.py` are new; `run1|run2|run3/`
  and `additional/` were superseded (byte-identical content preserved under
  `fixtures/` and `exploratory/`). **No core WatCon module was modified.**
- **Reason:** Audit findings 1–9. The two that lost data without erroring were
  bare-letter variety (117/1152 records) and the absence of any signal that
  ConSurf omits non-standard residues, which makes POS a non-dense index.
  Research into official ConSurf docs added three more: the confidence interval
  is Bayesian-only, a second grades dialect exists, and ConSurf also runs on
  nucleotides.
- **Tests:** `python -m pytest WatCon/tests -q` -> **217 passed**
  (parser 108, adversarial 58, dialects 34, cross-check 17).
  Validation sweep over all 10 real grades files -> **1508 records, 0 failures,
  0 unparsed-variety warnings, 173 fully-conserved records retained**.
  B-factor cross-check -> **794/794 residues agree** across 3 bundles.
- **Limits:** Maximum-likelihood runs, the ConSurf-DB dialect, nucleotide runs,
  insertion codes and unusual chain IDs are implemented from documentation and
  covered only by **synthetic** tests — no real file of those kinds is held.
  `test_general.py`/`test_static.py` still fail to import without MDAnalysis
  (pre-existing); they are now skipped at collection so the ConSurf suite runs.
- **Decision:** "Structurally present but unscored" is a *query result*
  (`lookup_structural` -> `SCORED` / `NO_SCORE_IN_FILE` / `NOT_IN_FILE`), not a
  record field — the residue has no row at all, so a field would have been
  another overload. Rejected: adding an `unscored` flag to `ConSurfRecord`.
  Also decided an unparseable ATOM field is a hard error rather than a warning,
  because leaving it as `None` made it indistinguishable from a legitimately
  unmapped residue — silent corruption of the join key.

---

## 2026-09-06 — Phase 5a: clear the defect backlog

- **What:** Fixed three real defects carried since the audit, and made one test
  cwd-independent. The suite is now fully green for the first time.
- **Files:**
  - `WatCon/generate_static_networks.py:1442` — `max_neighbors=max_neigbhbors`
    -> `max_neighbors=max_neighbors`. `max_neigbhbors` was never defined, so
    `extract_objects(directed=True)` raised `NameError` on every call. Static
    **directed** networks had therefore never executed in any released version.
  - `WatCon/generate_static_networks.py:617` — `np.linalg.norm(water1)` ->
    `np.linalg.norm(prot_prot)`. In that branch the numerator is
    `dot(prot_water, prot_prot)` and `water1` is undefined; the cosine was being
    normalised by the wrong (nonexistent) vector. **:655 and :702 use `water1`
    correctly** — it is assigned two lines above each — and were left alone.
  - `WatCon/tests/test_inputs.py` — added the missing `import sys`, which made
    the module's single test raise `NameError`.
  - `WatCon/tests/test_general.py` — `initialize_network('water_dir', ...)`
    resolved the path against the cwd, so the test passed only when pytest ran
    from `WatCon/tests/`. Now resolved relative to `__file__`.
- **Verified:** The directed static path was *executed*, not just compiled:
  `initialize_network('WatCon/tests/water_dir', network_type='water-water',
  include_hydrogens=True)` builds a directed graph (8 nodes, 6 edges,
  `is_directed() == True`), against 8 nodes / 6 edges undirected.
- **Tests:** `python -m pytest WatCon/tests -q` -> **368 passed, 0 failed**
  (was 367 passed / 1 failed). Confirmed the `test_general` failure predated
  these edits by re-running it on a stashed tree.
- **Corrections to the plan, found while checking rather than assuming:**
  - `sequence_processing.py` was reported to define `perform_structure_alignment`
    and `msa_with_modeller` **twice**. It does not. Lines 563-656 sit inside a
    `'''...'''` block — commented-out old code that `grep '^def '` matched as if
    it were live. `ast.parse` confirms each function is defined exactly once, and
    the live `perform_structure_alignment` is the good one that returns the
    rotation/translation matrices `align_with_waters` consumes. **No change made.**
  - `residue_analysis.get_all_water_distances` really is broken: at line 124 it
    unpacks two values from `get_per_residue_interactions`, which returns one
    (`residue_dict`), and then indexes `interaction_data[selection]['Water-Protein'][1]`
    — a structure that function has never produced. It is written against an API
    that no longer exists and is called from nowhere. Repairing it would mean
    inventing a specification, so it is **left untouched and recorded here** as
    known-dead rather than silently "fixed" into something untested.
- **Limits:** The directed path is now reachable and produces a directed graph on
  a small water-only box. That is one execution, not validation of its
  hydrogen-bond geometry; the angle criteria in that path remain unexercised by
  any test.

---

## 2026-09-06 - Phase 5b: verify the ConSurf file describes the structure

- **What:** The join key is `(chain, resid, icode)`. It encodes nothing about
  *which amino acid* a score belongs to, so a ConSurf run numbered differently
  from the structure attaches a score to every residue, reports 100% coverage,
  and is wrong everywhere. Nothing downstream can see that. This adds the only
  check that can: compare the residue ConSurf recorded against the residue the
  structure actually contains.
- **Why now:** With two hand-checked structures it did not matter. The barnase
  study applies **one** ConSurf run to **21 depositions** by different groups,
  where numbering conventions vary. The check is what makes that trustworthy.
- **Files:**
  - `WatCon/evolutionary.py` - `ResidueConservation` gains `amino_acid` (SEQ,
    one-letter) and `residue_name` (ATOM, three-letter). `from_record` was
    dropping both; the data was parsed and discarded. `CoverageReport` gains
    `identity_matched`, `identity_mismatched`, `mismatches`, `identity_rate`,
    `describe_identity()`. `ConservationMap._compare_identity` and the check in
    `coverage()`. New `enforce_identity()` and `DEFAULT_IDENTITY_THRESHOLD`.
  - `WatCon/generate_static_networks.py`, `WatCon/generate_dynamic_networks.py`
    - call `enforce_identity` before the metric is recorded, and report
    `identity_matched` / `identity_mismatched` / `identity_rate`.
  - `WatCon/tests/test_identity_check.py` - new, 14 tests.
- **A rate, not a boolean.** Two different things cause mismatches. A point
  mutant genuinely differs from the sequence ConSurf aligned - correct data
  about a real difference. A numbering offset differs almost everywhere. The
  threshold (0.95) was chosen from measurement on real barnase, not taste:

  | case | identity rate |
  |---|---|
  | correct numbering | 100% (108/108) |
  | one point mutation | 99.1% (107/108) |
  | numbering shifted by one residue | **2.8%** (3/107) |

  0.95 sits in a wide empty gap between those regimes.
- **Degrades rather than lies.** Objects carrying no residue name (bare
  coordinates) report `identity_rate is None` and "not checked" - never "all
  mismatched". Identity is counted **per residue**, not per atom, so a large
  residue cannot outvote a small one.
- **Tests:** `python -m pytest WatCon/tests -q` -> **382 passed, 0 failed**
  (was 368). Verified on the real pipeline: 1BRS reads 18.9% coverage with 100%
  identity - ConSurf covered one of six chains, and was correct about that one.
  Coverage and identity are orthogonal and are now reported separately.
- **Limits:** The check compares only residues ConSurf **matched**. It cannot
  detect a run that covers the wrong region entirely while agreeing on the
  residues it does touch; low coverage remains the signal for that, and is
  reported but still not enforced. The 0.95 threshold is calibrated on one
  protein - a structure with many genuine mutations could trip it, which is why
  it is a parameter and the error lists the disagreeing residues.

---

## 2026-09-06 - Phase 5c: superposition without MODELLER

- **What:** New `WatCon/superpose.py`. Kabsch rigid-body superposition in numpy,
  producing exactly the `{'Rot': [...], 'Trans': [...]}` structure
  `perform_structure_alignment` returns, so it drops straight into
  `align_with_waters` - which already moves waters along with their protein.
- **Why:** WatCon could only superpose through MODELLER, which is licensed and
  cannot be pip-installed. The barnase study needs 21 structures in one frame.
  Nothing downstream of the transforms needed changing; only their source did.
- **Files:** `WatCon/superpose.py` (new), `WatCon/tests/test_superpose.py` (new,
  18 tests). **No existing module was modified** - this adds an alternative
  source of transforms rather than replacing the MODELLER path.
- **Verified on real data:** the three crystallographically independent barnase
  copies in 1BRS (chains A, B, C).

  | copy | RMSD before | RMSD after | CA matched |
  |---|---|---|---|
  | B onto A | 39.15 A | **0.31 A** | 108 |
  | C onto A | 35.67 A | **0.26 A** | 108 |

- **Two bugs found by testing rather than reasoning:**
  1. **The correspondence was keyed on the chain letter.** Copies A, B and C
     therefore shared *zero* residues and looked un-superposable. A chain letter
     is a label the depositor chose, not part of a residue's identity across
     structures - the same protein is chain A in one entry and chain B in
     another. Now keyed on `(resid, icode)`, with an explicit error if more than
     one chain is in play. Confirmed the regression tests catch it: reverting
     the fix fails 8 of 18.
  2. **`align_with_waters` always treats the sorted-first file as the
     reference**, indexing `rotation_matrices[i - 1]`. A different reference
     would pair every structure with the wrong transform and still produce
     something that looked like a superposition. `as_dict()` now refuses when
     the reference does not sort first, and says what to do instead.
- **Reflection correction is tested by name.** Without it the SVD returns an
  improper rotation (determinant -1) for a mirrored target, superposing a
  protein onto its own mirror image - not a rigid motion, and silently the wrong
  handedness.
- **Tests:** `python -m pytest WatCon/tests -q` -> **400 passed, 0 failed**
  (was 382).
- **Limits:** Correspondence is by **residue number**, so this handles one
  protein - different crystal forms, mutants, ligand complexes. It is not a
  sequence aligner and will not superpose homologues numbered differently; the
  MSA path remains the right tool for those. Mutated positions still contribute
  their CA, deliberately: a side-chain substitution barely moves the backbone,
  and excluding them would bias the fit toward unmutated regions.

---

## 2026-09-06 - Phase 5d: family scaffold, and a silent wrong-file bug

- **What:** Pooling conservation from *separate* ConSurf runs onto shared MSA
  columns, so extending to a real protein family means supplying ConSurf files
  rather than writing code. Plus a genuine bug found while validating it.
- **Files:** `WatCon/evolutionary.py` - `ColumnConservation`,
  `conservation_by_msa_column`, `family_summary`, `CROSS_RUN_AGREEMENT`, and a
  fix to `find_consurf_file`. `WatCon/tests/test_family_scaffold.py` (new, 14
  tests); two regression tests added to `test_evolutionary.py`.

### The bug: a longer name beaten by its own prefix

`find_consurf_file("P00648_50", ...)` returned **`P00648_150.grades.txt`**.

The lookup tried the loose token (everything before the first underscore)
*before* the full stem. The bare prefix `P00648` matched both files, and the
alphabetically first one won. The structure silently received a different run's
conservation, and nothing downstream could detect it - coverage stays at 100%,
identity stays at 100% (same protein), every score looks plausible.

Found because two "independent" ConSurf runs pooled to *identical* values at all
107 positions, which is not what two runs at different MSA depths should do.

Fixed by trying the **most specific token first**, and **refusing** when a token
matches more than one file rather than picking arbitrarily.

### Measured cross-run agreement

With the lookup fixed, the two real barnase runs are actually distinguishable,
and confirm the figure quoted in the scaffold:

| pair | score rho | grade rho | grade changed |
|---|---|---|---|
| P00648 depth 150 vs depth 50 | **0.955** | **0.947** | **42%** |
| 1BRS run vs P00648 depth 150 | 0.380 | 0.365 | 62% |
| 1BRS run vs P00648 depth 50 | 0.367 | 0.360 | 66% |

The second and third rows are the more sobering ones: the *same protein*, run
from a different starting structure, agrees at only rho ~= 0.37. Recorded in
`CROSS_RUN_AGREEMENT` as the honest floor - any family-level effect smaller than
this sits inside ConSurf's own reproducibility noise.

### Design consequences

- **Grades are never averaged.** A grade is an ordinal per-run percentile bin;
  its mean is not a grade. `ColumnConservation` exposes `max_grade`, not a mean.
- **`spread` is reported alongside every pooled value**, because scores are
  z-normalised *within* each run and are only relatively comparable.
- **`unanimous_conserved`** (every member independently called grade >= 8) is the
  one family-level claim that never compares values across runs, only each run's
  own verdict. It is the recommended statistic.
- **Pooling without an MSA is refused.** Falling back to residue numbers would
  align position 40 of one protein to position 40 of another.

- **Tests:** `python -m pytest WatCon/tests -q` -> **416 passed, 0 failed**
  (was 400).
- **Limits:** The scaffold is exercised on two runs of *one* sequence, which is
  the easiest possible case and still only 0.955 agreement. It has **not** been
  run on a real multi-protein family, because we hold no ConSurf data for one.
  The join and the arithmetic are tested; the biology is not.

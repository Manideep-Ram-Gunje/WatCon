# ConSurf Integration — Specification

Standing reference for `WatCon/consurf`. Edited in place as understanding changes.
For *what changed and when*, see [CONSURF_CHANGELOG.md](CONSURF_CHANGELOG.md).

**Status:** Phases 1-4 complete. The evolutionary layer is now consumed end to
end: conservation reaches WatCon's cross-structure conserved water sites, is
projected for viewing, and is reported as a table.

**Earlier status (kept for history):** Phases 1-3 complete. ConSurf results supplied by the user are parsed,
joined to WatCon residues by `(chain, resid, icode)`, and aggregated onto waters.
No combined score exists yet (Phase 5). WatCon does **not** contact the ConSurf
server -- see section 11.

---

## 1. Layout

```
WatCon/residue_index.py    residue identity + MSA-column mapping (stdlib core,
                           MDAnalysis only in one lazy adapter)
WatCon/evolutionary.py     ConSurf conservation joined to WatCon residues/waters
                           (stdlib + consurf + residue_index; no MDAnalysis)

WatCon/consurf/            source — standard library only, no WatCon imports
  errors.py                ConSurfParseError, WarningCode, ConSurfWarning, MalformedLine
  model.py                 dataclasses (below)
  dialects.py              row grammars, symbol tables, ATOM identifier parsing
  parser.py                parse_consurf()
  crosscheck.py            validation against ConSurf's annotated PDB
  validate.py              CLI:  python -m WatCon.consurf.validate

WatCon/tests/
  conftest.py                    fixture paths, bundle reader, optional-dep skips
  test_consurf_parser.py         REAL-output validation
  test_consurf_dialects.py       SYNTHETIC, specification-derived
  test_consurf_adversarial.py    SYNTHETIC, malformed input
  test_consurf_crosscheck.py     REAL-output, B-factor oracle
  test_residue_index.py          REAL-output, residue identity + MSA mapping
  test_core_mapping.py           REAL-output, live residue -> MSA in the builders
  test_evolutionary.py           REAL-output, ConSurf -> residues/waters/graph

WatCon/data/consurf/
  fixtures/                5 canonical grades files (4 distinct datasets)
  exploratory/             full ConSurf bundles (.tar.gz) and MSAs
```

The `consurf` package must stay dependency-free. Core WatCon needs MDAnalysis and
modeller; neither is installed in every environment, so this boundary is the only
reason the parser is testable on its own.

## 2. Fixtures

| File | Dataset | Notable |
|---|---|---|
| `1BRS_A_150.grades.txt` | 1BRS chain A, depth 150 | baseline; 2 leading unmapped |
| `1BRS_A_150.crlf.grades.txt` | same content, CRLF | line-ending twin |
| `P00648_150.grades.txt` | P00648, depth 150 | UniProt offset +47; **interior** unmapped POS 51 |
| `P00648_50.grades.txt` | P00648, depth 50 | shallow MSA; AlphaFold model |
| `7O7W_A_45.grades.txt` | 7O7W chain A, depth 45 | **negative/zero resnums, POS gap, `X` in variety, 41 fully conserved** |

4 distinct datasets, 573 records. The `exploratory/` bundles contain the same
grades files plus the annotated PDBs used by the cross-check.

## 3. Format contract

Confidence tiers: **[G]** guaranteed by ConSurf docs · **[O]** observed in our
files · **[T]** tolerated conservatively, unconfirmed.

| Field | Tier | Contract |
|---|---|---|
| POS | [G] | SEQRES-derived position. **Not guaranteed dense** — ConSurf may omit non-standard residues. Never assume `POS == pdb_number`. |
| SEQ | [G] | One-letter code; amino acid **or nucleotide** (A/C/G/T; U→T). |
| ATOM / 3LATOM | [G] | `NAME:NUM:CHAIN` (webserver) or `NAMENUM:CHAIN` (ConSurf-DB), or `-` / `___N:C` when unmapped. |
| SCORE | [G] | Normalised rate, mean 0 sd 1. **Negative = conserved.** Per-run normalised → **not comparable across runs**. |
| COLOR | [G] | Grade 1–9, **9 = conserved**, optional trailing `*`. |
| CONFIDENCE INTERVAL | [G] | **Bayesian only — absent under maximum likelihood.** 25th/75th percentiles plus both bound grades. |
| B/E, F/S | [G] | `b`/`e`, `f`/`s`. Absent entirely in the ConSurf-DB dialect; blank when unassigned. |
| MSA DATA | [G] | `present/total`; `total` constant per file. |
| RESIDUE VARIETY | [G] | Composition of the MSA column. |

Detail rules:

- **Bare-letter variety means 100%** [O→G]. Proven against `msa_aa_variety_percentage.csv`.
- **`<1%`** is a real sub-1% value, distinct from unknown [G].
- **Percentages are truncated and need not sum to 100** [O]. 134/249 positions in
  7O7W sum below 99%. Never renormalise silently.
- **Residue numbers may be negative or zero** [O] — expression tags (`HIS:-5:A` … `PRO:0:A`).
- **Insertion codes** [T] — round-trip as number + code; not present in our files.
- **Chain IDs are opaque strings** [G/T] — never assume length 1 or uppercase.
- **Non-standard residues** may be omitted from the table entirely [G] *or* appear
  as `X` in the variety [O]. Both occur in 7O7W.
- **Grade-layer thresholds differ per run** [G] — 3 distinct sets across our 4 datasets.

### Two dialects

| | webserver (ours) | consurfdb / gradesPE |
|---|---|---|
| Fields | 10 tab-separated | 9 |
| Identifier | `THR:-2:A` | `MET1:A` |
| Confidence | one merged field | two fields |
| B/E, F/S | present | absent |
| Missing density | `-` | `___1:A` |

Both are read by the same whitespace-tolerant grammar — the merged-vs-split
distinction vanishes once a row is tokenised. Only the identifier form and the
optional columns genuinely differ.

## 4. Data model

```
PdbResidueId(chain, seq, icode, name3, raw)     .key -> (chain, seq, icode)
ResidueFrequency(residue, percent, less_than_one_percent, is_full)
Confidence(score_lower, score_upper, grade_at_lower, grade_at_upper)
MsaSupport(present, total)                      .fraction
GradeLayer(grade, score_from, score_to)
ConSurfProvenance(source, dialect, method, alphabet, msa_total, line_ending)
ConSurfRecord(...)                              .is_fully_conserved, .warning_codes
ConSurfParseResult(...)                         .by_position() .by_pdb_residue()
                                                .mapped_records() .unmapped_records()
                                                .unverified_records() .covered_range()
                                                .lookup_structural()
```

Four overloads deliberately removed:

| Was | Problem | Now |
|---|---|---|
| `residue_variety == {}` | "100% conserved" vs "no data" | real entry with `is_full=True` |
| `pdb_residue_number: str` | number + icode in one string | `PdbResidueId.seq` / `.icode` |
| `confidence_grade_lower/upper` | names imply `lo <= hi`; grades invert scores | `grade_at_lower` / `grade_at_upper` |
| `atom is None` | "no mapping" vs "unparseable mapping" | `pdb_residue is None` vs `mapping_unverified` |

### The structural query direction

WatCon starts from a structure and asks for conservation. The grades file is
organised the other way, so `lookup_structural(chain, seq, icode)` returns three
states rather than a bare `None`:

- `SCORED` — ConSurf scored it.
- `NO_SCORE_IN_FILE` — inside the covered range but absent, e.g. 7O7W's chromophore `PIA:68:A`.
- `NOT_IN_FILE` — outside the covered range, e.g. chain B of 1BRS.

`position_gaps` carries the parser-side evidence (7O7W: `[(65, 69)]`).

## 5. Error handling

**Principle: never guess a residue mapping.** A mapping that cannot be verified
is an error, not a default.

Hard errors — missing file; no header; zero records; duplicate POS; duplicate
`(chain, seq, icode)`; inconsistent MSA counts; grade outside 1–9; inverted
confidence interval; **unparseable ATOM field that is not `-`**.

Warnings — unusual SEQ symbol (`X`,`B`,`Z`,`J`,`O`); low confidence `*`; shallow
MSA; wide interval; variety sum below 90%; PDB-number discontinuity;
non-positive residue number; unparsed variety token; mixed line endings; missing
grade layers; no confidence interval (ML run).

Tolerated silently — blank B/E and F/S; absent CI; LF vs CRLF; either dialect;
percent truncation; whitespace- or tab-delimited rows.

| | `strict=True` | `strict=False` |
|---|---|---|
| Hard error | raises `ConSurfParseError` | row → `malformed_lines`; file-level → warning |
| Return | result, or raise | always a result |

A record with `mapping_unverified` is **excluded from `by_pdb_residue()`**, so a
hit in the join table is always a mapping ConSurf actually stated.

## 6. Cross-validation

`crosscheck.check_against_annotated_pdb()` compares every mapped record's grade
against the B-factor column of `*_ATOMS_section_With_ConSurf.pdb`, which ConSurf
writes itself and leaves blank for unscored residues.

**Do not use `msa_aa_variety_percentage.csv` as an oracle.** Its `pos` and
`ConSurf grade` columns are indexed on grades POS (249/249 agreement), but its
amino-acid composition columns are shifted by one row from the point where
ConSurf omitted a residue (composition agrees at offset +1: 172/248 vs 78/249).
The grades file is the internally consistent artifact — its SEQ residue appears
in its own variety for 249/249 records.

## 7. Known limitations

Verified against real ConSurf output: webserver dialect, Bayesian method,
protein alphabet, LF and CRLF, negative/zero residue numbers, interior unmapped
residues, omitted non-standard residues, `X` in variety, MSA depths 45/50/150.

**Supported but NOT verified against real output** — implemented from
documentation, covered only by synthetic tests:

- maximum-likelihood runs (no confidence interval)
- the ConSurf-DB / `gradesPE` dialect
- nucleotide (DNA/RNA) runs
- insertion codes and non-single-letter chain identifiers

Obtaining one real ML-mode file and one real ConSurf-DB file would close this gap.

## 8. Residue identity — `WatCon/residue_index.py`

Replaces `msa_indices[atm.resid - 1]` with a validated lookup:

```
(chain, resid, icode) -> ordinal (sequence order) -> msa_indices[ordinal]
```

| Symbol | Role |
|---|---|
| `StructureResidue` | one residue; `resid` may be negative or zero |
| `residues_from_pdb_file` | stdlib PDB walk — the reference ordering |
| `residues_from_universe` | MDAnalysis adapter, lazily imported |
| `build_fasta_sequence` | sequence from a residue walk |
| `check_sequence_consistency` | walk vs an existing FASTA — catches drift |
| `ResidueIndex.build` | **raises on length mismatch or duplicate identity** |
| `ResidueIndex.msa_column` | lookup; `None` when unknown, never a wrong value |
| `legacy_indexing_report` | quantifies what the old arithmetic does per structure |

### Why the old arithmetic survived

Measured with `legacy_indexing_report` on real structures:

| Structure | Residues | Chains | Legacy correct | Silent wrap-around |
|---|---|---|---|---|
| 1AKI | 129 | A | **129/129** | 0 |
| P00648 (AlphaFold model) | 157 | A | **157/157** | 0 |
| 1BRS chain A | 108 | A | 0/108 | 0 |
| 1BRS all chains | 588 | A–F | 0/588 | 0 |
| **7O7W** | 225→236 | A | **0/225** | **6** |

Dense 1-origin single-chain structures — 1AKI and the AlphaFold model — give the
right answer, which is why nothing ever failed. Everything else is silently wrong.
On 7O7W, `HIS:-5:A` resolved to column 220 instead of 1.

Note the two coordinate systems must not be conflated: P00648's notorious +47
offset is between **ConSurf POS and the PDB residue number**, not inside the PDB.
The structure itself is dense 1..157.

### Alternate conformers

The walk keeps the **first** conformer of each residue whatever its altloc label.
An earlier draft whitelisted `{"", "A", "1"}` and silently dropped 11 residues of
7O7W whose only conformers are labelled `B`/`C`. Caught by the cross-layer check
below.

### Cross-layer agreement

Structure walk vs ConSurf mapped identities, chain A:

| Structure | ConSurf mapped | Walk | Disagreement |
|---|---|---|---|
| 7O7W_A | 236 | 236 | 0 |
| 1BRS_A | 108 | 108 | 0 |
| P00648 | 106 | 157 | 51 walk-only (model region ConSurf did not cover) |

### Known limitation

A residue is recognised by an atom named exactly `CA`. 7O7W's chromophore
`PIA:68:A` is a fused tripeptide with `CA1`/`CA2`/`CA3` and cannot be picked up
even with `include_hetatm` plus a custom name. ConSurf omits it too, so the two
stay consistent — but anything needing such a residue must find it another way.

## 9. Phase 2 — live mapping (done)

`msa_indices[atm.resid - 1]` is gone from the live path. Residues are now looked
up by identity, `(chain, resid, icode)`.

| Change | Where |
|---|---|
| MODELLER imported lazily | `sequence_processing.py` — was `from modeller import *` at module scope, which made the whole network chain unimportable without a licensed package |
| `pdb_to_fastas` uses the shared walk | `sequence_processing.py:~290` |
| `OtherAtom` carries `chain` + `icode` | `generate_static_networks.py:~104`, `generate_dynamic_networks.py:~107` |
| `ResidueIndex` built per structure | `_build_residue_index` in both builders |
| Identity lookup replaces arithmetic | `extract_objects`, `extract_objects_per_frame` |
| 7 duplicate graph-node lookups removed | now read `molecule.msa_resid` |
| `convert_msa_to_individual` exact-match | `sequence_processing.py:~414` |

### Why it mattered

| Structure | Residues | Chains | Old arithmetic correct | Silent wrap-around |
|---|---|---|---|---|
| 1AKI | 129 | 1 | 129/129 | 0 |
| P00648 (model) | 157 | 1 | 157/157 | 0 |
| **7O7W** | 236 | 1 | **0/236** | **6** |
| **1BRS** | 588 | **6** | **0/588** | 0 |

Dense, 1-origin, single-chain structures were right by luck. Everything else was
silently wrong.

### Behaviour changes (deliberate)

1. **`pdb_to_fastas` output changes for multi-chain and modified-residue PDBs.**
   The old predicate `('ATOM' in line) and ('CA' in line)` matched `HETATM` by
   substring, matched a stray `CA` anywhere on the line, ignored alternate
   conformers, and merged all chains. It also **crashed** on this repository's
   own test PDB: `KeyError: 'RY '`, from a `REMARK 500` line. Single-chain output
   is byte-identical.
2. **`convert_msa_to_individual` returns `None` instead of a neighbour** when the
   reference column is a gap. It previously used `np.argmin` to pick the nearest
   column.
3. **Dynamic no-MSA fallback.** `msa_indices = residues` (all residues, waters
   included) is replaced by an explicit residue-number mapping via
   `fallback_to_resid=True`.
4. **`generate_dynamic_networks` non-active-region oxygen path** previously wrote
   `MSA=None` while computing the index — a copy-paste bug. It now carries the
   value.
5. **`msa_classification/*.csv` produced before this change are invalid** and
   should be regenerated; their `MSA_Resid` column came from the old arithmetic.

### Still unverified

`msa_with_modeller` and `perform_structure_alignment` need licensed MODELLER,
which is not installable here (the PyPI `modeller` is an unrelated package).
They are untouched and untested. Everything else runs without MODELLER, because
`generate_msa_alignment` only invokes it when the alignment file is absent.

## 10. Deferred to later phases

- ConSurf metadata on `OtherAtom` / graph nodes (Phase 3)
- Water/network propagation (Phase 4); scoring (Phase 5)
- `residue_analysis.get_per_residue_interactions:45` — `msa=True` branch still
  unreachable; `WaterNetwork.get_per_residue_interactions` does not forward it
- Pre-existing crashers not on the Phase-2 path: `generate_static_networks.py`
  `get_density` (`S` before assignment for `selection != 'all'`), the directed
  `angle_criteria` branch (`water1` undefined), and
  `residue_analysis.get_all_water_distances` (arity mismatch)
- `test_inputs.py` fails on a missing `import sys` (pre-existing, tracked);
  `test_general.py` only passes when pytest runs from `WatCon/tests/`

---

## 11. Phase 3 — ConSurf attached to residues and waters (done)

### What was added

`WatCon/evolutionary.py`:

| Symbol | Role |
|---|---|
| `ResidueConservation` | frozen per-residue record: score, grade, low_confidence, MSA support, CI, b/e, f/s, ConSurf POS, source |
| `WaterConservation` | per-water aggregate: min_score, mean_score, max_grade, n_residues, n_low_confidence, n_unscored |
| `ConservationMap` | `(chain, resid, icode) -> ResidueConservation`; `for_residue`, `status`, `coverage` |
| `aggregate_water` | the aggregation; returns `None` when nothing is scored |
| `find_consurf_file` / `load_conservation` | discovery and loading of **supplied** results |

Wired in: `OtherAtom.evolutionary` (one shared immutable object per residue),
`WaterMolecule.evolutionary`, `WaterNetwork.annotate_water_conservation()`,
graph attributes `evo_score` / `evo_grade` / `evo_min_score` / `evo_n_residues`,
`initialize_network(consurf_directory=, consurf_chain_map=, consurf_strict=)`,
three `parse_inputs` keys, three appended CSV columns, and an
`evolutionary_coverage` entry in `metrics`.

### Scientific decisions and the evidence behind them

**Score is primary; grade is a label.** ConSurf's `score` is already normalised
within a run to mean 0, sd 1 (measured: mean -0.000, sd 0.995-0.998 across all
four datasets). It is a within-protein z-score of evolutionary rate -- comparable
as a rank inside its own protein, **not** as an absolute rate between proteins.
`grade` bins that score using per-run thresholds (layer widths 3.89-4.71 here).

Running the same protein twice at different MSA depths (P00648, 150 vs 50
sequences, 107 shared positions):

| | Spearman | mean abs diff | identical |
|---|---|---|---|
| score | **0.955** | 0.240 sd | — |
| grade | 0.947 | 0.505 | 57% |

Grade is the *less* stable representation -- binning amplifies small score moves
across bin edges. Both are stored; only score should be computed with.

At the conserved end, agreement is much better: for grade >= 8, 26/34 positions
identical, 8 differ by one, none by two. `ResidueConservation.is_conserved`
(grade >= 8) is labelled on that basis.

**Sign convention.** More negative score = more conserved; grade 9 = most
conserved. The two run opposite ways.

**Low confidence is flagged, never filtered.** Excluding ConSurf's low-confidence
positions moved cross-run agreement from 0.505 to 0.495 mean grade difference --
it costs data and buys almost nothing. Every residue keeps its flag and MSA
support; aggregates count low-confidence contributors.

**Water aggregation is multi-statistic and unweighted.** Contacts are
de-duplicated to residues first (`OtherAtom` is per-atom, so a residue reached
through two atoms would otherwise count twice). No single water-level number is
invented -- which summary is biologically right is still open, and committing now
would force a re-run to change it. Distance weighting is unavailable: the
connection tuples carry no distance.

**`None` means "no data", never "not conserved".** A water with no scored contact
gets `evolutionary = None`, not zero.

### Keeping the two conservations apart

WatCon already uses "conservation" for **structural water conservation**
(`find_conserved_networks.find_commonality`). Everything evolutionary is prefixed
`evo_` / `evolutionary_`. `find_conserved_networks` was **not touched**. The two
never share a dict key, CSV column, or plot axis.

### Behaviour notes

- The classification CSV gained `Evo_Score,Evo_Grade,Evo_LowConf` **appended at
  the end**; `plot_interactions_from_angles` and `identify_clustered_angles` read
  it by column position, so nothing before them moved. Unscored residues write
  `NA,NA,NA`.
- `consurf_directory=None` (the default) disables the feature entirely and leaves
  every residue unscored -- existing runs are unaffected.
- Fixed on the way: `residue_analysis.get_per_residue_interactions` used `or`
  chains and `if target_res:`, so residue number **0 was falsy** and silently
  dropped (7O7W's `PRO:0:A`). Now uses explicit `is not None`.

## 12. Deferred: automatic ConSurf submission (NOT implemented)

WatCon **only consumes ConSurf results you supply**. It never contacts the
ConSurf web server, submits a structure, or downloads results. Run ConSurf
yourself and place the `*_consurf_grades.txt` files in `consurf_directory`.

When we do add a runtime/API hook, the seam is
**`WatCon.evolutionary.find_consurf_file`** -- the single place that decides which
ConSurf file belongs to a structure. A future hook should:

1. submit the structure (or its sequence) to ConSurf and poll for completion;
2. write the returned `*_consurf_grades.txt` into `consurf_directory`;
3. return control to `load_conservation`, which is unchanged.

Everything downstream of `load_conservation` is already agnostic to how the file
arrived. Points needing decisions when that work starts: credentials and rate
limiting, caching and cache invalidation, run parameters (homolog search vs
supplied MSA, MSA depth, Bayesian vs ML) recorded as provenance, and offline /
air-gapped behaviour.

## 13. Still deferred to later phases

- **Phase 4** — network/site-level rollups; distance-weighted aggregation, which
  needs the connection tuple extended (it is already heterogeneous: undirected
  WAT-PROT has 6 elements, directed has 5, and `classify_waters` indexes
  `connection[5]`).
- **Phase 5** — whether a combined structural + evolutionary score should exist
  at all, and cross-protein renormalisation.
- **Phase 6** — end-to-end validation and regression across the full corpus.
- Unrelated pre-existing defects: `get_density` (`S` before assignment for
  `selection != 'all'`), the directed `angle_criteria` branch (`water1`
  undefined), `get_all_water_distances` (arity mismatch), `test_inputs.py`
  (missing `import sys`), `test_general.py` (cwd-dependent relative path).

---

## 14. Phase 4 — the consumer end (done)

### The gap this closed

After Phase 3, conservation was attached to residues and waters but **nothing
read it**: `find_conserved_networks.py` had 0 references, `visualize_structures.py`
had 0. The plumbing existed with no faucet.

### What was added

| Symbol | Module | Purpose |
|---|---|---|
| `water_residue_contacts` | `evolutionary.py` | the single contact walk, shared by both builders and the cluster join |
| `aggregate_site` | `evolutionary.py` | conservation of residues lining any water set (e.g. the active region) |
| `ClusterConservation` | `evolutionary.py` | one conserved water site: structural occupancy **and** evolutionary conservation |
| `conservation_of_clusters` | `evolutionary.py` | **the cross-structure join** |
| `write_conservation_report` | `evolutionary.py` | one CSV row per conserved water site |
| `project_clusters_by_conservation` | `visualize_structures.py` | cluster PDB, evolutionary conservation in B-factors |
| `pymol_project_evolutionary` | `visualize_structures.py` | residues coloured on ConSurf's own 1-9 scale |
| `conservation_report` | `WatCon.py` | post-analysis option wiring it together |

### The join

```
cluster centre --(within dist_cutoff)--> waters, per structure
                                              |
                                        contacting residues
                                              |
                        de-duplicated ACROSS structures
                                              |
                                      ClusterConservation
```

A residue lining a site in five structures contributes **once**.

### Deliberate design choices

- **Never blended.** The report carries `occupancy` (structural) and `evo_*`
  (evolutionary) as separate columns. Whether they correlate is the research
  question; the tool supplies the inputs and stops.
- **Two separate projection files**, never one combined B-factor. Load both and
  compare.
- **`None` means no data.** A site with no scored lining residue reports `NA`,
  never 0 — 0 would read as "not conserved".
- **Unweighted**, consistent with the measurement that distance weighting sits
  below ConSurf's own reproducibility noise.

### Bugs fixed in this phase

1. **`get_density(selection != 'all')` raised `UnboundLocalError`** —
   `generate_static_networks.py:1018` and `generate_dynamic_networks.py:1094`
   read `S.edges` before `S` was assigned. One line each. This blocked every
   region-restricted analysis. *(Correction to an earlier note: this affected one
   method, not six — the five siblings were already correct.)*
2. **B-factor written one column too far right** in `project_clusters`. The PDB
   spec puts tempFactor at columns 61-66; WatCon emitted an extra space, so a
   spec-compliant reader truncated `0.25` to `0.2`. Fixed in both the existing
   and the new projection so the two agree.
3. **numpy integers leaking into residue keys.** `atom.resid` from MDAnalysis is
   `np.int64`; it hashes like a Python int so lookups worked, but keys serialised
   as `np.int64(70)`. Now coerced.

### Validated

355 tests passing (up from 318): 24 new for the cluster join, 13 for the
projections. Real structures used throughout — 7O7W (negative resids, omitted
chromophore), 1BRS (six chains), 1AKI (no ConSurf data).

### NOT validated — the biology

**Mechanical correctness is verified; the scientific question is not answered.**
Testing whether conserved water sites *are* lined by conserved residues needs a
real protein family with ConSurf data for every member. We hold four ConSurf
datasets covering two proteins. The machinery is ready; the data is not.

---

## 15. Phase 5 — identity, superposition, family scaffold, and the study (done)

The integration was complete but had never been run on more than two structures,
and the question it exists to answer was unanswered. This phase closed both.

### 15.1 Identity verification — the last correctness gap

`ConservationMap` matched residues on `(chain, resid, icode)` only. Nothing in
that key says *which amino acid* a score belongs to, so a ConSurf run numbered
differently from the structure attaches a score to every residue, reports 100%
coverage, and is wrong everywhere. `ResidueConservation` was already parsing the
amino acid from both the SEQ and ATOM columns and discarding it.

Now kept as `amino_acid` / `residue_name` and compared in `coverage()`.

It is a **rate, not a pass/fail**, because a point mutant *should* disagree at
the mutated position. Measured on real barnase:

| case | identity rate |
|---|---|
| correct numbering | 100% (108/108) |
| one point mutation | 99.1% (107/108) |
| numbering shifted by one residue | **2.8%** (3/107) |

`DEFAULT_IDENTITY_THRESHOLD = 0.95` sits in that gap. `enforce_identity()`
raises under `consurf_strict`, warns otherwise, and names the disagreeing
residues. Objects with no residue name report "not checked", never
"all mismatched"; counting is per residue, not per atom.

**It fired on real data.** In the barnase study a mis-numbered structure (2F56,
offset +2) reached the builder and was stopped with *"only 2 of 106 compared
residues agree (1.9%)"*.

### 15.2 `WatCon/superpose.py` — alignment without MODELLER

`align_with_waters` already applies transforms to every atom including waters;
it only needed a source of matrices that is not MODELLER. `kabsch()` supplies
them in numpy, in exactly the `{'Rot': [...], 'Trans': [...]}` shape
`perform_structure_alignment` returns, so the two are interchangeable. No
existing module changed.

Correspondence is by **residue number**, so this covers one protein — different
crystal forms, mutants, complexes. It is not a sequence aligner.

Two traps, both with named tests: the reflection correction (without it the SVD
returns an improper rotation and superposes onto a mirror image), and keying the
correspondence on the chain letter (which made the three barnase copies in 1BRS
share zero residues).

### 15.3 Family scaffold — and a silent wrong-file bug

`conservation_by_msa_column` pools conservation from *separate* ConSurf runs at
shared MSA columns, so a real family needs ConSurf files, not code.

Validating it exposed a genuine bug: `find_consurf_file("P00648_50", ...)`
returned **`P00648_150.grades.txt`**. The lookup tried the loose token (before
the first underscore) *before* the full stem, so the prefix `P00648` matched both
and the alphabetically first won. Now the most specific token is tried first, and
an ambiguous match is refused rather than resolved arbitrarily.

**Measured cross-run agreement**, the honest floor for any family claim:

| pair | score ρ | grade ρ | grade changed |
|---|---|---|---|
| P00648 depth 150 vs depth 50 | 0.955 | 0.947 | 42% |
| 1BRS run vs P00648 depth 150 | **0.380** | 0.365 | 62% |
| 1BRS run vs P00648 depth 50 | **0.367** | 0.360 | 66% |

The same protein, run from a different starting structure, agrees at only
ρ ≈ 0.37. Recorded in `CROSS_RUN_AGREEMENT`. Grades are never averaged (ordinal
per-run bins); `spread` accompanies every pooled value; the recommended statistic
is `unanimous_conserved`, which uses each run's own verdict rather than comparing
separately normalised scores.

### 15.4 The study — `experiments/barnase_waters/`

22 barnase structures at ≤2.0 Å, one ConSurf run, 2753 waters, 594 sites, 462
with conservation.

> **Spearman ρ = −0.375**, p = 6.6×10⁻¹⁷, against a permutation null centred on
> +0.033 (z = −4.47). Stable at ≤1.8 Å (−0.332) and wild-type-only (−0.338).

The pre-registered statistic (`evo_min_score`) had to be corrected **by its own
control**: its shuffled null centres on −0.226, because occupancy correlates with
the *number* of lining residues (+0.555) and the minimum of a larger set is
mechanically lower (−0.553). The mean over lining residues has no such
dependence. Both are reported.

**The declared confound is unresolved**: buried residues are both more conserved
(median −0.694 vs +0.180) and more likely to hold ordered water. No causal claim
is made.

See `experiments/barnase_waters/FINDINGS.md` and `PREREGISTRATION.md`.

### 15.5 Defect backlog cleared

| defect | effect |
|---|---|
| `max_neigbhbors` typo | static **directed** networks had never executed in any released version |
| `norm(water1)` at static:617 | cosine normalised by an undefined vector (`:655`/`:702` were correct) |
| `test_inputs.py` missing `import sys` | its only test raised `NameError` |
| `test_general.py` cwd-dependent path | passed only when pytest ran from `WatCon/tests/` |

Two reported defects were **dropped after checking**: `sequence_processing.py`
does *not* define functions twice (lines 563–656 are inside a `'''...'''` block
that `grep '^def '` matched), and `get_all_water_distances` is genuinely broken
but dead and written against an API that no longer exists, so repairing it would
mean inventing a specification.

### 15.6 What remains unverified

- The family scaffold is exercised on two runs of **one** sequence. No real
  multi-protein family has been run — the join is tested, the biology is not.
- The study is **one protein**, correlational, with burial unresolved.
- The directed static path now executes and produces a directed graph, but its
  hydrogen-bond geometry is still covered by no test.
- `2F56`/`2F5M` were excluded rather than renumbered; the renumbered sensitivity
  check has not been run.

---

## 16. What you can do with this

Verified against the running code, not recalled.

### Supported today

| capability | entry point |
|---|---|
| ConSurf parsing, both dialects, cross-checked against ConSurf's annotated PDB | `watcon validate`, `consurf/parser.py` |
| **Numbering cannot silently mis-map** -- the guarantee published WatCon lacks | `enforce_identity`, `ConservationMap.coverage` |
| Conservation on residues, on waters, **and on the network graph** | protein nodes carry `evo_grade`/`evo_score`; water nodes carry `evo_min_score`/`evo_n_residues` |
| Which residues each water contacts | `water_residue_contacts` |
| Conserved-site to conservation join, at family scale | `conservation_of_clusters` |
| Summary over an arbitrary region or water set | `aggregate_site`, `aggregate_water` |
| Pooling separate ConSurf runs at shared MSA columns | `conservation_by_msa_column`, `family_summary` |
| One row per conserved site, CSV | `write_conservation_report` |
| Dataset preparation without a MODELLER licence | `watcon prepare`, `superpose.kabsch` |
| A PyMOL session that opens | `watcon view` |

The graph point is the one with the most reach: because conservation sits on
NetworkX node attributes, **any graph analysis can be filtered or weighted by it
without new code** -- shortest paths, centrality, subgraph extraction.

### Buildable on what exists

* **Water-mediated bridges between conserved residues.** The graph already holds
  conservation on nodes and water-protein edges; ask which conserved residue
  pairs are connected through a water.
* **Variant or engineering triage.** Given candidate mutation sites, which line a
  conserved water position? The residue -> water -> conservation lookup is there.
* **Conformational-state comparison.** WatCon already separates open/closed or
  apo/holo sets; conservation says which state's waters sit on conserved
  residues.
* **Family-level conserved-water maps.** `conservation_by_msa_column` is built
  and unit-tested; it needs one ConSurf run per family member.
* **Burial-matched analysis.** The sharpest open question from the benchmark --
  whether conservation separates sites beyond acting as a burial proxy.

### What the benchmark showed, so nobody builds on a false premise

Conservation **does not** improve prediction of where water sits: dAUC +0.0004
(barnase) and +0.0001/-0.0002 (PTP1B), all inside their own permutation nulls,
against a pre-registered bar of 0.02. It is *redundant* with occupancy, not
uninformative -- alone it scores AUC 0.65 and 0.57.

What it does do is separate sites that occupancy alone conflates: the two
measures disagree at 36-42% of sites, and sites that recur *without* conservation
are overwhelmingly surface. Burial and conservation are correlated, so part of
that separation is a burial proxy. See `experiments/benchmark/FINDINGS.md`.


# Benchmark: PTP1B at scale, and the family of fifteen

The evidence behind the claims in the top-level README that are not barnase.
Read [FINDINGS.md](FINDINGS.md) first: what holds up, what does not, and what
must not be claimed.

| Document | What |
|---|---|
| [FINDINGS.md](FINDINGS.md) | every result, each with its limits |
| [PHASE_A_MAPPING.md](PHASE_A_MAPPING.md) | the residue-mapping defect, measured against sequence-derived ground truth |
| [PREREGISTRATION.md](PREREGISTRATION.md) | hypothesis and method, written before the data were looked at |
| [data/consurf_raw/MANIFEST.md](data/consurf_raw/MANIFEST.md) | which ConSurf job number is which protein, and the settings all fifteen shared |

## What is here

* `results/` — the numbers the findings are computed from: benchmark JSON,
  mapping and disagreement CSVs. Run logs are excluded as noise.
* `scripts/` — fetch, prepare, benchmark, compare, and the two builders
  (`build_family15.py`, `build_ensemble_fixtures.py`).
* `data/ptp1b_selected.txt` — the 253 PTP1B entries, and
  `ptp1b_resolutions.json` — how they were chosen.

## What is not here, and why

**Structures and ConSurf bundles are not committed.** The full working tree is
about **1.1 GB**: 287 PTP1B entries, their prepared copies, the family
structures, the Zenodo MD system, and fifteen ConSurf result bundles at roughly
32 MB each. None of it is ours to redistribute and all of it is retrievable —
RCSB serves the structures, Zenodo the MD system
([10.5281/zenodo.15213225](https://doi.org/10.5281/zenodo.15213225)), and
ConSurf the runs, whose settings are recorded in the manifest above.

So what is committed is what a reader cannot regenerate: the findings, the
method, and the numbers.

## Reproducing the family of fifteen

```bash
pip install -e ".[test]"
python experiments/benchmark/scripts/build_family15.py
watcon family \
    --members experiments/benchmark/ptp_family15/members15.tsv \
    --alignment WatCon/data/examples/ptp_family/ptp_family_alignment.pir
```

Expect `344 columns, 229 covered by all 15 proteins, 58 of those unanimously
conserved`, and 556 occupied water sites of which 117 are lined by a unanimously
conserved column.

The ConSurf grades files the members list points at are **in the repository**,
under `WatCon/data/consurf/fixtures/` — only the 32 MB bundles they came from
are not. Structures are fetched with **native numbering**, which matters: the
authors' own copies are renumbered from 1 over resolved residues, so crystal
gaps collapse, no constant offset recovers the sequence, and `enforce_identity`
correctly refuses them (FINDINGS.md, section 6).

To repeat the comparison with the pseudophosphatase excluded, use
`members14_no_pseudophosphatase.tsv`. It should give **58** unanimous columns
too — that equality is the finding.

# ConSurf exploratory data

Complete ConSurf result bundles and MSAs, kept as evidence. **Not read by the
test suite** except `test_consurf_crosscheck.py`, which streams the annotated
PDBs out of the archives without extracting them.

## Result bundles

| Archive | Job | Protein | MSA depth | Homolog search | Note |
|---|---|---|---|---|---|
| `1788241750_ConSurf.tar.gz` | 1788241750 | 7O7W chain A | 45 | yes | crystal structure |
| `1788241750_ConSurf (1).tar.gz` | 1788241750 | 7O7W chain A | 45 | yes | **duplicate download** |
| `1788241883_ConSurf.tar.gz` | 1788241883 | P00648 | 50 | yes | AlphaFold model |
| `1788241897_ConSurf.tar.gz` | 1788241897 | 1BRS chain A | 150 | **no** | supplied MSA |
| `1788242051_ConSurf.tar.gz` | 1788242051 | 1BRS chain A | 150 | **no** | supplied MSA |

Bundles without a homolog search lack `sequences_found_hmmer.tar.gz`,
`query_cdhit.tar.gz`, `filtered_homolougues.tar.gz` and
`query_rejected_homologues.tar.gz`. The grades format is identical either way.

Five archives contain only **three distinct** grades files: the two 1788241750
bundles are identical, and 1788241897 / 1788242051 are identical to each other.

## Useful members

- `*_consurf_grades.txt` — the file the parser reads.
- `*_ATOMS_section_With_ConSurf.pdb` — grade written into the B-factor column,
  blank where unscored. Used as the cross-check oracle.
- `msa_fasta.aln` — the MSA. For 7O7W this is what shows the query has 250
  ungapped residues with `X` at position 79, against 249 grades records.
- `TheTree.txt` — Newick tree; query is named `Input_seq_SEQRES_A`.
- `msa_aa_variety_percentage.csv` — **do not use as an oracle.** Its composition
  columns are shifted one row relative to its own `pos`/`grade` columns. See
  `docs/CONSURF_INTEGRATION.md` section 6.

## Loose MSAs

Eight `msa_fasta*.aln` files from ConSurf runs where only the alignment was kept.
At least five distinct proteins, including a depth series (50/50/100/200) on one
sequence. Not parsed by anything; retained as format evidence.

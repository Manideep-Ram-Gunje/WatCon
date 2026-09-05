# ConSurf fixtures

Canonical grades files read by the test suite. Five files, **four distinct
datasets** — `1BRS_A_150.crlf` is a byte-identical CRLF twin of `1BRS_A_150`,
kept deliberately so line-ending handling stays covered.

| File | Source bundle | Protein | MSA depth | Line endings |
|---|---|---|---|---|
| `1BRS_A_150.grades.txt` | job 1788241897 / 1788242051 | 1BRS chain A (barnase) | 150 | LF |
| `1BRS_A_150.crlf.grades.txt` | same content, CRLF copy | 1BRS chain A | 150 | CRLF |
| `P00648_150.grades.txt` | earlier run | P00648 (AlphaFold model) | 150 | CRLF |
| `P00648_50.grades.txt` | job 1788241883 | P00648 (AlphaFold model) | 50 | LF |
| `7O7W_A_45.grades.txt` | job 1788241750 | 7O7W chain A (rsEGFP2) | 45 | LF |

What each one exercises:

- **1BRS_A_150** — baseline; 2 leading unmapped residues; 3 fully conserved.
- **P00648_150 / _50** — UniProt numbering (POS 1 -> residue 48, offset +47);
  **interior** unmapped residue at POS 51; predicted rather than experimental
  structure; two MSA depths of the same protein.
- **7O7W_A_45** — the hard one: **negative and zero residue numbers**
  (`HIS:-5:A` … `PRO:0:A`), a **POS discontinuity** at 65->69 where ConSurf
  omitted the chromophore `PIA:68:A`, `X` in the residue variety, 41 fully
  conserved positions, and a shallow 45-sequence MSA.

Full result bundles, including the annotated PDBs used for the B-factor
cross-check, are in `../exploratory/`.

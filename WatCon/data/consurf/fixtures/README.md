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

## Protein tyrosine phosphatase family

Fifteen runs, one per member, chain A, all webserver / Bayesian / 150 homologues,
all LF. Used by `tests/test_ptp_family_runs.py`, `tests/test_family_fifteen.py`
and the family path.

| Grades file | CA extract | PDB | UniProt | Gene | Query state |
|---|---|---|---|---|---|
| `1AAX_A.grades.txt` | `1AAX_A.consurf_ca.pdb` | 1AAX | P18031 | PTPN1 | C215S trap |
| `4GRZ_A.grades.txt` | `4GRZ_A.consurf_ca.pdb` | 4GRZ | P29350 | PTPN6 | C453S trap |
| `1ZC0_A.grades.txt` | `1ZC0_A.consurf_ca.pdb` | 1ZC0 | P35236 | PTPN7 | wild type |
| `5HDE_A.grades.txt` | `5HDE_A.consurf_ca.pdb` | 5HDE | Q05209 | PTPN12 | CSP231 phospho-Cys |
| `3BRH_A.grades.txt` | `3BRH_A.consurf_ca.pdb` | 3BRH | Q9Y2R2 | PTPN22 | C227S + D195A trap |
| `4S0G_A.grades.txt` | `4S0G_A.consurf_ca.pdb` | 4S0G | P26045 | PTPN3 | C842S trap |
| `2I75_A.grades.txt` | `2I75_A.consurf_ca.pdb` | 2I75 | P29074 | PTPN4 | wild type |
| `8SLS_A.grades.txt` | `8SLS_A.consurf_ca.pdb` | 8SLS | P54829 | PTPN5 | wild type |
| `6KZQ_A.grades.txt` | `6KZQ_A.consurf_ca.pdb` | 6KZQ | P43378 | PTPN9 | C515A trap |
| `3ZM1_A.grades.txt` | `3ZM1_A.consurf_ca.pdb` | 3ZM1 | Q06124 | PTPN11 | wild type |
| `1WCH_A.grades.txt` | `1WCH_A.consurf_ca.pdb` | 1WCH | Q12923 | PTPN13 | wild type |
| `6IWD_A.grades.txt` | `6IWD_A.consurf_ca.pdb` | 6IWD | Q15678 | PTPN14 | wild type |
| `4GFU_A.grades.txt` | `4GFU_A.consurf_ca.pdb` | 4GFU | Q99952 | PTPN18 | C229S trap |
| `8GVV_A.grades.txt` | `8GVV_A.consurf_ca.pdb` | 8GVV | Q16825 | PTPN21 | C1108S trap |
| `2QEP_A.grades.txt` | `2QEP_A.consurf_ca.pdb` | 2QEP | Q92932 | PTPRN2 | pseudophosphatase |

The `*.consurf_ca.pdb` files keep only the CA line of ConSurf's
`*_ATOMS_section_With_ConSurf.pdb` (the grade is written on every atom, so the
CA carries it). They reproduce the full file's cross-check exactly — identical
checked/agree/mismatch/missing/blank counts, 4289/4289 over the fifteen —
at ~25 kB instead of 245–507 kB. The full bundles (~33 MB each) live outside this repository in
`experiments/benchmark/data/consurf_raw/`.

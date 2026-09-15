# PTP family example data

Real data for the family path and its tests.

| File | What | Source |
|---|---|---|
| `5HDE_A_csp_site.pdb` | 5HDE chain A trimmed to 12 A around CSP231, with waters | RCSB 5HDE |
| `<ID>_A_ca.pdb` | CA atoms of chain A, native numbering, coordinates unchanged | RCSB 2F71, 8U1E, 4GRZ, 4HJP, 1ZC0, 3O4U, 5HDE, 5J8R, 3BRH, 3OLR |
| `ptp_family_alignment.pir` | the ten corresponding rows of the authors' family alignment, verbatim | Brownless, Harrison-Rawn & Kamerlin, *JACS Au* 2025; Zenodo 10.5281/zenodo.15213225 (CC-BY-4.0), `PTPs_combined/alignment.txt` |

The alignment rows correspond to the authors' renumbered structures, not these
native-numbered ones; the family path maps each structure to its row **by
sequence**, never by position or residue number.

Two defects of that alignment are kept deliberately, because the mapper must
catch them: the 5HDE row omits CSP231 (the catalytic phosphocysteine), and the
3O4U row slides E124, D125 and E241 across disordered loops.

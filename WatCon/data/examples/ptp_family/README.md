# PTP family example data

Real data for the family path and its tests.

| File | What | Source |
|---|---|---|
| `5HDE_A_csp_site.pdb` | 5HDE chain A trimmed to 12 A around CSP231, with waters | RCSB 5HDE |
| `<ID>_A_ca.pdb` | CA atoms of chain A, native numbering, coordinates unchanged | RCSB, all twenty-four entries listed below |
| `<ID>_A_site.pdb` | the closed structure of each protein, trimmed to 12 A around its own catalytic nucleophile, waters kept, ligands dropped | RCSB 2F71, 4GRZ, 1ZC0, 3BRH (5HDE's equivalent is `5HDE_A_csp_site.pdb`) |
| `ptp_family_alignment.pir` | the authors' family alignment, all twenty-four rows, verbatim | Brownless, Harrison-Rawn & Kamerlin, *JACS Au* 2025; Zenodo 10.5281/zenodo.15213225 (CC-BY-4.0), `PTPs_combined/alignment.txt` |

The alignment rows correspond to the authors' renumbered structures, not these
native-numbered ones; the family path maps each structure to its row **by
sequence**, never by position or residue number.

Two defects of that alignment are kept deliberately, because the mapper must
catch them: the 5HDE row omits CSP231 (the catalytic phosphocysteine), and the
3O4U row slides E124, D125 and E241 across disordered loops.

## The fifteen proteins

Fifteen human protein tyrosine phosphatases, twenty-four structures, one ConSurf
run each. Every run reproduces its alignment row exactly, residue for residue.

| Protein | ConSurf run | Structures | Nucleophile as modelled |
|---|---|---|---|
| PTPN1 | 1AAX | 2F71, 8U1E | Cys215 |
| PTPN3 | 4S0G | 4S0G, 2B49 | C842S |
| PTPN4 | 2I75 | 2I75 | Cys852 |
| PTPN5 | 8SLS | 8SLS | Cys472 |
| PTPN6 | 4GRZ | 4GRZ, 4HJP | C453S / Cys |
| PTPN7 | 1ZC0 | 1ZC0, 3O4U | Cys270 |
| PTPN9 | 6KZQ | 6KZQ, 4GE6 | C515A |
| PTPN11 | 3ZM1 | 3ZM1 | Cys459 |
| PTPN12 | 5HDE | 5HDE, 5J8R | CSP231 (phospho-Cys) |
| PTPN13 | 1WCH | 1WCH | Cys2408 |
| PTPN14 | 6IWD | 6IWD | Cys1121 |
| PTPN18 | 4GFU | 4GFU, 2OC3 | C229S |
| PTPN21 | 8GVV | 8GVV, 8GWH | C1108S |
| PTPN22 | 3BRH | 3BRH, 3OLR | C227S |
| PTPRN2 | 2QEP | 2QEP | Cys945, in a dead P-loop |

Ten of the twenty-four are catalytically dead: seven C->S traps, one C->A
(PTPN9), 5HDE's phosphocysteine reaction intermediate, and PTPRN2.

**PTPRN2 is a pseudophosphatase.** Its P-loop reads `CSDGAGR` where the family
reads `CSAGIGR`, which is why IA-2beta has no activity. It is kept in the family
because removing it changes the conservation result not at all -- 58
unanimously conserved columns either way -- and because it is the clearest case
in the set of a column that is *constrained* without being *identical*. The
authors' dataset labels this entry PTPN2; the entry is PTPRN2 (Q92932), and real
PTPN2 is P17706.

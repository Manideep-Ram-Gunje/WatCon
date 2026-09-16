# ConSurf result bundles

Fifteen ConSurf runs, one per protein tyrosine phosphatase. Kept here
rather than in the WatCon-ConSurf repository because the bundles total
~470 MB; only the `*_consurf_grades.txt` and a CA extract of the
annotated PDB are copied into the repository as fixtures.

ConSurf names its downloads after an opaque job number, so the mapping
is recorded here. Every run used identical settings, which is what makes
them poolable: HMMER, 1 iteration, E-value 1e-4, UNIREF-90; automatic
selection of 150 sequences sampling the homologue list, max 95% / min 35%
identity; MAFFT-L-INS-i; Bayesian with the best-fitting model.

| Job | PDB | Chain | Gene | UniProt | scored | mapped | MSA | method | size |
|---|---|---|---|---|---|---|---|---|---|
| 1789014979 | 1AAX | A | PTPN1 | P18031 | 321 | 297 | 150 | bayesian | 31 MB |
| 1789492637 | 4GRZ | A | PTPN6 | P29350 | 288 | 282 | 150 | bayesian | 34 MB |
| 1789492776 | 1ZC0 | A | PTPN7 | P35236 | 309 | 286 | 150 | bayesian | 33 MB |
| 1789492820 | 5HDE | A | PTPN12 | Q05209 | 307 | 300 | 150 | bayesian | 33 MB |
| 1789492864 | 3BRH | A | PTPN22 | Q9Y2R2 | 310 | 296 | 150 | bayesian | 33 MB |
| 1789530895 | 4S0G | A | PTPN3 | P26045 | 306 | 279 | 150 | bayesian | 32 MB |
| 1789530908 | 2I75 | A | PTPN4 | P29074 | 320 | 273 | 150 | bayesian | 33 MB |
| 1789530924 | 8SLS | A | PTPN5 | P54829 | 282 | 282 | 150 | bayesian | 32 MB |
| 1789530944 | 6KZQ | A | PTPN9 | P43378 | 307 | 296 | 150 | bayesian | 32 MB |
| 1789530967 | 3ZM1 | A | PTPN11 | Q06124 | 284 | 259 | 150 | bayesian | 34 MB |
| 1789530984 | 1WCH | A | PTPN13 | Q12923 | 315 | 308 | 150 | bayesian | 32 MB |
| 1789530998 | 6IWD | A | PTPN14 | Q15678 | 303 | 289 | 150 | bayesian | 30 MB |
| 1789531015 | 4GFU | A | PTPN18 | Q99952 | 297 | 268 | 150 | bayesian | 32 MB |
| 1789531030 | 8GVV | A | PTPN21 | Q16825 | 299 | 286 | 150 | bayesian | 29 MB |
| 1789531046 | 2QEP | A | PTPRN2 | Q92932 | 304 | 288 | 150 | bayesian | 33 MB |

`scored` counts rows in the grades table; `mapped` counts those ConSurf
could place on the structure. The difference is residues in the query
sequence that the crystal did not resolve.

Verified 2026-09-16: all fifteen parse in strict mode, and every mapped grade
agrees with the grade ConSurf itself wrote into the B-factor column of
`*_ATOMS_section_With_ConSurf.pdb`.

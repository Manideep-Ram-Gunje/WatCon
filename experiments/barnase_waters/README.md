# Barnase: are conserved water sites lined by conserved residues?

The study behind the headline result in the top-level README:

> Spearman rho = -0.375 (p = 6.6e-17, n = 462 sites) against a permutation null
> centred on +0.033, over 22 barnase crystal structures and one ConSurf run.

Read [FINDINGS.md](FINDINGS.md) for the result, the method and the caveats, and
[PREREGISTRATION.md](PREREGISTRATION.md) for what was committed to before the
data were looked at.

## What is here

| Path | What |
|---|---|
| `FINDINGS.md` | the result, with its limits stated |
| `PREREGISTRATION.md` | hypothesis, method and stopping rule, written first |
| `data/` | the 22 chosen entries and how they were chosen |
| `consurf/` | the ConSurf grades files, one per structure |
| `results/` | the outputs the findings are computed from |
| `scripts/` | fetch, prepare, analyse, and a runner for all three |

## What is not here, and why

`structures/` (7.6 MB of raw PDB files) and `prepared/` (1.8 MB of superposed
copies) are **not committed**: both are regenerated exactly by the scripts, from
the entry list in `data/`. Keeping them would multiply the repository size for
files RCSB already serves.

## Reproducing it

```bash
pip install -e ".[test]"
python experiments/barnase_waters/scripts/run_study.py
```

That fetches the structures, prepares them, and rewrites `results/`. The scripts
use no absolute paths and take about two minutes over a normal connection.

Conservation is a property of the **sequence**, so the 22 grades files in
`consurf/` are the same run named once per structure — that is the file-matching
convention WatCon uses, not 22 separate submissions.

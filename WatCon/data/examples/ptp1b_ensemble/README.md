# Ensemble and MD fixtures

Real coordinates for the two paths that crystal structures cannot exercise:
reading a **trajectory**, and building a **directed** network from hydrogens.

Rebuild with `experiments/benchmark/scripts/build_ensemble_fixtures.py`.

| File | What | Source |
|---|---|---|
| `ptp1b_ensemble_ragged.pdb` | three superposed PTP1B structures, each trimmed to 12 Å around the catalytic centre, as a multi-model PDB | RCSB 1AAX, 1C83, 1C88 |
| `ptp1b_ensemble_constant.pdb` | the same three, reduced to atoms every model holds | as above |
| `md_active_site_h.pdb` | the active site of the Zenodo MD system, with its hydrogens | Zenodo 10.5281/zenodo.15213225, `MD/starting_structures/WT_PTP1B_Apo_Closed.gro` |

## Why there are two ensembles

`ptp1b_ensemble_ragged.pdb` is what an ensemble of crystal structures actually
is: 397, 395 and 396 atoms, because each entry resolves a different set of
waters. **It cannot be read as a trajectory**, and it is here to prove that the
dynamic path says so. MDAnalysis builds the topology from the first model alone,
so a ragged file is accepted, reports a frame count, and fails only once a later
frame is read — which is why `structure_io.require_constant_atom_count` refuses
it at the door and points at the static path instead.

`ptp1b_ensemble_constant.pdb` is the same three structures made readable: the
protein atoms all three resolve, plus the **9 water oxygens nearest the catalytic
centre** in each. Every coordinate is a real crystallographic position, but the
fixed water count is a *stated selection*, not a property of the structures. It
exists so the dynamic path can be tested on real water positions. **No
scientific claim rests on it** — a real crystal ensemble shares no common set of
waters, and belongs in the static path, which reads each structure separately.

## The MD fixture and its numbering

`md_active_site_h.pdb` is the only real hydrogen-bearing PTP1B system available,
so the directed path needs it. Amber ff14SB, waters named `WAT`, protonation
states carried in the residue names (`HID`, `HIE`, `HIP`, `CYM`, `GLH`).

Its numbering is the MD system's own: **resid 214 here is PTP1B's Cys215**. The
whole system matches the 1AAX ConSurf run at offset **+1** with identity
**1.000** over 295 residues, every other offset below 0.09. That offset is
stated wherever the fixture is used and never applied silently.

## A note on cost

The network itself is cheap. On the full MD system (4,851 waters) building it
takes **1.4 s**, while `shortest_path` and `characteristic_path_length` — both
on when `analysis_conditions='all'` — take **819 s** and **45 s**, because they
run all-pairs shortest paths over the whole graph. The tests here switch those
two off; so should anyone running a real trajectory.

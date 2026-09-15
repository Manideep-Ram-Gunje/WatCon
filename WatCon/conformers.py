"""Alternate conformers: drop a position only when another is more populated.

Crystal structures model some atoms in two or more alternate positions, each
labelled (A, B, C, ...) with a fractional occupancy. MDAnalysis keeps every
position as a separate atom, so WatCon's network builders turned one
hydrogen-bonding atom into two nodes and one water contact into two edges.

What was measured before choosing a rule
----------------------------------------
* Duplicate water-protein edges: 1.0% (2F71), 4.3% (8U1E), and about 70% in the
  PanDDA ensemble models of PTP1B (5QFP, 5QFS, 5QE2, 5QF8).
* **Exact occupancy ties are the norm, not an edge case.** 92% of residues with
  alternate positions in the 253 prepared PTP1B structures (21,719 of 23,521),
  and 95% in the bundled barnase set, have no single most-populated position --
  typically 0.50/0.50, or four positions at 0.25 in the PanDDA ensembles. In the
  ten raw PTP family structures ties are 28%.
* So a "highest occupancy, ties to the first-listed" rule is decided almost
  entirely by the tie-break, which is arbitrary. On barnase it discarded real
  contacts: a grade-9 site lined only by Asp54 lost it, because its water touched
  conformer B at 0.50 (165 -> 164 sites with conservation, 57 -> 55 conserved).

The rule
--------
Within a residue, compare conformers by mean occupancy. **Keep every conformer
tied for the highest; drop only those strictly lower.** A state is never
discarded without evidence that another is more populated. Atoms with no
alternate label are always kept.

Costs, stated rather than hidden
--------------------------------
* Equally populated duplicates remain, so WatCon's own edge counts stay inflated
  for structures with many tied positions -- the PanDDA ensembles above all.
* A WatCon water node holds one oxygen. A water with two *tied* positions keeps
  the builder's existing behaviour -- the first surviving oxygen -- because two
  positions cannot share one node. 470 such waters in the PTP1B set, none in
  barnase.

This module does not import MDAnalysis; it operates on an AtomGroup it is given.
"""

from __future__ import annotations

from collections import OrderedDict
from typing import Dict, List, Tuple

__all__ = [
    "OCCUPANCY_TOLERANCE",
    "conformer_choices",
    "has_alternate_conformers",
    "highest_occupancy_atoms",
]

#: Mean occupancies closer than this are tied. Deposited occupancies carry two
#: decimals, so 0.50 and 0.50 are equal and 0.49 and 0.51 are not.
OCCUPANCY_TOLERANCE = 1e-6

ResidueKey = Tuple[str, int, str, str]


def _labels(atoms):
    labels = getattr(atoms, "altLocs", None)
    if labels is None:
        return None
    return ["" if label is None else str(label).strip() for label in labels]


def has_alternate_conformers(atoms) -> bool:
    """True only when the topology carries at least one alternate-location label.

    The network builders apply :func:`highest_occupancy_atoms` only then. An MD
    topology has no altlocs, and leaving it untouched keeps its per-frame
    ``updating=True`` selections exactly as they were.
    """
    labels = _labels(atoms)
    return bool(labels) and any(labels)


def _residue_keys(atoms, labels) -> List[ResidueKey]:
    chains = getattr(atoms, "chainIDs", None)
    if chains is None:
        chains = getattr(atoms, "segids", None)
    if chains is None:
        chains = [""] * len(labels)
    icodes = getattr(atoms, "icodes", None)
    if icodes is None:
        icodes = [""] * len(labels)
    resids, resnames = atoms.resids, atoms.resnames
    return [(str(chains[i]), int(resids[i]), str(icodes[i] or ""), str(resnames[i]))
            for i in range(len(labels))]


def conformer_choices(atoms) -> Dict[ResidueKey, Tuple[str, ...]]:
    """``{(chain, resid, icode, resname): kept labels}`` for residues with altlocs.

    Kept labels are every label whose mean occupancy is tied for the highest, in
    the order they first appear in the file. Usually one label; several when
    they are equally populated.
    """
    labels = _labels(atoms)
    if labels is None or not any(labels):
        return {}

    occupancies = getattr(atoms, "occupancies", None)
    keys = _residue_keys(atoms, labels)

    groups: "OrderedDict[ResidueKey, OrderedDict[str, List[float]]]" = OrderedDict()
    for i, label in enumerate(labels):
        if not label:
            continue
        occupancy = float(occupancies[i]) if occupancies is not None else 1.0
        groups.setdefault(keys[i], OrderedDict()).setdefault(label, []).append(occupancy)

    choices = {}
    for key, by_label in groups.items():
        mean = {label: sum(v) / len(v) for label, v in by_label.items()}
        top = max(mean.values())
        choices[key] = tuple(label for label in by_label
                             if top - mean[label] <= OCCUPANCY_TOLERANCE)
    return choices


def highest_occupancy_atoms(atoms):
    """``atoms`` without the alternate positions that are strictly less populated.

    Returns ``atoms`` unchanged when the topology carries no alternate-location
    information (most non-PDB formats), so it is always safe to apply.
    """
    labels = _labels(atoms)
    if labels is None or not any(labels):
        return atoms

    choices = conformer_choices(atoms)
    keys = _residue_keys(atoms, labels)
    keep = [i for i, label in enumerate(labels)
            if not label or label in choices.get(keys[i], ())]
    return atoms[keep]

"""Turn a folder of raw PDB files into a set WatCon can analyse together.

WatCon's input file starts at ``structure_directory: aligned_with_waters`` --
it assumes someone has already picked the right chain out of each structure,
superposed them all into one frame, and carried the waters along.  That step was
manual, needed MODELLER, and is where most of the real difficulty lives.

This module does it:

  1. **Find the right chain.**  Structures of one protein are routinely deposited
     as complexes (barnase with barstar) or fusions, and the protein of interest
     is not reliably chain 'A'.  Chains are scored against a reference.
  2. **Extract that chain plus its waters**, and relabel the chain, so one
     ConSurf file serves every structure with no chain_map.
  3. **Superpose** onto the reference with :mod:`WatCon.superpose`, moving the
     waters with the protein.
  4. **Reject loudly, with a reason.**  A structure that quietly vanishes from
     the dataset is worse than one that stops the run.

Scope: this matches residues by number, so it is for many structures of **one
protein** -- crystal forms, mutants, complexes.  Homologues numbered differently
need WatCon's MSA path.

Why residue number and not sequence position
--------------------------------------------
Scoring a chain by comparing sequence strings position-by-position is wrong in
both directions, and both errors are present in real PDB data.  Measured on the
22-structure barnase set:

* 1BSA starts at residue 4, 1B2S at 1 and 1RNB at 2, against a reference
  starting at 3.  Positionally they score 0.03-0.05 and are **rejected**, though
  by residue number they agree at 0.98-0.99.
* 2F56 and 2F5M start at residue 1 with the residue everyone else numbers 3.
  Positionally they score a perfect **1.00 and are accepted**, but their
  numbering is shifted by two, so every ConSurf score would land on the wrong
  residue.  They superpose at 6.05 A.

Residue number is also what the ConSurf join uses, so comparing the way the join
compares gets both cases right.
"""

from __future__ import annotations

import os
import shutil
import tempfile
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .residue_index import StructureResidue, residues_from_pdb_file
from .structure_io import as_pdb_directory, needs_conversion
from .superpose import kabsch, rmsd

__all__ = [
    "PreparationError",
    "PreparationReport",
    "StructureOutcome",
    "detect_offset",
    "extract_chain_with_waters",
    "find_target_chain",
    "identity_by_resid",
    "prepare_directory",
]

#: Minimum agreement with the reference for a chain to be accepted.  Generous,
#: because point mutants genuinely differ -- but far above what an unrelated
#: chain in a complex could reach (barstar against barnase scores ~0.03).
DEFAULT_MIN_IDENTITY = 0.80

#: Waters within this distance of any retained protein atom travel with it.
DEFAULT_WATER_CUTOFF = 5.0

#: Fewest shared residue numbers that can define a superposition worth trusting.
DEFAULT_MIN_SHARED = 30

_WATER_NAMES = ("HOH", "WAT", "SOL", "TIP3", "TIP4", "H2O")


class PreparationError(ValueError):
    """Raised when a dataset cannot be prepared at all."""


@dataclass
class StructureOutcome:
    """What happened to one structure."""

    pdb_id: str
    chain: Optional[str] = None
    identity: float = 0.0
    n_protein_atoms: int = 0
    n_waters: int = 0
    n_shared: int = 0
    rmsd: Optional[float] = None
    status: str = ""
    note: str = ""

    @property
    def ok(self) -> bool:
        return self.status == "ok"


@dataclass
class PreparationReport:
    """Outcome for a whole directory."""

    reference: str = ""
    reference_chain: str = ""
    outcomes: List[StructureOutcome] = field(default_factory=list)

    @property
    def prepared(self) -> List[StructureOutcome]:
        return [o for o in self.outcomes if o.ok]

    @property
    def rejected(self) -> List[StructureOutcome]:
        return [o for o in self.outcomes if not o.ok]

    @property
    def n_waters(self) -> int:
        return sum(o.n_waters for o in self.prepared)

    def describe(self) -> str:
        lines = ["reference: %s chain %s" % (self.reference, self.reference_chain)]
        for outcome in self.outcomes:
            if outcome.ok:
                lines.append(
                    "  %-8s chain %s  identity %.2f  %4d waters  RMSD %.2f A over %d CA"
                    % (outcome.pdb_id, outcome.chain, outcome.identity,
                       outcome.n_waters, outcome.rmsd, outcome.n_shared)
                )
            else:
                lines.append(
                    "  %-8s REJECTED  %s%s"
                    % (outcome.pdb_id, outcome.status,
                       "  [" + outcome.note + "]" if outcome.note else "")
                )
        lines.append(
            "prepared %d/%d structures, %d waters"
            % (len(self.prepared), len(self.outcomes), self.n_waters)
        )
        return "\n".join(lines)

    def write_csv(self, path) -> int:
        import csv

        columns = ["pdb_id", "chain", "identity", "n_protein_atoms", "n_waters",
                   "n_shared", "rmsd", "status", "note"]
        with open(path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            for outcome in self.outcomes:
                writer.writerow({
                    "pdb_id": outcome.pdb_id,
                    "chain": outcome.chain or "",
                    "identity": round(outcome.identity, 3),
                    "n_protein_atoms": outcome.n_protein_atoms,
                    "n_waters": outcome.n_waters,
                    "n_shared": outcome.n_shared,
                    "rmsd": "" if outcome.rmsd is None else round(outcome.rmsd, 3),
                    "status": outcome.status,
                    "note": outcome.note,
                })
        return len(self.outcomes)


# ---------------------------------------------------------------------------
# Comparing a chain to the reference
# ---------------------------------------------------------------------------

def identity_by_resid(
    residues: Sequence[StructureResidue], reference: Dict[int, str]
) -> float:
    """Fraction of the reference whose residue number carries the same residue.

    Keyed on residue number, not sequence position -- see the module docstring
    for the two real failures that distinction prevents.

    The denominator is the reference length, so a chain that overlaps only a
    fragment of the reference cannot score highly merely by agreeing on the few
    residues it happens to contain.
    """
    if not reference:
        return 0.0
    same = sum(
        1 for r in residues
        if r.resid in reference and r.one_letter == reference[r.resid]
    )
    return same / len(reference)


def detect_offset(
    residues: Sequence[StructureResidue], reference: Dict[int, str], span: int = 10
) -> Tuple[int, float]:
    """Best integer renumbering as ``(offset, identity)`` -- diagnostic only.

    Reported so a rejection can say *why*, instead of leaving a perfectly good
    structure looking like an unrelated protein.  Deliberately **never applied**:
    silently renumbering a deposited structure is the same helpful guess this
    package exists to stamp out.  If the offset is real, the user should say so.
    """
    best = (0, identity_by_resid(residues, reference))
    for offset in range(-span, span + 1):
        if offset == 0:
            continue
        shifted = [
            StructureResidue(r.chain, r.resid + offset, r.icode, r.resname,
                             r.one_letter)
            for r in residues
        ]
        score = identity_by_resid(shifted, reference)
        if score > best[1]:
            best = (offset, score)
    return best


def find_target_chain(
    path: str,
    reference: Dict[int, str],
    min_identity: float = DEFAULT_MIN_IDENTITY,
) -> Tuple[Optional[str], float, str]:
    """``(chain, identity, note)`` for the chain best matching the reference.

    ``chain`` is None when nothing clears ``min_identity``; ``note`` then carries
    a diagnosis when a simple renumbering would have matched.
    """
    residues = residues_from_pdb_file(path)
    if not residues:
        return None, 0.0, "no protein residues found"

    best_chain, best_score = None, 0.0
    for chain in sorted({r.chain for r in residues}):
        score = identity_by_resid([r for r in residues if r.chain == chain],
                                  reference)
        if score > best_score:
            best_chain, best_score = chain, score

    if best_score >= min_identity:
        return best_chain, best_score, ""

    for chain in sorted({r.chain for r in residues}):
        offset, score = detect_offset(
            [r for r in residues if r.chain == chain], reference
        )
        if offset and score >= min_identity:
            return None, best_score, (
                "chain %s matches at %.2f if renumbered by %+d -- renumber it "
                "yourself if that is correct; this is not applied automatically"
                % (chain, score, offset)
            )
    return None, best_score, ""


# ---------------------------------------------------------------------------
# Extraction and output
# ---------------------------------------------------------------------------

def _atom_lines(path: str):
    """(line, chain, resname, xyz) for every coordinate record of the first model."""
    rows = []
    with open(path, "r", errors="replace") as handle:
        for line in handle:
            if line.startswith("ENDMDL"):
                break                       # first model only
            if not line.startswith(("ATOM", "HETATM")) or len(line) < 54:
                continue
            try:
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            rows.append((line.rstrip("\n"), line[21], line[17:20].strip().upper(), xyz))
    return rows


def extract_chain_with_waters(
    path: str, chain: str, water_cutoff: float = DEFAULT_WATER_CUTOFF
):
    """``(protein, waters)`` line/coordinate pairs for one chain.

    Waters are selected by proximity rather than by chain label: crystallographic
    waters are frequently deposited under a chain of their own, or under the
    label of whichever molecule they were nearest during refinement, so
    filtering them by chain silently discards most of them.
    """
    protein, candidates = [], []
    for line, chain_id, resname, xyz in _atom_lines(path):
        if resname in _WATER_NAMES:
            candidates.append((line, xyz))
        elif chain_id == chain:
            protein.append((line, xyz))

    if not protein:
        return [], []

    coordinates = np.array([xyz for _, xyz in protein])
    waters = [
        (line, xyz) for line, xyz in candidates
        if np.min(np.linalg.norm(coordinates - np.array(xyz), axis=1)) <= water_cutoff
    ]
    return protein, waters


def _transform_and_relabel(lines, rotation, translation, chain_label):
    """Apply the rigid transform and rewrite the chain, in fixed PDB columns."""
    out = []
    for line, xyz in lines:
        moved = np.asarray(xyz) @ rotation.T + translation
        padded = line.ljust(80)
        out.append(
            padded[:21] + chain_label + padded[22:30]
            + "%8.3f%8.3f%8.3f" % (moved[0], moved[1], moved[2])
            + padded[54:80]
        )
    return out


def _ca_coordinates(lines):
    """``(resid, icode) -> CA coordinate`` from already-extracted lines."""
    result = {}
    for line, xyz in lines:
        if line[12:16].strip() != "CA":
            continue
        try:
            resid = int(line[22:26].strip())
        except ValueError:
            continue
        result.setdefault((resid, line[26:27].strip() or None), np.array(xyz))
    return result


# ---------------------------------------------------------------------------
# The whole directory
# ---------------------------------------------------------------------------

def prepare_directory(
    input_dir: str,
    out_dir: str,
    reference: Optional[str] = None,
    chain_label: str = "A",
    min_identity: float = DEFAULT_MIN_IDENTITY,
    water_cutoff: float = DEFAULT_WATER_CUTOFF,
    min_shared: int = DEFAULT_MIN_SHARED,
    verbose: bool = True,
) -> PreparationReport:
    """Prepare every structure in ``input_dir`` into a shared frame in ``out_dir``.

    Accepts PDB, mmCIF and gzipped forms of either. mmCIF is converted to PDB
    text first (see :mod:`WatCon.structure_io`) because RCSB no longer issues
    PDB files for large or recent entries -- 31 of the 287 PTP1B structures we
    selected have no PDB file at all -- and neither this module nor MDAnalysis
    reads mmCIF.

    The conversion happens in a temporary directory that is removed afterwards,
    so nothing downstream sees anything but PDB and ``input_dir`` is untouched.

    Every parameter is passed through unchanged; see
    :func:`_prepare_pdb_directory` for the full description.
    """
    if not os.path.isdir(input_dir):
        raise PreparationError("no such directory: %r" % (input_dir,))

    if not any(needs_conversion(name) or name.lower().endswith(".gz")
               for name in os.listdir(input_dir)):
        return _prepare_pdb_directory(
            input_dir, out_dir, reference, chain_label, min_identity,
            water_cutoff, min_shared, verbose,
        )

    work = tempfile.mkdtemp(prefix="watcon_convert_")
    try:
        readable, _converted, _failed = as_pdb_directory(input_dir, work, verbose)
        return _prepare_pdb_directory(
            readable, out_dir, reference, chain_label, min_identity,
            water_cutoff, min_shared, verbose,
        )
    finally:
        shutil.rmtree(work, ignore_errors=True)


def _prepare_pdb_directory(
    input_dir: str,
    out_dir: str,
    reference: Optional[str] = None,
    chain_label: str = "A",
    min_identity: float = DEFAULT_MIN_IDENTITY,
    water_cutoff: float = DEFAULT_WATER_CUTOFF,
    min_shared: int = DEFAULT_MIN_SHARED,
    verbose: bool = True,
) -> PreparationReport:
    """Prepare every PDB in ``input_dir`` into a shared frame in ``out_dir``.

    Parameters
    ----------
    input_dir : str
        Folder of raw ``.pdb`` files.
    out_dir : str
        Destination.  **Cleared, not merged** -- see the note below.
    reference : str, optional
        Filename or PDB id to superpose onto.  Defaults to the first in sorted
        order.  Its largest chain defines the reference sequence.
    chain_label : str
        Chain letter every output is relabelled to, so one ConSurf file serves
        the whole set without a chain map.
    min_identity, water_cutoff, min_shared
        See the module-level defaults.

    Returns
    -------
    PreparationReport

    Note
    ----
    ``out_dir`` is emptied first.  Merging would let a structure rejected on this
    run survive from a previous one -- which happened during the barnase study,
    where a stale copy of a mis-numbered structure reached the analysis and was
    caught only by :func:`WatCon.evolutionary.enforce_identity`.
    """
    if not os.path.isdir(input_dir):
        raise PreparationError("no such directory: %r" % (input_dir,))

    names = sorted(f for f in os.listdir(input_dir) if f.lower().endswith(".pdb"))
    if not names:
        raise PreparationError("no .pdb files in %r" % (input_dir,))

    if reference is None:
        reference_file = names[0]
    else:
        candidates = [reference, reference + ".pdb", os.path.basename(reference)]
        reference_file = next((c for c in candidates if c in names), None)
        if reference_file is None:
            raise PreparationError(
                "reference %r is not among the structures (%s)"
                % (reference, ", ".join(names))
            )

    reference_path = os.path.join(input_dir, reference_file)
    reference_residues = residues_from_pdb_file(reference_path)
    if not reference_residues:
        raise PreparationError("no protein residues in reference %r" % reference_file)

    # Largest chain, ties broken alphabetically.  The sort matters: `max` over a
    # set returns whichever equal-length chain iteration happens to reach first,
    # and set order for strings varies with PYTHONHASHSEED -- so without it, a
    # structure with several identical copies (1A2P has three) could pick a
    # different reference chain on different runs and shift every coordinate.
    reference_chain = max(
        sorted({r.chain for r in reference_residues}),
        key=lambda c: len([r for r in reference_residues if r.chain == c]),
    )
    reference_map = {
        r.resid: r.one_letter for r in reference_residues if r.chain == reference_chain
    }

    report = PreparationReport(
        reference=reference_file, reference_chain=reference_chain
    )
    if verbose:
        print("Preparing %d structures (reference %s chain %s, %d residues)"
              % (len(names), reference_file, reference_chain, len(reference_map)))

    reference_protein, _ = extract_chain_with_waters(
        reference_path, reference_chain, water_cutoff
    )
    reference_cas = _ca_coordinates(reference_protein)

    # Cleared, never merged -- see the docstring.
    if os.path.isdir(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)

    for name in names:
        pdb_id = os.path.splitext(name)[0]
        path = os.path.join(input_dir, name)
        outcome = StructureOutcome(pdb_id=pdb_id)

        chain, score, note = find_target_chain(path, reference_map, min_identity)
        outcome.identity, outcome.note = score, note

        if chain is None:
            outcome.status = "identity %.2f by residue number, below %.2f" % (
                score, min_identity)
            report.outcomes.append(outcome)
            if verbose:
                print("  %-8s REJECTED (%s)%s" % (
                    pdb_id, outcome.status, "  [" + note + "]" if note else ""))
            continue

        protein, waters = extract_chain_with_waters(path, chain, water_cutoff)
        cas = _ca_coordinates(protein)
        shared = sorted(set(reference_cas) & set(cas))

        if len(shared) < min_shared:
            outcome.chain = chain
            outcome.status = "only %d residue numbers shared with the reference" % len(shared)
            report.outcomes.append(outcome)
            if verbose:
                print("  %-8s REJECTED (%s)" % (pdb_id, outcome.status))
            continue

        rotation, translation = kabsch(
            np.array([cas[k] for k in shared]),
            np.array([reference_cas[k] for k in shared]),
        )
        fit = rmsd(
            np.array([cas[k] for k in shared]) @ rotation.T + translation,
            np.array([reference_cas[k] for k in shared]),
        )

        lines = _transform_and_relabel(protein, rotation, translation, chain_label)
        lines += _transform_and_relabel(waters, rotation, translation, chain_label)
        with open(os.path.join(out_dir, pdb_id + ".pdb"), "w") as handle:
            handle.write("\n".join(lines) + "\nEND\n")

        outcome.chain = chain
        outcome.n_protein_atoms = len(protein)
        outcome.n_waters = len(waters)
        outcome.n_shared = len(shared)
        outcome.rmsd = fit
        outcome.status = "ok"
        report.outcomes.append(outcome)
        if verbose:
            print("  %-8s chain %s  identity %.2f  %4d waters  RMSD %.2f A over %d CA"
                  % (pdb_id, chain, score, len(waters), fit, len(shared)))

    if not report.prepared:
        raise PreparationError(
            "no structure could be prepared. Every one was rejected -- check that "
            "they are the same protein and that the reference is right.\n"
            + report.describe()
        )

    if verbose:
        print("Prepared %d/%d structures, %d waters -> %s"
              % (len(report.prepared), len(report.outcomes), report.n_waters, out_dir))
    return report

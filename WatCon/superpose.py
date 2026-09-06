"""Rigid-body superposition without MODELLER.

WatCon's existing alignment path runs through
:func:`WatCon.sequence_processing.perform_structure_alignment`, which needs
MODELLER -- licensed software that cannot be pip-installed.  Everything
downstream of it is fine:
:func:`WatCon.sequence_processing.align_with_waters` already applies rotation
and translation matrices to *every* atom, waters included, and can write out a
chosen chain plus its non-protein atoms.  It just needs the matrices.

This module produces them with numpy, so a user without a MODELLER licence can
superpose structures of the same protein and cluster their waters together.

Scope, stated plainly: this matches residues by **residue number**, so it is for
structures of the same protein -- different crystal forms, mutants, ligand
complexes.  It is not a sequence aligner and will not superpose homologues whose
numbering differs.  For that, WatCon's MSA path remains the right tool.

The output convention matches ``perform_structure_alignment`` exactly -- a
``{'Rot': [...], 'Trans': [...]}`` dict holding one entry per *non-reference*
structure, in sorted filename order -- so the two are interchangeable as inputs
to ``align_with_waters``.
"""

from __future__ import annotations

import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .residue_index import (
    _ATOMNAME,
    _CHAIN,
    _ICODE,
    _RESNAME,
    _RESSEQ,
    ResidueKey,
    _resolve_table,
)

__all__ = [
    "SuperpositionResult",
    "ca_coordinates",
    "kabsch",
    "rmsd",
    "superpose_structures",
]


class SuperpositionError(ValueError):
    """Raised when two structures cannot be meaningfully superposed."""


# ---------------------------------------------------------------------------
# Reading CA coordinates
# ---------------------------------------------------------------------------

def ca_coordinates(
    path: str,
    chain: Optional[str] = None,
    custom_residues: Optional[Dict[str, str]] = None,
) -> Dict[ResidueKey, np.ndarray]:
    """``(chain, resid, icode) -> CA coordinate`` for one structure.

    Follows exactly the conventions
    :func:`WatCon.residue_index.residues_from_pdb_file` uses -- fixed PDB
    columns rather than substring matching, an atom named exactly ``CA``, and
    the first alternate conformer of a residue rather than a whitelist of altloc
    labels.  The column slices are imported from that module rather than
    restated here, so the two cannot drift apart.
    """
    table = _resolve_table(custom_residues)
    coordinates: Dict[ResidueKey, np.ndarray] = {}

    with open(path, "r", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if len(line) < 54:
                continue
            if line[_ATOMNAME].strip() != "CA":
                continue
            if line[_RESNAME].strip().upper() not in table:
                continue

            chain_id = line[_CHAIN]
            if chain is not None and chain_id != chain:
                continue

            try:
                resid = int(line[_RESSEQ].strip())
                xyz = (
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                )
            except ValueError:
                continue

            key = (chain_id, resid, line[_ICODE].strip() or None)
            if key in coordinates:      # first conformer wins, as in residue_index
                continue
            coordinates[key] = np.array(xyz, dtype=float)

    return coordinates


# ---------------------------------------------------------------------------
# The algorithm
# ---------------------------------------------------------------------------

def kabsch(mobile: np.ndarray, target: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Optimal rigid transform taking ``mobile`` onto ``target``.

    Returns ``(R, t)`` such that ``mobile @ R.T + t`` is the least-squares best
    fit to ``target`` -- the same convention
    :func:`WatCon.sequence_processing.align_with_waters` applies.

    Both arrays are ``(N, 3)`` and must be in corresponding order.  The
    reflection correction matters: without it the SVD can return an improper
    rotation (determinant -1), which superposes the coordinates onto a *mirror
    image* of the target.  For a protein that is not a rigid motion, and it
    would silently produce a structure of the wrong handedness.
    """
    mobile = np.asarray(mobile, dtype=float)
    target = np.asarray(target, dtype=float)
    if mobile.shape != target.shape:
        raise SuperpositionError(
            "mobile and target must have the same shape, got %s and %s"
            % (mobile.shape, target.shape)
        )
    if mobile.ndim != 2 or mobile.shape[1] != 3:
        raise SuperpositionError("expected (N, 3) coordinates, got %s" % (mobile.shape,))
    if len(mobile) < 3:
        raise SuperpositionError(
            "at least 3 corresponding atoms are needed to define a rotation, got %d"
            % len(mobile)
        )

    mobile_centre = mobile.mean(axis=0)
    target_centre = target.mean(axis=0)

    covariance = (mobile - mobile_centre).T @ (target - target_centre)
    u, _, vt = np.linalg.svd(covariance)

    # Reflection correction -- see the docstring.
    correction = np.diag([1.0, 1.0, np.sign(np.linalg.det(vt.T @ u.T))])
    rotation = vt.T @ correction @ u.T

    translation = target_centre - rotation @ mobile_centre
    return rotation, translation


def rmsd(mobile: np.ndarray, target: np.ndarray) -> float:
    """Root-mean-square deviation between two already-superposed point sets."""
    mobile = np.asarray(mobile, dtype=float)
    target = np.asarray(target, dtype=float)
    return float(np.sqrt(np.mean(np.sum((mobile - target) ** 2, axis=1))))


# ---------------------------------------------------------------------------
# Whole structures
# ---------------------------------------------------------------------------

class SuperpositionResult:
    """Transforms for a set of structures, plus what they were derived from.

    ``rotations``/``translations`` are ordered to match
    ``perform_structure_alignment``: one entry per structure *after* the
    reference, in the same sorted order, so they can be handed straight to
    ``align_with_waters``.
    """

    def __init__(self, reference: str, names: Sequence[str]):
        self.reference = reference
        self.names: List[str] = list(names)
        self.rotations: List[np.ndarray] = []
        self.translations: List[np.ndarray] = []
        #: name -> (RMSD after superposition, number of atoms matched)
        self.quality: Dict[str, Tuple[float, int]] = {}

    def as_dict(self) -> Dict[str, List[np.ndarray]]:
        """The ``{'Rot': [...], 'Trans': [...]}`` shape WatCon already uses.

        Safe to hand straight to
        :func:`WatCon.sequence_processing.align_with_waters`.

        That function sorts the directory and *always* treats the first file as
        the reference, applying ``rotation_matrices[i - 1]`` to the i-th.  So a
        reference that is not alphabetically first would pair every structure
        with the wrong transform -- and the output would still look like a
        superposition, just a wrong one.  Refused here rather than discovered
        later in a water cluster that does not exist.

        Note ``align_with_waters`` lists the *whole* directory, so keep it free
        of non-PDB files or its ordering will not match this one.
        """
        first = sorted(self.names)[0]
        if self.reference != first:
            raise SuperpositionError(
                "reference is %r but %r sorts first. align_with_waters assumes "
                "the sorted-first structure is the reference, so these "
                "transforms would be applied to the wrong files. Either use %r "
                "as the reference, or apply .rotations/.translations yourself "
                "in the order given by .names."
                % (self.reference, first, first)
            )
        return {"Rot": list(self.rotations), "Trans": list(self.translations)}

    def describe(self) -> str:
        lines = ["reference: %s" % self.reference]
        for name in self.names:
            if name == self.reference:
                continue
            value, count = self.quality.get(name, (float("nan"), 0))
            lines.append("  %-24s RMSD %5.2f A over %4d CA" % (name, value, count))
        return "\n".join(lines)


def _correspondence(
    path: str, name: str, chain: Optional[str]
) -> Dict[Tuple[int, Optional[str]], np.ndarray]:
    """CA coordinates keyed by ``(resid, icode)`` -- deliberately *without* the chain.

    A chain letter is a label a depositor chose, not part of a residue's
    identity across structures.  The same protein is chain A in one entry and
    chain B in another, and the three barnase copies inside 1BRS are A, B and C.
    Keying on the chain would make those share zero residues and appear
    un-superposable, which is how this was first found.

    Dropping the chain is only safe once at most one chain is in play, so a
    collision is reported rather than silently resolved.
    """
    coordinates = ca_coordinates(path, chain=chain)
    reduced: Dict[Tuple[int, Optional[str]], np.ndarray] = {}
    for (chain_id, resid, icode), xyz in coordinates.items():
        key = (resid, icode)
        if key in reduced:
            raise SuperpositionError(
                "%s contains residue %d%s in more than one chain. Pass chain= to "
                "say which copy to superpose on; without it the correspondence "
                "between structures is ambiguous."
                % (name, resid, icode or "")
            )
        reduced[key] = xyz
    return reduced


def superpose_structures(
    pdb_dir,
    reference: Optional[str] = None,
    chain: Optional[str] = "A",
    min_atoms: int = 3,
    max_rmsd: Optional[float] = None,
) -> SuperpositionResult:
    """Superpose every structure in ``pdb_dir`` onto one reference.

    Parameters
    ----------
    pdb_dir : str or list of str
        Directory of PDB files, or an explicit list of paths.
    reference : str, optional
        Filename (not path) to align onto.  Defaults to the first in sorted
        order, matching ``perform_structure_alignment``.
    chain : str or None
        Chain to superpose on.  None uses every chain present.
    min_atoms : int
        Refuse a structure sharing fewer than this many residue numbers with the
        reference.  Guards against "aligned" on a handful of incidental matches.
    max_rmsd : float, optional
        If given, refuse a structure whose post-fit RMSD exceeds it.

    Returns
    -------
    SuperpositionResult

    Raises
    ------
    SuperpositionError
        If the reference is missing, or a structure shares too few residues.
        Failing loudly is deliberate: a silently skipped structure would drop
        out of the water clustering without appearing anywhere in the results.

    Notes
    -----
    Correspondence is by residue number, so **this assumes one protein**.  Only
    residues present in both structures are used, which handles differing
    disorder at the termini.  Mutated positions still contribute their CA -- a
    side-chain substitution does not move the backbone enough to matter, and
    excluding them would bias the fit toward the unmutated regions.
    """
    if isinstance(pdb_dir, (list, tuple)):
        paths = {os.path.basename(p): p for p in pdb_dir}
    else:
        paths = {
            f: os.path.join(pdb_dir, f)
            for f in os.listdir(pdb_dir)
            if f.lower().endswith(".pdb")
        }

    names = sorted(paths)
    if not names:
        raise SuperpositionError("no PDB files found in %r" % (pdb_dir,))

    if reference is None:
        reference = names[0]
    elif reference not in paths:
        raise SuperpositionError(
            "reference %r not among the structures (%s)"
            % (reference, ", ".join(names))
        )

    # Put the reference first, so the transform lists line up with the order
    # align_with_waters walks the directory in.
    ordered = [reference] + [n for n in names if n != reference]

    result = SuperpositionResult(reference, ordered)
    ref_coords = _correspondence(paths[reference], reference, chain)
    if not ref_coords:
        raise SuperpositionError(
            "no CA atoms found in reference %r (chain=%r)" % (reference, chain)
        )
    result.quality[reference] = (0.0, len(ref_coords))

    for name in ordered[1:]:
        coords = _correspondence(paths[name], name, chain)
        shared = sorted(set(ref_coords) & set(coords))
        if len(shared) < min_atoms:
            raise SuperpositionError(
                "%s shares only %d residue numbers with reference %s "
                "(need %d). Are these the same protein, numbered the same way?"
                % (name, len(shared), reference, min_atoms)
            )

        mobile = np.array([coords[k] for k in shared])
        target = np.array([ref_coords[k] for k in shared])
        rotation, translation = kabsch(mobile, target)

        fitted = mobile @ rotation.T + translation
        value = rmsd(fitted, target)
        if max_rmsd is not None and value > max_rmsd:
            raise SuperpositionError(
                "%s superposes onto %s at RMSD %.2f A, above the %.2f A limit"
                % (name, reference, value, max_rmsd)
            )

        result.rotations.append(rotation)
        result.translations.append(translation)
        result.quality[name] = (value, len(shared))

    return result

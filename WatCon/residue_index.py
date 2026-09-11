"""Explicit residue identity and MSA-column mapping.

Why this module exists
----------------------

``sequence_processing.generate_msa_alignment`` returns ``msa_indices``, a
**positional** list: ``msa_indices[k]`` is the alignment column of the *k-th
non-gap residue* of a sequence.  Several places in WatCon index it with an
absolute PDB residue number instead::

    msa_resid = msa_indexing[atm.resid - 1]

That is only correct when residue numbering is 1-origin, gapless and
single-chain.  Real structures break all three, and the failure is silent:

* 1BRS chain A starts at residue 3   -> off by 2
* P00648 is UniProt-numbered from 48 -> off by 47
* 7O7W carries an expression tag numbered -5 .. 0, and ``msa_indices[-6]``
  wraps around to the *end* of the list without raising

This module replaces the arithmetic with an explicit, validated lookup::

    (chain, resid, icode)  ->  ordinal (sequence order)  ->  MSA column

The second reason it exists is consistency.  The FASTA that feeds the alignment
and the residue list that consumes it are currently produced by two independent
walks -- raw PDB text parsing in ``sequence_processing.pdb_to_fastas`` versus an
MDAnalysis selection in the network builders -- and nothing checks that they
agree.  :meth:`ResidueIndex.build` refuses to construct a mapping when the two
lengths differ, turning a silent shift into an immediate error.

Dependencies
------------

The core of this module is standard library only, so it is testable without
MDAnalysis.  :func:`residues_from_universe` is a thin adapter that imports
MDAnalysis lazily, and is the only part that needs it.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

__all__ = [
    "AMINO_ACID_3TO1",
    "ResidueIndexError",
    "ResidueKey",
    "StructureResidue",
    "ResidueIndex",
    "SequenceConsistencyReport",
    "residues_from_pdb_file",
    "residues_from_universe",
    "atom_identity",
    "build_fasta_sequence",
    "check_sequence_consistency",
    "legacy_indexing_report",
]

#: ``(chain, resid, icode)`` -- the identity used everywhere in this module and
#: the key ConSurf records are joined on.
ResidueKey = Tuple[str, int, Optional[str]]


class ResidueIndexError(ValueError):
    """Raised when residue identity or the MSA mapping cannot be trusted."""


# ---------------------------------------------------------------------------
# Residue name table
# ---------------------------------------------------------------------------
#
# Kept identical to the table in ``sequence_processing.pdb_to_fastas`` so the
# two produce the same sequence.  It lives here rather than there because this
# module is the lower-level one and, unlike sequence_processing, does not import
# modeller at module scope.

AMINO_ACID_3TO1: Dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
    "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "HIP": "H", "HID": "H", "HIE": "H", "HISD": "H",
    "HISE": "H", "HISP": "H", "AS4": "D", "ASH": "D",
    "GL4": "E", "GLH": "E", "ARN": "R", "LYN": "K",
    "CYX": "C", "CYM": "C", "CSP": "C", "SEP": "S", "ASX": "D",
    "HSD": "H", "HSP": "H", "HSE": "H", "MSE": "M",
}


def _resolve_table(custom_residues: Optional[Dict[str, str]]) -> Dict[str, str]:
    table = dict(AMINO_ACID_3TO1)
    if custom_residues:
        table.update({k.upper(): v for k, v in custom_residues.items()})
    return table


# ---------------------------------------------------------------------------
# One residue
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class StructureResidue:
    """One protein residue, with its identity kept explicit.

    ``resid`` may be negative or zero -- expression tags are numbered backwards
    from the mature start -- so never use it as a list index.
    """

    chain: str
    resid: int
    icode: Optional[str]
    resname: str
    one_letter: str

    @property
    def key(self) -> ResidueKey:
        return (self.chain, self.resid, self.icode)

    def __str__(self) -> str:
        return f"{self.resname}:{self.resid}{self.icode or ''}:{self.chain}"


# ---------------------------------------------------------------------------
# Choosing which files in a directory are structures
# ---------------------------------------------------------------------------

#: Extensions treated as structures. Deliberately explicit: the previous rule
#: was "everything in the directory that does not contain 'swp'", which handed
#: READMEs, subdirectories and stray output files to the structure reader.
STRUCTURE_SUFFIXES = (
    ".pdb", ".ent", ".cif", ".mmcif", ".gro", ".pdbqt", ".xyz", ".mol2",
    ".prmtop", ".parm7", ".psf", ".top", ".tpr",
)


class NoStructuresFound(ValueError):
    """Raised when a directory contains nothing that can be read."""


def list_structure_files(directory, suffixes=STRUCTURE_SUFFIXES):
    """``(structure filenames, skipped filenames)`` for one directory.

    Anything that is not a regular file with a known structure extension is
    skipped and reported rather than passed to the reader. A user who leaves a
    README, a results folder or a ``.DS_Store`` beside their structures
    previously got ``ValueError: 'TXT' isn't a valid topology format`` from deep
    inside MDAnalysis, which names neither the file nor the fix.

    Compressed structures (``.pdb.gz``) are recognised by their inner extension.
    """
    import os

    structures, skipped = [], []
    for name in sorted(os.listdir(directory)):
        path = os.path.join(directory, name)
        if not os.path.isfile(path):
            skipped.append(name)
            continue
        stem = name[:-3] if name.lower().endswith(".gz") else name
        if stem.lower().endswith(suffixes):
            structures.append(name)
        else:
            skipped.append(name)
    return structures, skipped


def require_structure_files(directory, suffixes=STRUCTURE_SUFFIXES):
    """Like :func:`list_structure_files`, but refuses to return nothing.

    An empty result used to surface much later as
    ``ValueError: not enough values to unpack (expected 2, got 0)``.
    """
    import os

    if not os.path.isdir(directory):
        # %s, not %r: repr of a Windows path doubles every backslash, so
        # the message shows a path the user cannot paste back.
        raise NoStructuresFound(
            "No such directory: %s. This should hold the structure files to "
            "analyse." % (directory,)
        )
    structures, skipped = list_structure_files(directory, suffixes)
    if not structures:
        detail = ""
        if skipped:
            shown = ", ".join(skipped[:6])
            detail = (" It contains %d item(s) that are not structures (%s%s)."
                      % (len(skipped), shown, ", ..." if len(skipped) > 6 else ""))
        raise NoStructuresFound(
            "No structure files in %s.%s Recognised extensions: %s."
            % (directory, detail, ", ".join(suffixes))
        )
    return structures, skipped


# ---------------------------------------------------------------------------
# Reading residues from a structure
# ---------------------------------------------------------------------------

# PDB column slices (0-based, end-exclusive) per the PDB format specification.
_ALTLOC = slice(16, 17)
_RESNAME = slice(17, 20)
_CHAIN = slice(21, 22)
_RESSEQ = slice(22, 26)
_ICODE = slice(26, 27)
_ATOMNAME = slice(12, 16)

# Alternate conformers are handled by keeping the FIRST occurrence of each
# residue key, whatever its altloc label.  An earlier version whitelisted
# {"", "A", "1"} and silently dropped 11 residues of 7O7W whose only conformers
# are labelled B and C -- the same class of silent loss this module exists to
# remove.  Never assume a residue has a blank or "A" conformer.


def residues_from_pdb_file(
    path: Union[str, Path],
    chain: Optional[str] = None,
    custom_residues: Optional[Dict[str, str]] = None,
    include_hetatm: bool = False,
) -> List[StructureResidue]:
    """Ordered protein residues read straight from a PDB file.

    One residue per alpha-carbon, in file order.  This is the reference walk:
    the FASTA and the MSA index must both be derived from it.

    Unlike the historical text scan in ``sequence_processing.pdb_to_fastas``
    (``('ATOM' in line) and ('CA' in line)``), this parses fixed PDB columns, so
    it does not match ``HETATM`` records by substring, does not match a stray
    ``CA`` elsewhere on the line, collapses alternate conformers to one residue,
    and can restrict to one chain.

    Parameters
    ----------
    path:
        PDB file to read.
    chain:
        Keep only this chain.  ``None`` keeps every chain, in file order.
    custom_residues:
        Extra three-letter to one-letter mappings.
    include_hetatm:
        Also consider ``HETATM`` records, for modified residues that are part of
        the chain (``MSE``, for example).  Off by default.

    Returns
    -------
    list of StructureResidue

    Limitation
    ----------
    A residue is recognised by an atom named exactly ``CA``.  Fused or heavily
    modified groups whose alpha carbons carry other names cannot be picked up
    this way -- 7O7W's chromophore ``PIA:68:A`` has ``CA1``/``CA2``/``CA3`` and
    is skipped even with ``include_hetatm`` and a custom name.  ConSurf omits
    that residue too, so the two remain consistent, but anything that needs such
    a residue must find it by another route.
    """
    table = _resolve_table(custom_residues)
    prefixes = ("ATOM  ", "HETATM") if include_hetatm else ("ATOM  ",)

    residues: List[StructureResidue] = []
    seen: set = set()

    text = Path(path).read_text(encoding="utf-8", errors="replace")

    for line in text.splitlines():
        if not line.startswith(prefixes):
            continue
        if len(line) < 27:
            continue
        if line[_ATOMNAME].strip() != "CA":
            continue

        resname = line[_RESNAME].strip().upper()
        if resname not in table:
            continue

        chain_id = line[_CHAIN]
        if chain is not None and chain_id != chain:
            continue

        try:
            resid = int(line[_RESSEQ].strip())
        except ValueError:
            continue

        icode = line[_ICODE].strip() or None
        key = (chain_id, resid, icode)
        if key in seen:
            continue
        seen.add(key)

        residues.append(
            StructureResidue(
                chain=chain_id,
                resid=resid,
                icode=icode,
                resname=resname,
                one_letter=table[resname],
            )
        )

    return residues


def residues_from_universe(
    universe,
    chain: Optional[str] = None,
    custom_residues: Optional[Dict[str, str]] = None,
    selection: str = "protein",
) -> List[StructureResidue]:
    """Ordered protein residues from an MDAnalysis Universe.

    The MDAnalysis counterpart of :func:`residues_from_pdb_file`, kept
    deliberately parallel so the two produce the same ordering.  MDAnalysis is
    imported lazily, so importing this module never requires it.

    Note
    ----
    Not exercised by the test suite in environments without MDAnalysis.  Verify
    against :func:`residues_from_pdb_file` on the same structure before relying
    on it -- :func:`check_sequence_consistency` is there for exactly that.
    """
    table = _resolve_table(custom_residues)

    residues: List[StructureResidue] = []
    for residue in universe.select_atoms(selection).residues:
        resname = str(residue.resname).strip().upper()
        if resname not in table:
            continue

        chain_id = _universe_chain_id(residue)
        if chain is not None and chain_id != chain:
            continue

        icode = str(getattr(residue, "icode", "") or "").strip() or None

        residues.append(
            StructureResidue(
                chain=chain_id,
                resid=int(residue.resid),
                icode=icode,
                resname=resname,
                one_letter=table[resname],
            )
        )

    return residues


def _universe_chain_id(residue) -> str:
    """Best available chain identifier for an MDAnalysis residue.

    PDB topologies expose ``chainID``; others only ``segid``.  Falls back to the
    empty string, which matches how a blank PDB chain column is read.
    """
    for attribute in ("chainID", "segid"):
        value = getattr(residue, attribute, None)
        if value is None:
            continue
        text = str(value).strip()
        if text and text != "SYSTEM":
            return text
    return ""


def atom_identity(atom) -> Tuple[str, Optional[str]]:
    """Chain identifier and insertion code for an MDAnalysis atom.

    Returns ``(chain, icode)``.  Together with ``atom.resid`` these form the
    residue identity used by :class:`ResidueIndex`.  Kept here so the network
    builders and the residue walk agree on how identity is read out of
    MDAnalysis -- if they diverged, lookups would miss silently.
    """
    residue = getattr(atom, "residue", atom)
    chain = _universe_chain_id(residue)
    icode = str(getattr(residue, "icode", "") or "").strip() or None
    return chain, icode


# ---------------------------------------------------------------------------
# Sequence
# ---------------------------------------------------------------------------

def build_fasta_sequence(residues: Sequence[StructureResidue]) -> str:
    """One-letter sequence for an ordered residue list."""
    return "".join(residue.one_letter for residue in residues)


@dataclass
class SequenceConsistencyReport:
    """Comparison between a residue walk and an existing FASTA sequence."""

    matches: bool
    structure_length: int
    fasta_length: int
    first_mismatch: Optional[int] = None
    structure_sequence: str = ""
    fasta_sequence: str = ""

    def describe(self) -> str:
        if self.matches:
            return f"consistent ({self.structure_length} residues)"
        if self.structure_length != self.fasta_length:
            return (
                f"length mismatch: structure has {self.structure_length} "
                f"residues, FASTA has {self.fasta_length}"
            )
        return (
            f"sequences differ from ordinal {self.first_mismatch} "
            f"(structure {self.structure_sequence[self.first_mismatch]!r} vs "
            f"FASTA {self.fasta_sequence[self.first_mismatch]!r})"
        )


def check_sequence_consistency(
    residues: Sequence[StructureResidue],
    fasta_sequence: str,
) -> SequenceConsistencyReport:
    """Check a residue walk against the FASTA that fed the alignment.

    This is the guard against the two walks drifting apart.  If they disagree,
    every MSA column derived from the alignment is shifted, and nothing
    downstream would notice.
    """
    structure_sequence = build_fasta_sequence(residues)
    clean = "".join(fasta_sequence.split()).upper()

    if len(structure_sequence) != len(clean):
        return SequenceConsistencyReport(
            matches=False,
            structure_length=len(structure_sequence),
            fasta_length=len(clean),
            structure_sequence=structure_sequence,
            fasta_sequence=clean,
        )

    for ordinal, (left, right) in enumerate(zip(structure_sequence, clean)):
        if left != right:
            return SequenceConsistencyReport(
                matches=False,
                structure_length=len(structure_sequence),
                fasta_length=len(clean),
                first_mismatch=ordinal,
                structure_sequence=structure_sequence,
                fasta_sequence=clean,
            )

    return SequenceConsistencyReport(
        matches=True,
        structure_length=len(structure_sequence),
        fasta_length=len(clean),
        structure_sequence=structure_sequence,
        fasta_sequence=clean,
    )


# ---------------------------------------------------------------------------
# The index
# ---------------------------------------------------------------------------

class ResidueIndex:
    """Validated mapping from residue identity to MSA alignment column.

    Build it with :meth:`build`, then look up by identity::

        index = ResidueIndex.build(residues, msa_indices)
        column = index.msa_column("A", -5)

    Lookups never guess.  An unknown residue returns ``None`` from
    :meth:`msa_column`, or raises from :meth:`require_msa_column`.
    """

    def __init__(
        self,
        residues: Sequence[StructureResidue],
        msa_indices: Optional[Sequence[int]] = None,
        source: Optional[str] = None,
    ) -> None:
        self._residues: List[StructureResidue] = list(residues)
        self._msa_indices: Optional[List[int]] = (
            None if msa_indices is None else [int(i) for i in msa_indices]
        )
        self.source = source

        self._ordinal_by_key: Dict[ResidueKey, int] = {}
        for ordinal, residue in enumerate(self._residues):
            self._ordinal_by_key[residue.key] = ordinal

    # -- construction ------------------------------------------------------

    @classmethod
    def build(
        cls,
        residues: Sequence[StructureResidue],
        msa_indices: Optional[Sequence[int]] = None,
        source: Optional[str] = None,
    ) -> "ResidueIndex":
        """Build an index, refusing to produce one that would be wrong.

        Raises
        ------
        ResidueIndexError
            If two residues share an identity, or if ``msa_indices`` is a
            different length from ``residues``.  The length check is the whole
            point: a mismatch means the FASTA that produced the alignment and
            this residue walk describe different things, so every column would
            be shifted.
        """
        keys = [residue.key for residue in residues]
        if len(keys) != len(set(keys)):
            duplicates = sorted({k for k in keys if keys.count(k) > 1})
            raise ResidueIndexError(
                f"Duplicate residue identities in structure: {duplicates[:5]}"
            )

        if msa_indices is not None and len(msa_indices) != len(residues):
            raise ResidueIndexError(
                f"Sequence/structure length mismatch: {len(residues)} residues "
                f"but {len(msa_indices)} MSA indices"
                + (f" (source: {source})" if source else "")
                + ". The FASTA used for the alignment does not describe this "
                "residue selection, so every MSA column would be shifted."
            )

        return cls(residues, msa_indices, source=source)

    @classmethod
    def from_pdb_file(
        cls,
        path: Union[str, Path],
        msa_indices: Optional[Sequence[int]] = None,
        chain: Optional[str] = None,
        custom_residues: Optional[Dict[str, str]] = None,
    ) -> "ResidueIndex":
        """Convenience: read a PDB and build the index in one step."""
        residues = residues_from_pdb_file(
            path, chain=chain, custom_residues=custom_residues
        )
        return cls.build(residues, msa_indices, source=str(path))

    # -- introspection -----------------------------------------------------

    def __len__(self) -> int:
        return len(self._residues)

    def __contains__(self, key: object) -> bool:
        return key in self._ordinal_by_key

    @property
    def residues(self) -> List[StructureResidue]:
        return list(self._residues)

    @property
    def msa_indices(self) -> Optional[List[int]]:
        return None if self._msa_indices is None else list(self._msa_indices)

    @property
    def has_msa(self) -> bool:
        return self._msa_indices is not None

    def chains(self) -> List[str]:
        seen: List[str] = []
        for residue in self._residues:
            if residue.chain not in seen:
                seen.append(residue.chain)
        return seen

    def sequence(self) -> str:
        return build_fasta_sequence(self._residues)

    # -- lookup ------------------------------------------------------------

    def ordinal(
        self, chain: str, resid: int, icode: Optional[str] = None
    ) -> Optional[int]:
        """Zero-based position of a residue in sequence order, or ``None``."""
        return self._ordinal_by_key.get((chain, resid, icode))

    def msa_column(
        self, chain: str, resid: int, icode: Optional[str] = None
    ) -> Optional[int]:
        """MSA alignment column for a residue, or ``None`` if not mappable.

        ``None`` means "this parser cannot say", never "column zero".  Callers
        that must have an answer should use :meth:`require_msa_column`.
        """
        if self._msa_indices is None:
            return None
        ordinal = self.ordinal(chain, resid, icode)
        if ordinal is None:
            return None
        return self._msa_indices[ordinal]

    def require_msa_column(
        self, chain: str, resid: int, icode: Optional[str] = None
    ) -> int:
        """Like :meth:`msa_column`, but raises rather than returning ``None``."""
        if self._msa_indices is None:
            raise ResidueIndexError("This index carries no MSA indices")
        column = self.msa_column(chain, resid, icode)
        if column is None:
            raise ResidueIndexError(
                f"No MSA column for residue ({chain!r}, {resid}, {icode!r})"
            )
        return column

    def residue_at(self, ordinal: int) -> StructureResidue:
        return self._residues[ordinal]

    def keys(self) -> List[ResidueKey]:
        return [residue.key for residue in self._residues]

    # -- diagnostics -------------------------------------------------------

    def legacy_indexing_report(self) -> "LegacyIndexingReport":
        """Quantify what the old ``msa_indices[resid - 1]`` arithmetic would do."""
        return legacy_indexing_report(self._residues, self._msa_indices)


@dataclass
class LegacyIndexingReport:
    """How badly ``msa_indices[resid - 1]`` misbehaves on a given structure.

    Retained as a diagnostic so the effect of the fix can be demonstrated per
    structure rather than asserted in the abstract.
    """

    total: int
    correct: int
    wrong: int
    out_of_range: int
    negative_wraparound: int
    examples: List[Tuple[ResidueKey, Optional[int], Optional[int]]]

    @property
    def is_safe(self) -> bool:
        """True when the old arithmetic happened to give the right answer."""
        return self.wrong == 0 and self.out_of_range == 0

    def describe(self) -> str:
        return (
            f"{self.correct}/{self.total} correct, {self.wrong} wrong, "
            f"{self.out_of_range} out of range, "
            f"{self.negative_wraparound} silent negative wrap-around"
        )


def legacy_indexing_report(
    residues: Sequence[StructureResidue],
    msa_indices: Optional[Sequence[int]],
) -> LegacyIndexingReport:
    """Compare correct ordinal lookup against ``msa_indices[resid - 1]``.

    ``negative_wraparound`` counts the dangerous cases: a residue numbered zero
    or below makes the old expression index from the end of the list, which
    Python resolves happily and silently.
    """
    total = len(residues)
    correct = wrong = out_of_range = wraparound = 0
    examples: List[Tuple[ResidueKey, Optional[int], Optional[int]]] = []

    for ordinal, residue in enumerate(residues):
        expected = None if msa_indices is None else msa_indices[ordinal]

        position = residue.resid - 1
        if msa_indices is None:
            legacy: Optional[int] = None
        elif position < 0:
            if -position <= len(msa_indices):
                legacy = msa_indices[position]  # silently wraps to the end
                wraparound += 1
            else:
                legacy = None
                out_of_range += 1
        elif position >= len(msa_indices):
            legacy = None
            out_of_range += 1
        else:
            legacy = msa_indices[position]

        if legacy is not None and legacy == expected:
            correct += 1
        elif legacy is None:
            pass  # already counted as out_of_range
        else:
            wrong += 1
            if len(examples) < 5:
                examples.append((residue.key, expected, legacy))

    return LegacyIndexingReport(
        total=total,
        correct=correct,
        wrong=wrong,
        out_of_range=out_of_range,
        negative_wraparound=wraparound,
        examples=examples,
    )

"""Validate a parsed grades file against ConSurf's own annotated PDB.

Every ConSurf result bundle ships ``*_ATOMS_section_With_ConSurf.pdb``, in which
the conservation **grade** is written into the B-factor column and left blank for
residues ConSurf did not score.  That gives an oracle produced by ConSurf itself,
independent of the grades table, for the exact thing WatCon depends on: the
residue-identity mapping.

Note that the sibling ``msa_aa_variety_percentage.csv`` is **not** a usable
oracle.  Its ``pos`` and ``ConSurf grade`` columns are indexed on the grades
POS, but its amino-acid composition columns are shifted by one row wherever
ConSurf omitted a residue.  Do not cross-check composition against it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

from .model import ConSurfParseResult, ResidueKey

# PDB column slices (0-based, end-exclusive) per the PDB format specification.
_RESNAME = slice(17, 20)
_CHAIN = slice(21, 22)
_RESSEQ = slice(22, 26)
_ICODE = slice(26, 27)
_TEMPFACTOR = slice(60, 66)


@dataclass
class CrossCheckResult:
    """Outcome of comparing grades against an annotated PDB."""

    checked: int = 0
    agreements: int = 0
    mismatches: List[Tuple[ResidueKey, int, Optional[float]]] = field(
        default_factory=list
    )
    missing_from_pdb: List[ResidueKey] = field(default_factory=list)
    blank_in_pdb: List[ResidueKey] = field(default_factory=list)
    unscored_residues: List[ResidueKey] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        """True when every comparable residue agreed."""
        return not self.mismatches and not self.missing_from_pdb

    def summary(self) -> str:
        return (
            f"checked={self.checked} agree={self.agreements} "
            f"mismatch={len(self.mismatches)} "
            f"missing_from_pdb={len(self.missing_from_pdb)} "
            f"blank_in_pdb={len(self.blank_in_pdb)} "
            f"unscored_in_pdb={len(self.unscored_residues)}"
        )


def read_annotated_pdb(
    path: Union[str, Path],
) -> Dict[ResidueKey, Optional[float]]:
    """Map ``(chain, resseq, icode)`` to the B-factor ConSurf wrote.

    The value is ``None`` where the column is blank, which is how ConSurf marks
    a residue it did not score (modified residues, waters, and any part of the
    model outside the analysed region).
    """
    residues: Dict[ResidueKey, Optional[float]] = {}

    for raw in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if not raw.startswith(("ATOM", "HETATM")):
            continue
        if len(raw) < 27:
            continue

        try:
            resseq = int(raw[_RESSEQ].strip())
        except ValueError:
            continue

        icode = raw[_ICODE].strip() or None
        key = (raw[_CHAIN], resseq, icode)

        if key in residues:
            continue

        field_text = raw[_TEMPFACTOR].strip() if len(raw) >= 66 else ""
        try:
            residues[key] = float(field_text) if field_text else None
        except ValueError:
            residues[key] = None

    return residues


def check_against_annotated_pdb(
    result: ConSurfParseResult,
    pdb_path: Union[str, Path],
) -> CrossCheckResult:
    """Compare every mapped record's grade with the PDB B-factor column.

    Returns a report rather than raising, so callers can decide how strict to
    be.  Residues whose B-factor column is blank are counted separately rather
    than treated as mismatches: some bundles (notably predicted-structure runs)
    leave the column empty throughout.
    """
    pdb_residues = read_annotated_pdb(pdb_path)
    report = CrossCheckResult()

    by_key = result.by_pdb_residue()

    for key, record in by_key.items():
        if key not in pdb_residues:
            report.missing_from_pdb.append(key)
            continue

        b_factor = pdb_residues[key]
        if b_factor is None:
            report.blank_in_pdb.append(key)
            continue

        report.checked += 1
        if int(b_factor) == record.grade:
            report.agreements += 1
        else:
            report.mismatches.append((key, record.grade, b_factor))

    # Residues the structure contains but the grades table never scored.  This
    # is the population that motivates ConSurfParseResult.lookup_structural.
    for key, b_factor in pdb_residues.items():
        if key not in by_key and b_factor is None:
            report.unscored_residues.append(key)

    return report

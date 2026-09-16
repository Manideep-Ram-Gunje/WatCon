"""Evolutionary conservation pooled across a protein family.

One ConSurf run per protein, one or more structures per protein, one shared
alignment. This module is the orchestration the family scaffold in
:mod:`WatCon.evolutionary` was built for and never had:

1. every structure is checked against its own protein's ConSurf run
   (:func:`WatCon.evolutionary.enforce_identity`), so a numbering mismatch stops
   the run instead of attaching scores to the wrong residues;
2. every structure is placed on the alignment by sequence, and structures of the
   same protein are cross-checked (:mod:`WatCon.alignment`), which catches rows
   that are themselves misaligned;
3. conservation is pooled onto alignment columns
   (:func:`WatCon.evolutionary.conservation_by_msa_column`).

What the pooled result can and cannot say
-----------------------------------------
Scores are z-normalised within each ConSurf run, so a pooled column reports each
member's score and grade side by side, never an average grade.
``unanimous_conserved`` -- every run independently grading the column 8 or 9 --
is the most defensible family statement because it compares verdicts, not
separately normalised values.

Measured on five PTPs (PTPN1, 6, 7, 12, 22): their homologue sets are almost
disjoint (Jaccard <= 0.01), yet their scores correlate at rho = 0.91-0.96 over
shared columns. The runs are independent samples and they agree, so unanimity
across them is corroboration, not five copies of one alignment.

Nothing here reads waters; site-level family analysis builds on this.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from .alignment import (
    DEFAULT_MIN_IDENTITY,
    AlignmentError,
    ConsensusColumns,
    RowMapping,
    consensus_columns,
    map_structure_to_row,
    read_alignment,
    residue_index_for,
)
from .consurf import parse_consurf
from .evolutionary import (
    ColumnConservation,
    ConservationMap,
    conservation_by_msa_column,
    enforce_identity,
    family_summary,
)
from .residue_index import StructureResidue, residues_from_pdb_file

__all__ = [
    "FamilyStructure",
    "read_members",
    "FamilyProtein",
    "StructureReport",
    "FamilyConservation",
    "build_family_conservation",
]


@dataclass
class FamilyStructure:
    """One structure of a family member, numbered as its ConSurf run is."""

    pdb_id: str
    path: str
    chain: str = "A"
    #: Alignment row for this structure. Default: the one row whose name starts
    #: with ``pdb_id``.
    row_name: Optional[str] = None


@dataclass
class FamilyProtein:
    """One protein: its ConSurf run and its structures."""

    name: str
    consurf_path: str
    structures: List[FamilyStructure]
    #: Which structure's residues represent the protein when pooling. Default:
    #: the first structure.
    reference: Optional[str] = None


@dataclass
class StructureReport:
    """What happened to one structure on its way into the family."""

    protein: str
    pdb_id: str
    n_residues: int
    identity: Optional[float]
    #: (residue key, structure residue, ConSurf residue) -- engineered mutations,
    #: modified residues, or genuine sequence differences.
    identity_mismatches: List[Tuple[tuple, str, str]]
    coverage_fraction: float
    row_name: str
    row_identity: float
    #: Residues the alignment row does not place.
    unplaced: List[Tuple[int, Optional[str]]]


@dataclass
class FamilyConservation:
    """Conservation pooled across a family, with everything needed to audit it."""

    proteins: List[str]
    structures: List[StructureReport]
    consensus: Dict[str, ConsensusColumns]
    columns: Dict[int, ColumnConservation]
    reference_residues: Dict[str, List[StructureResidue]] = field(default_factory=dict)

    def column_for(self, protein: str, resid: int, icode: Optional[str] = None) -> Optional[int]:
        """Alignment column of one protein's residue, or None if it has none."""
        return self.consensus[protein].column(resid, icode)

    def at(self, protein: str, resid: int, icode: Optional[str] = None) -> Optional[ColumnConservation]:
        """Pooled conservation at the column holding one protein's residue."""
        column = self.column_for(protein, resid, icode)
        return None if column is None else self.columns.get(column)

    def conflicts(self) -> Dict[str, List[Tuple[int, Optional[str]]]]:
        """Residues excluded because structures of one protein disagreed on their column."""
        return {name: sorted(c.conflicts, key=lambda k: (k[0], k[1] or ""))
                for name, c in self.consensus.items() if c.conflicts}

    def summary(self) -> dict:
        """Headline counts, including the columns every protein covers."""
        base = family_summary(self.columns)
        everyone = [c for c in self.columns.values() if c.n_members == len(self.proteins)]
        spreads = sorted(c.spread for c in everyone if c.spread is not None)
        base.update({
            "n_proteins": len(self.proteins),
            "n_columns_all_proteins": len(everyone),
            "n_unanimous_all_proteins": sum(1 for c in everyone if c.unanimous_conserved),
            "median_spread_all_proteins": spreads[len(spreads) // 2] if spreads else None,
        })
        return base

    def pairwise_spearman(self) -> Dict[Tuple[str, str], Tuple[int, float]]:
        """``{(protein a, protein b): (shared columns, Spearman rho of scores)}``."""
        from itertools import combinations

        from scipy.stats import spearmanr

        out = {}
        for a, b in combinations(self.proteins, 2):
            shared = [c for c in self.columns.values() if a in c.scores and b in c.scores]
            if len(shared) < 3:
                continue
            rho = spearmanr([c.scores[a] for c in shared], [c.scores[b] for c in shared]).correlation
            out[(a, b)] = (len(shared), float(rho))
        return out


def _row_for(rows: Dict[str, str], structure: FamilyStructure) -> Tuple[str, str]:
    if structure.row_name is not None:
        if structure.row_name not in rows:
            raise AlignmentError("alignment has no row named %r" % structure.row_name)
        return structure.row_name, rows[structure.row_name]
    matches = [name for name in rows if name.upper().startswith(structure.pdb_id.upper())]
    if len(matches) != 1:
        raise AlignmentError(
            "%s matches %d alignment rows (%s); pass row_name to say which"
            % (structure.pdb_id, len(matches), ", ".join(matches) or "none"))
    return matches[0], rows[matches[0]]


def build_family_conservation(
    proteins: Sequence[FamilyProtein],
    alignment_path: str,
    strict: bool = True,
    min_row_identity: float = DEFAULT_MIN_IDENTITY,
) -> FamilyConservation:
    """Check, place and pool a family. See the module docstring.

    Raises
    ------
    ConservationError
        Under ``strict``, if any structure's residues are not the residues its
        protein's ConSurf run describes.
    AlignmentError
        If a structure does not match its alignment row.
    """
    proteins = list(proteins)
    if len({p.name for p in proteins}) != len(proteins):
        raise ValueError("protein names must be unique")

    rows = read_alignment(alignment_path)
    reports: List[StructureReport] = []
    consensus: Dict[str, ConsensusColumns] = {}
    members = []
    reference_residues: Dict[str, List[StructureResidue]] = {}

    for protein in proteins:
        if not protein.structures:
            raise ValueError("%s has no structures" % protein.name)
        conservation = ConservationMap.build(parse_consurf(protein.consurf_path, strict=True))

        mappings: List[RowMapping] = []
        residues_by_id: Dict[str, List[StructureResidue]] = {}
        for structure in protein.structures:
            residues = residues_from_pdb_file(structure.path, chain=structure.chain)
            if not residues:
                raise ValueError("%s chain %s has no residues" % (structure.path, structure.chain))
            residues_by_id[structure.pdb_id] = residues

            coverage = conservation.coverage(residues)
            identity = enforce_identity(coverage, label="%s (%s)" % (structure.pdb_id, protein.name),
                                        strict=strict)

            row_name, row = _row_for(rows, structure)
            mapping = map_structure_to_row(residues, row, label=structure.pdb_id,
                                           row_name=row_name, min_identity=min_row_identity)
            mappings.append(mapping)

            reports.append(StructureReport(
                protein=protein.name, pdb_id=structure.pdb_id, n_residues=len(residues),
                identity=identity, identity_mismatches=list(coverage.mismatches),
                coverage_fraction=coverage.fraction, row_name=row_name,
                row_identity=mapping.identity, unplaced=list(mapping.unmapped),
            ))

        protein_consensus = consensus_columns(protein.name, mappings)
        consensus[protein.name] = protein_consensus

        reference_id = protein.reference or protein.structures[0].pdb_id
        if reference_id not in residues_by_id:
            raise ValueError("%s: reference %r is not one of its structures" % (protein.name, reference_id))
        reference_residues[protein.name] = residues_by_id[reference_id]
        members.append((protein.name, conservation,
                        residue_index_for(residues_by_id[reference_id], protein_consensus)))

    return FamilyConservation(
        proteins=[p.name for p in proteins],
        structures=reports,
        consensus=consensus,
        columns=conservation_by_msa_column(members),
        reference_residues=reference_residues,
    )


# ---------------------------------------------------------------------------
# Describing a family in a file
# ---------------------------------------------------------------------------

def read_members(path: str) -> List[FamilyProtein]:
    """Read a tab-separated members file into :class:`FamilyProtein` objects.

    One line per protein::

        # protein   structures directory   ConSurf grades file   [reference]
        PTPN1       prepared/PTPN1         consurf/1AAX_A.grades.txt   2F71

    Blank lines and lines starting with ``#`` are ignored. The reference is the
    structure whose residues represent the protein when pooling; it defaults to
    the first structure in sorted order.

    Tab-separated rather than a ``name:dir:file`` string because Windows paths
    contain colons, and a format that breaks on ``C:/structures`` is no format
    at all.
    """
    from .residue_index import list_structure_files

    if not os.path.isfile(path):
        raise ValueError("no such members file: %s" % path)

    proteins: List[FamilyProtein] = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = [f.strip() for f in line.split("	") if f.strip()]
            if len(fields) < 3:
                raise ValueError(
                    "%s line %d: expected at least three tab-separated fields "
                    "(protein, structures directory, ConSurf file), found %d: %r"
                    % (path, number, len(fields), line))
            name, directory, consurf = fields[0], fields[1], fields[2]
            reference = fields[3] if len(fields) > 3 else None

            if not os.path.isdir(directory):
                raise ValueError("%s line %d: no such directory: %s" % (path, number, directory))
            if not os.path.isfile(consurf):
                raise ValueError("%s line %d: no such ConSurf file: %s" % (path, number, consurf))

            names, _skipped = list_structure_files(directory)
            if not names:
                raise ValueError("%s line %d: no structures in %s" % (path, number, directory))
            structures = [FamilyStructure(os.path.splitext(f)[0], os.path.join(directory, f))
                          for f in names]
            ids = [s.pdb_id for s in structures]
            if reference is not None and reference not in ids:
                raise ValueError("%s line %d: reference %r is not among %s"
                                 % (path, number, reference, ", ".join(ids)))
            proteins.append(FamilyProtein(name=name, consurf_path=consurf,
                                          structures=structures, reference=reference))

    if not proteins:
        raise ValueError("no members listed in %s" % path)
    return proteins

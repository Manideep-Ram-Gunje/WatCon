"""Map structures of a protein family onto a shared alignment -- by sequence, checked.

Pooling conservation across a family (``evolutionary.conservation_by_msa_column``)
needs, for every residue of every structure, the alignment column it occupies.
WatCon's historical route, ``sequence_processing.generate_msa_alignment``,
returns one column per non-gap letter of a structure's row, in order, and trusts
that the row describes the structure exactly. Real family data breaks that trust
in two ways, both present in the authors' own published PTP alignment:

* **A residue missing from a row.** 5HDE's catalytic nucleophile is deposited as
  CSP231, a HETATM phosphocysteine; its row omits it. Every column after it is
  then assigned to the wrong residue by a positional walk, or -- if the residue
  walk also skipped the HETATM -- the catalytic residue silently vanishes.
* **A row that is itself misaligned.** 3O4U's WPD loop is disordered (native
  W234-P235, then E241). The aligner slid E241 left across the gap into the
  column every other PTP fills with the WPD aspartate, and did the same to
  E124/D125 at another disordered loop. Mapping the structure to that row "by
  sequence" reproduces the slide faithfully; nothing about the row looks wrong.

This module handles both without trusting either:

1. :func:`map_structure_to_row` aligns a structure's own residue sequence to its
   row and refuses below an identity threshold. Residues the row lacks get no
   column, rather than a neighbour's.
2. :func:`consensus_columns` merges the rows of several structures **of one
   protein** by native residue number. Two structures of one protein must put
   residue 241 in the same column; where their rows disagree, the residue is
   reported as a conflict and excluded from pooling. That is what catches the
   3O4U slide -- 1ZC0's ordered loop puts E241 where it belongs. Where only one
   row places a residue, that column is used and its source recorded, which is
   how 5HDE's CSP231 gets a column from 5J8R's row.

Nothing here renumbers a structure, and nothing pools across proteins; that is
:func:`WatCon.evolutionary.conservation_by_msa_column`'s job, fed by
:func:`residue_index_for`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .residue_index import ResidueIndex, StructureResidue

__all__ = [
    "AlignmentError",
    "RowMapping",
    "ConsensusColumns",
    "read_alignment",
    "map_structure_to_row",
    "consensus_columns",
    "residue_index_for",
    "DEFAULT_MIN_IDENTITY",
]

#: Fraction of a structure's residues that must align to its row. Guards against
#: a short local match standing in for the whole chain.
MIN_COVERAGE = 0.80

#: Structure-to-row identity below which a mapping is refused. A structure and
#: its own alignment row describe the same chain and should agree completely;
#: measured on the ten PTP structures against the authors' rows: 1.000 for all.
#: 0.95 leaves room for a handful of genuine differences, not for a wrong row.
DEFAULT_MIN_IDENTITY = 0.95

#: (resid, icode). Chain is deliberately not part of the key: consensus is taken
#: across structures of one protein, whose chain labels are the depositor's
#: choice -- the same reasoning as ``superpose._correspondence``.
PositionKey = Tuple[int, Optional[str]]


class AlignmentError(ValueError):
    """Raised when a structure cannot be placed on an alignment safely."""


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------

def read_alignment(path) -> Dict[str, str]:
    """``{row name: gapped sequence}`` from a PIR or FASTA alignment.

    PIR (MODELLER's format, the one WatCon writes) is recognised by ``>P1;`` or
    ``>F1;`` headers, whose following description line is skipped and whose
    sequence ends with ``*``. Every row must have the same length, or the file is
    not an alignment and is refused.
    """
    # Parsed line by line, with a header being a line that STARTS with '>'.
    # Splitting the text on '>' is wrong for PIR: MODELLER writes description
    # lines such as "<unknown name> - 2F71_PTPN1_closed_aligned", and the '>' in
    # "<unknown name>" would start a bogus record named '-'.
    records: List[Tuple[str, List[str]]] = []
    skip_description = False
    for raw in open(path, encoding="utf-8", errors="replace"):
        line = raw.rstrip("\n")
        if line.startswith(">"):
            header = line[1:].strip()
            if header[:3] in ("P1;", "F1;"):
                name = header[3:].strip()
                skip_description = True          # PIR: next line describes, not sequence
            else:
                name = header.split()[0] if header.split() else ""
                skip_description = False
            if not name:
                raise AlignmentError("a sequence in %s has no name" % path)
            records.append((name, []))
            continue
        if skip_description:
            skip_description = False
            continue
        if records:
            records[-1][1].append(line.strip())

    if not records:
        raise AlignmentError("no sequences in %s" % path)

    rows: Dict[str, str] = {}
    for name, body in records:
        if name in rows:
            raise AlignmentError("row %r appears twice in %s" % (name, path))
        rows[name] = "".join(body).replace("*", "").upper()

    lengths = {len(s) for s in rows.values()}
    if len(lengths) != 1:
        raise AlignmentError(
            "%s is not an alignment: rows have lengths %s" % (path, sorted(lengths)))
    return rows


# ---------------------------------------------------------------------------
# One structure onto its own row
# ---------------------------------------------------------------------------

@dataclass
class RowMapping:
    """Where each residue of one structure sits in the alignment.

    Columns are **1-based**, matching ``generate_msa_alignment``.
    """

    label: str
    row_name: str
    columns: Dict[PositionKey, int]
    identity: float
    n_aligned: int
    #: Which pairwise alignment placed the structure: "global" for a full chain,
    #: "local" for a fragment such as a trimmed active site.
    mode: str = "global"
    #: Residues of the structure with no letter in the row.
    unmapped: List[PositionKey] = field(default_factory=list)
    #: (position, structure letter, row letter) where an aligned pair differs.
    mismatches: List[Tuple[PositionKey, str, str]] = field(default_factory=list)


def _aligner(mode: str = "global"):
    from Bio.Align import PairwiseAligner

    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5
    return aligner


def map_structure_to_row(
    residues: Sequence[StructureResidue],
    row: str,
    label: str = "structure",
    row_name: str = "",
    min_identity: float = DEFAULT_MIN_IDENTITY,
) -> RowMapping:
    """Align one structure's residues to its alignment row, by sequence.

    Raises
    ------
    AlignmentError
        If fewer than ``min_identity`` of the aligned pairs agree -- the row does
        not describe this structure, and any column assigned would be a guess.
    """
    residues = list(residues)
    if not residues:
        raise AlignmentError("%s has no residues to map" % label)

    row_columns = [i + 1 for i, letter in enumerate(row) if letter != "-"]
    row_letters = "".join(row[c - 1] for c in row_columns)
    structure_letters = "".join(r.one_letter for r in residues)

    # Global first, then local, choosing whichever places MORE residues
    # correctly -- the count of identical pairs, not the ratio. Ratio alone
    # rewards a short local match: on one trimmed fixture local alignment
    # matched a single 18-residue segment at 94% and beat a global alignment
    # covering the whole chain. Coverage is then required outright below.
    best = None
    for mode in ("global", "local"):
        alignment = _aligner(mode).align(structure_letters, row_letters)[0]

        columns: Dict[PositionKey, int] = {}
        mismatches = []
        identical = aligned = 0
        for (s_start, s_end), (r_start, r_end) in zip(*alignment.aligned):
            for i, j in zip(range(s_start, s_end), range(r_start, r_end)):
                residue = residues[i]
                key = (residue.resid, residue.icode)
                columns[key] = row_columns[j]
                aligned += 1
                if structure_letters[i] == row_letters[j]:
                    identical += 1
                else:
                    mismatches.append((key, structure_letters[i], row_letters[j]))

        identity = identical / aligned if aligned else 0.0
        if best is None or identical > best[4]:
            best = (identity, columns, mismatches, aligned, identical, mode)
        if identity >= min_identity and aligned >= MIN_COVERAGE * len(residues):
            break

    identity, columns, mismatches, aligned, identical, mode = best
    if aligned < MIN_COVERAGE * len(residues):
        raise AlignmentError(
            "%s could not be placed on alignment row %r: only %d of its %d residues "
            "align to it (need %.0f%%). A structure trimmed to a pocket is several "
            "disconnected segments, which no pairwise alignment can place against a "
            "full-length row -- supply the whole chain."
            % (label, row_name, aligned, len(residues), 100 * MIN_COVERAGE))
    if identity < min_identity:
        # A structure much shorter than its row is usually a trimmed one, and a
        # pocket is several disconnected segments: every residue aligns, but not
        # where it belongs. Say that, rather than only "different sequence".
        hint = ""
        if len(residues) < MIN_COVERAGE * len(row_letters):
            hint = (" %s has %d residues against the row's %d: if it was trimmed to a "
                    "pocket or a domain, supply the whole chain -- disconnected "
                    "segments cannot be placed by pairwise alignment."
                    % (label, len(residues), len(row_letters)))
        raise AlignmentError(
            "%s does not match alignment row %r: %d of %d aligned residues agree "
            "(%.1f%%, need %.0f%%), under both global and local alignment.%s"
            % (label, row_name, identical, aligned, 100 * identity, 100 * min_identity,
               hint or " The row describes a different sequence."))

    unmapped = [(r.resid, r.icode) for r in residues if (r.resid, r.icode) not in columns]
    return RowMapping(label=label, row_name=row_name, columns=columns,
                      identity=identity, n_aligned=aligned, mode=mode,
                      unmapped=unmapped, mismatches=mismatches)


# ---------------------------------------------------------------------------
# Several structures of one protein
# ---------------------------------------------------------------------------

@dataclass
class ConsensusColumns:
    """One column per native residue position of a protein, agreed across rows."""

    protein: str
    columns: Dict[PositionKey, int]
    #: Positions whose rows disagree: {position: {structure label: column}}.
    #: Excluded from ``columns`` -- a residue with two columns has none.
    conflicts: Dict[PositionKey, Dict[str, int]]
    #: Positions placed by only one structure's row: {position: label}.
    single_source: Dict[PositionKey, str]
    labels: List[str]

    def column(self, resid: int, icode: Optional[str] = None) -> Optional[int]:
        """Alignment column for one residue, or None if it has none.

        None is returned rather than guessed: a residue the structures of this
        protein disagree about is excluded, and so is one whose alignment row
        omits it entirely.
        """
        return self.columns.get((resid, icode))


def consensus_columns(protein: str, mappings: Iterable[RowMapping]) -> ConsensusColumns:
    """Merge the row mappings of several structures of **one** protein.

    Positions are native residue numbers, which structures of one protein share
    (verified for every PTP pair: offset 0, identity >= 0.996). Two rows placing
    the same residue in different columns cannot both be right; neither is
    guessed between.
    """
    mappings = list(mappings)
    if not mappings:
        raise AlignmentError("no row mappings given for %s" % protein)

    seen: Dict[PositionKey, Dict[str, int]] = {}
    for mapping in mappings:
        for position, column in mapping.columns.items():
            seen.setdefault(position, {})[mapping.label] = column

    columns, conflicts, single = {}, {}, {}
    for position, by_label in seen.items():
        distinct = set(by_label.values())
        if len(distinct) > 1:
            conflicts[position] = dict(by_label)
            continue
        columns[position] = next(iter(distinct))
        if len(by_label) == 1:
            single[position] = next(iter(by_label))

    # Two different residues of one protein must not share a column either.
    owner: Dict[int, PositionKey] = {}
    for position, column in sorted(columns.items(), key=lambda kv: (kv[0][0], kv[0][1] or "")):
        if column in owner:
            first = owner[column]
            conflicts.setdefault(first, {"column": column})
            conflicts.setdefault(position, {"column": column})
        else:
            owner[column] = position
    for position in conflicts:
        columns.pop(position, None)
        single.pop(position, None)

    return ConsensusColumns(protein=protein, columns=columns, conflicts=conflicts,
                            single_source=single, labels=[m.label for m in mappings])


def residue_index_for(
    residues: Sequence[StructureResidue],
    consensus: ConsensusColumns,
    source: Optional[str] = None,
) -> ResidueIndex:
    """A :class:`ResidueIndex` whose MSA columns come from the consensus.

    Residues without an agreed column -- absent from every row, or in conflict --
    map to ``None``, which ``conservation_by_msa_column`` skips.
    """
    residues = list(residues)
    return ResidueIndex.build(
        residues,
        [consensus.column(r.resid, r.icode) for r in residues],
        source=source,
    )

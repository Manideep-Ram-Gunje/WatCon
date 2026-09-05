"""Data model for parsed ConSurf output.

Design rule: no field carries two meanings.  In particular

* ``residue_variety`` never uses an empty dict to mean "100% conserved" -- a
  fully conserved position has a real entry with ``is_full=True``.
* ``PdbResidueId`` keeps the residue number and the insertion code apart.
* ``Confidence`` grade bounds are named for the *score* bound they belong to,
  because ConSurf grades run opposite to scores.
* "no structural mapping" (``pdb_residue is None``) is distinct from "mapping
  present but unparseable" (``mapping_unverified``) and from "residue exists in
  the structure but ConSurf never scored it"
  (see :meth:`ConSurfParseResult.lookup_structural`).
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from enum import Enum
from typing import Dict, List, Optional, Tuple

from .errors import ConSurfWarning, MalformedLine

# Key used to join a ConSurf record onto a structure.
ResidueKey = Tuple[str, int, Optional[str]]


@dataclass(frozen=True)
class PdbResidueId:
    """Structural identity of a residue, as reported in the ATOM/3LATOM column."""

    chain: str
    seq: int
    icode: Optional[str]
    name3: str
    raw: str

    @property
    def key(self) -> ResidueKey:
        """The (chain, seq, icode) tuple used for joining."""
        return (self.chain, self.seq, self.icode)

    def __str__(self) -> str:
        suffix = self.icode or ""
        return f"{self.name3}:{self.seq}{suffix}:{self.chain}"


@dataclass(frozen=True)
class ResidueFrequency:
    """Frequency of one residue in a ConSurf MSA column."""

    residue: str
    percent: Optional[float]
    less_than_one_percent: bool = False
    is_full: bool = False

    def __post_init__(self) -> None:
        if self.less_than_one_percent and self.percent is not None:
            raise ValueError("'<1%' entries must carry percent=None")
        if self.is_full and self.percent != 100.0:
            raise ValueError("is_full entries must carry percent=100.0")


@dataclass(frozen=True)
class Confidence:
    """Bayesian confidence interval.  Absent entirely under maximum likelihood."""

    score_lower: float
    score_upper: float
    grade_at_lower: int
    grade_at_upper: int

    @property
    def grade_span(self) -> int:
        """Absolute difference between the two bound grades.

        ConSurf marks a position unreliable when this difference exceeds 3.
        Note that ``grade_at_lower >= grade_at_upper`` normally holds, because a
        lower score means a higher (more conserved) grade.
        """
        return abs(self.grade_at_lower - self.grade_at_upper)


@dataclass(frozen=True)
class MsaSupport:
    """Non-gapped sequence support for one alignment column."""

    present: int
    total: int

    @property
    def fraction(self) -> float:
        if self.total == 0:
            return 0.0
        return self.present / self.total


@dataclass(frozen=True)
class GradeLayer:
    """One row of the per-run score-to-grade threshold table in the preamble."""

    grade: int
    score_from: float
    score_to: float


class Dialect(str, Enum):
    WEBSERVER = "webserver"
    CONSURFDB = "consurfdb"


class Method(str, Enum):
    BAYESIAN = "bayesian"
    MAXIMUM_LIKELIHOOD = "maximum_likelihood"


class Alphabet(str, Enum):
    PROTEIN = "protein"
    NUCLEOTIDE = "nucleotide"


class LineEnding(str, Enum):
    LF = "lf"
    CRLF = "crlf"
    MIXED = "mixed"
    UNKNOWN = "unknown"


@dataclass(frozen=True)
class ConSurfProvenance:
    """How this file was produced, inferred from its own content."""

    source: str
    dialect: Dialect
    method: Method
    alphabet: Alphabet
    msa_total: Optional[int]
    line_ending: LineEnding


@dataclass
class ConSurfRecord:
    """Everything ConSurf reports for one scored position."""

    # Sequence level.  NOTE: `position` is ConSurf's SEQRES-derived POS.  It is
    # NOT guaranteed to be a dense index into the query sequence -- ConSurf may
    # omit non-standard residues entirely.
    position: int
    sequence_residue: str

    # Structural mapping.  None when ConSurf reported '-'.
    pdb_residue: Optional[PdbResidueId]

    # Conservation
    score: float
    grade: int
    low_confidence: bool

    # Reliability (None under maximum likelihood)
    confidence: Optional[Confidence]

    # Optional ConSurf annotations
    buried_exposed: Optional[str]
    functional_structural: Optional[str]

    # MSA support
    msa: MsaSupport

    # Composition
    residue_variety: Dict[str, ResidueFrequency] = field(default_factory=dict)
    residue_variety_raw: str = ""

    # Traceability
    raw_line: str = ""
    line_number: int = 0
    warnings: List[ConSurfWarning] = field(default_factory=list)

    # True when an ATOM field was present but could not be parsed.  Such records
    # are deliberately excluded from by_pdb_residue().
    mapping_unverified: bool = False

    @property
    def is_structurally_mapped(self) -> bool:
        return self.pdb_residue is not None

    @property
    def msa_fraction(self) -> float:
        return self.msa.fraction

    @property
    def is_fully_conserved(self) -> bool:
        """True when a single residue accounts for 100% of the MSA column.

        Independent of notation: ConSurf may write either a bare letter or an
        explicit ``100%``, and both mean the same thing.  ``ResidueFrequency.is_full``
        records which notation was used; this property records the fact.
        """
        return any(f.percent == 100.0 for f in self.residue_variety.values())

    @property
    def warning_codes(self) -> List[str]:
        return [w.code for w in self.warnings]

    def to_dict(self) -> dict:
        """JSON-friendly representation."""
        result = asdict(self)
        result["pdb_residue"] = (
            None if self.pdb_residue is None else asdict(self.pdb_residue)
        )
        result["confidence"] = (
            None if self.confidence is None else asdict(self.confidence)
        )
        result["msa"] = asdict(self.msa)
        result["residue_variety"] = {
            aa: {
                "percent": value.percent,
                "less_than_one_percent": value.less_than_one_percent,
                "is_full": value.is_full,
            }
            for aa, value in self.residue_variety.items()
        }
        result["warnings"] = [
            {"code": w.code, "message": w.message, "line_number": w.line_number}
            for w in self.warnings
        ]
        return result


class LookupStatus(str, Enum):
    """Outcome of asking a ConSurf file about a residue seen in a structure."""

    SCORED = "scored"
    NO_SCORE_IN_FILE = "no_score_in_file"
    NOT_IN_FILE = "not_in_file"


@dataclass(frozen=True)
class ConSurfLookup:
    """Three-valued answer to a structural query.

    SCORED
        ConSurf scored this residue; ``record`` is set.
    NO_SCORE_IN_FILE
        The residue falls inside the range ConSurf covered for this chain but
        has no record -- for example a modified residue ConSurf omitted.
    NOT_IN_FILE
        Outside the covered range entirely.
    """

    status: LookupStatus
    record: Optional[ConSurfRecord] = None

    def __bool__(self) -> bool:
        return self.status is LookupStatus.SCORED


@dataclass
class ConSurfParseResult:
    """Complete result of parsing one ConSurf grades file."""

    records: List[ConSurfRecord]
    provenance: ConSurfProvenance
    grade_layers: List[GradeLayer] = field(default_factory=list)
    position_gaps: List[Tuple[int, int]] = field(default_factory=list)
    warnings: List[ConSurfWarning] = field(default_factory=list)
    malformed_lines: List[MalformedLine] = field(default_factory=list)

    # -- convenience -------------------------------------------------------

    @property
    def source(self) -> str:
        return self.provenance.source

    @property
    def warning_codes(self) -> List[str]:
        return [w.code for w in self.warnings]

    def by_position(self) -> Dict[int, ConSurfRecord]:
        """Records keyed by ConSurf POS (SEQRES-derived, possibly non-dense)."""
        return {record.position: record for record in self.records}

    def by_pdb_residue(self) -> Dict[ResidueKey, ConSurfRecord]:
        """Records keyed by (chain, seq, icode).

        Records whose ATOM field could not be parsed are excluded, so a hit in
        this mapping is always a mapping ConSurf actually stated.
        """
        return {
            record.pdb_residue.key: record
            for record in self.records
            if record.pdb_residue is not None
        }

    def mapped_records(self) -> List[ConSurfRecord]:
        return [r for r in self.records if r.pdb_residue is not None]

    def unmapped_records(self) -> List[ConSurfRecord]:
        """Records ConSurf explicitly reported as having no ATOM record ('-')."""
        return [
            r
            for r in self.records
            if r.pdb_residue is None and not r.mapping_unverified
        ]

    def unverified_records(self) -> List[ConSurfRecord]:
        """Records whose ATOM field was present but unparseable (tolerant mode)."""
        return [r for r in self.records if r.mapping_unverified]

    def chains(self) -> List[str]:
        seen = []
        for record in self.records:
            if record.pdb_residue is not None:
                chain = record.pdb_residue.chain
                if chain not in seen:
                    seen.append(chain)
        return seen

    def covered_range(self, chain: str) -> Optional[Tuple[int, int]]:
        """Lowest and highest residue number ConSurf mapped on ``chain``."""
        nums = [
            r.pdb_residue.seq
            for r in self.records
            if r.pdb_residue is not None and r.pdb_residue.chain == chain
        ]
        if not nums:
            return None
        return (min(nums), max(nums))

    def lookup_structural(
        self,
        chain: str,
        seq: int,
        icode: Optional[str] = None,
    ) -> ConSurfLookup:
        """Ask what this file knows about a residue observed in a structure.

        This is the query direction WatCon needs: it starts from the structure,
        not from the ConSurf sequence.  It never guesses -- an unmatched residue
        inside the covered range is reported as NO_SCORE_IN_FILE rather than
        being silently treated as absent.
        """
        record = self.by_pdb_residue().get((chain, seq, icode))
        if record is not None:
            return ConSurfLookup(LookupStatus.SCORED, record)

        bounds = self.covered_range(chain)
        if bounds is not None and bounds[0] <= seq <= bounds[1]:
            return ConSurfLookup(LookupStatus.NO_SCORE_IN_FILE)

        return ConSurfLookup(LookupStatus.NOT_IN_FILE)

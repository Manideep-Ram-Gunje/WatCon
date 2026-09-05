"""Typed errors and warning codes for the ConSurf parser.

Warnings are structured rather than bare strings so that tests and callers can
assert on a stable ``code`` instead of matching English prose.
"""

from __future__ import annotations

from dataclasses import dataclass


class ConSurfParseError(ValueError):
    """Raised when parsing encounters content that cannot be trusted.

    Used for conditions that would otherwise produce a silently wrong result --
    notably an unverifiable ATOM/structural mapping, duplicate identities, and
    values outside ConSurf's documented ranges.
    """


class WarningCode:
    """Stable identifiers for non-fatal parser observations."""

    # Reliability, as defined by ConSurf itself
    LOW_CONFIDENCE = "low_confidence"
    SHALLOW_MSA = "shallow_msa"
    WIDE_CONFIDENCE_INTERVAL = "wide_confidence_interval"

    # Symbol handling
    UNUSUAL_SEQ_SYMBOL = "unusual_seq_symbol"

    # Residue variety
    UNPARSED_VARIETY_ITEM = "unparsed_variety_item"
    VARIETY_PERCENT_SUM_LOW = "variety_percent_sum_low"
    DUPLICATE_VARIETY_RESIDUE = "duplicate_variety_residue"

    # Structural mapping
    NONPOSITIVE_RESIDUE_NUMBER = "nonpositive_residue_number"
    PDB_NUMBER_DISCONTINUITY = "pdb_number_discontinuity"
    MAPPING_UNVERIFIED = "mapping_unverified"

    # Sequence position
    NONPOSITIVE_POSITION = "nonpositive_position"
    NON_ASCENDING_POSITIONS = "non_ascending_positions"

    # File-level
    NO_CONFIDENCE_INTERVAL = "no_confidence_interval"
    MISSING_GRADE_LAYERS = "missing_grade_layers"
    MIXED_LINE_ENDINGS = "mixed_line_endings"
    MALFORMED_ROW = "malformed_row"


@dataclass(frozen=True)
class ConSurfWarning:
    """A non-fatal observation attached to a record or to the whole result."""

    code: str
    message: str
    line_number: int = 0

    def __str__(self) -> str:
        return self.message


@dataclass(frozen=True)
class MalformedLine:
    """A data row that could not be parsed, retained for inspection."""

    line_number: int
    raw: str
    reason: str

    def __str__(self) -> str:
        return f"line {self.line_number}: {self.raw} :: {self.reason}"

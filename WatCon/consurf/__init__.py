"""ConSurf evolutionary-conservation parsing for WatCon.

Standard library only -- this package must stay importable without MDAnalysis
or modeller so that it can be tested independently of core WatCon.

Typical use::

    from WatCon.consurf import parse_consurf

    result = parse_consurf("7O7W_A_consurf_grades.txt")
    hit = result.lookup_structural(chain="A", seq=-5)
    if hit:
        print(hit.record.grade, hit.record.score)
"""

from .crosscheck import (
    CrossCheckResult,
    check_against_annotated_pdb,
    read_annotated_pdb,
)
from .errors import (
    ConSurfParseError,
    ConSurfWarning,
    MalformedLine,
    WarningCode,
)
from .model import (
    Alphabet,
    Confidence,
    ConSurfLookup,
    ConSurfParseResult,
    ConSurfProvenance,
    ConSurfRecord,
    Dialect,
    GradeLayer,
    LineEnding,
    LookupStatus,
    Method,
    MsaSupport,
    PdbResidueId,
    ResidueFrequency,
    ResidueKey,
)
from .parser import parse_consurf

__all__ = [
    # parsing
    "parse_consurf",
    # errors and warnings
    "ConSurfParseError",
    "ConSurfWarning",
    "MalformedLine",
    "WarningCode",
    # model
    "ConSurfParseResult",
    "ConSurfRecord",
    "ConSurfProvenance",
    "ConSurfLookup",
    "LookupStatus",
    "PdbResidueId",
    "ResidueFrequency",
    "ResidueKey",
    "Confidence",
    "MsaSupport",
    "GradeLayer",
    "Dialect",
    "Method",
    "Alphabet",
    "LineEnding",
    # cross-validation
    "check_against_annotated_pdb",
    "read_annotated_pdb",
    "CrossCheckResult",
]

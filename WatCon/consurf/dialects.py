"""Row grammars and symbol tables for the ConSurf grades format.

Two dialects are supported.

**webserver** -- what the ConSurf web server emits today.  Ten tab-separated
fields; the residue identifier is ``THR:-2:A``; the confidence interval and its
bound grades share one field; B/E and F/S columns are present.

**consurfdb** -- the ConSurf-DB / ``gradesPE`` layout.  Nine fields; the
identifier column is called ``3LATOM`` and is written ``MET1:A``; the confidence
interval and its colours are separate fields; there are no B/E or F/S columns.

Both are read by the same whitespace-tolerant regexes, because the difference
between "one merged field" and "two adjacent fields" vanishes once the row is
tokenised on whitespace.  The parts that genuinely differ -- the identifier form
and the presence of B/E and F/S -- are handled by :func:`parse_atom_field` and by
the optional group in the row patterns.

The confidence interval is **Bayesian only**.  ConSurf's documentation states
that no confidence interval is given under maximum likelihood, so a row grammar
without those columns is tried as a fallback.
"""

from __future__ import annotations

import re
from typing import Optional, Tuple

from .model import Dialect, PdbResidueId

# --------------------------------------------------------------------------
# Symbol tables
# --------------------------------------------------------------------------

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

# ConSurf also runs on nucleotides; U is replaced by T when the MSA is built.
NUCLEOTIDES = set("ACGTU")

# Ambiguity / non-standard single-letter codes that legitimately appear.
#   X unknown, B Asx, Z Glx, J Xle, U selenocysteine, O pyrrolysine
AMBIGUITY_CODES = set("XBZJUO")

# Alignment and placeholder characters.
SPECIAL_SYMBOLS = set("-*?")

ACCEPTED_SEQ_SYMBOLS = (
    STANDARD_AA | NUCLEOTIDES | AMBIGUITY_CODES | SPECIAL_SYMBOLS
)

# Residues that mark a position as protein rather than nucleotide.
PROTEIN_ONLY_AA = STANDARD_AA - NUCLEOTIDES

# Codes that are accepted but worth reporting.  U is excluded because it is a
# legitimate RNA base as well as selenocysteine, and warning on it would be pure
# noise for nucleotide runs.
NOTEWORTHY_SEQ_SYMBOLS = AMBIGUITY_CODES - NUCLEOTIDES


# --------------------------------------------------------------------------
# Row grammars
# --------------------------------------------------------------------------

_NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)"

# Shared leading columns: POS SEQ ATOM SCORE COLOR
_HEAD = rf"""
    ^\s*
    (?P<pos>\d+)                    \s+
    (?P<seq>[A-Za-z*?-])            \s+
    (?P<atom>\S+)                   \s+
    (?P<score>{_NUMBER})            \s+
    (?P<color>[1-9]\*?)             \s+
"""

# Shared trailing columns: [B/E] [F/S] MSA-DATA RESIDUE-VARIETY
#
# The optional group is what absorbs the B/E and F/S columns.  In the
# consurfdb dialect they simply are not there and the group matches empty.
_TAIL = r"""
    (?P<optional>(?:[befsBEFS]\s+){0,2})
    (?P<msa>\d+\s*/\s*\d+)
    \s*
    (?P<variety>.*?)
    \s*$
"""

# Bayesian: "lower, upper  gradeLower,gradeUpper"
_CONFIDENCE = rf"""
    (?P<ci_lower>{_NUMBER})   \s*,\s*
    (?P<ci_upper>{_NUMBER})   \s+
    (?P<ci_grade_lower>[1-9]) \s*,\s*
    (?P<ci_grade_upper>[1-9]) \s+
"""

ROW_WITH_CONFIDENCE = re.compile(
    _HEAD + _CONFIDENCE + _TAIL,
    re.VERBOSE,
)

# Maximum-likelihood runs omit the confidence columns entirely.
ROW_WITHOUT_CONFIDENCE = re.compile(
    _HEAD + _TAIL,
    re.VERBOSE,
)

# A line that could plausibly be a data row (used to skip preamble prose).
DATA_ROW_HINT = re.compile(r"^\s*\d+\s+\S")

# Preamble: "from -1.269 to -0.987 the grade is 9"
GRADE_LAYER_RE = re.compile(
    rf"^\s*from\s+(?P<from>{_NUMBER})\s+to\s+(?P<to>{_NUMBER})"
    r"\s+the\s+grade\s+is\s+(?P<grade>[1-9])\s*$",
    re.IGNORECASE,
)


# --------------------------------------------------------------------------
# ATOM / 3LATOM identifier
# --------------------------------------------------------------------------

# webserver form:  NAME:NUMBER[ICODE]:CHAIN   e.g. THR:-2:A, HIS:100A:A
_ATOM_COLON_RE = re.compile(
    r"^(?P<name>[A-Za-z][A-Za-z0-9]{0,2})"
    r":(?P<number>-?\d+)(?P<icode>[A-Za-z])?"
    r":(?P<chain>[^:]+)$"
)

# consurfdb form: NAMENUMBER[ICODE]:CHAIN     e.g. MET1:A, LEU-5:A
# The name must be three non-digit characters so the boundary is unambiguous.
_ATOM_CONCAT_RE = re.compile(
    r"^(?P<name>[A-Za-z_]{3})"
    r"(?P<number>-?\d+)(?P<icode>[A-Za-z])?"
    r":(?P<chain>[^:]+)$"
)

# Values that mean "this residue has no ATOM record".
_UNMAPPED_TOKENS = {"-", "--", "---"}

# consurfdb writes missing density as ___<pos>:<chain>
_DB_MISSING_RE = re.compile(r"^_{3}-?\d+(?:[A-Za-z])?:[^:]+$")


def is_unmapped_atom(raw: str) -> bool:
    """True when the ATOM field states that no structural mapping exists."""
    token = raw.strip()
    if token in _UNMAPPED_TOKENS:
        return True
    return bool(_DB_MISSING_RE.match(token))


def parse_atom_field(raw: str) -> Optional[PdbResidueId]:
    """Parse an ATOM/3LATOM identifier.

    Returns ``None`` when the field cannot be interpreted.  The caller decides
    whether that is fatal -- this function never guesses a mapping.

    Handles negative and zero residue numbers (expression tags) and insertion
    codes, and treats the chain identifier as an opaque string.
    """
    token = raw.strip()

    for pattern in (_ATOM_COLON_RE, _ATOM_CONCAT_RE):
        match = pattern.match(token)
        if match is None:
            continue
        return PdbResidueId(
            chain=match.group("chain"),
            seq=int(match.group("number")),
            icode=match.group("icode"),
            name3=match.group("name").upper(),
            raw=token,
        )

    return None


# --------------------------------------------------------------------------
# Header handling
# --------------------------------------------------------------------------

def is_header_line(stripped: str) -> bool:
    """Recognise the grades table header row."""
    if not stripped.upper().startswith("POS"):
        return False
    upper = stripped.upper()
    if "SEQ" not in upper:
        return False
    return ("MSA DATA" in upper) or ("3LATOM" in upper) or ("ATOM" in upper)


def detect_dialect(header: Optional[str], sample_atom: Optional[str]) -> Dialect:
    """Infer the dialect from the header row, falling back to identifier shape."""
    if header is not None and "3LATOM" in header.upper():
        return Dialect.CONSURFDB

    if sample_atom:
        token = sample_atom.strip()
        if _ATOM_COLON_RE.match(token):
            return Dialect.WEBSERVER
        if _ATOM_CONCAT_RE.match(token) or _DB_MISSING_RE.match(token):
            return Dialect.CONSURFDB

    return Dialect.WEBSERVER


# --------------------------------------------------------------------------
# Residue variety
# --------------------------------------------------------------------------

# "A 20%", "G <1%", or a bare "M" meaning 100%.
VARIETY_ITEM_RE = re.compile(
    r"^(?P<aa>[A-Za-z*?])"
    r"(?:\s+(?P<pct><1|\d+(?:\.\d+)?)%)?$"
)


def split_optional_flags(raw: str) -> Tuple[Optional[str], Optional[str]]:
    """Split the merged B/E and F/S capture into its two components."""
    flags = re.sub(r"\s+", "", raw or "").lower()
    buried_exposed = next((c for c in flags if c in "be"), None)
    functional_structural = next((c for c in flags if c in "fs"), None)
    return buried_exposed, functional_structural

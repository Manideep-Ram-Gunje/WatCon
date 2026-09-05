"""Robust parser for ConSurf ``*_consurf_grades.txt`` output.

This module is intentionally independent of WatCon and of every third-party
package: it imports only the standard library.  That is deliberate.  The rest of
WatCon needs MDAnalysis and modeller, so keeping this package dependency-free is
what makes it testable on its own.

What the parser guarantees
--------------------------

* It never invents a structural mapping.  An ATOM field that is present but
  unparseable is an error in strict mode, and in tolerant mode the record is
  flagged ``mapping_unverified`` and excluded from :meth:`by_pdb_residue`.
* It records how the file was produced (dialect, Bayesian vs maximum
  likelihood, protein vs nucleotide) rather than assuming.
* It reports the per-run score-to-grade thresholds, without which ``score``
  values from different ConSurf runs are not comparable.
* It reports discontinuities in the PDB numbering, which are the only in-file
  signal that ConSurf omitted a residue (for example a modified residue such as
  a chromophore) and that ``position`` is therefore not a dense sequence index.

Sign convention
---------------

``score`` is a normalised rate: **more negative means more conserved**.
``grade`` runs the other way: **9 is most conserved, 1 most variable**.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Optional, TextIO, Tuple, Union

from .dialects import (
    ACCEPTED_SEQ_SYMBOLS,
    DATA_ROW_HINT,
    GRADE_LAYER_RE,
    NOTEWORTHY_SEQ_SYMBOLS,
    NUCLEOTIDES,
    PROTEIN_ONLY_AA,
    ROW_WITH_CONFIDENCE,
    ROW_WITHOUT_CONFIDENCE,
    SPECIAL_SYMBOLS,
    VARIETY_ITEM_RE,
    detect_dialect,
    is_header_line,
    is_unmapped_atom,
    parse_atom_field,
    split_optional_flags,
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
    ConSurfParseResult,
    ConSurfProvenance,
    ConSurfRecord,
    Dialect,
    GradeLayer,
    LineEnding,
    Method,
    MsaSupport,
    ResidueFrequency,
)

__all__ = ["parse_consurf"]

# ConSurf's own reliability thresholds, quoted from the file preamble:
# "If the difference between the colors of the CONFIDENCE INTERVAL COLORS is
#  more than 3 or the msa number ... is less than 6, there is insufficient data"
_SHALLOW_MSA_THRESHOLD = 6
_WIDE_INTERVAL_THRESHOLD = 3

# Percentages are truncated by ConSurf and legitimately fall short of 100.
# Only a large shortfall is worth reporting.
_VARIETY_SUM_WARN_BELOW = 90.0


# ---------------------------------------------------------------------------
# Residue variety
# ---------------------------------------------------------------------------

def _parse_variety(
    raw: str,
    line_number: int,
    warnings: List[ConSurfWarning],
) -> Dict[str, ResidueFrequency]:
    """Parse the RESIDUE VARIETY column.

    A bare residue letter with no percentage means 100%.  ConSurf writes fully
    conserved positions that way, and dropping them -- as an earlier version of
    this parser did -- silently discarded the composition of exactly the most
    conserved positions.
    """
    raw = raw.strip()
    if not raw:
        return {}

    result: Dict[str, ResidueFrequency] = {}

    for item in raw.split(","):
        item = item.strip()
        if not item:
            continue

        match = VARIETY_ITEM_RE.fullmatch(item)
        if not match:
            warnings.append(
                ConSurfWarning(
                    WarningCode.UNPARSED_VARIETY_ITEM,
                    f"Unparsed residue-variety item: {item!r}",
                    line_number,
                )
            )
            continue

        aa = match.group("aa").upper()
        pct = match.group("pct")

        if aa in result:
            # ConSurf should never list a residue twice.  Overwriting would
            # discard the first value with no trace, which is exactly the kind
            # of silent loss this parser exists to avoid.
            warnings.append(
                ConSurfWarning(
                    WarningCode.DUPLICATE_VARIETY_RESIDUE,
                    f"Residue {aa!r} listed more than once in RESIDUE VARIETY; "
                    f"keeping the first value ({result[aa].percent})",
                    line_number,
                )
            )
            continue

        if pct is None:
            # Bare letter: ConSurf omits the percentage at 100%.
            result[aa] = ResidueFrequency(
                residue=aa, percent=100.0, less_than_one_percent=False, is_full=True
            )
        elif pct == "<1":
            result[aa] = ResidueFrequency(
                residue=aa, percent=None, less_than_one_percent=True
            )
        else:
            result[aa] = ResidueFrequency(
                residue=aa, percent=float(pct), less_than_one_percent=False
            )

    return result


# ---------------------------------------------------------------------------
# One data row
# ---------------------------------------------------------------------------

def _match_row(line: str) -> Tuple[Optional[re.Match], bool]:
    """Match a data row, returning ``(match, has_confidence)``.

    The Bayesian grammar is tried first.  A file produced with maximum
    likelihood has no confidence columns at all, so the shorter grammar is the
    fallback rather than an error.
    """
    match = ROW_WITH_CONFIDENCE.match(line)
    if match is not None:
        return match, True

    match = ROW_WITHOUT_CONFIDENCE.match(line)
    if match is not None:
        return match, False

    return None, False


def _parse_record(line: str, line_number: int, strict: bool) -> ConSurfRecord:
    match, has_confidence = _match_row(line)
    if match is None:
        raise ConSurfParseError(f"Line {line_number}: malformed data row")

    data = match.groupdict()
    warnings: List[ConSurfWarning] = []

    position = int(data["pos"])
    sequence_residue = data["seq"].upper()

    if (
        sequence_residue not in ACCEPTED_SEQ_SYMBOLS
        or sequence_residue in NOTEWORTHY_SEQ_SYMBOLS
    ):
        # Tolerated, never fatal: ConSurf legitimately emits ambiguity codes,
        # and one unfamiliar symbol must not cost us the whole file.
        warnings.append(
            ConSurfWarning(
                WarningCode.UNUSUAL_SEQ_SYMBOL,
                f"Line {line_number}: unusual SEQ value {sequence_residue!r}",
                line_number,
            )
        )

    # -- structural mapping ------------------------------------------------
    atom_raw = data["atom"]
    pdb_residue = None
    mapping_unverified = False

    if not is_unmapped_atom(atom_raw):
        pdb_residue = parse_atom_field(atom_raw)

        if pdb_residue is None:
            message = (
                f"Line {line_number}: unrecognised ATOM identifier {atom_raw!r}; "
                "refusing to guess a structural mapping"
            )
            if strict:
                raise ConSurfParseError(message)
            mapping_unverified = True
            warnings.append(
                ConSurfWarning(WarningCode.MAPPING_UNVERIFIED, message, line_number)
            )

        elif pdb_residue.seq <= 0:
            # Real and legitimate (expression tags), but a downstream consumer
            # that indexes a list by residue number will wrap around silently.
            warnings.append(
                ConSurfWarning(
                    WarningCode.NONPOSITIVE_RESIDUE_NUMBER,
                    f"Line {line_number}: non-positive PDB residue number "
                    f"{pdb_residue.seq} ({atom_raw!r})",
                    line_number,
                )
            )

    # -- conservation ------------------------------------------------------
    score = float(data["score"])
    grade_token = data["color"]
    grade = int(grade_token.rstrip("*"))
    low_confidence = grade_token.endswith("*")

    if not 1 <= grade <= 9:
        raise ConSurfParseError(f"Line {line_number}: COLOR grade outside 1-9")

    if low_confidence:
        warnings.append(
            ConSurfWarning(
                WarningCode.LOW_CONFIDENCE,
                "ConSurf marks this residue as low-confidence (*)",
                line_number,
            )
        )

    # -- confidence interval (Bayesian only) -------------------------------
    confidence = None
    if has_confidence:
        confidence = Confidence(
            score_lower=float(data["ci_lower"]),
            score_upper=float(data["ci_upper"]),
            grade_at_lower=int(data["ci_grade_lower"]),
            grade_at_upper=int(data["ci_grade_upper"]),
        )
        if confidence.score_lower > confidence.score_upper:
            raise ConSurfParseError(
                f"Line {line_number}: confidence interval bounds are inverted"
            )
        if confidence.grade_span > _WIDE_INTERVAL_THRESHOLD:
            warnings.append(
                ConSurfWarning(
                    WarningCode.WIDE_CONFIDENCE_INTERVAL,
                    "Confidence-interval grade span exceeds 3 grades",
                    line_number,
                )
            )

    # -- optional annotations ----------------------------------------------
    buried_exposed, functional_structural = split_optional_flags(data.get("optional"))

    # -- MSA support -------------------------------------------------------
    msa_parts = data["msa"].replace(" ", "").split("/")
    msa = MsaSupport(present=int(msa_parts[0]), total=int(msa_parts[1]))

    if msa.total <= 0 or msa.present < 0 or msa.present > msa.total:
        raise ConSurfParseError(
            f"Line {line_number}: invalid MSA DATA {data['msa']!r}"
        )

    if msa.present < _SHALLOW_MSA_THRESHOLD:
        warnings.append(
            ConSurfWarning(
                WarningCode.SHALLOW_MSA,
                "Fewer than 6 non-gapped sequences support this position",
                line_number,
            )
        )

    # -- residue variety ---------------------------------------------------
    residue_variety_raw = data["variety"].strip()
    residue_variety = _parse_variety(residue_variety_raw, line_number, warnings)

    known = sum(
        f.percent for f in residue_variety.values() if f.percent is not None
    )
    implied = sum(1 for f in residue_variety.values() if f.less_than_one_percent)
    if residue_variety and (known + implied) < _VARIETY_SUM_WARN_BELOW:
        warnings.append(
            ConSurfWarning(
                WarningCode.VARIETY_PERCENT_SUM_LOW,
                f"Residue-variety percentages sum to {known:.0f}%, "
                "well below 100%",
                line_number,
            )
        )

    return ConSurfRecord(
        position=position,
        sequence_residue=sequence_residue,
        pdb_residue=pdb_residue,
        score=score,
        grade=grade,
        low_confidence=low_confidence,
        confidence=confidence,
        buried_exposed=buried_exposed,
        functional_structural=functional_structural,
        msa=msa,
        residue_variety=residue_variety,
        residue_variety_raw=residue_variety_raw,
        raw_line=line.rstrip("\r\n"),
        line_number=line_number,
        warnings=warnings,
        mapping_unverified=mapping_unverified,
    )


# ---------------------------------------------------------------------------
# Input handling
# ---------------------------------------------------------------------------

def _read_source(source: Union[str, Path, TextIO]) -> Tuple[str, str, LineEnding]:
    """Return ``(text, name, line_ending)`` without normalising newlines."""
    if hasattr(source, "read"):
        text = source.read()
        name = getattr(source, "name", "<stream>")
    else:
        path = Path(source)
        if not path.exists():
            raise ConSurfParseError(f"ConSurf file does not exist: {path}")
        if not path.is_file():
            raise ConSurfParseError(f"ConSurf path is not a file: {path}")
        text = path.read_bytes().decode("utf-8", errors="replace")
        name = str(path)

    crlf = text.count("\r\n")
    bare_lf = text.count("\n") - crlf
    if crlf and bare_lf:
        ending = LineEnding.MIXED
    elif crlf:
        ending = LineEnding.CRLF
    elif bare_lf:
        ending = LineEnding.LF
    else:
        ending = LineEnding.UNKNOWN

    return text, name, ending


#: A/C/G/T are valid amino-acid codes too, so a handful of rows is never
#: enough to call a file "nucleotide".  Require a real run's worth of evidence.
_MIN_ROWS_FOR_NUCLEOTIDE_CALL = 10
_MIN_DISTINCT_BASES = 3


def _infer_alphabet(records: List[ConSurfRecord]) -> Alphabet:
    """Protein unless the evidence for a nucleotide run is unambiguous.

    Protein is the safe default: every nucleotide symbol is also a valid
    amino-acid code, so only a sizeable file using nothing but bases counts.
    """
    symbols = {r.sequence_residue for r in records}

    if symbols & PROTEIN_ONLY_AA:
        return Alphabet.PROTEIN
    if len(records) < _MIN_ROWS_FOR_NUCLEOTIDE_CALL:
        return Alphabet.PROTEIN
    if len(symbols & NUCLEOTIDES) < _MIN_DISTINCT_BASES:
        return Alphabet.PROTEIN
    if symbols <= (NUCLEOTIDES | SPECIAL_SYMBOLS):
        return Alphabet.NUCLEOTIDE
    return Alphabet.PROTEIN


def _find_position_gaps(records: List[ConSurfRecord]) -> List[Tuple[int, int]]:
    """Discontinuities in PDB numbering among consecutive mapped records.

    A gap is the fingerprint of a residue ConSurf omitted from the table -- for
    example the chromophore of a fluorescent protein.  It is the only in-file
    evidence that POS is not a dense sequence index.
    """
    gaps: List[Tuple[int, int]] = []
    previous: Dict[str, int] = {}

    for record in records:
        if record.pdb_residue is None:
            continue
        chain = record.pdb_residue.chain
        current = record.pdb_residue.seq
        if chain in previous and current != previous[chain] + 1:
            gaps.append((previous[chain], current))
        previous[chain] = current

    return gaps


# ---------------------------------------------------------------------------
# Public parser
# ---------------------------------------------------------------------------

def parse_consurf(
    source: Union[str, Path, TextIO],
    *,
    strict: bool = True,
) -> ConSurfParseResult:
    """Parse a ConSurf grades file.

    Parameters
    ----------
    source:
        File path or an already-open text stream.
    strict:
        ``True`` -- anything that would yield untrustworthy data raises
        :class:`~WatCon.consurf.errors.ConSurfParseError`.

        ``False`` -- malformed rows are skipped and recorded in
        ``result.malformed_lines``; file-level problems become warnings.

    Returns
    -------
    ConSurfParseResult
    """
    text, name, line_ending = _read_source(source)

    records: List[ConSurfRecord] = []
    warnings: List[ConSurfWarning] = []
    malformed_lines: List[MalformedLine] = []
    grade_layers: List[GradeLayer] = []

    header_line: Optional[str] = None
    first_atom: Optional[str] = None
    confidence_seen = False
    confidence_absent_seen = False

    for line_number, line in enumerate(text.splitlines(), start=1):
        stripped = line.strip()
        if not stripped:
            continue

        # Per-run score-to-grade thresholds, needed to compare scores between
        # runs.  Must be tested before the data-row check: these lines start
        # with "from", not a digit, so ordering is safe either way.
        layer = GRADE_LAYER_RE.match(stripped)
        if layer is not None:
            grade_layers.append(
                GradeLayer(
                    grade=int(layer.group("grade")),
                    score_from=float(layer.group("from")),
                    score_to=float(layer.group("to")),
                )
            )
            continue

        if header_line is None and is_header_line(stripped):
            header_line = stripped
            continue

        # Everything that is not plausibly a data row is preamble prose.
        if not DATA_ROW_HINT.match(line):
            continue

        try:
            record = _parse_record(line, line_number, strict=strict)
        except ConSurfParseError as exc:
            if strict:
                raise
            malformed_lines.append(
                MalformedLine(line_number, line.rstrip(), str(exc))
            )
            continue

        if record.confidence is None:
            confidence_absent_seen = True
        else:
            confidence_seen = True

        if first_atom is None and record.pdb_residue is not None:
            first_atom = record.pdb_residue.raw

        records.append(record)

    # -- file-level validation --------------------------------------------

    def _fail(code: str, message: str) -> None:
        if strict:
            raise ConSurfParseError(message)
        warnings.append(ConSurfWarning(code, message))

    if header_line is None:
        _fail(WarningCode.MALFORMED_ROW, "ConSurf grades header was not found")

    if not records:
        _fail(WarningCode.MALFORMED_ROW, "No ConSurf residue records were parsed")

    positions = [r.position for r in records]
    if len(positions) != len(set(positions)):
        _fail(WarningCode.MALFORMED_ROW, "Duplicate residue positions detected")

    # ConSurf POS is 1-based and written in ascending order.  Neither is worth
    # failing over, but a violation means the file is not what we think it is,
    # and position_gaps assumes file order is ascending.
    nonpositive = [p for p in positions if p < 1]
    if nonpositive:
        warnings.append(
            ConSurfWarning(
                WarningCode.NONPOSITIVE_POSITION,
                f"POS values below 1 detected ({nonpositive[:5]}); ConSurf POS "
                "is 1-based",
            )
        )

    if any(b <= a for a, b in zip(positions, positions[1:])):
        warnings.append(
            ConSurfWarning(
                WarningCode.NON_ASCENDING_POSITIONS,
                "POS values are not strictly ascending; the file may be "
                "concatenated or reordered, and reported position_gaps are "
                "unreliable",
            )
        )

    keys = [r.pdb_residue.key for r in records if r.pdb_residue is not None]
    if len(keys) != len(set(keys)):
        duplicates = sorted({k for k in keys if keys.count(k) > 1})
        _fail(
            WarningCode.MALFORMED_ROW,
            f"Duplicate PDB residue identities detected: {duplicates}",
        )

    # -- provenance --------------------------------------------------------

    method = Method.BAYESIAN if confidence_seen else Method.MAXIMUM_LIKELIHOOD
    if not confidence_seen and confidence_absent_seen:
        warnings.append(
            ConSurfWarning(
                WarningCode.NO_CONFIDENCE_INTERVAL,
                "No confidence interval present; treating this as a "
                "maximum-likelihood run",
            )
        )

    msa_totals = {r.msa.total for r in records}
    msa_total = msa_totals.pop() if len(msa_totals) == 1 else None

    provenance = ConSurfProvenance(
        source=name,
        dialect=detect_dialect(header_line, first_atom),
        method=method,
        alphabet=_infer_alphabet(records),
        msa_total=msa_total,
        line_ending=line_ending,
    )

    if not grade_layers:
        warnings.append(
            ConSurfWarning(
                WarningCode.MISSING_GRADE_LAYERS,
                "No score-to-grade threshold table found; scores from this file "
                "cannot be compared with scores from other ConSurf runs",
            )
        )

    if line_ending is LineEnding.MIXED:
        warnings.append(
            ConSurfWarning(
                WarningCode.MIXED_LINE_ENDINGS, "File mixes LF and CRLF line endings"
            )
        )

    position_gaps = _find_position_gaps(records)
    for previous, current in position_gaps:
        warnings.append(
            ConSurfWarning(
                WarningCode.PDB_NUMBER_DISCONTINUITY,
                f"PDB numbering jumps {previous} -> {current}; ConSurf omitted "
                "one or more residues, so POS is not a dense sequence index",
            )
        )

    return ConSurfParseResult(
        records=records,
        provenance=provenance,
        grade_layers=sorted(grade_layers, key=lambda layer: layer.grade),
        position_gaps=position_gaps,
        warnings=warnings,
        malformed_lines=malformed_lines,
    )

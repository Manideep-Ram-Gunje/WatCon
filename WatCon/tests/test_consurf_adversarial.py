"""SYNTHETIC adversarial tests: malformed input and error-handling contract.

These are constructed inputs, not real ConSurf output.  They pin down the
promise that matters most for the WatCon integration:

    the parser never silently guesses a residue mapping.

Anything that would produce untrustworthy structural identity is a hard error in
strict mode, and in tolerant mode is quarantined -- never quietly folded into the
join table.
"""

from __future__ import annotations

import io

import pytest

from WatCon.consurf import ConSurfParseError, WarningCode, parse_consurf

from .test_consurf_dialects import PREAMBLE, WEBSERVER_HEADER, webserver, ws_row


def raw(*rows: str) -> io.StringIO:
    """A file with the header but no preamble grade-layer table."""
    return io.StringIO(WEBSERVER_HEADER + "".join(rows))


# ==========================================================================
# Structural mapping: never guess
# ==========================================================================

BAD_ATOM_FIELDS = [
    "WEIRD_THING",
    "HIS:notanumber:A",
    "HIS::A",
    "HIS:12:",
    "12345",
    "HIS-12-A",
]


@pytest.mark.parametrize("atom", BAD_ATOM_FIELDS)
def test_unparseable_atom_is_fatal_in_strict_mode(atom):
    with pytest.raises(ConSurfParseError, match="refusing to guess"):
        parse_consurf(webserver(ws_row(1, "H", atom)), strict=True)


@pytest.mark.parametrize("atom", BAD_ATOM_FIELDS)
def test_unparseable_atom_is_quarantined_in_tolerant_mode(atom):
    result = parse_consurf(webserver(ws_row(1, "H", atom)), strict=False)

    assert len(result.records) == 1
    record = result.records[0]

    assert record.mapping_unverified is True
    assert record.pdb_residue is None
    assert WarningCode.MAPPING_UNVERIFIED in record.warning_codes

    # The crucial guarantee: it never reaches the join table, and it is not
    # confused with a residue ConSurf legitimately reported as unmapped.
    assert result.by_pdb_residue() == {}
    assert result.unmapped_records() == []
    assert result.unverified_records() == [record]


def test_explicit_dash_is_unmapped_not_unverified():
    result = parse_consurf(webserver(ws_row(1, "H", "-")), strict=True)
    record = result.records[0]
    assert record.pdb_residue is None
    assert record.mapping_unverified is False
    assert result.unmapped_records() == [record]


# ==========================================================================
# Identity uniqueness
# ==========================================================================

def test_duplicate_positions_are_fatal():
    with pytest.raises(ConSurfParseError, match="Duplicate residue positions"):
        parse_consurf(
            webserver(ws_row(1, "G", "GLY:1:A"), ws_row(1, "A", "ALA:2:A")),
            strict=True,
        )


def test_duplicate_pdb_identity_is_fatal():
    """Two records claiming the same residue would corrupt the join key."""
    with pytest.raises(ConSurfParseError, match="Duplicate PDB residue identities"):
        parse_consurf(
            webserver(ws_row(1, "G", "GLY:7:A"), ws_row(2, "A", "ALA:7:A")),
            strict=True,
        )


def test_duplicate_identity_downgrades_to_warning_in_tolerant_mode():
    result = parse_consurf(
        webserver(ws_row(1, "G", "GLY:7:A"), ws_row(2, "A", "ALA:7:A")),
        strict=False,
    )
    assert WarningCode.MALFORMED_ROW in result.warning_codes
    assert len(result.records) == 2


def test_same_number_different_chain_is_not_a_duplicate():
    result = parse_consurf(
        webserver(ws_row(1, "G", "GLY:7:A"), ws_row(2, "A", "ALA:7:B")),
        strict=True,
    )
    assert set(result.by_pdb_residue()) == {("A", 7, None), ("B", 7, None)}


# ==========================================================================
# Malformed numeric fields
# ==========================================================================

MALFORMED_ROWS = {
    "bad score": "1\tG\tGLY:1:A\tBAD\t7\t-0.8, -0.2\t9,6\tb\t\t20/20\tG 100%\n",
    "grade zero": "1\tG\tGLY:1:A\t-0.5\t0\t-0.8, -0.2\t9,6\tb\t\t20/20\tG 100%\n",
    "grade ten": "1\tG\tGLY:1:A\t-0.5\t10\t-0.8, -0.2\t9,6\tb\t\t20/20\tG 100%\n",
    "msa zero total": "1\tG\tGLY:1:A\t-0.5\t7\t-0.8, -0.2\t9,6\tb\t\t0/0\tG 100%\n",
    "msa present>total": "1\tG\tGLY:1:A\t-0.5\t7\t-0.8, -0.2\t9,6\tb\t\t30/20\tG 100%\n",
    "msa not a ratio": "1\tG\tGLY:1:A\t-0.5\t7\t-0.8, -0.2\t9,6\tb\t\tmany\tG 100%\n",
    "no msa field": "1\tG\tGLY:1:A\t-0.5\t7\t-0.8, -0.2\t9,6\tb\t\tG 100%\n",
    "truncated row": "1\tG\tGLY:1:A\n",
}


@pytest.mark.parametrize("row", MALFORMED_ROWS.values(), ids=list(MALFORMED_ROWS))
def test_malformed_row_is_fatal_in_strict_mode(row):
    with pytest.raises(ConSurfParseError):
        parse_consurf(io.StringIO(PREAMBLE + WEBSERVER_HEADER + row), strict=True)


@pytest.mark.parametrize("row", MALFORMED_ROWS.values(), ids=list(MALFORMED_ROWS))
def test_malformed_row_is_recorded_in_tolerant_mode(row):
    good = ws_row(2, "A", "ALA:2:A")
    result = parse_consurf(
        io.StringIO(PREAMBLE + WEBSERVER_HEADER + row + good), strict=False
    )

    assert [r.position for r in result.records] == [2]
    assert len(result.malformed_lines) == 1

    entry = result.malformed_lines[0]
    assert entry.line_number > 0
    assert entry.raw.strip()
    assert entry.reason


def test_inverted_confidence_interval_is_fatal():
    row = "1\tG\tGLY:1:A\t-0.5\t7\t0.9, -0.9\t9,6\tb\t\t20/20\tG 100%\n"
    with pytest.raises(ConSurfParseError, match="inverted"):
        parse_consurf(io.StringIO(PREAMBLE + WEBSERVER_HEADER + row), strict=True)


# ==========================================================================
# Headers
# ==========================================================================

def test_missing_header_is_fatal():
    with pytest.raises(ConSurfParseError, match="header"):
        parse_consurf(io.StringIO(ws_row(1, "G", "GLY:1:A")), strict=True)


def test_wrong_header_is_fatal():
    content = "NAME\tVALUE\n" + ws_row(1, "G", "GLY:1:A")
    with pytest.raises(ConSurfParseError, match="header"):
        parse_consurf(io.StringIO(content), strict=True)


def test_empty_file_is_fatal():
    with pytest.raises(ConSurfParseError):
        parse_consurf(io.StringIO(""), strict=True)


def test_header_only_file_is_fatal():
    with pytest.raises(ConSurfParseError, match="No ConSurf residue records"):
        parse_consurf(io.StringIO(WEBSERVER_HEADER), strict=True)


def test_missing_header_is_a_warning_in_tolerant_mode():
    result = parse_consurf(io.StringIO(ws_row(1, "G", "GLY:1:A")), strict=False)
    assert len(result.records) == 1
    assert WarningCode.MALFORMED_ROW in result.warning_codes


def test_missing_file_is_fatal(tmp_path):
    with pytest.raises(ConSurfParseError, match="does not exist"):
        parse_consurf(tmp_path / "nope.txt", strict=True)


def test_directory_is_fatal(tmp_path):
    with pytest.raises(ConSurfParseError, match="not a file"):
        parse_consurf(tmp_path, strict=True)


# ==========================================================================
# Grade-layer table
# ==========================================================================

def test_missing_grade_layers_warns_but_parses():
    """Without thresholds, scores cannot be compared across runs -- say so."""
    result = parse_consurf(raw(ws_row(1, "G", "GLY:1:A")), strict=True)
    assert result.grade_layers == []
    assert WarningCode.MISSING_GRADE_LAYERS in result.warning_codes


def test_grade_layer_lines_are_not_mistaken_for_data():
    result = parse_consurf(webserver(ws_row(1, "G", "GLY:1:A")), strict=True)
    assert len(result.records) == 1
    assert len(result.grade_layers) == 2


# ==========================================================================
# Preamble and trailer prose must never be parsed as data
# ==========================================================================

def test_trailing_footnote_is_ignored():
    footnote = (
        "\n*Below the confidence cut-off - The calculations for this site were "
        "performed on less than 6 non-gaped homologue sequences,\n"
        "or the confidence interval for the estimated score is equal to- or "
        "larger than- 4 color grades.\n"
    )
    content = PREAMBLE + WEBSERVER_HEADER + ws_row(1, "G", "GLY:1:A") + footnote
    result = parse_consurf(io.StringIO(content), strict=True)
    assert len(result.records) == 1
    assert result.malformed_lines == []


def test_blank_lines_are_ignored():
    content = (
        PREAMBLE
        + WEBSERVER_HEADER
        + "\n\n"
        + ws_row(1, "G", "GLY:1:A")
        + "\n   \n"
        + ws_row(2, "A", "ALA:2:A")
    )
    result = parse_consurf(io.StringIO(content), strict=True)
    assert [r.position for r in result.records] == [1, 2]


# ==========================================================================
# Whitespace and line endings
# ==========================================================================

def test_space_separated_rows_parse():
    """ConSurf documentation renders rows space-aligned; accept both."""
    content = (
        PREAMBLE
        + "POS SEQ ATOM SCORE COLOR CONFIDENCE INTERVAL B/E F/S MSA DATA "
        "RESIDUE VARIETY\n"
        + "  1  G  GLY:1:A  -0.500   7   -0.800, -0.200  9,6   b      20/20  G 100%\n"
    )
    result = parse_consurf(io.StringIO(content), strict=True)
    record = result.records[0]
    assert record.pdb_residue.key == ("A", 1, None)
    assert record.buried_exposed == "b"
    assert record.residue_variety["G"].percent == 100.0


def test_crlf_input_parses():
    content = (PREAMBLE + WEBSERVER_HEADER + ws_row(1, "G", "GLY:1:A")).replace(
        "\n", "\r\n"
    )
    result = parse_consurf(io.StringIO(content), strict=True)
    assert len(result.records) == 1
    assert "\r" not in result.records[0].raw_line


def test_mixed_line_endings_warn():
    content = PREAMBLE + WEBSERVER_HEADER + ws_row(1, "G", "GLY:1:A").replace(
        "\n", "\r\n"
    )
    result = parse_consurf(io.StringIO(content), strict=True)
    assert WarningCode.MIXED_LINE_ENDINGS in result.warning_codes


# ==========================================================================
# Low-confidence marker
# ==========================================================================

def test_low_confidence_marker_is_parsed():
    result = parse_consurf(
        webserver(ws_row(1, "G", "GLY:1:A", color="7*")), strict=True
    )
    record = result.records[0]
    assert record.grade == 7
    assert record.low_confidence is True
    assert WarningCode.LOW_CONFIDENCE in record.warning_codes


def test_absent_low_confidence_marker():
    result = parse_consurf(webserver(ws_row(1, "G", "GLY:1:A", color="7")), strict=True)
    assert result.records[0].low_confidence is False


def test_wide_confidence_interval_warns():
    result = parse_consurf(
        webserver(ws_row(1, "G", "GLY:1:A", confidence="-0.900, 0.900\t9,1")),
        strict=True,
    )
    record = result.records[0]
    assert record.confidence.grade_span == 8
    assert WarningCode.WIDE_CONFIDENCE_INTERVAL in record.warning_codes


def test_shallow_msa_warns():
    result = parse_consurf(
        webserver(ws_row(1, "G", "GLY:1:A", msa="3/20")), strict=True
    )
    assert WarningCode.SHALLOW_MSA in result.records[0].warning_codes


# ==========================================================================
# Residue variety edge cases
# ==========================================================================

def test_unparseable_variety_item_warns_and_is_skipped():
    result = parse_consurf(
        webserver(ws_row(1, "G", "GLY:1:A", variety="G 90%, ???nonsense")),
        strict=True,
    )
    record = result.records[0]
    assert set(record.residue_variety) == {"G"}
    assert WarningCode.UNPARSED_VARIETY_ITEM in record.warning_codes


def test_very_low_percentage_sum_warns():
    result = parse_consurf(
        webserver(ws_row(1, "G", "GLY:1:A", variety="G 40%")), strict=True
    )
    assert WarningCode.VARIETY_PERCENT_SUM_LOW in result.records[0].warning_codes


def test_truncation_does_not_warn():
    """98% is ordinary ConSurf rounding, not a problem."""
    result = parse_consurf(
        webserver(ws_row(1, "G", "GLY:1:A", variety="G 64%, A 20%, L 14%")),
        strict=True,
    )
    assert (
        WarningCode.VARIETY_PERCENT_SUM_LOW
        not in result.records[0].warning_codes
    )


def test_fractional_percentages_parse():
    result = parse_consurf(
        webserver(ws_row(1, "G", "GLY:1:A", variety="G 97.778%, A 2.222%")),
        strict=True,
    )
    variety = result.records[0].residue_variety
    assert variety["G"].percent == pytest.approx(97.778)
    assert variety["A"].percent == pytest.approx(2.222)


# ==========================================================================
# Model invariants
# ==========================================================================

def test_residue_frequency_rejects_contradictory_state():
    from WatCon.consurf import ResidueFrequency

    with pytest.raises(ValueError):
        ResidueFrequency(residue="G", percent=5.0, less_than_one_percent=True)

    with pytest.raises(ValueError):
        ResidueFrequency(residue="G", percent=50.0, is_full=True)


def test_lookup_on_empty_chain_is_not_in_file():
    from WatCon.consurf import LookupStatus

    result = parse_consurf(webserver(ws_row(1, "G", "GLY:1:A")), strict=True)
    assert result.lookup_structural("Z", 1).status is LookupStatus.NOT_IN_FILE
    assert result.covered_range("Z") is None

"""SYNTHETIC, specification-derived tests for ConSurf format variation.

IMPORTANT -- read before trusting these.

Everything in this module is constructed from ConSurf's published documentation
and from the community ConSurf-DB parsers, **not** from ConSurf output we hold.
We currently possess only Bayesian, protein, web-server-dialect runs.  So these
tests prove that the parser behaves as specified on formats we believe ConSurf
can emit; they do not prove that our reading of those formats is correct.

Paths covered here that remain UNVERIFIED against real output:

* maximum-likelihood runs (no CONFIDENCE INTERVAL columns)
* the ConSurf-DB / gradesPE dialect (``MET1:A``, no B/E or F/S)
* nucleotide (DNA/RNA) runs
* insertion codes and unusual chain identifiers

Real-output validation lives in ``test_consurf_parser.py``.
"""

from __future__ import annotations

import io

import pytest

from WatCon.consurf import (
    Alphabet,
    Dialect,
    Method,
    WarningCode,
    parse_consurf,
)

from .conftest import FIXTURES

# --------------------------------------------------------------------------
# Builders
# --------------------------------------------------------------------------

PREAMBLE = (
    "\t Amino Acid Conservation Scores\n"
    "\t=======================================\n"
    "\n"
    "The layers for assigning grades are as follows.\n"
    "from -1.000 to -0.500 the grade is 9\n"
    "from -0.500 to  0.000 the grade is 8\n"
)

WEBSERVER_HEADER = (
    "POS\tSEQ\tATOM\tSCORE\tCOLOR\tCONFIDENCE INTERVAL\t"
    "B/E\tF/S\tMSA DATA\tRESIDUE VARIETY\n"
)

CONSURFDB_HEADER = (
    "POS\tSEQ\t3LATOM\tSCORE\tCOLOR\tCONFIDENCE INTERVAL\t"
    "CONFIDENCE INTERVAL COLORS\tMSA DATA\tRESIDUE VARIETY\n"
)


def webserver(*rows: str) -> io.StringIO:
    return io.StringIO(PREAMBLE + WEBSERVER_HEADER + "".join(rows))


def consurfdb(*rows: str) -> io.StringIO:
    return io.StringIO(PREAMBLE + CONSURFDB_HEADER + "".join(rows))


def ws_row(
    pos: int,
    seq: str,
    atom: str,
    score: str = "-0.500",
    color: str = "7",
    confidence: str = "-0.800, -0.200\t9,6",
    be: str = "b",
    fs: str = "",
    msa: str = "20/20",
    variety: str = "",
) -> str:
    """A web-server dialect row.  Pass confidence='' for a maximum-likelihood row."""
    variety = variety or f"{seq} 100%"
    parts = [str(pos), seq, atom, score, color]
    if confidence:
        parts.append(confidence)
    parts += [be, fs, msa, variety]
    return "\t".join(parts) + "\n"


def db_row(
    pos: int,
    seq: str,
    atom: str,
    score: str = "-0.500",
    color: str = "7",
    msa: str = "20/20",
    variety: str = "",
) -> str:
    """A ConSurf-DB row: no B/E or F/S, confidence split across two fields."""
    variety = variety or f"{seq} 100%"
    return (
        "\t".join(
            [str(pos), seq, atom, score, color, "-0.800, -0.200", "9,6", msa, variety]
        )
        + "\n"
    )


# ==========================================================================
# Maximum likelihood: no confidence interval at all
# ==========================================================================

def test_maximum_likelihood_row_parses():
    """ConSurf gives no confidence interval under maximum likelihood."""
    result = parse_consurf(
        webserver(ws_row(1, "G", "ALA:1:A", confidence="")), strict=True
    )

    assert len(result.records) == 1
    assert result.records[0].confidence is None
    assert result.provenance.method is Method.MAXIMUM_LIKELIHOOD
    assert WarningCode.NO_CONFIDENCE_INTERVAL in result.warning_codes


def test_maximum_likelihood_without_optional_columns():
    row = "\t".join(["1", "G", "ALA:1:A", "-0.500", "7", "20/20", "G 100%"]) + "\n"
    result = parse_consurf(webserver(row), strict=True)

    record = result.records[0]
    assert record.confidence is None
    assert record.buried_exposed is None
    assert record.functional_structural is None
    assert result.provenance.method is Method.MAXIMUM_LIKELIHOOD


def test_bayesian_is_still_detected_alongside():
    result = parse_consurf(webserver(ws_row(1, "G", "ALA:1:A")), strict=True)
    assert result.provenance.method is Method.BAYESIAN
    assert result.records[0].confidence is not None
    assert WarningCode.NO_CONFIDENCE_INTERVAL not in result.warning_codes


# ==========================================================================
# ConSurf-DB dialect
# ==========================================================================

def test_consurfdb_dialect_is_detected_and_parsed():
    result = parse_consurf(consurfdb(db_row(1, "M", "MET1:A")), strict=True)

    assert result.provenance.dialect is Dialect.CONSURFDB
    record = result.records[0]
    assert record.pdb_residue.name3 == "MET"
    assert record.pdb_residue.seq == 1
    assert record.pdb_residue.chain == "A"
    assert record.confidence is not None
    # This dialect has no B/E or F/S columns at all.
    assert record.buried_exposed is None
    assert record.functional_structural is None


def test_consurfdb_negative_residue_number():
    result = parse_consurf(consurfdb(db_row(1, "L", "LEU-5:A")), strict=True)
    residue = result.records[0].pdb_residue
    assert (residue.name3, residue.seq, residue.chain) == ("LEU", -5, "A")


def test_consurfdb_missing_density_marker():
    """ConSurf-DB writes missing density as ___<pos>:<chain>."""
    result = parse_consurf(consurfdb(db_row(1, "M", "___1:A")), strict=True)

    record = result.records[0]
    assert record.pdb_residue is None
    assert record.mapping_unverified is False
    assert record in result.unmapped_records()


def test_webserver_dialect_is_the_default():
    result = parse_consurf(webserver(ws_row(1, "M", "MET:1:A")), strict=True)
    assert result.provenance.dialect is Dialect.WEBSERVER


def test_dialects_do_not_cross_match():
    """A colon-form identifier must not be read as the concatenated form."""
    from WatCon.consurf.dialects import parse_atom_field

    colon = parse_atom_field("THR:-2:A")
    assert (colon.name3, colon.seq, colon.chain, colon.icode) == ("THR", -2, "A", None)

    concat = parse_atom_field("MET1:A")
    assert (concat.name3, concat.seq, concat.chain, concat.icode) == (
        "MET",
        1,
        "A",
        None,
    )


# ==========================================================================
# Insertion codes and chain identifiers
# ==========================================================================

@pytest.mark.parametrize(
    "atom,expected",
    [
        ("HIS:100A:A", ("HIS", 100, "A", "A")),
        ("HIS:100:A", ("HIS", 100, None, "A")),
        ("HIS:-3B:A", ("HIS", -3, "B", "A")),
    ],
)
def test_insertion_codes_webserver(atom, expected):
    result = parse_consurf(webserver(ws_row(1, "H", atom)), strict=True)
    residue = result.records[0].pdb_residue
    assert (residue.name3, residue.seq, residue.icode, residue.chain) == expected


def test_insertion_code_is_part_of_the_join_key():
    result = parse_consurf(
        webserver(
            ws_row(1, "H", "HIS:100:A"),
            ws_row(2, "H", "HIS:100A:A"),
        ),
        strict=True,
    )
    keys = set(result.by_pdb_residue())
    assert keys == {("A", 100, None), ("A", 100, "A")}
    assert result.lookup_structural("A", 100).record.position == 1
    assert result.lookup_structural("A", 100, "A").record.position == 2


@pytest.mark.parametrize("chain", ["A", "a", "1", "AAA", "Z9"])
def test_chain_identifier_is_opaque(chain):
    """Never assume a chain is a single uppercase letter."""
    result = parse_consurf(webserver(ws_row(1, "H", f"HIS:10:{chain}")), strict=True)
    assert result.records[0].pdb_residue.chain == chain
    assert result.lookup_structural(chain, 10).record is not None


# ==========================================================================
# Sequence alphabet
# ==========================================================================

@pytest.mark.parametrize("symbol", ["X", "B", "Z", "J", "O"])
def test_ambiguity_codes_warn_but_never_abort(symbol):
    """A single unusual symbol must not cost us the file."""
    result = parse_consurf(
        webserver(ws_row(1, symbol, "CRO:66:A", variety=f"{symbol} 100%")),
        strict=True,
    )
    record = result.records[0]
    assert record.sequence_residue == symbol
    assert WarningCode.UNUSUAL_SEQ_SYMBOL in record.warning_codes


def test_standard_amino_acids_do_not_warn():
    result = parse_consurf(webserver(ws_row(1, "W", "TRP:5:A")), strict=True)
    assert WarningCode.UNUSUAL_SEQ_SYMBOL not in result.records[0].warning_codes


def test_nucleotide_run_is_detected():
    bases = "ACGTACGTACGT"
    rows = [
        ws_row(i, base, f"D{base}:{i}:A", variety=f"{base} 100%")
        for i, base in enumerate(bases, start=1)
    ]
    result = parse_consurf(webserver(*rows), strict=True)
    assert result.provenance.alphabet is Alphabet.NUCLEOTIDE


def test_small_files_are_not_mistaken_for_nucleotide_runs():
    """A/C/G/T are valid amino acids, so a few rows prove nothing."""
    rows = [ws_row(i, base, f"ALA:{i}:A") for i, base in enumerate("ACGT", start=1)]
    result = parse_consurf(webserver(*rows), strict=True)
    assert result.provenance.alphabet is Alphabet.PROTEIN


def test_protein_run_with_only_acgt_like_residues_stays_protein():
    rows = [ws_row(i, "W", f"TRP:{i}:A") for i in range(1, 12)]
    result = parse_consurf(webserver(*rows), strict=True)
    assert result.provenance.alphabet is Alphabet.PROTEIN


# ==========================================================================
# Optional column handling
# ==========================================================================

@pytest.mark.parametrize(
    "be,fs,expected",
    [
        ("b", "", ("b", None)),
        ("e", "", ("e", None)),
        ("e", "f", ("e", "f")),
        ("b", "s", ("b", "s")),
        ("", "", (None, None)),
    ],
)
def test_optional_flag_combinations(be, fs, expected):
    result = parse_consurf(
        webserver(ws_row(1, "H", "HIS:1:A", be=be, fs=fs)), strict=True
    )
    record = result.records[0]
    assert (record.buried_exposed, record.functional_structural) == expected


# ==========================================================================
# Residue variety syntax
# ==========================================================================

def test_bare_letter_and_explicit_100_agree():
    bare = parse_consurf(
        webserver(ws_row(1, "G", "GLY:1:A", variety="G")), strict=True
    ).records[0]
    explicit = parse_consurf(
        webserver(ws_row(1, "G", "GLY:1:A", variety="G 100%")), strict=True
    ).records[0]

    assert bare.residue_variety["G"].percent == 100.0
    assert explicit.residue_variety["G"].percent == 100.0
    # Only the bare form carries the "ConSurf omitted the percentage" marker.
    assert bare.residue_variety["G"].is_full is True
    assert bare.is_fully_conserved and explicit.is_fully_conserved


def test_mixed_variety_forms():
    result = parse_consurf(
        webserver(ws_row(1, "G", "GLY:1:A", variety="G 80%, A 19%, C <1%")),
        strict=True,
    )
    variety = result.records[0].residue_variety

    assert variety["G"].percent == 80.0
    assert variety["A"].percent == 19.0
    assert variety["C"].percent is None
    assert variety["C"].less_than_one_percent is True
    assert result.records[0].is_fully_conserved is False


def test_empty_variety_is_allowed():
    result = parse_consurf(
        webserver(ws_row(1, "G", "GLY:1:A", variety=" ")), strict=True
    )
    assert result.records[0].residue_variety == {}
    assert result.records[0].is_fully_conserved is False


# ===========================================================================
# The fixtures' line endings are test data, and git will destroy them
# ===========================================================================

def test_fixture_line_endings_survive_checkout():
    """Each fixture must contain the line endings its name claims.

    ``consurf/parser.py`` detects line endings from raw bytes, so
    ``test_line_ending_detection`` depends on exactly what git checked out --
    and git normalises line endings unless told not to.

    Without the ``-text`` rule in ``.gitattributes`` this breaks in both
    directions and never on the machine that committed it: the CRLF fixture is
    stored as LF and arrives as LF on Linux and macOS, while ``core.autocrlf``
    converts the LF fixtures to CRLF on a fresh Windows clone. The parser test
    then fails as "LineEnding.LF is not LineEnding.CRLF", which describes the
    symptom and not the cause.

    This checks the bytes directly so the failure says what to fix.
    """
    import glob
    import os

    fixtures = sorted(glob.glob(os.path.join(str(FIXTURES), "*.grades.txt")))
    assert fixtures, "no grades fixtures found"

    for path in fixtures:
        name = os.path.basename(path)
        data = open(path, "rb").read()
        crlf = data.count(b"\r\n")
        bare_lf = data.count(b"\n") - crlf

        if ".crlf." in name:
            assert crlf and not bare_lf, (
                "%s should be CRLF throughout but has CRLF=%d, bare LF=%d. "
                "git has normalised it -- check the '-text' rule for "
                "WatCon/data/consurf/fixtures/*.grades.txt in .gitattributes."
                % (name, crlf, bare_lf)
            )
        else:
            assert bare_lf and not crlf, (
                "%s should be LF throughout but has CRLF=%d, bare LF=%d. "
                "core.autocrlf has converted it on checkout -- check the "
                "'-text' rule for WatCon/data/consurf/fixtures/*.grades.txt "
                "in .gitattributes." % (name, crlf, bare_lf)
            )

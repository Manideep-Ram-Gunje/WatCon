"""Parser tests driven by REAL ConSurf output.

Every assertion in this module is anchored to a value observed in an actual
ConSurf run held under ``WatCon/data/consurf/fixtures``.  Specification-derived
and adversarial cases live in ``test_consurf_dialects.py`` and
``test_consurf_adversarial.py`` respectively, so that real-output validation
stays clearly separated from synthetic validation.
"""

from __future__ import annotations

import pytest

from WatCon.consurf import (
    Alphabet,
    Dialect,
    LineEnding,
    LookupStatus,
    Method,
    WarningCode,
    parse_consurf,
)

from .conftest import (
    ALL_FIXTURES,
    BRS_150,
    BRS_150_CRLF,
    O7W_45,
    P00648_50,
    P00648_150,
    REAL_FIXTURES,
)


def ids(paths):
    return [p.name for p in paths]


# ===========================================================================
# Baseline: every real file parses cleanly in strict mode
# ===========================================================================

@pytest.mark.parametrize("path", ALL_FIXTURES, ids=ids(ALL_FIXTURES))
def test_real_file_parses_strict(path):
    result = parse_consurf(path, strict=True)

    assert result.records
    assert result.malformed_lines == []
    assert all(record.position >= 1 for record in result.records)
    assert all(1 <= record.grade <= 9 for record in result.records)
    assert all(isinstance(record.score, float) for record in result.records)


@pytest.mark.parametrize("path", ALL_FIXTURES, ids=ids(ALL_FIXTURES))
def test_positions_are_unique(path):
    result = parse_consurf(path)
    positions = [record.position for record in result.records]
    assert len(positions) == len(set(positions))


@pytest.mark.parametrize("path", ALL_FIXTURES, ids=ids(ALL_FIXTURES))
def test_pdb_identities_are_unique(path):
    """The join key WatCon will use must be unique within a file."""
    result = parse_consurf(path)
    keys = [r.pdb_residue.key for r in result.mapped_records()]
    assert len(keys) == len(set(keys))


@pytest.mark.parametrize(
    "path,expected",
    [(BRS_150, 110), (P00648_150, 107), (P00648_50, 107), (O7W_45, 249)],
    ids=ids(REAL_FIXTURES),
)
def test_record_counts(path, expected):
    assert len(parse_consurf(path).records) == expected


# ===========================================================================
# Residue variety -- the bare-letter 100% defect
# ===========================================================================

@pytest.mark.parametrize("path", ALL_FIXTURES, ids=ids(ALL_FIXTURES))
def test_no_unparsed_variety_warnings(path):
    """Regression guard: bare-letter entries used to be silently discarded."""
    result = parse_consurf(path)
    offenders = [
        (r.position, r.residue_variety_raw)
        for r in result.records
        if WarningCode.UNPARSED_VARIETY_ITEM in r.warning_codes
    ]
    assert offenders == []


@pytest.mark.parametrize("path", ALL_FIXTURES, ids=ids(ALL_FIXTURES))
def test_every_record_has_residue_variety(path):
    """No real record should end up with an empty composition."""
    result = parse_consurf(path)
    empty = [r.position for r in result.records if not r.residue_variety]
    assert empty == []


def test_bare_letter_is_full_conservation():
    """7O7W POS 1 has RESIDUE VARIETY 'M', which means M at 100%."""
    record = parse_consurf(O7W_45).by_position()[1]

    assert record.residue_variety_raw == "M"
    assert set(record.residue_variety) == {"M"}

    entry = record.residue_variety["M"]
    assert entry.percent == 100.0
    assert entry.is_full is True
    assert entry.less_than_one_percent is False
    assert record.is_fully_conserved is True


@pytest.mark.parametrize(
    "path,expected",
    [(BRS_150, 3), (P00648_150, 8), (P00648_50, 15), (O7W_45, 41)],
    ids=ids(REAL_FIXTURES),
)
def test_fully_conserved_counts(path, expected):
    """These are exactly the records the old parser lost (117 of 1152)."""
    result = parse_consurf(path)
    assert sum(r.is_fully_conserved for r in result.records) == expected


def test_less_than_one_percent_from_real_data():
    result = parse_consurf(BRS_150)
    entries = [
        f
        for r in result.records
        for f in r.residue_variety.values()
        if f.less_than_one_percent
    ]
    assert entries, "1BRS is expected to contain '<1%' entries"
    assert all(f.percent is None for f in entries)
    assert all(f.is_full is False for f in entries)


def test_non_standard_residue_in_variety():
    """7O7W POS 106 carries 'X 2%' -- a non-standard residue in the MSA column."""
    record = parse_consurf(O7W_45).by_position()[106]
    assert record.residue_variety["X"].percent == 2.0


def test_percentages_need_not_sum_to_100():
    """ConSurf truncates; 7O7W POS 79 sums to 98%.  Never renormalise silently."""
    record = parse_consurf(O7W_45).by_position()[79]
    total = sum(f.percent for f in record.residue_variety.values() if f.percent)
    assert total == pytest.approx(98.0)


@pytest.mark.parametrize("path", ALL_FIXTURES, ids=ids(ALL_FIXTURES))
def test_seq_residue_appears_in_its_own_variety(path):
    """Self-consistency invariant: holds 249/249 for 7O7W.

    This is what shows the grades table -- not the sibling variety CSV -- is the
    internally correct artifact.
    """
    result = parse_consurf(path)
    bad = [
        r.position
        for r in result.records
        if r.sequence_residue not in r.residue_variety
    ]
    assert bad == []


# ===========================================================================
# Structural identity
# ===========================================================================

def test_negative_and_zero_residue_numbers():
    """7O7W carries an expression tag numbered HIS:-5:A .. PRO:0:A."""
    by_position = parse_consurf(O7W_45).by_position()

    first_tag = by_position[8].pdb_residue
    assert (first_tag.name3, first_tag.seq, first_tag.chain) == ("HIS", -5, "A")
    assert first_tag.icode is None

    zero = by_position[13].pdb_residue
    assert (zero.name3, zero.seq, zero.chain) == ("PRO", 0, "A")

    assert WarningCode.NONPOSITIVE_RESIDUE_NUMBER in by_position[8].warning_codes
    assert WarningCode.NONPOSITIVE_RESIDUE_NUMBER in by_position[13].warning_codes


@pytest.mark.parametrize(
    "path,position,expected",
    [
        (BRS_150, 3, ("VAL", 3, "A")),
        (BRS_150, 110, ("ARG", 110, "A")),
        # UniProt-numbered run: POS 1 maps to residue 48, an offset of +47.
        (P00648_150, 1, ("ALA", 48, "A")),
        (P00648_150, 107, ("PHE", 153, "A")),
        (P00648_50, 1, ("ALA", 48, "A")),
        (O7W_45, 8, ("HIS", -5, "A")),
        (O7W_45, 14, ("MET", 1, "A")),
    ],
)
def test_known_position_to_pdb_mappings(path, position, expected):
    """POS is never assumed equal to the PDB residue number."""
    residue = parse_consurf(path).by_position()[position].pdb_residue
    assert (residue.name3, residue.seq, residue.chain) == expected


def test_leading_unmapped_residues():
    result = parse_consurf(BRS_150)
    unmapped = [r.position for r in result.unmapped_records()]
    assert unmapped == [1, 2]
    assert all(r.pdb_residue is None for r in result.unmapped_records())
    assert result.unverified_records() == []


def test_interior_unmapped_residue():
    """P00648 POS 51 is unmapped in the MIDDLE of the chain, not at a terminus."""
    result = parse_consurf(P00648_150)
    unmapped = [r.position for r in result.unmapped_records()]
    assert unmapped == [51]

    positions = [r.position for r in result.records]
    assert min(positions) < 51 < max(positions)


def test_mapped_records_may_be_non_contiguous():
    result = parse_consurf(P00648_150)
    mapped = [r.position for r in result.mapped_records()]
    assert 50 in mapped and 52 in mapped and 51 not in mapped


def test_by_pdb_residue_round_trip():
    result = parse_consurf(O7W_45)
    lookup = result.by_pdb_residue()
    assert len(lookup) == len(result.mapped_records())
    for key, record in lookup.items():
        assert record.pdb_residue.key == key


# ===========================================================================
# POS discontinuity -- the omitted-residue signal
# ===========================================================================

def test_position_gap_is_detected():
    """7O7W's chromophore (PIA:68:A) is absent, leaving a 65 -> 69 jump."""
    result = parse_consurf(O7W_45)
    assert result.position_gaps == [(65, 69)]
    assert WarningCode.PDB_NUMBER_DISCONTINUITY in result.warning_codes


@pytest.mark.parametrize(
    "path", [BRS_150, P00648_150, P00648_50], ids=["1BRS", "P00648_150", "P00648_50"]
)
def test_no_spurious_position_gaps(path):
    assert parse_consurf(path).position_gaps == []


def test_lookup_structural_three_states():
    """The query direction WatCon needs: structure -> conservation."""
    result = parse_consurf(O7W_45)

    scored = result.lookup_structural("A", -5)
    assert scored.status is LookupStatus.SCORED
    assert scored.record.pdb_residue.name3 == "HIS"
    assert bool(scored) is True

    # The chromophore is present in the structure but ConSurf never scored it.
    omitted = result.lookup_structural("A", 68)
    assert omitted.status is LookupStatus.NO_SCORE_IN_FILE
    assert omitted.record is None
    assert bool(omitted) is False

    outside = result.lookup_structural("A", 9999)
    assert outside.status is LookupStatus.NOT_IN_FILE

    # 1BRS was run on chain A only; chain B exists in the structure.
    other_chain = parse_consurf(BRS_150).lookup_structural("B", 10)
    assert other_chain.status is LookupStatus.NOT_IN_FILE


# ===========================================================================
# Score / grade semantics
# ===========================================================================

@pytest.mark.parametrize("path", REAL_FIXTURES, ids=ids(REAL_FIXTURES))
def test_score_sign_convention(path):
    """Negative score == conserved; grade 9 == conserved.  The scales invert."""
    records = parse_consurf(path).records
    assert min(records, key=lambda r: r.score).grade == 9
    assert max(records, key=lambda r: r.score).grade == 1


@pytest.mark.parametrize("path", ALL_FIXTURES, ids=ids(ALL_FIXTURES))
def test_grade_layers_captured(path):
    result = parse_consurf(path)
    assert len(result.grade_layers) == 9
    assert [layer.grade for layer in result.grade_layers] == list(range(1, 10))
    assert WarningCode.MISSING_GRADE_LAYERS not in result.warning_codes


def test_grade_layers_differ_between_runs():
    """Scores are per-run normalised, so cross-run comparison needs these."""
    o7w = max(layer.score_to for layer in parse_consurf(O7W_45).grade_layers)
    brs = max(layer.score_to for layer in parse_consurf(BRS_150).grade_layers)
    assert o7w == pytest.approx(3.395)
    assert brs == pytest.approx(2.709)
    assert o7w != brs


@pytest.mark.parametrize("path", ALL_FIXTURES, ids=ids(ALL_FIXTURES))
def test_grade_layers_bracket_observed_scores(path):
    result = parse_consurf(path)
    low = min(layer.score_from for layer in result.grade_layers)
    high = max(layer.score_to for layer in result.grade_layers)
    for record in result.records:
        assert low <= record.score <= high


# ===========================================================================
# Confidence interval
# ===========================================================================

@pytest.mark.parametrize("path", ALL_FIXTURES, ids=ids(ALL_FIXTURES))
def test_confidence_intervals_present_and_ordered(path):
    """All four real datasets are Bayesian runs."""
    result = parse_consurf(path)
    assert result.provenance.method is Method.BAYESIAN

    for record in result.records:
        confidence = record.confidence
        assert confidence is not None
        assert confidence.score_lower <= confidence.score_upper
        assert 1 <= confidence.grade_at_lower <= 9
        assert 1 <= confidence.grade_at_upper <= 9


def test_confidence_grade_bounds_are_not_ordered_low_to_high():
    """Grades invert scores, so grade_at_lower >= grade_at_upper is normal.

    7O7W POS 11 has interval 2.038..3.510 with grades 2,1.
    """
    confidence = parse_consurf(O7W_45).by_position()[11].confidence
    assert confidence.grade_at_lower == 2
    assert confidence.grade_at_upper == 1
    assert confidence.grade_at_lower > confidence.grade_at_upper
    assert confidence.grade_span == 1


# ===========================================================================
# MSA support and reliability
# ===========================================================================

@pytest.mark.parametrize(
    "path,total",
    [(BRS_150, 150), (P00648_150, 150), (P00648_50, 50), (O7W_45, 45)],
    ids=ids(REAL_FIXTURES),
)
def test_msa_totals(path, total):
    result = parse_consurf(path)
    assert result.provenance.msa_total == total
    for record in result.records:
        assert record.msa.total == total
        assert 0 <= record.msa.present <= total
        assert 0.0 <= record.msa_fraction <= 1.0


@pytest.mark.parametrize("path", ALL_FIXTURES, ids=ids(ALL_FIXTURES))
def test_low_confidence_records_are_flagged(path):
    for record in parse_consurf(path).records:
        if record.low_confidence:
            assert WarningCode.LOW_CONFIDENCE in record.warning_codes


def test_shallow_msa_warning():
    """7O7W is a shallow run (45 sequences) and has 6 low-confidence positions."""
    result = parse_consurf(O7W_45)
    assert sum(r.low_confidence for r in result.records) == 6
    shallow = [
        r for r in result.records if WarningCode.SHALLOW_MSA in r.warning_codes
    ]
    assert shallow
    assert all(r.msa.present < 6 for r in shallow)


# ===========================================================================
# Optional B/E and F/S columns
# ===========================================================================

@pytest.mark.parametrize("path", ALL_FIXTURES, ids=ids(ALL_FIXTURES))
def test_optional_flag_domains(path):
    result = parse_consurf(path)
    for record in result.records:
        assert record.buried_exposed in (None, "b", "e")
        assert record.functional_structural in (None, "f", "s")

    # Every real row carries B/E; F/S only on highly conserved positions.
    assert all(r.buried_exposed is not None for r in result.records)
    assert any(r.functional_structural is None for r in result.records)


def test_blank_functional_structural_is_none():
    record = parse_consurf(O7W_45).by_position()[11]
    assert record.buried_exposed == "e"
    assert record.functional_structural is None


# ===========================================================================
# Provenance
# ===========================================================================

@pytest.mark.parametrize("path", ALL_FIXTURES, ids=ids(ALL_FIXTURES))
def test_provenance_of_real_files(path):
    provenance = parse_consurf(path).provenance
    assert provenance.dialect is Dialect.WEBSERVER
    assert provenance.method is Method.BAYESIAN
    assert provenance.alphabet is Alphabet.PROTEIN
    assert provenance.source


def test_line_ending_detection():
    assert parse_consurf(BRS_150).provenance.line_ending is LineEnding.LF
    assert parse_consurf(BRS_150_CRLF).provenance.line_ending is LineEnding.CRLF


def test_crlf_and_lf_parse_identically():
    """The two 1BRS fixtures differ only in line endings."""
    lf = parse_consurf(BRS_150)
    crlf = parse_consurf(BRS_150_CRLF)

    assert len(lf.records) == len(crlf.records)
    for left, right in zip(lf.records, crlf.records):
        assert left.position == right.position
        assert left.score == right.score
        assert left.grade == right.grade
        assert left.residue_variety == right.residue_variety
        assert (left.pdb_residue is None) == (right.pdb_residue is None)
        if left.pdb_residue is not None:
            assert left.pdb_residue.key == right.pdb_residue.key


def test_chains_reported():
    assert parse_consurf(O7W_45).chains() == ["A"]
    assert parse_consurf(BRS_150).covered_range("A") == (3, 110)
    assert parse_consurf(BRS_150).covered_range("B") is None


# ===========================================================================
# Serialisation
# ===========================================================================

@pytest.mark.parametrize("path", REAL_FIXTURES, ids=ids(REAL_FIXTURES))
def test_records_serialise(path):
    import json

    for record in parse_consurf(path).records:
        data = record.to_dict()
        assert data["position"] == record.position
        assert data["grade"] == record.grade
        assert data["score"] == record.score
        json.dumps(data)  # must be JSON-clean

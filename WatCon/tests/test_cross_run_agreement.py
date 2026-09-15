"""How well independent ConSurf runs agree -- measured with residues paired correctly.

For a long time this project documented that "two runs of the same barnase
sequence from different starting structures agree at only rho ~= 0.37", and
treated that as the noise floor for any family-level claim.

It was an artefact. The 1BRS run is numbered from residue 3 and the P00648 run
from residue 48, and the two had been paired by ConSurf position. Half the
position-paired residues were different amino acids. Paired by sequence, the
same runs agree at 0.97 (depth 150) and 0.94 (depth 50).

These tests pin both the correct figure and the artefact, so the mistake cannot
come back quietly: pairing runs of differently numbered structures must go
through sequence, never through position.
"""

from __future__ import annotations

import pytest

pytest.importorskip("Bio", reason="sequence pairing needs Biopython")
from scipy.stats import spearmanr

from WatCon.consurf import parse_consurf
from WatCon.evolutionary import CROSS_RUN_AGREEMENT

from .conftest import BRS_150, P00648_50, P00648_150


def _records(path):
    return parse_consurf(path, strict=True).records


def _pair_by_sequence(first, second):
    from Bio.Align import PairwiseAligner

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5
    a = "".join(r.sequence_residue for r in first)
    b = "".join(r.sequence_residue for r in second)
    alignment = aligner.align(a, b)[0]
    return [(first[i], second[j])
            for (a0, a1), (b0, b1) in zip(*alignment.aligned)
            for i, j in zip(range(a0, a1), range(b0, b1))]


def _pair_by_position(first, second):
    by_position = {r.position: r for r in second}
    return [(r, by_position[r.position]) for r in first if r.position in by_position]


def _rho(pairs, attribute="score"):
    return spearmanr([getattr(x, attribute) for x, _ in pairs],
                     [getattr(y, attribute) for _, y in pairs]).correlation


def _identity(pairs):
    return sum(x.sequence_residue == y.sequence_residue for x, y in pairs) / len(pairs)


def test_same_sequence_two_depths_matches_the_recorded_constant():
    pairs = _pair_by_sequence(_records(P00648_150), _records(P00648_50))
    assert len(pairs) == 107
    assert _identity(pairs) == 1.0
    assert _rho(pairs) == pytest.approx(CROSS_RUN_AGREEMENT["score_spearman"], abs=0.001)


@pytest.mark.parametrize("other, expected", [(P00648_150, 0.969), (P00648_50, 0.941)])
def test_different_starting_structures_agree_closely_when_paired_by_sequence(other, expected):
    pairs = _pair_by_sequence(_records(BRS_150), _records(other))
    assert len(pairs) == 106
    assert _identity(pairs) == 1.0
    assert _rho(pairs) == pytest.approx(expected, abs=0.0015)


@pytest.mark.parametrize("other", [P00648_150, P00648_50])
def test_pairing_by_position_produced_the_old_037_floor(other):
    """The artefact, pinned: half the position-paired residues differ."""
    pairs = _pair_by_position(_records(BRS_150), _records(other))
    assert _identity(pairs) == pytest.approx(0.50, abs=0.02)
    assert _rho(pairs) < 0.40

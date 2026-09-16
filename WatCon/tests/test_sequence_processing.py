"""Sequence helpers in :mod:`WatCon.sequence_processing`.

This module was an empty placeholder. What it covers now is ``seq_similarity``,
which was both untested and the package's last use of ``Bio.pairwise2`` -- an
API Biopython has deprecated and intends to remove. It is now
``PairwiseAligner`` configured identically, and these values were taken from the
old implementation so the swap cannot have changed an answer silently.
"""

from __future__ import annotations

import pytest

from WatCon.sequence_processing import seq_similarity

# ===========================================================================
# seq_similarity
#
# Untested until now, and the package's last use of Bio.pairwise2, which
# Biopython has deprecated and intends to remove. It is now PairwiseAligner
# configured identically -- global, match 1, mismatch 0, no gap penalty -- and
# these pin the behaviour so the swap cannot have changed the answer silently.
# ===========================================================================

def test_identical_sequences_are_completely_similar():
    assert seq_similarity("ACDEFGHIKL", "ACDEFGHIKL") == 1.0


def test_sequences_sharing_nothing_are_not_similar():
    assert seq_similarity("WWWWWW", "YYYYYY") == 0.0


def test_a_single_substitution_and_a_single_deletion():
    """Exact values, carried over from the Bio.pairwise2 implementation."""
    assert seq_similarity("ACDEFGHIKL", "ACDEFGHIKM") == pytest.approx(9 / 11)
    assert seq_similarity("ACDEFGHIKL", "ACDEGHIKL") == pytest.approx(0.9)
    assert seq_similarity("MKVLAT", "MKVLATTT") == pytest.approx(0.75)


def test_similarity_is_symmetric():
    assert seq_similarity("MKVLAT", "MKVLATTT") == seq_similarity("MKVLATTT", "MKVLAT")

"""Phase A: residue -> alignment-column mapping, published WatCon vs ours.

What is being measured
----------------------
Published WatCon assigns a residue's alignment column with

    MSA_index = MSA_indices[molecule.resid-1]

at `generate_static_networks.py` lines 714, 798 and 816 of the released version
(commit 89e992d).

`MSA_indices` comes from `sequence_processing.generate_msa_alignment` and is
**positional**: element *k* is the alignment column of the *k*-th residue of that
structure's sequence. So the correct index for the residue at ordinal *k* is
`MSA_indices[k]`, and `MSA_indices[resid-1]` is right only when
``resid == k + 1`` -- that is, only when the structure's residue numbering starts
at 1 and contains no gaps.

Crystal structures routinely violate both: disordered termini, expression tags
numbered from zero or negative, and construct numbering that starts elsewhere.

Every residue therefore falls into one of four outcomes:

  correct          resid - 1 == k
  wrong, in range  0 <= resid-1 < len, but != k  -> a real column, silently wrong
  wrong, wrapped   resid - 1 < 0                 -> Python indexes from the end
  out of range     resid - 1 >= len              -> IndexError

Only the last is visible to a user. The other two produce a plausible number.

Ours routes the same question through `WatCon.residue_index.ResidueIndex`, which
pairs residues with columns by ordinal and refuses to build when the two lengths
disagree -- so it is correct by construction, and this script asserts that rather
than assuming it.

Usage
-----
    python scripts/compare_mapping.py <directory of .pdb files> [--label NAME]
"""

from __future__ import annotations

import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
WATCON = os.path.join(os.path.dirname(os.path.dirname(ROOT)),
                      "tools", "WatCon_ConSurf")
sys.path.insert(0, WATCON)

from WatCon.residue_index import residues_from_pdb_file   # noqa: E402

CORRECT, WRONG_IN_RANGE, WRONG_WRAPPED, OUT_OF_RANGE = (
    "correct", "wrong_in_range", "wrong_wrapped", "out_of_range")


def classify(residues):
    """Outcome of ``MSA_indices[resid-1]`` for each residue, in order.

    ``len(MSA_indices)`` equals the number of residues in the structure, because
    the list is built one entry per residue of that structure's sequence.
    """
    n = len(residues)
    outcomes = []
    for ordinal, residue in enumerate(residues):
        index = residue.resid - 1
        if index == ordinal:
            outcomes.append(CORRECT)
        elif index < 0:
            outcomes.append(WRONG_WRAPPED)
        elif index >= n:
            outcomes.append(OUT_OF_RANGE)
        else:
            outcomes.append(WRONG_IN_RANGE)
    return outcomes


def ours_is_correct(path, residues):
    """Assert our mapping is right, rather than assuming it.

    ``ResidueIndex`` is handed the same positional list; asking it for each
    residue's column must return that residue's own ordinal entry.
    """
    from WatCon.residue_index import ResidueIndex

    columns = list(range(len(residues)))          # stand-in alignment columns
    index = ResidueIndex(residues, columns)
    wrong = 0
    for ordinal, residue in enumerate(residues):
        got = index.msa_column(residue.chain, residue.resid, residue.icode)
        if got != columns[ordinal]:
            wrong += 1
    return wrong


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory")
    parser.add_argument("--label", default=None)
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()

    label = args.label or os.path.basename(args.directory.rstrip("/\\"))
    files = sorted(f for f in os.listdir(args.directory) if f.lower().endswith(".pdb"))
    if args.limit:
        files = files[:args.limit]
    if not files:
        print("no .pdb files in %s" % args.directory)
        return 1

    print("=" * 76)
    print("Residue -> alignment column: published WatCon (89e992d) vs ours")
    print("dataset: %s   %d structures" % (label, len(files)))
    print("=" * 76)

    rows = []
    totals = {CORRECT: 0, WRONG_IN_RANGE: 0, WRONG_WRAPPED: 0, OUT_OF_RANGE: 0}
    ours_wrong_total = 0
    perfect_structures = 0

    for name in files:
        path = os.path.join(args.directory, name)
        residues = residues_from_pdb_file(path)
        if not residues:
            continue
        outcomes = classify(residues)
        counts = {k: outcomes.count(k) for k in totals}
        for k in totals:
            totals[k] += counts[k]
        if counts[CORRECT] == len(residues):
            perfect_structures += 1

        ours_wrong = ours_is_correct(path, residues)
        ours_wrong_total += ours_wrong

        rows.append({
            "structure": os.path.splitext(name)[0],
            "n_residues": len(residues),
            "first_resid": residues[0].resid,
            "chains": len({r.chain for r in residues}),
            "published_correct": counts[CORRECT],
            "published_wrong_in_range": counts[WRONG_IN_RANGE],
            "published_wrong_wrapped": counts[WRONG_WRAPPED],
            "published_out_of_range": counts[OUT_OF_RANGE],
            "published_fraction_correct": round(counts[CORRECT] / len(residues), 4),
            "ours_wrong": ours_wrong,
        })

    total_residues = sum(r["n_residues"] for r in rows)
    print()
    print("PUBLISHED WatCon")
    print("  residues mapped correctly     %7d / %d  (%.1f%%)"
          % (totals[CORRECT], total_residues, 100 * totals[CORRECT] / total_residues))
    print("  silently wrong, in range      %7d" % totals[WRONG_IN_RANGE])
    print("  silently wrong, negative wrap %7d" % totals[WRONG_WRAPPED])
    print("  would raise IndexError        %7d" % totals[OUT_OF_RANGE])
    print("  structures mapped perfectly   %7d / %d" % (perfect_structures, len(rows)))
    print()
    print("OURS")
    print("  residues mapped incorrectly   %7d / %d" % (ours_wrong_total, total_residues))

    print()
    print("Mechanism: correct only where residue numbering starts at 1 with no gaps.")
    by_first = {}
    for r in rows:
        by_first.setdefault(r["first_resid"], []).append(r["published_fraction_correct"])
    print("  %-14s %-12s %s" % ("first resid", "structures", "mean fraction correct"))
    for first in sorted(by_first)[:12]:
        v = by_first[first]
        print("  %-14s %-12d %.3f" % (first, len(v), sum(v) / len(v)))

    out = os.path.join(ROOT, "results", "mapping_%s.csv" % label)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print()
    print("Wrote %s" % out)

    # A worked example, so the numbers can be checked by hand.
    example = next((r for r in rows if r["published_fraction_correct"] < 1.0), rows[0])
    path = os.path.join(args.directory, example["structure"] + ".pdb")
    residues = residues_from_pdb_file(path)
    print()
    print("Worked example -- %s (first residue numbered %d):"
          % (example["structure"], residues[0].resid))
    print("  %-8s %-8s %-14s %-14s %s"
          % ("ordinal", "resid", "correct slot", "published slot", "same?"))
    for ordinal in range(min(6, len(residues))):
        r = residues[ordinal]
        print("  %-8d %-8d %-14d %-14d %s"
              % (ordinal, r.resid, ordinal, r.resid - 1,
                 "yes" if r.resid - 1 == ordinal else "NO"))
    return 0


if __name__ == "__main__":
    sys.exit(main())

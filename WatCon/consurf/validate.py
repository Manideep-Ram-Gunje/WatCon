"""Diagnostics CLI for ConSurf grades files.

    python -m WatCon.consurf.validate FILE [FILE ...] [--tolerant] [--json]
    python -m WatCon.consurf.validate FILE --pdb ANNOTATED.pdb

With no arguments it validates every fixture under ``WatCon/data/consurf/fixtures``.

``--pdb`` cross-checks the parsed grades against ConSurf's own annotated PDB,
in which the grade is written into the B-factor column.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path
from typing import List, Optional

from .crosscheck import check_against_annotated_pdb
from .parser import parse_consurf


def _fixture_dir() -> Path:
    return Path(__file__).resolve().parents[1] / "data" / "consurf" / "fixtures"


def summarise(path: Path, strict: bool = True) -> dict:
    """Return a JSON-friendly summary of one grades file."""
    result = parse_consurf(path, strict=strict)
    codes = Counter(
        code for record in result.records for code in record.warning_codes
    )
    return {
        "file": path.name,
        "dialect": result.provenance.dialect.value,
        "method": result.provenance.method.value,
        "alphabet": result.provenance.alphabet.value,
        "line_ending": result.provenance.line_ending.value,
        "msa_total": result.provenance.msa_total,
        "records": len(result.records),
        "mapped": len(result.mapped_records()),
        "unmapped": len(result.unmapped_records()),
        "unverified": len(result.unverified_records()),
        "fully_conserved": sum(r.is_fully_conserved for r in result.records),
        "low_confidence": sum(r.low_confidence for r in result.records),
        "grade_layers": len(result.grade_layers),
        "score_range": [
            min(r.score for r in result.records),
            max(r.score for r in result.records),
        ],
        "chains": result.chains(),
        "position_gaps": result.position_gaps,
        "file_warnings": result.warning_codes,
        "record_warning_counts": dict(codes),
        "malformed_lines": len(result.malformed_lines),
    }


def _print_report(summary: dict) -> None:
    print("=" * 68)
    print(f"  {summary['file']}")
    print("=" * 68)
    keys = [
        "dialect",
        "method",
        "alphabet",
        "line_ending",
        "msa_total",
        "records",
        "mapped",
        "unmapped",
        "unverified",
        "fully_conserved",
        "low_confidence",
        "grade_layers",
        "score_range",
        "chains",
        "position_gaps",
        "malformed_lines",
    ]
    for key in keys:
        print(f"  {key:<20} {summary[key]}")
    if summary["file_warnings"]:
        print(f"  {'file_warnings':<20} {summary['file_warnings']}")
    if summary["record_warning_counts"]:
        print(f"  {'record_warnings':<20} {summary['record_warning_counts']}")
    print()


def main(argv: Optional[List[str]] = None) -> int:
    cli = argparse.ArgumentParser(
        prog="python -m WatCon.consurf.validate",
        description="Validate and summarise ConSurf grades files.",
    )
    cli.add_argument("files", type=Path, nargs="*", help="grades file(s)")
    cli.add_argument(
        "--tolerant",
        action="store_true",
        help="Skip malformed rows instead of failing.",
    )
    cli.add_argument("--json", action="store_true", help="Emit JSON.")
    cli.add_argument(
        "--pdb",
        type=Path,
        default=None,
        help="Cross-check against a ConSurf annotated PDB (B-factor = grade).",
    )
    args = cli.parse_args(argv)

    paths = args.files or sorted(_fixture_dir().glob("*.grades.txt"))
    if not paths:
        print("No grades files found.", file=sys.stderr)
        return 1

    summaries = [summarise(path, strict=not args.tolerant) for path in paths]

    if args.json:
        print(json.dumps(summaries, indent=2))
    else:
        for summary in summaries:
            _print_report(summary)

    if args.pdb is not None:
        if len(paths) != 1:
            print("--pdb requires exactly one grades file.", file=sys.stderr)
            return 2
        result = parse_consurf(paths[0], strict=not args.tolerant)
        report = check_against_annotated_pdb(result, args.pdb)
        print(f"cross-check vs {args.pdb.name}: {report.summary()}")
        if not report.ok:
            print(f"  mismatches: {report.mismatches[:10]}", file=sys.stderr)
            return 3

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

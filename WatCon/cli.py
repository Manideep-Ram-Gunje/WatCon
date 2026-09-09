"""The ``watcon`` command.

    watcon prepare  --input-dir raw/ --out-dir prepared/
    watcon run      --input input.txt [--analysis analysis.txt]
    watcon validate --consurf FILE...
    watcon demo

``run`` delegates to the existing :mod:`WatCon.WatCon` entry points, so
``python WatCon/WatCon.py --input ...`` keeps working exactly as before; this
only puts a discoverable name in front of it and adds the steps that were
previously left to the user.

Imports are deliberately deferred into each subcommand.  ``watcon --help`` and
``watcon validate`` should not need MDAnalysis, matplotlib and scikit-learn to
be importable -- the ConSurf parser is pure standard library, and a user
diagnosing a broken grades file is often the one whose environment is broken.
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import List, Optional


def _version() -> str:
    try:
        from ._version import __version__
        return __version__
    except Exception:                       # noqa: BLE001 - version is cosmetic
        return "unknown"


# ---------------------------------------------------------------------------
# prepare
# ---------------------------------------------------------------------------

def cmd_prepare(args) -> int:
    from .prepare import PreparationError, prepare_directory

    try:
        report = prepare_directory(
            args.input_dir,
            args.out_dir,
            reference=args.reference,
            chain_label=args.chain_label,
            min_identity=args.min_identity,
            water_cutoff=args.water_cutoff,
        )
    except PreparationError as error:
        print("error: %s" % error, file=sys.stderr)
        return 1

    if args.report:
        report.write_csv(args.report)
        print("Wrote %s" % args.report)

    # Rejections are the point of this step: say so again at the end, where it
    # will not scroll past.
    if report.rejected:
        print()
        print("%d structure(s) were NOT prepared:" % len(report.rejected))
        for outcome in report.rejected:
            print("  %-8s %s" % (outcome.pdb_id, outcome.status))
            if outcome.note:
                print("           %s" % outcome.note)
    return 0


# ---------------------------------------------------------------------------
# run
# ---------------------------------------------------------------------------

def cmd_run(args) -> int:
    import pickle

    from .WatCon import parse_analysis, parse_inputs, run_watcon, run_watcon_postanalysis

    if args.input is None and args.analysis is None:
        print("error: give --input, --analysis, or both", file=sys.stderr)
        return 1

    if args.input is not None:
        structure_type, kwargs = parse_inputs(args.input)
        results = run_watcon(structure_type, kwargs)
        os.makedirs("watcon_output", exist_ok=True)
        out = os.path.join("watcon_output", args.name + ".pkl")
        with open(out, "wb") as handle:
            pickle.dump(results, handle)
        print("Wrote %s" % out)

    if args.analysis is not None:
        run_watcon_postanalysis(**parse_analysis(args.analysis))
    return 0


# ---------------------------------------------------------------------------
# validate
# ---------------------------------------------------------------------------

def cmd_validate(args) -> int:
    from .consurf.validate import main as validate_main

    argv: List[str] = list(args.consurf)
    if args.tolerant:
        argv.append("--tolerant")
    if args.json:
        argv.append("--json")
    if args.pdb:
        argv += ["--pdb", args.pdb]
    return validate_main(argv)


# ---------------------------------------------------------------------------
# demo
# ---------------------------------------------------------------------------

def cmd_demo(args) -> int:
    from .demo import run_demo
    return run_demo(args.out_dir, keep=args.keep)


# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="watcon",
        description=(
            "Conserved water network analysis with ConSurf evolutionary "
            "conservation. Extends WatCon (Brownless & Kamerlin, JACS Au 2025)."
        ),
    )
    parser.add_argument("--version", action="version",
                        version="watcon-consurf %s" % _version())
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")

    # -- prepare ------------------------------------------------------------
    prepare = subparsers.add_parser(
        "prepare",
        help="align a folder of structures into one frame, keeping their waters",
        description=(
            "Pick the right chain out of each structure, superpose them all onto "
            "a reference, and carry the waters along. Structures that do not "
            "match are rejected with a reason, never dropped silently."
        ),
    )
    prepare.add_argument("--input-dir", required=True, help="folder of raw .pdb files")
    prepare.add_argument("--out-dir", required=True,
                         help="destination (CLEARED, not merged)")
    prepare.add_argument("--reference", default=None,
                         help="PDB id or filename to superpose onto "
                              "(default: first in sorted order)")
    prepare.add_argument("--chain-label", default="A",
                         help="chain letter for every output (default: A)")
    prepare.add_argument("--min-identity", type=float, default=0.80,
                         help="minimum agreement with the reference (default: 0.80)")
    prepare.add_argument("--water-cutoff", type=float, default=5.0,
                         help="keep waters within this distance, Angstrom "
                              "(default: 5.0)")
    prepare.add_argument("--report", default=None, help="write a per-structure CSV")
    prepare.set_defaults(func=cmd_prepare)

    # -- run ----------------------------------------------------------------
    run = subparsers.add_parser(
        "run", help="run WatCon from an input file",
        description="Identical to `python WatCon/WatCon.py`, which still works.",
    )
    run.add_argument("--input", default=None, help="network-building input file")
    run.add_argument("--analysis", default=None, help="post-analysis input file")
    run.add_argument("--name", default="results", help="output basename")
    run.set_defaults(func=cmd_run)

    # -- validate -----------------------------------------------------------
    validate = subparsers.add_parser(
        "validate", help="check ConSurf grades files before using them",
        description=(
            "Parse ConSurf grades files and report what they contain. Run this "
            "before a real analysis: a file that parses is not necessarily one "
            "that describes your structure."
        ),
    )
    validate.add_argument("--consurf", nargs="+", required=True, metavar="FILE",
                          help="ConSurf grades file(s)")
    validate.add_argument("--pdb", default=None,
                          help="cross-check against ConSurf's annotated PDB, "
                               "in which the grade is written into the B-factors")
    validate.add_argument("--tolerant", action="store_true",
                          help="warn instead of failing on unusual records")
    validate.add_argument("--json", action="store_true", help="machine-readable output")
    validate.set_defaults(func=cmd_validate)

    # -- demo ---------------------------------------------------------------
    demo = subparsers.add_parser(
        "demo", help="run the whole pipeline on the bundled barnase example",
        description=(
            "Prepare, build networks, cluster, join to ConSurf conservation and "
            "write a report plus PyMOL projections, using data shipped with the "
            "package. Needs no network access."
        ),
    )
    demo.add_argument("--out-dir", default="watcon_demo", help="output directory")
    demo.add_argument("--keep", action="store_true",
                      help="reuse an existing output directory instead of clearing it")
    demo.set_defaults(func=cmd_demo)

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if getattr(args, "func", None) is None:
        parser.print_help()
        return 1
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())

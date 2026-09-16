"""Assemble the fifteen-protein PTP family and write its members files.

Fetches the twenty-four RCSB entries, lays them out one directory per protein,
and writes the two members files `watcon family` reads -- the full fifteen and
the fourteen without the pseudophosphatase, so the comparison in FINDINGS.md can
be reproduced.

Structures are fetched rather than committed: they are RCSB's to serve, and the
set is 9 MB. **Native numbering matters** -- the authors' own copies are
renumbered from 1 over resolved residues, which collapses crystal gaps, so no
constant offset recovers the sequence and `enforce_identity` refuses them. See
FINDINGS.md section 6.

Run from the repository root::

    python experiments/benchmark/scripts/build_family15.py
    watcon family \\
        --members experiments/benchmark/ptp_family15/members15.tsv \\
        --alignment WatCon/data/examples/ptp_family/ptp_family_alignment.pir
"""

from __future__ import annotations

import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BENCHMARK = os.path.dirname(HERE)
REPO = os.path.dirname(os.path.dirname(BENCHMARK))

OUT = os.path.join(BENCHMARK, "ptp_family15")
STRUCTURES = os.path.join(OUT, "structures")
RAW = os.path.join(OUT, "raw")
FIXTURES = os.path.join("WatCon", "data", "consurf", "fixtures")

#: protein -> (ConSurf run, reference structure, structures)
#:
#: One run per protein, all chain A, all webserver/Bayesian/150 homologues --
#: identical settings are what make them poolable. 2QEP is labelled PTPN2 in the
#: authors' dataset; the entry is PTPRN2 (Q92932) and real PTPN2 is P17706.
MEMBERS = [
    ("PTPN1",  "1AAX", "2F71", ["2F71", "8U1E"]),
    ("PTPN3",  "4S0G", "4S0G", ["4S0G", "2B49"]),
    ("PTPN4",  "2I75", "2I75", ["2I75"]),
    ("PTPN5",  "8SLS", "8SLS", ["8SLS"]),
    ("PTPN6",  "4GRZ", "4GRZ", ["4GRZ", "4HJP"]),
    ("PTPN7",  "1ZC0", "1ZC0", ["1ZC0", "3O4U"]),
    ("PTPN9",  "6KZQ", "6KZQ", ["6KZQ", "4GE6"]),
    ("PTPN11", "3ZM1", "3ZM1", ["3ZM1"]),
    ("PTPN12", "5HDE", "5HDE", ["5HDE", "5J8R"]),
    ("PTPN13", "1WCH", "1WCH", ["1WCH"]),
    ("PTPN14", "6IWD", "6IWD", ["6IWD"]),
    ("PTPN18", "4GFU", "4GFU", ["4GFU", "2OC3"]),
    ("PTPN21", "8GVV", "8GVV", ["8GVV", "8GWH"]),
    ("PTPN22", "3BRH", "3BRH", ["3BRH", "3OLR"]),
    ("PTPRN2", "2QEP", "2QEP", ["2QEP"]),
]

#: Kept out of the fourteen-protein run. PTPRN2 is catalytically dead: its
#: P-loop reads CSDGAGR where the family reads CSAGIGR. Excluding it changes the
#: conservation result not at all -- 58 unanimous columns either way -- which is
#: the point of running it both ways.
PSEUDOPHOSPHATASE = "PTPRN2"


def main():
    from WatCon.fetch import fetch_structures

    ids = sorted({pdb_id for _n, _r, _ref, group in MEMBERS for pdb_id in group})
    os.makedirs(RAW, exist_ok=True)
    print("Fetching %d entries from RCSB (native numbering)..." % len(ids))
    fetch_structures(ids, RAW)

    os.makedirs(STRUCTURES, exist_ok=True)
    rows_all, rows_without = [], []
    header = "# protein\tstructures directory\tConSurf grades\treference"
    rows_all.append(header)
    rows_without.append(header)

    for name, run, reference, group in MEMBERS:
        directory = os.path.join(STRUCTURES, name)
        os.makedirs(directory, exist_ok=True)
        for pdb_id in group:
            source = os.path.join(RAW, "%s.pdb" % pdb_id)
            if not os.path.isfile(source):
                sys.exit("missing %s -- did the fetch succeed?" % source)
            shutil.copyfile(source, os.path.join(directory, "%s.pdb" % pdb_id))

        grades = os.path.join(REPO, FIXTURES, "%s_A.grades.txt" % run)
        row = "\t".join([name, directory, grades, reference])
        rows_all.append(row)
        if name != PSEUDOPHOSPHATASE:
            rows_without.append(row)

    for path, rows, label in (
        (os.path.join(OUT, "members15.tsv"), rows_all, "fifteen proteins"),
        (os.path.join(OUT, "members14_no_pseudophosphatase.tsv"), rows_without,
         "fourteen, PTPRN2 excluded"),
    ):
        with open(path, "w", encoding="utf-8", newline="\n") as handle:
            handle.write("\n".join(rows) + "\n")
        print("wrote %s (%s)" % (os.path.relpath(path, REPO), label))

    print()
    print("Now run:")
    print("  watcon family --members %s \\"
          % os.path.relpath(os.path.join(OUT, "members15.tsv"), REPO))
    print("      --alignment WatCon/data/examples/ptp_family/ptp_family_alignment.pir")


if __name__ == "__main__":
    main()

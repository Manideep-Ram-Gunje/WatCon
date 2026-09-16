"""Fetch the barnase crystal structures for the conserved-water study.

ConSurf conservation is a property of the *sequence*, so the single ConSurf run
we hold for barnase (1BRS chain A, 150 sequences) describes every barnase
structure in the PDB.  That is what makes this study possible without any new
ConSurf data.

Selection is deliberately mechanical and recorded, so the set can be rebuilt:

  * UniProt P00648 (barnase, Bacillus amyloliquefaciens)
  * X-ray diffraction only -- water positions are the observable here, and
    NMR/predicted models do not have them
  * resolution <= RESOLUTION_LIMIT

Writes:
  data/entries.csv       every P00648 X-ray entry with its resolution
  data/selected.txt      the ones passing the cut, one PDB id per line
  structures/XXXX.pdb    the coordinate files (gitignored; rerun to restore)
"""

from __future__ import annotations

import csv
import json
import os
import sys
import time
import urllib.parse
import urllib.request

UNIPROT = "P00648"
RESOLUTION_LIMIT = 2.0

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
STRUCTURES = os.path.join(ROOT, "structures")

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
ENTRY = "https://data.rcsb.org/rest/v1/core/entry/"
DOWNLOAD = "https://files.rcsb.org/download/"


def search_entries(accession: str) -> list:
    """Every X-ray entry whose polymer maps to this UniProt accession."""
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": (
                            "rcsb_polymer_entity_container_identifiers"
                            ".reference_sequence_identifiers.database_accession"
                        ),
                        "operator": "exact_match",
                        "value": accession,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "exptl.method",
                        "operator": "exact_match",
                        "value": "X-RAY DIFFRACTION",
                    },
                },
            ],
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": 500},
            "results_content_type": ["experimental"],
        },
    }
    url = SEARCH + "?json=" + urllib.parse.quote(json.dumps(query))
    with urllib.request.urlopen(url, timeout=60) as handle:
        payload = json.load(handle)
    return sorted(hit["identifier"] for hit in payload.get("result_set", []))


def entry_metadata(pdb_id: str) -> dict:
    with urllib.request.urlopen(ENTRY + pdb_id, timeout=45) as handle:
        entry = json.load(handle)

    info = entry.get("rcsb_entry_info") or {}
    combined = info.get("resolution_combined") or []
    return {
        "pdb_id": pdb_id,
        "resolution": combined[0] if combined else None,
        "title": (entry.get("struct") or {}).get("title", "").strip(),
        "deposited_waters": info.get("deposited_solvent_atom_count"),
        "polymer_entities": info.get("polymer_entity_count_protein"),
    }


def download(pdb_id: str) -> bool:
    target = os.path.join(STRUCTURES, pdb_id + ".pdb")
    if os.path.exists(target) and os.path.getsize(target) > 0:
        return True
    try:
        with urllib.request.urlopen(DOWNLOAD + pdb_id + ".pdb", timeout=90) as handle:
            content = handle.read()
    except Exception as error:                       # noqa: BLE001
        print("  %s FAILED: %s" % (pdb_id, error))
        return False
    with open(target, "wb") as out:
        out.write(content)
    return True


def main() -> int:
    os.makedirs(DATA, exist_ok=True)
    os.makedirs(STRUCTURES, exist_ok=True)

    print("Searching RCSB for %s X-ray entries..." % UNIPROT)
    ids = search_entries(UNIPROT)
    print("  found %d" % len(ids))

    rows = []
    for i, pdb_id in enumerate(ids, 1):
        try:
            rows.append(entry_metadata(pdb_id))
        except Exception as error:                   # noqa: BLE001
            print("  %s metadata FAILED: %s" % (pdb_id, error))
            rows.append({"pdb_id": pdb_id, "resolution": None, "title": "",
                         "deposited_waters": None, "polymer_entities": None})
        if i % 10 == 0:
            print("  metadata %d/%d" % (i, len(ids)))
        time.sleep(0.05)

    path = os.path.join(DATA, "entries.csv")
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["pdb_id", "resolution", "deposited_waters",
                        "polymer_entities", "title"],
        )
        writer.writeheader()
        for row in sorted(rows, key=lambda r: (r["resolution"] is None,
                                               r["resolution"] or 0)):
            writer.writerow(row)
    print("Wrote %s (%d entries)" % (path, len(rows)))

    selected = sorted(
        r["pdb_id"] for r in rows
        if r["resolution"] is not None and r["resolution"] <= RESOLUTION_LIMIT
    )
    with open(os.path.join(DATA, "selected.txt"), "w") as handle:
        handle.write("\n".join(selected) + "\n")
    print("Selected %d entries at <= %.1f A" % (len(selected), RESOLUTION_LIMIT))

    ok = 0
    for i, pdb_id in enumerate(selected, 1):
        if download(pdb_id):
            ok += 1
        if i % 5 == 0:
            print("  downloaded %d/%d" % (i, len(selected)))
    print("Downloaded %d/%d structures into %s" % (ok, len(selected), STRUCTURES))
    return 0 if ok == len(selected) else 1


if __name__ == "__main__":
    sys.exit(main())

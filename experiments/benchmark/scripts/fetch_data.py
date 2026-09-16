"""Fetch the data for the WatCon comparison.

Two sources, both re-runnable and both cited:

  * The WatCon authors' own published dataset -- Zenodo 10.5281/zenodo.15213225,
    CC-BY-4.0, from Brownless, Harrison-Rawn & Kamerlin, JACS Au 2025. This is
    what makes the correctness comparison fair: their structures, their input
    files, their alignments, none of it chosen by us.

  * PTP1B crystal structures from RCSB, for the held-out benchmark. PTP1B is
    UniProt P18031; 1AAX chain A is the authors' own MSA reference.

Writes into data/, which is gitignored -- everything here is re-obtainable.
"""

from __future__ import annotations

import json
import os
import sys
import time
import urllib.parse
import urllib.request
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")

ZENODO = "https://zenodo.org/api/records/15213225/files/archive.zip/content"
UNIPROT = "P18031"          # PTP1B
RESOLUTION_LIMIT = 2.0
UA = {"User-Agent": "watcon-consurf-benchmark (research)"}


def fetch(url, target, tries=3):
    """Download with retries; Zenodo drops connections on long reads."""
    if os.path.exists(target) and os.path.getsize(target) > 0:
        print("  already have %s (%.1f MB)" % (os.path.basename(target),
                                               os.path.getsize(target) / 1e6))
        return True
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=300) as response:
                data = response.read()
            with open(target, "wb") as handle:
                handle.write(data)
            print("  fetched %s (%.1f MB)" % (os.path.basename(target), len(data) / 1e6))
            return True
        except Exception as error:                       # noqa: BLE001
            print("    attempt %d failed: %s" % (attempt + 1, type(error).__name__))
            time.sleep(3)
    return False


def get_zenodo():
    print("Zenodo 10.5281/zenodo.15213225 (WatCon authors' dataset, CC-BY-4.0)")
    archive = os.path.join(DATA, "watcon_zenodo.zip")
    if not fetch(ZENODO, archive):
        print("  FAILED -- Phase A cannot run without it")
        return False

    out = os.path.join(DATA, "zenodo")
    if not os.path.isdir(out):
        with zipfile.ZipFile(archive) as z:
            z.extractall(DATA)
        print("  unpacked to %s" % out)

    # The structure sets are themselves zipped inside the archive.
    for root, _dirs, files in os.walk(DATA):
        for name in files:
            if not name.endswith(".zip") or name == "watcon_zenodo.zip":
                continue
            path = os.path.join(root, name)
            target = os.path.join(root, name[:-4])
            if os.path.isdir(target):
                continue
            try:
                with zipfile.ZipFile(path) as z:
                    z.extractall(target)
                n = sum(len(f) for _, _, f in os.walk(target))
                print("  unpacked %-24s -> %d files" % (name, n))
            except zipfile.BadZipFile:
                print("  %s is not a readable zip" % name)
    return True


def search_rcsb(accession):
    query = {
        "query": {"type": "group", "logical_operator": "and", "nodes": [
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": ("rcsb_polymer_entity_container_identifiers"
                              ".reference_sequence_identifiers.database_accession"),
                "operator": "exact_match", "value": accession}},
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "exptl.method", "operator": "exact_match",
                "value": "X-RAY DIFFRACTION"}},
        ]},
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": 1000},
                            "results_content_type": ["experimental"]},
    }
    url = ("https://search.rcsb.org/rcsbsearch/v2/query?json="
           + urllib.parse.quote(json.dumps(query)))
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=90) as r:
        payload = json.load(r)
    return sorted(h["identifier"] for h in payload.get("result_set", []))


def resolution(pdb_id):
    url = "https://data.rcsb.org/rest/v1/core/entry/" + pdb_id
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=45) as r:
            entry = json.load(r)
        combined = (entry.get("rcsb_entry_info") or {}).get("resolution_combined") or []
        return combined[0] if combined else None
    except Exception:                                    # noqa: BLE001
        return None


def get_ptp1b():
    print()
    print("PTP1B (UniProt %s) from RCSB" % UNIPROT)
    ids = search_rcsb(UNIPROT)
    print("  %d X-ray entries" % len(ids))

    cache = os.path.join(DATA, "ptp1b_resolutions.json")
    known = {}
    if os.path.exists(cache):
        known = json.load(open(cache))

    for i, pdb_id in enumerate(ids, 1):
        if pdb_id not in known:
            known[pdb_id] = resolution(pdb_id)
            time.sleep(0.03)
        if i % 50 == 0:
            print("    resolutions %d/%d" % (i, len(ids)))
            json.dump(known, open(cache, "w"))
    json.dump(known, open(cache, "w"))

    selected = sorted(p for p in ids
                      if known.get(p) is not None and known[p] <= RESOLUTION_LIMIT)
    print("  %d at <= %.1f A" % (len(selected), RESOLUTION_LIMIT))

    with open(os.path.join(DATA, "ptp1b_selected.txt"), "w") as handle:
        handle.write("\n".join(selected) + "\n")

    out = os.path.join(DATA, "ptp1b_structures")
    os.makedirs(out, exist_ok=True)
    ok = 0
    for i, pdb_id in enumerate(selected, 1):
        target = os.path.join(out, pdb_id + ".pdb")
        if os.path.exists(target) and os.path.getsize(target) > 0:
            ok += 1
            continue
        try:
            url = "https://files.rcsb.org/download/%s.pdb" % pdb_id
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120) as r:
                content = r.read()
            with open(target, "wb") as handle:
                handle.write(content)
            ok += 1
        except Exception as error:                       # noqa: BLE001
            print("    %s failed: %s" % (pdb_id, type(error).__name__))
        if i % 25 == 0:
            print("    downloaded %d/%d" % (i, len(selected)))
    print("  %d/%d structures in %s" % (ok, len(selected), out))
    return ok > 0


def main() -> int:
    os.makedirs(DATA, exist_ok=True)
    zen = get_zenodo()
    ptp = get_ptp1b()
    print()
    print("Zenodo: %s   PTP1B: %s" % ("ok" if zen else "FAILED",
                                      "ok" if ptp else "FAILED"))
    return 0 if (zen and ptp) else 1


if __name__ == "__main__":
    sys.exit(main())

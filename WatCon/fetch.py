"""Download structures from RCSB by PDB id.

So that getting started is ``watcon fetch --ids 1AAX 2F71`` rather than a manual
trip to the website, and so the PyMOL plugin can offer a box you type ids into.

**PDB format is tried first, mmCIF second.** Not the other way round: the PDB
file is smaller and needs no conversion, but it does not always exist. Measured
against the real archive -- ``7GSA`` returns 404 for ``.pdb`` and 560 kB for
``.cif``, and 31 of the 287 PTP1B entries we selected behave the same way. A
fetcher that only asked for ``.pdb`` would silently lose a tenth of that set,
which is exactly what happened the first time.

Nothing here contacts ConSurf. Conservation data is yours to supply -- see
:mod:`WatCon.evolutionary`.
"""

from __future__ import annotations

import os
import urllib.error
import urllib.request

#: Where files come from. One entry per format, tried in order.
RCSB_URL = "https://files.rcsb.org/download/%s.%s"

#: Formats to try, best first.
FORMATS = ("pdb", "cif")

#: RCSB asks that automated clients identify themselves.
USER_AGENT = {"User-Agent": "WatCon-ConSurf (structural biology research tool)"}

DEFAULT_TIMEOUT = 60


class FetchError(RuntimeError):
    """Raised only when nothing at all could be downloaded."""


def _looks_like_pdb_id(text):
    """Four characters, starting with a digit -- the PDB's own rule.

    Extended 12-character ids (``pdb_0000xxxx``) are accepted too, since RCSB
    has begun issuing them and will rely on them once four characters run out.
    """
    text = text.strip()
    if len(text) == 4:
        return text[0].isdigit() and text.isalnum()
    return len(text) == 12 and text.lower().startswith("pdb_")


def fetch_structure(pdb_id, out_dir, timeout=DEFAULT_TIMEOUT, overwrite=False):
    """Download one entry. Returns the path written, or raises :class:`FetchError`.

    An existing non-empty file is kept unless ``overwrite``, so re-running is
    cheap and works offline for anything already fetched.
    """
    pdb_id = pdb_id.strip().upper()
    if not _looks_like_pdb_id(pdb_id):
        raise FetchError(
            "%r is not a PDB id. Ids are four characters beginning with a "
            "digit, such as 1AAX." % (pdb_id,)
        )

    os.makedirs(out_dir, exist_ok=True)

    if not overwrite:
        for suffix in FORMATS:
            existing = os.path.join(out_dir, "%s.%s" % (pdb_id, suffix))
            if os.path.exists(existing) and os.path.getsize(existing) > 0:
                return existing

    reasons = []
    for suffix in FORMATS:
        url = RCSB_URL % (pdb_id, suffix)
        try:
            request = urllib.request.Request(url, headers=USER_AGENT)
            with urllib.request.urlopen(request, timeout=timeout) as response:
                content = response.read()
        except urllib.error.HTTPError as error:
            # 404 on .pdb is the normal case for a large or recent entry, not a
            # failure -- it is why .cif is tried next.
            reasons.append("%s: HTTP %s" % (suffix, error.code))
            continue
        except Exception as error:          # noqa: BLE001 - urllib raises broadly
            reasons.append("%s: %s" % (suffix, type(error).__name__))
            continue

        if not content:
            reasons.append("%s: empty response" % suffix)
            continue

        target = os.path.join(out_dir, "%s.%s" % (pdb_id, suffix))
        with open(target, "wb") as handle:
            handle.write(content)
        return target

    raise FetchError(
        "Could not download %s from RCSB (%s). Check the id and your network "
        "connection; obsolete entries are served from "
        "https://files.rcsb.org/download/ only while they remain current."
        % (pdb_id, "; ".join(reasons))
    )


def fetch_structures(ids, out_dir, timeout=DEFAULT_TIMEOUT, overwrite=False,
                     verbose=True):
    """Download several entries. Returns ``(paths, failures)``.

    One bad id does not stop the rest -- ``failures`` is ``[(id, reason), ...]``
    and the caller decides what that means. Raises :class:`FetchError` only if
    **nothing** was downloaded, because a run that produces an empty directory
    should not proceed to report "no structure files found" as if the user had
    pointed at the wrong place.
    """
    ids = [str(i).strip() for i in ids if str(i).strip()]
    if not ids:
        raise FetchError("No PDB ids given.")

    paths, failures = [], []
    for index, pdb_id in enumerate(ids, 1):
        try:
            path = fetch_structure(pdb_id, out_dir, timeout=timeout,
                                   overwrite=overwrite)
        except FetchError as error:
            failures.append((pdb_id, str(error)))
            if verbose:
                print("  %-6s FAILED" % pdb_id.upper())
            continue
        paths.append(path)
        if verbose:
            print("  %-6s %s (%.0f kB)"
                  % (pdb_id.upper(), os.path.splitext(path)[1].lstrip("."),
                     os.path.getsize(path) / 1024.0))
        if verbose and index % 25 == 0:
            print("  ... %d/%d" % (index, len(ids)))

    if not paths:
        raise FetchError(
            "None of the %d requested structures could be downloaded. First "
            "reason: %s" % (len(ids), failures[0][1] if failures else "unknown")
        )

    if verbose:
        print("Fetched %d/%d structures into %s" % (len(paths), len(ids), out_dir))
        if failures:
            print("Could not fetch: %s"
                  % ", ".join(pdb_id for pdb_id, _ in failures))
    return paths, failures

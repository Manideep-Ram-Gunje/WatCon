"""Reading the structure formats the PDB actually serves today.

WatCon's preparation step parses ``ATOM``/``HETATM`` text directly rather than
going through MDAnalysis, which is fast and dependency-free but means it only
understands PDB format. That is no longer enough:

    $ watcon prepare --input-dir raw/ ...
    ValueError: 'CIF' isn't a valid topology format

RCSB has stopped issuing PDB-format files for large and recent entries. This is
not hypothetical -- fetching the PTP1B benchmark set, **31 of 287 selected
entries had no PDB file at all**, among them every ``7GS*``/``7GT*`` fragment
structure. ``7GSA`` returns 404 for ``.pdb`` and 560 kB for ``.cif``.

MDAnalysis 2.10 does not read mmCIF either (it knows MMTF, which RCSB has since
retired), so the conversion has to happen before either reader sees the file.
``gemmi`` does it in-process, and is the same library the PDB's own tooling uses.

Everything here converts *to* PDB text and then gets out of the way, so the rest
of WatCon is unchanged and cannot acquire a second, subtly different parser.
"""

from __future__ import annotations

import gzip
import os
import shutil

#: Extensions that must be converted before WatCon can read them.
CONVERTIBLE_SUFFIXES = (".cif", ".mmcif")

#: Extensions WatCon reads directly, once any gzip wrapper is removed.
NATIVE_SUFFIXES = (".pdb", ".ent")


class ConversionError(ValueError):
    """Raised when a structure file cannot be turned into readable PDB text."""


def _strip_gz(name):
    """``('1abc.cif', True)`` for ``1abc.cif.gz``; unchanged and False otherwise."""
    if name.lower().endswith(".gz"):
        return name[:-3], True
    return name, False


def needs_conversion(name):
    """True if this filename is a format WatCon cannot parse as-is."""
    stem, _ = _strip_gz(os.path.basename(name))
    return stem.lower().endswith(CONVERTIBLE_SUFFIXES)


def is_readable_structure(name):
    """True if WatCon can use this file, with or without conversion."""
    stem, _ = _strip_gz(os.path.basename(name))
    return stem.lower().endswith(NATIVE_SUFFIXES + CONVERTIBLE_SUFFIXES)


def _require_gemmi():
    try:
        import gemmi
    except ImportError as error:            # pragma: no cover - trivial branch
        raise ConversionError(
            "Reading mmCIF needs the 'gemmi' package, which is not installed. "
            "Install it with:  pip install gemmi\n"
            "(WatCon lists gemmi as a dependency, so this normally means the "
            "environment is only partly installed.)"
        ) from error
    return gemmi


def convert_to_pdb(path, out_path):
    """Write ``path`` out as PDB text at ``out_path``. Returns ``out_path``.

    Handles mmCIF and gzipped inputs. A file that is already plain PDB is copied
    rather than round-tripped through gemmi -- reformatting a file WatCon can
    already read would only risk changing it.

    Raises
    ------
    ConversionError
        With the offending filename in the message. The common real cause is an
        entry too large for PDB format at all (more than 62 chains or 99,999
        atoms), which is precisely why such entries are mmCIF-only; the message
        says so rather than reporting gemmi's internal complaint.
    """
    name = os.path.basename(path)
    stem, gzipped = _strip_gz(name)

    if not needs_conversion(name):
        if gzipped:
            with gzip.open(path, "rb") as source, open(out_path, "wb") as target:
                shutil.copyfileobj(source, target)
        else:
            shutil.copyfile(path, out_path)
        return out_path

    gemmi = _require_gemmi()

    # gemmi reads gzip itself, by extension.
    try:
        structure = gemmi.read_structure(path)
    except Exception as error:              # noqa: BLE001 - gemmi raises broadly
        raise ConversionError(
            "Could not read %s as mmCIF: %s: %s"
            % (name, type(error).__name__, error)
        ) from error

    # setup_entities fills in the entity/subchain bookkeeping that a PDB record
    # needs and an mmCIF file may not carry explicitly.
    structure.setup_entities()

    try:
        structure.write_pdb(out_path)
    except Exception as error:              # noqa: BLE001
        raise ConversionError(
            "%s cannot be expressed in PDB format (%s). Entries exceeding "
            "62 chains or 99,999 atoms are mmCIF-only by design; WatCon cannot "
            "prepare them. Use a smaller entry, or extract the chain of "
            "interest first." % (name, type(error).__name__)
        ) from error
    return out_path


def pdb_name_for(name):
    """The filename ``name`` will have once converted: ``7GSA.cif.gz`` -> ``7GSA.pdb``.

    Only the format suffix is replaced. A name carrying its own dots keeps them,
    so ``1A2P.run1.cif`` becomes ``1A2P.run1.pdb`` and still matches its ConSurf
    file -- the same truncation bug fixed in the network builder.
    """
    stem, _ = _strip_gz(os.path.basename(name))
    root, _ = os.path.splitext(stem)
    return root + ".pdb"


def as_pdb_directory(input_dir, work_dir, verbose=True):
    """Give back a directory of PDB files, converting whatever needs it.

    Returns ``(directory, converted, failed)``:

    * ``directory`` -- ``input_dir`` itself when nothing needed converting, so
      the common case copies nothing;
    * ``converted`` -- ``[(original name, new name), ...]``;
    * ``failed`` -- ``[(name, reason), ...]``, one entry per file that could not
      be converted. Conversion failures do **not** abort the run: one unreadable
      entry in a folder of three hundred should not cost you the other 299.
      The caller reports them.
    """
    from WatCon.residue_index import require_structure_files

    names, _skipped = require_structure_files(input_dir)
    to_convert = [n for n in names if needs_conversion(n) or n.lower().endswith(".gz")]
    if not to_convert:
        return input_dir, [], []

    os.makedirs(work_dir, exist_ok=True)
    converted, failed = [], []
    for name in names:
        target = os.path.join(work_dir, pdb_name_for(name))
        try:
            convert_to_pdb(os.path.join(input_dir, name), target)
        except ConversionError as error:
            failed.append((name, str(error)))
            continue
        if name != os.path.basename(target):
            converted.append((name, os.path.basename(target)))

    if verbose and converted:
        print("Converted %d file(s) to PDB format: %s"
              % (len(converted), ", ".join(a for a, _ in converted[:6])
                 + (", ..." if len(converted) > 6 else "")))
    if verbose and failed:
        print("Could not convert %d file(s):" % len(failed))
        for name, reason in failed[:6]:
            print("  %s: %s" % (name, reason.splitlines()[0]))

    if not os.listdir(work_dir):
        raise ConversionError(
            "Every structure in %s failed to convert. First reason: %s"
            % (input_dir, failed[0][1] if failed else "unknown")
        )
    return work_dir, converted, failed

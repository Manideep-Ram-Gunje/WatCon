"""Build the CI fixtures for WatCon's trajectory and directed paths.

Three files, all real coordinates, all small enough for CI:

``ptp1b_ensemble_ragged.pdb``
    Three superposed PTP1B crystal structures around the catalytic centre.
    The models hold different numbers of atoms, because each entry resolves a
    different set of waters. This is what a crystal ensemble actually looks
    like, and it is the input the dynamic path must refuse.

``ptp1b_ensemble_constant.pdb``
    The same three, made readable as a trajectory: the protein atoms all three
    resolve, plus the N water oxygens nearest the catalytic centre in each,
    N being the smallest any of them offers. Every coordinate is a real
    crystallographic position; the fixed water count is a stated selection, so
    this fixture tests the code path and makes no claim about water occupancy.

``md_active_site_h.pdb``
    The active site of the Zenodo MD system, which is where the hydrogens are.

**Needs data this repository does not carry**: the 253 prepared PTP1B
structures and the Zenodo MD system. It is committed so the fixtures' provenance
is inspectable -- so a reader can see exactly how they were cut -- not so that
they can be rebuilt from a clean clone. See the README beside it.

    python experiments/benchmark/scripts/build_ensemble_fixtures.py
"""

from __future__ import annotations

import os
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import MDAnalysis as mda

HERE = os.path.dirname(os.path.abspath(__file__))
BENCHMARK = os.path.dirname(HERE)
RESEARCH = os.path.dirname(os.path.dirname(BENCHMARK))

PREPARED = os.path.join(BENCHMARK, "prepared", "ptp1b")
GRO = os.path.join(BENCHMARK, "data", "zenodo", "MD", "starting_structures",
                   "WT_PTP1B_Apo_Closed.gro")
DEST = os.path.join(RESEARCH, "tools", "WatCon_ConSurf",
                    "WatCon", "data", "examples", "ptp1b_ensemble")

#: Three of the 253 prepared structures. Superposed already, so the models share
#: a frame; 1AAX is the C215S trap, which is why the centre falls back to CA.
MEMBERS = ("1AAX", "1C83", "1C88")
RADIUS = 12.0
MD_RADIUS = 10.0


def _atom_lines(path):
    return [line for line in open(path, encoding="utf-8", errors="replace")
            if line.startswith(("ATOM", "HETATM"))]


def _coords(lines):
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in lines])


def _is_water(line):
    return line[17:20].strip() in ("HOH", "WAT")


def _key(line):
    return line[21], line[22:27].strip(), line[12:16].strip()


def trimmed_models():
    """Each member trimmed to RADIUS around the catalytic centre of the first."""
    centre = None
    out = []
    for pdb_id in MEMBERS:
        lines = _atom_lines(os.path.join(PREPARED, "%s.pdb" % pdb_id))
        coords = _coords(lines)
        if centre is None:
            found = [i for i, l in enumerate(lines)
                     if l[22:26].strip() == "215" and l[12:16].strip() in ("SG", "OG")]
            if not found:
                found = [i for i, l in enumerate(lines)
                         if l[22:26].strip() == "215" and l[12:16].strip() == "CA"]
            centre = coords[found[0]]
        keep = np.linalg.norm(coords - centre, axis=1) <= RADIUS
        out.append((pdb_id, [l for l, k in zip(lines, keep) if k]))
    return out, centre


def write_ragged(models, path):
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(
            "REMARK   Three superposed PTP1B crystal structures, each trimmed to %g A\n"
            "REMARK   around the catalytic centre. The models hold different numbers of\n"
            "REMARK   atoms because each entry resolves a different set of waters. This is\n"
            "REMARK   what an ensemble of crystal structures looks like, and it is why one\n"
            "REMARK   cannot be read as a trajectory -- use the static path for such a set.\n"
            % RADIUS)
        for i, (pdb_id, _lines) in enumerate(models, 1):
            handle.write("REMARK   MODEL %d is RCSB %s\n" % (i, pdb_id))
        for i, (_pdb_id, lines) in enumerate(models, 1):
            handle.write("MODEL     %4d\n" % i)
            handle.writelines(lines)
            handle.write("ENDMDL\n")
        handle.write("END\n")
    return [len(lines) for _pdb_id, lines in models]


def write_constant(models, centre, path):
    common = None
    for _pdb_id, lines in models:
        keys = {_key(l) for l in lines if not _is_water(l)}
        common = keys if common is None else (common & keys)

    waters = {}
    for pdb_id, lines in models:
        found = [l for l in lines if _is_water(l) and l[12:16].strip().startswith("O")]
        found.sort(key=lambda l: float(np.linalg.norm(_coords([l])[0] - centre)))
        waters[pdb_id] = found
    n_waters = min(len(v) for v in waters.values())

    counts = []
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(
            "REMARK   The same three structures as the ragged ensemble, made readable as a\n"
            "REMARK   trajectory: the protein atoms all three resolve, plus the %d water\n"
            "REMARK   oxygens nearest the catalytic centre in each. Every coordinate is a\n"
            "REMARK   real crystallographic position; the fixed water count is a stated\n"
            "REMARK   selection, not a property of the structures, so this file tests the\n"
            "REMARK   code path and makes no claim about water occupancy.\n" % n_waters)
        for i, (pdb_id, _lines) in enumerate(models, 1):
            handle.write("REMARK   MODEL %d is RCSB %s\n" % (i, pdb_id))
        for i, (pdb_id, lines) in enumerate(models, 1):
            kept, seen = [], set()
            for line in lines:
                key = _key(line)
                if key in common and key not in seen:
                    seen.add(key)
                    kept.append(line)
            kept.sort(key=lambda l: (l[21], int(l[22:26]), l[12:16]))
            kept.extend(waters[pdb_id][:n_waters])
            counts.append(len(kept))
            handle.write("MODEL     %4d\n" % i)
            handle.writelines(kept)
            handle.write("ENDMDL\n")
        handle.write("END\n")
    return counts, n_waters


def write_md_site(path):
    universe = mda.Universe(GRO)
    near = universe.select_atoms(
        "(around %g (resid 214 and name SG)) or (resid 214)" % MD_RADIUS)
    site = universe.select_atoms("same residue as group s", s=near)
    site.write(path)
    body = open(path, encoding="utf-8").read()
    header = (
        "REMARK   Active site of the Zenodo MD system WT_PTP1B_Apo_Closed.gro: every\n"
        "REMARK   residue with an atom within %g A of the nucleophile. Amber ff14SB, so\n"
        "REMARK   hydrogens are present and waters are named WAT -- this is the only real\n"
        "REMARK   hydrogen-bearing PTP1B system available, and the directed path needs it.\n"
        "REMARK   Numbering is the MD system's own: resid 214 here is PTP1B Cys215.\n"
        % MD_RADIUS)
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(header + body)
    return len(site), len(site.select_atoms("resname WAT").residues), \
        len(site.select_atoms("name H*"))


def main():
    os.makedirs(DEST, exist_ok=True)
    models, centre = trimmed_models()

    ragged = os.path.join(DEST, "ptp1b_ensemble_ragged.pdb")
    counts = write_ragged(models, ragged)
    print("ragged    %6d B  model atom counts %s" % (os.path.getsize(ragged), counts))

    constant = os.path.join(DEST, "ptp1b_ensemble_constant.pdb")
    counts, n_waters = write_constant(models, centre, constant)
    print("constant  %6d B  model atom counts %s  (%d waters each)"
          % (os.path.getsize(constant), counts, n_waters))

    md = os.path.join(DEST, "md_active_site_h.pdb")
    atoms, waters, hydrogens = write_md_site(md)
    print("md site   %6d B  %d atoms, %d waters, %d hydrogens"
          % (os.path.getsize(md), atoms, waters, hydrogens))


if __name__ == "__main__":
    main()

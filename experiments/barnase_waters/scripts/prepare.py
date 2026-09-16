"""Put every barnase structure into one coordinate frame, with its waters.

Four steps, each of which can reject a structure loudly rather than quietly
degrading the dataset:

  1. Find the barnase chain.  Several entries are barnase-barstar complexes or
     fusions, so the chain is identified by sequence identity to the reference
     barnase, not by assuming chain 'A'.
  2. Keep that chain plus the waters nearest to it, and relabel the chain 'A'.
     Relabelling means one ConSurf file serves every structure with no
     chain_map: the ConSurf run is on barnase, and after this so is chain A
     everywhere.
  3. Superpose onto the reference with Kabsch (WatCon.superpose), carrying the
     waters along.
  4. Report RMSD and rejections.

Writes prepared/XXXX.pdb and results/preparation.csv.
"""

from __future__ import annotations

import csv
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
# ROOT is <repo>/experiments/barnase_waters, so two levels up is the workspace.
WATCON = os.path.join(
    os.path.dirname(os.path.dirname(ROOT)), "tools", "WatCon_ConSurf",
)
sys.path.insert(0, WATCON)

from WatCon.residue_index import (                       # noqa: E402
    build_fasta_sequence,
    residues_from_pdb_file,
)
from WatCon.superpose import kabsch, rmsd                # noqa: E402

STRUCTURES = os.path.join(ROOT, "structures")
PREPARED = os.path.join(ROOT, "prepared")
RESULTS = os.path.join(ROOT, "results")

#: Highest-resolution entry whose barnase is a plain single chain.  2C4B (1.30 A)
#: is a fusion protein and 6PQK (1.20 A) a disulfide-engineered complex, so the
#: reference is the best ordinary barnase: the 1.5 A wild-type structure.
REFERENCE = "1A2P"

#: Minimum identity to the reference barnase sequence for a chain to count.
#: Generous, because many entries are point mutants -- but far above what an
#: unrelated chain such as barstar (89 residues, no homology) could reach.
MIN_IDENTITY = 0.80

#: Waters within this distance of any barnase atom are kept with the chain.
WATER_CUTOFF = 5.0


def read_atoms(path):
    """All ATOM/HETATM lines, as (line, chain, resname, xyz)."""
    rows = []
    with open(path, "r", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")) or len(line) < 54:
                continue
            # Only the first model of a multi-model file.
            try:
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            rows.append((line.rstrip("\n"), line[21], line[17:20].strip(), xyz))
    return rows


def identity(residues, reference_by_resid) -> float:
    """Fraction of shared residue NUMBERS carrying the same amino acid.

    Keyed on residue number rather than position in the string, because residue
    number is what the ConSurf join uses.  Comparing strings positionally is
    wrong twice over, and both errors are present in this dataset:

      * 1BSA starts at residue 4, 1B2S at 1, 1RNB at 2, against the reference's
        3.  Positionally their sequences look unrelated (identity 0.03-0.05) and
        all three were wrongly rejected -- yet by residue number they agree
        perfectly.
      * 2F56 and 2F5M start at residue 1 with the residue everyone else numbers
        3.  Positionally they score a perfect 1.00 and were accepted, but their
        numbering is shifted by two, so every ConSurf score would land on the
        wrong residue.  They superposed at 6.05 A, which is what first showed it.

    Comparing the way the join compares gets both cases right.
    """
    shared = [r for r in residues if r.resid in reference_by_resid]
    if not shared:
        return 0.0
    same = sum(1 for r in shared if r.one_letter == reference_by_resid[r.resid])
    # Denominator is the reference length, so a chain matching only a fragment
    # cannot score highly just by agreeing on the few residues it has.
    return same / len(reference_by_resid)


def detect_offset(residues, reference_by_resid, span=6):
    """Best integer renumbering, as (offset, identity) -- diagnostic only.

    Reported so a rejection says *why*, rather than leaving a real barnase
    structure looking like an unrelated chain.  Never applied: silently
    renumbering a deposited structure is the same class of helpful guess this
    project has spent its time removing.
    """
    best = (0, identity(residues, reference_by_resid))
    for offset in range(-span, span + 1):
        if offset == 0:
            continue
        shifted = [type(r)(r.chain, r.resid + offset, r.icode, r.resname,
                           r.one_letter) for r in residues]
        score = identity(shifted, reference_by_resid)
        if score > best[1]:
            best = (offset, score)
    return best


def find_barnase_chain(path: str, reference_by_resid):
    """(chain, identity, offset_note) for the chain most like barnase."""
    best_chain, best_score = None, 0.0
    residues = residues_from_pdb_file(path)
    for chain in sorted({r.chain for r in residues}):
        score = identity([r for r in residues if r.chain == chain],
                         reference_by_resid)
        if score > best_score:
            best_chain, best_score = chain, score

    if best_score >= MIN_IDENTITY:
        return best_chain, best_score, ""

    # Explain the rejection: is this barnase with shifted numbering?
    note = ""
    for chain in sorted({r.chain for r in residues}):
        offset, score = detect_offset(
            [r for r in residues if r.chain == chain], reference_by_resid
        )
        if offset and score >= MIN_IDENTITY:
            note = ("chain %s matches at %.2f if renumbered by %+d"
                    % (chain, score, offset))
            break
    return None, best_score, note


def extract(path: str, chain: str):
    """(protein_lines, water_lines) for one chain, waters within the cutoff."""
    rows = read_atoms(path)
    protein, waters = [], []
    for line, chain_id, resname, xyz in rows:
        if resname == "HOH" or resname == "WAT":
            waters.append((line, xyz))
        elif chain_id == chain:
            protein.append((line, xyz))

    if not protein:
        return [], []

    coords = np.array([xyz for _, xyz in protein])
    kept = []
    for line, xyz in waters:
        distance = np.min(np.linalg.norm(coords - np.array(xyz), axis=1))
        if distance <= WATER_CUTOFF:
            kept.append((line, xyz))
    return protein, kept


def relabel_and_transform(lines, rotation, translation):
    """Rewrite chain to 'A' and apply the rigid transform, in PDB columns."""
    out = []
    for line, xyz in lines:
        moved = np.asarray(xyz) @ rotation.T + translation
        padded = line.ljust(80)
        out.append(
            padded[:21] + "A" + padded[22:30]
            + "%8.3f%8.3f%8.3f" % (moved[0], moved[1], moved[2])
            + padded[54:80]
        )
    return out


def ca_map(lines):
    """(resid, icode) -> CA coordinate, from already-extracted lines."""
    result = {}
    for line, xyz in lines:
        if line[12:16].strip() != "CA":
            continue
        try:
            resid = int(line[22:26].strip())
        except ValueError:
            continue
        key = (resid, line[26:27].strip() or None)
        result.setdefault(key, np.array(xyz))
    return result


def main() -> int:
    # Clear rather than merge.  A structure rejected on this run must not be
    # left behind from a previous one: that is exactly what happened when the
    # sequence comparison was fixed, and a stale 2F56.pdb survived into the
    # study, where the ConSurf identity check caught it independently.
    if os.path.isdir(PREPARED):
        shutil.rmtree(PREPARED)
    os.makedirs(PREPARED)
    os.makedirs(RESULTS, exist_ok=True)

    ids = sorted(
        f[:-4] for f in os.listdir(STRUCTURES) if f.lower().endswith(".pdb")
    )
    print("Preparing %d structures (reference %s)" % (len(ids), REFERENCE))

    ref_path = os.path.join(STRUCTURES, REFERENCE + ".pdb")
    ref_residues = residues_from_pdb_file(ref_path)
    ref_chains = sorted({r.chain for r in ref_residues})
    ref_chain = max(
        ref_chains,
        key=lambda c: len([r for r in ref_residues if r.chain == c]),
    )
    reference_by_resid = {
        r.resid: r.one_letter for r in ref_residues if r.chain == ref_chain
    }
    print("  reference chain %s, %d residues (resid %d..%d)"
          % (ref_chain, len(reference_by_resid),
             min(reference_by_resid), max(reference_by_resid)))

    ref_protein, ref_waters = extract(ref_path, ref_chain)
    ref_cas = ca_map(ref_protein)

    rows = []
    for pdb_id in ids:
        path = os.path.join(STRUCTURES, pdb_id + ".pdb")
        chain, score, note = find_barnase_chain(path, reference_by_resid)
        row = {
            "pdb_id": pdb_id, "chain": chain or "", "seq_identity": round(score, 3),
            "n_protein_atoms": 0, "n_waters": 0, "n_ca_matched": 0,
            "rmsd": "", "status": "",
        }

        if chain is None:
            row["status"] = ("rejected: identity %.2f by residue number%s"
                             % (score, "; " + note if note else ""))
            rows.append(row)
            print("  %s  REJECTED (identity %.2f)%s"
                  % (pdb_id, score, "  [" + note + "]" if note else ""))
            continue

        protein, waters = extract(path, chain)
        cas = ca_map(protein)
        shared = sorted(set(ref_cas) & set(cas))

        if len(shared) < 30:
            row["status"] = "rejected: only %d CA shared with reference" % len(shared)
            rows.append(row)
            print("  %s  REJECTED (%d CA shared)" % (pdb_id, len(shared)))
            continue

        mobile = np.array([cas[k] for k in shared])
        target = np.array([ref_cas[k] for k in shared])
        rotation, translation = kabsch(mobile, target)
        fit = rmsd(mobile @ rotation.T + translation, target)

        lines = relabel_and_transform(protein, rotation, translation)
        lines += relabel_and_transform(waters, rotation, translation)
        with open(os.path.join(PREPARED, pdb_id + ".pdb"), "w") as handle:
            handle.write("\n".join(lines) + "\nEND\n")

        row.update(
            n_protein_atoms=len(protein), n_waters=len(waters),
            n_ca_matched=len(shared), rmsd=round(fit, 3), status="ok",
        )
        rows.append(row)
        print("  %s  chain %s  id %.2f  %3d waters  RMSD %.2f A over %d CA"
              % (pdb_id, chain, score, len(waters), fit, len(shared)))

    path = os.path.join(RESULTS, "preparation.csv")
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    ok = [r for r in rows if r["status"] == "ok"]
    print()
    print("Prepared %d/%d structures, %d waters total"
          % (len(ok), len(rows), sum(r["n_waters"] for r in ok)))
    print("Wrote %s" % path)
    return 0


if __name__ == "__main__":
    sys.exit(main())

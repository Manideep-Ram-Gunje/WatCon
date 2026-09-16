"""A PyMOL session for a protein family: the structures in one frame, and the
water sites they share.

The single-protein equivalent is :mod:`WatCon.scene`, and the same rules apply
here, each learned the hard way:

* the ``.pml`` **loads what it colours**, so opening it shows something;
* no emitted line contains a semicolon -- PyMOL splits on ``;`` even inside a
  ``#`` comment, and ran a prose comment as Python;
* a group carries its own enabled flag, so disabling its members is not enough;
* paths are written with forward slashes, which PyMOL reads on every platform.

What the picture says
---------------------
A site is drawn larger the more proteins of the family hold a water there, and
coloured by whether the residues lining it are conserved in **every** member --
the one family-level claim that never compares separately normalised scores.
Occupancy and conservation stay in separate columns of the PDB (occupancy and
B-factor), exactly as in the single-protein session.
"""

from __future__ import annotations

import os
from typing import List, Optional

from .family_sites import FamilySites

__all__ = ["FAMILY_COLOURS", "write_family_session"]

#: Distinguishable, colour-blind-safe enough, and stable per position so two
#: runs of the same family colour the same protein the same way.
FAMILY_COLOURS = [
    (0.12, 0.47, 0.71), (0.89, 0.47, 0.13), (0.17, 0.63, 0.35),
    (0.84, 0.24, 0.31), (0.58, 0.40, 0.74), (0.55, 0.34, 0.29),
    (0.89, 0.47, 0.76), (0.50, 0.50, 0.50), (0.74, 0.74, 0.13),
    (0.09, 0.75, 0.81),
]

#: Sphere radius for a site held by one protein, and by every protein.
SITE_RADIUS_RANGE = (0.25, 0.75)


def _posix(path) -> str:
    return str(path).replace("\\", "/")


def _write_sites_pdb(path: str, sites, n_proteins: int) -> None:
    """Sites as waters: occupancy = proteins holding one, B-factor = best grade.

    B-factor 0 means no ConSurf-scored residue lines the site, which is not the
    same as low conservation -- stated in the REMARKs, as elsewhere.
    """
    with open(path, "w", newline="\n") as handle:
        handle.write("REMARK   WatCon + ConSurf family water sites.\n")
        handle.write("REMARK   occupancy = fraction of the family's proteins holding\n")
        handle.write("REMARK               a water here (0-1).\n")
        handle.write("REMARK   B-factor  = highest ConSurf grade among the residues\n")
        handle.write("REMARK               lining it, in any protein. 0 means NO DATA.\n")
        for serial, site in enumerate(sites, 1):
            grade = site.max_grade or 0
            fraction = site.n_proteins_occupied / n_proteins if n_proteins else 0.0
            handle.write(
                "ATOM  %5d  O   HOH S%4d    %8.3f%8.3f%8.3f%6.2f%6.2f           O\n"
                % (serial, min(site.site_id, 9999), site.centre[0], site.centre[1],
                   site.centre[2], fraction, float(grade)))
        handle.write("END\n")


def write_family_session(
    sites: FamilySites,
    out_dir: Optional[str] = None,
    session_name: str = "watcon_family.pml",
) -> str:
    """Write the family session. Returns the path to the ``.pml``.

    Needs ``sites`` to carry a ``superposed_dir``, which
    :func:`WatCon.family_sites.build_family_sites` fills in -- the structures
    must be in one frame before anything is drawn together.
    """
    out_dir = out_dir or sites.superposed_dir
    if not out_dir:
        raise ValueError("no output directory, and the sites carry no superposed_dir")
    os.makedirs(out_dir, exist_ok=True)

    proteins: List[str] = []
    for fit in sites.fits:
        if fit.protein not in proteins:
            proteins.append(fit.protein)

    all_path = os.path.join(out_dir, "family_sites_all.pdb")
    shared_path = os.path.join(out_dir, "family_sites_shared.pdb")
    everywhere = [s for s in sites.sites if s.n_proteins_occupied == len(proteins)]
    _write_sites_pdb(all_path, sites.sites, len(proteins))
    _write_sites_pdb(shared_path, everywhere, len(proteins))

    low, high = SITE_RADIUS_RANGE
    lines = [
        "# WatCon + ConSurf -- a protein family in one frame, and the water sites",
        "# it shares. Open with:   pymol %s" % session_name,
        "#",
        "# Frame reference: %s" % sites.reference,
        "# %d proteins, %d structures, %d sites, %d held by every protein"
        % (len(proteins), len(sites.fits), len(sites.sites), len(everywhere)),
        "",
        "bg white",
    ]

    for index, protein in enumerate(proteins):
        r, g, b = FAMILY_COLOURS[index % len(FAMILY_COLOURS)]
        lines.append("set_color family_%d, [%.3f, %.3f, %.3f]" % (index, r, g, b))
    lines.append("")

    for fit in sorted(sites.fits, key=lambda f: (proteins.index(f.protein), f.pdb_id)):
        path = os.path.join(sites.superposed_dir or out_dir, fit.pdb_id + ".pdb")
        lines.append("load %s, %s" % (_posix(path), fit.pdb_id))
    lines += [
        "load %s, sites" % _posix(all_path),
        "load %s, sites_shared" % _posix(shared_path),
        "",
        "hide everything",
    ]

    for fit in sites.fits:
        index = proteins.index(fit.protein)
        lines.append("show cartoon, %s" % fit.pdb_id)
        lines.append("color family_%d, %s" % (index, fit.pdb_id))
    lines += [
        "set cartoon_transparency, 0.6",
        "",
        "# Every site, sized by how many proteins hold a water there. Off by",
        "# default because there are %d of them." % len(sites.sites),
        "show spheres, sites",
        "color grey70, sites",
        "alter sites, vdw=%.2f+%.2f*q" % (low, high - low),
        "set sphere_scale, 1.0, sites",
        "rebuild sites",
        "disable sites",
        "",
        "# The sites every protein of the family holds. That is the picture.",
        "show spheres, sites_shared",
        "color firebrick, sites_shared",
        "alter sites_shared, vdw=%.2f" % high,
        "set sphere_scale, 1.0, sites_shared",
        "rebuild sites_shared",
        "",
        "orient %s" % sites.reference,
        "zoom %s, 4" % sites.reference,
        "",
        "print('')",
        "print('  WatCon + ConSurf, family view')",
    ]
    for index, protein in enumerate(proteins):
        members = ", ".join(f.pdb_id for f in sites.fits if f.protein == protein)
        lines.append("print('  %-10s %s')" % (protein, members))
    lines += [
        "print('')",
        "print('  firebrick    %d site(s) where every protein holds a water')" % len(everywhere),
        "print('  enable sites all %d sites, larger where more proteins hold one')"
        % len(sites.sites),
        "print('  a site B-factor of 0 means NO conservation data, not low conservation')",
        "print('')",
    ]

    pml = os.path.join(out_dir, session_name)
    with open(pml, "w", newline="\n") as handle:
        handle.write("\n".join(lines) + "\n")
    return pml

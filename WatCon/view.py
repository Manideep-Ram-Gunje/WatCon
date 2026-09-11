"""``watcon view`` -- build a PyMOL session you can actually open.

The projection functions write correct data but not a usable session. A ``.pml``
of colour commands with no ``load`` opens an empty window and colours nothing --
every command succeeds, so there is no error to notice:

    $ pymol EVO_RESIDUES.pml
    objects loaded: []
    atoms visible: 0

And a cluster PDB carries ConSurf grades in its B-factor column but arrives grey.
Both outputs are right, and neither shows you anything.

This assembles one self-contained script:

    protein           cartoon, coloured on ConSurf's own 1-9 scale
    sites_conserved   sites lined by a highly conserved residue -- shown
    sites             every occupied site, coloured by grade -- available

so that ``pymol watcon_view.pml`` opens on the picture the tool exists to make.
"""

from __future__ import annotations

import os
import shutil
import tempfile
import warnings
from pathlib import Path

#: Grade at or above which a site is called highly conserved. 8 because that is
#: the most reproducible end of a ConSurf run -- comparing two runs of the same
#: sequence, grade >= 8 was far more stable than the variable end.
HIGHLY_CONSERVED = 8


def build_session(prepared, consurf, out_dir="watcon_view", reference=None,
                  site_radius=1.5, min_cluster_samples=2, verbose=True):
    """Build the conserved-water picture and the script that shows it.

    Parameters
    ----------
    prepared : str
        Directory of prepared structures -- output of ``watcon prepare``.
    consurf : str
        A ConSurf ``*_consurf_grades.txt`` file. One run describes every
        structure of the same protein, because conservation is a property of the
        sequence.
    out_dir : str
        Where to write the session.
    reference : str, optional
        Which structure to draw the protein from. Defaults to the first in
        sorted order. Affects only what is *drawn* -- the sites always come from
        every structure.
    site_radius : float
        A water occupies a site within this distance, in Angstrom.
    min_cluster_samples : int
        Minimum waters forming a site.

    Returns
    -------
    str
        Path to the ``.pml`` to open.
    """
    import numpy as np

    from .evolutionary import conservation_of_clusters, write_conservation_report
    from .find_conserved_networks import cluster_coordinates_only
    from .generate_static_networks import initialize_network
    from .visualize_structures import CONSURF_GRADE_COLORS

    warnings.filterwarnings("ignore")
    out = Path(out_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)

    names = sorted(f[:-4] for f in os.listdir(prepared) if f.endswith(".pdb"))
    if not names:
        raise ValueError("no .pdb files in %r -- run `watcon prepare` first" % prepared)
    reference = reference or names[0]
    if reference not in names:
        raise ValueError("reference %r is not in %s (have: %s)"
                         % (reference, prepared, ", ".join(names[:6])))

    # One ConSurf run serves every structure; stage it under each name so the
    # per-structure lookup finds it.
    staging = tempfile.mkdtemp()
    for name in names:
        shutil.copyfile(consurf, os.path.join(staging, name + "_consurf_grades.txt"))

    if verbose:
        print("Building networks for %d structure(s)..." % len(names))
    metrics, networks, _, built = initialize_network(
        prepared, network_type="water-protein", msa_indexing=False,
        classify_water=False, return_network=True, num_workers=1,
        max_distance=3.3, water_name="HOH",
        consurf_directory=staging, consurf_strict=True,
    )

    coordinates = np.array([w.O.coordinates for net in networks
                            for w in net.water_molecules])
    if verbose:
        print("Clustering %d waters..." % len(coordinates))
    _, centers = cluster_coordinates_only(
        coordinates, cluster="hdbscan", min_samples=min_cluster_samples, eps=0.0,
        source=names)
    clusters = conservation_of_clusters(networks, centers, dist_cutoff=site_radius)

    write_conservation_report(clusters, out / "conservation.csv")

    # --- the two objects the session shows --------------------------------
    ordered = (centers.items() if hasattr(centers, "items") else enumerate(centers))
    sites, conserved_sites = [], []
    for cluster_id, centre in ordered:
        record = clusters.get(cluster_id)
        if record is None or record.occupancy == 0:
            continue
        grade = record.max_grade if record.has_conservation else 0
        entry = (cluster_id, centre, grade or 0)
        sites.append(entry)
        if grade and grade >= HIGHLY_CONSERVED:
            conserved_sites.append(entry)

    def write_sites(path, rows):
        with open(path, "w") as handle:
            handle.write("REMARK   Conserved water sites. B-factor = ConSurf grade\n")
            handle.write("REMARK   of the most conserved residue lining the site.\n")
            handle.write("REMARK   B-factor 0 means NO DATA, not 'not conserved'.\n")
            for serial, (cid, centre, grade) in enumerate(rows, 1):
                handle.write(
                    "ATOM  %5d  O   HOH S%4d    %8.3f%8.3f%8.3f  1.00 %6.2f           O\n"
                    % (serial, min(cid, 9999), centre[0], centre[1], centre[2],
                       float(grade)))
            handle.write("END\n")

    all_path = out / "sites_all.pdb"
    top_path = out / "sites_conserved.pdb"
    write_sites(all_path, sites)
    write_sites(top_path, conserved_sites)

    structure = Path(prepared).resolve() / (reference + ".pdb")

    # --- residue colouring, one entry per residue not per atom ------------
    reference_network = networks[built.index(reference)]
    by_residue = {}
    for atom in reference_network.protein_atoms:
        if atom.evolutionary is not None:
            by_residue[(atom.chain, int(atom.resid))] = atom.evolutionary.grade

    by_grade = {}
    for (chain, resid), grade in by_residue.items():
        by_grade.setdefault(grade, []).append((chain, resid))

    def posix(path):
        return str(path).replace("\\", "/")

    pml = out / "watcon_view.pml"
    lines = [
        "# WatCon + ConSurf -- conserved water sites on a conservation-coloured protein.",
        "# Open with:   pymol %s" % pml.name,
        "#",
        "# Protein : %s" % structure.name,
        "# Sites   : %d occupied, %d lined by a grade>=%d residue"
        % (len(sites), len(conserved_sites), HIGHLY_CONSERVED),
        "",
        "bg white",
        "load %s, protein" % posix(structure),
        "load %s, sites" % posix(all_path),
    ]
    if conserved_sites:
        lines.append("load %s, sites_conserved" % posix(top_path))
    lines += [
        "",
        "hide everything",
        "show cartoon, protein",
        "color grey70, protein",
        "",
        "# ConSurf's own scale: 1 cyan (variable) ... 9 maroon (conserved)",
    ]
    for grade in sorted(CONSURF_GRADE_COLORS):
        r, g, b = CONSURF_GRADE_COLORS[grade]
        lines.append("set_color consurf_%d, [%.3f, %.3f, %.3f]" % (grade, r, g, b))
    lines.append("")
    for grade in sorted(by_grade):
        for chain, resid in sorted(by_grade[grade]):
            chain_sel = " and chain %s" % chain if chain.strip() else ""
            lines.append("color consurf_%d, protein and resi %d%s"
                         % (grade, resid, chain_sel))
    lines += [
        "",
        "# Every occupied site, coloured by the grade in its B-factor.",
        "# Off by default: %d spheres bury the protein and the point." % len(sites),
        "# Turn it on with:   enable sites",
        "show spheres, sites",
        "set sphere_scale, 0.25, sites",
        "color grey60, sites",
    ]
    for grade in sorted(CONSURF_GRADE_COLORS):
        lines.append("color consurf_%d, sites and b > %.1f and b < %.1f"
                     % (grade, grade - 0.5, grade + 0.5))
    lines += ["disable sites", ""]

    if conserved_sites:
        lines += [
            "# The sites worth looking at first: lined by a residue ConSurf calls",
            "# highly conserved. That is what the tool is for, so it is what you",
            "# see on opening.",
            "show spheres, sites_conserved",
            "set sphere_scale, 0.40, sites_conserved",
            "color firebrick, sites_conserved",
            "",
        ]
    lines += [
        "set cartoon_transparency, 0.15, protein",
        "orient protein",
        "zoom protein, 4",
        "",
        "print('')",
        "print('  ConSurf grade  9 = most conserved (maroon), 1 = most variable (cyan)')",
        "print('  showing: protein + %d sites lined by a grade>=%d residue (firebrick)')"
        % (len(conserved_sites), HIGHLY_CONSERVED),
        "print('  enable sites    to add all %d occupied sites, coloured by grade')"
        % len(sites),
        "print('  a site B-factor of 0 means NO conservation data, not low conservation')",
        "print('')",
    ]
    with open(pml, "w") as handle:
        handle.write("\n".join(lines) + "\n")

    if verbose:
        print()
        print("Wrote %s" % out)
        print("  watcon_view.pml       open this")
        print("  sites_all.pdb         %d occupied sites, grade in B-factor" % len(sites))
        print("  sites_conserved.pdb   %d sites lined by a grade>=%d residue"
              % (len(conserved_sites), HIGHLY_CONSERVED))
        print("  conservation.csv      the numbers behind the picture")
        print()
        print("  pymol %s" % pml)
    return str(pml)

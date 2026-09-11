"""The conserved-water picture, as data rather than as a file.

``watcon view`` writes a ``.pml``; the PyMOL plugin drives a live session. Those
are two consumers of one thing, and if each built its own scene they would drift
until the file and the plugin disagreed about what the data shows -- a whole
class of bug this project has already been bitten by twice (the duplicated
contact walk, and the two pickle shapes).

So the scene is built once, here, as a :class:`Scene`: the files it wrote, the
per-site records behind it, and an ordered list of PyMOL commands. Writing those
commands to a file gives you ``watcon view``; feeding them to ``cmd.do`` gives
you the plugin. Neither contains any scene logic of its own, and this module
imports PyMOL not at all, so it is testable everywhere -- which matters, because
CI has no PyMOL.

What the picture shows, and why
-------------------------------
* **Protein**, cartoon, on ConSurf's own 1-9 scale. Grades 3-6 are genuinely
  near-white on that scale; that is ConSurf's convention, not a washed-out
  palette.
* **Residues with no ConSurf score are yellow**, again ConSurf's convention.
  They used to be ``grey70``, which is almost exactly the near-white of an
  average grade -- so "we have no data here" and "this is averagely conserved"
  looked alike. They are opposite claims.
* **Conserved sites**, shown on opening: water positions lined by a residue
  ConSurf calls highly conserved.
* **Every occupied site**, loaded but disabled -- 190 spheres bury the protein.
* **Sphere radius tracks occupancy**, so a site found in every structure reads
  as more than one found in a single structure.
* **The contacts**, off by default: sticks for the residues WatCon says line each
  conserved site, and dashes for the polar contacts among them.
"""

from __future__ import annotations

import os
import shutil
import tempfile
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

#: Grade at or above which a site is called highly conserved. 8 because that is
#: the most reproducible end of a ConSurf run -- comparing two runs of the same
#: sequence, grade >= 8 was far more stable than the variable end.
HIGHLY_CONSERVED = 8

#: ConSurf's colour for "insufficient data". Distinct from every grade colour,
#: deliberately: absence of evidence must not look like evidence of averageness.
NO_DATA_COLOR = (1.00, 1.00, 0.59)

#: Sphere radius in Angstrom for a site occupied in no structures ... all of them.
SITE_RADIUS_RANGE = (0.25, 0.70)

#: Distance cutoff, Angstrom, for drawing polar contacts in PyMOL.
CONTACT_CUTOFF = 3.6


@dataclass
class Site:
    """One conserved water site, with everything a table or a label needs."""

    cluster_id: int
    centre: Tuple[float, float, float]
    grade: int                      # 0 means NO DATA, not "not conserved"
    occupancy: int
    n_structures_occupied: int
    n_structures_total: int
    residues: Tuple = ()            # ((chain, resid, icode, grade or None), ...)

    @property
    def occupancy_fraction(self) -> float:
        if self.n_structures_total == 0:
            return 0.0
        return self.n_structures_occupied / self.n_structures_total

    @property
    def has_conservation(self) -> bool:
        return self.grade > 0

    # Deliberately no `conserved` property here. It would have to read the
    # module constant, which silently disagrees with a Scene built at a
    # different threshold -- exactly the bug this arrangement replaced. The
    # threshold belongs to the Scene, so the Scene does the filtering.

    def residue_label(self, limit: int = 6) -> str:
        """``A:27, A:54, ...`` -- what lines this site, for a table cell."""
        names = ["%s:%d" % (chain, resid) for chain, resid, _icode, _g
                 in self.residues[:limit]]
        return ", ".join(names) + (", ..." if len(self.residues) > limit else "")


@dataclass
class Scene:
    """A built picture: the files, the records, and how to display them."""

    out_dir: str
    pml_path: str
    structure_path: str
    sites_all_path: str
    sites_conserved_path: str
    report_path: str

    reference: str
    names: List[str] = field(default_factory=list)
    sites: List[Site] = field(default_factory=list)
    clusters: Dict = field(default_factory=dict)
    commands: List[str] = field(default_factory=list)
    n_clusters: int = 0
    highly_conserved: int = HIGHLY_CONSERVED

    @property
    def conserved_sites(self) -> List[Site]:
        return [s for s in self.sites if s.grade >= self.highly_conserved]

    @property
    def n_with_conservation(self) -> int:
        return sum(1 for s in self.sites if s.has_conservation)

    def summary(self) -> str:
        return ("%d clusters / %d occupied / %d with conservation / %d conserved "
                "(grade >= %d)"
                % (self.n_clusters, len(self.sites), self.n_with_conservation,
                   len(self.conserved_sites), self.highly_conserved))

    def write_pml(self, path: Optional[str] = None) -> str:
        """Write the commands as a self-contained ``.pml``. Returns its path."""
        path = path or self.pml_path
        with open(path, "w") as handle:
            handle.write("\n".join(self.commands) + "\n")
        return path


def _posix(path) -> str:
    """PyMOL reads forward slashes on every platform; Windows paths do not."""
    return str(path).replace("\\", "/")


def _write_sites(path, sites: Sequence[Site]) -> None:
    """Sites as a PDB: grade in the B-factor, occupancy fraction in occupancy.

    The occupancy column is what lets the radius track how often a site is
    actually seen, in a single PyMOL ``alter`` rather than one command per site.

    The column layout is exact here. The previous writer put the B-factor one
    column late; it happened to survive parsing because the trailing digit fell
    off the end of the field, which is the kind of thing that works until a
    value changes width.
    """
    with open(path, "w") as handle:
        handle.write("REMARK   WatCon + ConSurf conserved water sites.\n")
        handle.write("REMARK   B-factor  = ConSurf grade of the most conserved\n")
        handle.write("REMARK               residue lining the site (1-9).\n")
        handle.write("REMARK               0 means NO DATA, not 'not conserved'.\n")
        handle.write("REMARK   occupancy = fraction of structures in which a\n")
        handle.write("REMARK               water occupies this site (0-1).\n")
        for serial, site in enumerate(sites, 1):
            x, y, z = site.centre[0], site.centre[1], site.centre[2]
            handle.write(
                "ATOM  %5d  O   HOH S%4d    %8.3f%8.3f%8.3f%6.2f%6.2f           O\n"
                % (serial, min(site.cluster_id, 9999), x, y, z,
                   site.occupancy_fraction, float(site.grade)))
        handle.write("END\n")


def build_scene(prepared, consurf, out_dir="watcon_view", reference=None,
                site_radius=1.5, min_cluster_samples=2,
                max_distance=3.3, water_name="HOH",
                network_type="water-protein", num_workers=1, eps=0.0,
                highly_conserved=HIGHLY_CONSERVED,
                consurf_chain_map=None, progress=None, verbose=True) -> Scene:
    """Run the analysis and describe the picture. Returns a :class:`Scene`.

    Parameters
    ----------
    prepared : str
        Directory of prepared structures -- output of ``watcon prepare``.
    consurf : str
        A ConSurf ``*_consurf_grades.txt`` file. One run describes every
        structure of the same protein, because conservation is a property of the
        sequence, not of the crystal.
    out_dir : str
        Where to write the session.
    reference : str, optional
        Which structure to draw the protein from. Defaults to the first in
        sorted order. Affects only what is *drawn* -- sites always come from
        every structure.
    site_radius : float
        A water occupies a site within this distance, in Angstrom.
    min_cluster_samples : int
        Minimum waters forming a site.
    max_distance : float
        Hydrogen-bond distance cutoff for building the network, in Angstrom.
    water_name : str
        Residue name of water in these files. ``None`` accepts the usual set.
    network_type : str
        ``'water-protein'`` or ``'water-water'``.
    num_workers : int
        Parallel workers for network building.
    eps : float
        Clustering epsilon.
    highly_conserved : int
        Grade at or above which a site is called conserved.
    consurf_chain_map : dict, optional
        ``{consurf_chain: structure_chain}`` when the run's chain labels differ.
    progress : callable, optional
        Called as ``progress(fraction, message)``. The plugin uses it to drive a
        progress bar; the CLI leaves it None.
    verbose : bool
        Print as it goes.

    Notes
    -----
    This function does not import PyMOL and never will. It is the whole of the
    scene logic, so the ``.pml`` and the plugin cannot disagree.
    """
    import numpy as np

    from .evolutionary import conservation_of_clusters, write_conservation_report
    from .find_conserved_networks import cluster_coordinates_only
    from .generate_static_networks import initialize_network
    from .visualize_structures import CONSURF_GRADE_COLORS

    def step(fraction, message):
        if progress is not None:
            progress(fraction, message)
        if verbose:
            print(message)

    warnings.filterwarnings("ignore")
    out = Path(out_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)

    names = sorted(f[:-4] for f in os.listdir(prepared) if f.endswith(".pdb"))
    if not names:
        raise ValueError(
            "No prepared structures in %s. This should be the output of "
            "`watcon prepare`, which writes one .pdb per accepted structure."
            % prepared)
    reference = reference or names[0]
    if reference not in names:
        raise ValueError("reference %r is not in %s (have: %s)"
                         % (reference, prepared, ", ".join(names[:6])))

    # One ConSurf run serves every structure; stage it under each name so the
    # per-structure lookup finds it.
    staging = tempfile.mkdtemp(prefix="watcon_consurf_")
    try:
        for name in names:
            shutil.copyfile(consurf,
                            os.path.join(staging, name + "_consurf_grades.txt"))

        step(0.05, "Building networks for %d structure(s)..." % len(names))
        _metrics, networks, _centers, built = initialize_network(
            prepared, network_type=network_type, msa_indexing=False,
            classify_water=False, return_network=True, num_workers=num_workers,
            max_distance=max_distance, water_name=water_name,
            consurf_directory=staging, consurf_chain_map=consurf_chain_map,
            consurf_strict=True,
        )
    finally:
        shutil.rmtree(staging, ignore_errors=True)

    coordinates = np.array([w.O.coordinates for net in networks
                            for w in net.water_molecules])
    step(0.55, "Clustering %d waters..." % len(coordinates))
    _labels, centers = cluster_coordinates_only(
        coordinates, cluster="hdbscan", min_samples=min_cluster_samples,
        eps=eps, source=names)

    step(0.75, "Joining %d sites to conservation..." % len(centers))
    clusters = conservation_of_clusters(networks, centers,
                                        dist_cutoff=site_radius)
    report_path = out / "conservation.csv"
    write_conservation_report(clusters, report_path)

    # --- conservation per residue of the reference, for colouring and labels --
    reference_network = networks[built.index(reference)]
    grade_of_residue = {}
    for atom in reference_network.protein_atoms:
        if atom.evolutionary is not None:
            grade_of_residue[(atom.chain, int(atom.resid))] = atom.evolutionary.grade

    # --- the sites ---------------------------------------------------------
    step(0.85, "Describing the scene...")
    ordered = (centers.items() if hasattr(centers, "items")
               else enumerate(centers))
    sites: List[Site] = []
    for cluster_id, centre in ordered:
        record = clusters.get(cluster_id)
        if record is None or record.occupancy == 0:
            continue
        residues = []
        for key in record.residue_keys:
            chain, resid = key[0], int(key[1])
            icode = key[2] if len(key) > 2 else None
            residues.append((chain, resid, icode,
                             grade_of_residue.get((chain, resid))))
        residues.sort(key=lambda r: (r[0], r[1]))
        sites.append(Site(
            cluster_id=cluster_id,
            centre=(float(centre[0]), float(centre[1]), float(centre[2])),
            grade=int(record.max_grade) if record.has_conservation else 0,
            occupancy=record.occupancy,
            n_structures_occupied=record.n_structures_occupied,
            n_structures_total=record.n_structures_total,
            residues=tuple(residues),
        ))

    conserved = [s for s in sites if s.grade >= highly_conserved]

    all_path = out / "sites_all.pdb"
    top_path = out / "sites_conserved.pdb"
    _write_sites(all_path, sites)
    _write_sites(top_path, conserved)

    structure = Path(prepared).resolve() / (reference + ".pdb")

    scene = Scene(
        out_dir=str(out), pml_path=str(out / "watcon_view.pml"),
        structure_path=str(structure), sites_all_path=str(all_path),
        sites_conserved_path=str(top_path), report_path=str(report_path),
        reference=reference, names=list(names), sites=sites, clusters=clusters,
        n_clusters=len(centers), highly_conserved=highly_conserved,
    )
    scene.commands = _commands_for(scene, grade_of_residue,
                                   CONSURF_GRADE_COLORS)
    step(1.0, "Done.")
    return scene


def _commands_for(scene: Scene, grade_of_residue, palette) -> List[str]:
    """The ordered PyMOL commands that draw ``scene``. No PyMOL import."""
    conserved = scene.conserved_sites
    low, high = SITE_RADIUS_RANGE

    lines = [
        "# WatCon + ConSurf -- conserved water sites on a conservation-coloured protein.",
        "# Open with:   pymol %s" % os.path.basename(scene.pml_path),
        "#",
        "# Protein : %s" % os.path.basename(scene.structure_path),
        "# Sites   : %d occupied, %d lined by a grade>=%d residue"
        % (len(scene.sites), len(conserved), scene.highly_conserved),
        "",
        "bg white",
        "load %s, protein" % _posix(scene.structure_path),
        "load %s, sites" % _posix(scene.sites_all_path),
    ]
    if conserved:
        lines.append("load %s, sites_conserved" % _posix(scene.sites_conserved_path))

    lines += [
        "",
        "hide everything",
        "show cartoon, protein",
        "",
        "# ConSurf's own scale: 1 cyan (variable) ... 9 maroon (conserved).",
        "# Grades 3-6 really are near-white on this scale.",
    ]
    for grade in sorted(palette):
        r, g, b = palette[grade]
        lines.append("set_color consurf_%d, [%.3f, %.3f, %.3f]" % (grade, r, g, b))
    lines.append("set_color consurf_nodata, [%.3f, %.3f, %.3f]" % NO_DATA_COLOR)
    lines += [
        "",
        "# Yellow = no ConSurf score for this residue, ConSurf's own convention.",
        "# It was grey70, which is almost the near-white of an average grade --",
        "# so 'no data' and 'averagely conserved' looked alike. They are opposite",
        "# claims.",
        "color consurf_nodata, protein",
        "",
    ]

    by_grade = {}
    for (chain, resid), grade in grade_of_residue.items():
        by_grade.setdefault(grade, []).append((chain, resid))
    for grade in sorted(by_grade):
        for chain, resid in sorted(by_grade[grade]):
            chain_sel = " and chain %s" % chain if chain.strip() else ""
            lines.append("color consurf_%d, protein and resi %d%s"
                         % (grade, resid, chain_sel))

    lines += [
        "",
        "# Every occupied site, coloured by the grade in its B-factor.",
        "# Off by default: %d spheres bury the protein and the point." % len(scene.sites),
        "# Turn it on with:   enable sites",
        "show spheres, sites",
        "color consurf_nodata, sites",
    ]
    for grade in sorted(palette):
        lines.append("color consurf_%d, sites and b > %.1f and b < %.1f"
                     % (grade, grade - 0.5, grade + 0.5))
    lines += [
        "",
        "# Radius tracks occupancy: the fraction of structures holding a water",
        "# here is in the occupancy column, so one alter does every sphere.",
        "alter sites, vdw=%.2f+%.2f*q" % (low, high - low),
        "set sphere_scale, 1.0, sites",
        "rebuild sites",
        "disable sites",
        "",
    ]

    if conserved:
        lines += [
            "# The sites worth looking at first: lined by a residue ConSurf calls",
            "# highly conserved. That is what the tool is for, so it is what you",
            "# see on opening.",
            "show spheres, sites_conserved",
            "color firebrick, sites_conserved",
            "alter sites_conserved, vdw=%.2f+%.2f*q" % (low + 0.15, high - low),
            "set sphere_scale, 1.0, sites_conserved",
            "rebuild sites_conserved",
            "",
            "# Which residues line them, and the polar contacts among those.",
            "# WatCon chose the residues. PyMOL only decides which atom pair",
            "# each dash connects. Off by default -- it is dense.",
        ]
        lining = sorted({(chain, resid)
                         for site in conserved
                         for chain, resid, _icode, _g in site.residues})
        by_chain = {}
        for chain, resid in lining:
            by_chain.setdefault(chain, []).append(resid)
        parts = []
        for chain in sorted(by_chain):
            resids = "+".join(str(r) for r in sorted(by_chain[chain]))
            parts.append("(chain %s and resi %s)" % (chain, resids)
                         if chain.strip() else "(resi %s)" % resids)
        if parts:
            lines += [
                "select lining_residues, protein and (%s)" % " or ".join(parts),
                "show sticks, lining_residues and sidechain",
                "color grey40, lining_residues and elem C",
                "distance site_contacts, sites_conserved, lining_residues, "
                "%.1f, mode=2" % CONTACT_CUTOFF,
                "set dash_color, black, site_contacts",
                "set dash_gap, 0.3, site_contacts",
                "hide labels, site_contacts",
                "group WatCon_contacts, lining_residues site_contacts",
                # The group has its own enabled flag -- disabling the members
                # leaves the container switched on, so `enable WatCon_contacts`
                # would appear to do nothing the first time it is used.
                "disable WatCon_contacts",
                "disable lining_residues",
                "disable site_contacts",
                "deselect",
            ]
        lines += [
            "",
            "# Site numbers, off by default.",
            'label sites_conserved, "S%s" % resi',
            "set label_size, 14",
            "set label_color, black",
            "hide labels, sites_conserved",
            "",
        ]

    lines += [
        "orient protein",
        "zoom protein, 4",
        "",
        "print('')",
        "print('  WatCon + ConSurf')",
        "print('  protein colour  9 maroon = most conserved ... 1 cyan = most variable')",
        "print('                  yellow = NO ConSurf score for that residue')",
        "print('  showing         %d site(s) lined by a grade>=%d residue, in firebrick')"
        % (len(conserved), scene.highly_conserved),
        "print('  sphere size     larger = occupied in more structures')",
        "print('')",
        "print('  enable sites              all %d occupied sites, coloured by grade')"
        % len(scene.sites),
        "print('  enable WatCon_contacts    the residues lining them, and polar contacts')",
        "print('  show labels, sites_conserved   site numbers')",
        "print('')",
        "print('  a site B-factor of 0 means NO conservation data, not low conservation')",
        "print('')",
    ]
    return lines

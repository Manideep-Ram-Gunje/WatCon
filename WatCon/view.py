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
    WatCon_contacts   the residues lining them, and the polar contacts

so that ``pymol watcon_view.pml`` opens on the picture the tool exists to make.

The scene itself lives in :mod:`WatCon.scene`, because the PyMOL plugin draws
the same picture live and the two must not drift apart. This module is now just
"build the scene, write it to a file".
"""

from __future__ import annotations

from .scene import HIGHLY_CONSERVED, Scene, build_scene

__all__ = ["HIGHLY_CONSERVED", "Scene", "build_scene", "build_session"]


def build_session(prepared, consurf, out_dir="watcon_view", reference=None,
                  site_radius=1.5, min_cluster_samples=2, verbose=True,
                  **kwargs):
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
    **kwargs
        Passed to :func:`WatCon.scene.build_scene` -- ``max_distance``,
        ``water_name``, ``network_type``, ``num_workers``, ``eps``,
        ``highly_conserved``, ``consurf_chain_map``.

    Returns
    -------
    str
        Path to the ``.pml`` to open.
    """
    scene = build_scene(prepared, consurf, out_dir=out_dir, reference=reference,
                        site_radius=site_radius,
                        min_cluster_samples=min_cluster_samples,
                        verbose=verbose, **kwargs)
    path = scene.write_pml()

    if verbose:
        print()
        print("Wrote %s" % scene.out_dir)
        print("  watcon_view.pml       open this")
        print("  sites_all.pdb         %d occupied sites, grade in B-factor"
              % len(scene.sites))
        print("  sites_conserved.pdb   %d sites lined by a grade>=%d residue"
              % (len(scene.conserved_sites), scene.highly_conserved))
        print("  conservation.csv      the numbers behind the picture")
        print()
        print("  pymol %s" % path)
    return path

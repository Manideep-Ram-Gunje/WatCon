WatCon Package
==============

WatCon is organized in the following manner:

.. figure:: images/workflow.png
   :width: 800

   Image sourced from DOI: 10.1021/jacsau.5c00447 

Details on the contents of WatCon are outlined below.


Subpackages
-----------

.. toctree::
   :maxdepth: 4

   source/modules/WatCon.tests

Submodules
----------

.. toctree::
   :maxdepth: 2

   source/modules/WatCon
   source/modules/WatCon.sequence_processing
   source/modules/WatCon.generate_dynamic_networks
   source/modules/WatCon.generate_static_networks
   source/modules/WatCon.find_conserved_networks
   source/modules/WatCon.residue_analysis
   source/modules/WatCon.visualize_structures

Evolutionary conservation (ConSurf)
-----------------------------------

Added by the ConSurf extension.  ``evolutionary`` holds the join between water
sites and per-residue conservation; ``residue_index`` the ``(chain, resid,
icode)`` identity mapping everything is keyed on; ``superpose`` rigid-body
alignment without MODELLER; ``prepare`` the dataset preparation step; and the
``consurf`` subpackage the grades-file parser and its diagnostics.

.. toctree::
   :maxdepth: 2

   source/modules/WatCon.evolutionary
   source/modules/WatCon.residue_index
   source/modules/WatCon.superpose
   source/modules/WatCon.prepare
   source/modules/WatCon.structure_io
   source/modules/WatCon.fetch
   source/modules/WatCon.conformers
   source/modules/WatCon.cli
   source/modules/WatCon.demo
   source/modules/WatCon.consurf
   source/modules/WatCon.consurf.parser
   source/modules/WatCon.consurf.model
   source/modules/WatCon.consurf.dialects
   source/modules/WatCon.consurf.crosscheck
   source/modules/WatCon.consurf.validate
   source/modules/WatCon.consurf.errors


Protein families
----------------

One ConSurf run per protein, several structures each, one shared alignment.
``alignment`` places every structure on its row **by sequence** and cross-checks
structures of the same protein against each other; ``family`` pools conservation
onto alignment columns; ``family_sites`` superposes the family and finds the
water sites they share; ``family_scene`` draws the result.

.. toctree::
   :maxdepth: 2

   source/modules/WatCon.alignment
   source/modules/WatCon.family
   source/modules/WatCon.family_sites
   source/modules/WatCon.family_scene

Viewing results
---------------

``scene`` turns an analysis into plain data -- sites, their conservation and the
commands that would draw them -- so the command line and the PyMOL plugin cannot
disagree about what a result looks like.  ``view`` writes that scene as a
self-contained ``.pml``; ``pymol_plugin`` runs the identical file inside PyMOL.

.. toctree::
   :maxdepth: 2

   source/modules/WatCon.scene
   source/modules/WatCon.view
   source/modules/WatCon.pymol_plugin

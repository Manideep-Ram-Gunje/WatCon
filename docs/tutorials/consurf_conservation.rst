Evolutionary Conservation of Water Sites
=========================================

WatCon finds water sites that recur across many structures of a protein.
ConSurf scores how conserved each residue is across evolution. This tutorial
joins the two and asks whether they agree.

The worked example is the study in ``experiments/barnase_waters/``: 22 barnase
crystal structures, one ConSurf run. Everything here runs from the command line.

.. contents:: :local:


The idea that makes this cheap
------------------------------

**Conservation is a property of the sequence, not of a crystal.** One ConSurf
run therefore describes *every* structure of that protein — all 22 barnase
entries were covered by a single run. See :doc:`../consurf_data` for how to
obtain one.


Try it first
------------

.. code-block:: bash

   watcon demo

Six barnase structures ship with the package. The whole pipeline runs offline in
about ten seconds, and the output directory contains everything the rest of this
tutorial explains.


Step 1 — collect structures
---------------------------

Any PDB files of the same protein. For barnase, all X-ray entries for UniProt
P00648 at ≤2.0 Å — 24 of 47. ``experiments/barnase_waters/scripts/fetch_structures.py``
does this from RCSB.

Resolution matters: water positions are the observable, and they are only as
good as the map.


Step 2 — prepare them
---------------------

.. code-block:: bash

   watcon prepare --input-dir structures/ --out-dir prepared/ \
                  --reference 1A2P --report preparation.csv

This finds the right chain in each structure, superposes everything onto the
reference, carries the waters along, and relabels chains so one ConSurf file
serves them all::

   1A2P     chain A  identity 1.00   148 waters  RMSD 0.00 A over 108 CA
   1BRN     chain L  identity 1.00   120 waters  RMSD 0.45 A over 108 CA
   1BRS     chain A  identity 1.00   174 waters  RMSD 0.46 A over 108 CA

Note ``1BRN`` — barnase is chain **L** there. Structures are routinely deposited
as complexes, and assuming chain 'A' quietly analyses the wrong protein.

Read the rejections
~~~~~~~~~~~~~~~~~~~

::

   2F56     REJECTED (identity 0.02 by residue number, below 0.80)
            [chain A matches at 1.00 if renumbered by +2 -- renumber it
             yourself if that is correct; this is not applied automatically]

2F56 numbers the residue everyone else calls 3 as 1. Every ConSurf score would
land on the wrong residue. WatCon diagnoses it and refuses; it does not renumber
your file for you.

.. note::

   Chains are compared by **residue number**, not by sequence position. On this
   real dataset the positional alternative failed both ways: it rejected three
   good structures (1BSA, 1B2S, 1RNB, which start at residues 4, 1 and 2) and
   accepted the two mis-numbered ones, which then superposed at 6.05 Å.

No MODELLER needed — ``WatCon.superpose`` does this in numpy.


Step 3 — validate your ConSurf file
-----------------------------------

.. code-block:: bash

   watcon validate --consurf consurf/1A2P_consurf_grades.txt

Do this before the real run. See :doc:`../consurf_data`.


Step 4 — build networks with conservation attached
--------------------------------------------------

In your input file:

.. code-block:: text

   structure_directory: prepared
   network_type: water-protein
   max_distance: 3.3
   cluster_coordinates: on
   clustering_method: hdbscan
   min_cluster_samples: 2

   consurf_directory: consurf
   consurf_strict: on

.. code-block:: bash

   watcon run --input input.txt --name barnase

``min_cluster_samples: 2`` is deliberate — low-occupancy sites are the baseline
you compare conserved sites *against*, so excluding them deletes your control
group.

Every structure's coverage and identity is reported::

   entry      waters   coverage   identity
   1A2P          148       100%       100%
   1BSA           81       100%        99%

Those measure different things. **Coverage** is how much of the structure the
run reached — 1BRS reads 18.9% because ConSurf covered one of its six chains.
**Identity** is whether the residues it did reach are the ones present. A run
can cover 100% of a structure and still describe a different protein.


Step 5 — join sites to conservation
-----------------------------------

In your analysis file:

.. code-block:: text

   input_directory: watcon_output
   cluster_concatenated: on
   conservation_report: on
   conservation_dist_cutoff: 1.5

.. code-block:: bash

   watcon run --analysis analysis.txt

This writes one row per conserved water site:

=====================  ==========================================
column                  meaning
=====================  ==========================================
``occupancy``           waters found at this site
``n_structures_occupied``  structures in which it is occupied
``occupancy_fraction``  **structural** conservation
``evo_min_score``       most conserved lining residue
``evo_mean_score``      mean over lining residues
``evo_max_grade``       highest ConSurf grade, 1–9
``evo_n_residues``      residues lining the site
``evo_n_low_confidence``  lining residues ConSurf flagged
``evo_n_unscored``      lining residues with no score
=====================  ==========================================

.. important::

   Structural and evolutionary conservation stay in **separate columns**. Whether
   they agree is the scientific question, so WatCon will not blend them into one
   number. ``NA`` means no data — never "not conserved".


Step 6 — look at it
-------------------

.. code-block:: bash

   pymol cluster_pdbs/*_evolutionary.pdb pymol_projections/*.pml

The cluster PDB carries the ConSurf grade in the B-factor column; the ``.pml``
colours residues on ConSurf's own 1–9 scale. WatCon's existing
``project_clusters`` writes *structural* conservation into the same column of a
*different* file, so you can load both and compare.


What the barnase study found
----------------------------

462 sites carried conservation data.

  Spearman ρ = **−0.375**, p = 6.6×10⁻¹⁷, against a permutation null centred on
  +0.033 (z = −4.47). Stable at ≤1.8 Å (−0.332) and wild-type only (−0.338).

Recurring water sites *are* lined by more conserved residues. Scores are more
negative when more conserved, hence the negative correlation.

The statistic had to be corrected
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The pre-registered choice was ``evo_min_score``, which gives a more emphatic
−0.501. Its shuffled-conservation null centres on **−0.226**, not zero:

- Spearman(occupancy, number of lining residues) = **+0.555**
- Spearman(number of lining residues, min score) = **−0.553**

The minimum of a larger set is mechanically lower, and busier sites have more
lining residues — so much of −0.50 was the statistic measuring its own set size.
``evo_mean_score`` has no such dependence and its null does centre on zero.

.. tip::

   Always run a shuffled-conservation null. Permute conservation across residues,
   keep the geometry fixed, recompute. If the null does not centre on zero, your
   statistic is measuring something structural before any biology enters.

What this does not show
~~~~~~~~~~~~~~~~~~~~~~~

Buried residues are both more conserved (median −0.694 vs +0.180) and more
likely to hold ordered water. This design cannot separate "conserved residues
hold water" from "buried positions are conserved *and* hold water", so **no
causal claim is made**. It is one protein, and it is correlational.


Beyond one protein
------------------

For a true family, each member needs its own ConSurf run, and residue 40 of one
protein is not residue 40 of another — the correspondence runs through the MSA:

.. code-block:: python

   from WatCon.evolutionary import conservation_by_msa_column, family_summary

   columns = conservation_by_msa_column([
       ("member1", map1, index1),
       ("member2", map2, index2),
   ])

.. caution::

   ConSurf scores are z-normalised **within each run**, so pooled values are
   only relatively comparable, and grades are per-run bins that must never be
   averaged. Two runs of the *same* barnase sequence from different starting
   structures agree at only ρ ≈ 0.37 — that is the floor. Prefer
   ``unanimous_conserved``, which uses each run's own verdict instead of
   comparing separately normalised numbers.

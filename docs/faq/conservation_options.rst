Conservation Options
--------------------

Every setting the ConSurf extension adds, and what each one refuses to guess.

Conservation is a property of the **sequence**, so one ConSurf run describes
every structure of that protein.  See :doc:`../consurf_data` for how to obtain a
run, and :doc:`../tutorials/consurf_conservation` for a worked example.


Input file settings
~~~~~~~~~~~~~~~~~~~

These go in the input file (``input_static.txt`` / ``input_dynamic.txt``):

.. code-block:: text

    consurf_directory: consurf/     ; directory of *_consurf_grades.txt files
    consurf_strict: on              ; missing or unparseable file is an error
    consurf_chain_map: None         ; 'A:A,B:C' ConSurf-chain -> structure-chain

``consurf_directory``
    ``None`` disables conservation entirely.  Otherwise WatCon looks for a
    grades file whose name matches each structure.  Nothing here contacts the
    ConSurf server; the files are supplied by you.

``consurf_strict``
    ``on`` stops the run when a file is missing or unreadable.  ``off`` warns
    and continues with every residue unscored -- useful when only some of your
    structures have runs.

``consurf_chain_map``
    Needed when the ConSurf run and the structure label chains differently.
    A run on chain A joined to a file with no chain column at all -- a
    GROMACS ``.gro``, for instance -- needs ``A:``.

And in the analysis file:

.. code-block:: text

    conservation_report: on         ; write the site-by-site CSV


What WatCon checks for you
~~~~~~~~~~~~~~~~~~~~~~~~~~

The join key is ``(chain, resid, icode)``, which says nothing about *which amino
acid* a score belongs to.  A run numbered differently from your structure would
attach a score to every residue, report 100% coverage, and be wrong everywhere
-- invisibly.  So two things are verified.

**Residue identity.**  Every matched residue is compared with the amino acid
ConSurf scored there.  Below 95% agreement the run stops:

.. code-block:: text

    ConSurf data does not describe 2F56.pdb: only 2 of 106 compared residues
    agree (1.9%, threshold 95%). ... A3: structure has ASN, ConSurf has VAL

It is a *rate*, not a pass/fail, because a point mutant genuinely differs at one
position (99.1%) while a numbering offset collapses to 2.8%.  Force-field
protonation names (``HID``, ``HIE``, ``HIP``, ``CYM``, ``GLH`` …) are normalised
first, so they do not read as differences; chemical modifications such as
phosphocysteine ``CSP`` deliberately still do.

**Coverage.**  How much of the structure the run described.  Low coverage is
legitimate -- a run covering one chain of several -- so it is reported rather
than enforced.  Coverage of *zero* is never legitimate, and warns:

.. code-block:: text

    Warning: ConSurf data was supplied for system.gro but matched none of its
    residues, so nothing is scored. The structure has chain(s) (blank) and the
    ConSurf run describes chain(s) A. If those differ, pass consurf_chain_map ...


Checking a file before you rely on it
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    watcon validate --consurf my_run_consurf_grades.txt

This reports the dialect, the calculation method, the MSA depth, how many
residues were scored and how many could be placed on a structure, and any
malformed lines -- before the run rather than during it.


What you get back
~~~~~~~~~~~~~~~~~

A CSV with one row per conserved water site: its structural occupancy beside the
evolutionary conservation of the residues lining it, as **two separate
columns**.  Deciding how they relate is the science, so the tool does not blend
them into one number.


Across a protein family
~~~~~~~~~~~~~~~~~~~~~~~

``watcon family`` takes one ConSurf run per protein and a shared alignment, and
pools conservation onto alignment columns.  Because scores are z-normalised
within each run, separate runs are only relatively comparable -- so the family
result reports each member's grade side by side and counts *verdicts*
(``unanimous_conserved``: every run grading a column 8 or 9) rather than
averaging values that were never on a common scale.

.. code-block:: bash

    watcon family --members members.tsv --alignment alignment.pir

``members.tsv`` is tab-separated, one protein per line:

.. code-block:: text

    # protein   structures directory   ConSurf grades file   [reference]
    PTPN1       prepared/PTPN1         consurf/1AAX_A.grades.txt   2F71

Tab-separated rather than colon-separated because Windows paths contain colons.

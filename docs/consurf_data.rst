Getting ConSurf Data
====================

This is the one step WatCon cannot do for you, so it is worth doing carefully.

**There is no public ConSurf API.** The ConSurf server has no documented
programmatic interface, its terms restrict automated use, and ConSurf-DB was
unreachable when this was written. Scripting the web form would be fragile and
discourteous, so WatCon deliberately does **not** submit anything. You run
ConSurf yourself, once, and drop the result in a folder.

That is less painful than it sounds, because of one fact:

.. important::

   **Conservation is a property of the sequence, not of a crystal.**

   One ConSurf run therefore describes *every* structure of that protein. The
   22-structure barnase study in this repository used exactly **one** ConSurf
   run. You do not need one per PDB file.

You need a new run only for a genuinely different protein — a different member
of a protein family.


1. Submit your sequence
-----------------------

Go to https://consurf.tau.ac.il/ and choose the protein pathway.

- **Upload a structure** if you have one. ConSurf then returns a grades file
  whose ``ATOM`` column carries real PDB residue numbers and chain — which is
  what WatCon joins on, so this is the better option.
- **Paste a sequence** if you have no structure. This still works; the mapping
  falls back to the one-letter ``SEQ`` column.

Defaults are reasonable. Two choices matter:

- **Number of sequences in the MSA.** More is better. Our two runs of the same
  barnase sequence at depths 150 and 50 agree at Spearman ρ = 0.955 on the score
  but the *grade* differs at 42% of positions. Shallow alignments are noisy.
- **Bayesian vs maximum likelihood.** Bayesian gives confidence intervals, which
  WatCon records and reports. Maximum likelihood does not; WatCon handles both.


2. Keep the right file
----------------------

ConSurf returns several files. The one you want is the **grades** file — the
one containing a table with ``POS``, ``SEQ``, ``ATOM``, ``SCORE``, ``COLOR``
columns::

   POS  SEQ  ATOM      SCORE    COLOR   CONFIDENCE INTERVAL ...
     3   V   VAL:3:A   -0.729      7    -1.017,-0.523 ...

Typically named ``consurf_grades.txt`` or ``<job>_consurf_grades.txt``. The
others — the annotated PDB, the MSA, the tree, the coloured images — are not
needed, though the annotated PDB is useful for cross-checking (see below).


3. Name it after your structure
-------------------------------

WatCon matches ConSurf files to structures by name, using the same convention it
already uses for FASTA files. For structures ``1A2P.pdb`` and ``1BRS.pdb``::

   consurf/
     1A2P_consurf_grades.txt
     1BRS_consurf_grades.txt

Because one run covers every structure of a protein, those are simply **copies
of the same file** under different names. ``watcon demo`` does exactly this.

.. warning::

   **Each structure must match exactly one file.** The lookup tries the full
   name first, then the part before the first underscore, and **refuses** when a
   token matches more than one file rather than picking arbitrarily.

   This is not hypothetical. An earlier version tried the loose token first, so
   ``P00648_50`` silently resolved to ``P00648_150.grades.txt`` — a different
   run's conservation, attached with no error anywhere.


4. Validate it before you rely on it
------------------------------------

.. code-block:: bash

   watcon validate --consurf consurf/1A2P_consurf_grades.txt

Reports the dialect, method, MSA depth, and how many records carry a structural
mapping::

     dialect              webserver
     method               bayesian
     msa_total            150
     records              110
     mapped               108
     unmapped             2

``unmapped`` records are normal: ConSurf scores the whole sequence, and residues
unresolved in the crystal have no ``ATOM`` entry.

If you kept ConSurf's annotated PDB, cross-check against it. ConSurf writes the
grade into the B-factor column, so this verifies the parse against ConSurf's own
output:

.. code-block:: bash

   watcon validate --consurf grades.txt --pdb annotated.pdb


5. Point WatCon at it
---------------------

In your input file:

.. code-block:: text

   consurf_directory: consurf     ; folder of *_consurf_grades.txt files
   consurf_strict: on             ; a missing or unparseable file is an error
   consurf_chain_map: None        ; optional 'A:B' ConSurf-chain -> structure-chain

and in your analysis file:

.. code-block:: text

   conservation_report: on
   conservation_dist_cutoff: 1.5

Leave ``consurf_strict: on``. Tolerant mode continues without conservation,
which is occasionally what you want and is otherwise an excellent way to produce
an empty report and not notice.


What WatCon checks for you
--------------------------

The join key is ``(chain, resid, icode)``. Nothing in it records *which amino
acid* a score belongs to, so a ConSurf run numbered differently from your
structure would attach a score to every residue, report 100% coverage, and be
wrong everywhere — with nothing downstream able to detect it.

So the amino acid ConSurf recorded is compared against the one your structure
actually contains::

   ConSurf data does not describe 2F56.pdb: only 2 of 106 compared residues
   agree (1.9%, threshold 95%). This usually means the run and the structure
   disagree about residue numbering, in which case every score attached would
   be wrong. First disagreements: A3: structure has ASN, ConSurf has VAL ...

That is a real message from the barnase study. PDB entry 2F56 numbers the
residue everyone else calls 3 as 1.

It is a **rate**, not a pass/fail, because two different things cause
disagreement:

============================  =============
case                          identity rate
============================  =============
correct numbering             100%
one point mutation            99.1%
numbering shifted by one      **2.8%**
============================  =============

A point mutant *should* disagree at the mutated position — that is correct data.
A numbering offset disagrees nearly everywhere. The 95% default sits in the wide
gap between them, and is configurable.


Reading the numbers
-------------------

Two conventions catch people out, so WatCon states them in every output:

- **SCORE: more negative = more conserved.** It is z-normalised *within each
  run* (mean ≈ 0, sd ≈ 1), so it is a relative measure, never absolute.
- **GRADE: 9 = most conserved, 1 = most variable.** The opposite sense to the
  score. Grades are per-run percentile bins, so they are ordinal — WatCon
  reports ``max_grade`` and never averages them.

.. caution::

   **Do not compare raw scores between different ConSurf runs.** Both are
   normalised within their own run. Two runs of the same barnase sequence from
   different starting structures agree at only ρ ≈ 0.37. For pooling across a
   family, use :func:`WatCon.evolutionary.conservation_by_msa_column` and prefer
   ``unanimous_conserved``, which uses each run's own verdict rather than
   comparing separately normalised numbers.


If you need many runs
---------------------

For a real protein family you need one ConSurf run per member, and the manual
route becomes tedious. Two options:

- **ConSurf-DB** (https://consurfdb.tau.ac.il/) holds precomputed results for
  many PDB entries. It was unreachable during this work, but is worth trying.
- **Standalone ConSurf** (https://github.com/Rostlab/ConSurf) runs locally. It
  needs sequence databases and setup — roughly a day — but removes the
  bottleneck permanently and makes runs reproducible.

If an official API appears, the single place to add it is
:func:`WatCon.evolutionary.find_consurf_file`, which is the only function that
decides which ConSurf file belongs to a structure. Everything downstream would
be unaffected.

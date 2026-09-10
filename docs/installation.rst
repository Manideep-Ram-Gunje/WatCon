Installation Guide
==================

Quick install (recommended)
---------------------------

.. code-block:: bash

   pip install "git+https://github.com/Manideep-Ram-Gunje/WatCon.git@consurf-integration"

.. warning::

   The ``@consurf-integration`` suffix is required. This fork's ``main`` is a
   clean mirror of upstream WatCon, which does **not** contain the ConSurf
   integration; an install URL without the branch gives you plain WatCon.

Then check it:

.. code-block:: bash

   watcon --version
   watcon demo          # full pipeline on bundled example data, offline

Every runtime dependency is installed automatically. Python >= 3.10.

.. note::

   **MODELLER is optional.** It is needed only for MSA-based structural
   alignment across a protein *family*. For many structures of one protein --
   crystal forms, mutants, complexes -- ``watcon prepare`` superposes them with
   :mod:`WatCon.superpose`, which is plain numpy and needs no licence.

   Install MODELLER only if you need family-level alignment. It is licensed,
   conda-only, and cannot be pip-installed; follow the conda instructions below.


Conda install (for MODELLER, or a fully pinned environment)
-----------------------------------------------------------

1. Clone the WatCon Repository
------------------------------

Clone the WatCon repository and change your working directory into WatCon

.. code-block:: bash

   git clone https://github.com/kamerlinlab/WatCon.git

   cd WatCon


2. Create a WatCon Conda Environment
------------------------------------

.. code-block:: bash

   conda env create -f WatCon.yaml


3. Add Modeller license
-----------------------

After you create the WatCon conda environment, you will receive this message:

.. code-block:: txt

   Edit /anaconda3/envs/WatCon/lib/modeller-10.7/modlib/modeller/config.py
   and replace XXXX with your Modeller license key
   (or set the KEY_MODELLER environment variable before running 'conda install').

This message is prompted because certain features of WatCon are dependent on the Modeller package by the Sali Lab. This package requires a license which can be obtained `here <https://salilab.org/modeller/>`_. Once you receive the license, simply edit the modeller config file with the license key.


Once completed, a WatCon environment should have been created so that you can freely use the code. With input files, call WatCon on the command line by

.. code-block:: console

   $ python -m WatCon.WatCon --input input.txt --name name_of_system

or simply import WatCon as a python package directly:

.. code-block:: python

   import WatCon.sequence_processing

   sequence_processing.pdb_to_fastas('clean_pdbs/structure1.pdb', 'fasta', 'PTP1B')


.. note::
   WatCon is compatible with Python versions >3.9 and <3.12. Attempting to use WatCon with Python versions outside of this range may result in unexpected behavior. 

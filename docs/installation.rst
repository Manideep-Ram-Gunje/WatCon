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


Installing PyMOL (needed only for the plugin)
---------------------------------------------

Everything except the interactive plugin works without PyMOL: the whole command
line, every analysis, the CSV outputs, and the ``.pml`` session files. You only
need PyMOL to *open* those sessions and to use the plugin.

.. note::

   Sessions written by ``watcon view`` and ``watcon family`` are plain text. You
   can generate them on a machine with no PyMOL and open them on one that has
   it.

**Recommended: conda-forge.** This is the reliable route on all three operating
systems, and the one CI uses.

.. code-block:: bash

   conda install -c conda-forge pymol-open-source

**Windows.** Either the conda command above, or the installer from
`pymol.org <https://pymol.org/>`_. If you use the installer, install WatCon into
*PyMOL's own* Python so the plugin can import it — see the check below.

**macOS and Linux.** The conda command above. On Linux your distribution may
also package it (``apt install pymol``), but a distribution PyMOL often cannot
see a pip-installed WatCon, which breaks the plugin.

**Does it work?** Both of these must succeed *in the same Python*:

.. code-block:: bash

   python -c "import pymol; print(pymol.__file__)"
   python -c "import WatCon; print(WatCon.__version__)"

If the first works and the second does not, PyMOL has its own Python and WatCon
is not in it. Install WatCon there:

.. code-block:: bash

   # substitute the interpreter PyMOL actually uses
   /path/to/pymol/python -m pip install "git+https://github.com/Manideep-Ram-Gunje/WatCon.git@consurf-integration"

Then install the plugin and restart PyMOL:

.. code-block:: bash

   watcon plugin --install

It appears under **Plugin → WatCon + ConSurf**. ``watcon plugin --uninstall``
removes it; ``--force`` overwrites an existing copy.


Setting up a clean machine
--------------------------

The order to do things on a computer that has none of this, with a check after
each step so you find problems where they happen.

.. list-table::
   :header-rows: 1
   :widths: 5 45 50

   * - #
     - Do
     - Check
   * - 1
     - Install Python 3.10 or newer
     - ``python --version``
   * - 2
     - Install the package (the quick install above)
     - ``watcon --version`` prints ``watcon-consurf 0.9.0``
   * - 3
     - Run the bundled demo
     - ``watcon demo`` ends with ``57 scored site(s) ...``
   * - 4
     - *Optional:* install PyMOL
     - ``python -c "import pymol"`` is silent
   * - 5
     - *Optional:* install the plugin
     - ``watcon plugin --install``, restart, the menu entry exists
   * - 6
     - *Optional:* MODELLER, for family MSA alignment
     - only if you need it; see below

Steps 1–3 need no network beyond the install itself and take about five minutes.
If step 3 works you have a functioning installation, whatever happens later.

.. warning::

   ``watcon --version`` printing ``1+unknown`` means the package was built
   without git metadata — usually from a downloaded zip rather than a clone. The
   tool works, but install from the git URL above to get a real version number.

``docs/MANUAL_TESTING.md`` walks the whole thing with expected output at every
step, including what to do when something fails.


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

.. code-block:: text

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

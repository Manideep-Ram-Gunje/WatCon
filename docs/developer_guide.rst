Developer Guide
===============

For anyone reading, changing or extending this codebase. If you only want to
*use* the tool, :doc:`getting_started` is the page you want.

Two companions to this one: ``docs/ARCHITECTURE.md`` explains how the pieces fit
together and the four rules behind most design decisions, and
``docs/CODE_MAP.md`` tells you which file owns which capability. Read
ARCHITECTURE first; it is short.


Getting set up
--------------

.. code-block:: bash

   git clone -b consurf-integration https://github.com/Manideep-Ram-Gunje/WatCon.git
   cd WatCon
   python -m pip install -e ".[test]"
   python -m pytest WatCon/tests -q

You should see **795 passed, 2 skipped**. The two skips are opt-in tests that
contact RCSB; run them with ``WATCON_NETWORK_TESTS=1``.

Python 3.10 or newer. PyMOL is optional and only needed for the plugin — see
:doc:`installation`. MODELLER is optional and only needed for MSA-based
alignment across a family; ``WatCon.superpose`` covers the single-protein case
in numpy.


Running the tests
-----------------

.. code-block:: bash

   python -m pytest WatCon/tests -q              # everything
   python -m pytest WatCon/tests/test_family.py  # one module
   python -m pytest WatCon/tests -q -k identity  # by name

**Two numbers, and why they differ.** From a repository checkout you get 795
tests. From an installed wheel — ``pytest --pyargs WatCon.tests`` — you get 569,
because ten modules need raw ConSurf bundles and test structures that a
distribution deliberately does not ship. The suite prints a warning naming those
ten rather than skipping them silently. Both numbers are correct; neither is a
failure.

CI runs fourteen jobs: twelve across Linux, macOS and Windows on Python
3.10–3.13, one that installs PyMOL from conda-forge so the plugin's 28 tests
actually execute, and one that builds a wheel and runs the shipped suite from
outside the source tree.


How the tests are written
-------------------------

Three conventions worth matching.

**Names are sentences.** ``test_a_member_whose_run_matches_nothing_is_refused``
rather than ``test_family_error_3``. A failing test should tell you what broke
without opening the file.

**Modules that exist for a defect explain it.** Where a test module was written
because something went wrong, its docstring says what went wrong and why the
test prevents a recurrence. ``test_chain_identity.py`` and
``test_family_fifteen.py`` are the clearest examples, and are a reasonable way
to learn the failure modes.

**Fixtures are real, and some are deliberately awkward.** The barnase example
has the protein as chain **L** in one entry and as half a complex in another.
5HDE's catalytic residue is a phosphocysteine the published alignment omits.
``ptp1b_ensemble_ragged.pdb`` exists *so that* the tool can be shown refusing
it. Do not tidy these up — the mess is the point.

Fixtures live in ``WatCon/data/examples/`` and ``WatCon/data/consurf/fixtures/``.
Keep them small enough for CI; the CA-only extracts reproduce the full bundles'
cross-check exactly at a twentieth of the size.


Adding things
-------------

**A CLI subcommand.** Write ``cmd_yours(args)`` in ``WatCon/cli.py``, register a
subparser in ``build_parser()``, and set ``func=cmd_yours`` as its default. Add
a dispatch test to ``tests/test_cli.py`` — ``test_every_subcommand_has_something_to_run``
will catch a subcommand wired into the parser but bound to nothing.

**A ConSurf dialect.** Extend ``detect_dialect()`` in
``WatCon/consurf/dialects.py`` and add a real grades file under
``data/consurf/fixtures/``. The parser has no network access and no guessing:
if a dialect cannot be detected it raises rather than assuming one.

**A modified residue.** Add it to ``MODIFIED_RESIDUES`` in
``WatCon/residue_index.py``. If it is the *same* amino acid in a different
protonation state — Amber's ``HID``, CHARMM's ``HSE`` — it belongs in
``PROTONATION_VARIANTS`` instead, which is normalised away before identity
comparison. A chemically modified residue such as ``CSP`` must stay
distinguishable, because reporting ``CYS`` against ``CSP`` is correct behaviour.

``docs/CODE_MAP.md`` has a fuller "if you want to change X, touch Y" table.


Things to be careful about
--------------------------

.. warning::

   **Never join on array position.** Everything is keyed on
   ``(chain, resid, icode)``. Position-based indexing is the original defect
   this fork exists to fix, and it fails silently — it returns a real value for
   the wrong residue.

.. warning::

   **Do not let the two structure readers drift apart.** WatCon reads structures
   both as text (``residues_from_pdb_file``) and through MDAnalysis
   (``atom_identity``). The ConSurf join is keyed on the chain both produce. When
   they disagreed, every lookup missed and nothing said so.

Other traps, each of which was a real defect:

* ``mol.atoms`` on a water residue returns the whole residue, including the
  alternate conformers you just filtered out.
* A ``;`` inside a ``#`` comment in a ``.pml`` file still splits the line —
  PyMOL runs the prose after it as Python.
* ``shortest_path`` and ``characteristic_path_length`` run all-pairs shortest
  paths and are on by default. They were 819 s and 45 s of a 957 s run whose
  network took 1.4 s to build.


The changelog is part of the work
---------------------------------

``docs/CONSURF_CHANGELOG.md`` records every phase in a fixed shape: what
changed, which files, why, what the tests said, what limits remain, and which
alternative was rejected. It is append-only and newest last.

This is not ceremony. Several entries exist because a later session needed to
know why something was done, and two claims have been publicly withdrawn after
re-measurement — that history is in there rather than quietly overwritten. If
you change behaviour, add an entry.


Reporting problems
------------------

Issues and pull requests for **this fork**:
https://github.com/Manideep-Ram-Gunje/WatCon/issues

For the original WatCon — the water-network analysis, clustering and PyMOL
projections, which are the original authors' work — use
https://github.com/kamerlinlab/WatCon/ instead. ``docs/CODE_MAP.md`` lists which
files are inherited and which were added here, so you can tell where a problem
belongs.

Everything is GPL-3.0, inherited from WatCon. If you use it, cite the original
paper: Brownless, Harrison-Rawn and Kamerlin, *JACS Au* 2025,
`10.1021/jacsau.5c00447 <https://pubs.acs.org/doi/10.1021/jacsau.5c00447>`_.

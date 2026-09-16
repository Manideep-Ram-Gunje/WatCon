"""The test suite's own assumptions about where it is running.

The suite ships inside the wheel, so ``pytest --pyargs WatCon.tests`` is
something a user can type. Some modules need data a distribution deliberately
does not carry -- the 5 MB of raw ConSurf bundles under
``data/consurf/exploratory`` (excluded by ``MANIFEST.in``, and asserted absent
by the CI wheel check), and the raw structures under ``tests/inputs`` and
``tests/water_dir``. Run from an installed copy those modules produced **89
collection errors and 41 failures**, which reads as a broken package rather
than as "these tests need the repository".

``conftest.py`` now declines to collect them there and says so. These tests
keep that machinery honest: the lists must name modules and directories that
really exist, or the guard would quietly stop guarding after a rename.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from . import conftest

HERE = Path(conftest.__file__).resolve().parent

#: True when the suite is running from a repository checkout, where every input
#: is present. The checks below that assert data *exists* only make sense there;
#: asserting them from an installed copy would recreate the very failures this
#: module exists to prevent.
IN_CHECKOUT = not conftest.collect_ignore
needs_checkout = pytest.mark.skipif(
    not IN_CHECKOUT, reason="running from an installed copy; repository data absent")


@pytest.mark.parametrize("attribute", [
    "_REQUIRES_MDANALYSIS",
    "_REQUIRES_EXPLORATORY",
    "_REQUIRES_INPUTS",
    "_REQUIRES_WATER_DIR",
])
def test_every_guarded_module_exists(attribute):
    """A renamed module would silently fall out of its guard."""
    for name in getattr(conftest, attribute):
        assert (HERE / name).is_file(), "%s lists a missing module: %s" % (attribute, name)


def test_no_module_is_guarded_twice_in_the_ignore_list():
    assert len(conftest.collect_ignore) == len(set(conftest.collect_ignore))


@needs_checkout
def test_the_repository_checkout_collects_everything():
    """This suite only passes in full from a checkout, so here it must be empty.

    If this fails, the run is from an installed copy -- in which case the
    modules named in the warning are genuinely not covered.
    """
    assert conftest.collect_ignore == [], (
        "not collecting: %s" % ", ".join(conftest.collect_ignore))


@needs_checkout
@pytest.mark.parametrize("directory", ["inputs", "water_dir"])
def test_the_data_directories_the_guards_name_exist(directory):
    assert (HERE / directory).is_dir()


@needs_checkout
def test_the_exploratory_bundles_are_present_in_a_checkout():
    assert conftest.EXPLORATORY.is_dir()
    assert any(conftest.EXPLORATORY.glob("*.tar.gz"))

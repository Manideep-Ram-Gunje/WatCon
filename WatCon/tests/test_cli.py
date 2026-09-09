"""The ``watcon`` command.

A console script is the first thing a new user touches, so these check the
promises the packaging makes: that the subcommands exist, that ``--help`` and
``validate`` work without the heavy scientific stack being importable, and that
the old ``python WatCon/WatCon.py`` entry point still exists.
"""

from __future__ import annotations

import os
import subprocess
import sys

import pytest

from WatCon.cli import build_parser, main

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FIXTURE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "data", "consurf", "fixtures", "1BRS_A_150.grades.txt",
)


# ===========================================================================
# Shape of the interface
# ===========================================================================

def test_every_advertised_subcommand_exists():
    parser = build_parser()
    actions = [a for a in parser._actions if hasattr(a, "choices") and a.choices]
    commands = set()
    for action in actions:
        commands.update(action.choices)
    assert {"prepare", "run", "validate", "demo"} <= commands


def test_help_exits_cleanly():
    with pytest.raises(SystemExit) as excinfo:
        main(["--help"])
    assert excinfo.value.code == 0


def test_version_is_reported():
    with pytest.raises(SystemExit) as excinfo:
        main(["--version"])
    assert excinfo.value.code == 0


def test_no_command_prints_help_and_fails():
    assert main([]) == 1


def test_prepare_requires_its_directories():
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["prepare"])


# ===========================================================================
# The commands do something
# ===========================================================================

def test_validate_reads_a_real_grades_file(capsys):
    assert main(["validate", "--consurf", FIXTURE]) == 0
    assert "1BRS" in capsys.readouterr().out


def test_validate_reports_a_bad_file_without_crashing(tmp_path, capsys):
    bad = tmp_path / "broken_consurf_grades.txt"
    bad.write_text("this is not a ConSurf grades file\n")
    code = main(["validate", "--consurf", str(bad)])
    assert code != 0
    capsys.readouterr()


def test_prepare_runs_end_to_end(tmp_path):
    examples = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "data", "examples", "barnase", "structures",
    )
    if not os.path.isdir(examples):
        pytest.skip("bundled example data not present")

    out = tmp_path / "prepared"
    report = tmp_path / "prep.csv"
    code = main([
        "prepare", "--input-dir", examples, "--out-dir", str(out),
        "--reference", "1A2P", "--report", str(report),
    ])
    assert code == 0
    assert len(os.listdir(out)) == 6
    assert report.exists()


def test_prepare_fails_cleanly_on_a_missing_directory(tmp_path, capsys):
    code = main(["prepare", "--input-dir", str(tmp_path / "nope"),
                 "--out-dir", str(tmp_path / "out")])
    assert code == 1
    assert "no such directory" in capsys.readouterr().err


def test_run_without_any_input_is_an_error(capsys):
    assert main(["run"]) == 1
    assert "give --input" in capsys.readouterr().err


# ===========================================================================
# Promises the packaging makes
# ===========================================================================

def test_help_does_not_need_the_scientific_stack():
    """`watcon --help` must work in a half-broken environment.

    Imports are deferred into each subcommand for this reason: a user whose
    MDAnalysis install is broken still needs the tool to tell them what it can
    do, and `validate` to diagnose a grades file.
    """
    code = (
        "import sys, types;"
        "blocked = ('MDAnalysis','matplotlib','sklearn','networkx','joblib','pandas');"
        "sys.meta_path.insert(0, type('B',(object,),{"
        "'find_module': lambda self, name, path=None:"
        " self if name.split('.')[0] in blocked else None,"
        "'load_module': lambda self, name: (_ for _ in ()).throw(ImportError(name))"
        "})());"
        "from WatCon.cli import build_parser;"
        "build_parser().parse_args(['--help'])"
    )
    result = subprocess.run([sys.executable, "-c", code], cwd=REPO,
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "prepare" in result.stdout


def test_the_original_entry_point_still_works():
    """`python WatCon/WatCon.py --input ...` must not have been broken."""
    result = subprocess.run(
        [sys.executable, os.path.join("WatCon", "WatCon.py"), "--help"],
        cwd=REPO, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "--input" in result.stdout


def test_the_console_script_target_is_importable():
    """pyproject declares `watcon = WatCon.cli:main`; it must resolve."""
    import importlib

    module = importlib.import_module("WatCon.cli")
    assert callable(getattr(module, "main", None))

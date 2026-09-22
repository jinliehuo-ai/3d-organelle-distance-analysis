"""Tests for the public-release check.

This file is itself published, so it must not contain any of the private terms or
private paths it exercises. Real identifiers come from the denylist named by
$PUBLIC_RELEASE_DENY_FILE; fixture strings are assembled at run time so the literals
never appear in the source.
"""
import os
from pathlib import Path

import pytest

from scripts.qa.check_public_release import load_deny_terms, main

REPO = Path(__file__).resolve().parents[1]


def test_release_tree_is_data_free():
    assert main(REPO, quiet=True) == 0


@pytest.mark.skipif(not os.environ.get("PUBLIC_RELEASE_DENY_FILE"),
                    reason="set PUBLIC_RELEASE_DENY_FILE to the private denylist kept outside this repo")
def test_release_tree_has_no_private_terms():
    """Gate a release on the project's real identifiers without naming them here."""
    assert main(REPO, deny_file=os.environ["PUBLIC_RELEASE_DENY_FILE"], quiet=True) == 0


def test_denied_term_is_detected(tmp_path, capsys):
    """A denied term in ordinary source text must fail, not only a data-like file name."""
    term = "Fixture" + "-7"
    (tmp_path / "notes.md").write_text(f"Measured in the {term} clone.\n", encoding="utf-8")
    deny = tmp_path / "denylist.txt"
    deny.write_text(term + "\n", encoding="utf-8")
    assert main(tmp_path, deny_file=deny, quiet=True) == 1
    assert term in capsys.readouterr().out


def test_denied_term_matches_case_insensitively(tmp_path):
    term = "Fixture" + "-7"
    (tmp_path / "notes.md").write_text(term.lower() + "\n", encoding="utf-8")
    deny = tmp_path / "denylist.txt"
    deny.write_text(term.upper() + "\n", encoding="utf-8")
    assert main(tmp_path, deny_file=deny, quiet=True) == 1


def test_private_path_is_detected(tmp_path):
    fake_home = "/" + "Users" + "/someone/data"
    (tmp_path / "script.py").write_text(f'BASE = "{fake_home}"\n', encoding="utf-8")
    assert main(tmp_path, quiet=True) == 1


def test_home_expression_is_detected(tmp_path):
    (tmp_path / "script.py").write_text("BASE = Path." + "home()\n", encoding="utf-8")
    assert main(tmp_path, quiet=True) == 1


def test_data_like_file_is_detected(tmp_path):
    (tmp_path / "results.csv").write_text("a,b\n1,2\n", encoding="utf-8")
    assert main(tmp_path, quiet=True) == 1


def test_clean_tree_passes(tmp_path):
    (tmp_path / "readme.md").write_text("Nothing private here.\n", encoding="utf-8")
    assert main(tmp_path, quiet=True) == 0


def test_deny_file_comments_and_blanks_ignored(tmp_path):
    deny = tmp_path / "d.txt"
    deny.write_text("# comment\n\n  TERM  # trailing\n", encoding="utf-8")
    assert load_deny_terms(deny) == ["TERM"]


def test_missing_deny_file_fails_loudly(tmp_path):
    with pytest.raises(SystemExit):
        load_deny_terms(tmp_path / "absent.txt")

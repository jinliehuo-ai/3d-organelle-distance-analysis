"""Fail if a source tree looks like it still carries data, private paths or private terms.

Three layers, because the first two alone let identifiers through in source comments
and string literals:

1. Data-like files, by extension and by size.
2. Private machine paths, including home-directory expressions.
3. Private terms supplied by --deny-file.

Layer 3 exists because the terms worth blocking - unpublished cell line or clone
names, internal batch identifiers, collaborator names - are themselves confidential.
Keep that file OUTSIDE the repository and pass it in, or set PUBLIC_RELEASE_DENY_FILE.
A denylist committed to the repository it protects publishes exactly what it hides.

Exit status is 0 when clean and 1 on any finding, so CI can gate a release on it.
"""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path


FORBIDDEN_SUFFIXES = {".nd2", ".tif", ".tiff", ".xls", ".xlsx", ".csv", ".png", ".pdf", ".svg"}
MAX_FILE_BYTES = 2_000_000
PRIVATE_PATTERNS = (
    re.compile(r"/(?:Users|Volumes|private|home)/"),
    re.compile(r"[A-Za-z]:\\"),
    re.compile(r"Path\.home\(\)"),
    re.compile(r"os\.path\.expanduser"),
    # Built from a code point so this file does not contain the literal it hunts for.
    re.compile(r"(?<![\w.])" + chr(126) + "/"),
)


def load_deny_terms(deny_file):
    """One term per line; blank lines and # comments ignored. Matched case-insensitively."""
    if deny_file is None:
        return []
    path = Path(deny_file)
    if not path.is_file():
        raise SystemExit(f"deny file not found: {path}")
    terms = []
    for line in path.read_text(encoding="utf-8").splitlines():
        term = line.split("#", 1)[0].strip()
        if term:
            terms.append(term)
    return terms


def tracked_files(root):
    """Files git would publish, or None when this is not a usable git work tree."""
    if not (root / ".git").exists():
        return None
    import subprocess
    try:
        run = subprocess.run(["git", "-C", str(root), "ls-files", "-z"],
                             capture_output=True, text=True, check=True)
    except (OSError, subprocess.CalledProcessError):
        return None
    return [root / name for name in run.stdout.split("\0") if name]


def iter_text_files(root):
    """Prefer git-tracked files: untracked build output such as __pycache__ is never published,
    and scanning it produces findings nobody can act on."""
    tracked = tracked_files(root)
    if tracked is not None:
        for path in sorted(tracked):
            if path.is_file():
                yield path
        return
    for path in sorted(root.rglob("*")):
        if not path.is_file() or ".git" in path.parts or "__pycache__" in path.parts:
            continue
        yield path


def main(root: Path, deny_file=None, quiet: bool = False) -> int:
    deny_terms = load_deny_terms(deny_file)
    deny_res = [(t, re.compile(re.escape(t), re.IGNORECASE)) for t in deny_terms]
    findings = []

    for path in iter_text_files(root):
        rel = path.relative_to(root)
        if path.suffix.lower() in FORBIDDEN_SUFFIXES:
            findings.append(f"data-like file: {rel}")
            continue
        if path.stat().st_size > MAX_FILE_BYTES:
            findings.append(f"unexpectedly large file: {rel}")
            continue
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        for pattern in PRIVATE_PATTERNS:
            for match in pattern.finditer(text):
                line = text.count("\n", 0, match.start()) + 1
                findings.append(f"private-path pattern {pattern.pattern!r}: {rel}:{line}")
                break
        for term, pattern in deny_res:
            for match in pattern.finditer(text):
                line = text.count("\n", 0, match.start()) + 1
                findings.append(f"denied term {term!r}: {rel}:{line}")
                break

    if not deny_terms and not quiet:
        print("NOTE: no --deny-file given, so private terms in source text were not checked.")
    if findings:
        print("PUBLIC RELEASE CHECK: FAIL")
        print("\n".join(findings))
        return 1
    print(f"PUBLIC RELEASE CHECK: PASS ({len(deny_terms)} denied terms checked)")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("root", type=Path)
    parser.add_argument("--deny-file", default=os.environ.get("PUBLIC_RELEASE_DENY_FILE"),
                        help="File of private terms to reject, one per line. Keep it outside "
                             "this repository. Defaults to $PUBLIC_RELEASE_DENY_FILE.")
    parser.add_argument("--quiet", action="store_true", help="Suppress the no-denylist note.")
    args = parser.parse_args()
    raise SystemExit(main(args.root, args.deny_file, args.quiet))

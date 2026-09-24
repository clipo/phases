#!/usr/bin/env python3
"""Compare the visible text of each tracked figure SVG against HEAD.

Matplotlib emits an SVG comment carrying the source string before each glyph
run, so the multiset of those comments is a faithful proxy for the text a reader
sees. Regenerating a figure from a script whose labels have drifted from the
committed artifact shows up here as a text difference, which is the failure this
guards against: commits d4f127c and 467bfc6 both corrected labels in the SVGs
without touching the code that writes them, so any regeneration silently
reverted them ("Cultural $F_{ST}$" back to "Cultural F_ST", and the
aggregated-conformity mimic back to "aggregated signaling").

Run it after regenerating figures and before committing them. A reported
difference is not automatically wrong: it is either a label regression to fix in
the generating script, or an intended change, in which case the manuscript's
wording for that label should be checked in the same pass.

Usage: .venv/bin/python scripts/check_figure_labels.py
Exit 1 if any tracked figure's visible text differs from HEAD.
"""
from __future__ import annotations

import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COMMENT = re.compile(r"<!--\s*(.*?)\s*-->", re.S)
# Matplotlib's own provenance comments, which carry a timestamp and so differ on
# every run regardless of content.
BOILERPLATE = ("Matplotlib", "Created with")


def text_of(svg: str) -> Counter:
    return Counter(
        m.group(1)
        for m in COMMENT.finditer(svg)
        if not m.group(1).startswith(BOILERPLATE)
    )


def head_version(rel: str) -> str | None:
    r = subprocess.run(
        ["git", "show", f"HEAD:{rel}"], cwd=ROOT, capture_output=True, text=True
    )
    return r.stdout if r.returncode == 0 else None


def main() -> int:
    tracked = subprocess.run(
        ["git", "ls-files", "figures/*.svg"], cwd=ROOT, capture_output=True, text=True
    ).stdout.split()
    checked = bad = 0
    for rel in sorted(tracked):
        p = ROOT / rel
        old = head_version(rel)
        if not p.exists() or old is None:
            continue
        checked += 1
        a, b = text_of(old), text_of(p.read_text(errors="replace"))
        if a == b:
            continue
        bad += 1
        print(f"\n{rel}")
        for s in sorted((a - b).keys()):
            print(f"    committed only:   {s!r}")
        for s in sorted((b - a).keys()):
            print(f"    regenerated only: {s!r}")
    print(f"\n{checked} tracked figures compared, {bad} with changed visible text")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())

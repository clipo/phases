"""90_placebo_divisions.py - is the phases' ABC posterior a property of the phase lines or of the map?

Analysis 84 gives P(copying boundary at the phase lines) = 0.79 against a prior
of 0.75. If an arbitrary division with the phases' group sizes gives the same or
a higher value, that lean is a property of dividing this map (the shorter-range
variation the copying model misses), not of the phase lines. This wrapper runs
84 with --placebo for three seeds (each replaces the phases with a same-size
division around random centers; 84 writes
`output/findings/phases_as_groups_placebo<seed>.md` and its runs CSV) and
collects the posteriors into one table.

Output: output/findings/placebo_divisions.md
Usage: .venv/bin/python analyses/90_placebo_divisions.py
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT_MD = ROOT / "output" / "findings" / "placebo_divisions.md"
SEEDS = (11, 12, 13)


def read(path: Path) -> dict:
    t = path.read_text(encoding="utf-8")
    grab = lambda pat: float(re.search(pat, t).group(1))
    return dict(ari=grab(r"adjusted Rand with the phases ([0-9.]+)") if "PLACEBO" in t else 1.0,
                boundary=grab(r"P\(some copying boundary at phase lines, leak < 1\) = \*\*([0-9.]+)\*\*"),
                local=grab(r"P\(local innovation, share > 0\) = \*\*([0-9.]+)\*\*"))


def main() -> int:
    for s in SEEDS:
        subprocess.run([sys.executable, str(ROOT / "analyses" / "84_phases_as_groups.py"), "--placebo", str(s)],
                       check=True, cwd=ROOT)
    rows = [("the published phases", read(ROOT / "output" / "findings" / "phases_as_groups.md"))]
    rows += [(f"placebo, seed {s}", read(ROOT / "output" / "findings" / f"phases_as_groups_placebo{s}.md")) for s in SEEDS]
    L = ["# Is the phases' ABC posterior a property of the phase lines or of the map?", "",
         "Produced by `analyses/90_placebo_divisions.py`, which runs `analyses/84_phases_as_groups.py --placebo <seed>` "
         "for each seed: the phases are replaced by a division with their group sizes around random centers, and the "
         "whole grid, posterior and all, is rerun on it (300 runs per cell, 5 percent acceptance, prior 0.75 on both).", "",
         "| groups | adjusted Rand with the phases | P(copying boundary at the lines) | P(local innovation) |",
         "|---|---|---|---|"]
    for name, r in rows:
        L.append(f"| {name} | {r['ari']:.3f} | {r['boundary']:.2f} | {r['local']:.2f} |")
    L += ["", "If the placebo rows match or exceed the phases, the lean toward a copying boundary is what any "
          "division of this map produces, not evidence about the phase lines.", ""]
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

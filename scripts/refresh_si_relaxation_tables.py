#!/usr/bin/env python3
"""Regenerate the four tables in the SI's "Relaxing four assumptions" section.

Those tables were typed into SUPPLEMENTAL_TEXT.md by hand from a run on the
29-assemblage basin, and they stayed there through the move to the
43-assemblage phase set: on 2026-09-21 the innovation table still showed a
boundary strength of 0.25 at a 1.1-fold shortfall with diversity matched, while
the main text, re-measured, said 1.9-fold and "still fails". Forty rows of
numbers nobody regenerates is how that happens, so the rows are now written
from the files the analyses produce.

Sources, all for the basin:
  unequal populations   output/unequal_populations.json       (64, 300 per cell)
  transport geography   output/other_departures.json, dep. A  (65, 250 per cell)
  accumulation spans    output/other_departures.json, dep. B  (65, 250 per cell)
  innovation boundary   output/other_departures_c_highrep.json (65 --reps 1200),
                        because boundary strength 0.25 sits on the diversity
                        tolerance edge and 250 draws call it the wrong way; the
                        two strongest settings, which the fine grid omits, come
                        from the default run and are labelled as such.

Each table is located by the bold heading above it and replaced whole. The
script asserts it found exactly one table under each heading before writing, so
a restructured SI fails loudly instead of being half-edited (rule 11).

Usage:
    python scripts/refresh_si_relaxation_tables.py [--check]
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SI = ROOT / "docs" / "manuscript" / "SUPPLEMENTAL_TEXT.md"
OUT = ROOT / "output"


def _row(setting, r):
    short = r["observed"] / max(r["fst_median"], 1e-9)
    return (f"| {setting} | {r['fst_median']:.4f} | [{r['fst_lo']:.4f}, {r['fst_hi']:.4f}] | "
            f"{short:.1f}x | {'matched' if r['diversity_matched'] else 'broken'} |")


def tables() -> dict[str, list[str]]:
    dep = [r for r in json.loads((OUT / "other_departures.json").read_text())
           if r["region"] == "basin"]
    hi = [r for r in json.loads((OUT / "other_departures_c_highrep.json").read_text())
          if r["region"] == "basin"]
    uneq = [r for r in json.loads((OUT / "unequal_populations.json").read_text())
            if r["region"] == "basin"]
    head = ["| setting | median $F_{ST}$ | 95 percent | shortfall | diversity |",
            "|---|---:|---|---:|:---:|"]
    t = {}
    t["**Unequal site populations**"] = (
        ["| CV | harmonic N | median $F_{ST}$ | 95 percent | shortfall | diversity |",
         "|---:|---:|---:|---|---:|:---:|"]
        + [f"| {r['cv']:.1f} | {r['harmonic_n']:.0f} | {r['fst_median']:.4f} | "
           f"[{r['fst_lo']:.4f}, {r['fst_hi']:.4f}] | "
           f"{r['observed'] / max(r['fst_median'], 1e-9):.1f}x | "
           f"{'matched' if r['diversity_matched'] else 'broken'} |" for r in uneq])
    t["**Transport geography**"] = head + [_row(r["setting"], r) for r in dep
                                           if r["departure"] == "A"]
    t["**Unequal accumulation spans**"] = head + [_row(r["setting"], r) for r in dep
                                                  if r["departure"] == "B"]
    fine = {r["setting"] for r in hi}
    extra = [r for r in dep if r["departure"] == "C" and r["setting"] not in fine
             and float(r["setting"].split()[-1]) > 0.4]
    t["**An innovation boundary**"] = (
        head + [_row(r["setting"], r) for r in hi]
        + [_row(r["setting"] + " (250 realizations)", r) for r in extra])
    return t


def main() -> int:
    check = "--check" in sys.argv
    s = SI.read_text(encoding="utf-8")
    new = s
    for heading, lines in tables().items():
        # heading paragraph, blank line, then a contiguous block of table rows
        pat = re.compile(r"(" + re.escape(heading) + r"[^\n]*\n\n)((?:\|[^\n]*\n)+)")
        found = pat.findall(new)
        assert len(found) == 1, f"expected one table under {heading}, found {len(found)}"
        new = pat.sub(lambda m: m.group(1) + "\n".join(lines) + "\n", new, count=1)
    if new == s:
        print("SI relaxation tables are current.")
        return 0
    if check:
        print("SI relaxation tables are STALE; run without --check.", file=sys.stderr)
        return 1
    SI.write_text(new, encoding="utf-8")
    print("rewrote the four relaxation tables in", SI.name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

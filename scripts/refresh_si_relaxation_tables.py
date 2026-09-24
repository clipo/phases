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


def _calibration_tables() -> dict[str, list[str]]:
    """Tables S3 and S4 from analysis 47's summary.json (added 2026-09-22)."""
    import pandas as pd
    s = json.loads((OUT / "revision_2026_09" / "summary.json").read_text())
    t: dict[str, list[str]] = {}
    label = {"pooled": "pooled profile", "uniform": "uniform"}
    rows = ["| Set | Innovation model | Matched cells (of 297; six-seed screen / fifty-seed re-check) | "
            "Selected N, innovation, mixing | Observed H_S / richness / H_T | Achieved H_S / richness / H_T |",
            "|---|---|---:|---|---|---|"]
    for c in s["calibration"]:
        if c["region"] != "basin":
            continue
        o, a = c["observed"], c["achieved"]
        flag = "" if c["matched"] else f" (unmatched; {c['selection']})"
        rows.append(f"| St. Francis basin | {label[c['model']]} | {c['n_matched']} / {c['n_matched_stage2']}{flag} | "
                    f"{c['n_ind']:,}, {c['innovation']:g}, {c['mixing']:g} | "
                    f"{o['hs']:.3f} / {o['rich']:.2f} / {o['ht']:.3f} | "
                    f"{a['hs']['median']:.3f} / {a['rich']['median']:.2f} / {a['ht']['median']:.3f} |")
    t["**Table S3.**"] = rows
    base = next(r for r in s["comparisons"] if r["region"] == "basin" and r["model"] == "pooled"
                and r["sampling"] == "time_transgressive" and r["metric"] == "spatial_fst")
    cal = next(c for c in s["calibration"] if c["region"] == "basin" and c["model"] == "pooled")
    obs = s["observed"]["basin"]
    rows = ["| Set | Case | Median [95%] | Frac. reaching obs. | H_S | Richness |", "|---|---|---|---:|---:|---:|",
            f"| Basin (observed {obs['spatial_fst']:.4f}; H_S {obs['hs']:.2f}, richness {obs['rich']:.1f}) | baseline | "
            f"{base['median']:.4f} [{base['lo']:.4f}, {base['hi']:.4f}] | {base['p_upper']:.3f} | "
            f"{cal['achieved']['hs']['median']:.2f} | {cal['achieved']['rich']['median']:.1f} |"]
    names = [("longer_burnin", "burn-in 2,400"), ("monomorphic_burnin", "monomorphic start, burn-in"),
             ("monomorphic_no_burnin", "monomorphic start, no burn-in"), ("innovation_halved", "innovation halved"),
             ("innovation_doubled", "innovation doubled"), ("mixing_halved", "mixing halved"), ("mixing_doubled", "mixing doubled")]
    sens = {r["case"]: r for r in s["sensitivity"] if r["region"] == "basin" and r["model"] == "pooled"
            and r["metric"] == "spatial_fst"}
    for key, name in names:
        r = sens[key]
        rows.append(f"| | {name} | {r['median']:.4f} [{r['lo']:.4f}, {r['hi']:.4f}] | {r['p_upper']:.3f} | "
                    f"{r['hs']:.2f} | {r['rich']:.1f} |")
    for leak in (0.5, 0.1, 0.03):
        r = next(x for x in s["boundary_grid"] if x["region"] == "basin" and x["model"] == "pooled"
                 and x["length"] == 24.0 and x["leak"] == leak and x["metric"] == "spatial_fst")
        rows.append(f"| | boundary, multiplier {leak:g} (24 km) | {r['median']:.4f} [{r['lo']:.4f}, {r['hi']:.4f}] | "
                    f"{r['p_upper']:.3f} | | |")
    t["**Table S4.**"] = rows
    return t


def main() -> int:
    check = "--check" in sys.argv
    s = SI.read_text(encoding="utf-8")
    new = s
    for heading, lines in {**tables(), **_calibration_tables()}.items():
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
    print("rewrote the four relaxation tables and Tables S3 and S4 in", SI.name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

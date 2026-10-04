#!/usr/bin/env python3
"""Cross-check the statistics drawn INSIDE each figure against its caption.

`check_figure_labels.py` answers a narrower question: does this figure still
match what its own script produces. A figure can pass that and still be wrong,
because the script and the manuscript can drift apart. figS5_dynamic did exactly
that: it passed regeneration while displaying a drift of +0.023 against a text
that said -0.093, the two differing in sign on a resolved result, and nothing
caught it until the numbers were compared by eye.

WHAT IT COMPARES, AND WHY NOT THE CAPTION ALONE. A first version diffed each
figure's numbers against its own caption. That was the wrong comparison twice
over: it drowned in noise, because captions legitimately cite model parameters a
plot never prints and plots carry axis ticks no caption repeats, and it would
have MISSED the figS6 case it was written for, whose printed drift contradicted
the supplement's body text rather than its caption.

So the test is: does every statistic a figure prints appear SOMEWHERE in the
manuscript? A figure showing a number the paper never states is the signature of
a figure left behind by a correction. figS6 printing +0.023 after the text moved
to -0.093 would fail this, which is the case that motivated it.

Axis ticks are excluded by construction: only values with three or more decimals
or an explicit percentage count, since round tick values are scaffolding and the
paper's statistics are not. A flagged value is still a worklist entry rather than
a defect, so this always exits 0.

TWO KNOWN LIMITS, so the output is not over-read. Finer axis ticks still slip
through: Figure 9's 0.0025/0.0075/0.0125 are scale, not claims. And matching a
bare value against the whole manuscript can clear a figure by coincidence, when
some other quantity happens to share the number: Figure 8's 0.041 was briefly
cleared that way by an unrelated interval bound added to Figure 11's caption.
The tool narrows where to look; it does not decide.

Usage: .venv/bin/python scripts/check_figure_claims.py [figure_stem ...]
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / "docs" / "manuscript" / "MAIN_TEXT.md"
SUPP = ROOT / "docs" / "manuscript" / "SUPPLEMENTAL_TEXT.md"
FIGS = ROOT / "figures"

COMMENT = re.compile(r"<!--\s*(.*?)\s*-->", re.S)
BOILER = ("Matplotlib", "Created with")
# A statistic, not an axis tick: two or more decimals, or a percentage.
# Three or more decimals, or an explicit percentage. Two decimals admitted
# far too many axis ticks (0.05, 0.25) to be usable.
STAT = re.compile(r"[-+\u2212]?\d+\.\d{3,}|\d+(?:\.\d+)?\s*(?:%|percent)")


def norm_num(s: str) -> str:
    s = s.replace("−", "-").replace(" ", "").rstrip("%")
    s = re.sub(r"percent$", "", s)
    try:
        return f"{float(s):+.4g}"
    except ValueError:
        return s


def figure_stats(stem: str) -> set[str]:
    svg = FIGS / f"{stem}.svg"
    if not svg.exists():
        return set()
    text = svg.read_text(errors="replace")
    labels = [m.group(1) for m in COMMENT.finditer(text)
              if not m.group(1).startswith(BOILER)]
    joined = " ".join(labels)
    return {norm_num(m.group(0)) for m in STAT.finditer(joined)}


def captions() -> dict[str, tuple[str, str]]:
    """stem -> (label, caption text), from both manuscript files."""
    out: dict[str, tuple[str, str]] = {}
    for path in (MAIN, SUPP):
        if not path.exists():
            continue
        t = path.read_text()
        for m in re.finditer(r"!\[(.*?)\]\(([^)]*?/)?([A-Za-z0-9_]+)\.(?:png|svg|pdf)\)",
                             t, re.S):
            body, stem = m.group(1), m.group(3)
            lab = re.match(r"\*\*(Figure [^.]*)\.\*\*", body)
            out[stem] = (lab.group(1) if lab else stem, body)
    return out


def manuscript_stats() -> set[str]:
    """Every statistic stated anywhere in the manuscript, both files."""
    out: set[str] = set()
    for path in (MAIN, SUPP):
        if path.exists():
            out |= {norm_num(m.group(0)) for m in STAT.finditer(path.read_text())}
    return out


def main() -> int:
    caps = captions()
    stated = manuscript_stats()
    wanted = sys.argv[1:] or sorted(caps)
    flagged = 0
    print(f"{len(caps)} captioned figures; "
          f"{len(stated)} distinct statistics stated in the manuscript\n")
    for stem in wanted:
        if stem not in caps:
            print(f"{stem}: no caption in the manuscript")
            continue
        label = caps[stem][0]
        fig = figure_stats(stem)
        if not fig:
            continue
        orphan = sorted(fig - stated)
        if orphan:
            flagged += 1
            print(f"[REVIEW] {label}  ({stem})")
            print(f"    prints {len(fig)} statistic(s); "
                  f"{len(orphan)} appear nowhere in the manuscript:")
            print(f"    {', '.join(orphan)}")
            print()
    if not flagged:
        print("Every statistic printed in a figure is also stated in the "
              "manuscript.")
    else:
        print(f"{flagged} figure(s) print a value the manuscript never states.")
        print("A value a figure shows and the paper never mentions is how a "
              "figure left behind by a correction announces itself. Check each "
              "against the text before dismissing it.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

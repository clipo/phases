#!/usr/bin/env python3
"""Check that every signed number in the manuscripts has the sign the outputs give it.

WHY THIS EXISTS. A flipped sign is the error a magnitude check cannot see. This
project has had it several times: Figure S5 printed a drift of +0.023 against a
text that said -0.093; the seriation axis's correlation with latitude was
quoted with both signs before its orientation was pinned; and on 2026-10-02 the
boundary excess went from +15.9 to -12.3 on a one-pair bin, where the question
"is the sign in the text the sign in the output" had to be asked by hand.
`check_figure_claims.py` normalizes a value before comparing and does not ask
it. This does, and nothing else.

IT COVERS HALF THE PROBLEM. This compares the text with the outputs. A sign
flipped inside a calculation reaches the outputs and the text together and
passes here. That half is covered by `tests/signatures/test_sign_conventions.py`,
which feeds each signed statistic a constructed record whose direction is
known and asserts the sign.

WHAT IT COMPARES. Every number written with an explicit sign in MAIN_TEXT.md
and SUPPLEMENTAL_TEXT.md ("+7.3", "-0.20", and the typographic minus) is looked
up in the text outputs (output/findings, output/*.md, output/revision_2026_09),
by magnitude at the precision the manuscript prints, so -0.2045 in an output
matches "-0.20" in the text.

Magnitude alone is not enough: a first version searched all outputs at once
and cleared three deliberately flipped signs, because some unrelated output
always prints 7.3 or 0.63 with the other sign. So a value is judged IN
CONTEXT. Its neighbors are the other numbers within CONTEXT_CHARS of it in the
manuscript. Every place an output prints the magnitude is scored by how many
of those neighbors appear within WINDOW lines of it (numbers with more
decimals count for more, since 0.0121 identifies a passage and 15 does not),
and the best-scoring places are taken as where the sentence got the value.
Each value is then one of:

  ok          the best-matching places print this magnitude with this sign
  FLIPPED?    they print it ONLY with the other sign
  alone       the outputs print the magnitude, but never near a neighbor, so
              its sign could not be checked in context
  unmatched   the outputs do not print the magnitude at all (a value from a
              CSV, a literature value, or a derived number)

An output value written without a sign counts as positive.

THE VERIFIED LIST. The neighbor match is a heuristic, and it flags some values
whose source is a JSON or CSV field the text outputs do not print beside the
passage's other numbers. A flag is cleared only by a person or a session
opening the source and confirming the sign, and recording that in
`scripts/check_signs_verified.tsv` (manuscript file, token, a phrase from the
same line, the source and field, the date). An entry is keyed on the token and
the phrase, so it stops applying as soon as the number or its sentence
changes, and an entry that no longer matches any flag is reported as stale.

LIMITS, so the result is not over-read. It cannot see a sign carried by words
("falls", "below", "less than") or by a column whose sign is implied; a value
quoted with no other number near it is only weakly checked; and it knows
nothing about which output a sentence means beyond the neighbors. A FLIPPED?
line is a place to look, and a clean run is not proof. It exits 1 when any
value is flagged, so the closing checks of a rerun stop on it.

Usage: .venv/bin/python scripts/check_signs.py [--show-all] [file ...]
"""
from __future__ import annotations

import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPTS = [ROOT / "docs" / "manuscript" / "MAIN_TEXT.md",
               ROOT / "docs" / "manuscript" / "SUPPLEMENTAL_TEXT.md"]
VERIFIED = ROOT / "scripts" / "check_signs_verified.tsv"
OUTPUT_GLOBS = ["output/findings/*.md", "output/findings/*.json", "output/*.md",
                "output/revision_2026_09/*.md", "output/revision_2026_09/*.json"]

MINUS = ("-", "\N{MINUS SIGN}")     # hyphen-minus and the typographic minus
SIGN_CLASS = r"[+\-\N{MINUS SIGN}]"
CONTEXT_CHARS = 200                 # how far, in the manuscript, a neighbor may sit
WINDOW = 6                          # output lines either side that count as near
# A signed number in prose: the sign must start a token, so "10-15", "k-means",
# "e-5" and the citation form "[-@key" are not read as negatives.
SIGNED = re.compile(r"(?<![\w.,)\]])(" + SIGN_CLASS + r")(\d+(?:\.\d+)?)(?![\d@]|\.\d)")
# Any number, with whatever sign character precedes it.
ANY = re.compile(r"(?<![\w.])(" + SIGN_CLASS + r"?)(\d+\.\d+|\d+)(?![\d]|\.\d)")


def _informative(tok: str) -> bool:
    """Small integers are everywhere and identify nothing."""
    return "." in tok or float(tok) >= 10


def _decimals(s: str) -> int:
    return len(s.split(".")[1]) if "." in s else 0


def _same(value: float, target: float, decimals: int) -> bool:
    """Does an output value print as `target` at `decimals` places? Half-way
    cases count either way (-0.435 is "-0.44" in one text and "-0.43" in
    another), since which way the manuscript rounded is not this check's
    business."""
    return abs(abs(value) - abs(target)) <= 0.5 * 10 ** -decimals + 1e-12


def signed_tokens(line: str):
    """[(sign '+' or '-', magnitude string, match)] for the signed numbers in a line."""
    return [("-" if m.group(1) in MINUS else "+", m.group(2), m)
            for m in SIGNED.finditer(line) if float(m.group(2)) != 0]


def load_outputs(globs=OUTPUT_GLOBS, root=ROOT):
    """[(file name, [[signed value, ...] per line])] for every output file."""
    out = []
    for pattern in globs:
        for path in sorted(root.glob(pattern)):
            try:
                text = path.read_text(errors="replace")
            except OSError:
                continue
            lines = []
            for line in text.split("\n"):
                vals = [(-float(n) if sg in MINUS else float(n)) for sg, n in ANY.findall(line)]
                if vals:
                    lines.append(vals)
            out.append((path.name, lines))
    return out


def _weight(decimals: int) -> float:
    """How much a shared number says about the source: 15 little, 0.0121 a lot."""
    return {0: 0.25, 1: 0.5, 2: 1.0, 3: 3.0}.get(decimals, 5.0)


def judge(sign: str, mag: str, neighbors: list[str], outputs):
    """Return (verdict, detail) for one signed manuscript value.

    `neighbors` are the other informative numbers near it in the manuscript,
    as (sign or None, magnitude string). A neighbor matches an output number at
    the neighbor's own precision; a neighbor the manuscript writes WITH a sign
    must also match in sign, which is what tells the passage's source from a
    table that merely shares its magnitudes; and a whole-number percentage
    also matches the fraction (13 percent ~ 0.13).
    """
    d = _decimals(mag)
    want = float(mag)
    nb = sorted({(sg or "", float(n), _decimals(n)) for sg, n in neighbors})

    def shared(values):
        vals = list(values)
        got = 0.0
        for sg, target, nd in nb:
            pool = vals if not sg else [v for v in vals if (v < 0) == (sg == "-")]
            if any(_same(v, target, nd) for v in pool) or (
                    nd == 0 and any(_same(v, target / 100, 2) for v in pool)):
                got += _weight(nd) * (2.0 if sg else 1.0)
        return got

    top, signs, where, found = 0.0, set(), set(), False
    for name, lines in outputs:
        for i, values in enumerate(lines):
            here = {("-" if v < 0 else "+") for v in values if _same(v, want, d)}
            if not here:
                continue
            found = True
            near = [v for row in lines[max(0, i - WINDOW):i + WINDOW + 1] for v in row
                    if not _same(v, want, d)]
            k = shared(near)
            if k > top:
                top, signs, where = k, set(here), {name}
            elif k == top and k > 0:
                signs |= here
                where.add(name)
    if not found:
        return "unmatched", ""
    if top <= 0:
        return "alone", ""
    if sign in signs:
        return "ok", ""
    return "FLIPPED?", (f"where the outputs print it nearest the passage's other numbers "
                        f"({', '.join(sorted(where)[:3])}), it has only the other sign")


def check_text(text: str, outputs):
    """[(verdict, line number, token, context, detail)] for every signed value in `text`."""
    results = []
    for lineno, line in enumerate(text.split("\n"), 1):
        for sign, mag, m in signed_tokens(line):
            lo, hi = max(0, m.start() - CONTEXT_CHARS), m.end() + CONTEXT_CHARS
            neighbors = [(("-" if a.group(1) in MINUS else "+") if a.group(1) else None, a.group(2))
                         for a in ANY.finditer(line[lo:hi])
                         if _informative(a.group(2)) and lo + a.start(2) != m.start(2)]
            verdict, why = judge(sign, mag, neighbors, outputs)
            results.append((verdict, lineno, m.group(0),
                            line[max(0, m.start() - 45):m.end() + 25], why))
    return results


def load_verified(path=VERIFIED):
    """[(manuscript file name, token, phrase, source)] from the verified list."""
    if not path.exists():
        return []
    rows = []
    for line in path.read_text().split("\n"):
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            raise SystemExit(f"{path.name}: expected file, token, phrase, source, date; got {line!r}")
        rows.append(tuple(x.strip() for x in parts[:4]))
    return rows


def main() -> int:
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    show = "--show-all" in sys.argv
    files = [Path(a) for a in args] or MANUSCRIPTS
    outputs = load_outputs()
    verified = load_verified()
    used = set()
    counts = defaultdict(int)
    rows = []
    for path in files:
        text_lines = path.read_text().split("\n")
        for verdict, lineno, tok, context, why in check_text("\n".join(text_lines), outputs):
            if verdict == "FLIPPED?":
                norm = tok.replace("\N{MINUS SIGN}", "-")
                for k, (vf, vt, phrase, source) in enumerate(verified):
                    if vf == path.name and vt == norm and phrase in text_lines[lineno - 1]:
                        verdict, why = "verified", source
                        used.add(k)
                        break
            counts[verdict] += 1
            rows.append((verdict, path.name, lineno, tok, context, why))
    print(f"{sum(counts.values())} signed values in {', '.join(p.name for p in files)}: "
          f"{counts['ok']} carry the sign the outputs print beside their neighbors, "
          f"{counts['verified']} flagged and cleared by the verified list, "
          f"{counts['FLIPPED?']} flagged, {counts['alone']} found only without a neighbor, "
          f"{counts['unmatched']} not in the outputs.")
    for verdict, name, lineno, tok, context, why in rows:
        if verdict == "FLIPPED?":
            print(f"FLIPPED? {name}:{lineno}  {tok}  ...{context}...")
            print(f"          {why}")
    if files == MANUSCRIPTS:
        for k, (vf, vt, phrase, source) in enumerate(verified):
            if k not in used:
                print(f"stale     {VERIFIED.name}: {vf} {vt} \"{phrase}\" no longer matches a flagged value; remove it")
    if show:
        for verdict, name, lineno, tok, context, why in rows:
            if verdict in ("alone", "unmatched", "verified"):
                print(f"{verdict:9s} {name}:{lineno}  {tok}  ...{context}..." + (f"  [{why}]" if verdict == "verified" else ""))
    elif counts["alone"] + counts["unmatched"] + counts["verified"]:
        print("(run with --show-all to list the values cleared by hand or not checkable in context)")
    if not counts["FLIPPED?"]:
        print("No unverified signed value contradicts the sign the outputs print for it.")
    return 1 if counts["FLIPPED?"] else 0


if __name__ == "__main__":
    raise SystemExit(main())

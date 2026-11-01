#!/usr/bin/env python3
"""Which files under data/ does the code actually read?

Answers the question the public-release packaging asks: ship the data the
analyses consume and nothing else. Written to be re-run, not trusted once
(rule 1); `--check` fails when the tracked set drifts from what the code reads.

Method. Every quoted filename in the Python and R sources is collected, along
with brace patterns such as f"basin_members_{which}.txt", which are globbed.
Each literal is then matched against the real tree under data/ by basename, so
a path assembled from module constants (DATA / "raw" / "x.xlsx") resolves
without having to interpret the constants themselves. Directory literals, such
as a shapefile directory handed to geopandas, pull in their contents.

The bias is deliberate and one-directional: a name that appears in a source
file counts as read even if the call sits behind a branch that never runs.
Over-inclusion leaves a spare file in the release; under-inclusion ships a
repository that cannot reproduce its own figures.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
SOURCE_DIRS = ["analyses", "src", "tests", "scripts"]
SOURCE_SUFFIXES = {".py", ".R", ".r", ".stan"}
# Extensions that denote data rather than code or prose.
DATA_SUFFIXES = {
    ".xlsx", ".xls", ".xlsm", ".csv", ".tsv", ".txt", ".14c", ".json",
    ".shp", ".shx", ".dbf", ".prj", ".cpg", ".sbn", ".sbx", ".xml",
    ".npz", ".nc", ".rds", ".geojson", ".gpkg", ".zip",
}
# A quoted string that looks like a filename, with an optional {placeholder}.
LITERAL = re.compile(r"""['"]([^'"\n]*?\.[A-Za-z0-9_]{1,8})['"]""")
DIR_LITERAL = re.compile(r"""['"]([A-Za-z0-9_][A-Za-z0-9_ ./-]*/)['"]""")


def source_files() -> list[Path]:
    out: list[Path] = []
    for d in SOURCE_DIRS:
        p = ROOT / d
        if not p.exists():
            continue
        out += [f for f in p.rglob("*") if f.suffix in SOURCE_SUFFIXES
                and "__pycache__" not in f.parts
                # This file quotes data paths to describe them, which would
                # otherwise make it a reader of everything it documents.
                and f.resolve() != Path(__file__).resolve()]
    return sorted(out)


def data_index() -> dict[str, list[Path]]:
    """Basename -> every file under data/ carrying it."""
    idx: dict[str, list[Path]] = {}
    for f in DATA.rglob("*"):
        if f.is_file():
            idx.setdefault(f.name, []).append(f)
    return idx


def resolve(literal: str, idx: dict[str, list[Path]]) -> list[Path]:
    """Files under data/ that a source literal names."""
    name = literal.split("/")[-1]
    if "{" in name:  # f-string pattern: glob the fixed parts
        pattern = re.sub(r"\{[^}]*\}", "*", name)
        hits: list[Path] = []
        for known, paths in idx.items():
            if Path(known).match(pattern):
                hits += paths
        return hits
    return idx.get(name, [])


def scan() -> tuple[dict[Path, set[str]], dict[str, set[str]]]:
    """(data file -> sources that name it, unresolved literal -> sources)."""
    idx = data_index()
    used: dict[Path, set[str]] = {}
    unresolved: dict[str, set[str]] = {}
    for src in source_files():
        text = src.read_text(encoding="utf-8", errors="replace")
        rel = str(src.relative_to(ROOT))
        for lit in set(LITERAL.findall(text)):
            if Path(lit).suffix.lower() not in DATA_SUFFIXES:
                continue
            hits = resolve(lit, idx)
            if hits:
                for h in hits:
                    used.setdefault(h, set()).add(rel)
            elif "/" in lit and lit.split("/")[0] in {"data", "../data"}:
                unresolved.setdefault(lit, set()).add(rel)
        # A directory handed to a reader (shapefile dirs, cache dirs).
        for lit in set(DIR_LITERAL.findall(text)):
            if "data/" not in lit:
                continue
            sub = lit.split("data/")[-1].rstrip("/")
            # A bare "data/" names the root, not a directory being read from.
            d = DATA / sub if sub else None
            if d is not None and d.is_dir():
                for f in d.rglob("*"):
                    if f.is_file():
                        used.setdefault(f, set()).add(rel)
    return used, unresolved


def add_shapefile_companions(used: dict[Path, set[str]]) -> None:
    """A .shp is unreadable without its sidecars, which no source names.

    geopandas opens LMVstates.shp and silently needs LMVstates.dbf/.shx/.prj
    beside it, so a literal-scan alone would ship a shapefile that cannot be
    opened. Sidecars are credited to the same readers as their .shp.
    """
    SIDECARS = {".dbf", ".shx", ".prj", ".cpg", ".sbn", ".sbx", ".qix", ".xml"}
    for shp in [f for f in list(used) if f.suffix.lower() == ".shp"]:
        for sib in shp.parent.iterdir():
            if (sib.is_file() and sib.stem.lower() == shp.stem.lower()
                    and sib.suffix.lower() in SIDECARS):
                used.setdefault(sib, set()).update(used[shp])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="exit 1 if any tracked file under data/ is unread")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    used, unresolved = scan()
    add_shapefile_companions(used)
    all_data = {f for f in DATA.rglob("*") if f.is_file()}
    unused = sorted(all_data - set(used))

    if not args.quiet:
        print(f"# Data files read by code: {len(used)} of {len(all_data)} under data/\n")
        for f in sorted(used):
            readers = ", ".join(sorted(used[f])[:3])
            more = f" (+{len(used[f]) - 3} more)" if len(used[f]) > 3 else ""
            print(f"{f.relative_to(ROOT)}\t{f.stat().st_size:>10,} B\t{readers}{more}")
        if unused:
            print(f"\n# Not read by any source file ({len(unused)}):")
            for f in unused:
                print(f"{f.relative_to(ROOT)}\t{f.stat().st_size:>10,} B")
        if unresolved:
            print(f"\n# Literals naming data/ that match no file ({len(unresolved)}):")
            for lit, srcs in sorted(unresolved.items()):
                print(f"{lit}\t{', '.join(sorted(srcs))}")
    if args.check and unused:
        print(f"\nFAIL: {len(unused)} data files are not read by any source file.",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

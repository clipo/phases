#!/usr/bin/env python3
"""Refresh the public release repo (clipo/phases) from this working tree.

The release carries the code, the data that code reads, the computed outputs,
the figures, and the manuscript. It does not carry working documents: the
review logs, the decision registers, the postmortems, CLAUDE.md or STATUS.md
stay private.

What goes in is derived, not listed by hand (rule 15). Code, outputs, figures
and the manuscript come from `git ls-files`, so an untracked scratch file is
never published by accident. Data comes from
`scripts/list_data_dependencies.py`, so the release holds the files the
analyses actually open and nothing else.

Within the managed paths the sync is authoritative: a file the release has and
this tree does not is deleted, because a stale analysis left behind in a public
repository is worse than a missing one. Everything outside them, the licenses,
citation metadata, container recipe, CI and data/README.md, belongs to the
release repo and is left untouched.

Usage:
    python scripts/sync_public_repo.py --dest ../phases [--apply]

Without --apply it prints the plan and changes nothing.
"""
from __future__ import annotations

import argparse
import filecmp
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Tracked files under these prefixes are published.
MANAGED_PREFIXES = [
    "analyses/", "src/", "tests/", "scripts/", "stan/", "R/",
    "output/", "figures/", "docs/manuscript/",
]
# Individually published tracked files at the repo root.
MANAGED_ROOT_FILES = ["pyproject.toml", "MANIFEST.md"]
# Never published, whatever git says. `sync_drive.sh` is the authors' own
# working tool: it names a private Google Drive folder and does nothing for
# reproduction, so it is a working document in the same sense as CLAUDE.md.
# `rerun_pipeline.sh` is excluded for a different reason: it reads its step
# order from `output/rerun/order.txt`, which output/* gitignores and the
# release therefore does not carry, so a reader who ran it would get a script
# that cannot start. The release's own `run_all.sh` is the reproduction path.
EXCLUDE_NAMES = {"CLAUDE.md", "STATUS.md", "sync_drive.sh", "rerun_pipeline.sh"}
# Left to the release repo: licensing, citation, container, CI, its own README.
UNMANAGED = {
    "README.md", "LICENSE", "LICENSE-data", "CITATION.cff", ".zenodo.json",
    "Dockerfile", ".dockerignore", ".gitignore", "environment.yml",
    "conda-lock.yml", "requirements.txt", "requirements-lock.txt",
    "postBuild", "run_all.sh", "CORRECTIONS.md", "data/README.md",
}


def tracked_files() -> list[str]:
    out = subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True,
                         text=True, check=True).stdout.split("\n")
    return [p for p in out if p]


def data_files() -> list[str]:
    """The data the code reads, from the dependency scanner."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "ldd", ROOT / "scripts" / "list_data_dependencies.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    used, _ = mod.scan()
    mod.add_shapefile_companions(used)
    return sorted(str(f.relative_to(ROOT)) for f in used)


def figure_renderings() -> list[str]:
    """Figure files the release needs that git here does not track.

    `figures/*.png` is gitignored in the working repo because the PNGs are
    regenerated on every build, but the manuscript markdown embeds the PNGs,
    so a release without them renders with eight broken images. The .tif
    exports are excluded: they are print-resolution siblings, tens of MB each,
    and nothing in the release references them.
    """
    out = []
    for f in sorted((ROOT / "figures").iterdir()):
        if f.is_file() and f.suffix.lower() in {".png", ".pdf", ".svg"}:
            out.append(str(f.relative_to(ROOT)))
    return out


def result_summaries() -> list[str]:
    """Computed results the release publishes, which git here does not track.

    `output/*` is gitignored in the working repo (only output/findings/** and
    the frozen baseline are tracked) because the summaries are rewritten on
    every run. The release publishes them so a reader can check a reported
    number without rerunning the pipeline. Only text is published: the .npz,
    .nc and .rds caches beside them run to hundreds of MB.
    """
    out = []
    for f in sorted((ROOT / "output").iterdir()):
        if f.is_file() and f.suffix.lower() in {".md", ".csv"}:
            out.append(str(f.relative_to(ROOT)))
    return out


def release_set() -> set[str]:
    keep = set(figure_renderings()) | set(result_summaries())
    for p in tracked_files():
        if Path(p).name in EXCLUDE_NAMES or p in UNMANAGED:
            continue
        if any(p.startswith(pre) for pre in MANAGED_PREFIXES) or p in MANAGED_ROOT_FILES:
            keep.add(p)
    keep |= set(data_files())
    keep.discard("data/processed/.gitkeep")
    return keep


def managed(path: str) -> bool:
    """Is this release path one the sync owns (and may therefore delete)?"""
    if path in UNMANAGED or path.startswith(".git/"):
        return False
    return (any(path.startswith(pre) for pre in MANAGED_PREFIXES)
            or path in MANAGED_ROOT_FILES
            or path.startswith("data/"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dest", required=True, help="path to the release clone")
    ap.add_argument("--apply", action="store_true", help="write the changes")
    args = ap.parse_args()

    dest = Path(args.dest).resolve()
    if not (dest / ".git").is_dir():
        print(f"{dest} is not a git clone", file=sys.stderr)
        return 2

    want = release_set()
    missing = [p for p in sorted(want) if not (ROOT / p).is_file()]
    if missing:
        print("source files named for release but absent:", file=sys.stderr)
        for p in missing:
            print(f"  {p}", file=sys.stderr)
        return 2

    have = {p for p in subprocess.run(["git", "ls-files"], cwd=dest,
                                      capture_output=True, text=True,
                                      check=True).stdout.split("\n") if p}

    add = sorted(want - have)
    # output/ is additive: a published result whose summary is not in this
    # working tree may still be valid, produced by a run whose artifacts were
    # cleaned. Removing one is a judgment call, so the sync lists them for
    # review instead of deleting them.
    delete = sorted(p for p in have - want
                    if managed(p) and not p.startswith("output/"))
    stale_outputs = sorted(p for p in have - want if p.startswith("output/"))
    keep_untouched = sorted(p for p in have - want if not managed(p))
    update = sorted(p for p in (want & have)
                    if not filecmp.cmp(ROOT / p, dest / p, shallow=False))

    print(f"release set: {len(want)} files")
    print(f"  add:    {len(add)}")
    print(f"  update: {len(update)}")
    print(f"  delete: {len(delete)}")
    print(f"  left alone in the release (unmanaged): {len(keep_untouched)}")
    print(f"  published outputs not regenerated here (kept, review): "
          f"{len(stale_outputs)}")
    for p in stale_outputs:
        print(f"STALE? {p}")
    for label, items in (("ADD", add), ("UPDATE", update), ("DELETE", delete)):
        for p in items:
            print(f"{label:6s} {p}")

    if not args.apply:
        print("\n(dry run; pass --apply to write)")
        return 0

    for p in add + update:
        target = dest / p
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / p, target)
    for p in delete:
        (dest / p).unlink(missing_ok=True)
    print(f"\nwrote {len(add) + len(update)} files, removed {len(delete)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

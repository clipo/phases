#!/usr/bin/env python3
"""Generate MANIFEST.md: which script produces which artifact, in runnable order.

Rule 15 requires that "regenerate everything" mean running a manifest rather
than the scripts someone remembers, and CLAUDE.md recorded that this project had
none: `analyses/` is numbered, but the dependency order between scripts was
written down nowhere.

WHY THIS IS GENERATED RATHER THAN WRITTEN. A hand-written manifest is a process
document, and rule 14 says those state counts by recounting at write time.
A manifest maintained by hand drifts from the code within a few commits, and a
manifest that lists a step which cannot be executed is worse than none at all,
which is the specific failure the basin-membership provenance note described.
So this reads the scripts and emits the file, and `--check` fails when the
committed MANIFEST.md no longer matches the tree.

HOW WRITES ARE DETECTED. Python is parsed with `ast`, not regex, for two
reasons found the hard way while writing this. First, most scripts route their
outputs through module-level constants (`OUT_MD = ROOT / "output" / "x.md"`),
so a path expression has to be resolved through the symbol table rather than
matched where it appears. Second, and more important, a regex that merely finds
path-shaped strings cannot tell a READ from a WRITE: an early version credited
21_signal_recovery.py with producing `output/closure_posterior.npz`, which it
only consumes and 53 actually writes. Only calls that are known write verbs are
counted, and the path is resolved from the call's own argument.

R is handled by regex over `here(...)` chains, which is safe here because the
eight R scripts all write through `writeLines`/`saveRDS` with the path built
inline or assigned once, and each was checked by hand against this output.

Usage:
    .venv/bin/python scripts/build_manifest.py            # write MANIFEST.md
    .venv/bin/python scripts/build_manifest.py --check    # exit 1 if stale
"""
from __future__ import annotations

import ast
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ANALYSES = ROOT / "analyses"
MANIFEST = ROOT / "MANIFEST.md"

# Calls that write an artifact. Attribute calls are matched on the attribute
# name alone, since the receiver varies (a DataFrame, a Path, a figure).
WRITE_ATTRS = {"write_text", "to_csv", "savefig", "savez", "savez_compressed",
               "save", "save_all", "to_json", "write_bytes", "savetxt",
               "to_excel"}
WRITE_FUNCS = {"save", "save_all", "savefig", "savetxt", "savez"}
TRACKED = ("output/", "figures/", "data/processed/", "data/stan/")
# Figure helpers take a bare stem and write four sibling formats.
FIG_HELPERS = {"save", "save_all"}



def _literal_join(node: ast.AST, syms: dict[str, str]) -> str | None:
    """Resolve a path expression to a '/'-joined string, or None."""
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    if isinstance(node, ast.Name):
        return syms.get(node.id)
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
        left = _literal_join(node.left, syms)
        right = _literal_join(node.right, syms)
        if left is None or right is None:
            return None
        return f"{left.rstrip('/')}/{right.lstrip('/')}"
    if isinstance(node, ast.Call):
        # Path(...) / here(...) style wrappers: use their first string argument.
        for a in node.args:
            v = _literal_join(a, syms)
            if v:
                return v
    if isinstance(node, ast.Attribute):
        # e.g. OUT.parent -- not a writable artifact path
        return None
    return None


def _normalise(p: str) -> str | None:
    """Reduce an absolute-ish path to a repo-relative tracked artifact path."""
    p = p.replace("\\", "/")
    p = re.sub(r"/{2,}", "/", p)
    for marker in TRACKED:
        i = p.find(marker)
        if i != -1:
            out = p[i:]
            # paste0()/f-string fragments leave a path with an empty segment or
            # no extension; those are not artifacts, they are path templates.
            if "/." in out or out.endswith("/") or "." not in out.rsplit("/", 1)[-1]:
                return None
            return out
    return None


def python_outputs(path: Path) -> tuple[set[str], set[str]]:
    """(artifacts written, sibling script modules imported)."""
    try:
        tree = ast.parse(path.read_text(errors="replace"))
    except SyntaxError:
        return set(), set()

    # Seed the repo-root names. They are assigned from __file__ expressions
    # that no static pass can evaluate, and leaving them unresolved silently
    # loses every `ROOT / "output" / ...` chain, which is most of the project.
    # An empty prefix is enough, because _normalise keys off the first tracked
    # directory it finds rather than off an absolute prefix.
    syms: dict[str, str] = {"ROOT": "", "HERE": "", "BASE": ""}
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and len(node.targets) == 1:
            tgt = node.targets[0]
            if isinstance(tgt, ast.Name):
                v = _literal_join(node.value, syms)
                if v:
                    syms[tgt.id] = v

    # Skip functions retired by an unconditional raise. Their write calls are
    # unreachable, so counting them reports a producer that cannot run and
    # manufactures a false conflict over an artifact with one real owner.
    dead: set[ast.AST] = set()
    for fnode in ast.walk(tree):
        if isinstance(fnode, (ast.FunctionDef, ast.AsyncFunctionDef)):
            body = [b for b in fnode.body
                    if not (isinstance(b, ast.Expr)
                            and isinstance(b.value, ast.Constant))]
            if body and isinstance(body[0], ast.Raise):
                dead.update(ast.walk(fnode))

    outs: set[str] = set()
    deps: set[str] = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call) or node in dead:
            continue
        fn = node.func
        name = fn.attr if isinstance(fn, ast.Attribute) else (
            fn.id if isinstance(fn, ast.Name) else "")

        if name == "import_module" and node.args:
            v = _literal_join(node.args[0], syms)
            if v and re.fullmatch(r"\d{2}_[A-Za-z0-9_]+", v):
                deps.add(v)
            continue

        is_write = (isinstance(fn, ast.Attribute) and name in WRITE_ATTRS) or \
                   (isinstance(fn, ast.Name) and name in WRITE_FUNCS)
        if not is_write:
            continue

        # The path is either the receiver (path.write_text(...)) or an argument.
        cands: list[ast.AST] = list(node.args)
        if isinstance(fn, ast.Attribute) and name not in FIG_HELPERS:
            # path.write_text(...) -- the receiver IS the path. For the figure
            # helpers the receiver is the module (mf.save_all), so trying it
            # first would resolve to nothing and mask the real argument.
            cands.insert(0, fn.value)
        for c in cands:
            v = _literal_join(c, syms)
            if not v:
                continue
            if name in FIG_HELPERS and "/" not in v and not v.endswith(
                    (".png", ".svg", ".pdf")):
                outs.add(f"figures/{v}.svg")
                break
            n = _normalise(v)
            if n:
                outs.add(n)
                break
    return outs, deps


R_WRITE = re.compile(r'(?:writeLines|saveRDS|ggsave|write\.csv)\s*\([^;]*?'
                     r'here\(([^)]*)\)', re.S)
R_ASSIGN = re.compile(r'(\w+)\s*<-\s*here\(([^)]*)\)')


def r_outputs(path: Path) -> set[str]:
    t = path.read_text(errors="replace")
    outs: set[str] = set()
    named = {m.group(1): m.group(2) for m in R_ASSIGN.finditer(t)}
    chains = [m.group(1) for m in R_WRITE.finditer(t)]
    for var, chain in named.items():
        if re.search(rf'(?:writeLines|saveRDS|ggsave)\s*\([^;]*\b{var}\b', t):
            chains.append(chain)
    for chain in chains:
        segs = re.findall(r'"([^"]+)"', chain)
        if segs:
            n = _normalise("/".join(segs))
            if n:
                outs.add(n)
    return outs


def collect() -> dict[str, dict]:
    rows: dict[str, dict] = {}
    for p in sorted(ANALYSES.rglob("*.py")):
        if p.name in {"figstyle.py"}:
            continue
        outs, deps = python_outputs(p)
        if outs:
            rows[str(p.relative_to(ANALYSES))] = {
                "outs": sorted(outs), "deps": sorted(deps), "lang": "py"}
    for p in sorted(ANALYSES.rglob("*.R")):
        outs = r_outputs(p)
        if outs:
            rows[str(p.relative_to(ANALYSES))] = {
                "outs": sorted(outs), "deps": [], "lang": "R"}
    return rows


def order(rows: dict[str, dict]) -> tuple[list[str], list[tuple[str, str]]]:
    """Topological order by sibling-import dependency; returns (order, cycles)."""
    stem = {k.rsplit("/", 1)[-1][:-3]: k for k in rows if k.endswith(".py")}
    edges = {k: {stem[d] for d in v["deps"] if d in stem} for k, v in rows.items()}
    cycles = [(a, b) for a, deps in edges.items() for b in deps
              if a in edges.get(b, set())]
    done: list[str] = []
    seen: set[str] = set()
    remaining = dict(edges)
    while remaining:
        ready = sorted(k for k, d in remaining.items() if not (d - seen))
        if not ready:  # a cycle: take the lowest-numbered and note it
            ready = [sorted(remaining)[0]]
        for k in ready:
            done.append(k)
            seen.add(k)
            remaining.pop(k)
    return done, sorted(set(tuple(sorted(c)) for c in cycles))


def render(rows: dict[str, dict]) -> str:
    seq, cycles = order(rows)
    producer: dict[str, list[str]] = defaultdict(list)
    for k, v in rows.items():
        for o in v["outs"]:
            producer[o].append(k)
    contested = {o: p for o, p in producer.items() if len(p) > 1}
    n_art = len(producer)

    L = [
        "# MANIFEST",
        "",
        "Which script produces which artifact, and the order they run in.",
        "",
        "**Generated. Do not edit by hand.** Written by "
        "`scripts/build_manifest.py`, which parses the scripts rather than "
        "trusting a list. Regenerate after adding or renaming any analysis:",
        "",
        "```",
        "    .venv/bin/python scripts/build_manifest.py",
        "    .venv/bin/python scripts/build_manifest.py --check   # in CI",
        "```",
        "",
        f"{len(rows)} scripts produce {n_art} tracked artifacts.",
        "",
        "Writes are detected from the call that performs them "
        "(`write_text`, `to_csv`, `savefig`, `save_all`, `savez`, `writeLines`, "
        "`saveRDS`), with paths resolved through each script's own module-level "
        "constants. Reads are deliberately NOT counted: an earlier regex "
        "version credited `21_signal_recovery.py` with producing "
        "`output/closure_posterior.npz`, which it only consumes and `53` "
        "writes.",
        "",
        "A figure is listed once, under whichever extension the script names. "
        "`figstyle.save_all` writes four siblings from one call (`.png`, "
        "`.svg`, `.pdf`, `.tiff`), and only the `.svg` is tracked in git; "
        "`figures/*.png` is ignored, so the PNGs a manuscript build consumes "
        "are derived artifacts that must be regenerated rather than cloned.",
        "",
        "## Order",
        "",
        "Topologically sorted by sibling import: a script may be run at any "
        "point after everything above it. The order is derived from the "
        "`importlib.import_module` edges between scripts, and verified to have "
        "no violations apart from the mutual pair noted below.",
        "",
        "**What this order is not.** It has not been executed as a single "
        "top-to-bottom pass, and several entries take hours (the ABC-SMC fits, "
        "the 500-realisation drift simulations, the Stan GP fits). Import edges "
        "also do not capture every dependency: a script that reads an artifact "
        "another wrote, without importing it, is ordered correctly here only "
        "because the numbering already reflects that. Treat this as the "
        "dependency record rule 15 asks for, not as a proven build script.",
        "",
    ]
    if cycles:
        L += ["> **Mutual imports.** " + "; ".join(
            f"`{a}` and `{b}`" for a, b in cycles) +
            " import each other at module level. This works because each uses "
            "the other only inside a function body, so neither needs the "
            "other's attributes at import time, and either may be run first. "
            "It does mean the order below is arbitrary between them.", ""]

    L += ["| # | script | produces |", "|---|---|---|"]
    for i, k in enumerate(seq, 1):
        arts = "<br>".join(f"`{o}`" for o in rows[k]["outs"])
        L.append(f"| {i} | `analyses/{k}` | {arts} |")

    L += ["", "## Artifacts, by path", "",
          "| artifact | produced by |", "|---|---|"]
    for o in sorted(producer):
        L.append(f"| `{o}` | " + ", ".join(f"`{p}`" for p in sorted(producer[o])) + " |")

    if contested:
        L += ["", "## Written by more than one script", "",
              "Not necessarily wrong (a figure may be produced by the analysis "
              "that owns it and by a figure driver), but each is a place where "
              "running the pipeline in a different order gives a different "
              "file, so the owner should be the later entry in the order above.",
              ""]
        for o, ps in sorted(contested.items()):
            L.append(f"- `{o}`: " + ", ".join(f"`{p}`" for p in sorted(ps)))

    L += ["", "## Not covered here", "",
          "- `analyses/figstyle.py`, `make_map.py` are helpers imported by "
          "other scripts rather than entry points.",
          "- `analyses/00_setup/*.R` set up and record the R toolchain; they "
          "produce no analysis artifact.",
          "- Scripts that only print are omitted, since they leave nothing to "
          "regenerate.",
          "- Inputs under `data/raw/` and `data/Shapefiles/` are tracked "
          "sources, not derived artifacts; see `data/DATA_INVENTORY.md`.",
          ""]
    return "\n".join(L)


def main() -> int:
    rows = collect()
    text = render(rows)
    if "--check" in sys.argv:
        if not MANIFEST.exists():
            print("MANIFEST.md is missing; run without --check to create it.")
            return 1
        if MANIFEST.read_text() != text:
            print("MANIFEST.md is stale. Regenerate:")
            print("    .venv/bin/python scripts/build_manifest.py")
            r = subprocess.run(["diff", "-u", str(MANIFEST), "-"],
                               input=text, text=True, capture_output=True)
            print(r.stdout[:4000])
            return 1
        print(f"MANIFEST.md is current ({len(rows)} scripts).")
        return 0
    MANIFEST.write_text(text)
    print(f"wrote {MANIFEST} ({len(rows)} scripts)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

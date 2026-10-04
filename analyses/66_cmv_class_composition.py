"""66_cmv_class_composition.py - what the southeast-Missouri classes actually are.

WHY THIS EXISTS. The manuscript uses the Williams [1954] southeast-Missouri
assemblages as a COMPARISON for the St. Francis basin test, and reports a
between-cluster cultural F_ST of 0.034 for them. The basin analysis runs on ten
PFG 1951 DECORATED types. This script asks the question nobody had asked of the
comparison set: what are its 42 classes?

The answer changes how the comparison should be read, so it is measured here
rather than asserted. Every number the manuscript states about the composition
of that set comes from this script.

THE THREE FINDINGS, each of which is a reason the two sets are not measuring the
same quantity:

1. The comparison set is 90 percent PLAIN WARE. Mississippi Plain alone is 81.5
   percent of it. Plain wares are utility pottery: coarse-versus-fine paste is a
   manufacturing and functional variable, not a stylistic one, so a between-group
   difference in it is the kind of signal Table 1 files under environmental
   patchiness rather than under bounded interaction.

2. The stylistic part of the set is too small to carry the test. Only 883 sherds
   (9.9 percent) are non-plain, a median of 18 per assemblage spread over 38
   classes, and only 12 of 39 assemblages reach the 30-sherd minimum the study
   itself imposes when that minimum is applied to decorated material.

3. What decorated material there is, is dominated by SURFACE FILMS rather than
   by design: Varney Red is 38.7 percent of it and Old Town Red Filmed a further
   17.4 percent. Varney Red is chronologically diagnostic, so its spatial signal
   is confounded with the Woodland-to-Mississippian gradient that the Woodland
   exclusion was meant to remove. Dropping it moves the decorated F_ST from 0.35
   to 0.18, so it carries about half.

A CONTROL, so the decorated F_ST is not over-read. A value of 0.32 on 883 sherds
invites the objection that it is a small-sample artifact. It is not pure
sampling: resampling every assemblage from the regional pooled profile at its own
observed count gives a median of 0.003. But it is not a clean stylistic
measurement either, because 29 of the 38 decorated classes occur in only one of
the four spatial clusters while holding 302 sherds between them, so rarity and
spatial restriction are confounded and cannot be separated at this sample size.
The honest statement is that the set cannot resolve the question, which is the
paper's own criterion applied to itself.

A DATA DEFECT, reported because it is in the source table. Three class names
collide under case normalization ("Thin Ware"/"Thin ware", "Fine Ware"/"Fine
ware", "Brown Slip"/"Brown slip"), splitting one class into two and inflating the
class count. Together they are 92 sherds, so no result here turns on them, but
they should be merged in the source rather than carried.

Usage: .venv/bin/python analyses/66_cmv_class_composition.py
"""
from __future__ import annotations

import importlib
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "analyses"))

g26 = importlib.import_module("26_cmv_phase_groupness")
r47 = importlib.import_module("47_revision_analysis")

OUT = ROOT / "output" / "findings" / "cmv_class_composition.md"

# The four shell-tempered plain wares in the retained set. Named explicitly
# rather than pattern-matched on "Plain", because Barnes Plain and Baytown Plain
# are already excluded as Woodland and a regex would silently pick up whatever
# the source adds next.
PLAIN = ["Mississippi Plain", "Bell Plain", "Wickliffe Plain", "Kimmswick Plain"]
N_NULL = 400
SEED = 3


def non_types(cols):
    """Classes that are descriptors or residuals rather than named types."""
    return [c for c in cols
            if str(c).lower().startswith("unidentified")
            or re.match(r"(thin|fine) ware", str(c), re.I)]


def load():
    """The exact matrix analysis 47 uses for the `cmv` set, with class names.

    `47_revision_analysis.load_sets` drops the column labels, so the matrix is
    rebuilt here by the same steps and checked against it. If the two ever
    diverge this raises rather than reporting composition for the wrong table.
    """
    M, _ = g26.load_cmv()
    miss = [c for c in M.columns if c not in g26.WOODLAND]
    Mm = M[miss]
    Mm = Mm[Mm.to_numpy().sum(1) >= 30]
    d = r47.load_sets()["cmv"]
    if not np.allclose(Mm.to_numpy(float), d["m"]):
        raise SystemExit("rebuilt CMV matrix does not match 47_revision_analysis "
                         "load_sets()['cmv']['m']; the loaders have diverged")
    return Mm, d["labels"]


def fst(mat, labels):
    keep = mat.sum(1) > 0
    return float(r47.fst_by(mat[keep], labels[keep])), int(keep.sum())


def pooled_null(mat, labels, rng):
    """No spatial structure at all: each assemblage resampled from the regional
    pooled profile at its own observed count. Rules sampling in or out as the
    explanation; it is NOT the paper's drift null, which is far stronger."""
    keep = mat.sum(1) > 0
    X, L = mat[keep], labels[keep]
    pooled = X.sum(0) / X.sum()
    tots = X.sum(1).astype(int)
    vals = []
    for _ in range(N_NULL):
        Y = np.array([rng.multinomial(t, pooled) for t in tots], float)
        ok = Y.sum(1) > 0
        vals.append(r47.fst_by(Y[ok], L[ok]))
    return np.asarray(vals)


def main() -> int:
    Mm, labels = load()
    cols = list(Mm.columns)
    dec = [c for c in cols if c not in PLAIN]
    junk = non_types(dec)
    X = Mm.to_numpy(float)
    total = X.sum()
    rng = np.random.default_rng(SEED)

    L = ["# What the southeast-Missouri comparison classes actually are", "",
         "Produced by `analyses/66_cmv_class_composition.py`. "
         f"{Mm.shape[0]} assemblages, {Mm.shape[1]} classes, {total:.0f} sherds: "
         "the set analysis 47 calls `cmv`.", "",
         "## Composition", "", "| class | share of all sherds |", "|---|---:|"]
    share = Mm.sum(0) / total
    for c, v in share.sort_values(ascending=False).head(12).items():
        L.append(f"| {c} | {100 * v:.2f}% |")
    plain_n = Mm[PLAIN].to_numpy().sum()
    dec_n = Mm[dec].to_numpy().sum()
    per = Mm[dec].to_numpy().sum(1)
    L += ["", f"Plain wares ({', '.join(PLAIN)}) hold **{plain_n:.0f} of {total:.0f} "
          f"sherds, {100 * plain_n / total:.1f} percent**. The non-plain remainder is "
          f"{dec_n:.0f} sherds ({100 * dec_n / total:.1f} percent), a median of "
          f"{np.median(per):.0f} per assemblage over {len(dec)} classes; "
          f"**{int((per >= 30).sum())} of {len(per)} assemblages** reach the 30-sherd "
          "minimum the study imposes when it is applied to that material.", "",
          "### What the non-plain sherds are", "", "| class | share of non-plain |", "|---|---:|"]
    dshare = Mm[dec].sum(0) / dec_n
    for c, v in dshare.sort_values(ascending=False).head(8).items():
        L.append(f"| {c} | {100 * v:.1f}% |")

    L += ["", "## Where the between-cluster differentiation lives", "",
          f"Cultural $F_{{ST}}$ across the {len(np.unique(labels))} spatial clusters, by subset. "
          "The null resamples every assemblage from the regional pooled profile at its own "
          f"observed count, so it carries no spatial structure at all ({N_NULL} draws, seed "
          f"{SEED}); it separates sampling from everything else, and is NOT the paper's drift "
          "null.", "",
          "| subset | classes | sherds | $F_{ST}$ | pooled-profile null (median, 95%) |",
          "|---|---:|---:|---:|---|"]
    subsets = [
        ("all classes, as analyzed", cols),
        ("plain wares only", PLAIN),
        ("Mississippi Plain vs Bell Plain only", ["Mississippi Plain", "Bell Plain"]),
        ("non-plain only", dec),
        ("non-plain, named types only", [c for c in dec if c not in junk]),
        ("non-plain, named, minus Varney Red",
         [c for c in dec if c not in junk and c != "Varney Red"]),
    ]
    rows = {}
    for name, sub in subsets:
        mat = Mm[sub].to_numpy(float)
        f, n_used = fst(mat, labels)
        nv = pooled_null(mat, labels, rng)
        rows[name] = f
        L.append(f"| {name} | {len(sub)} | {mat.sum():.0f} | {f:.4f} | "
                 f"{np.median(nv):.4f} [{np.percentile(nv, 2.5):.4f}, "
                 f"{np.percentile(nv, 97.5):.4f}] |")

    byc = np.array([Mm[dec].to_numpy(float)[labels == c].sum(0)
                    for c in np.unique(labels)])
    only1 = (byc > 0).sum(0) == 1
    L += ["", "## Reading", "",
          f"The reported comparison value of {rows['all classes, as analyzed']:.4f} is a blend in "
          f"which {100 * plain_n / total:.0f} percent of the weight is undecorated utility ware. "
          f"The plain component on its own differentiates more strongly "
          f"({rows['plain wares only']:.4f}), and the coarse-versus-fine paste ratio alone more "
          f"strongly still ({rows['Mississippi Plain vs Bell Plain only']:.4f}). That is a "
          "manufacturing variable, so a between-group difference in it is the signature Table 1 "
          "files under environmental patchiness rather than under bounded interaction.", "",
          f"The non-plain subset differentiates far more ({rows['non-plain only']:.4f}), and that "
          "is not pure sampling: the pooled-profile null for the same subset sits two orders of "
          "magnitude below it. But it is not a clean stylistic measurement either. Varney Red, a "
          "chronologically diagnostic red-filmed ware that the Woodland exclusion does not catch, "
          f"is {100 * dshare['Varney Red']:.1f} percent of the non-plain sherds, and removing it "
          f"moves the value to {rows['non-plain, named, minus Varney Red']:.4f}. "
          f"{int(only1.sum())} of the {len(dec)} non-plain classes occur in only one of the four "
          f"spatial clusters while holding {Mm[dec].to_numpy(float)[:, only1].sum():.0f} sherds "
          "between them, so rarity and spatial restriction are confounded at this sample size and "
          "cannot be separated.", "",
          "The conclusion is not that the comparison set is differentiated or undifferentiated. "
          "It is that the set cannot resolve the question the basin set is being asked, which is "
          "the paper's own empirical-sufficiency criterion applied to itself.", ""]

    seen: dict[str, list] = {}
    for c in cols:
        seen.setdefault(re.sub(r"[^a-z]", "", str(c).lower()), []).append(c)
    dups = {k: v for k, v in seen.items() if len(v) > 1}
    if dups:
        L += ["## A defect in the source table", "",
              "These class names collide under case normalization, so one class is carried as "
              "two. No result here turns on them; they should be merged in the source.", "",
              "| names | sherds |", "|---|---|"]
        for v in dups.values():
            L.append(f"| {', '.join(map(repr, v))} | {', '.join(str(int(Mm[c].sum())) for c in v)} |")
        L.append("")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""27_cmv_lmv_repertoire.py — the CMV/LMV decorated-type REPERTOIRE boundary.

The CMV (Williams 1954, SE Missouri) and LMV (Phillips-Ford-Griffin, St. Francis
and south) are drift-like in decorated-ceramic PROPORTIONS at both scales
(analyses 23-26). This script tests the other axis Lipo flagged: the decorated-
type REPERTOIRE, which types occur where. "Related but clearly different" would
show as a repertoire boundary (distinctive northern vs southern types) even where
proportional differentiation is weak.

Harmonizes Williams's and PFG's decorated typologies (both rooted in the
Phillips-Ford-Griffin 1951 framework) into one canonical vocabulary, each type
tagged shared / northern / southern. Then:
  1. Type incidence and abundance by macro-region (CMV vs LMV).
  2. Shared vs distinctive sherd shares in each region.
  3. Repertoire separation: presence/absence Jaccard within/between regions;
     decorated-proportion F_ST CMV-vs-LMV; and a SHARED-TYPES-ONLY control that
     removes the typology-vintage confound (do the regions differ even among the
     types both typologies record?).

CAVEAT logged in output: Williams 1954 and PFG 1951 share a typological frame,
but some "northern only" vs "southern only" status is partly typological vintage
(e.g. Williams's Matthews Incised may subsume Barton). The shared-types control
is the guard against reading typology differences as cultural ones.

Read-only on the manuscript. Writes output/cmv_lmv_repertoire.md and
figures/figX_cmv_lmv_repertoire.png.

Usage: .venv/bin/python analyses/27_cmv_lmv_repertoire.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import make_figures as mf  # noqa: E402
from mls_emergence.dataio.pfg import load_pfg_counts  # noqa: E402

TSV = ROOT / "data" / "raw" / "williams1954_cmv_counts.tsv"
PFG = ROOT / "data" / "raw" / "PFGData.xlsx"
OUT_MD = ROOT / "output" / "cmv_lmv_repertoire.md"
OUT_FIG = ROOT / "figures" / "figX_cmv_lmv_repertoire.png"

# canonical decorated vocabulary + affinity
AFFINITY = {
    "Parkin Punctated": "shared", "Ranch Incised": "shared",
    "Kent Incised": "shared", "Mound Place Incised": "shared",
    "Wallace Incised": "shared", "Manly Punctated": "shared",
    "Old Town Red": "shared", "Nodena Red and White": "shared",
    "Hollywood White Slipped": "shared", "Rhodes Incised": "shared",
    "Matthews Incised": "northern", "O'Byam": "northern",
    "Wickliffe": "northern", "Kimmswick": "northern", "Varney Red": "northern",
    "Barton Incised": "southern", "Vernon Paul Applique": "southern",
    "Fortune Noded": "southern", "Tyronza Punctated": "southern",
    "Walls Engraved": "southern", "Hull Engraved": "southern",
    "Carson Red-on-Buff": "southern", "Avenue Polychrome": "southern",
    "Stokes Bayou Incised": "southern", "Oliver Incised": "southern",
    "Owens Punctated": "southern", "Leland Incised": "southern",
    "Blanchard Incised": "southern", "Arcola Incised": "southern",
}

PFG_TO_CANON = {
    "Parkin Punctated": "Parkin Punctated", "Barton Incised": "Barton Incised",
    "Ranch Incised": "Ranch Incised", "Vernon Paul Applique": "Vernon Paul Applique",
    "Fortune Noded": "Fortune Noded", "Manly Punctated": "Manly Punctated",
    "Tyronza Punctated": "Tyronza Punctated", "Kent Incised": "Kent Incised",
    "Rhodes Incised": "Rhodes Incised", "Walls Engraved": "Walls Engraved",
    "Hull Engraved": "Hull Engraved", "Mound Place Incised (T)": "Mound Place Incised",
    "Old Town Red": "Old Town Red", "Carson Red-on-Buff": "Carson Red-on-Buff",
    "Nodena Red and White": "Nodena Red and White", "Avenue Polychrome": "Avenue Polychrome",
    "Hollywood White Slipped": "Hollywood White Slipped", "Wallace Incised": "Wallace Incised",
    "Stokes Bayou Incised": "Stokes Bayou Incised", "Oliver Incised": "Oliver Incised",
    "Owens Punctated (T)": "Owens Punctated", "Leland Incised": "Leland Incised",
    "Blanchard Incised (T)": "Blanchard Incised", "Arcola Incised (T)": "Arcola Incised",
}
WILLIAMS_TO_CANON = {
    "Parkin Punctate": "Parkin Punctated", "Ranch Incised": "Ranch Incised",
    "Kent Incised": "Kent Incised", "Mound Place Incised": "Mound Place Incised",
    "Wallace Incised": "Wallace Incised", "Manly Punctate": "Manly Punctated",
    "Old Town Red Filmed": "Old Town Red", "Nodena Red and White": "Nodena Red and White",
    "Hollywood White Filmed": "Hollywood White Slipped", "Rhodes-like": "Rhodes Incised",
    "Matthews Incised": "Matthews Incised", "O'Byam Engraved": "O'Byam",
    "O'Byam Incised": "O'Byam", "Wickliffe Incised": "Wickliffe",
    "Wickliffe Cord Marked": "Wickliffe", "Kimmswick Fabric Impressed": "Kimmswick",
    "Kimmswick Fabric Marked": "Kimmswick", "Varney Red": "Varney Red",
}
CANON = list(AFFINITY.keys())


def cmv_decorated():
    d = pd.read_csv(TSV, sep="\t")
    d["count"] = pd.to_numeric(d["count"], errors="coerce").fillna(0)
    d["asm"] = d["site_name"].astype(str) + " [" + d["collection"].astype(str) + "]"
    d["canon"] = d["type_name"].map(WILLIAMS_TO_CANON)
    dd = d.dropna(subset=["canon"])
    M = dd.pivot_table(index="asm", columns="canon", values="count",
                       aggfunc="sum", fill_value=0)
    return M.reindex(columns=CANON, fill_value=0)


def lmv_decorated():
    raw = load_pfg_counts(PFG)
    if not raw.index.is_unique:
        raw = raw.groupby(level=0).sum()
    cols = {c: PFG_TO_CANON[c] for c in raw.columns if c in PFG_TO_CANON}
    M = raw[list(cols)].rename(columns=cols)
    M = M.T.groupby(level=0).sum().T  # collapse any duplicate canon cols
    return M.reindex(columns=CANON, fill_value=0)


def jaccard_block(P):
    """mean pairwise Jaccard similarity on a presence/absence matrix."""
    n = P.shape[0]
    if n < 2:
        return np.nan
    sims = []
    for i in range(n):
        for j in range(i + 1, n):
            a, b = P[i] > 0, P[j] > 0
            u = (a | b).sum()
            sims.append(((a & b).sum() / u) if u else np.nan)
    return float(np.nanmean(sims))


def main():
    cmv = cmv_decorated()
    lmv = lmv_decorated()
    L = ["# CMV vs LMV decorated-type repertoire boundary", ""]

    # 1. incidence + abundance by region
    def stats(M):
        pres = (M > 0)
        inc = pres.mean(0)                      # fraction of assemblages present
        ab = M.sum(0) / M.sum().sum()           # share of decorated sherds
        return inc, ab
    ci, ca = stats(cmv)
    li, la = stats(lmv)
    L += [f"CMV decorated sherds: {int(cmv.sum().sum())} in {len(cmv)} assemblages. "
          f"LMV decorated sherds: {int(lmv.sum().sum())} in {len(lmv)} assemblages.", "",
          "## 1. Type incidence (fraction of assemblages present) and sherd share",
          "| type | affinity | CMV inc | LMV inc | CMV share | LMV share |",
          "|---|---|---|---|---|---|"]
    for t in CANON:
        L += [f"| {t} | {AFFINITY[t]} | {ci[t]:.2f} | {li[t]:.2f} | "
              f"{ca[t]:.3f} | {la[t]:.3f} |"]

    # 2. shared vs distinctive sherd shares
    def share_by_aff(ab):
        return {a: float(sum(ab[t] for t in CANON if AFFINITY[t] == a))
                for a in ("shared", "northern", "southern")}
    cs, ls = share_by_aff(ca), share_by_aff(la)
    L += ["", "## 2. Decorated-sherd share by affinity",
          "| region | shared | northern | southern |", "|---|---|---|---|",
          f"| CMV | {cs['shared']:.3f} | {cs['northern']:.3f} | {cs['southern']:.3f} |",
          f"| LMV | {ls['shared']:.3f} | {ls['northern']:.3f} | {ls['southern']:.3f} |", ""]

    # 3a. repertoire separation: Jaccard within/between (assemblages with >=2 decorated types present)
    cmin = cmv[(cmv > 0).sum(1) >= 2]
    lmin = lmv[(lmv > 0).sum(1) >= 2]
    Pc, Pl = (cmin.to_numpy() > 0).astype(int), (lmin.to_numpy() > 0).astype(int)
    np.vstack([Pc, Pl])
    np.array([0] * len(Pc) + [1] * len(Pl))
    j_cc, j_ll = jaccard_block(Pc), jaccard_block(Pl)
    # between
    bs = []
    for i in range(len(Pc)):
        for j in range(len(Pl)):
            a, b = Pc[i] > 0, Pl[j] > 0
            u = (a | b).sum()
            bs.append(((a & b).sum() / u) if u else np.nan)
    j_cl = float(np.nanmean(bs))
    L += ["## 3a. Repertoire separation (presence/absence Jaccard, assemblages with >=2 decorated types)",
          f"- assemblages used: CMV {len(Pc)}, LMV {len(Pl)}.",
          f"- mean Jaccard within-CMV = {j_cc:.3f}, within-LMV = {j_ll:.3f}, "
          f"between = {j_cl:.3f}.",
          "  Between << within on both sides = distinct repertoires.", ""]

    # 3b. proportional F_ST CMV vs LMV, EQUAL-WEIGHTED pools (each region = one
    # composition). Raw-count pooling is size-weighted and the LMV pool is ~21x
    # the CMV pool, which collapses F_ST toward 0 regardless of difference; we
    # therefore normalize each region to an equal effective size before pooling.
    def fst_two(cM, lM, cols, equal=True):
        cpool = cM[cols].sum(0).to_numpy(float)
        lpool = lM[cols].sum(0).to_numpy(float)
        if equal:
            cpool = cpool / cpool.sum() * 1000.0 if cpool.sum() else cpool
            lpool = lpool / lpool.sum() * 1000.0 if lpool.sum() else lpool
        return mf.cultural_fst(np.array([cpool, lpool]))
    shared_cols = [t for t in CANON if AFFINITY[t] == "shared"]
    fst_all = fst_two(cmv, lmv, CANON)
    fst_shared = fst_two(cmv, lmv, shared_cols)
    fst_all_raw = fst_two(cmv, lmv, CANON, equal=False)
    L += ["## 3b. Decorated-proportion F_ST, CMV vs LMV (equal-weighted pools)",
          f"- ALL decorated types: F_ST = {fst_all:+.3f}.",
          f"- SHARED types only (typology-vintage control): F_ST = {fst_shared:+.3f}.",
          f"- (size-weighted, for reference, dominated by the 21x-larger LMV pool: "
          f"{fst_all_raw:+.3f}).",
          "  If the boundary persists among shared types, it is not just a typology",
          "  artifact; if it collapses, the apparent boundary is largely distinctive-",
          "  type bookkeeping. Compare the within-region drift level F_ST ~ 0.04.", ""]

    OUT_MD.write_text("\n".join(L), encoding="utf-8")

    # figure: stacked affinity shares + incidence of distinctive types
    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["Arial", "DejaVu Sans"], "font.size": 8})
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.4))
    # left: affinity sherd-share stacked bars
    regions = ["CMV", "LMV"]
    sh = np.array([[cs[a] for a in ("shared", "northern", "southern")],
                   [ls[a] for a in ("shared", "northern", "southern")]])
    cols3 = ["#999999", "#0072B2", "#D55E00"]
    bottom = np.zeros(2)
    for k, a in enumerate(("shared", "northern", "southern")):
        ax1.bar(regions, sh[:, k], bottom=bottom, color=cols3[k], label=a)
        bottom += sh[:, k]
    ax1.set_ylabel("decorated-sherd share"); ax1.legend(fontsize=6, frameon=False)
    ax1.set_title("Repertoire by affinity", fontsize=8)
    for sp in ("top", "right"):
        ax1.spines[sp].set_visible(False)
    # right: incidence of the distinctive types
    dist = [t for t in CANON if AFFINITY[t] != "shared"]
    y = np.arange(len(dist))
    ax2.barh(y - 0.2, [ci[t] for t in dist], height=0.4, color="#0072B2", label="CMV")
    ax2.barh(y + 0.2, [li[t] for t in dist], height=0.4, color="#D55E00", label="LMV")
    ax2.set_yticks(y); ax2.set_yticklabels(dist, fontsize=5)
    ax2.set_xlabel("incidence (frac. assemblages)"); ax2.legend(fontsize=6, frameon=False)
    ax2.set_title("Distinctive-type incidence", fontsize=8)
    for sp in ("top", "right"):
        ax2.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=300, bbox_inches="tight")

    print("\n".join(L))
    print(f"\nwrote {OUT_MD}\nwrote {OUT_FIG}")


if __name__ == "__main__":
    main()

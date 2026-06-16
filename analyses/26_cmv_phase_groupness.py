"""26_cmv_phase_groupness.py — do the CMV (SE Missouri) phases form real groups?

Frequency-level test on Williams 1954's southeast-Missouri survey assemblages
(data/raw/williams1954_cmv_counts.tsv: 50 assemblages, 28 sites, 6 physiographic
regions). This is the CMV-scale analog of the within-basin LMV test (analyses
23-25) and a frequency-level extension of Fox 1998 (who used cluster/MDS on a
subset and found the phases fail as groups).

Asks:
  1. Do Williams's physiographic regions (proxy for his phases) cohere as groups?
     within-region vs between-region Brainerd-Robinson similarity (permutation
     test), and cultural F_ST by region. Compare F_ST to the LMV within-basin
     value (~0.04).
  2. Does the CMV seriate into one coherent order or fragment? CA1 axis + IDSS
     seriation-group structure.
  3. The decoration contrast Lipo flagged: plain vs Woodland vs painted/slipped
     vs Mississippian-decorated fractions by region, and the distinctive NORTHERN
     decorated types (Wickliffe/O'Byam/Matthews/Kimmswick) the LMV lacks.

Read-only on the manuscript. Writes output/cmv_phase_groupness.md and
figures/figX_cmv_mds.png.

Usage: .venv/bin/python analyses/26_cmv_phase_groupness.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import networkx as nx

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import make_figures as mf  # noqa: E402
from mls_emergence.signatures.assortativity import similarity_matrix  # noqa: E402

TSV = ROOT / "data" / "raw" / "williams1954_cmv_counts.tsv"
OUT_MD = ROOT / "output" / "cmv_phase_groupness.md"
OUT_FIG = ROOT / "figures" / "figX_cmv_mds.png"
MIN_TOTAL = 30

# --- type categories (Williams 1954 taxonomy) ------------------------------- #
# Woodland / Baytown tradition (clay/sand tempered, pre-Mississippian)
WOODLAND = {
    "Barnes Plain", "Barnes Cord Marked", "Barnes Fabric Impressed",
    "Barnes Pinched", "Barnes Punctate", "Baytown Plain",
    "Mulberry Creek Cord Marked", "Withers Fabric Impressed",
    "Wheeler Check Stamped", "Larto Red Filmed", "Cord marked with cord design",
    "Grit Tempered Plain", "Unidentified Plain Grit", "Decorated Sand Tempered",
    "Unidentified Clay Tempered Incised", "Unidentified Clay Tempered Punctate",
    "Unidentified Clay Tempered Cord Marked Incised", "Unidentified Cord Marked",
}
# Mississippian shell-tempered PLAIN
PLAIN_SHELL = {
    "Mississippi Plain", "Bell Plain", "Kimmswick Plain", "Wickliffe Plain",
    "Polished Shell Ware", "Thin ware", "Thin Ware", "Thin ware (shell tempered)",
    "Fine ware", "Fine Ware",
}
# Painted / slipped (red, white, negative, brown slip)
PAINTED = {
    "Varney Red", "Old Town Red Filmed", "Nodena Red and White",
    "Hollywood White Filmed", "Angel Negative Painted", "Polished Brown",
    "Brown Slip", "Brown slip", "Unidentified Incised Red Filmed",
}
# Distinctive NORTHERN decorated types largely absent in the LMV St. Francis set
NORTHERN_DEC = {
    "Wickliffe Incised", "Wickliffe Cord Marked", "O'Byam Engraved",
    "O'Byam Incised", "Matthews Incised", "Kimmswick Fabric Impressed",
    "Kimmswick Fabric Marked",
}
# Decorated Mississippian shared with / typical of the LMV
SHARED_DEC = {
    "Parkin Punctate", "Ranch Incised", "Kent Incised", "Mound Place Incised",
    "Wallace Incised", "Manly Punctate", "Rhodes-like", "L'eau Noir Incised",
}


def load_cmv():
    d = pd.read_csv(TSV, sep="\t")
    d["count"] = pd.to_numeric(d["count"], errors="coerce").fillna(0)
    d["asm"] = (d["site_name"].astype(str) + " [" +
                d["collection"].astype(str) + "]")
    region = d.groupby("asm")["region"].first()
    M = d.pivot_table(index="asm", columns="type_name", values="count",
                      aggfunc="sum", fill_value=0)
    tot = M.sum(1)
    keep = tot[tot >= MIN_TOTAL].index
    M = M.loc[keep]
    region = region.loc[keep]
    return M, region


def communities(counts, seed=0):
    S = similarity_matrix(counts)
    n = S.shape[0]
    G = nx.Graph(); G.add_nodes_from(range(n))
    for i in range(n):
        for j in range(i + 1, n):
            if S[i, j] > 0:
                G.add_edge(i, j, weight=S[i, j])
    comms = nx.community.greedy_modularity_communities(G, weight="weight")
    lab = np.full(n, -1)
    for c, mem in enumerate(comms):
        for m in mem:
            lab[m] = c
    return lab, float(nx.community.modularity(G, comms, weight="weight")), len(comms)


def region_coherence(S, region_labels, n_perm=4999, seed=0):
    """within-region minus between-region mean BR; permutation p."""
    n = S.shape[0]
    iu = np.triu_indices(n, 1)
    same = (region_labels[iu[0]] == region_labels[iu[1]])
    s = S[iu]
    obs = s[same].mean() - s[~same].mean()
    rng = np.random.default_rng(seed)
    cnt = 0
    for _ in range(n_perm):
        p = rng.permutation(n)
        rl = region_labels[p]
        sm = (rl[iu[0]] == rl[iu[1]])
        if sm.sum() == 0 or (~sm).sum() == 0:
            continue
        if (s[sm].mean() - s[~sm].mean()) >= obs:
            cnt += 1
    return float(obs), float((cnt + 1) / (n_perm + 1))


def frac(M, typeset):
    cols = [c for c in M.columns if c in typeset]
    return M[cols].sum(1) / M.sum(1) if cols else pd.Series(0.0, index=M.index)


def mds2(counts):
    S = similarity_matrix(counts)
    D = (200.0 - S) ** 2
    n = D.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ D @ J
    w, V = np.linalg.eigh(B)
    idx = np.argsort(-w)[:2]
    coords = V[:, idx] * np.sqrt(np.maximum(w[idx], 0))
    # Deterministic sign convention (eigh sign is arbitrary): orient each axis so
    # its largest-magnitude coordinate is positive, for cross-platform reproducibility.
    for j in range(coords.shape[1]):
        if coords[np.argmax(np.abs(coords[:, j])), j] < 0:
            coords[:, j] = -coords[:, j]
    return coords


def main():
    M, region = load_cmv()
    counts = M.to_numpy(float)
    n = counts.shape[0]
    reg = region.to_numpy()
    reg_codes = pd.Categorical(reg).codes
    S = similarity_matrix(counts)

    L = ["# Do the CMV (SE Missouri) phases form real groups? (Williams 1954)",
         "",
         f"{n} assemblages with >= {MIN_TOTAL} sherds, {len(set(reg))} physiographic "
         "regions (proxy for Williams's phases). Frequency-level test, the CMV-scale",
         "analog of the LMV within-basin analysis.", ""]

    # 1. region coherence
    obs, p = region_coherence(S, reg_codes)
    # cultural F_ST by region (pooled region count vectors)
    region_pools = np.array([counts[reg_codes == c].sum(0)
                             for c in np.unique(reg_codes)
                             if counts[reg_codes == c].sum() > 0])
    fst_region = mf.cultural_fst(region_pools)
    lab, Q, ncom = communities(counts)
    L += ["## 1. Region/phase coherence",
          f"- Within-region minus between-region BR similarity: {obs:+.1f} BR units "
          f"(permutation p = {p:.4f}).",
          f"- Cultural F_ST by region: {fst_region:+.3f}  "
          f"(LMV within-basin F_ST ~ +0.04 for comparison).",
          f"- Data-driven communities (greedy modularity): {ncom}, Q = {Q:.3f}.",
          "  Do detected communities line up with regions? cross-tab below.", ""]
    ct = pd.crosstab(pd.Series(reg, name="region"),
                     pd.Series(lab, name="community"))
    L += ["```", ct.to_string(), "```", ""]

    # 2. seriation coherence
    ca1, _, _ = mf.correspondence_axis(counts)
    sg = mf.seriation_groups(counts, cont=0.10)
    nbridge = sum(1 for i in range(n) if len(sg["membership"].get(i, [])) > 1)
    L += ["## 2. Seriation coherence (IDSS, cont=0.10)",
          f"- n_groups = {sg['n_groups']}, bridges = {nbridge}/{n}. "
          "Many groups / few bridges = fragments rather than one coherent order.", ""]

    # 3. decoration contrast
    cats = {"Woodland": WOODLAND, "plain shell": PLAIN_SHELL,
            "painted/slipped": PAINTED, "northern decorated": NORTHERN_DEC,
            "shared decorated": SHARED_DEC}
    L += ["## 3. Decoration composition by region (sherd-weighted)",
          "| region | n | Woodland | plain shell | painted | northern dec | shared dec |",
          "|---|---|---|---|---|---|---|"]
    for rg in sorted(set(reg)):
        sub = M[region.to_numpy() == rg]
        tot = sub.sum().sum()
        row = [rg, str((region.to_numpy() == rg).sum())]
        for cat, ts in cats.items():
            cols = [c for c in sub.columns if c in ts]
            row.append(f"{(sub[cols].sum().sum()/tot):.3f}" if cols else "0.000")
        L += ["| " + " | ".join(row) + " |"]
    # overall CMV
    tot_all = M.sum().sum()
    over = ["**CMV overall**", str(n)]
    for cat, ts in cats.items():
        cols = [c for c in M.columns if c in ts]
        over.append(f"{(M[cols].sum().sum()/tot_all):.3f}" if cols else "0.000")
    L += ["| " + " | ".join(over) + " |", "",
          "Note: 'northern decorated' (Wickliffe/O'Byam/Matthews/Kimmswick) are "
          "the CMV-diagnostic types largely absent from the LMV St. Francis basin.", ""]

    OUT_MD.write_text("\n".join(L), encoding="utf-8")

    # figure: MDS colored by region
    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["Arial", "DejaVu Sans"],
                         "font.size": 8})
    xy = mds2(counts)
    fig, ax = plt.subplots(figsize=(5.5, 5.0))
    palette = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
    for c, rg in enumerate(sorted(set(reg))):
        m = reg == rg
        ax.scatter(xy[m, 0], xy[m, 1], s=30, c=palette[c % len(palette)],
                   edgecolor="white", linewidth=0.4, label=rg)
    ax.legend(fontsize=6, frameon=False, loc="best")
    ax.set_xlabel("MDS axis 1"); ax.set_ylabel("MDS axis 2")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.savefig(OUT_FIG, dpi=300, bbox_inches="tight")

    print("\n".join(L))
    print(f"\nwrote {OUT_MD}\nwrote {OUT_FIG}")


if __name__ == "__main__":
    main()

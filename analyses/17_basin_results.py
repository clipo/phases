"""17_basin_results.py — consolidated Results numbers for the drainage basin.

Recomputes the headline manuscript quantities on the drainage-defined basin
(n given by analyses/16): the four-signature ordinal trajectory (with the CA axis
oriented by radiocarbon so increasing = later), the convergence-score trend, the
empirical four-signature correlation matrix, the IDSS group structure, and the
within-basin mound-height ranking. Writes output/basin_results.md.

Usage: .venv/bin/python analyses/17_basin_results.py
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import make_figures as mf  # noqa: E402
grid = importlib.import_module("12_sensitivity_grid")
ch = importlib.import_module("11_chronology_14c")
from mls_emergence.dataio.pfg import load_pfg_counts  # noqa: E402
from mls_emergence.dataio.settlement import load_lmv, join_pfg_to_lmv, normalize_grid  # noqa: E402

OUT = ROOT / "output" / "basin_results.md"
PARKIN_BROAD = "11-N-1"


def oriented_ca(counts):
    M = counts.to_numpy(float)
    ca1, _, _ = mf.correspondence_axis(M)
    ca = pd.Series(ca1, index=counts.index)
    # orient by pooled 14C median age of matching proveniences
    df = ch.parse_dates()
    grid_ad = np.arange(1000, 1750)
    cur_norm = {mf.norm_name(a): a for a in counts.index}
    anch = {}
    for prov, g in df.groupby("prov"):
        pn = mf.norm_name(prov)
        hit = cur_norm.get(pn) or next((v for k, v in cur_norm.items() if pn and (pn in k or k in pn)), None)
        if hit is not None:
            anch[hit] = ch.median_ad(g, grid_ad)
    if len(anch) >= 3:
        a_ca = [ca[h] for h in anch]
        a_age = [anch[h] for h in anch]
        rho, _ = spearmanr(a_ca, a_age)
        if np.isfinite(rho) and rho < 0:
            ca = -ca
    return ca, len(anch)


def main():
    counts, coords = mf._load_curated()
    n = counts.shape[0]
    M = counts.to_numpy(float)
    ca, n_anch = oriented_ca(counts)

    L = [f"# Drainage-basin Results (n = {n} curated assemblages)", "",
         "CA axis oriented by radiocarbon so increasing = later.", ""]

    # per-bin signatures with oriented bins
    cdf = coords.dropna(); have = list(cdf.index)
    cc = cdf[["Latitude", "Longitude"]].to_numpy(float); cc_c = cc - cc.mean(0)
    sil = {kk: mf.silhouette_mean(cc_c, mf._kmeans_labels(cc_c, kk, seed=7)) for kk in range(2, 7)}
    k = max(sil, key=sil.get)
    cl = mf._kmeans_labels(cc_c, k, seed=7); cluster_of = dict(zip(have, cl))
    ca_have = ca.reindex(have); counts_have = counts.reindex(have)
    bins = pd.qcut(ca_have, 6, labels=False, duplicates="drop")

    # seriation membership per assemblage
    sg = mf.seriation_groups(M, cont=0.10); memb = sg["membership"]
    idx = list(counts.index)
    memb_count = {idx[i]: len(memb.get(i, [])) for i in range(len(idx))}
    g2a = {}
    for ai, gl in memb.items():
        for gi in gl:
            g2a.setdefault(gi, []).append(ai)
    n_groups = sg["n_groups"]; max_sz = max(len(v) for v in g2a.values())
    nbridge = sum(1 for i in range(n) if len(memb.get(i, [])) > 1)
    mc = sorted([(idx[i], len(memb.get(i, []))) for i in range(n)], key=lambda x: -x[1])
    p_rank = next((r + 1 for r, (nm, _) in enumerate(mc) if nm == "Parkin"), None)
    p_memb = next((v for nm, v in mc if nm == "Parkin"), None)

    rows = {}
    for b in sorted(pd.Series(bins).dropna().unique()):
        ids = [i for i in have if not pd.isna(bins[i]) and bins[i] == b]
        sc = counts_have.loc[ids].to_numpy(float); scl = np.array([cluster_of[i] for i in ids])
        sco = cc_c[[have.index(i) for i in ids]]
        nd = []
        for c in np.unique(scl):
            p = sc[scl == c].sum(0)
            if p.sum() < 2 or (p > 0).sum() < 2:
                continue
            tf, te = mf.theta_f(p), mf.theta_e(p)
            if np.isfinite(tf) and te > 0:
                nd.append(abs(1 - tf / te))
        rep = np.unique(scl)
        rows[int(b)] = {
            "neutral": np.mean(nd) if nd else np.nan,
            "seriation": np.mean([memb_count[i] for i in ids]),
            "fst": mf.cultural_fst(np.array([sc[scl == c].sum(0) for c in rep])) if len(rep) >= 2 else np.nan,
            "spatial": mf.boundary_excess(sc, sco, seed=7) if len(ids) >= 4 else np.nan,
        }
    P = pd.DataFrame(rows).T.sort_index()

    def tr(col):
        s = P[col].dropna()
        return spearmanr(s.index.to_numpy(float), s.values) if len(s) >= 3 else (np.nan, np.nan)

    L += ["## Four-signature trajectory (k = %d clusters, 6 bins)" % k, "",
          "| signature | Spearman rho | p |", "|---|---|---|"]
    for c in ["neutral", "seriation", "fst", "spatial"]:
        r, p = tr(c)
        L.append(f"| {c} | {r:+.2f} | {p:.3f} |")
    # convergence score on the 3 continuous signatures
    cont3 = P[["neutral", "fst", "spatial"]]
    z = (cont3 - cont3.mean()) / cont3.std(ddof=0).replace(0, np.nan)
    conv = z.mean(axis=1).dropna()
    cr, cp = spearmanr(conv.index.to_numpy(float), conv.values) if len(conv) >= 3 else (np.nan, np.nan)
    cslope = np.polyfit(conv.index.to_numpy(float), conv.values, 1)[0] if len(conv) >= 2 else np.nan
    L += ["", f"- Convergence score (3 continuous signatures): slope = {cslope:+.3f}, "
          f"Spearman rho = {cr:+.2f}, p = {cp:.3f}.",
          "- Seriation fragmentation (per-assemblage membership vs CA): see table above.", ""]

    # correlation matrix (4 signatures across bins)
    C = P[["neutral", "seriation", "fst", "spatial"]].corr(method="spearman")
    offdiag = C.values[np.triu_indices(4, 1)]
    L += ["## Empirical 4-signature correlation (across bins)",
          f"- Mean |r| off-diagonal = {np.nanmean(np.abs(offdiag)):.2f} (simulation genuine 0.88).",
          f"- neutral-seriation r = {C.loc['neutral','seriation']:+.2f}; "
          f"fst-spatial r = {C.loc['fst','spatial']:+.2f}.", ""]

    L += ["## IDSS group structure (cont = 0.1)",
          f"- n_groups = {n_groups}, max size = {max_sz}, bridges = {nbridge}/{n}.",
          f"- Parkin: rank {p_rank}/{n}, membership {p_memb}.",
          f"- 14C anchors matching curated assemblages: {n_anch}.", ""]

    # settlement: mound-height ranking on the broad basin
    broad = load_pfg_counts(ROOT / "data" / "raw" / "PFGData_sherds.csv")
    if not broad.index.is_unique:
        broad = broad.groupby(level=0).sum()
    lmv = load_lmv(ROOT / "data" / "LMVData_locations.csv")
    joined, _ = join_pfg_to_lmv(broad, lmv)
    bm = joined.dropna(subset=["Easting", "Northing"]).copy()
    members = mf._basin_members("broad")
    bm = bm[[str(i) in members for i in bm.index]].copy()
    lmv2 = pd.read_csv(ROOT / "data" / "LMVData-22March2006.csv").dropna(subset=["Number"])
    lmv2["_k"] = lmv2["Number"].astype(str).map(normalize_grid)
    lmv2 = lmv2.drop_duplicates("_k").set_index("_k")
    ext = lmv2.reindex(pd.Index([normalize_grid(str(i)) for i in bm.index]))
    ext.index = bm.index
    ht = pd.to_numeric(ext["Max Mound Height (ft)"], errors="coerce").dropna()
    ht = ht[ht > 0].sort_values(ascending=False)
    p_ht_rank = int((ht.index == PARKIN_BROAD).argmax() + 1) if PARKIN_BROAD in ht.index else None
    primacy = float(ht.iloc[0] / ht.iloc[1]) if len(ht) >= 2 else np.nan
    L += ["## Settlement (broad basin)",
          f"- Broad basin sites: {len(bm)}.",
          f"- Mound-height ranking: {len(ht)} sites with height > 0; Parkin "
          f"{'rank ' + str(p_ht_rank) + '/' + str(len(ht)) if p_ht_rank else 'unranked'} "
          f"at {ht.get(PARKIN_BROAD, float('nan')):.0f} ft; tallest/second = {primacy:.2f}.", ""]

    OUT.write_text("\n".join(L), encoding="utf-8")
    print(f"wrote {OUT}")
    print("\n".join(L))


if __name__ == "__main__":
    main()

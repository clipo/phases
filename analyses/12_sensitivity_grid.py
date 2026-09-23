"""12_sensitivity_grid.py — robustness of the no-convergence verdict.

Re-runs the basin convergence test across the grid of analyst choices the
reviewers flagged (latitude cut x number of ordinal bins x number of spatial
clusters k), and tabulates, for each cell, the ordinal Spearman trend of each
continuous signature (neutral departure, cultural F_ST, spatial boundary) and
whether the cell shows CONVERGENCE (all three rising together, rho > +0.3).
If no cell converges, the negative result is robust to the garden of forking
paths. The IDSS continuity-threshold sensitivity is separately handled in
analyses/09 (group counts 51/124/683 at cont 0.05/0.1/0.2; bridge topology
stable), so cont is not re-gridded here (it does not enter the convergence score).

Writes output/sensitivity_grid.md. Read-only on the manuscript.

Usage: .venv/bin/python analyses/12_sensitivity_grid.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from mls_emergence.dataio.coords import read_assemblage_xy
from mls_emergence.dataio.matrix import read_analysis_matrix  # noqa: E402
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import make_figures as mf  # noqa: E402

DATA = ROOT / "data" / "raw"
OUT = ROOT / "output" / "sensitivity_grid.md"
DECORATED = mf.DECORATED_TYPES
RISE = 0.30  # Spearman threshold for "rising"


def load_curated_full():
    """Curated decorated set WITHOUT the basin latitude filter (whole-LMV, 55)."""
    cur = read_analysis_matrix().dropna(
        subset=["Assemblages"])
    cur["Assemblages"] = cur["Assemblages"].astype(str).str.strip()
    cur = cur.drop_duplicates(subset=["Assemblages"], keep="first").set_index("Assemblages")
    cols = [c for c in DECORATED if c in cur.columns]
    counts = cur[cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    counts = counts[counts.sum(axis=1) > 0]
    xy = read_assemblage_xy(DATA / "mainfort-pfg-cplXY.txt")
    xy["Assemblages"] = xy["Assemblages"].astype(str).str.strip()
    xy = xy.drop_duplicates(subset=["Assemblages"], keep="first").set_index("Assemblages")
    coords = xy.reindex(counts.index)[["Latitude", "Longitude"]].apply(pd.to_numeric, errors="coerce")
    return counts, coords


def per_bin_signatures(counts, coords, n_bins, k_fixed=None):
    M = counts.to_numpy(float)
    ca1, _, _ = mf.correspondence_axis(M)
    cdf = coords.dropna()
    have = list(cdf.index)
    cc = cdf[["Latitude", "Longitude"]].to_numpy(float)
    cc_c = cc - cc.mean(0)
    if k_fixed:
        k = min(k_fixed, len(have) - 1)
    else:
        sil = {kk: mf.silhouette_mean(cc_c, mf._kmeans_labels(cc_c, kk, seed=7)) for kk in range(2, 7)}
        k = max(sil, key=sil.get)
    cl = mf._kmeans_labels(cc_c, k, seed=7)
    cluster_of = dict(zip(have, cl))
    ca = pd.Series(ca1, index=counts.index)
    ca_have = ca.reindex(have)
    counts_have = counts.reindex(have)
    bins = pd.qcut(ca_have, n_bins, labels=False, duplicates="drop")
    rows = {}
    for b in sorted(pd.Series(bins).dropna().unique()):
        ids = [i for i in have if not pd.isna(bins[i]) and bins[i] == b]
        sc = counts_have.loc[ids].to_numpy(float)
        scl = np.array([cluster_of[i] for i in ids])
        sco = cc_c[[have.index(i) for i in ids]]
        nd = []
        for c in np.unique(scl):
            p = sc[scl == c].sum(0)
            if p.sum() < 2 or (p > 0).sum() < 2:
                continue
            tf, te = mf.theta_f(p), mf.theta_e(p)
            if np.isfinite(tf) and te > 0:
                nd.append(abs(1 - tf / te))
        nd = np.mean(nd) if nd else np.nan
        rep = np.unique(scl)
        fst = mf.cultural_fst(np.array([sc[scl == c].sum(0) for c in rep])) if len(rep) >= 2 else np.nan
        be = mf.boundary_excess(sc, sco, seed=7) if len(ids) >= 4 else np.nan
        rows[b] = {"neutral": nd, "fst": fst, "spatial": be}
    return pd.DataFrame(rows).T.sort_index(), k


def trend(series):
    s = series.dropna()
    if len(s) < 3:
        return np.nan
    r, _ = spearmanr(s.index.to_numpy(float), s.values)
    return r


def main():
    counts_full, coords_full = load_curated_full()
    lat = coords_full["Latitude"]

    L = ["# Sensitivity of the no-convergence verdict (forking-paths grid)", "",
         "For each cell: ordinal Spearman trend of neutral departure / cultural F_ST / "
         "spatial boundary along the CA axis, and the CONVERGENCE verdict (all three "
         f"rising, rho > +{RISE:.1f}). k = number of spatial clusters (auto = silhouette).", "",
         "| lat cut | n | bins | k | neutral | F_ST | spatial | converges? |",
         "|---|---|---|---|---|---|---|---|"]
    any_conv = False
    for latcut in (34.0, 34.5, 35.0):
        keep = lat.index[(lat >= latcut).fillna(False).to_numpy()]
        c_sub = counts_full.loc[keep]
        co_sub = coords_full.loc[keep]
        n = len(c_sub)
        for nb in (4, 6, 8):
            for kspec in ("auto", 3, 4):
                kf = None if kspec == "auto" else kspec
                try:
                    panel, kused = per_bin_signatures(c_sub, co_sub, nb, k_fixed=kf)
                except Exception as e:
                    L.append(f"| {latcut} | {n} | {nb} | {kspec} | err | err | err | {e.__class__.__name__} |")
                    continue
                rn, rf, rs = trend(panel["neutral"]), trend(panel["fst"]), trend(panel["spatial"])
                conv = all((x is not np.nan and np.isfinite(x) and x > RISE) for x in (rn, rf, rs))
                any_conv = any_conv or conv
                kshow = f"{kused}" if kspec == "auto" else f"{kspec}"
                # An undefined trend is a finding about the design (too few bins
                # with two or more clusters in them, or too few assemblages per
                # bin), so it is written out as that and never as "nan".
                show = lambda v: f"{v:+.2f}" if np.isfinite(v) else "undefined"
                L.append(f"| {latcut} | {n} | {nb} | {kshow} | {show(rn)} | {show(rf)} | "
                         f"{show(rs)} | {'YES' if conv else 'no'} |")
    # The reading is computed from the table, not written in advance: an earlier
    # version asserted "consistently negative" neutral trends while the table
    # beside it showed them positive (2026-09-22).
    n_conv = sum(1 for l in L if l.endswith("| YES |"))
    n_cells = sum(1 for l in L if l.startswith("| ") and l.rstrip().endswith(("| YES |", "| no |")))
    neu = [float(l.split("|")[5]) for l in L if l.startswith("| ") and l.split("|")[5].strip()[:1] in "+-"]
    L += ["",
          f"**Convergence (all three continuous signatures rising, rho > {RISE:+.1f}) appears in "
          f"{n_conv} of {n_cells} cells.** Neutral-departure trends run {min(neu):+.2f} to "
          f"{max(neu):+.2f} across cells. These are raw, unrarefied trends: assemblage size rises "
          "along the seriation, and the record-matched recovery experiment shows only cultural "
          "F_ST is interpretable at this resolution, so a cell that converges here is not "
          "evidence of closure. The grid shows how far raw trends move with the analyst's "
          "choices of latitude cut, bin count and cluster number."]
    OUT.write_text("\n".join(L), encoding="utf-8")
    print(f"wrote {OUT}")
    print("\n".join(L))


if __name__ == "__main__":
    main()

"""28_macro_boundary.py — the CMV/LMV macro transect: boundary or smooth cline?

Combines the georeferenced CMV (Williams 1954, SE Missouri) and LMV (Phillips-
Ford-Griffin, St. Francis and south) decorated-ceramic assemblages into one
north-south transect and asks whether decorated-type similarity STEPS DOWN at the
CMV/LMV boundary (a real interaction/transmission boundary) or declines smoothly
with distance (isolation-by-distance). Tests:

  1. Map of the full transect (the spatial gap between the blocks).
  2. Distance-decay of decorated Brainerd-Robinson similarity across the transect.
  3. Boundary excess: within-region minus between-region similarity at MATCHED
     geographic distance (CMV-vs-LMV split). >0 beyond IBD = a step.
  4. Spatial-drift null (Lipo et al. 2021) on the combined geography: does the
     observed step exceed neutral drift on this layout?

CRITICAL CAVEAT (see output): the CMV/LMV decorated difference is confounded with
time (early-middle CMV types vs late LMV types; analyses 27). The spatial tests
here cannot separate a synchronic spatial boundary from early-CMV/late-LMV
succession. Reported as spatial PATTERN, not as a demonstrated synchronic boundary.

Read-only on the manuscript. Writes output/macro_boundary.md and
figures/figS6_macro_transect.png.

Usage: .venv/bin/python analyses/28_macro_boundary.py
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
rep = importlib.import_module("27_cmv_lmv_repertoire")
demo = importlib.import_module("25_drift_vs_groups_demo")
sd = importlib.import_module("23_phases_vs_spatial_drift")
from mls_emergence.signatures.assortativity import similarity_matrix, mantel  # noqa: E402

OUT_MD = ROOT / "output" / "macro_boundary.md"
OUT_FIG = ROOT / "figures" / "figS6_macro_transect.png"
MIN_DEC = 10
C_CMV = "0.62"
C_LMV = "0.28"


def _norm_grid(s):
    return s.strip().upper().replace(" ", "") if isinstance(s, str) else ""


def grid_coords():
    """LMS grid -> lat/lon from the digitized Phillips/Ford/PFG point files."""
    frames = []
    for f in ["PHILLIPS.SHP", "GEOFILE.SHP", "PFGOUT.SHP"]:
        g = gpd.read_file(ROOT / "data" / "CMV" / f).to_crs(4326)
        g["grid"] = g["NAME1_"].map(_norm_grid)
        g["lat"] = g.geometry.y; g["lon"] = g.geometry.x
        frames.append(g[["grid", "lat", "lon"]])
    G = pd.concat(frames, ignore_index=True)
    G = G[(G.lat.between(33, 38)) & (G.lon.between(-92, -88))]
    G = G[G.grid != ""].drop_duplicates("grid").set_index("grid")
    return G


def geo_km(coords):
    return sd.geo_km(coords)


def main():
    G = grid_coords()

    # --- LMV decorated (PFG), geolocated via grid -------------------------- #
    lmv = rep.lmv_decorated()                 # assemblage x canonical decorated
    lmv_grid = [_norm_grid(str(i)) for i in lmv.index]
    lat = np.array([G.loc[g, "lat"] if g in G.index else np.nan for g in lmv_grid])
    lon = np.array([G.loc[g, "lon"] if g in G.index else np.nan for g in lmv_grid])
    lmv = lmv.assign(lat=lat, lon=lon, region="LMV")
    lmv = lmv.dropna(subset=["lat"])

    # --- CMV decorated (Williams), coords from processed file -------------- #
    cmv = rep.cmv_decorated()
    cc = pd.read_csv(ROOT / "data" / "processed" / "williams1954_cmv_coords.tsv", sep="\t")
    # cmv index is "site [collection]"; map by site name
    cc["key"] = cc["site_name"].astype(str)
    site_of = {asm: asm.split(" [")[0] for asm in cmv.index}
    latlon = cc.set_index("key")[["lat", "lon"]]
    clat = [latlon.loc[site_of[a], "lat"] if site_of[a] in latlon.index else np.nan
            for a in cmv.index]
    clon = [latlon.loc[site_of[a], "lon"] if site_of[a] in latlon.index else np.nan
            for a in cmv.index]
    cmv = cmv.assign(lat=clat, lon=clon, region="CMV").dropna(subset=["lat"])

    canon = rep.CANON
    keep_cols = canon + ["lat", "lon", "region"]
    both = pd.concat([cmv[keep_cols], lmv[keep_cols]], ignore_index=True)
    dec_tot = both[canon].sum(1)
    both = both[dec_tot >= MIN_DEC].reset_index(drop=True)
    counts = both[canon].to_numpy(float)
    coords = both[["lat", "lon"]].to_numpy(float)
    region = both["region"].to_numpy()
    n = len(both)
    n_cmv = int((region == "CMV").sum()); n_lmv = int((region == "LMV").sum())

    geo = geo_km(coords)
    S = similarity_matrix(counts)

    # spatial gap between blocks
    cmv_lat = coords[region == "CMV", 0]; lmv_lat = coords[region == "LMV", 0]
    gap_km = (cmv_lat.min() - lmv_lat.max()) * 111.32

    L = ["# CMV/LMV macro transect: boundary or smooth cline?",
         "",
         f"{n} decorated assemblages (>= {MIN_DEC} decorated sherds): "
         f"CMV {n_cmv}, LMV {n_lmv}. Decorated types harmonized (analyses 27).",
         f"Latitudinal gap between the southernmost CMV and northernmost LMV "
         f"assemblage in this set: ~{gap_km:.0f} km.", ""]

    # distance-decay
    r_dd, p_dd = mantel(S, geo, n_perm=4999, seed=1)
    L += ["## Distance-decay across the transect",
          f"- Mantel r (decorated BR similarity vs geographic distance) = "
          f"{r_dd:+.3f} (p = {p_dd:.4f}).", ""]

    # boundary excess: CMV vs LMV, distance-controlled
    lab = (region == "LMV").astype(int)
    be = demo_boundary(counts, geo, lab)
    iu = np.triu_indices(n, 1)
    same = lab[iu[0]] == lab[iu[1]]
    raw = float(S[iu][same].mean() - S[iu][~same].mean())
    L += ["## Boundary excess (CMV vs LMV, distance-controlled)",
          f"- raw within-minus-between region BR = {raw:+.1f}.",
          f"- distance-controlled boundary excess = {be:+.1f} BR units.",
          "  >0 beyond IBD = a step in decorated similarity at the boundary.", ""]

    # spatial-drift null on combined geography (proportional structure only)
    totals = counts.sum(1)
    K = counts.shape[1]
    Wd = demo.drift_weights(coords, length_km=24.0)
    rng_stats = []
    for s in range(60):
        M = demo.simulate(Wd, K=K, total_per_node=totals, m_between=0.02,
                          steps=700, seed=9000 + s)
        if M.sum() == 0:
            continue
        Ssim = similarity_matrix(M)
        rsim, _ = mantel(Ssim, geo, n_perm=1, seed=s)
        besim = demo_boundary(M, geo, lab)
        rng_stats.append((rsim, besim))
    arr = np.array(rng_stats)
    dd_lo, dd_hi = np.percentile(arr[:, 0], [2.5, 97.5])
    be_lo, be_hi = np.percentile(arr[:, 1], [2.5, 97.5])
    L += ["## Spatial-drift null on the combined geography (60 sims)",
          f"- distance-decay r: observed {r_dd:+.3f}, null 95% [{dd_lo:+.3f}, {dd_hi:+.3f}] "
          f"({'inside' if dd_lo<=r_dd<=dd_hi else 'OUTSIDE'}).",
          f"- boundary excess: observed {be:+.1f}, null 95% [{be_lo:+.1f}, {be_hi:+.1f}] "
          f"({'inside' if be_lo<=be<=be_hi else 'OUTSIDE'}).",
          "  Note: the drift null randomizes type labels, so it cannot reproduce "
          "the mutually-exclusive distinctive repertoires; it tests only whether the "
          "PROPORTIONAL step exceeds neutral spatial drift.", ""]

    L += ["## CRITICAL CAVEAT",
          "The CMV/LMV decorated difference is confounded with time (early-middle CMV",
          "types vs late LMV types; analyses 27). These spatial statistics describe",
          "PATTERN; they do not separate a synchronic spatial boundary from early-CMV/",
          "late-LMV succession. A synchronic test needs chronological control per",
          "assemblage (independent CMV radiocarbon).", ""]

    OUT_MD.write_text("\n".join(L), encoding="utf-8")

    # --- figure: map + distance-decay ------------------------------------- #
    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["Arial", "DejaVu Sans"], "font.size": 8})
    fig, (axm, axd) = plt.subplots(1, 2, figsize=(7.2, 4.4),
                                   gridspec_kw={"width_ratios": [1.0, 1.0]})
    for rg, c in [("CMV", C_CMV), ("LMV", C_LMV)]:
        m = region == rg
        axm.scatter(coords[m, 1], coords[m, 0], s=18, c=c, edgecolor="white",
                    linewidth=0.3, label=f"{rg} (n={int(m.sum())})")
    axm.set_xlabel("longitude"); axm.set_ylabel("latitude")
    axm.legend(fontsize=6, frameon=False, loc="lower left")
    axm.set_title("Decorated assemblages, CMV vs LMV", fontsize=8)
    axm.set_aspect(1.25)

    d = geo[iu]; s = S[iu]
    axd.scatter(d[~same], s[~same], s=4, c="#999999", alpha=0.4, label="between-region")
    axd.scatter(d[same], s[same], s=4, c="#333333", alpha=0.5, label="within-region")
    axd.set_xlabel("geographic distance (km)"); axd.set_ylabel("decorated BR similarity")
    axd.legend(fontsize=6, frameon=False)
    axd.set_title(f"Distance-decay (Mantel r={r_dd:+.2f})", fontsize=8)
    for sp in ("top", "right"):
        axd.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=200, bbox_inches="tight")

    print("\n".join(L))
    print(f"\nwrote {OUT_MD}\nwrote {OUT_FIG}")


def demo_boundary(counts, geo_dist, labels, n_bins=4):
    """distance-controlled within-minus-between BR for supplied labels."""
    S = similarity_matrix(counts)
    n = S.shape[0]
    iu = np.triu_indices(n, 1)
    d = geo_dist[iu]; s = S[iu]
    same = labels[iu[0]] == labels[iu[1]]
    if same.sum() == 0 or (~same).sum() == 0:
        return 0.0
    if d.max() == d.min():
        return float(s[same].mean() - s[~same].mean())
    edges = np.linspace(d.min(), d.max() + 1e-9, n_bins + 1)
    gaps = []
    for b in range(n_bins):
        inb = (d >= edges[b]) & (d < edges[b + 1])
        w = inb & same; bt = inb & ~same
        if w.sum() and bt.sum():
            gaps.append(float(s[w].mean() - s[bt].mean()))
    return float(np.mean(gaps)) if gaps else float(s[same].mean() - s[~same].mean())


if __name__ == "__main__":
    main()

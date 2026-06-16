"""37_lmv_drift_groups.py - how many phase-like groups does drift produce across
the wider lower Mississippi Valley, and how strong are they?

Runs the forward neutral-drift model of analyses/33 over all 55 Mainfort-PFG
decorated assemblages of the lower Mississippi Valley, on the river-network
geography, with no boundary imposed (500 realizations). For each realization it
detects spatially coherent groups (greedy-modularity communities on the
Brainerd-Robinson similarity graph) and records, separately for the
contemporaneous (pure-space) and time-transgressive samples, the number of
groups and the between-group cultural F_ST. The observed data and the seven named
culture-historical phases (Mainfort 1996) are marked for comparison.

The point is descriptive, not a claim of "emergence": drift structured by
geography and time already produces a few coherent, phase-like groups, fewer
than the named scheme draws, at about the observed level of differentiation.
Figure B (analyses/35) then asks whether one named phase, Parkin, exceeds drift.

Per-run results cache to output/lmv_drift_groups_runs.csv and the example
realization to output/lmv_drift_groups_example.csv; delete those to force a rerun.

Writes output/lmv_drift_groups.md and figures/fig_lmv_drift_groups.png.

Usage: PYTHONPATH=src python3 analyses/37_lmv_drift_groups.py
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))
sys.path.insert(0, str(ROOT / "src"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import geopandas as gpd  # noqa: E402
import make_figures as mf  # noqa: E402  (house style + DECORATED_TYPES)
import make_map as mm  # noqa: E402  (river basemap + river-network distance)
m33 = importlib.import_module("33_time_aware_emergence")
m35 = importlib.import_module("35_basin_pullout")
m36 = importlib.import_module("36_canonical_phase_map")

OUT_MD = ROOT / "output" / "lmv_drift_groups.md"
OUT_RUNS = ROOT / "output" / "lmv_drift_groups_runs.csv"
OUT_EX = ROOT / "output" / "lmv_drift_groups_example.csv"
OUT_FIG = ROOT / "figures" / "fig_lmv_drift_groups.png"
N_CONS = 500
N_NAMED = len(m36.PHASES)   # seven named culture-historical phases

# Distinct grays and markers for the example-realization groups (grayscale-safe).
GRP_GRAY = ["0.15", "0.55", "0.78", "0.40", "0.90", "0.66"]
GRP_MARK = ["o", "s", "^", "D", "v", "P"]


def consensus(coords, ranks, totals, dist_km):
    """For each realization, the number of groups and between-group F_ST under
    both contemporaneous (pure-space) and time-transgressive sampling. Returns
    arrays plus the labels of one example contemporaneous realization (seed 0)."""
    nc_co, nc_tt, fst_co, fst_tt = [], [], [], []
    ex_lab_co = None
    for s in range(N_CONS):
        rng = np.random.default_rng(5000 + s)
        _, freqs = m33.drift_field_over_time(coords, seed=s, dist_km=dist_km)
        Mco = m33.sample_contemporaneous(freqs, totals, rng)
        Mtt = m33.sample_time_transgressive(freqs, ranks, totals, rng)
        lab_co, n_co, f_co = m33.fst_communities(Mco, seed=s)
        _, n_tt, f_tt = m33.fst_communities(Mtt, seed=s)
        nc_co.append(n_co); nc_tt.append(n_tt); fst_co.append(f_co); fst_tt.append(f_tt)
        if s == 0:
            ex_lab_co = lab_co
        if (s + 1) % 100 == 0:
            print(f"  ... {s + 1}/{N_CONS} runs", flush=True)
    return (np.array(nc_co, float), np.array(nc_tt, float),
            np.array(fst_co), np.array(fst_tt), ex_lab_co)


def main():
    counts_df, coords_df = m35.load_full()
    names = [str(x) for x in counts_df.index]
    counts = counts_df.to_numpy(float)
    coords = coords_df[["Latitude", "Longitude"]].to_numpy(float)
    n = counts.shape[0]
    totals = counts.sum(1).astype(int)
    lon, lat = coords[:, 1], coords[:, 0]
    gpts = gpd.GeoSeries(gpd.points_from_xy(lon, lat),
                         crs="EPSG:4326").to_crs("EPSG:26915")
    E, Nm = gpts.x.to_numpy(), gpts.y.to_numpy()

    ca1, _, _ = mf.correspondence_axis(counts)
    order = np.argsort(ca1)
    ranks = np.empty(n)
    ranks[order] = np.linspace(0, 1, n)

    # Observed structure (real data clustered the same way).
    _, obs_ncom, obs_fst = m33.fst_communities(counts)

    if OUT_RUNS.exists() and OUT_EX.exists():
        runs = pd.read_csv(OUT_RUNS)
        ex = pd.read_csv(OUT_EX)
        nc_co, nc_tt = runs["nc_co"].to_numpy(), runs["nc_tt"].to_numpy()
        fst_co, fst_tt = runs["fst_co"].to_numpy(), runs["fst_tt"].to_numpy()
        ex_lab_co = ex["ex_lab_co"].to_numpy()
        print(f"loaded cached runs from {OUT_RUNS.name}, {OUT_EX.name}")
    else:
        river_km, _, _, rinfo = mm.river_distance_matrix(coords)
        print(f"river network: largest component {rinfo['largest_component']} nodes, "
              f"max access {rinfo['max_access_km']:.1f} km")
        nc_co, nc_tt, fst_co, fst_tt, ex_lab_co = consensus(
            coords, ranks, totals, river_km)
        pd.DataFrame({"seed": np.arange(N_CONS), "nc_co": nc_co, "nc_tt": nc_tt,
                      "fst_co": fst_co, "fst_tt": fst_tt}).to_csv(OUT_RUNS, index=False)
        pd.DataFrame({"name": names, "E": E, "N": Nm,
                      "ex_lab_co": ex_lab_co}).to_csv(OUT_EX, index=False)

    # ---- figure (grayscale): two distribution panels, no map ----
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.0, 3.2))

    # Panel A: distribution of the number of groups, contemporaneous vs time-tr.
    bins = np.arange(1.5, max(nc_co.max(), nc_tt.max()) + 1.5, 1.0)
    axA.hist(nc_co, bins=bins, color="0.6", edgecolor="0.2", linewidth=0.4,
             label="space (contemp.)")
    axA.hist(nc_tt, bins=bins, color="0.25", edgecolor="black", linewidth=0.4,
             alpha=0.8, label="space + time")
    axA.axvline(N_NAMED, ls="--", c="0.1", lw=1.4,
                label=f"named phases ({N_NAMED})")
    axA.axvline(obs_ncom, ls=":", c="0.1", lw=1.4, label=f"observed ({obs_ncom})")
    axA.set_xlabel("number of groups", fontsize=8)
    axA.set_ylabel("drift runs", fontsize=8)
    axA.set_title("A. Groups under drift", fontsize=8.5)
    axA.legend(fontsize=6.0, loc="upper right", framealpha=0.9)
    axA.tick_params(labelsize=7)

    # Panel B: between-group F_ST distribution against the observed level.
    fco, ftt = fst_co[np.isfinite(fst_co)], fst_tt[np.isfinite(fst_tt)]
    axB.hist(fco, bins=20, color="0.6", edgecolor="0.2", linewidth=0.3,
             label="space (contemp.)")
    axB.hist(ftt, bins=20, color="0.25", edgecolor="black", linewidth=0.3,
             alpha=0.7, label="space + time")
    axB.axvline(obs_fst, ls="--", c="0.1", lw=1.6, label=f"observed {obs_fst:.3f}")
    axB.set_xlabel("between-group $F_{ST}$", fontsize=8)
    axB.set_ylabel("drift runs", fontsize=8)
    axB.set_title("B. Differentiation", fontsize=8.5)
    axB.legend(fontsize=6.0, loc="upper right", framealpha=0.9)
    axB.tick_params(labelsize=7)

    mf.save_all(fig, OUT_FIG)
    plt.close(fig)

    # ---- summary ----
    L = [
        f"# Phase-like groups under drift across the lower Mississippi Valley (n = {n})",
        "",
        f"Forward neutral drift on the river-network geography of all {n} Mainfort-PFG "
        f"decorated LMV assemblages, no boundary imposed, {N_CONS} realizations. Groups "
        f"are greedy-modularity communities on the Brainerd-Robinson similarity graph.",
        "",
        f"- Number of groups, contemporaneous (pure-space) sampling: mean {nc_co.mean():.2f} "
        f"(range {int(nc_co.min())}-{int(nc_co.max())}).",
        f"- Number of groups, time-transgressive sampling: mean {nc_tt.mean():.2f} "
        f"(range {int(nc_tt.min())}-{int(nc_tt.max())}); the early-to-late repertoire "
        f"turnover collapses the count toward two time-defined groups.",
        f"- Named culture-historical phases in this set: {N_NAMED}. Observed data cluster "
        f"into {obs_ncom} groups.",
        f"- Between-group cultural F_ST: observed {obs_fst:.3f}; contemporaneous drift "
        f"mean {fco.mean():.3f} (95% {np.percentile(fco, 2.5):.3f} to "
        f"{np.percentile(fco, 97.5):.3f}); time-transgressive drift mean {ftt.mean():.3f}.",
        "",
        "Interpretation: drift structured by geography already produces a few spatially "
        "coherent, phase-like groups, fewer than the seven named phases the regional "
        "scheme draws, at about the observed level of differentiation. The named phases "
        "subdivide the valley more finely than the drift structure supports.",
        "",
        f"Figure: {OUT_FIG.relative_to(ROOT)}",
    ]
    OUT_MD.parent.mkdir(exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))
    print(f"\nwrote {OUT_FIG}")


if __name__ == "__main__":
    main()

"""35_basin_pullout.py — does the Parkin phase emerge as a coherent community
from neutral drift on the wider LMV assemblage set?

Runs the time-aware, river-network drift model of analyses/33 over ALL 55
Mainfort-PFG decorated assemblages of the lower Mississippi Valley, with no
boundary imposed anywhere. Assemblages are labeled with their culture-historical
phase after Mainfort (1996; see analyses/36). Across many realizations it records,
for each assemblage, the probability that it falls in the same emergent community
as Parkin, then asks whether the Parkin-phase assemblages separate from the
other-phase (Nodena, Kent, Walls, Tipton, Jones Bayou, Parchman) assemblages, and
whether the observed Parkin-vs-others between-group F_ST exceeds the drift null.
A leave-one-out jackknife (relabel each Parkin assemblage into the comparison
group, test against a matched per-drop null) reports how much of that contrast
rests on individual assemblages. If the separation is no more than drift produces,
or collapses when single assemblages are removed, the phase distinction is a
drift-and-hydrology effect rather than a social boundary.

Per-run results cache to output/basin_pullout_runs.csv and the per-assemblage
probabilities to output/basin_pullout_prob.csv; delete those to force a rerun.

Writes output/basin_pullout.md and figures/figS8_basin_pullout.png.

Usage: PYTHONPATH=src python3 analyses/35_basin_pullout.py
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
from matplotlib.colors import Normalize, LinearSegmentedColormap  # noqa: E402
from matplotlib.cm import ScalarMappable  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
import matplotlib.patheffects as pe  # noqa: E402
import pandas as pd  # noqa: E402
import geopandas as gpd  # noqa: E402
import make_figures as mf  # noqa: E402  (house style + DECORATED_TYPES, basin members)
import make_map as mm  # noqa: E402  (river basemap + river-network distance)
m33 = importlib.import_module("33_time_aware_emergence")
m36 = importlib.import_module("36_canonical_phase_map")
from scipy.stats import mannwhitneyu  # noqa: E402
from mls_emergence.signatures.variance import cultural_fst  # noqa: E402

DATA = ROOT / "data"
OUT_MD = ROOT / "output" / "basin_pullout.md"
OUT_RUNS = ROOT / "output" / "basin_pullout_runs.csv"
OUT_PROB = ROOT / "output" / "basin_pullout_prob.csv"
OUT_FIG = ROOT / "figures" / "figS8_basin_pullout.png"
N_CONS = 500


def load_full():
    """All Mainfort-PFG decorated assemblages with coordinates, WITHOUT the
    drainage-basin restriction (the inverse of make_figures._load_curated)."""
    cur = pd.read_csv(DATA / "raw" / "mainfort-pfg-cpl.csv").dropna(subset=["Assemblages"])
    cur["Assemblages"] = cur["Assemblages"].astype(str).str.strip()
    cur = cur.drop_duplicates(subset=["Assemblages"], keep="first").set_index("Assemblages")
    type_cols = [c for c in mf.DECORATED_TYPES if c in cur.columns]
    counts = cur[type_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    counts = counts[counts.sum(axis=1) > 0]
    xy = pd.read_csv(DATA / "raw" / "mainfort-pfg-cplXY.txt", sep="\t")
    xy["Assemblages"] = xy["Assemblages"].astype(str).str.strip()
    xy = xy.drop_duplicates(subset=["Assemblages"], keep="first").set_index("Assemblages")
    coords = xy.reindex(counts.index)[["Latitude", "Longitude"]].apply(
        pd.to_numeric, errors="coerce")
    ok = coords.notna().all(axis=1)
    return counts[ok], coords[ok]


def basin_off_fst(M, is_parkin):
    """Between-group cultural F_ST for the fixed basin vs off-drainage partition."""
    groups = np.array([M[is_parkin].sum(0), M[~is_parkin].sum(0)])
    return float(cultural_fst(groups))


def consensus(coords, ranks, totals, pk, is_parkin, dist_km):
    n = coords.shape[0]
    same = np.zeros(n)
    nc, prec, rec, fst_bo = [], [], [], []
    # Leave-one-out null: for each Parkin assemblage, a null F_ST distribution
    # with that assemblage relabeled into the comparison group, on the same drift
    # realizations, so each jackknife drop is tested against its matched null.
    pk_pos = np.where(is_parkin)[0]
    jk_masks = []
    for j in pk_pos:
        mj = is_parkin.copy()
        mj[j] = False
        jk_masks.append(mj)
    jk_null = np.empty((N_CONS, len(pk_pos)))
    for s in range(N_CONS):
        rng = np.random.default_rng(3000 + s)
        _, freqs = m33.drift_field_over_time(coords, seed=s, dist_km=dist_km)
        M = m33.sample_time_transgressive(freqs, ranks, totals, rng)
        lab, ncom, _ = m33.fst_communities(M, seed=s)
        comm = lab == lab[pk]                 # Parkin's emergent community
        same += comm.astype(float)
        inter = (comm & is_parkin).sum()
        prec.append(inter / comm.sum())       # of Parkin's community, fraction in basin
        rec.append(inter / is_parkin.sum())    # of basin, fraction in Parkin's community
        nc.append(ncom)
        fst_bo.append(basin_off_fst(M, is_parkin))   # drift null for basin-vs-neighbors F_ST
        for jj, mj in enumerate(jk_masks):
            jk_null[s, jj] = basin_off_fst(M, mj)
        if (s + 1) % 100 == 0:
            print(f"  ... {s + 1}/{N_CONS} runs", flush=True)
    return (same / N_CONS, np.array(nc, float), np.array(prec), np.array(rec),
            np.array(fst_bo), pk_pos, jk_null)


def main():
    counts_df, coords_df = load_full()
    names = [str(x) for x in counts_df.index]
    counts = counts_df.to_numpy(float)
    coords = coords_df[["Latitude", "Longitude"]].to_numpy(float)
    n = counts.shape[0]
    totals = counts.sum(1).astype(int)
    is_parkin = m36.assign_phases(names, coords) == "Parkin"  # Mainfort Parkin phase
    pk = m33.parkin_index(names)
    lon, lat = coords[:, 1], coords[:, 0]

    ca1, _, _ = mf.correspondence_axis(counts)
    order = np.argsort(ca1)
    ranks = np.empty(n)
    ranks[order] = np.linspace(0, 1, n)

    gpts = gpd.GeoSeries(gpd.points_from_xy(lon, lat),
                         crs="EPSG:4326").to_crs("EPSG:26915")
    E, Nm = gpts.x.to_numpy(), gpts.y.to_numpy()

    obs_fst_bo = basin_off_fst(counts, is_parkin)   # observed basin-vs-neighbors F_ST
    pk_pos = np.where(is_parkin)[0]
    pk_names = [names[j] for j in pk_pos]
    jk_cols = [f"jk{j}" for j in range(len(pk_pos))]

    cache_ok = (OUT_RUNS.exists() and OUT_PROB.exists()
                and all(c in pd.read_csv(OUT_RUNS, nrows=1).columns for c in jk_cols))
    if cache_ok:
        runs = pd.read_csv(OUT_RUNS)
        prob = pd.read_csv(OUT_PROB)
        nc, prec, rec = runs["nc"].to_numpy(), runs["prec"].to_numpy(), runs["rec"].to_numpy()
        fst_bo = runs["fst_bo"].to_numpy()
        jk_null = runs[jk_cols].to_numpy()
        P = prob["P"].to_numpy()
        print(f"loaded cached consensus from {OUT_RUNS.name}, {OUT_PROB.name}")
    else:
        river_km, _, _, rinfo = mm.river_distance_matrix(coords)
        print(f"river network: largest component {rinfo['largest_component']} nodes, "
              f"max access {rinfo['max_access_km']:.1f} km, "
              f"{rinfo['n_unreachable']} unreachable pairs")
        P, nc, prec, rec, fst_bo, pk_pos, jk_null = consensus(
            coords, ranks, totals, pk, is_parkin, river_km)
        runs_df = pd.DataFrame({"seed": np.arange(N_CONS), "nc": nc, "prec": prec,
                                "rec": rec, "fst_bo": fst_bo})
        for jj, c in enumerate(jk_cols):
            runs_df[c] = jk_null[:, jj]
        runs_df.to_csv(OUT_RUNS, index=False)
        pd.DataFrame({"name": names, "E": E, "N": Nm, "lon": lon, "lat": lat,
                      "is_parkin": is_parkin.astype(int), "P": P}).to_csv(OUT_PROB, index=False)

    # ---- leave-one-out jackknife on the observed Parkin-vs-others F_ST ----
    # Relabel each Parkin assemblage into the comparison group; test the reduced
    # observed F_ST against its matched per-drop null (jk_null column).
    jk = []
    for jj, j in enumerate(pk_pos):
        mj = is_parkin.copy()
        mj[j] = False
        obs_j = basin_off_fst(counts, mj)
        col = jk_null[:, jj]
        col = col[np.isfinite(col)]
        jk.append((pk_names[jj], obs_j, float(np.mean(col >= obs_j)),
                   float(np.percentile(col, 97.5))))
    jk.sort(key=lambda r: r[1])
    n_jk_inside = sum(1 for r in jk if r[1] <= r[3])   # drops falling inside their null

    # ---- pull-out metrics ----
    other = np.arange(n) != pk
    P_basin = P[is_parkin & other]
    P_off = P[~is_parkin]
    u, p_mw = mannwhitneyu(P_basin, P_off, alternative="greater")
    auc = u / (len(P_basin) * len(P_off))   # P(basin member ranks above off-basin)
    nullf = fst_bo[np.isfinite(fst_bo)]
    null_p = float(np.mean(nullf >= obs_fst_bo))   # drift runs at or above observed

    # ---- figure ----
    fig = plt.figure(figsize=(7.4, 6.6))
    gsf = fig.add_gridspec(3, 2, width_ratios=[1.7, 1.0], height_ratios=[1, 1, 1],
                           hspace=0.62, wspace=0.36)
    axA = fig.add_subplot(gsf[:, 0])
    axB = fig.add_subplot(gsf[0, 1])
    axC = fig.add_subplot(gsf[1, 1])
    axD = fig.add_subplot(gsf[2, 1])

    margin = 12_000.0
    ext = (E.min() - margin, E.max() + margin, Nm.min() - margin, Nm.max() + margin)
    axA.add_patch(Rectangle((ext[0], ext[2]), ext[1] - ext[0], ext[3] - ext[2],
                            facecolor="0.93", edgecolor="none", zorder=-5))
    mm.basin_basemap(axA, ext, geology=False, grayscale=True)
    norm = Normalize(0.0, 1.0)
    cmap = LinearSegmentedColormap.from_list(
        "greys_t", plt.get_cmap("Greys")(np.linspace(0.20, 1.0, 256)))
    off = ~is_parkin
    bas = is_parkin & other
    axA.scatter(E[off], Nm[off], c=P[off], cmap=cmap, norm=norm, marker="s", s=40,
                edgecolor="black", linewidth=0.5, zorder=10, label="other phases")
    axA.scatter(E[bas], Nm[bas], c=P[bas], cmap=cmap, norm=norm, marker="o", s=66,
                edgecolor="black", linewidth=1.0, zorder=11, label="Parkin phase")
    axA.scatter([E[pk]], [Nm[pk]], marker="*", s=320, c="white",
                edgecolor="black", linewidth=0.8, zorder=12)
    axA.annotate("Parkin", (E[pk], Nm[pk]), fontsize=7, zorder=13,
                 xytext=(7, 5), textcoords="offset points",
                 path_effects=[pe.withStroke(linewidth=2, foreground="white")])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=axA,
                      orientation="horizontal", fraction=0.05, pad=0.03)
    cb.set_label("P(shares Parkin's drift-detected group)", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    axA.legend(fontsize=6.5, loc="upper right", framealpha=0.9)
    axA.set_title("A. Does the Parkin phase pull out of the wider LMV set?\n"
                  f"({N_CONS} river-network drift runs, no boundary imposed)", fontsize=8.5)

    axB.boxplot([P_basin, P_off], positions=[0, 1], widths=0.6,
                medianprops=dict(color="0.15"))
    axB.scatter(np.zeros(len(P_basin)), P_basin, s=14, c="0.15", alpha=0.55, zorder=3)
    axB.scatter(np.ones(len(P_off)), P_off, s=14, c="0.55", alpha=0.55, zorder=3)
    axB.set_xticks([0, 1])
    axB.set_xticklabels([f"Parkin\n(n={len(P_basin)})", f"others\n(n={len(P_off)})"],
                        fontsize=7.5)
    axB.set_ylabel("P(shares Parkin's community)", fontsize=8)
    axB.set_ylim(0, 1)
    axB.set_title("B. Parkin vs other phases", fontsize=8.5)
    axB.tick_params(labelsize=7)
    axB.text(0.04, 0.04, f"AUC = {auc:.2f}", transform=axB.transAxes, fontsize=8,
             bbox=dict(boxstyle="round", fc="white", ec="0.6", alpha=0.9))

    # Panel C: observed basin-vs-neighbors F_ST against the drift null
    axC.hist(nullf, bins=20, color="0.72", edgecolor="0.3", linewidth=0.3)
    axC.axvline(obs_fst_bo, ls="--", c="0.1", lw=1.6,
                label=f"observed {obs_fst_bo:.3f}\n({null_p*100:.1f}% of null ≥)")
    axC.set_xlabel("Parkin vs others $F_{ST}$", fontsize=8)
    axC.set_ylabel("drift runs", fontsize=8)
    axC.set_title("C. Observed vs drift", fontsize=8.5)
    axC.legend(fontsize=6.0, loc="upper right")
    axC.tick_params(labelsize=7)

    # Panel D: leave-one-out jackknife. Each row is the observed F_ST with one
    # Parkin assemblage relabeled into the comparison group, against the 2.5-97.5
    # band of its matched per-drop null; the top row keeps all eight. Open points
    # fall inside the band (drift-consistent); filled points stay outside it.
    rows = [("all 8", obs_fst_bo, nullf)]
    for jj, j in enumerate(pk_pos):
        mj = is_parkin.copy()
        mj[j] = False
        col = jk_null[:, jj]
        rows.append((pk_names[jj], basin_off_fst(counts, mj), col[np.isfinite(col)]))
    rows = [rows[0]] + sorted(rows[1:], key=lambda r: r[1], reverse=True)
    yy = np.arange(len(rows))[::-1]
    for y, (lab, ob, col) in zip(yy, rows):
        lo, hi = np.percentile(col, 2.5), np.percentile(col, 97.5)
        axD.plot([lo, hi], [y, y], color="0.6", lw=2.0, solid_capstyle="butt", zorder=2)
        inside = ob <= hi
        axD.scatter([ob], [y], s=26, zorder=3,
                    facecolor=("white" if inside else "0.1"),
                    edgecolor="0.1", linewidth=1.0)
    axD.set_yticks(yy)
    axD.set_yticklabels([r[0] for r in rows], fontsize=6.0)
    axD.set_xlabel("Parkin vs others $F_{ST}$", fontsize=8)
    axD.set_title("D. Leave-one-out (drop one Parkin site)", fontsize=8.5)
    axD.tick_params(labelsize=7)
    axD.margins(y=0.08)

    mf.save_all(fig, OUT_FIG)
    plt.close(fig)

    # ---- summary ----
    L = [
        f"# Does the Parkin phase pull out of the wider LMV set under drift? (n = {n})",
        "",
        f"All {n} Mainfort-PFG decorated LMV assemblages ({int(is_parkin.sum())} "
        f"Parkin-phase, {int((~is_parkin).sum())} in other phases, after Mainfort 1996); "
        f"neutral time-transgressive drift on the river network, no boundary imposed, "
        f"{N_CONS} realizations. Focal node: Parkin.",
        "",
        "## Pull-out of the Parkin phase",
        f"- Mean P(shares Parkin's community): Parkin phase {P_basin.mean():.2f} "
        f"(range {P_basin.min():.2f}-{P_basin.max():.2f}) vs other phases "
        f"{P_off.mean():.2f} (range {P_off.min():.2f}-{P_off.max():.2f}).",
        f"- Separation of Parkin from other phases by P: AUC = {auc:.2f} "
        f"(Mann-Whitney p = {p_mw:.1e}). This separation is what drift on the river "
        f"network predicts: Parkin-phase assemblages lie near Parkin and so co-occur "
        f"in its community more than distant assemblages do.",
        f"- Per run, Parkin's emergent community is {prec.mean()*100:.0f}% Parkin-phase "
        f"(precision) and captures {rec.mean()*100:.0f}% of Parkin-phase assemblages "
        f"(recall); mean {nc.mean():.1f} communities (range {int(nc.min())}-{int(nc.max())}).",
        f"- Observed Parkin-vs-others F_ST = {obs_fst_bo:.3f} against a drift null of "
        f"{nullf.mean():.3f} (95% {np.percentile(nullf, 2.5):.3f} to "
        f"{np.percentile(nullf, 97.5):.3f}, max {nullf.max():.3f}); {null_p*100:.1f}% of "
        f"drift runs reach or exceed the observed value (observed at the "
        f"{100*np.mean(nullf < obs_fst_bo):.1f}th percentile of the null).",
        f"- Leave-one-out jackknife (matched per-drop null): relabeling single Parkin "
        f"assemblages into the comparison group moves the observed F_ST to "
        f"{jk[0][1]:.3f}-{jk[-1][1]:.3f}; {len(jk) - n_jk_inside} of {len(jk)} single "
        f"drops stay outside their matched drift null, so the contrast is not driven by any "
        f"one assemblage. Dropping {jk[0][0]} gives the lowest value, F_ST = {jk[0][1]:.3f}.",
        "",
        "Interpretation: with no boundary imposed, drift on the river network already "
        "makes the Parkin-phase assemblages share Parkin's emergent community more often "
        "than the other-phase (Nodena, Kent, Walls, Tipton, Jones Bayou, Parchman) "
        "assemblages do. The pull-out (AUC) is therefore a drift-and-hydrology effect: "
        "nearby assemblages co-occur, and the co-membership is graded rather than "
        "all-or-nothing. The observed between-group F_ST sits at the extreme upper tail of "
        "the drift null and survives most single-site deletions, but it is tiny and the "
        "null's own tail reaches past it (the maximum drift value exceeds the observed). "
        "Time-averaging, which is not modeled in the null and tends to deflate between-group "
        "variance, and the small Parkin sample (n = 8) further limit the claim. We read the "
        "data as consistent with drift structured by geography, with at most a weak hint of "
        "structure beyond drift at Parkin rather than a demonstrated social boundary.",
        "",
        f"Figure: {OUT_FIG.relative_to(ROOT)}",
    ]
    OUT_MD.parent.mkdir(exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))
    print(f"\nwrote {OUT_FIG}")


if __name__ == "__main__":
    main()

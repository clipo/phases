"""Manuscript figures from the calibrated spatial pipeline (analysis 47).

Regenerates, in the July layouts and grayscale style, the figures that
analyses 25, 31, 33, 34, 35, and 37 produced from the earlier uncalibrated
drift model: figS8_within_region (Supplemental Figure S9), figS3 (basin
positive control), figS4 (basin co-membership), and figS5 (the calibration
sweep). Reads only output/revision_2026_09; run analysis 47 first.

The two wider-valley figures (fig8_lmv_drift_groups, fig9_parkin_pullout) are
no longer cited by the manuscript and are produced only under --valley; see
the note in main().

Usage: .venv/bin/python analyses/50_revision_figures.py [--valley]
"""
from pathlib import Path
import importlib
import json
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analyses"))
sys.path.insert(0, str(ROOT / "src"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Patch
import geopandas as gpd

import make_map as mm
from figstyle import panel_label, save
from mls_emergence.transmission.spatial import sample_record
rev = importlib.import_module("47_revision_analysis")
demo = importlib.import_module("25_drift_vs_groups_demo")

OUT = rev.OUT
# `output/` is gitignored, so a clean clone has none of analysis 47's products.
# The small ones are tracked with the revision they belong to; the large
# per-realization tables are not. Reading through this fallback means the
# figures that need only the small files can be regenerated from a fresh clone,
# and the ones that need the big tables say so instead of crashing on a missing
# path (rule 15: a manifest step that cannot be executed is worse than none).
RESULTS_FALLBACK = (ROOT / "docs" / "manuscript" / "revisions" / "2026-09-09"
                    / "results")
MAIN = "pooled"           # main baseline; "uniform" is the alternative null
ALT = "uniform"
N_NAMED = 7               # named late Mississippian phases in the 55-assemblage set
GRAY_MAIN, GRAY_ALT, GRAY_GROUP = "0.25", "0.62", "0.82"
CMAP = LinearSegmentedColormap.from_list(
    "greys_t", plt.get_cmap("Greys")(np.linspace(0.20, 1.0, 256)))


def _find(name):
    """The analysis-47 product `name`, from output/ or from the tracked copy."""
    for d in (OUT, RESULTS_FALLBACK):
        if (d / name).exists():
            return d / name
    return None


def _csv(name):
    p = _find(name)
    return None if p is None else pd.read_csv(p)


def load():
    sets = rev.load_sets()
    sp = _find("summary.json")
    if sp is None:
        raise SystemExit("no summary.json in output/revision_2026_09 or the "
                         "tracked results; run analyses/47_revision_analysis.py")
    s = json.loads(sp.read_text())
    return (sets, s, _csv("baseline.csv"), _csv("boundary_grid.csv"),
            _csv("co_membership.csv"), _csv("calibrated_rates.csv"))


def utm(coords):
    """Project (lat, lon) rows to UTM 15N easting/northing for the basemap."""
    pts = gpd.GeoSeries(gpd.points_from_xy(coords[:, 1], coords[:, 0]), crs="EPSG:4326").to_crs("EPSG:26915")
    return pts.x.to_numpy(), pts.y.to_numpy()


def sel(base, region, model, sampling):
    return base[(base.region == region) & (base.model == model) & (base.sampling == sampling)]


def rate_of(rates, region, model):
    r = rates[(rates.region == region) & (rates.model == model)].iloc[0]
    return dict(innovation=float(r.innovation), mixing=float(r.mixing), n_ind=int(r.n_ind))


# --------------------------------------------------------------------------- #
def fig8(sets, s, base):
    obs = s["observed"]["valley"]
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.0, 3.2))
    co, tt = sel(base, "valley", MAIN, "contemporaneous"), sel(base, "valley", MAIN, "time_transgressive")
    bins = np.arange(1.5, max(co.ncom.max(), tt.ncom.max()) + 1.5, 1.0)
    axA.hist(co.ncom, bins=bins, color="0.6", edgecolor="0.2", linewidth=0.4, label="space (contemp.)")
    axA.hist(tt.ncom, bins=bins, color="0.25", edgecolor="black", linewidth=0.4, alpha=0.8, label="space + time")
    axA.axvline(N_NAMED, ls="--", c="0.1", lw=1.4, label=f"named phases ({N_NAMED})")
    axA.axvline(obs["ncom"], ls=":", c="0.1", lw=1.4, label=f"observed ({obs['ncom']})")
    axA.set_xlabel("number of groups"); axA.set_ylabel("drift runs")
    axA.legend(fontsize=6, loc="upper center", frameon=False)
    panel_label(axA, "A")
    axB.hist(co.community_fst, bins=20, color="0.6", edgecolor="0.2", linewidth=0.3, label="space (contemp.)")
    axB.hist(tt.community_fst, bins=20, color="0.25", edgecolor="black", linewidth=0.3, alpha=0.7, label="space + time")
    axB.axvline(obs["community_fst"], ls="--", c="0.1", lw=1.6, label=f"observed {obs['community_fst']:.3f}")
    axB.set_xlabel("between-group $F_{ST}$"); axB.set_ylabel("drift runs")
    axB.legend(fontsize=6, loc="upper right", frameon=False)
    panel_label(axB, "B")
    fig.tight_layout(); save(fig, "fig8_lmv_drift_groups")
    return dict(nc_co=co.ncom.describe().to_dict(), nc_tt=tt.ncom.describe().to_dict(),
                fst_co=np.nanpercentile(co.community_fst, [2.5, 50, 97.5]).tolist(),
                fst_tt=np.nanpercentile(tt.community_fst, [2.5, 50, 97.5]).tolist(),
                p_co=rev.tail(obs["community_fst"], co.community_fst),
                p_tt=rev.tail(obs["community_fst"], tt.community_fst))


def fig9(sets, s, base, co):
    data = sets["valley"]; names = data["names"]; is_parkin = data["parkin"]
    pk = names.index("Parkin"); other = np.arange(len(names)) != pk
    E, Nm = utm(data["coords"])
    P = co[(co.region == "valley") & (co.model == MAIN)].set_index("site").loc[names, "probability"].to_numpy()
    P_basin, P_off = P[is_parkin & other], P[~is_parkin]
    u, p_mw = mannwhitneyu(P_basin, P_off, alternative="greater")
    auc = u / (len(P_basin) * len(P_off))
    null = sel(base, "valley", MAIN, "time_transgressive").parkin_fst.to_numpy()
    obs = s["observed"]["valley"]["parkin_fst"]
    p_obs = rev.tail(obs, null)
    loo = [r for r in s["leave_one_out"] if r["model"] == MAIN]

    fig = plt.figure(figsize=(7.4, 6.6))
    gsf = fig.add_gridspec(3, 2, width_ratios=[1.7, 1.0], hspace=0.62, wspace=0.36)
    axA = fig.add_subplot(gsf[:, 0]); axB = fig.add_subplot(gsf[0, 1])
    axC = fig.add_subplot(gsf[1, 1]); axD = fig.add_subplot(gsf[2, 1])

    margin = 12_000.0
    ext = (E.min() - margin, E.max() + margin, Nm.min() - margin, Nm.max() + margin)
    mm.basin_basemap(axA, ext, geology=False, grayscale=True, show_counties=False, show_states=False, draw_rivers=False)
    mm.draw_hydrorivers(axA, ext, max_ord=6, main_ord=3, grayscale=False, zorder=4.4)
    norm = Normalize(0.0, 1.0)
    off, bas = ~is_parkin, is_parkin & other
    axA.scatter(E[off], Nm[off], c=P[off], cmap=CMAP, norm=norm, marker="s", s=40,
                edgecolor="black", linewidth=0.5, zorder=10, label="other phases")
    axA.scatter(E[bas], Nm[bas], c=P[bas], cmap=CMAP, norm=norm, marker="o", s=66,
                edgecolor="black", linewidth=1.0, zorder=11, label="Parkin phase")
    axA.scatter([E[pk]], [Nm[pk]], marker="*", s=320, c="white", edgecolor="black", linewidth=0.8, zorder=12)
    axA.annotate("Parkin", (E[pk], Nm[pk]), fontsize=7, zorder=13, xytext=(7, 5), textcoords="offset points",
                 path_effects=[pe.withStroke(linewidth=2, foreground="white")])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=CMAP), ax=axA, orientation="horizontal", fraction=0.05, pad=0.03)
    cb.set_label("P(shares Parkin's drift-detected group)", fontsize=8); cb.ax.tick_params(labelsize=7)
    axA.legend(fontsize=6.5, loc="upper right", framealpha=0.9)
    panel_label(axA, "A")

    axB.boxplot([P_basin, P_off], positions=[0, 1], widths=0.6, medianprops=dict(color="0.15"))
    axB.scatter(np.zeros(len(P_basin)), P_basin, s=14, c="0.15", alpha=0.55, zorder=3)
    axB.scatter(np.ones(len(P_off)), P_off, s=14, c="0.55", alpha=0.55, zorder=3)
    axB.set_xticks([0, 1]); axB.set_xticklabels([f"Parkin\n(n={len(P_basin)})", f"others\n(n={len(P_off)})"], fontsize=7.5)
    axB.set_ylabel("P(shares Parkin's community)", fontsize=8); axB.set_ylim(0, 1); axB.tick_params(labelsize=7)
    axB.text(0.04, 0.04, f"AUC = {auc:.2f}", transform=axB.transAxes, fontsize=8,
             bbox=dict(boxstyle="round", fc="white", ec="0.6", alpha=0.9))
    panel_label(axB, "B")

    axC.hist(null, bins=20, color="0.72", edgecolor="0.3", linewidth=0.3)
    axC.axvline(obs, ls="--", c="0.1", lw=1.6, label=f"observed {obs:.3f}\n(p = {p_obs:.3f})")
    axC.set_xlabel("Parkin vs others $F_{ST}$", fontsize=8); axC.set_ylabel("drift runs", fontsize=8)
    axC.legend(fontsize=6.0, loc="upper right", frameon=False); axC.tick_params(labelsize=7)
    panel_label(axC, "C")

    rows = [("all 8", obs, np.nanpercentile(null, 2.5), np.nanpercentile(null, 97.5))]
    rows += sorted([(r["site"].replace("_", " "), r["observed"], r["lo"], r["hi"]) for r in loo],
                   key=lambda r: r[1], reverse=True)
    yy = np.arange(len(rows))[::-1]
    for y, (lab, ob, lo, hi) in zip(yy, rows):
        axD.plot([lo, hi], [y, y], color="0.6", lw=2.0, solid_capstyle="butt", zorder=2)
        inside = ob <= hi
        axD.scatter([ob], [y], s=26, zorder=3, facecolor=("white" if inside else "0.1"), edgecolor="0.1", linewidth=1.0)
    axD.set_yticks(yy); axD.set_yticklabels([r[0] for r in rows], fontsize=6.0)
    axD.set_xlabel("Parkin vs others $F_{ST}$", fontsize=8); axD.tick_params(labelsize=7); axD.margins(y=0.08)
    panel_label(axD, "D")
    save(fig, "fig9_parkin_pullout")
    n_outside = sum(1 for r in rows[1:] if r[1] > r[3])
    return dict(auc=float(auc), p_mw=float(p_mw), P_basin_mean=float(P_basin.mean()), P_off_mean=float(P_off.mean()),
                obs=obs, null_pct=np.nanpercentile(null, [2.5, 50, 97.5]).tolist(), p_upper=p_obs,
                loo_outside=n_outside, loo_total=len(rows) - 1)


def figS9(sets, s, base, grid):
    fig = plt.figure(figsize=(7.0, 5.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 0.85], hspace=0.42, wspace=0.28)
    marks = ["o", "s", "^", "D"]; grays = ["0.15", "0.55", "0.8", "0.35"]
    regions = [("basin", "LMV St. Francis basin"), ("cmv", "CMV southeast Missouri (Mississippian)")]
    for j, (name, label) in enumerate(regions):
        data = sets[name]; ax = fig.add_subplot(gs[0, j])
        xy = demo.mds2(data["m"].astype(float))
        for c in np.unique(data["labels"]):
            m = data["labels"] == c
            ax.scatter(xy[m, 0], xy[m, 1], s=24, c=grays[c % 4], marker=marks[c % 4], edgecolor="white", linewidth=0.4)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel("MDS 1", fontsize=7, labelpad=2); ax.set_ylabel("MDS 2", fontsize=7, labelpad=2)
        ax.text(0.02, 0.98, f"{label}\n(n={len(data['m'])}, {len(np.unique(data['labels']))} spatial clusters, "
                f"$F_{{ST}}$={s['observed'][name]['spatial_fst']:+.3f})", transform=ax.transAxes, fontsize=7, va="top")
        for sp in ax.spines.values():
            sp.set_edgecolor("#bbbbbb")
        panel_label(ax, "AB"[j])
    axb = fig.add_subplot(gs[1, :])
    summary = {}
    for (name, label), y in zip(regions, [1, 0]):
        obs = s["observed"][name]["spatial_fst"]
        drift = sel(base, name, MAIN, "time_transgressive").spatial_fst
        drift_alt = sel(base, name, ALT, "time_transgressive").spatial_fst
        grp = grid[(grid.region == name) & (grid.model == MAIN) & (grid.length == 24) & (grid.leak == 0.03)].spatial_fst
        for a, off, color, lw in [(drift, 0.18, GRAY_MAIN, 6), (drift_alt, 0.0, GRAY_ALT, 6), (grp, -0.18, GRAY_GROUP, 6)]:
            lo, hi = np.nanpercentile(a, [2.5, 97.5])
            axb.plot([lo, hi], [y + off, y + off], color=color, lw=lw, solid_capstyle="butt")
        axb.plot([obs], [y], "D", ms=8, color="black", zorder=5)
        axb.text(obs, y + 0.34, f"observed {obs:+.3f}", ha="center", fontsize=6.5)
        summary[name] = dict(obs=obs, drift=np.nanpercentile(drift, [2.5, 97.5]).tolist(), p_drift=rev.tail(obs, drift),
                             drift_alt=np.nanpercentile(drift_alt, [2.5, 97.5]).tolist(), p_alt=rev.tail(obs, drift_alt),
                             groups=np.nanpercentile(grp, [2.5, 97.5]).tolist(), p_groups=rev.tail(obs, grp))
    axb.set_yticks([1, 0]); axb.set_yticklabels(["LMV basin", "CMV (Miss.)"], fontsize=8); axb.set_ylim(-0.6, 1.7)
    axb.set_xlabel("between-cluster cultural $F_{ST}$")
    from matplotlib.ticker import MultipleLocator as _ML
    axb.xaxis.set_major_locator(_ML(0.01))   # round ticks only (check_figure_claims)
    axb.legend(handles=[Patch(color=GRAY_MAIN, label="calibrated drift, pooled-profile innovation (95%)"),
                        Patch(color=GRAY_ALT, label="calibrated drift, uniform innovation (95%)"),
                        Patch(color=GRAY_GROUP, label="bounded groups, multiplier 0.03 (95%)"),
                        plt.Line2D([], [], marker="D", color="black", ls="", label="observed")],
               fontsize=6.5, frameon=False, loc="upper right")
    for sp in ("top", "right", "left"):
        axb.spines[sp].set_visible(False)
    panel_label(axb, "C")
    save(fig, "figS8_within_region")
    return summary


def figS3(sets, s, base, grid, rates):
    data = sets["basin"]; obs = s["observed"]["basin"]
    keys = [("dd_r", "distance-decay r"), ("Q", "modularity Q"), ("be_spatial", "boundary excess (BR)"), ("spatial_fst", "cultural $F_{ST}$")]
    drift = sel(base, "basin", MAIN, "time_transgressive")
    grp = grid[(grid.region == "basin") & (grid.model == MAIN) & (grid.length == 24) & (grid.leak == 0.03)]
    fig = plt.figure(figsize=(7.0, 4.6))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.0, 1.1], hspace=0.5, wspace=0.55)
    summary = {}
    for j, (key, lab) in enumerate(keys):
        ax = fig.add_subplot(gs[0, j])
        ax.hist(drift[key], bins=18, color=GRAY_MAIN, alpha=0.55, density=True, label="spatial drift")
        ax.hist(grp[key], bins=12, color=GRAY_GROUP, alpha=0.8, density=True, label="bounded groups")
        ax.axvline(obs[key], color="black", lw=1.4)
        ax.set_xlabel(lab); ax.set_yticks([])
        # Round ticks only: the default locator picks 0.005 steps on these
        # ranges, which scripts/check_figure_claims.py reads as statistics
        # the text never states (its documented "finer ticks slip" limit).
        from matplotlib.ticker import MultipleLocator as _ML
        if key in ("spatial_fst", "Q"):
            ax.xaxis.set_major_locator(_ML(0.01))
        panel_label(ax, "ABCD"[j])
        summary[key] = dict(obs=obs[key], drift=np.nanpercentile(drift[key], [2.5, 97.5]).tolist(),
                            groups=np.nanpercentile(grp[key], [2.5, 97.5]).tolist(),
                            in_drift=bool(np.nanpercentile(drift[key], 2.5) <= obs[key] <= np.nanpercentile(drift[key], 97.5)),
                            in_groups=bool(np.nanpercentile(grp[key], 2.5) <= obs[key] <= np.nanpercentile(grp[key], 97.5)))
    fig.axes[0].legend(loc="upper left", fontsize=6, frameon=False)
    r = rate_of(rates, "basin", MAIN); totals = data["m"].sum(1)
    Md = sample_record(rev.simulate(data, 5001, MAIN, **r), data["ranks"], totals, np.random.default_rng(5001))
    Mg = sample_record(rev.simulate(data, 5001, MAIN, leak=0.03, **r), data["ranks"], totals, np.random.default_rng(5001))
    marks = ["o", "s", "^", "D"]; grays = ["0.15", "0.55", "0.8", "0.35"]
    for j, (ttl, M) in enumerate([("observed", data["m"].astype(float)), ("one spatial-drift realization", Md),
                                  ("one bounded-groups realization", Mg)]):
        ax = fig.add_subplot(gs[1, j])
        xy = demo.mds2(M)
        for c in np.unique(data["labels"]):
            m = data["labels"] == c
            ax.scatter(xy[m, 0], xy[m, 1], s=18, c=grays[c % 4], marker=marks[c % 4], edgecolor="white", linewidth=0.4)
        if ttl == "observed":
            pk = data["names"].index("Parkin")
            ax.scatter(xy[pk, 0], xy[pk, 1], marker="*", s=120, c="white", edgecolor="black", linewidth=0.8, zorder=5)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel("MDS 1", fontsize=7, labelpad=2); ax.set_ylabel("MDS 2", fontsize=7, labelpad=2)
        ax.text(0.02, 0.98, ttl, transform=ax.transAxes, fontsize=7, va="top")
        panel_label(ax, "EFG"[j])
    save(fig, "figS2_drift_vs_groups")
    return summary


def figS4(sets, s, base, co):
    data = sets["basin"]; names = data["names"]; pk = names.index("Parkin")
    E, Nm = utm(data["coords"])
    P = co[(co.region == "basin") & (co.model == MAIN)].set_index("site").loc[names, "probability"].to_numpy()
    fig = plt.figure(figsize=(7.2, 3.6))
    gsf = fig.add_gridspec(1, 2, width_ratios=[1.7, 1.0], wspace=0.30)
    axA = fig.add_subplot(gsf[0, 0]); axB = fig.add_subplot(gsf[0, 1])
    margin = 8_000.0
    ext = (E.min() - margin, E.max() + margin, Nm.min() - margin, Nm.max() + margin)
    mm.basin_basemap(axA, ext, geology=False, grayscale=True, show_counties=False, show_states=False, draw_rivers=False)
    mm.draw_hydrorivers(axA, ext, max_ord=6, main_ord=3, grayscale=False, zorder=4.4)
    norm = Normalize(0.0, 1.0); nonpk = np.arange(len(names)) != pk
    axA.scatter(E[nonpk], Nm[nonpk], c=P[nonpk], cmap=CMAP, norm=norm, s=54, edgecolor="black", linewidth=0.6, zorder=10)
    axA.scatter([E[pk]], [Nm[pk]], marker="*", s=320, c="white", edgecolor="black", linewidth=0.8, zorder=12)
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=CMAP), ax=axA, orientation="horizontal", fraction=0.05, pad=0.03)
    cb.set_label("P(shares Parkin's drift-detected group)", fontsize=8); cb.ax.tick_params(labelsize=7)
    panel_label(axA, "A")
    d_pk = data["d"][pk][nonpk]
    from scipy.stats import spearmanr
    rho = spearmanr(d_pk, P[nonpk]).correlation
    obs = s["observed"]["basin"]["community_fst"]
    x = 0
    for model, color in [(MAIN, GRAY_MAIN), (ALT, GRAY_ALT)]:
        for kind, mk in [("contemporaneous", "o"), ("time_transgressive", "s")]:
            a = sel(base, "basin", model, kind).community_fst
            lo, med, hi = np.nanpercentile(a, [2.5, 50, 97.5])
            axB.errorbar(x, med, yerr=[[med - lo], [hi - med]], fmt=mk, color=color, ms=5, capsize=3, lw=1,
                         label=f"{rev.MODEL_LABEL[model]}, {kind.replace('_', ' ')}")
            x += 1
    axB.axhline(obs, ls="--", color="black", lw=1.2, label=f"observed {obs:.3f}")
    axB.set_xticks(range(4)); axB.set_xticklabels(["pooled\nspace", "pooled\nspace+time", "uniform\nspace", "uniform\nspace+time"], fontsize=6.5)
    axB.set_ylabel("between-group $F_{ST}$", fontsize=8); axB.tick_params(labelsize=7)
    axB.legend(fontsize=5.5, frameon=False, loc="upper right")
    panel_label(axB, "B")
    save(fig, "figS3_emergent_phases")
    return dict(rho_distance=float(rho), P_range=[float(P[nonpk].min()), float(P[nonpk].max())], obs=obs,
                p_co=rev.tail(obs, sel(base, "basin", MAIN, "contemporaneous").community_fst),
                p_tt=rev.tail(obs, sel(base, "basin", MAIN, "time_transgressive").community_fst))


def figS5(sets, s, base):
    """Calibration sweep: every grid cell's partition F_ST against its diversity misfit."""
    calib = _csv("calibration.csv")
    # Two panels, not three. The 55-assemblage wider-valley set was dropped from
    # the paper: it mixes deposits of different ages and very different sample
    # sizes, so a partition F_ST measured across it confounds chronology with
    # spatial process. The St. Francis basin is the test and southeast Missouri
    # is the comparison.
    fig, axes = plt.subplots(1, 2, figsize=(5.0, 2.9), sharey=False)
    for ax, (region, metric, label) in zip(axes, [("basin", "spatial_fst", "St. Francis basin"),
                                                  ("cmv", "spatial_fst", "SE Missouri")]):
        obs = s["observed"][region][metric]
        ax.axhline(obs, ls="--", color="black", lw=1.2, label="observed")
        for model, color, mk in [(MAIN, GRAY_MAIN, "o"), (ALT, GRAY_ALT, "s")]:
            f = calib[(calib.region == region) & (calib.model == model)]
            misfit = np.sqrt(f.loss)
            ax.scatter(misfit[~f.matched], f.fst_median[~f.matched], s=9, marker=mk, facecolor="none",
                       edgecolor=color, linewidth=0.6, label=f"{rev.MODEL_LABEL[model]}, unmatched")
            # Matched cells are drawn at their stage-2 (fifty-seed) median, the
            # value the selection is made on; unmatched cells have only the screen.
            y2 = f.stage2_fst_median.where(f.stage2_fst_median.notna(), f.fst_median)
            ax.scatter(misfit[f.matched], y2[f.matched], s=16, marker=mk, color=color,
                       label=f"{rev.MODEL_LABEL[model]}, diversity-matched (50 seeds)")
        ax.axvline(np.sqrt(3) * rev.CONFIG["tolerance"], color="0.7", lw=0.8)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("diversity misfit (root summed squared relative error)", fontsize=7)
        ax.set_ylabel(r"partition $F_{ST}$ (median; 6 seeds unmatched, 50 matched)", fontsize=7)
        ax.text(0.03, 0.03, label, transform=ax.transAxes, fontsize=7.5)
        panel_label(ax, "AB"[list(axes).index(ax)])
    axes[0].legend(fontsize=5, frameon=False, loc="upper right")
    fig.tight_layout(); save(fig, "figS4_emergence_robustness")
    out = {}
    for region, metric in [("basin", "spatial_fst"), ("cmv", "spatial_fst")]:
        obs = s["observed"][region][metric]
        for model in (MAIN, ALT):
            f = calib[(calib.region == region) & (calib.model == model) & calib.matched]
            f2 = f[f.stage2_matched.astype(bool)]
            out[f"{region}/{model}"] = dict(n_matched=int(len(f)), n_matched_stage2=int(len(f2)),
                                             max_fst_stage2=float(f2.stage2_fst_median.max()) if len(f2) else None,
                                             n_reaching_observed_stage2=int((f2.stage2_fst_median >= obs).sum()),
                                             n_covering_observed_stage2=int(((f2.stage2_lo <= obs) & (obs <= f2.stage2_hi)).sum()))
    return out


def main():
    sets, s, base, grid, co, rates = load()
    out = {}
    # Each figure states the per-realization tables it needs, and is skipped
    # with a message when they are absent, rather than failing the whole run.
    todo = [("figS9", lambda: figS9(sets, s, base, grid), (base, grid)),
            ("figS3", lambda: figS3(sets, s, base, grid, rates), (base, grid, rates)),
            ("figS4", lambda: figS4(sets, s, base, co), (base, co)),
            ("figS5", lambda: figS5(sets, s, base), ())]
    for name, fn, needs in todo:
        if any(x is None for x in needs):
            print(f"skipping {name}: needs analysis-47 tables not present "
                  f"(run analyses/47_revision_analysis.py)")
            continue
        out[name] = fn()
    # The two wider-valley figures are off by default. The 55-assemblage valley
    # set was dropped from the paper because it mixes deposits of different ages
    # and very different sample sizes, so a partition F_ST measured across it
    # confounds chronology with spatial process. The code is kept, behind a
    # flag, because the set itself is real and the decision is editorial: a
    # future analysis that controls for age could use it. Nothing in the
    # manuscript cites either figure.
    if "--valley" in sys.argv:
        out.update(fig8=fig8(sets, s, base), fig9=fig9(sets, s, base, co))
    # MERGE, do not clobber. When some figures were skipped for want of the
    # per-realization tables, overwriting this file would silently delete the
    # values of the figures that did not run this time.
    OUT.mkdir(parents=True, exist_ok=True)
    fv = OUT / "figure_values.json"
    prev = json.loads(fv.read_text()) if fv.exists() else {}
    if not prev and (RESULTS_FALLBACK / "figure_values.json").exists():
        prev = json.loads((RESULTS_FALLBACK / "figure_values.json").read_text())
    # Merge, but only for figures this script can still produce. A key it no
    # longer emits at all is an orphan: nothing recomputes it, so it survives
    # every rerun and reads as current. `fig11` sat here at the 29-assemblage
    # F_ST for a day after the set changed to 30, and was quoted as though it
    # were fresh.
    producible = {name for name, _, _ in todo} | {"fig8", "fig9", "_written"}
    orphans = sorted(set(prev) - producible)
    for k in orphans:
        prev.pop(k)
    if orphans:
        print(f"dropped orphaned figure values, no longer produced here: {orphans}")
    prev.update(out)
    # Keep each block's own write time, so a value carried over from an earlier
    # run is legible as carried over rather than passing for fresh.
    from datetime import datetime
    stamps = dict(prev.get("_written", {}))
    stamps.update({k: datetime.now().isoformat(timespec="seconds") for k in out})
    prev["_written"] = {k: v for k, v in stamps.items() if k in prev}
    fv.write_text(json.dumps(prev, indent=2, default=float))
    print(json.dumps(out, indent=2, default=float))


if __name__ == "__main__":
    main()

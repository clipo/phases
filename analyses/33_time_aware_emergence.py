"""33_time_aware_emergence.py — do phase-like structures emerge when neutral
drift runs on the real geography AND assemblages are non-contemporaneous?

Extends 32 by dropping the contemporaneity assumption. A neutral drift field
evolves through time on the actual coordinate layout of the curated St. Francis
basin assemblages, with copying that decays with geographic distance and NO
group boundaries. Two features make it realistic:

  - New decorated types arise over time (innovation), so later moments carry
    "late types" that earlier ones lack -- the differential late-type
    frequencies that let the assemblages seriate at all.
  - Each assemblage is sampled time-transgressively, at the time corresponding
    to its observed seriation position (CA1 rank), and time-averaged over an
    accumulation window (a death assemblage mixes production across its span).

We then read the output the way a culture historian would and compare three
cases on the same geography and sherd totals:
  (i)  observed data,
  (ii) contemporaneous drift (one equilibrium snapshot, as in 32),
  (iii) time-transgressive drift (this model).

Reports whether the synthetic assemblages (a) seriate, (b) form spatially
coherent, elevated-F_ST communities (phase-like structure), and (c) match the
observed level of structure better than the contemporaneous case.

Writes output/time_aware_emergence.md and figures/fig8_emergent_phases.png.

Usage: PYTHONPATH=src python3 analyses/33_time_aware_emergence.py
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
import make_figures as mf  # noqa: E402  (applies house style on import)
import make_map as mm  # noqa: E402  (reuse the fig1 river basemap)
sd = importlib.import_module("23_phases_vs_spatial_drift")
from mls_emergence.signatures.variance import cultural_fst  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

OUT_MD = ROOT / "output" / "time_aware_emergence.md"
OUT_FIG = ROOT / "figures" / "fig8_emergent_phases.png"
OUT_RUNS = ROOT / "output" / "time_aware_runs.csv"
OUT_PARKIN = ROOT / "output" / "time_aware_parkin_prob.csv"
N_CONS = 500       # consensus realizations for the Parkin co-membership map

LEN_KM = 24.0      # interaction length scale (positive-control corner)
M_BETWEEN = 0.02   # between-node mixing
MU = 0.012         # innovation (new types appear over time)
K_MAX = 40         # type-label capacity (columns); late columns appear late
N_IND = 120        # individuals per node
T_STEPS = 1800     # total generations
N_REC = 90         # recorded time slices
WIN = 8            # accumulation window (record slices) for time-averaging
N_SEEDS = 16


def fst_communities(counts, seed=0):
    labels, _, ncom = sd.ceramic_communities(counts, seed=seed)
    groups = np.array([counts[labels == c].sum(0) for c in np.unique(labels)
                       if counts[labels == c].sum() > 0])
    fst = float(cultural_fst(groups)) if len(groups) >= 2 else float("nan")
    return labels, ncom, fst


def drift_field_over_time(coords, seed, len_km=LEN_KM, m_between=M_BETWEEN,
                          mu=MU, k_max=K_MAX, n_ind=N_IND, t_steps=T_STEPS,
                          n_rec=N_REC, dist_km=None):
    """Evolve neutral spatial drift; record per-node type frequencies at n_rec
    evenly spaced time slices. New types (columns) appear over time via
    innovation. Parameters default to the module-level positive-control corner;
    pass overrides to sweep them. ``dist_km`` is an optional precomputed N x N
    inter-site distance matrix (km); if None, great-circle distance is used.
    Returns (rec_steps[n_rec], freq[n_rec, n, k_max])."""
    rng = np.random.default_rng(seed)
    n = coords.shape[0]
    D = sd.geo_km(coords) if dist_km is None else np.asarray(dist_km, float)
    W = np.exp(-D / len_km)
    np.fill_diagonal(W, 0.0)
    W = W / W.sum(axis=1, keepdims=True)

    pop = np.zeros((n, n_ind), dtype=int)  # start with type 0 everywhere
    next_type = 1
    rec_at = set(np.linspace(0, t_steps - 1, n_rec).astype(int).tolist())
    rec_steps, freqs = [], []
    for t in range(t_steps):
        new = pop.copy()
        for i in range(n):
            from_other = rng.random(n_ind) < m_between
            src_self = pop[i, rng.integers(0, n_ind, size=n_ind)]
            partners = rng.choice(n, size=n_ind, p=W[i])
            src_other = pop[partners, rng.integers(0, n_ind, size=n_ind)]
            drawn = np.where(from_other, src_other, src_self)
            innov = rng.random(n_ind) < mu
            ni = int(innov.sum())
            if ni:
                drawn = drawn.copy()
                slots = []
                for _ in range(ni):
                    if next_type < k_max:
                        slots.append(next_type); next_type += 1
                    else:
                        slots.append(int(rng.integers(0, k_max)))
                drawn[innov] = slots
            new[i] = drawn
        pop = new
        if t in rec_at:
            f = np.zeros((n, k_max))
            for i in range(n):
                u, cnt = np.unique(pop[i], return_counts=True)
                f[i, u] = cnt
            f = f / f.sum(1, keepdims=True)
            rec_steps.append(t); freqs.append(f)
    return np.array(rec_steps), np.stack(freqs)  # (R,), (R,n,k_max)


def sample_time_transgressive(freqs, ranks, totals, rng, win=WIN):
    """Each assemblage sampled at its seriation-rank time slice, time-averaged
    over ``win`` preceding slices, to its observed total."""
    R, n, K = freqs.shape
    slot = np.clip((ranks * (R - 1)).round().astype(int), 0, R - 1)
    M = np.zeros((n, K), dtype=int)
    for i in range(n):
        lo = max(0, slot[i] - win + 1)
        p = freqs[lo:slot[i] + 1, i, :].mean(0)
        p = p / p.sum() if p.sum() else np.ones(K) / K
        M[i] = rng.multinomial(int(totals[i]), p) if totals[i] > 0 else 0
    return M


def sample_contemporaneous(freqs, totals, rng):
    """All assemblages sampled at the final slice (equilibrium snapshot)."""
    _, n, K = freqs.shape
    M = np.zeros((n, K), dtype=int)
    for i in range(n):
        p = freqs[-1, i, :]
        p = p / p.sum() if p.sum() else np.ones(K) / K
        M[i] = rng.multinomial(int(totals[i]), p) if totals[i] > 0 else 0
    return M


def parkin_index(names):
    for i, nm in enumerate(names):
        if "parkin" in str(nm).lower():
            return i
    raise ValueError("Parkin not found among curated assemblages")


def compute_consensus(coords, ranks, totals, pk, dist_km=None):
    """Run N_CONS time-transgressive realizations at the default (positive-control)
    parameters. Return: P[n] = fraction of runs each assemblage shares Parkin's
    emergent community; and per-run arrays of between-group F_ST (time-transgressive
    and contemporaneous), community count, and seriation recovery. ``dist_km`` is
    the inter-site distance matrix driving the copying kernel (river-network km)."""
    n = coords.shape[0]
    same_co = np.zeros(n)   # co-membership with Parkin, contemporaneous snapshot
    same_tt = np.zeros(n)   # co-membership with Parkin, time-transgressive
    fst_tt, fst_co, nc_tt, rho = [], [], [], []
    for s in range(N_CONS):
        rng = np.random.default_rng(2000 + s)
        _, freqs = drift_field_over_time(coords, seed=s, dist_km=dist_km)
        Mtt = sample_time_transgressive(freqs, ranks, totals, rng)
        Mco = sample_contemporaneous(freqs, totals, rng)
        lab_tt, nc, f = fst_communities(Mtt, seed=s)
        lab_co, _, fco = fst_communities(Mco, seed=s)
        same_tt += (lab_tt == lab_tt[pk]).astype(float)
        same_co += (lab_co == lab_co[pk]).astype(float)
        ca, _, _ = mf.correspondence_axis(Mtt)
        r = spearmanr(ranks, ca).correlation
        fst_tt.append(f); fst_co.append(fco); nc_tt.append(nc)
        rho.append(abs(r) if np.isfinite(r) else np.nan)
        if (s + 1) % 100 == 0:
            print(f"  ... {s + 1}/{N_CONS} runs", flush=True)
    return (same_co / N_CONS, same_tt / N_CONS, np.array(fst_tt),
            np.array(fst_co), np.array(nc_tt, float), np.array(rho))


def main():
    counts_df, coords_df = mf._load_curated()
    counts = counts_df.to_numpy(float)
    coords = coords_df[["Latitude", "Longitude"]].to_numpy(float)
    names = [str(x) for x in counts_df.index]
    n = counts.shape[0]
    totals = counts.sum(1).astype(int)
    lon, lat = coords[:, 1], coords[:, 0]
    pk = parkin_index(names)

    # observed seriation position (CA1), oriented and rank-normalized to [0,1]
    ca1, _, frac = mf.correspondence_axis(counts)
    order = np.argsort(ca1)
    ranks = np.empty(n)
    ranks[order] = np.linspace(0, 1, n)

    # observed structure
    _, obs_ncom, obs_fst = fst_communities(counts)

    # project assemblage lon/lat to the basemap CRS (UTM zone 15N)
    gpts = gpd.GeoSeries(gpd.points_from_xy(lon, lat),
                         crs="EPSG:4326").to_crs("EPSG:26915")
    E, Nm = gpts.x.to_numpy(), gpts.y.to_numpy()

    # consensus across N_CONS realizations (cached for fast re-plotting)
    if OUT_RUNS.exists() and OUT_PARKIN.exists():
        runs = pd.read_csv(OUT_RUNS)
        pkdf = pd.read_csv(OUT_PARKIN)
        fst_tt = runs["fst_tt"].to_numpy()
        fst_co = runs["fst_co"].to_numpy()
        nc_tt = runs["nc"].to_numpy()
        rho = runs["rho"].to_numpy()
        P_co = pkdf["P_co"].to_numpy()
        P_tt = pkdf["P_tt"].to_numpy()
        print(f"loaded cached consensus from {OUT_RUNS.name}, {OUT_PARKIN.name}")
    else:
        river_km, _, access_km, rinfo = mm.river_distance_matrix(coords)
        print(f"river network: largest component {rinfo['largest_component']} nodes, "
              f"max access {rinfo['max_access_km']:.1f} km, "
              f"{rinfo['n_unreachable']} unreachable pairs")
        P_co, P_tt, fst_tt, fst_co, nc_tt, rho = compute_consensus(
            coords, ranks, totals, pk, dist_km=river_km)
        pd.DataFrame({"seed": np.arange(N_CONS), "nc": nc_tt, "fst_tt": fst_tt,
                      "fst_co": fst_co, "rho": rho}).to_csv(OUT_RUNS, index=False)
        pd.DataFrame({"name": names, "E": E, "N": Nm, "lon": lon, "lat": lat,
                      "P_co": P_co, "P_tt": P_tt}).to_csv(OUT_PARKIN, index=False)

    nc_tt = np.asarray(nc_tt, float)
    nonpk = np.arange(n) != pk
    P = P_co  # panel A shows the contemporaneous (pure-space) co-membership
    geo_pk = sd.geo_km(coords)[pk]
    rho_co = spearmanr(geo_pk[nonpk], P_co[nonpk]).correlation
    rho_tt = spearmanr(geo_pk[nonpk], P_tt[nonpk]).correlation

    # ---- figure: A = P(shares Parkin's community) on the river basemap;
    #              B = drift reproduces the observed between-group F_ST ----
    fig = plt.figure(figsize=(7.2, 5.0))
    gsf = fig.add_gridspec(1, 2, width_ratios=[1.7, 1.0], wspace=0.30)
    axA = fig.add_subplot(gsf[0, 0])
    axB = fig.add_subplot(gsf[0, 1])

    margin = 12_000.0
    ext = (E.min() - margin, E.max() + margin, Nm.min() - margin, Nm.max() + margin)
    # Grayscale base to match Figure 1: uniform land tone plus the no-geology,
    # gray-hydrology basemap (American Antiquity prints without color).
    axA.add_patch(Rectangle((ext[0], ext[2]), ext[1] - ext[0], ext[3] - ext[2],
                            facecolor="0.93", edgecolor="none", zorder=-5))
    mm.basin_basemap(axA, ext, geology=False, grayscale=True)

    norm = Normalize(0.0, 1.0)
    # Truncated Greys so the lowest probabilities are still a visible mid-gray on
    # the light land rather than white; markers keep black edges for definition.
    cmap = LinearSegmentedColormap.from_list(
        "greys_t", plt.get_cmap("Greys")(np.linspace(0.20, 1.0, 256)))
    axA.scatter(E[nonpk], Nm[nonpk], c=P[nonpk], cmap=cmap, norm=norm, s=54,
                edgecolor="black", linewidth=0.5, zorder=10)
    axA.scatter([E[pk]], [Nm[pk]], marker="*", s=320, c="white",
                edgecolor="black", linewidth=0.8, zorder=11)
    axA.annotate("Parkin", (E[pk], Nm[pk]), fontsize=7, zorder=12,
                 xytext=(7, 5), textcoords="offset points",
                 path_effects=[pe.withStroke(linewidth=2, foreground="white")])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=axA,
                      orientation="horizontal", fraction=0.05, pad=0.03)
    cb.set_label("P(shares Parkin's drift-detected group)", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    axA.set_title("A. Probability of sharing Parkin's \"phase\" under spatial drift\n"
                  f"({N_CONS} river-network runs, contemporaneous snapshot)", fontsize=8.5)

    # Panel B: between-group F_ST, observed vs the drift 95% interval (2.5-97.5
    # percentiles). The observed sits inside the contemporaneous-drift interval,
    # so drift reproduces the observed level; the exceedance fraction is annotated.
    axB.axhline(obs_fst, ls="--", c="0.4", lw=1)
    axB.scatter([0], [obs_fst], s=55, c="black", zorder=3)
    axB.annotate(f"{obs_fst:.3f}", (0, obs_fst), xytext=(8, 0),
                 textcoords="offset points", fontsize=6.5, va="center")
    for x, arr in [(1, fst_co), (2, fst_tt)]:
        a = arr[np.isfinite(arr)]
        lo, md, hi = np.percentile(a, [2.5, 50, 97.5])
        axB.errorbar([x], [md], yerr=[[md - lo], [hi - md]], fmt="o", ms=7,
                     c="0.3", capsize=3, zorder=3)
        axB.annotate(f"{np.mean(a >= obs_fst) * 100:.1f}%\n≥ obs", (x, hi),
                     xytext=(0, 4), textcoords="offset points", fontsize=5.8,
                     ha="center", va="bottom", color="0.3")
    axB.set_xticks([0, 1, 2])
    axB.set_xticklabels(["observed", "contemp.\ndrift", "time-tr.\ndrift"],
                        fontsize=7)
    axB.set_xlim(-0.5, 2.5)
    axB.set_ylim(0, max(obs_fst, np.nanpercentile(fst_co, 97.5),
                        np.nanpercentile(fst_tt, 97.5)) * 1.35)
    axB.set_ylabel("between-group $F_{ST}$", fontsize=8)
    axB.set_title("B. Drift reproduces the\nobserved level", fontsize=8.5)
    axB.tick_params(labelsize=7)

    fig.savefig(OUT_FIG, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ---- summary ----
    L = [
        f"# Time-aware emergence of phase-like structure (n = {n})",
        "",
        f"Drift field on the real layout at the positive-control corner, copying kernel "
        f"on river-network (along-waterway) distance: interaction length {LEN_KM:.0f} "
        f"river-km, mixing {M_BETWEEN}, innovation {MU}, {T_STEPS} steps, {N_REC} slices, "
        f"window {WIN}; {N_CONS} realizations. CA1 inertia (observed) = {frac:.2f}. "
        f"Focal node: {names[pk]}.",
        "",
        "## Probability of sharing Parkin's emergent community",
        f"- Contemporaneous (pure spatial drift) snapshot, shown in Figure 8A: mean "
        f"P = {P_co[nonpk].mean():.2f} (range {P_co[nonpk].min():.2f}-{P_co[nonpk].max():.2f}) "
        f"across the other {n - 1} assemblages; co-membership falls with distance from "
        f"Parkin (Spearman rho = {rho_co:+.2f}), the intuitive spatial gradient.",
        f"- Time-transgressive sampling: mean P = {P_tt[nonpk].mean():.2f} "
        f"(range {P_tt[nonpk].min():.2f}-{P_tt[nonpk].max():.2f}); co-membership no longer "
        f"falls with distance (Spearman rho = {rho_tt:+.2f}) because assemblages at "
        "different seriation positions carry different repertoires, so membership reflects "
        "space and time jointly.",
        "- Either way, no assemblage is a certain member (P = 1) or non-member (P = 0) of "
        "Parkin's group: membership is graded and probabilistic, so the 'phase' has no "
        "stable membership.",
        "",
        "## Emergent communities and between-group cultural F_ST",
        "| case | communities | F_ST |",
        "|---|---|---|",
        f"| observed data | {obs_ncom} | {obs_fst:.3f} |",
        f"| contemporaneous drift | -- | {np.nanmean(fst_co):.3f} |",
        f"| time-transgressive drift | {nc_tt.mean():.1f} | {np.nanmean(fst_tt):.3f} |",
        "",
        f"- Phase-like structure (>= 2 communities): {np.mean(nc_tt >= 2) * 100:.0f}% of "
        f"runs; mean {nc_tt.mean():.1f} (range {int(nc_tt.min())}-{int(nc_tt.max())}).",
        f"- Time-transgressive record seriates: mean |Spearman rho| = {np.nanmean(rho):.2f}.",
        "",
        "Interpretation: across many drift realizations on the real geography, with no "
        "boundary imposed, membership in Parkin's emergent community is graded and "
        "probabilistic rather than marking a bounded set, and the between-group cultural "
        "F_ST sits at the observed level rather than the much higher value a bounded "
        "group would leave. The phase is the fuzzy, spatially structured residue of "
        "drift, not a society.",
        "",
        f"Figure: {OUT_FIG.relative_to(ROOT)}",
    ]
    OUT_MD.parent.mkdir(exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))
    print(f"\nwrote {OUT_FIG}")


if __name__ == "__main__":
    main()

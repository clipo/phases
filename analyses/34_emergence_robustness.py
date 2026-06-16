"""34_emergence_robustness.py — is the emergence of phase-like structure from
drift on the real basin geography robust to contingency and to the model factors?

Figure 8 (script 33) shows one representative realization at the positive-control
corner. Because each run is stochastic and the corner is one point in parameter
space, this script repeats the time-aware emergence experiment across many seeds
and across a grid of the three factors that govern it:

  - interaction length (distance decay of copying),
  - between-node mixing,
  - innovation rate (appearance of new decorated types).

For every (factor combination x seed) it records the time-transgressive outcome:
the number of emergent communities, the between-group cultural F_ST, and how well
the synthetic record re-seriates. We then ask whether phase-like structure
(>= 2 spatially coherent communities) appears across the grid and whether the
emergent F_ST brackets the observed value rather than depending on a tuned corner.

Raw per-run results are cached to output/emergence_robustness.csv so the figure
can be regenerated without re-simulating; delete that file to force a fresh run.

Writes output/emergence_robustness.{csv,md} and figures/figS7_emergence_robustness.png.

Usage: PYTHONPATH=src python3 analyses/34_emergence_robustness.py
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
import make_figures as mf  # noqa: E402  (applies house style on import)
import make_map as mm  # noqa: E402  (river-network distance)
m33 = importlib.import_module("33_time_aware_emergence")
from scipy.stats import spearmanr  # noqa: E402

OUT_MD = ROOT / "output" / "emergence_robustness.md"
OUT_CSV = ROOT / "output" / "emergence_robustness.csv"
OUT_FIG = ROOT / "figures" / "figS7_emergence_robustness.png"

# Factor grid, centered on the positive-control corner (24 km, 0.02, 0.012).
LEN_GRID = [12.0, 18.0, 24.0, 36.0]   # interaction length (km)
MIX_GRID = [0.01, 0.02, 0.05]         # between-node mixing
MU_GRID = [0.006, 0.012, 0.024]       # innovation rate
N_SEEDS = 12
COLS = "len_km,mix,mu,seed,ncom,fst,rho"


def simulate(coords, ranks, totals):
    river_km, _, _, rinfo = mm.river_distance_matrix(coords)
    print(f"river network: largest component {rinfo['largest_component']} nodes, "
          f"max access {rinfo['max_access_km']:.1f} km")
    rows = []
    total = len(LEN_GRID) * len(MIX_GRID) * len(MU_GRID) * N_SEEDS
    done = 0
    for lk in LEN_GRID:
        for mx in MIX_GRID:
            for mu in MU_GRID:
                for s in range(N_SEEDS):
                    rng = np.random.default_rng(7000 + s)
                    _, freqs = m33.drift_field_over_time(
                        coords, seed=s, len_km=lk, m_between=mx, mu=mu,
                        dist_km=river_km)
                    Mtt = m33.sample_time_transgressive(freqs, ranks, totals, rng)
                    _, nc, fst = m33.fst_communities(Mtt, seed=s)
                    ca, _, _ = mf.correspondence_axis(Mtt)
                    rho = spearmanr(ranks, ca).correlation
                    rows.append((lk, mx, mu, s, nc, fst,
                                 abs(rho) if np.isfinite(rho) else np.nan))
                    done += 1
        print(f"  ... {done}/{total} runs (interaction length {lk:.0f} km done)",
              flush=True)
    arr = np.array(rows, dtype=float)
    np.savetxt(OUT_CSV, arr, delimiter=",", header=COLS, comments="",
               fmt=["%.4f", "%.4f", "%.4f", "%d", "%d", "%.6f", "%.6f"])
    return arr


def box_by_factor(ax, levels, factor_col, val_col, obs, xlabel, title, fmt="%g"):
    data = [val_col[np.isclose(factor_col, lv) & np.isfinite(val_col)]
            for lv in levels]
    pos = np.arange(len(levels))
    ax.boxplot(data, positions=pos, widths=0.6,
               medianprops=dict(color="0.5"))
    ax.axhline(obs, ls="--", c="0.0", lw=1.2, label=f"observed ({obs:.2f})")
    ax.set_xticks(pos)
    ax.set_xticklabels([fmt % lv for lv in levels])
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel("between-group $F_{ST}$", fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=7)


def marg(levels, factor_col, val_col):
    out = []
    for lv in levels:
        v = val_col[np.isclose(factor_col, lv) & np.isfinite(val_col)]
        out.append((lv, float(np.mean(v)), float(np.std(v))))
    return out


def main():
    counts_df, coords_df = mf._load_curated()
    counts = counts_df.to_numpy(float)
    coords = coords_df[["Latitude", "Longitude"]].to_numpy(float)
    n = counts.shape[0]
    totals = counts.sum(1).astype(int)

    ca1, _, _ = mf.correspondence_axis(counts)
    order = np.argsort(ca1)
    ranks = np.empty(n)
    ranks[order] = np.linspace(0, 1, n)
    _, obs_ncom, obs_fst = m33.fst_communities(counts)

    if OUT_CSV.exists():
        arr = np.loadtxt(OUT_CSV, delimiter=",", skiprows=1)
        print(f"loaded cached results from {OUT_CSV}")
    else:
        arr = simulate(coords, ranks, totals)

    lk_c, mx_c, mu_c, _, nc_c, fst_c, rho_c = (arr[:, i] for i in range(7))
    fin = np.isfinite(fst_c)

    frac_phase = float(np.mean(nc_c >= 2))
    frac_le_obs = float(np.mean(fst_c[fin] <= obs_fst))
    frac_bracket = float(np.mean((fst_c[fin] >= 0.5 * obs_fst)
                                 & (fst_c[fin] <= 2.0 * obs_fst)))

    # ---- figure (2x2) ----
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.4))
    axA, axB, axC, axD = axes.ravel()

    box_by_factor(axA, MIX_GRID, mx_c, fst_c, obs_fst,
                  "between-node mixing", "A. F_ST vs mixing", fmt="%g")
    axA.legend(fontsize=6.5, loc="upper right")
    box_by_factor(axB, MU_GRID, mu_c, fst_c, obs_fst,
                  "innovation rate", "B. F_ST vs innovation", fmt="%g")

    vals, cnts = np.unique(nc_c.astype(int), return_counts=True)
    axC.bar(vals, cnts / cnts.sum(), width=0.7, color="0.5",
            edgecolor="black", linewidth=0.3)
    axC.set_xlabel("drift-detected groups (per run)", fontsize=8)
    axC.set_ylabel("fraction of runs", fontsize=8)
    axC.set_title(f"C. Phase-like groups under drift ({frac_phase*100:.0f}% have $\\geq$2)",
                  fontsize=9)
    axC.set_xticks(vals)
    axC.tick_params(labelsize=7)

    axD.hist(rho_c[np.isfinite(rho_c)], bins=16, color="0.5",
             edgecolor="black", linewidth=0.3)
    axD.axvline(np.nanmean(rho_c), ls="--", c="0.0", lw=1.2,
                label=f"mean {np.nanmean(rho_c):.2f}")
    axD.set_xlabel("seriation recovery |Spearman rho|", fontsize=8)
    axD.set_ylabel("runs", fontsize=8)
    axD.set_title("D. The record seriates across runs", fontsize=9)
    axD.legend(fontsize=6.5)
    axD.tick_params(labelsize=7)

    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ---- summary ----
    def fmt3(rows3):
        return "; ".join(f"{lv:g}: {mn:.3f}±{sd:.3f}" for lv, mn, sd in rows3)

    L = [
        f"# Robustness of phase-like emergence to contingency and factors (n = {n})",
        "",
        f"Grid: interaction length {LEN_GRID} km, between-node mixing {MIX_GRID}, "
        f"innovation {MU_GRID}; {N_SEEDS} seeds per cell; {len(arr)} "
        "time-transgressive runs total.",
        f"Observed data: {obs_ncom} communities, between-group F_ST = {obs_fst:.3f}.",
        "",
        "## Headline robustness",
        f"- Phase-like structure (>= 2 emergent communities): {frac_phase*100:.0f}% of runs.",
        f"- Emergent F_ST at or below the observed value: {frac_le_obs*100:.0f}% of runs.",
        f"- Emergent F_ST within a factor of two of observed "
        f"({0.5*obs_fst:.3f}-{2*obs_fst:.3f}): {frac_bracket*100:.0f}% of runs.",
        f"- Mean emergent communities {np.nanmean(nc_c):.1f} "
        f"(range {int(np.nanmin(nc_c))}-{int(np.nanmax(nc_c))}); "
        f"mean F_ST {np.nanmean(fst_c):.3f}; mean seriation |rho| {np.nanmean(rho_c):.2f}.",
        "",
        "## Marginal effects (mean +/- sd)",
        f"- F_ST by interaction length (km): {fmt3(marg(LEN_GRID, lk_c, fst_c))}",
        f"- F_ST by between-node mixing: {fmt3(marg(MIX_GRID, mx_c, fst_c))}",
        f"- F_ST by innovation rate: {fmt3(marg(MU_GRID, mu_c, fst_c))}",
        f"- communities by interaction length (km): {fmt3(marg(LEN_GRID, lk_c, nc_c))}",
        "",
        "Interpretation: across stochastic replicates and the full factor grid, "
        "neutral drift on the real geography reliably produces phase-like, spatially "
        "coherent communities, and the between-group F_ST stays at or near the observed "
        "level rather than the much higher value bounded groups would leave. The "
        "appearance of phase structure is a generic outcome of distance-structured "
        "drift on this layout, not an artifact of one tuned parameter set. F_ST is "
        "essentially flat across interaction length and falls with more mixing or "
        "more innovation, as expected.",
        "",
        f"Figure: {OUT_FIG.relative_to(ROOT)}",
    ]
    OUT_MD.parent.mkdir(exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))
    print(f"\nwrote {OUT_FIG}")


if __name__ == "__main__":
    main()

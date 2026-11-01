"""Rim R1 readings writer (spec sections 4 to 7): statuses, geography statistics, pair association, sensitivities, priors.

Reads the cell CSVs written by rim_r1.py (`output/probe_latent_groups/rim_r1_<kind>_<point>.csv`) and the rim and decorated
tables, and writes `output/findings/rim_r1.md`. Three kinds of number are kept apart and named: the observed count-based
moment F; the generator band of the same statistic; posterior summaries of latent frequencies under the sampling model of
spec section 6 (direct conjugate draws; no MCMC). Descriptive throughout (rule 18): statuses below, inside, above or
unavailable; no threshold other than the band; no success level; nothing phrased as a choice among causes.
"""
from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT, ROOT / "analyses" / "per_type"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import readings as rd

from mls_emergence.rims import stats as st

OUTD = ROOT / "output" / "probe_latent_groups"
FINDINGS = ROOT / "output" / "findings" / "rim_r1.md"
ATTRIBUTES = ("beveled", "notched_exterior", "lip_flat", "paste_bell")
READ_BY_F = ("beveled", "notched_exterior", "lip_flat")
NOT_READ = ("beaded", "noded", "wide_nodes", "notched_top")
LIP_AMBIGUOUS = ("Beck", "Irby", "Starkley")
SMALL_SITE = "Woodlyn"
WINDOWS = ("T1", "8LR")
POINTS = ("primary", "secondary")
COMMON = 3
LABEL = {"beveled": "interior bevel", "notched_exterior": "exterior notching", "lip_flat": "lip flat", "paste_bell": "paste Bell (of N)"}


def observed_f(rims: pd.DataFrame, attr: str, sites=None) -> float:
    t = rims if sites is None else rims.loc[sites]
    k = t[attr].to_numpy(int); N = t["N"].to_numpy(int)
    return float(st.moment_f(np.column_stack([k, N - k]))[0])


def cell_status(cells: pd.DataFrame, observed19: float, observed18: float) -> dict:
    b19 = rd.band(cells["F19"]); b18 = rd.band(cells["F18"])
    return {"band19": b19, "band18": b18, "status19": rd.status3(observed19, b19), "status18": rd.status3(observed18, b18), "n_finite": b19["n_finite"]}


def lip_sensitivity(rims: pd.DataFrame) -> tuple[float, float, int]:
    """Range of the lip-flat F over every combination of -1, 0, +1 rims at the three sites whose file percentages do not sum to 100."""
    vals = []
    for deltas in itertools.product((-1, 0, 1), repeat=len(LIP_AMBIGUOUS)):
        t = rims.copy()
        for s, d in zip(LIP_AMBIGUOUS, deltas):
            if s in t.index:
                t.loc[s, "lip_flat"] = int(np.clip(t.loc[s, "lip_flat"] + d, 0, t.loc[s, "N"]))
        vals.append(observed_f(t, "lip_flat"))
    return float(min(vals)), float(max(vals)), len(vals)


def paste_frequencies(rims: pd.DataFrame) -> dict:
    named = rims["paste_bell"] + rims["paste_nf"]
    with np.errstate(invalid="ignore", divide="ignore"):
        of_named = (rims["paste_bell"] / named).to_numpy(float)
    return {"bell_of_N": (rims["paste_bell"] / rims["N"]).to_numpy(float), "bell_of_named": of_named,
            "unresolved_rims": int(rims["paste_other"].sum()), "sites_without_named": int((named == 0).sum())}


def margins(cells: pd.DataFrame, rims: pd.DataFrame, counts19=None, pfg19=None, classes=None) -> pd.DataFrame:
    """Per trait, point, window and site set: the band, the observed F, the status, and the margin observed minus the band's upper
    edge (positive above the band; the gap to the lower edge when below). For lip flat, whether the plus-or-minus-one
    reconstruction range [lo, hi] reaches the band edge the status was read against (crosses_lip_range)."""
    keep18 = [s for s in rims.index if s != SMALL_SITE]; lip_lo, lip_hi, _ = lip_sensitivity(rims); rows = []
    for (point, kind, attr, window), g in cells.groupby(["point", "kind", "attribute", "window"], dropna=False):
        if kind == "rim":
            items = [(attr, g["F19"], g["F18"], observed_f(rims, attr), observed_f(rims, attr, keep18))]
        else:
            tab = counts19 if kind == "decoration_current" else pfg19
            if tab is None:
                continue
            f19 = st.moment_f(tab.to_numpy(float)); f18 = st.moment_f(tab.loc[keep18].to_numpy(float))
            items = [(f"{classes[k]}, {kind}", g[f"F19_{k}"], g[f"F18_{k}"], f19[k], f18[k]) for k in range(COMMON)]
        for trait, v19, v18, o19, o18 in items:
            for n_sites, vals, obs in ((19, v19, o19), (18, v18, o18)):
                b = rd.band(vals); status = rd.status3(obs, b)
                margin = np.nan if status == "unavailable" else (obs - b["hi"] if status != "below" else obs - b["lo"])
                edge = b["hi"] if status == "above" else (b["lo"] if status == "below" else np.nan)
                crosses = bool(trait == "lip_flat" and np.isfinite(edge) and lip_lo <= edge <= lip_hi)
                rows.append({"point": point, "window": window, "trait": trait, "sites": n_sites, "lo": b["lo"], "hi": b["hi"], "observed": obs, "status": status,
                                 "margin": margin, "crosses_lip_range": crosses})
    return pd.DataFrame(rows)


def _fmt(x, nd=3):
    return "undefined" if x is None or not np.isfinite(x) else f"{x:.{nd}f}"


def _band(b, nd=3):
    return f"[{_fmt(b['lo'], nd)}, {_fmt(b['hi'], nd)}] n {b['n_finite']}" + ("" if b["eligible"] else " (ineligible)")


def _post(s: dict) -> str:
    return f"{_fmt(s['observed'], 2)}; posterior {_fmt(s['median'], 2)} [{_fmt(s['q05'], 2)}, {_fmt(s['q95'], 2)}] (MC error of the median {s['mc_error_median']:.3f})"


def _comp_posterior(fn, counts: np.ndarray, draws: int, seeds=(21, 22), alpha=0.5) -> dict:
    """Posterior summary of a statistic of one decorated class's shares under Dirichlet(counts + alpha) draws of the compositions."""
    obs = fn(counts / counts.sum(1, keepdims=True)); runs = []
    for s in seeds:
        C = st.composition_draws(counts, np.random.default_rng(s), draws, alpha)
        vals = np.array([fn(C[d]) for d in range(draws)], float); fin = vals[np.isfinite(vals)]
        runs.append({"median": float(np.median(fin)) if len(fin) else np.nan, "q05": float(np.quantile(fin, .05)) if len(fin) else np.nan,
                     "q95": float(np.quantile(fin, .95)) if len(fin) else np.nan})
    a, b = runs
    return {"observed": float(obs), **a, "mc_error_median": abs(a["median"] - b["median"])}


def write_findings(G: dict, cells: pd.DataFrame | None = None, out: Path | None = None, draws: int = 4000, recovery: bool = True) -> str:
    """Render the findings and write them to `out`, or to the committed path `output/findings/rim_r1.md` when `out` is None."""
    rims: pd.DataFrame = G["rims"]; sites = list(rims.index); keep18 = [s for s in sites if s != SMALL_SITE]
    counts19: pd.DataFrame = G["counts19"]; pfg19: pd.DataFrame = G["pfg19"]; classes = G["classes"]
    D = np.asarray(G["D_obs"], float); lat = np.asarray(G["lat"], float); lon = np.asarray(G["lon"], float)
    if cells is None:
        files = sorted(OUTD.glob("rim_r1_*.csv"))
        if not files:
            raise FileNotFoundError("no rim_r1 cell CSVs; run rim_r1.py first")
        cells = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
        dup = cells.duplicated(subset=["point", "kind", "attribute", "window", "record"], keep=False)
        if dup.any():
            raise ValueError("duplicate record keys across cell files")
    points = [p for p in POINTS if (cells.point == p).any()]
    per_cell = cells.groupby(["point", "kind", "attribute", "window"], dropna=False).size()
    n_rec = f"{per_cell.min()}" if per_cell.min() == per_cell.max() else f"{per_cell.min()} to {per_cell.max()}"
    L = ["# Rim attributes R1: readings beside decoration on the 19 rim sites", "",
         (f"Produced by `analyses/rims/rim_r1.py --readings` (git {G['prov'].get('git_rev', '')[:12]}) from the cell CSVs and the rim and "
          "decorated tables. Three kinds of number appear and are named: the observed count-based moment F (binomial sampling variance "
          "removed; `mls_emergence.signatures.variance.moment_f`); the generator band of the same statistic (10th to 90th percentile of "
          f"{n_rec} records, NumPy linear rule, endpoints inclusive, eligible at {rd.MIN_FINITE} finite values), with status below, inside, above "
          "or unavailable, of which only 'above' is read as more variation than the baseline; and posterior summaries of latent site "
          f"frequencies under Beta(k + 1/2, N - k + 1/2) per site and Dirichlet(counts + 1/2) per composition ({draws:,} direct draws, two "
          "seeds, the Monte Carlo error of each median shown), with Beta(1, 1) and Dirichlet(1) sensitivities. The yardstick is the "
          "fixed-target copying simulator with all 28 populations transmitting and the 19 rim sites observed, at the frozen step 0 "
          "time scales, two windows (T1: ends over the last L_d; 8LR: the last 8 L_R). Everything is conditional on the 19-site sample "
          f"(absent: {', '.join(G['unmatched'])}), on that yardstick and on that population model; {SMALL_SITE} (27 rims) is below "
          "Mainfort's 45-rim minimum and Mound Place (48) just above it, so every statistic is also given on the 18 sites without "
          f"{SMALL_SITE}. Nothing here is phrased as a choice among causes (spec section 1). Descriptive summaries under rule 18."), ""]
    # ---- A. statuses
    L += ["## A. Per-attribute F against the band, beside the common decorated classes on the same sites", ""]
    for point in points:
        for window in WINDOWS:
            L += [f"### {point} point, window {window}", "",
                  "| trait | observed F (19) | band (19) | status | observed F (18) | band (18) | status (18) |", "|---|---|---|---|---|---|---|"]
            for a in READ_BY_F + ("paste_bell",):
                c = cells[(cells.point == point) & (cells.kind == "rim") & (cells.attribute == a) & (cells.window == window)]
                if c.empty:
                    continue
                o19, o18 = observed_f(rims, a), observed_f(rims, a, keep18); s = cell_status(c, o19, o18)
                L.append(f"| {LABEL[a]}{' (read separately)' if a == 'paste_bell' else ''} | {_fmt(o19)} | {_band(s['band19'])} | **{s['status19']}** | {_fmt(o18)} | {_band(s['band18'])} | {s['status18']} |")
            for kind, tab, label in (("decoration_current", counts19, "current counts"), ("decoration_pfg", pfg19, "1940 PFG parts")):
                c = cells[(cells.point == point) & (cells.kind == kind) & (cells.window == window)]
                if c.empty:
                    continue
                f19 = st.moment_f(tab.to_numpy(float)); f18 = st.moment_f(tab.loc[keep18].to_numpy(float))
                for k in range(COMMON):
                    b19 = rd.band(c[f"F19_{k}"]); b18 = rd.band(c[f"F18_{k}"])
                    L.append(f"| {classes[k]}, {label} | {_fmt(f19[k])} | {_band(b19)} | **{rd.status3(f19[k], b19)}** | {_fmt(f18[k])} | {_band(b18)} | {rd.status3(f18[k], b18)} |")
            L.append("")
    lo, hi, ncomb = lip_sensitivity(rims)
    L += [(f"Lip form is a reconstructed field: the file's flat and round percentages sum to 99, 101 and 98 at {', '.join(LIP_AMBIGUOUS)}; over "
           f"the {ncomb} combinations of one rim more or fewer at those sites the lip-flat F ranges {lo:.3f} to {hi:.3f}."), ""]
    M = margins(cells, rims, counts19, pfg19, classes); M_rim = M[M.trait.isin(READ_BY_F + ("paste_bell",))]
    L += [("Smallest margin of the observed F to the band edge it was read against, over the points, windows and both site sets (observed minus "
          "the upper edge when above; four decimals, since a status at the third decimal can rest on the fourth):"), ""]
    for a in READ_BY_F + ("paste_bell",):
        g = M_rim[M_rim.trait == a].dropna(subset=["margin"])
        if g.empty:
            L.append(f"- {LABEL[a]}: no eligible cell."); continue
        r = g.loc[g.margin.abs().idxmin()]
        note = (f"; the lip reconstruction range {lo:.4f} to {hi:.4f} contains that band edge ({r.hi:.4f}), so the status at this cell is not "
                "resolved by the record" if r.crosses_lip_range else "")
        L.append(f"- {LABEL[a]}: {r.margin:+.4f} at the {r.point} point, window {r.window}, {int(r.sites)} sites (observed {r.observed:.4f}, band "
                 f"[{r.lo:.4f}, {r.hi:.4f}]); {int((g.status == 'above').sum())} of {len(g)} cells above{note}.")
    L.append("")
    # ---- B. geography
    L += ["## B. Association with geography", "",
          ("Spearman of |p_i - p_j| with straight-line distance over the 171 pairs (positive: pairs farther apart differ more), and of p "
           "with latitude (positive: higher to the north) and with longitude (positive: higher to the east). Observed value; posterior median "
           "[5%, 95%]; the Beta(1, 1) or Dirichlet(1) posterior median in the last column. Rows for the 18 sites give observed values only."), "",
          "| trait | rho pair difference vs distance | rho vs latitude | rho vs longitude | alpha 1 medians (lat, lon, dist) |", "|---|---|---|---|---|"]
    geo_rows = []
    for a in READ_BY_F + ("paste_bell",):
        k = rims[a].to_numpy(int); N = rims["N"].to_numpy(int)
        summ = {key: st.posterior_summary(lambda p, key=key: st.geography_stats(p, D, lat, lon)[key], k, N, draws, seeds=(11, 12)) for key in ("rho_pair_distance", "rho_latitude", "rho_longitude")}
        alt = {key: st.posterior_summary(lambda p, key=key: st.geography_stats(p, D, lat, lon)[key], k, N, draws, seeds=(13, 14), alpha=1.0)["median"] for key in ("rho_pair_distance", "rho_latitude", "rho_longitude")}
        geo_rows.append((LABEL[a], summ, alt))
        i18 = [sites.index(s) for s in keep18]
        g18 = st.geography_stats((k / N)[i18], D[np.ix_(i18, i18)], lat[i18], lon[i18])
        L.append(f"| {LABEL[a]} | {_post(summ['rho_pair_distance'])} | {_post(summ['rho_latitude'])} | {_post(summ['rho_longitude'])} | {alt['rho_latitude']:+.2f}, {alt['rho_longitude']:+.2f}, {alt['rho_pair_distance']:+.2f} |")
        L.append(f"| {LABEL[a]}, 18 sites | {_fmt(g18['rho_pair_distance'], 2)} | {_fmt(g18['rho_latitude'], 2)} | {_fmt(g18['rho_longitude'], 2)} | |")
    for tab, label in ((counts19, "current counts"), (pfg19, "1940 PFG parts")):
        Y = tab.to_numpy(float)
        for c in range(COMMON):
            summ = {key: _comp_posterior(lambda C, key=key, c=c: st.geography_stats(C[:, c], D, lat, lon)[key], Y, draws) for key in ("rho_pair_distance", "rho_latitude", "rho_longitude")}
            alt = {key: _comp_posterior(lambda C, key=key, c=c: st.geography_stats(C[:, c], D, lat, lon)[key], Y, draws, seeds=(23, 24), alpha=1.0)["median"] for key in ("rho_pair_distance", "rho_latitude", "rho_longitude")}
            L.append(f"| {classes[c]}, {label} | {_post(summ['rho_pair_distance'])} | {_post(summ['rho_latitude'])} | {_post(summ['rho_longitude'])} | {alt['rho_latitude']:+.2f}, {alt['rho_longitude']:+.2f}, {alt['rho_pair_distance']:+.2f} |")
    L.append("")
    # ---- C. association
    L += ["## C. Rim and decoration pair differences after linear adjustment on ranks for distance", "",
          ("Rim difference |p_i - p_j|; decoration difference the total variation 0.5 sum_c |p_ic - p_jc| over the ten classes (zero for identical "
           "compositions). The statistic is the Spearman correlation of the two after least-squares removal of ranked distance from each: "
           "association after linear adjustment on ranks, not 'at a given distance' (a nonlinear dependence on distance survives it). "
           "Positive: pairs that differ more in the attribute also differ more in decoration. Each rim attribute separately; no combined rim "
           "measure, since the marginal table carries no co-occurrence and the two collections' physical overlap is unknown."), "",
          "| trait | against current counts: plain, adjusted | against 1940 PFG parts: plain, adjusted | alpha 1 adjusted medians |", "|---|---|---|---|"]
    for a in READ_BY_F + ("paste_bell",):
        k = rims[a].to_numpy(int); N = rims["N"].to_numpy(int); cells_txt = []; alts = []
        for tab in (counts19, pfg19):
            Y = tab.to_numpy(int)
            s_plain = st.posterior_summary(lambda p, C: st.pair_association(p, C, D)["spearman"], k, N, draws, seeds=(31, 32), counts=Y)
            s_adj = st.posterior_summary(lambda p, C: st.pair_association(p, C, D)["partial_spearman"], k, N, draws, seeds=(33, 34), counts=Y)
            alt = st.posterior_summary(lambda p, C: st.pair_association(p, C, D)["partial_spearman"], k, N, draws, seeds=(35, 36), counts=Y, alpha=1.0, comp_alpha=1.0)["median"]
            cells_txt.append(f"{_post(s_plain)}; adjusted {_post(s_adj)}"); alts.append(f"{alt:+.2f}")
        L.append(f"| {LABEL[a]} | {cells_txt[0]} | {cells_txt[1]} | {', '.join(alts)} |")
    L.append("")
    # ---- D. paste
    pf = paste_frequencies(rims); n_named = len(rims) - pf["sites_without_named"]
    named_note = "" if pf["sites_without_named"] == 0 else f" over the {n_named} of {len(rims)} sites with any named paste"
    L += ["## D. Paste, separately", "",
          (f"Paste is a type-level distinction of the 1951 typology (Bell Plain against Neeley's Ferry Plain) and is not independent of how the "
           f"phases were drawn; the file's paste percentages leave {pf['unresolved_rims']} of {int(rims['N'].sum())} rims unresolved, which are not "
           "negative Bell observations. The tables above read Bell of N; here Bell is also given with the named pastes as denominator (two "
           f"denominators). Bell of N: mean site frequency {np.nanmean(pf['bell_of_N']):.3f}; Bell of named: {np.nanmean(pf['bell_of_named']):.3f}{named_note}."), ""]
    # ---- E. not read
    L += ["## E. Published attributes reported but not read", "", "| attribute | rims | sites with any (of 19) | pooled frequency |", "|---|---|---|---|"]
    for a in NOT_READ:
        L.append(f"| {a} | {int(rims[a].sum())} | {int((rims[a] > 0).sum())} | {rims[a].sum() / rims['N'].sum():.4f} |")
    L += ["", "A recorded scope choice (at least 100 rims and 12 sites with any), not a demonstrated resolution limit. Bevel degrees as shares of the beveled rims: "
          + ", ".join(f"{d.split('_')[1]} {int(rims[d].sum())}" for d in ("bevel_sharp", "bevel_moderate", "bevel_slight")) + f" of {int(rims['beveled'].sum())}.", ""]
    # ---- F. priors
    N = rims["N"].to_numpy(int); pp = st.prior_predictive_counts(N, 0.5)
    L += ["## F. The priors (rule 20)", "",
          (f"Beta(1/2, 1/2) on each site frequency: one pseudo-observation in total, weight 1/(N + 1) on the posterior mean ({1/(N.min()+1):.3f} at "
           f"{N.min()} rims, {1/(N.max()+1):.4f} at {N.max()}); direction toward 1/2. Observable prior predictive, Beta-Binomial(N, 1/2, 1/2) at each "
           "site's N: 5th, 50th and 95th percentiles of the predicted count against the observed count for the three read attributes:"), "",
          "| site | N | prior predictive count (5%, 50%, 95%) | bevel | notching | lip flat |", "|---|---|---|---|---|---|"]
    for i, s in enumerate(sites):
        L.append(f"| {s} | {N[i]} | {int(pp[i,0])}, {int(pp[i,1])}, {int(pp[i,2])} | {int(rims.loc[s, 'beveled'])} | {int(rims.loc[s, 'notched_exterior'])} | {int(rims.loc[s, 'lip_flat'])} |")
    obs_counts = rims[list(READ_BY_F)].to_numpy(int); inside = int(((obs_counts >= pp[:, [0]]) & (obs_counts <= pp[:, [2]])).sum())
    span = (pp[:, 2] - pp[:, 0]) / N
    L += ["", (f"The 5th to 95th percentile interval of the predicted count spans {100 * span.min():.0f} to {100 * span.max():.0f} percent of N across the "
               f"sites and contains {inside} of the {obs_counts.size} observed counts: the prior concentrates the predicted count nowhere, so it does not "
               "already resemble the data. The Beta(1, 1) sensitivity medians are in the last columns of sections B and C."), ""]
    if recovery:
        reps, rec_draws = 20, 1000; rec = st.recovery_program(D, lat, lon, N, replicates=reps, draws=rec_draws, seed=7)
        L += [(f"Recovery program: {reps} replicates at three true configurations on these {len(N)} sites' totals, {rec_draws:,} draws per fit, 90 percent "
               f"intervals; a coverage proportion at {reps} replicates has binomial standard deviation {np.sqrt(0.9 * 0.1 / reps):.3f} at a true rate of 0.9, "
               "so the table resolves gross failure, not the nominal rate."), "",
              "| configuration | true rho latitude | coverage of rho latitude | coverage of rho pair distance | sign recovered |", "|---|---|---|---|---|"]
        for name, r in rec.items():
            L.append(f"| {name} | {r['truth_rho_latitude']:+.2f} | {r['coverage']:.2f} | {r['coverage_pair_distance']:.2f} | {r['sign_recovered']:.2f} |")
        L.append("")
    # ---- G. limits
    L += ["## G. Limits stated in advance", "",
          ("19 sites of 28; rims from collections Mainfort names only in aggregate, never per site and never as the 1940 survey collections, so "
           "whether a site's rims and its decorated counts share material is unknown; one observer, which helps consistency but does not make "
           "recording independent across sites; totals 27 to 277, so per-site frequencies carry wide intervals and the F estimator's own variance "
           "at 19 sites is large; the yardstick is the fixed-target baseline (no regional drift) under one population model; interior beveling "
           "is spatially as well as temporally patterned by Mainfort's own account (p. 37), so no clock reading is made."), ""]
    text = "\n".join(L)
    if out is None:
        FINDINGS.parent.mkdir(parents=True, exist_ok=True); FINDINGS.write_text(text); print(f"wrote {FINDINGS}")
    else:
        Path(out).parent.mkdir(parents=True, exist_ok=True); Path(out).write_text(text); print(f"wrote {out}")
    return text

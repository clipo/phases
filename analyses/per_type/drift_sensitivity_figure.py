"""Readings of the regional-drift sensitivity run (regional-drift spec draft 3, section 5): coverage sets, separation, figure.

Reads every C3 and M per-record CSV written by per_class_table.py
(output/probe_latent_groups/per_type_table_C3_<point>.csv and
per_type_table_M_<point>_<class>_<cluster>.csv), refuses duplicate record
keys, and writes:

- output/findings/per_type_record_row.md: the record's own per-class row and
  the band and quantile conventions;
- output/findings/per_type_drift_sensitivity.md: per operating point the
  10-class by amplitude status table (inside / outside / unavailable), the
  sets C, R, H and their availability sets with the pre-stated wording, the
  separation table per planted class, cluster, m and F_W (n_out of
  n_out + n_in, undefined beside it, the neutral cell's own count outside
  its band, the planted median and its eligibility, paired changes), and the
  sentence of spec section 1 filled in for the primary point;
- figures/per_type_drift_sensitivity.png (7 in wide, 300 dpi, no title,
  sans-serif, Okabe-Ito): (a) common-class F against F_W with the neutral
  band and the record's values; (b) planted medians and the fraction outside
  the band at each m; (c) rare-class sites present with the record's values
  and the planted Walls_Engraved occupancy.

Descriptive summaries only (D28, D32; rule 18's carve-out). No threshold
other than the band, no test statistic, no success level.
--standin builds a pseudo C3 frame from the step 0 C1 run (primary point,
W = 8 L_R cell, labeled F_W = 0) and writes everything under --out-dir, so
the script can be exercised before the C3 and M data exist; nothing it
produces in that mode is a finding.
Usage: .venv/bin/python analyses/per_type/drift_sensitivity_figure.py [--standin --out-dir DIR]
"""
from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT, ROOT / "analyses", ROOT / "analyses" / "probe_latent_groups", ROOT / "analyses" / "per_type"):
    sys.path.insert(0, str(_p))

import readings as rd
from figstyle import OKABE_ITO_COLOR

from mls_emergence.dataio.basin import load_curated

OUTD = ROOT / "output" / "probe_latent_groups"
FINDD = ROOT / "output" / "findings"
FIGD = ROOT / "figures"
RECORD_ROW = FINDD / "per_type_record_row.md"
FINDINGS = FINDD / "per_type_drift_sensitivity.md"
FIGURE = FIGD / "per_type_drift_sensitivity.png"
COMMON = 3
GRID = (0.0, 0.01, 0.03, 0.1, 0.3)
SEED_BASE_TEXT = "2026101950 + 100000 point + 1000 a + r"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans", "Helvetica"], "font.size": 8,
                     "axes.spines.top": False, "axes.spines.right": False})


def record_row() -> tuple[dict, list[str]]:
    pct = importlib.import_module("per_class_table")
    from a0_predictive_summaries import spatial_clusters
    m47 = importlib.import_module("47_revision_analysis")
    counts, coords = load_curated(); labels = spatial_clusters(coords.to_numpy(float))
    return pct.per_class_rows(counts.to_numpy(int), labels, m47.diversity), list(counts.columns)


def record_row_text(real: dict, cls: list[str]) -> str:
    L = ["# The record's per-class row", "",
         ("Produced by `analyses/per_type/drift_sensitivity_figure.py` from the curated 28-assemblage matrix and the two-cluster "
          "spatial partition (`a0_predictive_summaries.spatial_clusters`). F is the site-weighted method-of-moments F with binomial "
          "sampling variance removed; sites is the number of the 28 sites at which the class is present. Bands in the sensitivity "
          f"findings are the 10th to 90th percentile of a cell's finite values under NumPy's linear quantile rule, endpoints inclusive, "
          f"eligible when at least {rd.MIN_FINITE} of the records give a finite value."), "",
         "| class | F | sites present | status | Spearman share vs size |", "|---|---|---|---|---|"]
    for c, name in enumerate(cls):
        L.append(f"| {name} | {real[f'F_{c}']:.4f} | {real[f'occ_{c}']} | {real[f'status_{c}']} | {real[f'rho_{c}']:+.2f} |")
    L += ["", (f"Site F_ST (sherd-weighted pool) {real['site_fst']:.4f}; site F_ST (site-mean pool) {real['site_fst_sitepool']:.4f}; "
               f"two-cluster F_ST {real['fst_clusters']:.4f}; richness {real['rich']:.2f}."), ""]
    return "\n".join(L)


def load_frames(outd: Path, standin: bool) -> tuple[pd.DataFrame, pd.DataFrame]:
    if standin:
        c1 = pd.read_csv(outd / "per_type_table_C1.csv")
        d = c1[(c1.point == "primary") & (c1.setting.astype(str) == "8")].copy()
        if d.empty:
            raise FileNotFoundError("the stand-in needs the C1 run's primary W = 8 L_R cell")
        d["generator"] = "C3"; d["F_W"] = 0.0; d["a"] = 0; d["kappa"] = 0.0; d["total_f_S"] = 0.0
        d["record"] = d["seed"] if "seed" in d else d["record"]
        for k in range(10):
            d[f"status_{k}"] = np.where(d[f"occ_{k}"] == 0, "absent", np.where(np.isfinite(d[f"F_{k}"]), "finite", "undefined"))
            d[f"lost_{k}"] = False; d[f"pfT_{k}"] = np.nan; d[f"pfB_{k}"] = np.nan
        d["reason"] = ""
        return rd.load_runs([d]), pd.DataFrame()
    c3 = sorted(outd.glob("per_type_table_C3_*.csv")); mm = sorted(outd.glob("per_type_table_M_*.csv"))
    if not c3:
        raise FileNotFoundError("no C3 run found; run per_class_table.py --generator C3 first")
    neutral = rd.load_runs([pd.read_csv(p) for p in c3])
    markers = rd.load_runs([pd.read_csv(p) for p in mm]) if mm else pd.DataFrame()
    return neutral, markers


def statuses_for_point(neutral: pd.DataFrame, real: dict, cls: list[str]) -> tuple[dict, dict]:
    """Per amplitude: per-class band and status; returns (bands[a][c], statuses[a] = {"common": [...], "rare": [...]})."""
    bands, st = {}, {}
    for a in GRID:
        cell = neutral[np.isclose(neutral.F_W, a)]
        if cell.empty:
            continue
        bands[a] = {}; sc = []; sr = []
        for c in range(len(cls)):
            col = f"F_{c}" if c < COMMON else f"occ_{c}"
            b = rd.band(cell[col]); bands[a][c] = b
            (sc if c < COMMON else sr).append(rd.status(real[col], b))
        st[a] = {"common": sc, "rare": sr}
    return bands, st


def status_table(bands: dict, st: dict, real: dict, cls: list[str], neutral: pd.DataFrame) -> list[str]:
    grid = [a for a in GRID if a in bands]
    L = ["| class | record | " + " | ".join(f"F_W {a:g}" for a in grid) + " |", "|---|---|" + "---|" * len(grid)]
    for c, name in enumerate(cls):
        col = f"F_{c}" if c < COMMON else f"occ_{c}"
        cells = []
        for a in grid:
            b = bands[a][c]; s = (st[a]["common"] + st[a]["rare"])[c]
            cell = neutral[np.isclose(neutral.F_W, a)]; acc = rd.accounting(cell, c)
            fmt = (lambda x: f"{x:.3f}") if c < COMMON else (lambda x: f"{x:.0f}")
            cells.append(f"{fmt(b['median'])} [{fmt(b['lo'])}, {fmt(b['hi'])}] n {b['n_finite']}: **{s}**"
                         + (f"; lost {acc['lost']}" if acc["lost"] else "") + (f"; absent {acc['absent']}" if acc["absent"] else "")
                         + (f"; fixed {acc['fixed']}" if acc["fixed"] else "") + (f"; undefined pool {acc['undefined_pool']}" if acc["undefined_pool"] else ""))
        L.append(f"| {name} | {real[col]:.3f} | " if c < COMMON else f"| {name} | {real[col]} | ")
        L[-1] += " | ".join(cells) + " |"
    return L


def separation_tables(neutral: pd.DataFrame, markers: pd.DataFrame, cls: list[str], point: str) -> tuple[list[str], dict]:
    """Per (class, cluster): rows per (F_W, m) with the separation counts; returns the lines and the first-m readings."""
    L, firsts = [], {}
    if markers.empty:
        return ["No marker run found for this point."], firsts
    mk = markers[markers.point == point]
    for (mclass, cluster), grp in mk.groupby(["marker_class", "cluster"], sort=True):
        c = cls.index(mclass); col = f"F_{c}" if c < COMMON else f"occ_{c}"; kind = "F" if c < COMMON else "sites present"
        size = int(grp.cluster_size.iloc[0])
        L += ["", f"### {mclass}, {cluster} cluster ({size} sites), statistic: {kind}", "",
              "| F_W | m | planted outside the neutral band | neutral outside its own band | planted median (eligible) | neutral band | paired: up / down / tied / unavailable |",
              "|---|---|---|---|---|---|---|"]
        for a in GRID:
            ncell = neutral[(neutral.point == point) & np.isclose(neutral.F_W, a)].sort_values("record")
            if ncell.empty:
                continue
            b = rd.band(ncell[col]); nsep = rd.separation(ncell[col], b); first = None; any_planted = False
            for m in sorted(grp.m.unique()):
                pcell = grp[np.isclose(grp.F_W, a) & np.isclose(grp.m, m)].sort_values("record")
                if pcell.empty:
                    continue
                any_planted = True
                if not np.array_equal(pcell.record.to_numpy(), ncell.record.to_numpy()):
                    raise ValueError(f"planted and neutral records do not pair at F_W {a}, m {m}")
                sep = rd.separation(pcell[col], b); pc = rd.paired_changes(pcell[col], ncell[col])
                med = (f"{sep['median']:.3f}" if c < COMMON else f"{sep['median']:.0f}") if np.isfinite(sep["median"]) else "undefined"
                med += " (eligible)" if sep["median_eligible"] else " (ineligible)"
                if first is None and sep["available"] and sep["median_eligible"] and rd.status(sep["median"], b) == "outside":
                    first = m
                lost_n = int(ncell[f"lost_{c}"].astype(bool).sum()) if f"lost_{c}" in ncell else 0
                note = f" (class lost from the target in {lost_n} of {len(ncell)} neutral records; paired direction uninformative)" if lost_n >= len(ncell) // 2 else ""
                L.append(f"| {a:g} | {m:g} | {sep['text']} | {nsep['text']} | {med} | "
                         + (f"[{b['lo']:.3f}, {b['hi']:.3f}]" if c < COMMON else f"[{b['lo']:.0f}, {b['hi']:.0f}]")
                         + (" n " + str(b["n_finite"]) + ("" if b["eligible"] else ", ineligible"))
                         + f" | {pc['increased']} / {pc['decreased']} / {pc['tied']} / {pc['unavailable']}{note} |")
            if not any_planted:
                firsts[(mclass, cluster, a)] = "not run"
            elif not b["eligible"]:
                firsts[(mclass, cluster, a)] = "unavailable"
            else:
                firsts[(mclass, cluster, a)] = f"{first:g}" if first is not None else "none on this grid"
        L.append("")
        L.append("First tested m whose eligible planted median lies outside the neutral band: "
                 + "; ".join(f"F_W {a:g}: {firsts[(mclass, cluster, a)]}" for a in GRID if (mclass, cluster, a) in firsts) + ".")
    return L, firsts


def sentence(point: str, sets: dict, firsts: dict, markers: pd.DataFrame, neutral: pd.DataFrame, cls: list[str]) -> str:
    n_rec = int(neutral[neutral.point == point].groupby("F_W").size().max())
    W = int(neutral[neutral.point == point].W.iloc[0])
    nf = []
    for a in GRID:
        cell = neutral[(neutral.point == point) & np.isclose(neutral.F_W, a)]
        if len(cell):
            nf += [rd.band(cell[f"F_{c}"])["n_finite"] for c in range(COMMON)]
    finite_clause = (f"every common-class band over {nf[0]} finite values" if nf and min(nf) == max(nf)
                     else f"common-class bands over {min(nf)} to {max(nf)} finite values of {n_rec}")
    s = (f"Under neutral copying with regional drift of the innovation target (Dirichlet recursion), at the {point} operating point, "
         f"window W = 8 L_R = {W} steps, {n_rec} records per cell ({finite_clause}), the record's common-class F values and rare-class "
         f"occupancies are {rd.coverage_wording(sets)}.")
    if not markers.empty:
        parts = []
        for (mclass, cluster), grp in markers[markers.point == point].groupby(["marker_class", "cluster"], sort=True):
            c = cls.index(mclass); col = f"F_{c}" if c < COMMON else f"occ_{c}"
            bits = []
            for a in GRID:
                ncell = neutral[(neutral.point == point) & np.isclose(neutral.F_W, a)]
                if ncell.empty:
                    continue
                b = rd.band(ncell[col]); ms = []
                for m in sorted(grp.m.unique()):
                    pcell = grp[np.isclose(grp.F_W, a) & np.isclose(grp.m, m)]
                    if pcell.empty:
                        continue
                    ms.append(f"m {m:g}: {rd.separation(pcell[col], b)['text']}")
                if not ms:
                    continue
                bits.append(f"F_W {a:g} ({'; '.join(ms)}; first tested m with the eligible median outside: {firsts.get((mclass, cluster, a), 'unavailable')})")
            parts.append(f"For {mclass} multiplied on the odds scale in the {cluster} cluster ({int(grp.cluster_size.iloc[0])} sites), planted "
                         f"records with a defined statistic falling outside the unmarked band: " + "; ".join(bits) + ".")
        s += " " + " ".join(parts)
    s += (" These comparisons do not identify the cause of the observed pattern and do not exclude marker patterns outside this "
          "planting design.")
    return s


def figure(neutral: pd.DataFrame, markers: pd.DataFrame, real: dict, cls: list[str], point: str, path: Path) -> None:
    n = neutral[neutral.point == point]; grid = [a for a in GRID if np.isclose(n.F_W, a).any()]
    x = np.arange(len(grid)); colors = OKABE_ITO_COLOR
    fig, axes = plt.subplots(3, 3, figsize=(7, 7.5))
    for c in range(COMMON):
        ax = axes[0, c]; bands = [rd.band(n[np.isclose(n.F_W, a)][f"F_{c}"]) for a in grid]
        ax.fill_between(x, [b["lo"] for b in bands], [b["hi"] for b in bands], color=colors[1], alpha=.3, lw=0, label="neutral 10th to 90th")
        ax.plot(x, [b["median"] for b in bands], color=colors[1], lw=1.2, marker="o", ms=2.5, label="neutral median")
        ax.vlines(x, [b["lo"] for b in bands], [b["hi"] for b in bands], color=colors[1], lw=0.8)
        ax.axhline(real[f"F_{c}"], color=colors[7], lw=1, ls="--", label="record")
        ax.set_xticks(x); ax.set_xticklabels([f"{a:g}" for a in grid]); ax.set_xlabel("F_W"); ax.set_ylabel(f"F, {cls[c]}", fontsize=7)
        ax.text(0.02, 0.95, "(a)" if c == 0 else "", transform=ax.transAxes, va="top")
        ax2 = axes[1, c]
        if not markers.empty:
            mk = markers[(markers.point == point) & (markers.marker_class == cls[c])]
            for j, (cluster, grp) in enumerate(mk.groupby("cluster", sort=True)):
                for i, m in enumerate(sorted(grp.m.unique())):
                    meds, frac = [], []
                    for a, b in zip(grid, bands):
                        pcell = grp[np.isclose(grp.F_W, a) & np.isclose(grp.m, m)]
                        sep = rd.separation(pcell[f"F_{c}"], b) if len(pcell) else {"median": np.nan, "available": False}
                        meds.append(sep["median"]); frac.append(sep["n_out"] / (sep["n_out"] + sep["n_in"]) if sep["available"] and (sep["n_out"] + sep["n_in"]) else np.nan)
                    ax.plot(x, meds, color=colors[2 + j], lw=0.8, ls=["-", "--", ":"][i % 3], label=f"{cluster}, m {m:g}" if c == 0 else None)
                    ax2.plot(x, frac, color=colors[2 + j], lw=0.8, ls=["-", "--", ":"][i % 3], marker="o", ms=2)
        if c == 0:
            ax.legend(fontsize=5.5, frameon=False, loc="upper left", bbox_to_anchor=(0, 0.92), ncol=2)
        ax2.set_ylim(-0.02, 1.02); ax2.set_xticks(x); ax2.set_xticklabels([f"{a:g}" for a in grid]); ax2.set_xlabel("F_W")
        ax2.set_ylabel(f"planted outside band, {cls[c]}", fontsize=7); ax2.text(0.02, 0.95, "(b)" if c == 0 else "", transform=ax2.transAxes, va="top")
    gs = axes[2, 0].get_gridspec()
    for ax in axes[2]:
        ax.remove()
    ax3 = fig.add_subplot(gs[2, :])
    for k, c in enumerate(range(COMMON, len(cls))):
        bands = [rd.band(n[np.isclose(n.F_W, a)][f"occ_{c}"]) for a in grid]
        ax3.plot(x + 0.04 * (k - 3), [b["median"] for b in bands], color=colors[k % 8], lw=0.9, marker="o", ms=2, label=cls[c])
        ax3.vlines(x + 0.04 * (k - 3), [b["lo"] for b in bands], [b["hi"] for b in bands], color=colors[k % 8], lw=0.6)
        ax3.scatter(x + 0.04 * (k - 3), [real[f"occ_{c}"]] * len(grid), marker="_", color=colors[k % 8], s=40, lw=1.2)
    if not markers.empty and "Walls_Engraved" in cls:
        c = cls.index("Walls_Engraved"); mk = markers[(markers.point == point) & (markers.marker_class == "Walls_Engraved")]
        for j, (cluster, grp) in enumerate(mk.groupby("cluster", sort=True)):
            pm = grp[np.isclose(grp.m, grp.m.max())]
            meds = [rd.band(pm[np.isclose(pm.F_W, a)][f"occ_{c}"])["median"] for a in grid]
            ax3.plot(x, meds, color=colors[2 + j], lw=0.8, ls="--", label=f"Walls_Engraved planted, {cluster}, m {grp.m.max():g}")
    ax3.set_xticks(x); ax3.set_xticklabels([f"{a:g}" for a in grid]); ax3.set_xlabel("F_W"); ax3.set_ylabel("sites present (of 28)")
    ax3.text(0.01, 0.95, "(c)", transform=ax3.transAxes, va="top"); ax3.legend(fontsize=6, frameon=False, ncol=3, loc="upper right")
    fig.tight_layout(); path.parent.mkdir(parents=True, exist_ok=True)
    if path == FIGURE:
        fig.savefig(FIGURE, dpi=300)
    else:
        fig.savefig(path, dpi=300)
    plt.close(fig)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--standin", action="store_true", help="exercise the script on the step 0 C1 cell; writes under --out-dir only")
    ap.add_argument("--out-dir", default=None, help="directory for every output (required with --standin)")
    a = ap.parse_args(argv)
    if a.standin and not a.out_dir:
        raise SystemExit("--standin needs --out-dir so no finding is overwritten")
    out = Path(a.out_dir) if a.out_dir else None
    if out:
        out.mkdir(parents=True, exist_ok=True)
    real, cls = record_row()
    rec_text = record_row_text(real, cls)
    if out:
        (out / RECORD_ROW.name).write_text(rec_text)
    else:
        RECORD_ROW.write_text(rec_text)
    neutral, markers = load_frames(OUTD, a.standin)
    L = ["# Regional-drift sensitivity: coverage and separation readings", "",
         ("Produced by `analyses/per_type/drift_sensitivity_figure.py` from the C3 and M per-record CSVs of `per_class_table.py` "
          "(regional-drift spec draft 3, section 5). Statuses: inside, outside or unavailable (band ineligible: fewer than "
          f"{rd.MIN_FINITE} finite values). Separation: planted records with a defined statistic outside the neutral band of the same "
          "class and amplitude, with the undefined count beside them; the neutral cell's own count outside its band is a description "
          "of the band, not a calibration. Descriptive summaries only; no threshold other than the band, no test statistic, no success "
          "level; nothing here identifies the cause of the record's pattern.")
         + (" STAND-IN MODE: the step 0 C1 cell relabeled F_W = 0; not a finding." if a.standin else ""), "",
         ("Planting one class rescales the other classes' shares in the marked rows and redraws the sherds, so the realized occupancies "
          "of unplanted rare classes in a marker cell differ from the paired neutral record's; only the planted class is read. "
          "The paired direction (up / down / tied) is marked uninformative where the class was lost from the target in at least half of "
          "the neutral records."), ""]
    first_point_sentence = None; primary_H = []
    for point in [p for p in ("primary", "secondary") if (neutral.point == p).any()]:
        n = neutral[neutral.point == point]
        bands, st = statuses_for_point(n, real, cls); sets = rd.coverage_sets({a_: st[a_] for a_ in st}, [a_ for a_ in GRID if a_ in st])
        L += [f"## {point} operating point", "", (f"W = 8 L_R = {int(n.W.iloc[0])}, S = {int(n.S.iloc[0])}, L_d = {int(n.L_d.iloc[0])}, "
              f"L_R = {int(n.L_R.iloc[0])}, burn-in = {int(n.burnin.iloc[0]) if 'burnin' in n else 'see run'}, records per cell "
              f"{int(n.groupby('F_W').size().max())}."), ""]
        if point == "secondary" and not markers.empty and (markers.point == "secondary").any():
            amps = sorted(markers[markers.point == "secondary"].F_W.unique())
            L += [(f"Marker cells at this point were run at F_W {', '.join(f'{v:g}' for v in amps)} only, the fallback of spec 3.4: "
                   f"{'H was empty on the primary grid' if not primary_H else 'the first tested F_W in H on the primary grid'}."), ""]
        L += ["### Coverage status, record against the neutral bands", ""]
        L += status_table(bands, st, real, cls, n)
        if point == "primary":
            primary_H = sets["H"]
        L += ["", "Sets over the tested grid: " + "; ".join(f"{k} = {rd._set(v)}" for k, v in sets.items() if k != "unavailable_cells") + ".",
              "", f"Reading: {rd.coverage_wording(sets)}.", "", "### Separation of planted markers"]
        sep_lines, firsts = separation_tables(neutral, markers, cls, point); L += sep_lines
        if first_point_sentence is None:
            first_point_sentence = sentence(point, sets, firsts, markers, neutral, cls)
            figure(neutral, markers, real, cls, point, (out / FIGURE.name) if out else FIGURE)
        L.append("")
    n0 = neutral[neutral.point == "primary"] if (neutral.point == "primary").any() else neutral
    L += ["## Figure", "", (f"`figures/per_type_drift_sensitivity.png` ({n0.point.iloc[0]} operating point; W = 8 L_R = {int(n0.W.iloc[0])}, "
          f"S = {int(n0.S.iloc[0])}, {int(n0.groupby('F_W').size().max())} records per cell, record-key seed base {SEED_BASE_TEXT}). "
          "(a) Per-class moment F of the three common classes against F_W: neutral median and 10th to 90th percentile band (NumPy "
          "linear quantile rule, endpoints inclusive), the record's value as a dashed line, and the planted medians for each cluster "
          "(colors) and multiplier (solid 1.5, dashed 2, dotted 3). (b) Fraction of planted records with a defined statistic outside "
          "the neutral band of the same class and amplitude. (c) Rare-class sites present: neutral median and band per class, the "
          "record's value as a short horizontal bar, and the planted Walls_Engraved occupancy at m = 3 in each cluster."), ""]
    L += ["## The sentence of spec section 1, filled in", "", "> " + first_point_sentence, ""]
    text = "\n".join(L)
    if out:
        (out / FINDINGS.name).write_text(text); print(f"wrote {out / FINDINGS.name}, {out / RECORD_ROW.name}, {out / FIGURE.name}")
    else:
        FINDINGS.write_text(text); print(f"wrote {FINDINGS}, {RECORD_ROW}, {FIGURE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

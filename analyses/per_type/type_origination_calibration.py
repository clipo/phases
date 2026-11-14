"""Calibration of generator C4 at theta = 0 against generator C3 on identical ranked readouts (exploratory).

Type-origination spec 2026-11-14, draft 1, section 5, "Calibration of the
ranked readouts"; Codex build audit 2026-11-14, finding F2. The run script
set C4 by rank beside C3 by named class and printed medians only. This script
makes the comparison the spec asks for: the same function
(`ranked_readouts`) on both generators, medians and 10th to 90th percentile
bands side by side.

C3 arm: the drift-extension records of the baseline design (span limit L_d,
dates independent of location) are regenerated from their registered seeds
(`drift_extension.record_counts`, root entropy 2026111400) at the five
amplitudes and two windows of the origination grid, and read by rank. They
are the records behind `output/findings/drift_extension_cells.csv`.
C4 arm: the theta = 0 cells of `output/findings/type_origination_cells.csv`.

The two arms differ by construction in three ways, and the comparison shows
what each does to the readouts:
  1. C3 redraws the regional target from a Dirichlet at every step, so a
     class's target share can become very small but is lost only at
     floating-point zero. C4 resamples a finite pool, so a class at count
     zero is lost exactly.
  2. C4 starts from the pooled shares rounded to the pool size by largest
     remainder, so a class allocated zero draws is absent from the start.
  3. The seed streams differ, so the two arms are independent sets of 50
     records. Two medians of 50 independent records differ by sampling alone.

Descriptive comparison of simulation bands (rule 18). One seed set per arm.

Usage: .venv-np254/bin/python analyses/per_type/type_origination_calibration.py
Outputs: output/findings/type_origination_calibration.md, output/findings/type_origination_calibration.csv
"""
from __future__ import annotations

import multiprocessing as mp
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT, ROOT / "analyses", ROOT / "analyses" / "probe_latent_groups", ROOT / "analyses" / "per_type"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import drift_extension as de
import readings as rd
import type_origination as to

from mls_emergence.signatures import ranked_readouts as rr

FINDD = ROOT / "output" / "findings"
MANIFEST_OUTPUTS = ["output/findings/type_origination_calibration.md", "output/findings/type_origination_calibration.csv"]
SPAN, COUPLING = "Ld", "c0"
G: dict = {}


def c3_task(task) -> list[dict]:
    """Ranked readouts of the 50 C3 records of one (point, window, amplitude) at the baseline design."""
    point, window, F_W = task; ctx = G["ctx"]; a = de.amplitude_index(F_W); design = (SPAN, window, COUPLING); rows = []
    for r in range(de.RECORDS):
        y = de.record_counts(ctx, point, window, a, r, [design])[design]
        rows.append({"point": point, "window": window, "F_W": F_W, "record": r, **rr.ranked_readouts(y, ctx["lat"], ctx["D"])})
    return rows


def c3_bands(per_record: pd.DataFrame, real: dict) -> pd.DataFrame:
    """C3 bands by rank, with the record's status against each band (the reading rule of both specs)."""
    rows = []
    for (point, window, F_W), cell in per_record.groupby(["point", "window", "F_W"], sort=False):
        if len(cell) != de.RECORDS:
            raise ValueError(f"C3 cell {(point, window, F_W)} has {len(cell)} records, expected {de.RECORDS}")
        for s in to.STATS:
            b = rd.band(to.stat_values(cell, s))
            rows.append({"point": point, "window": window, "F_W": F_W, "stat": s, "arm": "C3", **b,
                         "status": rd.status3(to.record_value(real, s), b)})
    return pd.DataFrame(rows)


def c4_bands(cells: pd.DataFrame) -> pd.DataFrame:
    q = cells[(cells["span"] == SPAN) & (cells["coupling"] == COUPLING) & (cells["theta"] == 0) & cells["stat"].isin(to.STATS)]
    if len(q) != 2 * len(to.WINDOWS) * len(to.AMPLITUDES) * len(to.STATS):
        raise ValueError(f"expected every theta = 0 baseline cell of both points in the C4 cell table; found {len(q)} rows")
    return q[["point", "window", "F_W", "stat", "lo", "median", "hi", "n_finite", "status"]].assign(arm="C4")


def _cell(row, nd: int) -> str:
    return f"{row['median']:.{nd}f} [{row['lo']:.{nd}f}, {row['hi']:.{nd}f}]"


def compare(both: pd.DataFrame) -> pd.DataFrame:
    """One row per (point, window, F_W, stat): both arms' bands and whether each median lies inside the other arm's band."""
    w = both.pivot(index=["point", "window", "F_W", "stat"], columns="arm", values=["lo", "median", "hi", "n_finite", "status"])
    out = pd.DataFrame({f"{v}_{arm}": w[(v, arm)] for arm in ("C3", "C4") for v in ("lo", "median", "hi", "n_finite", "status")})
    for v in ("lo", "median", "hi"):
        for arm in ("C3", "C4"):
            out[f"{v}_{arm}"] = out[f"{v}_{arm}"].astype(float)
    out["same_status"] = out["status_C3"] == out["status_C4"]
    out["c4_median_in_c3_band"] = (out["median_C4"] >= out["lo_C3"]) & (out["median_C4"] <= out["hi_C3"])
    out["c3_median_in_c4_band"] = (out["median_C3"] >= out["lo_C4"]) & (out["median_C3"] <= out["hi_C4"])
    return out.reset_index()


def main() -> int:
    de.require_numpy()
    ctx = de.context(); G["ctx"] = ctx
    tasks = [(p, w, a) for p in de.POINT_INDEX for w in to.WINDOWS for a in to.AMPLITUDES]
    with mp.get_context("fork").Pool(min(mp.cpu_count(), len(tasks))) as pool:
        per_record = pd.DataFrame([r for rows in pool.map(c3_task, tasks, chunksize=1) for r in rows])
    real = rr.ranked_readouts(ctx["y"], ctx["lat"], ctx["D"])
    cmp = compare(pd.concat([c3_bands(per_record, real), c4_bands(pd.read_csv(FINDD / "type_origination_cells.csv"))], ignore_index=True))
    cmp.to_csv(FINDD / "type_origination_calibration.csv", index=False)
    nd = {s: (0 if s.startswith("occ_") or s == "K_obs" else 3) for s in to.STATS}
    label = {"F_r1": "F rank 1", "F_r2": "F rank 2", "F_r3": "F rank 3", "K_obs": "richness", "beyond10": "share beyond rank 10",
             "D_dist": "D_dist", "A_lat_r1": "|lat| rank 1", "A_lat_r2": "|lat| rank 2", "A_lat_r3": "|lat| rank 3",
             **{f"occ_r{r}": f"sites, rank {r}" for r in range(4, 11)}}
    L = ["# Type origination: calibration of C4 at theta = 0 against C3, by rank (exploratory)", "",
         ("Produced by `.venv-np254/bin/python analyses/per_type/type_origination_calibration.py`; "
          f"NumPy {np.__version__}; data hash {ctx['data_hash']}. Both arms are read with the same function "
          "(`ranked_readouts`) at the baseline design (span limit L_d, dates independent of location). C3: the drift-extension "
          "records regenerated from their registered seeds. C4: the theta = 0 cells of `type_origination_cells.csv`. Each cell is "
          "median [10th, 90th percentile] over 50 records. The arms are independent sets of records, so medians differ by "
          "sampling as well as by construction. Descriptive bands from one seed set per arm, not posterior intervals."), "",
         ("By construction C3 loses a class from the regional target only at floating-point zero, and C4 loses it when its pool "
          "count reaches zero. C4 also starts from pooled shares rounded to the pool size, so a class allocated zero draws is absent "
          "from the start (listed per cell below)."), ""]
    for point in de.POINT_INDEX:
        L += [f"## {point.capitalize()} point", ""]
        for window in to.WINDOWS:
            L += [f"### Window {window}", "",
                  "| F_W | pool size | starting classes present | classes rounded to zero draws |", "|---|---|---|---|"]
            for a in to.AMPLITUDES:
                N, k0, gone = to.initial_pool(ctx, point, window, a)
                L.append(f"| {a:g} | {N} | {k0} | {', '.join(gone) if gone else 'none'} |")
            L += ["", "| F_W | statistic | C3 by rank | C4 theta 0 by rank | C4 median in C3 band | C3 median in C4 band |", "|---|---|---|---|---|---|"]
            q = cmp[(cmp["point"] == point) & (cmp["window"] == window)].set_index(["F_W", "stat"])
            for a in to.AMPLITUDES:
                for s in to.STATS:
                    r = q.loc[(a, s)]
                    c3 = {"median": r["median_C3"], "lo": r["lo_C3"], "hi": r["hi_C3"]}; c4 = {"median": r["median_C4"], "lo": r["lo_C4"], "hi": r["hi_C4"]}
                    L.append(f"| {a:g} | {label[s]} | {_cell(c3, nd[s])} | {_cell(c4, nd[s])} | "
                             f"{'yes' if r['c4_median_in_c3_band'] else 'no'} | {'yes' if r['c3_median_in_c4_band'] else 'no'} |")
            L += [""]
    L += ["## Overlap by statistic", "",
          ("Cells are the 20 combinations of point, window and amplitude. Each check asks whether one arm's median lies inside "
           "the other arm's 10th to 90th percentile band. The check is coarse: two bands shifted by less than half their width "
           "pass it, so passing shows overlap and does not show that the generators agree to any stated tolerance (Codex repair "
           "check 2026-11-14, RC-F2). No tolerance for a material difference has been set. The last two columns give the size "
           "of what was observed: the largest absolute difference between the two medians over the 20 cells, and the "
           "largest ratio of that difference to the C3 band's width, over the cells where the width is positive."), "",
          "| statistic | C4 median inside C3 band | C3 median inside C4 band | cells where neither holds | largest median difference | as a fraction of the C3 band width |",
          "|---|---|---|---|---|---|"]
    for s in to.STATS:
        m = cmp[cmp["stat"] == s]
        neither = m[~m["c4_median_in_c3_band"] & ~m["c3_median_in_c4_band"]]
        where = "; ".join(f"{r.point} {r.window} F_W {r.F_W:g}" for r in neither.itertuples()) or "none"
        diff = (m["median_C4"] - m["median_C3"]).abs(); width = m["hi_C3"] - m["lo_C3"]; pos = width > 0
        if (diff[~pos] > 0).any():
            raise ValueError(f"{s}: the medians differ in a cell whose C3 band has zero width")
        frac = f"{(diff[pos] / width[pos]).max():.2f}" if pos.any() else "not defined"
        if not pos.all():
            frac += f" ({int((~pos).sum())} cells have a C3 band of zero width and equal medians)"
        L.append(f"| {label[s]} | {int(m['c4_median_in_c3_band'].sum())} of {len(m)} | {int(m['c3_median_in_c4_band'].sum())} of {len(m)} | {where} | "
                 f"{diff.max():.{nd[s]}f} | {frac} |")
    L += ["", "## Status of the record under each generator", "",
          ("The author's definition of a material difference (2026-11-14): one that changes the record's status (inside, above "
           "or below the band) for a statistic in a cell, since status is the reading rule of both designs. The table counts, "
           "over the 20 cells, how often the record has the same status against the C3 band and against the C4 band, and lists "
           "the cells where it does not. Each band is the 10th to 90th percentile of 50 records from one seed set, and the two "
           "arms are independent, so a record value near a band edge can change status by sampling alone. The count is a "
           "description, not a test."), "",
          "| statistic | same status | cells where the status differs (C3 status, C4 status) |", "|---|---|---|"]
    for s in to.STATS:
        m = cmp[cmp["stat"] == s]; d = m[~m["same_status"]]
        where = "; ".join(f"{r.point} {r.window} F_W {r.F_W:g} ({r.status_C3}, {r.status_C4})" for r in d.itertuples()) or "none"
        L.append(f"| {label[s]} | {int(m['same_status'].sum())} of {len(m)} | {where} |")
    L.append(f"| all statistics | {int(cmp['same_status'].sum())} of {len(cmp)} | |")
    common = cmp[cmp["stat"].isin(to.F3)]; rare = cmp[cmp["stat"].isin(to.RARE)]
    L += ["", (f"Common-class F (ranks 1 to 3): same status in {int(common['same_status'].sum())} of {len(common)}. "
               f"Rare-class occupancy (ranks 4 to 10): same status in {int(rare['same_status'].sum())} of {len(rare)}.")]
    out = FINDD / "type_origination_calibration.md"; out.write_text("\n".join(L) + "\n"); print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

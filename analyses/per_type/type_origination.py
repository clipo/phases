"""Type origination during the occupation: generator C4 against the record (exploratory).

Design: docs/superpowers/specs/2026-11-14-type-origination-generator-design.md,
draft 1. Exploratory by the author's decision of 2026-11-14: the spec was not
audited, the build and the first run were audited by Codex after the run
(docs/audits/codex_audit_2026-11-14_build_ledger.md), one seed set.

Generator C4 (`mls_emergence.transmission.origination`): the site dynamics of
C3 with the regional innovation target given by a finite pool of N_R draws,
resampled after burn-in as a Wright-Fisher population in which each draw is a
new type with probability theta / (2 N_R). N_R = round(kappa) + 1 for the
kappa that C3 uses at the same amplitude F_W and window W, so drift per step
matches C3. theta = 0 is the calibration arm (C3 with exact loss).

Grid per operating point: 2 span limits (L_d, 1) x 2 windows (8, 32 L_R) x 3
date arrangements (c0, lat100N, lat100S) x 5 amplitudes x 5 values of theta =
300 cells of 50 records. The statement form of spec section 7 is read on all
seven rare ranks; the five-of-seven threshold is a separate post hoc summary
(Codex build audit 2026-11-14, F4). The matched calibration against C3 is in
`type_origination_calibration.py` (F2). The six designs of a (window, amplitude, theta,
record) read one simulated history. The grid was narrowed after seeing the
drift-extension run.

Readouts are by rank of pooled count (`ranked_readouts`): F and |latitude
rank correlation| of ranks 1 to 3, sites present for ranks 4 to 10, richness,
share of sherds beyond rank 10, D_dist. Each cell's band is the 10th to 90th
percentile over its records (eligible at 45 finite values); the record's
status is inside, above, below or unavailable. Descriptive neutral
comparisons of marginal coverage (rule 18), not posterior intervals.

Seeds: SeedSequence(2026111401, spawn_key=(point, window, amplitude, theta, record, child[, site])),
children: site dynamics 0, lengths 1, ends 2, sherds 3 (per site), coupling 4, pool 5.

Usage: .venv-np254/bin/python analyses/per_type/type_origination.py
Outputs: output/probe_latent_groups/type_origination_<point>.csv (gitignored),
  output/findings/type_origination_cells.csv, output/findings/type_origination.md
"""
from __future__ import annotations

import itertools
import multiprocessing as mp
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT, ROOT / "analyses", ROOT / "analyses" / "probe_latent_groups", ROOT / "analyses" / "per_type"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import drift_extension as de
import readings as rd

from mls_emergence.probe import spans
from mls_emergence.signatures import ranked_readouts as rr
from mls_emergence.transmission import origination as og
from mls_emergence.transmission.spatial import kappa_for_window_f

OUTD = ROOT / "output" / "probe_latent_groups"
FINDD = ROOT / "output" / "findings"
MANIFEST_OUTPUTS = ["output/probe_latent_groups/type_origination_{point}.csv",
                    "output/findings/type_origination_cells.csv",
                    "output/findings/type_origination.md"]
ROOT_ENTROPY = 2026111401
RECORDS = 50
AMPLITUDES = (0.1, 0.3, 0.5, 0.7, 0.9)
THETAS = (0.0, 0.3, 1.0, 3.0, 10.0)
WINDOWS = ("8LR", "32LR")
SPANS = ("Ld", "1")
COUPLINGS = ("c0", "lat100N", "lat100S")
DESIGNS = tuple(itertools.product(SPANS, COUPLINGS))
N_CELLS = len(SPANS) * len(WINDOWS) * len(COUPLINGS) * len(AMPLITUDES) * len(THETAS)      # 300 per operating point
CHILD = {"sites": 0, "lengths": 1, "ends": 2, "sherds": 3, "coupling": 4, "pool": 5}
F3 = ("F_r1", "F_r2", "F_r3")
RARE = tuple(f"occ_r{r}" for r in range(4, 11))
RICH = ("K_obs", "beyond10")
SPACE3 = ("D_dist", "A_lat_r1", "A_lat_r3")
STATS = F3 + RARE + RICH + SPACE3 + ("A_lat_r2",)
CELL = ["span", "window", "coupling", "F_W", "theta"]
G: dict = {}


def pool_size(point: str, window: str, F_W: float) -> int:
    """N_R = round(kappa) + 1 for C3's kappa at this amplitude and window, so the per-step drift variance matches C3."""
    return int(round(kappa_for_window_f(F_W, de.scales(point, window)["W"]))) + 1


def realized_f(point: str, window: str, F_W: float) -> float:
    """The F the integer pool accumulates over the window: 1 - (1 - 1/N_R)^W."""
    return 1.0 - (1.0 - 1.0 / pool_size(point, window, F_W)) ** de.scales(point, window)["W"]


def seed_key(point: str, window: str, a: int, h: int, r: int, child: str, site: int | None = None) -> np.random.SeedSequence:
    if point not in de.POINT_INDEX or window not in WINDOWS or child not in CHILD:
        raise ValueError(f"unknown point {point!r}, window {window!r} or stream {child!r}")
    if not (0 <= a < len(AMPLITUDES) and 0 <= h < len(THETAS) and 0 <= r < RECORDS):
        raise ValueError("amplitude, theta or record index outside the registered grid")
    if (child == "sherds") != (site is not None):
        raise ValueError("the sherds stream takes a site index and no other stream does")
    key = (de.POINT_INDEX[point], WINDOWS.index(window), a, h, r, CHILD[child]) + (() if site is None else (int(site),))
    return np.random.SeedSequence(ROOT_ENTROPY, spawn_key=key)


def record_designs(ctx: dict, point: str, window: str, a: int, h: int, r: int) -> dict:
    """Span ends and lengths of the six designs of one record, from base draws on fresh generators."""
    sc = de.scales(point, window); n = ctx["N"]
    rng = lambda child: np.random.default_rng(seed_key(point, window, a, h, r, child))
    v = rng("lengths").random(n); ends0 = spans.ends_window(rng("ends"), n, sc["S"], sc["W"]); u = rng("coupling").random(n)
    out = {}
    for span, coupling in DESIGNS:
        lengths = spans.lengths_from_uniforms(v, de.span_limit(span, sc["L_d"]))
        ends = ends0.copy() if coupling == "c0" else spans.ends_coupled(ends0, ctx["lat"], 1.0, coupling == "lat100N", u)
        out[(span, coupling)] = (ends, lengths)
    return out


def one_record(ctx: dict, point: str, window: str, a: int, h: int, r: int) -> list[dict]:
    """One simulated history observed under the six designs: one row of ranked readouts per design."""
    sc = de.scales(point, window); op = ctx["ops"][point]; designs = record_designs(ctx, point, window, a, h, r)
    N = pool_size(point, window, AMPLITUDES[a])
    sim = og.origination_profiles(ctx["weights"][point], n_ind=op["n_ind"], mixing=op["mixing"], innovation=op["innovation"],
                                  steps=sc["S"], burnin=sc["burnin"],
                                  seed=int(np.random.default_rng(seed_key(point, window, a, h, r, "sites")).integers(0, 2**31 - 1)),
                                  pool_seed=seed_key(point, window, a, h, r, "pool"), target=ctx["pooled"], pool_size=N,
                                  theta=THETAS[h], spans={f"{s}|{c}": v for (s, c), v in designs.items()}, max_slots=8192)
    rows = []
    for (span, coupling) in DESIGNS:
        name = f"{span}|{coupling}"; prof = sim["profiles"][name]; origin = sim["origin"][name]
        y = np.array([np.random.default_rng(seed_key(point, window, a, h, r, "sherds", i)).multinomial(int(n), q)
                      for i, (n, q) in enumerate(zip(ctx["totals"], prof))])
        top = rr.rank_order(y)[:rr.N_RANKED]
        rows.append({"point": point, "span": span, "window": window, "coupling": coupling, "F_W": AMPLITUDES[a], "theta": THETAS[h],
                     "record": r, "pool_size": N, **rr.ranked_readouts(y, ctx["lat"], ctx["D"]),
                     "new_in_top3": int((origin[top[:3]] >= 0).sum()), "new_in_top10": int((origin[top] >= 0).sum()),
                     **{k: sim["info"][k] for k in ("initial_present", "originated", "pool_types_end", "max_live_slots")}})
    return rows


def run_task(task) -> list[dict]:
    point, window, a, h = task
    return [row for r in range(RECORDS) for row in one_record(G["ctx"], point, window, a, h, r)]


def stat_values(df: pd.DataFrame, stat: str) -> np.ndarray:
    if stat.startswith("A_lat_"):
        return np.abs(pd.to_numeric(df["lat_" + stat[6:]], errors="coerce").to_numpy(float))
    return pd.to_numeric(df[stat], errors="coerce").to_numpy(float)


def record_value(real: dict, stat: str) -> float:
    return abs(float(real["lat_" + stat[6:]])) if stat.startswith("A_lat_") else float(real[stat])


def cell_table(per_record: pd.DataFrame, real: dict) -> pd.DataFrame:
    rows = []
    for key, cell in per_record.groupby(CELL, sort=False):
        if len(cell) != RECORDS:
            raise ValueError(f"cell {key} has {len(cell)} records, expected {RECORDS}")
        head = dict(zip(CELL, key)) | {"point": cell["point"].iloc[0], "pool_size": int(cell["pool_size"].iloc[0])}
        for s in STATS:
            b = rd.band(stat_values(cell, s)); x = record_value(real, s)
            rows.append({**head, "stat": s, "record_value": x, "lo": b["lo"], "median": b["median"], "hi": b["hi"],
                         "n_finite": b["n_finite"], "status": rd.status3(x, b)})
        for s in ("new_in_top3", "new_in_top10", "originated", "initial_present", "pool_types_end", "max_live_slots"):
            b = rd.band(cell[s]); rows.append({**head, "stat": s, "record_value": np.nan, "lo": b["lo"], "median": b["median"],
                                               "hi": b["hi"], "n_finite": b["n_finite"], "status": "not compared"})
    return pd.DataFrame(rows)


def sets(table: pd.DataFrame) -> pd.DataFrame:
    """Per cell: membership of L (level), R (all seven rare ranks), R5 (at least five), K (richness), S3 (space), and counts."""
    t = table[table["stat"].isin(STATS)]
    w = t.pivot(index=CELL, columns="stat", values="status")
    if len(w) != N_CELLS or w[list(STATS)].isna().any().any():
        raise ValueError(f"the sets are read on one operating point's {N_CELLS} cells with every statistic; the table has {len(w)}")
    ins = w == "inside"
    out = pd.DataFrame({"L": ins[list(F3)].all(axis=1), "n_F": ins[list(F3)].sum(axis=1), "n_rare": ins[list(RARE)].sum(axis=1),
                        "K": ins[list(RICH)].all(axis=1), "S3": ins[list(SPACE3)].all(axis=1)})
    out["R"] = out["n_rare"] == len(RARE); out["R5"] = out["n_rare"] >= 5
    return out


def statement(s: pd.DataFrame) -> tuple[str, str]:
    """The statement form of spec section 7, read on R (all seven rare ranks inside), and its text.

    (a) L and R together is empty; (b) it is non-empty and no such cell is
    also in K and S3; (c) some cell is in all four sets. The five-of-seven
    threshold R5 does not enter the form (Codex build audit 2026-11-14, F4):
    it is reported by `r5_summary` as a separate, post hoc summary.
    """
    lr = s["L"] & s["R"]; full = lr & s["K"] & s["S3"]
    if not lr.any():
        return "a", ("At no tested origination intensity do the level of variation and the occupancy of the rare classes both lie "
                     "inside their bands.")
    if not full.any():
        fails_k, fails_s = int((lr & ~s["K"]).sum()), int((lr & ~s["S3"]).sum())
        return "b", (f"The level of variation and the occupancy of all seven rare ranks both lie inside their bands at {int(lr.sum())} "
                     f"cells. None of those cells also has richness and the three spatial readouts inside their bands: richness is "
                     f"outside at {fails_k} of them and space at {fails_s}.")
    return "c", (f"The level of variation, the occupancy of all seven rare ranks, richness and the three spatial readouts each lie "
                 f"inside their own band at {int(full.sum())} cells.")


def r5_summary(s: pd.DataFrame) -> str:
    """The same reading with at least five of the seven rare ranks, a threshold chosen after the drift-extension run."""
    lr5 = s["L"] & s["R5"]
    return (f"With the threshold relaxed to at least five of the seven rare ranks (chosen after the drift-extension run, so not part "
            f"of the statement form): level and rare ranks together at {int(lr5.sum())} cells, of which {int((lr5 & s['K']).sum())} "
            f"also have richness inside, {int((lr5 & s['S3']).sum())} also have the three spatial readouts inside, and "
            f"{int((lr5 & s['K'] & s['S3']).sum())} have both. The cells are listed below.")


def initial_pool(ctx: dict, point: str, window: str, F_W: float) -> tuple[int, int, list[str]]:
    """Pool size, number of starting classes with at least one draw, and the classes rounded to zero draws."""
    N = pool_size(point, window, F_W); c0 = og.largest_remainder(ctx["pooled"], N)
    return N, int((c0 > 0).sum()), [str(ctx["classes"][k]) for k in np.flatnonzero(c0 == 0)]


def _b(row, nd=3, signed=False) -> str:
    f = (lambda x: f"{x:+.{nd}f}") if signed else (lambda x: f"{x:.{nd}f}")
    return f"{f(row['median'])} [{f(row['lo'])}, {f(row['hi'])}] {row['status']}"


def point_section(point: str, t: pd.DataFrame, real: dict) -> list[str]:
    s = sets(t); theta = s.index.get_level_values("theta"); L = [f"## {point.capitalize()} point", ""]
    L += ["Pool sizes N_R (realized F over the window in parentheses): " + "; ".join(
        f"{w} " + ", ".join(f"F_W {a:g}: {pool_size(point, w, a)} ({realized_f(point, w, a):.3f})" for a in AMPLITUDES) for w in WINDOWS) + ".", "",
          "### Starting classes in the pool", "",
          ("The pool starts from the record's pooled shares rounded to N_R draws by largest remainder. A class allocated zero draws "
           "is absent from the start and, since a lost type cannot return, is absent from every record of that cell at every theta. This "
           "is absence from the start, distinct from loss by drift during the window. Where a class is listed, a theta = 0 record "
           "has at most the stated number of classes, so the occupancy of the lowest ranks is zero by construction."), "",
          "| window | F_W | N_R | starting classes present | classes rounded to zero draws |", "|---|---|---|---|---|"]
    for w in WINDOWS:
        for a in AMPLITUDES:
            N, k0, gone = initial_pool(G["ctx"], point, w, a)
            L.append(f"| {w} | {a:g} | {N} | {k0} | {', '.join(gone) if gone else 'none'} |")
    L += ["", "### Sets by origination intensity", "",
          ("Cells of 60 per theta at which the record's summaries each lie inside their own band. L: F of ranks 1 to 3. R: sites "
           "present for all of ranks 4 to 10. R5: at least five of those seven (a threshold chosen after the drift-extension run). "
           "K: richness and the share beyond rank 10. S3: D_dist and the latitude trends of ranks 1 and 3."), "",
          "| theta | L | R | R5 | K | S3 | L and R | all of L, R, K, S3 | L and R5 | L and R5 and K | L and R5 and S3 | all of L, R5, K, S3 | histories behind L and R5 |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for h in THETAS:
        m = s[theta == h]; lr5 = m["L"] & m["R5"]
        hist = RECORDS * len({(c[1], c[3]) for c in m.index[lr5]})
        L.append(f"| {h:g} | {int(m['L'].sum())} | {int(m['R'].sum())} | {int(m['R5'].sum())} | {int(m['K'].sum())} | {int(m['S3'].sum())} | "
                 f"{int((m['L'] & m['R']).sum())} | {int((m['L'] & m['R'] & m['K'] & m['S3']).sum())} | {int(lr5.sum())} | {int((lr5 & m['K']).sum())} | {int((lr5 & m['S3']).sum())} | "
                 f"{int((lr5 & m['K'] & m['S3']).sum())} | {hist} |")
    form, text = statement(s)
    lr5 = s["L"] & s["R5"]
    L += ["", f"**Statement (form {form}).** {text}", "", r5_summary(s), "",
          ("Counts are of grid cells from one seed set. The six designs of a window, amplitude and theta share simulated histories, and a "
           "band excludes about one record in five by construction, so a small count can arise from the bands' own sampling."), ""]
    listed = s.index[lr5 | (s["L"] & s["S3"])]
    L += ["### Cells with level and at least five rare ranks, or level and space", ""]
    if len(listed) == 0:
        L += ["None.", ""]
    else:
        L += ["| span | window | dates | F_W | theta | rare ranks inside | richness: median [band] (record 10) | beyond rank 10: median | K | S3 | new types among top 3 / top 10 (median) |",
              "|---|---|---|---|---|---|---|---|---|---|---|"]
        for c in listed:
            q = t[(t[CELL] == pd.Series(dict(zip(CELL, c)))).all(axis=1)].set_index("stat")
            L.append(f"| {c[0]} | {c[1]} | {c[2]} | {c[3]:g} | {c[4]:g} | {int(s.loc[c, 'n_rare'])} | {q.loc['K_obs', 'median']:.0f} "
                     f"[{q.loc['K_obs', 'lo']:.0f}, {q.loc['K_obs', 'hi']:.0f}] | {q.loc['beyond10', 'median']:.3f} | "
                     f"{'yes' if s.loc[c, 'K'] else 'no'} | {'yes' if s.loc[c, 'S3'] else 'no'} | "
                     f"{q.loc['new_in_top3', 'median']:.0f} / {q.loc['new_in_top10', 'median']:.0f} |")
        L += [""]
    L += ["### Baseline design (span limit L_d, window 8 L_R, dates independent of location) by amplitude and theta", "",
          "Cells: median [10th, 90th percentile] status. Rare ranks: number of ranks 4 to 10 whose record occupancy is inside its band.", "",
          f"| F_W | theta | F rank 1 (record {real['F_r1']:.3f}) | F rank 2 (record {real['F_r2']:.3f}) | F rank 3 (record {real['F_r3']:.3f}) | "
          f"rare ranks inside | sites, rank 4 (record {real['occ_r4']}) | sites, rank 7 (record {real['occ_r7']}) | sites, rank 10 (record {real['occ_r10']}) | "
          "richness (record 10) | beyond rank 10 | types originated | new among top 10 |", "|" + "---|" * 13]
    for a in AMPLITUDES:
        for h in THETAS:
            c = ("Ld", "8LR", "c0", a, h); q = t[(t[CELL] == pd.Series(dict(zip(CELL, c)))).all(axis=1)].set_index("stat")
            L.append(f"| {a:g} | {h:g} | " + " | ".join(_b(q.loc[x]) for x in F3) + f" | {int(s.loc[c, 'n_rare'])} | "
                     + " | ".join(_b(q.loc[x], nd=0) for x in ("occ_r4", "occ_r7", "occ_r10")) + f" | {_b(q.loc['K_obs'], nd=0)} | "
                     f"{_b(q.loc['beyond10'])} | {q.loc['originated', 'median']:.0f} | {q.loc['new_in_top10', 'median']:.0f} |")
    L += ["", "### Calibration: theta = 0 beside generator C3", "",
          ("The comparison of C4 at theta = 0 with C3 on identical ranked readouts, with medians and bands for every statistic, is "
           "in `output/findings/type_origination_calibration.md` (`analyses/per_type/type_origination_calibration.py`).")]
    return L + [""]


def main() -> int:
    de.require_numpy()
    rev = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()
    ctx = de.context(); G["ctx"] = ctx; real = rr.ranked_readouts(ctx["y"], ctx["lat"], ctx["D"])
    OUTD.mkdir(parents=True, exist_ok=True); tables, secs = [], {}
    for point in de.POINT_INDEX:
        tasks = [(point, w, a, h) for w in WINDOWS for a in range(len(AMPLITUDES)) for h in range(len(THETAS))]
        t0 = time.time()
        with mp.get_context("fork").Pool(min(mp.cpu_count(), len(tasks))) as pool:
            rows = [r for rr_ in pool.map(run_task, tasks, chunksize=1) for r in rr_]
        secs[point] = time.time() - t0
        per_record = pd.DataFrame(rows); per_record.to_csv(OUTD / f"type_origination_{point}.csv", index=False)
        tables.append(cell_table(per_record, real))
    pd.concat(tables, ignore_index=True).to_csv(FINDD / "type_origination_cells.csv", index=False)
    L = ["# Type origination during the occupation (generator C4): exploratory", "",
         (f"Produced by `.venv-np254/bin/python analyses/per_type/type_origination.py` at commit {rev[:12]}; NumPy {np.__version__}; data hash "
          f"{ctx['data_hash']}; run time " + ", ".join(f"{p} point {secs[p]:.0f} s" for p in secs) + ". Design: "
          "`docs/superpowers/specs/2026-11-14-type-origination-generator-design.md`, draft 1. Exploratory by the author's decision of "
          "2026-11-14. The spec was not audited. The build and the first run were audited by Codex on 2026-11-14, after the run (`docs/audits/codex_audit_2026-11-14_build_ledger.md`). Each cell is 50 records from one seed set."), "",
         ("Generator C4 is the neutral copying model of C3 with the regional innovation target given by a finite pool, resampled at every "
          "stored step, in which each draw is a new type with probability theta / (2 N_R). theta = 0 has no origination. The Ewens "
          "estimate from the record's richness (ten classes in 14,101 sherds) is theta = 0.98, a reference point only, since "
          f"archaeological types are not neutral variants. Grid per point: {len(SPANS)} span limits x {len(WINDOWS)} windows x "
          f"{len(COUPLINGS)} date arrangements x {len(AMPLITUDES)} amplitudes x {len(THETAS)} values of theta = {N_CELLS} cells."), "",
         ("Readouts are by rank of pooled count, because a simulated new type has no counterpart among the record's named classes. "
          "Every comparison places the record's value beside the 10th to 90th percentile band of the cell's simulated records. These "
          "are descriptive bands, not posterior intervals, and every statement is of marginal coverage: each named summary lies inside "
          "its own band. Nothing here shows that types originated during the occupation at any intensity."), "",
         "## The record, by rank", "",
         f"F of ranks 1 to 3: {real['F_r1']:.3f}, {real['F_r2']:.3f}, {real['F_r3']:.3f}. Sites present for ranks 4 to 10: "
         + ", ".join(str(real[f"occ_r{r}"]) for r in range(4, 11)) + f". Richness {real['K_obs']}; share beyond rank 10: {real['beyond10']:.3f}. "
         f"D_dist {real['D_dist']:+.3f}; latitude rank correlations of ranks 1 to 3: " + ", ".join(f"{real[f'lat_r{r}']:+.3f}" for r in (1, 2, 3)) + ".", ""]
    for point, t in zip(de.POINT_INDEX, tables):
        L += point_section(point, t, real)
    out = FINDD / "type_origination.md"; out.write_text("\n".join(L)); print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Regional drift at larger amplitudes and other occupation timings (drift-extension spec 2026-11-14, draft 3).

Generator C3 (`drift_record` with a drifting innovation target) is used
unchanged. Per operating point the grid is 3 span-length limits x 3 windows x
7 couplings of span end to location = 63 observation designs, each at 8
amplitudes F_W: 504 cells, 50 records each. The 21 designs that share a window
read the same simulated histories (3 x 8 x 50 = 1,200 per point), so cells that
share a window are dependent.

Designs. Span limit L in (L_d, ceil(L_d / 4), 1), length 1 + floor(v L). Window
W in (8, 2, 32) L_R, with F_W the regional F accumulated over that W, so kappa
is recomputed per window. Coupling: c0 (ends as drawn, independent of
location); lat50N, lat100N, lat50S, lat100S (the same ends reassigned along
latitude with coupling 0.5 or 1, later in the north or in the south); stepL,
stepS (the two spatial clusters in disjoint halves of the window, the larger
cluster late or early).

Seeds. Every stream is a fresh generator of
SeedSequence(2026111400, spawn_key=(point, window, amplitude, record, child[, site])),
children history 0, lengths 1, ends 2, sherds 3 (one stream per site),
coupling 4, step 5. No generator is carried between designs: designs share
only the base draws of a record, so an observation is the same whether it is
generated alone or with others, in any order.

Readouts per record: the carried per-class statistics of
`per_class_table.per_class_rows`, losses from the target, profile F of the
target and of the span profiles, and the two plug-in spatial readouts of
`mls_emergence.signatures.spatial_readouts` (lat_rank per class, D_dist).
Each cell's band is the 10th to 90th percentile of the finite values, eligible
at 45 finite values of 50; the record's status against it is inside, above,
below or unavailable. These are descriptive neutral comparisons (rule 18), not
posterior intervals. Reading rules R1 to R5 are in drift_extension_readings.py.

Environment. Frozen to NumPy 2.5.4, the version of the committed C3 tables
(`require_numpy`). Before the run the script regenerates the committed C3 cell
at F_W = 0.3 of each point through `per_class_table.one_record_c3m` and raises
unless its three common-class bands equal the committed findings file to
three decimals (`reproduce_committed`): a check of the environment and of the
registered frozen time scales, added in the build beyond the spec.

Usage: .venv-np254/bin/python analyses/per_type/drift_extension.py
Outputs: output/probe_latent_groups/drift_extension_<point>.csv (one row per record and design; gitignored),
  output/findings/drift_extension_cells.csv (one row per cell and statistic), output/findings/drift_extension.md
"""
from __future__ import annotations

import hashlib
import importlib
import itertools
import json
import math
import multiprocessing as mp
import re
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

import drift_extension_readings as R
import per_class_table as pct
import readings as rd

from mls_emergence.dataio.basin import geo_km, load_curated
from mls_emergence.probe import generators as gen
from mls_emergence.probe import resemblance as rs
from mls_emergence.probe import spans
from mls_emergence.signatures import spatial_readouts as sr
from mls_emergence.transmission.spatial import drift_record, kappa_for_window_f

OUTD = ROOT / "output" / "probe_latent_groups"
FINDD = ROOT / "output" / "findings"
MANIFEST_OUTPUTS = ["output/probe_latent_groups/drift_extension_{point}.csv",
                    "output/findings/drift_extension_cells.csv",
                    "output/findings/drift_extension.md"]
LOCK = Path(__file__).with_name("drift_extension_requirements.lock")
POINTS_FILE = ROOT / "analyses" / "probe_latent_groups" / "operating_points.json"
NUMPY_REQUIRED = "2.5.4"
ROOT_ENTROPY = 2026111400
RECORDS = R.RECORDS_PER_CELL
POINT_INDEX = {"primary": 0, "secondary": 1}
# Frozen time scales, registered from line 5 of the committed findings files
# output/findings/per_type_generators_C3_<point>.md (spec 3.1). Not re-measured here.
FROZEN = {"primary": {"L_d": 182, "L_R": 292, "burnin": 5000}, "secondary": {"L_d": 12, "L_R": 40, "burnin": 5000}}
TS_FILES = {p: OUTD / f"time_scales_{p}.json" for p in POINT_INDEX}
AMPLITUDES = (0.0, 0.01, 0.03, 0.1, 0.3, 0.5, 0.7, 0.9)
EARLIER_AMPLITUDES = (0.0, 0.01, 0.03, 0.1, 0.3)          # the registered grid of the 2026-10-19 spec
SPANS = R.SPANS
WINDOWS = R.WINDOWS
WINDOW_MULTIPLE = {"8LR": 8, "2LR": 2, "32LR": 32}        # the window index of the seed key is the position in WINDOWS
COUPLINGS = R.COUPLINGS
LATITUDE_COUPLING = {"lat50N": (0.5, True), "lat100N": (1.0, True), "lat50S": (0.5, False), "lat100S": (1.0, False)}
STEP_COUPLING = {"stepL": "large_label", "stepS": "small_label"}   # which cluster is late
DESIGNS = tuple(itertools.product(SPANS, WINDOWS, COUPLINGS))
D0 = R.D0
CHILD = {"history": 0, "lengths": 1, "ends": 2, "sherds": 3, "coupling": 4, "step": 5}
COMMON = pct.COMMON
G: dict = {}


def require_numpy(version: str | None = None) -> None:
    """Raise unless NumPy is the registered version: 2.4.6 and 2.5.4 give different multinomial draws (spec 3.1)."""
    v = np.__version__ if version is None else version
    if v != NUMPY_REQUIRED:
        raise RuntimeError(f"this run is frozen to NumPy {NUMPY_REQUIRED}; found {v}. Use .venv-np254 "
                           f"({LOCK.relative_to(ROOT)}).")


def amplitude_index(F_W: float) -> int:
    """Position of F_W on the registered amplitude grid; anything else raises."""
    F = float(F_W)
    hits = [i for i, x in enumerate(AMPLITUDES) if np.isfinite(F) and abs(F - x) <= 1e-12]
    if not hits:
        raise ValueError(f"F_W {F_W!r} is not on the registered amplitude grid {AMPLITUDES}")
    return hits[0]


def scales(point: str, window: str) -> dict:
    """L_d, L_R, burn-in, the window W and the stored length S = W + L_d rounded up to 100 (to 10 below 1,000).

    The frozen values are constants. When the gitignored time-scale file of the
    point is present it must agree, or this raises.
    """
    if point not in FROZEN or window not in WINDOW_MULTIPLE:
        raise ValueError(f"unknown operating point {point!r} or window {window!r}")
    f = FROZEN[point]
    if TS_FILES[point].exists():
        d = json.loads(TS_FILES[point].read_text())
        found = {"L_d": int(d["L_d"]), "L_R": int(d["L_R"]), "burnin": int(d["burnin_used"])}
        if found != f:
            raise ValueError(f"{TS_FILES[point]} gives {found}, the registered frozen scales are {f}")
    W = WINDOW_MULTIPLE[window] * f["L_R"]
    unit = 100 if W + f["L_d"] >= 1000 else 10
    return {**f, "W": W, "S": int(math.ceil((W + f["L_d"]) / unit) * unit)}


def span_limit(span: str, L_d: int) -> int:
    """The span-length limit L of a design: L_d, ceil(L_d / 4) or 1."""
    if span not in SPANS:
        raise ValueError(f"unknown span limit {span!r}; the designs use {SPANS}")
    return {"Ld": int(L_d), "Ld/4": int(math.ceil(int(L_d) / 4)), "1": 1}[span]


def check_fit(W: int, L: int, S: int) -> None:
    """A span of the longest length ending at the first step of the window must start inside the stored record."""
    if W + L > S:
        raise ValueError(f"window {W} plus the longest span {L} exceeds the stored length {S}")


def kappa(point: str, window: str, a: int):
    """Dirichlet concentration for amplitude index a over this window's W: None at F_W = 0."""
    return kappa_for_window_f(AMPLITUDES[a], scales(point, window)["W"])


def seed_key(point: str, window: str, a: int, r: int, child: str, site: int | None = None) -> np.random.SeedSequence:
    """SeedSequence(2026111400, spawn_key=(point, window, amplitude, record, child[, site])) (spec 5)."""
    if point not in POINT_INDEX or window not in WINDOWS or child not in CHILD:
        raise ValueError(f"unknown point {point!r}, window {window!r} or stream {child!r}")
    if not 0 <= int(a) < len(AMPLITUDES) or not 0 <= int(r) < RECORDS:
        raise ValueError(f"amplitude index {a} or record {r} is outside the registered grid")
    key = (POINT_INDEX[point], WINDOWS.index(window), int(a), int(r), CHILD[child])
    if (child == "sherds") != (site is not None):
        raise ValueError("the sherds stream takes a site index and no other stream does")
    if site is not None:
        if not 0 <= int(site) < 10000:
            raise ValueError("site index out of range")
        key += (int(site),)
    return np.random.SeedSequence(ROOT_ENTROPY, spawn_key=key)


def base_draws(point: str, window: str, a: int, r: int, n_sites: int) -> dict:
    """The record's base draws, each from a fresh generator of its own key: the only randomness designs share."""
    sc = scales(point, window)
    rng = lambda child: np.random.default_rng(seed_key(point, window, a, r, child))
    return {"sim_seed": int(rng("history").integers(0, 2**31 - 1)),
            "v": rng("lengths").random(n_sites),
            "ends0": spans.ends_window(rng("ends"), n_sites, sc["S"], sc["W"]),
            "u": rng("coupling").random(n_sites),
            "w": rng("step").random(n_sites)}


def context() -> dict:
    """The record and everything computed from it once: counts, latitudes, distances, clusters, pooled shares, totals."""
    from a0_predictive_summaries import spatial_clusters
    m47 = importlib.import_module("47_revision_analysis")
    counts, coords = load_curated()
    if list(coords.columns) != ["Latitude", "Longitude"]:
        raise ValueError(f"expected coordinate columns Latitude, Longitude; found {list(coords.columns)}")
    xy = coords.to_numpy(float); y = counts.to_numpy(int); labels = np.asarray(spatial_clusters(xy))
    sizes = {int(g): int((labels == g).sum()) for g in np.unique(labels)}
    if len(sizes) != 2 or len(set(sizes.values())) != 2:
        raise ValueError(f"the step designs need two spatial clusters of different size, got {sizes}")
    pts = json.loads(POINTS_FILE.read_text()); D = geo_km(xy)
    ops = {p: {k: pts[p][k] for k in pct.OP_KEYS} for p in POINT_INDEX}
    h = hashlib.sha256(); h.update(np.ascontiguousarray(y).tobytes()); h.update("|".join(map(str, counts.index)).encode())
    h.update("|".join(map(str, counts.columns)).encode()); h.update(labels.astype(int).tobytes())
    return {"y": y, "classes": list(counts.columns), "sites": list(counts.index), "N": y.shape[0], "K": y.shape[1],
            "lat": xy[:, 0].copy(), "D": D, "labels": labels, "large_label": max(sizes, key=sizes.get),
            "small_label": min(sizes, key=sizes.get), "pooled": y.astype(float).sum(0) / y.sum(), "totals": y.sum(1),
            "diversity": m47.diversity, "ops": ops, "data_hash": h.hexdigest()[:16],
            "weights": {p: gen.family_weights(ops[p]["kernel"], ops[p]["length"], D) for p in POINT_INDEX}}


def simulate(ctx: dict, point: str, window: str, a: int, r: int):
    """The simulated history of a record key and the target snapshots at every stored step (the C3 generator, unchanged)."""
    sc = scales(point, window); op = ctx["ops"][point]
    sim_seed = int(np.random.default_rng(seed_key(point, window, a, r, "history")).integers(0, 2**31 - 1))
    return drift_record(ctx["weights"][point], k=ctx["K"], n_ind=op["n_ind"], innovation=op["innovation"], mixing=op["mixing"],
                        steps=sc["S"], n_records=sc["S"], burnin=sc["burnin"], seed=sim_seed, initial="uniform",
                        target=ctx["pooled"], target_kappa=kappa(point, window, a), return_targets=True)


def design_spans(ctx: dict, point: str, window: str, span: str, coupling: str, base: dict):
    """Span ends and lengths of one observation design, from the record's base draws alone.

    c0 reads only the window draw; the latitude designs read the window draw
    and the coupling uniforms; the step designs read only the step uniforms.
    """
    sc = scales(point, window); L = span_limit(span, sc["L_d"]); check_fit(sc["W"], L, sc["S"])
    lengths = spans.lengths_from_uniforms(base["v"], L)
    if coupling == "c0":
        ends = np.array(base["ends0"], copy=True)
    elif coupling in LATITUDE_COUPLING:
        c, north = LATITUDE_COUPLING[coupling]
        ends = spans.ends_coupled(base["ends0"], ctx["lat"], c, north, base["u"])
    elif coupling in STEP_COUPLING:
        ends = spans.ends_step(base["w"], ctx["labels"], ctx[STEP_COUPLING[coupling]], sc["S"], sc["W"])
    else:
        raise ValueError(f"unknown coupling {coupling!r}; the designs use {COUPLINGS}")
    return ends, lengths


def observe(ctx: dict, rec, targets, ends, lengths, point: str, window: str, a: int, r: int):
    """Counts for every site from its span profile, each site from a fresh generator of its own key, and the record's extras."""
    S = len(rec); prof = spans.span_profile(rec, ends, lengths)
    y = np.array([np.random.default_rng(seed_key(point, window, a, r, "sherds", i)).multinomial(int(n), q)
                  for i, (n, q) in enumerate(zip(ctx["totals"], prof))])
    return y, {"lost": targets[S - 1] == 0.0, "pfT": rs.profile_f(targets[ends]), "pfB": rs.profile_f(prof)}


def readouts(ctx: dict, y) -> dict:
    """Carried per-class statistics and the two spatial readouts of one counts matrix (the record or a simulated record)."""
    out = pct.per_class_rows(y, ctx["labels"], ctx["diversity"])
    lat, why = sr.lat_rank(y, ctx["lat"])
    out.update({f"lat_{k}": float(lat[k]) for k in range(len(lat))})
    out.update({f"lat_reason_{k}": why[k] for k in range(len(lat))})
    out["D_dist"], out["D_reason"] = sr.distance_decay_rank(y, ctx["D"])
    return out


def _generate(ctx: dict, point: str, window: str, a: int, r: int, designs):
    if any(d[1] != window for d in designs) or any(d not in DESIGNS for d in designs):
        raise ValueError("every design must be a registered design of this window")
    rec, targets = simulate(ctx, point, window, a, r)
    base = base_draws(point, window, a, r, ctx["N"])
    for d in designs:
        ends, lengths = design_spans(ctx, point, window, d[0], d[2], base)
        y, ex = observe(ctx, rec, targets, ends, lengths, point, window, a, r)
        yield d, y, ex, ends, lengths, base


def record_counts(ctx: dict, point: str, window: str, a: int, r: int, designs) -> dict:
    """Counts of one record under each named design of one window, keyed by design."""
    return {d: y for d, y, *_ in _generate(ctx, point, window, a, r, list(designs))}


def run_task(task) -> list[dict]:
    """All 50 records of one (point, window, amplitude): 21 designs each."""
    point, window, a = task; ctx = G["ctx"]; sc = scales(point, window); kap = kappa(point, window, a); rows = []
    designs = [d for d in DESIGNS if d[1] == window]
    for r in range(RECORDS):
        for d, y, ex, ends, lengths, base in _generate(ctx, point, window, a, r, designs):
            row = {"point": point, "span": d[0], "L": span_limit(d[0], sc["L_d"]), "window": window, "W": sc["W"], "S": sc["S"],
                   "coupling": d[2], "F_W": AMPLITUDES[a], "a": a, "kappa": kap if kap is not None else 0.0, "record": r,
                   "sim_seed": base["sim_seed"], "mean_end": float(ends.mean()), "mean_length": float(lengths.mean())}
            row.update(readouts(ctx, y))
            row.update({f"lost_{k}": bool(ex["lost"][k]) for k in range(ctx["K"])})
            row.update({f"pfT_{k}": float(ex["pfT"][k]) for k in range(ctx["K"])})
            row.update({f"pfB_{k}": float(ex["pfB"][k]) for k in range(ctx["K"])})
            rows.append(row)
    return rows


# ---- cells -----------------------------------------------------------------------------------------

def stat_names(K: int = 10) -> list[str]:
    """Statistics of the cell table: the 14 the rules read, then the signed latitude correlations and profile F (context)."""
    return (list(R.F3) + list(R.SP) + list(R.RARE) + [f"lat_rank_{k}" for k in range(COMMON)]
            + [f"pfT_{k}" for k in range(COMMON)] + [f"pfB_{k}" for k in range(COMMON)])


def _stat_source(stat: str):
    """(per-record column, transform, class index or None, compared against the record)."""
    m = re.fullmatch(r"(F|occ|A_lat|lat_rank|pfT|pfB)_(\d+)", stat)
    if stat == "D_dist":
        return "D_dist", None, None, True
    if not m:
        raise ValueError(f"unknown statistic {stat!r}")
    kind, k = m.group(1), int(m.group(2))
    col = {"F": f"F_{k}", "occ": f"occ_{k}", "A_lat": f"lat_{k}", "lat_rank": f"lat_{k}", "pfT": f"pfT_{k}", "pfB": f"pfB_{k}"}[kind]
    return col, (np.abs if kind == "A_lat" else None), k, kind not in ("pfT", "pfB")


def record_value(real: dict, stat: str) -> float:
    """The record's own value of a cell statistic (NaN for the profile F statistics, which the record does not have)."""
    col, f, _k, compared = _stat_source(stat)
    if not compared:
        return float("nan")
    x = float(real[col])
    return float(f(x)) if f is not None else x


def stat_row(cell: pd.DataFrame, stat: str, x: float) -> dict:
    """Band, finite count, losses and the record's status for one statistic of one cell."""
    col, f, k, compared = _stat_source(stat)
    v = pd.to_numeric(cell[col], errors="coerce").to_numpy(float)
    if f is not None:
        v = f(v)
    b = rd.band(v)
    return {"stat": stat, "record_value": x, "lo": b["lo"], "median": b["median"], "hi": b["hi"], "n_finite": b["n_finite"],
            "eligible": b["eligible"], "status": rd.status3(x, b) if compared else "not compared",
            "n_lost": int(cell[f"lost_{k}"].astype(bool).sum()) if k is not None and f"lost_{k}" in cell else float("nan")}


def cell_table(per_record: pd.DataFrame, real: dict) -> pd.DataFrame:
    """One row per cell and statistic for one operating point."""
    rows = []
    for (span, window, coupling, fw), cell in per_record.groupby(["span", "window", "coupling", "F_W"], sort=False):
        if len(cell) != RECORDS:
            raise ValueError(f"cell {(span, window, coupling, fw)} has {len(cell)} records, expected {RECORDS}")
        head = {"point": cell["point"].iloc[0], "span": span, "L": int(cell["L"].iloc[0]), "window": window, "W": int(cell["W"].iloc[0]),
                "S": int(cell["S"].iloc[0]), "coupling": coupling, "F_W": float(fw), "kappa": float(cell["kappa"].iloc[0])}
        rows += [{**head, **stat_row(cell, s, record_value(real, s))} for s in stat_names()]
    return pd.DataFrame(rows)


# ---- the committed C3 table: baseline comparison and the environment check ---------------------------

def committed_cells(point: str) -> dict:
    """The committed C3 cells of the three common classes, parsed from output/findings/per_type_generators_C3_<point>.md.

    Returns {F_W: {"F_k": {median, lo, hi, n_finite, status}, "lost": [count per class]}}.
    """
    path = FINDD / f"per_type_generators_C3_{point}.md"
    lines = [ln for ln in path.read_text().splitlines() if ln.startswith("|")]
    hdr = [c.strip() for c in lines[0].strip("|").split("|")]
    fcols = [i for i, c in enumerate(hdr) if c.startswith("F ") and c.endswith("n_finite status")]
    lost_col = hdr.index("lost from target (count per class)"); fw_col = hdr.index("F_W")
    if len(fcols) != COMMON:
        raise ValueError(f"{path} does not have {COMMON} common-class F columns")
    out = {}
    for ln in lines[2:]:
        c = [x.strip() for x in ln.strip("|").split("|")]
        if not c[0].startswith("C3 "):
            continue
        cell = {}
        for k, i in enumerate(fcols):
            m = re.fullmatch(r"(-?[\d.]+) \[(-?[\d.]+), (-?[\d.]+)\] (\d+) (\w+)( \(ineligible\))?", c[i])
            if not m:
                raise ValueError(f"{path}: cannot parse the cell {c[i]!r}")
            cell[f"F_{k}"] = {"median": float(m.group(1)), "lo": float(m.group(2)), "hi": float(m.group(3)),
                              "n_finite": int(m.group(4)), "status": m.group(5)}
        cell["lost"] = [int(x) for x in c[lost_col].split(",")]
        out[float(c[fw_col])] = cell
    return out


def reproduce_committed(ctx: dict, point: str, F_W: float = 0.3) -> dict:
    """Regenerate one committed C3 cell through the committed code path and compare it with the committed findings file.

    Uses `per_class_table.one_record_c3m` and `record_streams` (the earlier
    design's own streams), with the registered frozen scales. Agreement to the
    three decimals the committed file prints shows that this environment and
    these time scales reproduce the committed table. Raises otherwise.
    """
    sc = scales(point, "8LR"); a = pct.AMPLITUDE_GRID.index(F_W); want = committed_cells(point)[F_W]
    pct.G.update(N=ctx["N"], K=ctx["K"], D=ctx["D"], pooled=ctx["pooled"], totals=ctx["totals"], labels=ctx["labels"])
    g = {"W": sc["W"], "kappa": kappa_for_window_f(F_W, sc["W"])}
    F, lost = [], np.zeros(ctx["K"], int)
    for r in range(RECORDS):
        y, ex = pct.one_record_c3m(ctx["ops"][point], g, sc, pct.record_streams(point, a, r), None)
        F.append(pct.moment_f(y)[:COMMON]); lost += ex["lost"].astype(int)
    F = np.array(F); got = {}
    for k in range(COMMON):
        b = rd.band(F[:, k]); got[f"F_{k}"] = {"median": round(b["median"], 3), "lo": round(b["lo"], 3), "hi": round(b["hi"], 3), "n_finite": b["n_finite"]}
        w = want[f"F_{k}"]
        if any(abs(got[f"F_{k}"][q] - w[q]) > 5e-4 for q in ("median", "lo", "hi")) or b["n_finite"] != w["n_finite"]:
            raise RuntimeError(f"{point} point, F_W {F_W}: the committed C3 cell is not reproduced (class {k}: regenerated "
                               f"{got[f'F_{k}']}, committed {w}). The environment or the frozen scales differ from the committed run.")
    if lost.tolist() != want["lost"]:
        raise RuntimeError(f"{point} point, F_W {F_W}: losses from the target {lost.tolist()} differ from the committed {want['lost']}")
    return {"point": point, "F_W": F_W, "cells": got, "lost": lost.tolist()}


# ---- findings ------------------------------------------------------------------------------------------

def _f(x, nd=3) -> str:
    return "NaN" if x is None or not np.isfinite(x) else f"{x:.{nd}f}"


def _sg(x, nd=3) -> str:
    return "NaN" if x is None or not np.isfinite(x) else f"{x:+.{nd}f}"


def _cell_row(t: pd.DataFrame, cell: tuple, stat: str) -> pd.Series:
    m = (t["span"] == cell[0]) & (t["window"] == cell[1]) & (t["coupling"] == cell[2]) & np.isclose(t["F_W"], cell[3]) & (t["stat"] == stat)
    if int(m.sum()) != 1:
        raise ValueError(f"cell {cell}, statistic {stat}: expected one row, found {int(m.sum())}")
    return t[m].iloc[0]


def _band_text(row: pd.Series, signed: bool = False, nd: int = 3) -> str:
    fmt = _sg if signed else _f
    return f"{fmt(row['median'], nd)} [{fmt(row['lo'], nd)}, {fmt(row['hi'], nd)}] {int(row['n_finite'])} {row['status']}"


def _cells_listing(t: pd.DataFrame, cells: list, classes: list[str]) -> list[str]:
    """Listed cells with the R4 diagnostics beside each."""
    if not cells:
        return ["None."]
    hdr = ["span limit", "window", "coupling", "F_W", "rare classes inside (of 7)", "rare classes unavailable",
           "records lost from target: " + ", ".join(classes[:COMMON]), "finite counts: F x 3, D_dist, A_lat x 3"]
    L = ["| " + " | ".join(hdr) + " |", "|" + "---|" * len(hdr)]
    for c in cells:
        d = R.diagnostics(t, c)
        L.append("| " + " | ".join([c[0], c[1], c[2], f"{c[3]:g}", str(d["rare_inside"]), str(d["rare_unavailable"]),
                                    ", ".join(str(d["n_lost"][s]) for s in R.F3),
                                    ", ".join(str(d["n_finite"][s]) for s in R.F3 + R.SP)]) + " |")
    return L


def point_section(point: str, t: pd.DataFrame, ctx: dict, real: dict, repro: dict) -> list[str]:
    cls = ctx["classes"]; names = {s: cls[i] for i, s in enumerate(R.F3)}
    out = R.read(t); sc8 = scales(point, "8LR")
    L = [f"## {point.capitalize()} point", "",
         (f"Frozen scales: L_d {sc8['L_d']}, L_R {sc8['L_R']}, burn-in {sc8['burnin']}. Windows and stored lengths: "
          + "; ".join(f"{w} W = {scales(point, w)['W']}, S = {scales(point, w)['S']}" for w in WINDOWS)
          + ". Span limits: " + ", ".join(f"{s} = {span_limit(s, sc8['L_d'])}" for s in SPANS) + "."), "",
         (f"Environment check: the committed C3 cell at F_W {repro['F_W']:g} was regenerated through `per_class_table.one_record_c3m` "
          "and equals the committed findings file to three decimals for the three common classes ("
          + "; ".join(f"{cls[k]} {repro['cells'][f'F_{k}']['median']:.3f} [{repro['cells'][f'F_{k}']['lo']:.3f}, {repro['cells'][f'F_{k}']['hi']:.3f}]"
                      for k in range(COMMON)) + "), with the same losses from the target."), "",
         "### R1. Amplitude at the baseline design d0 (span limit L_d, window 8 L_R, coupling c0)", "",
         "Cells: median [10th, 90th percentile] finite count status. Lost: records of 50 in which the class's target share is exactly zero at the last stored step.", ""]
    hdr = ["F_W"] + [f"F {cls[k]} (record {_f(real[f'F_{k}'])})" for k in range(COMMON)] + ["lost: " + ", ".join(cls[:COMMON])] \
        + [f"D_dist (record {_sg(real['D_dist'])})"] + [f"A_lat {cls[k]} (record {_f(abs(real[f'lat_{k}']))})" for k in range(COMMON)]
    L += ["| " + " | ".join(hdr) + " |", "|" + "---|" * len(hdr)]
    for fw in AMPLITUDES:
        c = (*D0, fw)
        L.append("| " + " | ".join([f"{fw:g}"] + [_band_text(_cell_row(t, c, s)) for s in R.F3]
                                   + [", ".join(str(int(_cell_row(t, c, s)["n_lost"])) for s in R.F3)]
                                   + [_band_text(_cell_row(t, c, "D_dist"), signed=True)]
                                   + [_band_text(_cell_row(t, c, f"A_lat_{k}")) for k in range(COMMON)]) + " |")
    L += [""] + R.r1_text(out["R1"], names) + [""]
    n_el, n_l, n_es, n_s = (len(out[k]) for k in ("E_L", "L", "E_S", "S"))
    L += ["### R2. Level", "",
          (f"Cells evaluable for the three common-class F values (E_L): {n_el} of {R.N_CELLS}. Cells with all three inside their bands (L): "
           f"{n_l}, resting on {R.n_histories(out['L'])} distinct simulated histories. Of the {n_el - n_l} evaluable cells outside L, "
           f"{out['n_any_above']} have at least one value above its band and {out['n_below_only']} have none above and at least one below."), "",
          "Cells of L:", ""] + _cells_listing(t, out["L"], cls) + [""]
    L += ["### R3. Level and space", "",
          (f"Cells of L evaluable for the four spatial readouts (E_S): {n_es}. Cells with the four also inside their bands (S): {n_s}, "
           f"resting on {R.n_histories(out['S'])} distinct simulated histories. Cells of L not evaluable for space: "
           f"{len(out['L_not_evaluable_for_space'])}."), "", "Cells of S:", ""] + _cells_listing(t, out["S"], cls) + [""]
    L += ["### R5. By coupling", "", "| coupling | cells | E_L | L | E_S | S |", "|---|---|---|---|---|---|"]
    L += [f"| {r.coupling} | {r.cells} | {r.E_L} | {r.L} | {r.E_S} | {r.S} |" for r in out["R5"].itertuples()]
    L += ["", f"**Statement (form {out['form']}).** {out['statement']}", "",
          ("A band excludes about one record in five by construction, each cell carries up to fourteen statuses, and cells that share "
           "a window share simulated histories. The counts above are counts of grid cells, not of independent confirmations, and a "
           "small count can arise from the bands' own sampling."), ""]
    # Recorded expectations (spec 6): reported, not rules.
    c0 = t[(t["coupling"] == "c0") & (t["stat"] == "D_dist")]
    L += ["### Recorded expectations", "",
          (f"D_dist under c0 (72 cells): the record's value {_sg(real['D_dist'])} is above the band at {int((c0['status'] == 'above').sum())}, "
           f"inside at {int((c0['status'] == 'inside').sum())}, below at {int((c0['status'] == 'below').sum())} and unavailable at "
           f"{int((c0['status'] == 'unavailable').sum())}. Cell medians of simulated D_dist under c0 run from {_sg(c0['median'].min())} "
           f"to {_sg(c0['median'].max())}."), "",
          "Median simulated F by span limit, at window 8 L_R and coupling c0 (median over 50 records, finite count in parentheses):", ""]
    hdr = ["F_W"] + [f"{cls[k]}, L = {span_limit(s, sc8['L_d'])}" for k in range(COMMON) for s in SPANS]
    L += ["| " + " | ".join(hdr) + " |", "|" + "---|" * len(hdr)]
    for fw in AMPLITUDES:
        vals = []
        for k in range(COMMON):
            for s in SPANS:
                row = _cell_row(t, (s, "8LR", "c0", fw), f"F_{k}")
                vals.append(f"{_f(row['median'])} ({int(row['n_finite'])})")
        L.append("| " + " | ".join([f"{fw:g}"] + vals) + " |")
    L += ["", "### Baseline comparison with the committed C3 table", "",
          ("Design d0 at the five earlier amplitudes, this run beside the committed table "
           f"(`output/findings/per_type_generators_C3_{point}.md`). The streams differ by design, so no equality is expected. "
           "The committed table reports inside or outside; this run's above and below both correspond to its outside. Agreement "
           "does not show that conclusions are robust to the seed set."), ""]
    hdr = ["F_W", "class", "this run: median [band] finite, lost, status", "committed: median [band] finite, lost, status", "status agrees"]
    L += ["| " + " | ".join(hdr) + " |", "|" + "---|" * len(hdr)]
    com = committed_cells(point); differ = []
    for fw in EARLIER_AMPLITUDES:
        for k in range(COMMON):
            row = _cell_row(t, (*D0, fw), f"F_{k}"); w = com[fw][f"F_{k}"]
            mine = "inside" if row["status"] == "inside" else "unavailable" if row["status"] == "unavailable" else "outside"
            same = mine == w["status"]
            if not same:
                differ.append(f"F_W {fw:g}, {cls[k]}")
            L.append("| " + " | ".join([f"{fw:g}", cls[k],
                                        f"{_f(row['median'])} [{_f(row['lo'])}, {_f(row['hi'])}] {int(row['n_finite'])}, {int(row['n_lost'])}, {row['status']}",
                                        f"{w['median']:.3f} [{w['lo']:.3f}, {w['hi']:.3f}] {w['n_finite']}, {com[fw]['lost'][k]}, {w['status']}",
                                        "yes" if same else "NO"]) + " |")
    L += ["", "Statuses that differ between the two runs: " + ("; ".join(differ) if differ else "none") + ".", ""]
    return L


def main() -> int:
    require_numpy()
    if not LOCK.exists():
        raise SystemExit(f"the lock file {LOCK} is missing")
    lock_hash = hashlib.sha256(LOCK.read_bytes()).hexdigest()
    rev = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()
    ctx = context(); G["ctx"] = ctx; real = readouts(ctx, ctx["y"])
    t0 = time.time()
    repro = {p: reproduce_committed(ctx, p) for p in POINT_INDEX}
    t_repro = time.time() - t0
    OUTD.mkdir(parents=True, exist_ok=True); FINDD.mkdir(parents=True, exist_ok=True)
    tables, secs = [], {}
    for point in POINT_INDEX:
        tasks = [(point, w, a) for w in WINDOWS for a in range(len(AMPLITUDES))]
        t1 = time.time()
        with mp.get_context("fork").Pool(min(mp.cpu_count(), len(tasks))) as pool:
            rows = [r for rr in pool.map(run_task, tasks, chunksize=1) for r in rr]
        secs[point] = time.time() - t1
        per_record = pd.DataFrame(rows)
        if len(per_record) != R.N_CELLS * RECORDS:
            raise RuntimeError(f"{point}: {len(per_record)} record rows, expected {R.N_CELLS * RECORDS}")
        per_record.to_csv(OUTD / f"drift_extension_{point}.csv", index=False)
        tables.append(cell_table(per_record, real))
    cells = pd.concat(tables, ignore_index=True)
    cells.to_csv(FINDD / "drift_extension_cells.csv", index=False)
    cls = ctx["classes"]
    L = ["# Regional drift at larger amplitudes and other occupation timings", "",
         (f"Produced by `.venv-np254/bin/python analyses/per_type/drift_extension.py` at commit {rev[:12]}"
          + (" (with uncommitted changes to tracked files)" if dirty else "") + f". NumPy {np.__version__}; lock file "
          f"`{LOCK.relative_to(ROOT)}`, SHA-256 {lock_hash[:16]}; data hash {ctx['data_hash']}. Run time: environment check "
          f"{t_repro:.0f} s, " + ", ".join(f"{p} point {secs[p]:.0f} s" for p in POINT_INDEX) + ". "
          "Design: drift-extension spec of 2026-11-14, draft 3 (`docs/superpowers/specs/2026-11-14-drift-extension-timing-design.md`)."), "",
         (f"Generator C3 (a neutral copying model whose regional innovation target drifts), unchanged. Per operating point: "
          f"{len(SPANS)} span-length limits x {len(WINDOWS)} windows x {len(COUPLINGS)} couplings = {len(DESIGNS)} observation designs at "
          f"{len(AMPLITUDES)} amplitudes F_W = {', '.join(f'{a:g}' for a in AMPLITUDES)}: {R.N_CELLS} cells of {RECORDS} records. F_W is the "
          "regional F accumulated over the cell's own window. The 21 designs that share a window read the same simulated histories "
          f"({len(WINDOWS) * len(AMPLITUDES) * RECORDS} histories per point), so cells that share a window are dependent."), "",
         ("Every comparison is a neutral comparison in the main text's sense: the record's value of a statistic beside the 10th to 90th "
          "percentile band of the same statistic over the cell's simulated records (NumPy's linear quantile rule, inclusive endpoints; a "
          f"band is eligible at {rd.MIN_FINITE} finite values of {RECORDS}). These are descriptive bands, not posterior intervals. Statements "
          "are of marginal coverage: each named summary of the record lies inside its own band. No statement here says that any simulated "
          "record resembles the record in several summaries at once. Only a value above the band of a per-class F is read as more "
          "variation than the generator produces. Losses from the target are reported beside the cells and exclude none."), "",
         "## The record", "",
         ("Plug-in values on observed shares, computed by the functions applied to the simulated records. They are not the posterior "
          "quantities the main text reports for the record. Sign conventions: a positive latitude correlation means the class's share "
          "is higher to the north; a positive D_dist means sites farther apart differ more in composition."), "",
         "| statistic | " + " | ".join(cls[:COMMON]) + " |", "|---|" + "---|" * COMMON,
         "| F (site-weighted moment F) | " + " | ".join(_f(real[f"F_{k}"]) for k in range(COMMON)) + " |",
         "| latitude rank correlation (signed) | " + " | ".join(_sg(real[f"lat_{k}"]) for k in range(COMMON)) + " |",
         "| A_lat (absolute value) | " + " | ".join(_f(abs(real[f"lat_{k}"])) for k in range(COMMON)) + " |", "",
         f"D_dist (rank correlation of total-variation distance with distance over {ctx['N'] * (ctx['N'] - 1) // 2} pairs): {_sg(real['D_dist'])}. "
         "Sites present for the seven rare classes: " + ", ".join(f"{cls[k]} {int(real[f'occ_{k}'])}" for k in range(COMMON, ctx["K"])) + ".", ""]
    for point, t in zip(POINT_INDEX, tables):
        L += point_section(point, t, ctx, real, repro[point])
    L += ["## What these tables cannot give", "",
          ("Joint reproduction of the record by any simulated record; the sufficiency or necessity of a timing arrangement; whether the "
           "sites differ in date; type origination, which this generator does not have; identified parameters, since results hold at "
           "two operating points; robustness to the seed set, since each cell is 50 records from one set of streams; a test of "
           "boundaries; and the fitted spatial share of the A0 model, which is not computed on simulated records (spec section 10)."), "",
          "Per-cell values: `output/findings/drift_extension_cells.csv`. Per-record values: `output/probe_latent_groups/drift_extension_<point>.csv` (not committed).", ""]
    out = FINDD / "drift_extension.md"; out.write_text("\n".join(L)); print("\n".join(L)); print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

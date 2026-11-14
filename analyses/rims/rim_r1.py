"""Rim attributes R1: neutral yardstick cells for the rim attributes and the matched decoration reference, and the readings.

Spec: docs/superpowers/specs/2026-10-19-rim-attributes-r1-design.md (draft 2, approved 2026-10-19).

Cells (section 3). The fixed-target D22 copying simulator with ALL 28
basin populations transmitting (the same weights, kernel and operating
points as every earlier run) and only the 19 rim sites observed. Time
scales frozen at the step 0 measurement (per_class_table.frozen_scales:
primary L_d 182, L_R 292, S 2,600; secondary 12, 40, 340; burn-in 5,000);
two windows, T1 (ends over the last L_d) and 8LR (the last 8 L_R). Each
observed site draws its sample from its own stream.

- rim cells: two classes at the attribute's pooled frequency; each site's
  rim N is its sample size; per record the binary moment F over the 19
  sites and over the 18 without Woodlyn. Attributes read: interior bevel,
  exterior notching, lip flat, paste Bell (the last read separately).
- decoration reference cells: the ten classes at the record's pooled
  shares, each site's decorated sherd total (current counts) or its 1940
  PFG part's total; per record the three common classes' F over the 19
  and the 18.

Record key: SeedSequence((2026101970, point, kind, attribute, window,
record)), a tuple entropy that cannot coincide with the integer roots of
the earlier runs; children history 0, lengths 1, ends 2, sherds 3, the
sherd child spawned per observed site (spawn key (3, i)).

Readings (sections 4 to 7) are written by `--readings` from the cells'
CSVs and the rim and decorated tables: statuses below, inside, above or
unavailable per attribute, point and window, beside the decorated
classes' statuses on the same 19 sites; geography statistics and the
rim-decoration pair association as observed values and posterior
summaries (mls_emergence.rims.stats); the Beta(1, 1) and Dirichlet(1)
sensitivities, the lip plus or minus one rim sensitivity at Beck, Irby and
Starkley, the 18-site sensitivity, the observable prior predictives, and
the recovery program's result. Descriptive throughout (rule 18).

Usage:
  .venv/bin/python analyses/rims/rim_r1.py --points primary --kinds rim,decoration_current,decoration_pfg --records 50
  .venv/bin/python analyses/rims/rim_r1.py --points secondary --kinds rim,decoration_current,decoration_pfg --records 50
  .venv/bin/python analyses/rims/rim_r1.py --readings
Outputs: output/probe_latent_groups/rim_r1_<kind>_<point>.csv (gitignored), output/findings/rim_r1.md
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import multiprocessing as mp
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT, ROOT / "analyses", ROOT / "analyses" / "probe_latent_groups", ROOT / "analyses" / "per_type", ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))


from mls_emergence.dataio.basin import geo_km, load_curated
from mls_emergence.probe import generators as gen
from mls_emergence.probe import spans
from mls_emergence.rims import data as rdata
from mls_emergence.rims import stats as st
from mls_emergence.transmission.spatial import drift_record

OUTD = ROOT / "output" / "probe_latent_groups"
FINDD = ROOT / "output" / "findings"
FINDINGS = FINDD / "rim_r1.md"
POINTS = ROOT / "analyses" / "probe_latent_groups" / "operating_points.json"
MANIFEST_OUTPUTS = ["output/probe_latent_groups/rim_r1_{kind}_{point}.csv"]
SEED_BASE_R1 = 2026101970
POINT_INDEX = {"primary": 0, "secondary": 1}
KINDS = ("rim", "decoration_current", "decoration_pfg")
WINDOWS = ("T1", "8LR")
ATTRIBUTES = ("beveled", "notched_exterior", "lip_flat", "paste_bell")   # read by F; paste separately (spec 2)
NOT_READ = ("beaded", "noded", "wide_nodes", "notched_top")
LIP_AMBIGUOUS = ("Beck", "Irby", "Starkley")                              # file percentages sum to 99, 101, 98 (spec 2)
SMALL_SITE = "Woodlyn"                                                    # below Mainfort's 45-rim minimum (spec 2)
COMMON = 3
OP_KEYS = ("kernel", "length", "n_ind", "innovation", "mixing")
G: dict = {}


def scales_for(point: str) -> dict:
    pct = importlib.import_module("per_class_table")
    return pct.frozen_scales(point)


def window_steps(window: str, sc: dict) -> int:
    if window == "T1":
        return int(sc["L_d"])
    if window == "8LR":
        return 8 * int(sc["L_R"])
    raise ValueError(f"unknown window {window!r}; the spec runs T1 and 8LR")


def record_streams(point: str, kind: str, attribute_index: int, window_index: int, record: int) -> dict:
    """The record key (spec 3): a tuple entropy and four named children; site i's stream is spawn key (3, i)."""
    if point not in POINT_INDEX or kind not in KINDS or window_index not in (0, 1) or attribute_index < 0 or record < 0:
        raise ValueError("record key out of range")
    entropy = (SEED_BASE_R1, POINT_INDEX[point], KINDS.index(kind), int(attribute_index), int(window_index), int(record))
    children = np.random.SeedSequence(entropy).spawn(4)
    return {"entropy": entropy, "children": children, "sim_seed": int(np.random.default_rng(children[0]).integers(0, 2**31 - 1))}


def _site_rngs(kids: dict, n_observed: int):
    return [np.random.default_rng(np.random.SeedSequence(kids["entropy"], spawn_key=(3, i))) for i in range(n_observed)]


def _spans_and_profile(rec: np.ndarray, sc: dict, W: int, kids: dict):
    S, L_d = sc["S"], sc["L_d"]
    if W + L_d > S:
        raise ValueError(f"window {W} plus the longest span {L_d} exceeds the stored length {S}")
    rng_len, rng_end = (np.random.default_rng(k) for k in kids["children"][1:3])
    lengths = spans.draw_lengths(rng_len, G["N"], L_d)
    ends = spans.ends_window(rng_end, G["N"], S, W)
    prof = spans.span_profile(rec, ends, lengths)
    return lengths, ends, prof


def one_binary_record(op: dict, p_attr: float, sc: dict, W: int, kids: dict):
    """One rim cell record: a two-class simulator at pooled frequency p_attr over all populations, observed at the rim sites."""
    if not 0 < p_attr < 1:
        raise ValueError("the attribute's pooled frequency must lie strictly inside (0, 1)")
    w = gen.family_weights(op["kernel"], op["length"], G["D"])
    rec = drift_record(w, k=2, n_ind=op["n_ind"], innovation=op["innovation"], mixing=op["mixing"], steps=sc["S"], n_records=sc["S"],
                       burnin=sc["burnin"], seed=kids["sim_seed"], initial="uniform", target=np.array([p_attr, 1 - p_attr]))
    lengths, ends, prof = _spans_and_profile(rec, sc, W, kids)
    obs = G["observed"]; totals = np.asarray(G["rim_totals"], int)
    y = np.array([r.multinomial(int(n), q) for r, n, q in zip(_site_rngs(kids, len(obs)), totals, prof[obs])])
    return y, {"weights_shape": w.shape, "profile_shape": prof.shape, "lengths": lengths[obs], "ends": ends[obs]}


def one_reference_record(op: dict, pooled: np.ndarray, totals: np.ndarray, sc: dict, W: int, kids: dict):
    """One decoration reference record: the ten classes at the pooled shares, observed at the rim sites with the given sherd totals."""
    w = gen.family_weights(op["kernel"], op["length"], G["D"])
    rec = drift_record(w, k=G["K"], n_ind=op["n_ind"], innovation=op["innovation"], mixing=op["mixing"], steps=sc["S"], n_records=sc["S"],
                       burnin=sc["burnin"], seed=kids["sim_seed"], initial="uniform", target=np.asarray(pooled, float))
    lengths, ends, prof = _spans_and_profile(rec, sc, W, kids)
    obs = G["observed"]; totals = np.asarray(totals, int)
    y = np.array([r.multinomial(int(n), q) for r, n, q in zip(_site_rngs(kids, len(obs)), totals, prof[obs])])
    return y, {"weights_shape": w.shape, "profile_shape": prof.shape, "lengths": lengths[obs], "ends": ends[obs]}


def record_f_binary(y: np.ndarray) -> float:
    return float(st.moment_f(y)[0])


def record_f_common(y: np.ndarray) -> np.ndarray:
    return st.moment_f(y)[:COMMON]


def output_stem(kind: str, point: str) -> str:
    if kind not in KINDS or point not in POINT_INDEX:
        raise ValueError("unknown kind or point")
    return f"rim_r1_{kind}_{point}"


# ---------------------------------------------------------------------------
# Data setup and cells
# ---------------------------------------------------------------------------

def setup(records: int) -> dict:
    """Load the 28-site decorated matrix, the 19 rim sites, totals, coordinates and provenance into G."""
    cbc = importlib.import_module("composition_by_collection")
    counts, coords = load_curated(); xy = coords.to_numpy(float); D = geo_km(xy)
    rims, unmatched = rdata.basin_rim_table(list(counts.index))
    observed = np.array([list(counts.index).index(a) for a in rims.index])
    _pfg_counts, pfg, _add, _cls, _prov = cbc.decompose()
    y_real = counts.to_numpy(int)
    G.update(N=len(counts), K=counts.shape[1], D=D, classes=list(counts.columns), site_names=list(counts.index), observed=observed,
             rim_sites=list(rims.index), rims=rims, unmatched=unmatched, counts=counts, pfg=pfg,
             counts19=counts.loc[rims.index], pfg19=pfg.loc[rims.index],
             pooled=y_real.astype(float).sum(0) / y_real.sum(), pooled_pfg=pfg.to_numpy(float).sum(0) / pfg.to_numpy(float).sum(),
             rim_totals=rims["N"].to_numpy(int), sherd_totals=counts.loc[rims.index].to_numpy(int).sum(1),
             pfg_totals=pfg.loc[rims.index].to_numpy(int).sum(1), lat=coords["Latitude"].to_numpy(float)[observed],
             lon=coords["Longitude"].to_numpy(float)[observed], D_obs=D[np.ix_(observed, observed)], records=records)
    h = hashlib.sha256(); h.update(np.ascontiguousarray(y_real).tobytes()); h.update("|".join(counts.index).encode())
    num = rims.select_dtypes("number")  # the table also carries the file name and rough coordinates as text
    h.update(np.ascontiguousarray(num.to_numpy(float)).tobytes()); h.update("|".join(num.columns).encode()); h.update("|".join(rims.index).encode())
    sh = hashlib.sha256()
    for p in POINT_INDEX:
        sh.update((OUTD / f"time_scales_{p}.json").read_bytes())
    G["prov"] = {"git_rev": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip(),
                 "numpy_version": np.__version__, "bit_generator": "PCG64", "scales_hash": sh.hexdigest()[:16], "data_hash": h.hexdigest()[:16],
                 "initial": "uniform (every population at the pooled shares)", "population_model": "28 populations transmit, 19 observed"}
    return G


def attribute_frequency(attr: str) -> float:
    r = G["rims"]; return float(r[attr].sum() / r["N"].sum())


def run_cell(args):
    point, kind, attr, window, op = args
    sc = scales_for(point); W = window_steps(window, sc); a_ix = ATTRIBUTES.index(attr) if attr else 0; w_ix = WINDOWS.index(window)
    rows = []; obs_names = G["rim_sites"]; keep18 = np.array([s != SMALL_SITE for s in obs_names])
    for r in range(G["records"]):
        kids = record_streams(point, kind, a_ix, w_ix, r)
        base = {"point": point, "kind": kind, "attribute": attr or "", "window": window, "W": W, "record": r, "entropy": str(kids["entropy"]),
                "child_spawn_keys": ";".join(str(c.spawn_key) for c in kids["children"]), "sim_seed": kids["sim_seed"],
                **{k: op[k] for k in OP_KEYS}, "L_d": sc["L_d"], "L_R": sc["L_R"], "S": sc["S"], "burnin": sc["burnin"], **G["prov"]}
        if kind == "rim":
            p = attribute_frequency(attr)
            y, _ex = one_binary_record(op, p, sc, W, kids)
            rows.append({**base, "pooled_frequency": p, "F19": record_f_binary(y), "F18": record_f_binary(y[keep18]),
                         "sites_with_any": int((y[:, 0] > 0).sum())})
        else:
            pooled = G["pooled"] if kind == "decoration_current" else G["pooled_pfg"]
            totals = G["sherd_totals"] if kind == "decoration_current" else G["pfg_totals"]
            y, _ex = one_reference_record(op, pooled, totals, sc, W, kids)
            f19 = record_f_common(y); f18 = record_f_common(y[keep18])
            rows.append({**base, **{f"F19_{c}": float(f19[c]) for c in range(COMMON)}, **{f"F18_{c}": float(f18[c]) for c in range(COMMON)}})
    return rows


def run_cells(points: list[str], kinds: list[str], records: int) -> None:
    setup(records); pts = json.loads(POINTS.read_text())
    for point in points:
        op = {k: pts[point][k] for k in OP_KEYS}
        for kind in kinds:
            cells = [(point, kind, a, w, op) for a in ATTRIBUTES for w in WINDOWS] if kind == "rim" else [(point, kind, None, w, op) for w in WINDOWS]
            t0 = time.time()
            with mp.get_context("fork").Pool(min(28, mp.cpu_count(), len(cells))) as pool:
                rows = [r for rr in pool.map(run_cell, cells) for r in rr]
            df = pd.DataFrame(rows); OUTD.mkdir(parents=True, exist_ok=True)
            out = OUTD / f"{output_stem(kind, point)}.csv"; df.to_csv(out, index=False)
            print(f"wrote {out} ({len(df)} rows, {time.time() - t0:.0f} s)")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--points", default="primary")
    ap.add_argument("--kinds", default=",".join(KINDS))
    ap.add_argument("--records", type=int, default=50)
    ap.add_argument("--readings", action="store_true", help="write output/findings/rim_r1.md from the existing cell CSVs")
    a = ap.parse_args(argv)
    if a.readings:
        from rim_r1_readings import write_findings
        setup(a.records); write_findings(G); return 0
    run_cells([p for p in a.points.split(",") if p], [k for k in a.kinds.split(",") if k], a.records)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

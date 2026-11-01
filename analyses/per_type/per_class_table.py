"""Per-class table of neutral and rival generators against the record (D30, spec 6.5; regional-drift spec sections 3 to 5).

Generalized from analyses/probe_latent_groups/per_class_resemblance.py, which
remains the committed procedure of its own findings file. For one generator
and each of its settings: records through the probe's span sampling at the
D22 operating points; per record the site-weighted moment F per class
(binomial sampling variance removed), sites present, site-level F_ST
(sherd-weighted and site-mean pool), two-cluster F_ST, richness, and the
per-class Spearman of share with collection size. The findings table gives
the record's row first, then per setting the median and 10th to 90th
percentile across records. No selection is made.

Generators. Step 0 (D30), each point at its own measured time scales:
neutral (the D22 simulator as the probe runs it); C1 occupation spread (span
ends over a window W, in units of L_R; setting "Ld" is the T1 case); C2 (C1
with the record's totals assigned to spans in rank order of length); R1
misclassification (rate); R2 vessel clumps (mean clump size). Seeds for
these: generator g, setting index j, candidate c use
SeedSequence(2026101940 + 10000 g + 100 j + c).

Regional drift (D35, regional-drift spec draft 3), time scales FROZEN per
operating point at the step 0 values (frozen_scales): C3, the innovation
target redrawn from Dirichlet(kappa target) after every post-burn-in learner
draw, swept by F_W, the regional F accumulated over W = 8 L_R (setting: F_W
on the registered grid); M, a C3 record with one class's odds multiplied by
m in one spatial cluster at the deposit stage (setting: "F_W/m", with
--marker-class and --marker-cluster small|large). Record key (point, F_W
index a, record r): root SeedSequence(2026101950 + 100000 point + 1000 a + r),
children history (0), lengths (1), ends (2), sherds (3); a marker record
uses its neutral record's children exactly, and each site draws its sherds
from its own stream (the sherds child spawned per site), so a marked site
changes no other site's counts. One CSV and findings file per operating
point; C3 and M runs take exactly one --points value.

Exploratory design validation (rule 1), not evidence about the record.
Usage: .venv/bin/python analyses/per_type/per_class_table.py --generator C1 --settings Ld,1,2,4,8 --records 50 --points primary,secondary
       .venv/bin/python analyses/per_type/per_class_table.py --generator C3 --settings 0,0.01,0.03,0.1,0.3 --records 50 --points primary
       .venv/bin/python analyses/per_type/per_class_table.py --generator M --settings 0/1.5,0/2,0/3 --marker-class Barton/Kent/MPI --marker-cluster small --records 50 --points primary
Outputs: output/probe_latent_groups/per_type_table_<generator>.csv (step 0 generators),
  output/probe_latent_groups/per_type_table_C3_<point>.csv, per_type_table_M_<point>_<class>_<cluster>.csv,
  output/findings/per_type_generators_<generator>.md, ..._C3_<point>.md, ..._M_<point>_<class>_<cluster>.md
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import math
import multiprocessing as mp
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "analyses")); sys.path.insert(0, str(ROOT / "analyses" / "probe_latent_groups"))
sys.path.insert(0, str(ROOT / "analyses" / "per_type"))

import readings as rd
import resemblance_screen as screen

from mls_emergence.dataio.basin import geo_km, load_curated
from mls_emergence.probe import generators as gen
from mls_emergence.probe import resemblance as rs
from mls_emergence.probe import rivals, spans
from mls_emergence.signatures.variance import moment_f
from mls_emergence.transmission.spatial import drift_record, kappa_for_window_f

OUTD = ROOT / "output" / "probe_latent_groups"
FINDD = ROOT / "output" / "findings"
# Artifact families this script writes, named from its arguments (rule 15; read by scripts/build_manifest.py).
MANIFEST_OUTPUTS = ["output/probe_latent_groups/per_type_table_{generator}.csv",
                    "output/findings/per_type_generators_{generator}.md",
                    "output/probe_latent_groups/per_type_table_C3_{point}.csv",
                    "output/findings/per_type_generators_C3_{point}.md",
                    "output/probe_latent_groups/per_type_table_M_{point}_{class}_{cluster}.csv",
                    "output/findings/per_type_generators_M_{point}_{class}_{cluster}.md"]
POINTS = ROOT / "analyses" / "probe_latent_groups" / "operating_points.json"
TS = ROOT / "output" / "probe_latent_groups" / "time_scales_primary.json"
SEED_BASE = 2026101940
GENERATORS = ("neutral", "C1", "C2", "R1", "R2", "C3", "M")
DRIFT_GENERATORS = ("C3", "M")
SEED_BASE_C3 = 2026101950
POINT_INDEX = {"primary": 0, "secondary": 1}
TS_FILES = {p: OUTD / f"time_scales_{p}.json" for p in POINT_INDEX}
AMPLITUDE_GRID = (0.0, 0.01, 0.03, 0.1, 0.3)
OP_KEYS = ("kernel", "length", "n_ind", "innovation", "mixing")
COMMON = 3   # the first three classes are read by F; the rest by occupancy
G: dict = {}


def _diversity():
    return importlib.import_module("47_revision_analysis").diversity


def fmt_f(x) -> str:
    return "absent" if not np.isfinite(x) else f"{x:.3f}"


def fmt_stat(x, status: str) -> str:
    """A finite value to three decimals; otherwise its status word (absent, fixed, undefined)."""
    return f"{x:.3f}" if np.isfinite(x) else (status if status in rd.STATUSES else "undefined")


def per_class_rows(y: np.ndarray, labels: np.ndarray, diversity=None) -> dict:
    """Per-class F, occupancy and status, both site-level F_ST weightings, two-cluster F_ST, richness, share-size Spearman per class.

    A record whose pooled sample holds a single class is kept: its context
    summaries are NaN with reason "undefined_pool" (the single-class case of
    cultural_fst, matched by message; any other error propagates).
    """
    diversity = diversity or _diversity()
    y = np.asarray(y, float); K = y.shape[1]; n = y.sum(1); P = y / n[:, None]
    f = moment_f(y); occ = (y > 0).sum(0); d = diversity(y); reason = ""
    try:
        s = rs.summaries(y, labels, diversity)
    except ValueError as e:
        if "single class" not in str(e):
            raise
        s = {"fst_clusters": float("nan"), "rich": d["rich"], "hs": d["hs"], "ht": d["ht"]}; reason = "undefined_pool"
    pb = P.mean(0); ht_site = 1 - (pb ** 2).sum()
    site_fst = float(1 - d["hs"] / d["ht"]) if d["ht"] > 0 else float("nan")
    site_fst_sitepool = float(1 - d["hs"] / ht_site) if ht_site > 0 else float("nan")
    with np.errstate(invalid="ignore"):
        rho = [spearmanr(P[:, c], n)[0] if occ[c] > 1 else float("nan") for c in range(K)]
    status = ["absent" if occ[c] == 0 else "fixed" if np.all(P[:, c] == 1) else "finite" if np.isfinite(f[c]) else "undefined"
              for c in range(K)]
    return {**{f"F_{c}": f[c] for c in range(K)}, **{f"occ_{c}": int(occ[c]) for c in range(K)}, **{f"rho_{c}": rho[c] for c in range(K)},
            **{f"status_{c}": status[c] for c in range(K)},
            "site_fst": site_fst, "site_fst_sitepool": site_fst_sitepool,
            "fst_clusters": s["fst_clusters"], "rich": s["rich"], "hs": s["hs"], "ht": s["ht"], "reason": reason}


def frozen_scales(point: str) -> dict:
    """Time scales frozen at the step 0 measurement (spec 3.3): W = 8 L_R; S = W + L_d rounded up to a multiple of 100 (or 10 below 1,000)."""
    d = json.loads(TS_FILES[point].read_text())
    L_d, L_R, burnin = int(d["L_d"]), int(d["L_R"]), int(d["burnin_used"])
    if "doublings" not in d:
        raise ValueError(f"{TS_FILES[point]} does not record the measurement's doublings, so its acceptance is unknown")
    W = 8 * L_R; unit = 100 if W + L_d >= 1000 else 10
    S = int(math.ceil((W + L_d) / unit) * unit)
    return {"L_d": L_d, "L_R": L_R, "W": W, "S": S, "burnin": burnin, "accepted": int(d["doublings"]) <= 3}


def record_streams(point: str, a: int, r: int) -> dict:
    """The record key's root and named child streams (spec 4.3): history 0, lengths 1, ends 2, sherds 3 (spawned into one stream per site)."""
    if point not in POINT_INDEX or a < 0 or a >= 100 or r < 0 or r >= 1000:
        raise ValueError("record key out of range")
    root = SEED_BASE_C3 + 100000 * POINT_INDEX[point] + 1000 * int(a) + int(r)
    kids = np.random.SeedSequence(root).spawn(4)
    return {"root": root, "children": kids, "sim_seed": int(np.random.default_rng(kids[0]).integers(0, 2**31 - 1))}


def generator_spec(kind: str, setting: str, L_d: int, L_R: int) -> dict:
    """The switch per_class_table applies to one record; settings are strings from the CLI."""
    if kind not in GENERATORS:
        raise ValueError(f"unknown generator {kind!r}")
    g = {"kind": kind, "setting": setting, "W": None, "rate": None, "mean_clump": None, "totals_by_length": kind == "C2",
         "F_W": None, "a": None, "kappa": None, "m": None}
    if kind == "neutral" and setting != "neutral":
        raise ValueError("the neutral generator takes no setting")
    if kind in ("C1", "C2"):
        g["W"] = int(L_d) if setting == "Ld" else int(round(float(setting) * L_R))
    elif kind == "R1":
        g["rate"] = float(setting)
    elif kind == "R2":
        g["mean_clump"] = float(setting)
    elif kind in DRIFT_GENERATORS:
        g["W"] = 8 * int(L_R)
        if kind == "M":
            parts = setting.split("/")
            if len(parts) != 2 or not parts[0] or not parts[1]:
                raise ValueError("an M setting is written F_W/m")
            fw_s, m = parts[0], float(parts[1])
            if not np.isfinite(m) or m <= 0:
                raise ValueError("m must be finite and positive")
            g["m"] = m
        else:
            fw_s = setting
        F = float(fw_s)
        hits = [i for i, x in enumerate(AMPLITUDE_GRID) if np.isclose(F, x, rtol=0, atol=1e-12)]
        if not hits:
            raise ValueError(f"F_W {F} is not on the registered amplitude grid {AMPLITUDE_GRID}")
        g["F_W"] = float(AMPLITUDE_GRID[hits[0]]); g["a"] = hits[0]; g["kappa"] = kappa_for_window_f(g["F_W"], g["W"])
    return g


def one_record(w, op, g: dict, sc: dict, kids) -> np.ndarray:
    """One record under a step 0 generator g at operating point op with time scales sc; kids are five seed children."""
    L_d, S, burnin = sc["L_d"], sc["S"], sc["burnin"]
    rec = drift_record(w, k=G["K"], n_ind=op["n_ind"], innovation=op["innovation"], mixing=op["mixing"], steps=S, n_records=S,
                       burnin=burnin, seed=int(np.random.default_rng(kids[0]).integers(0, 2**31 - 1)), initial="uniform", target=G["pooled"])
    rng_len, rng_end, rng_sub, rng_riv = (np.random.default_rng(k) for k in kids[1:5])
    lengths = spans.draw_lengths(rng_len, G["N"], L_d)
    W = g["W"] if g["W"] is not None else L_d
    if W + int(lengths.max()) > S:
        raise ValueError(f"window {W} plus the longest span exceeds the stored length {S}")
    ends = spans.ends_window(rng_end, G["N"], S, W) if g["kind"] in ("C1", "C2") else spans.ends_t1(rng_end, G["N"], S, L_d)
    totals = rivals.totals_by_span_length(G["totals"], lengths) if g["totals_by_length"] else G["totals"]
    if g["kind"] == "R2":
        prof = spans.span_profile(rec, ends, lengths)
        y = rivals.clumped_counts(prof, totals, g["mean_clump"], rng_riv)
    else:
        y = spans.sample_spans(rec, ends, lengths, totals, rng_sub)
    if g["kind"] == "R1":
        y = rivals.misclassify(y, g["rate"], G["pooled"], rng_riv)
    return y


def one_record_c3m(op: dict, g: dict, sc: dict, kids: dict, marker: dict | None):
    """One C3 (marker None) or M record (spec 3, 4): returns counts and the extras (spans, targets at ends, profile F, loss flags).

    Both paths draw each site's sherds from that site's own stream (the
    sherds child of the record key spawned into one generator per site, in
    site order), from the span profile (planted for M). A planted site
    therefore cannot change any other site's draw: an M record equals its
    C3 record bit for bit at every unmarked site, and at m = 1 everywhere
    (tests in tests/per_type). Before 2026-10-19 both paths shared one stream
    in site order, which redrew the sites after a marked one.
    """
    W, S, L_d, burnin = g["W"], sc["S"], sc["L_d"], sc["burnin"]
    if W + L_d > S:
        raise ValueError(f"window {W} plus the longest span {L_d} exceeds the stored length {S}")
    w = gen.family_weights(op["kernel"], op["length"], G["D"])
    rec, targets = drift_record(w, k=G["K"], n_ind=op["n_ind"], innovation=op["innovation"], mixing=op["mixing"], steps=S, n_records=S,
                                burnin=burnin, seed=kids["sim_seed"], initial="uniform", target=G["pooled"], target_kappa=g["kappa"],
                                return_targets=True)
    rng_len, rng_end = (np.random.default_rng(k) for k in kids["children"][1:3])
    # Site i's stream is the sherds child's i-th spawn, built directly (spawn_key (3, i)) so that reusing a key never
    # advances a spawn counter: the same key gives the same per-site streams on every call.
    site_rngs = [np.random.default_rng(np.random.SeedSequence(kids["root"], spawn_key=(3, i))) for i in range(G["N"])]
    lengths = spans.draw_lengths(rng_len, G["N"], L_d)
    ends = spans.ends_window(rng_end, G["N"], S, W)
    prof = spans.span_profile(rec, ends, lengths)
    ex = {"ends": ends, "lengths": lengths, "targets_at_ends": targets[ends], "profile_f_target": rs.profile_f(targets[ends]),
          "profile_f_before": rs.profile_f(prof), "lost": targets[S - 1] == 0.0,
          "total_f_S": float(1 - (g["kappa"] / (g["kappa"] + 1)) ** (S - 1)) if g["kappa"] is not None else 0.0}
    drawn = prof
    if marker is not None:
        drawn = rivals.plant_marker(prof, G["labels"], marker["cluster"], marker["c"], g["m"])
        ex["profile_f_after"] = rs.profile_f(drawn)
    y = np.array([r.multinomial(int(n), q) for r, n, q in zip(site_rngs, G["totals"], drawn)])
    return y, ex


def output_stem(kind: str, point: str, marker_class: str | None = None, cluster: str | None = None) -> str:
    if kind == "C3":
        return f"per_type_table_C3_{point}"
    if kind == "M":
        if not marker_class or cluster not in ("small", "large"):
            raise ValueError("M outputs need the marker class and the cluster (small or large)")
        return f"per_type_table_M_{point}_{marker_class.replace('/', '-').replace(' ', '_')}_{cluster}"
    return f"per_type_table_{kind}"


def run_candidate(args):
    c, label, point, op, g = args
    if g["kind"] in DRIFT_GENERATORS:
        return run_candidate_c3m(c, label, point, op, g)
    ss = np.random.SeedSequence(SEED_BASE + 10000 * GENERATORS.index(g["kind"]) + 100 * g["setting_index"] + c)
    kids = ss.spawn(5 * G["records"] + 5)
    w = gen.family_weights(op["kernel"], op["length"], G["D"])
    sc = screen.cell_scales(w, op["n_ind"], op["innovation"], op["mixing"], kids[5 * G["records"]:])
    gg = dict(g); gg.update(generator_spec(g["kind"], g["setting"], sc["L_d"], sc["L_R"]))
    rows = []
    for j in range(G["records"]):
        y = one_record(w, op, gg, sc, kids[5 * j:5 * j + 5])
        rows.append({"candidate": c, "label": label, "point": point, "generator": g["kind"], "setting": g["setting"], "W": gg["W"],
                     "rate": gg["rate"], "mean_clump": gg["mean_clump"], "record": j, "L_d": sc["L_d"], "L_R": sc["L_R"], "S": sc["S"],
                     "scales_accepted": sc["accepted"], **per_class_rows(y, G["labels"], G["diversity"])})
    return rows


def run_candidate_c3m(c, label, point, op, g):
    sc = frozen_scales(point)
    gg = dict(g); gg.update(generator_spec(g["kind"], g["setting"], sc["L_d"], sc["L_R"]))
    marker = G.get("marker")
    rows = []
    for r in range(G["records"]):
        kids = record_streams(point, gg["a"], r)
        y, ex = one_record_c3m(op, gg, sc, kids, marker)
        row = {"candidate": c, "label": label, "point": point, "generator": g["kind"], "setting": g["setting"], "W": gg["W"],
               "F_W": gg["F_W"], "a": gg["a"], "kappa": gg["kappa"] if gg["kappa"] is not None else 0.0, "total_f_S": ex["total_f_S"],
               "record": r, "root": kids["root"], "child_spawn_keys": ";".join(str(k.spawn_key) for k in kids["children"]), "sim_seed": kids["sim_seed"],
               **{k: op[k] for k in OP_KEYS}, "L_d": sc["L_d"], "L_R": sc["L_R"], "S": sc["S"], "burnin": sc["burnin"],
               "scales_accepted": True, **G["prov"]}
        if marker is not None:
            row.update(marker_class=G["marker_class"], marker_c=marker["c"], cluster=G["marker_cluster"], cluster_label=marker["cluster"],
                       cluster_size=int((G["labels"] == marker["cluster"]).sum()), m=gg["m"])
        row.update(per_class_rows(y, G["labels"], G["diversity"]))
        row.update({f"lost_{k}": bool(ex["lost"][k]) for k in range(G["K"])})
        row.update({f"pfT_{k}": float(ex["profile_f_target"][k]) for k in range(G["K"])})
        row.update({f"pfB_{k}": float(ex["profile_f_before"][k]) for k in range(G["K"])})
        if marker is not None:
            row.update({f"pfA_{k}": float(ex["profile_f_after"][k]) for k in range(G["K"])})
        rows.append(row)
    return rows


def q(x, p):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return float(np.quantile(x, p)) if len(x) else float("nan")


def summarize(per_seed: pd.DataFrame, real: dict, cls: list[str], cands: list) -> list[str]:
    K = len(cls)
    hdr = (["candidate", "point", "L_d, S", "site F_ST", "site F_ST (site pool)", "2-cluster F_ST", "rich", "rho PP", "rho Painted"]
           + [f"F {cls[c]}" for c in range(COMMON)] + [f"sites {cls[c]}" for c in range(COMMON, K)] + ["covered of 10"])
    L = ["| " + " | ".join(hdr) + " |", "|" + "---|" * len(hdr)]
    rec = (["**record**", "", "", f"{real['site_fst']:.4f}", f"{real['site_fst_sitepool']:.4f}", f"{real['fst_clusters']:.4f}", f"{real['rich']:.2f}",
            f"{real['rho_0']:+.2f}", f"{real['rho_2']:+.2f}"] + [fmt_f(real[f"F_{c}"]) for c in range(COMMON)]
           + [str(real[f"occ_{c}"]) for c in range(COMMON, K)] + [""])
    L.append("| " + " | ".join(rec) + " |")
    for c, (label, point, _op, _g) in enumerate(cands):
        d = per_seed[per_seed.candidate == c]; covered = 0; cells = []
        for k in range(COMMON):
            lo, me, hi = q(d[f"F_{k}"], .1), q(d[f"F_{k}"], .5), q(d[f"F_{k}"], .9)
            hit = np.isfinite(lo) and lo <= real[f"F_{k}"] <= hi; covered += int(hit)
            cells.append(f"{fmt_f(me)} [{fmt_f(lo)}, {fmt_f(hi)}]" + (" c" if hit else ""))
        for k in range(COMMON, K):
            lo, me, hi = q(d[f"occ_{k}"], .1), q(d[f"occ_{k}"], .5), q(d[f"occ_{k}"], .9)
            hit = lo <= real[f"occ_{k}"] <= hi; covered += int(hit)
            cells.append(f"{me:.0f} [{lo:.0f}, {hi:.0f}]" + (" c" if hit else ""))
        def iv(col, d=d):
            return f"{q(d[col], .5):.4f} [{q(d[col], .1):.4f}, {q(d[col], .9):.4f}]"
        row = ([label, point, f"{int(d.L_d.iloc[0])}, {int(d.S.iloc[0])}" + ("" if d.scales_accepted.iloc[0] else " (not accepted)"),
                iv("site_fst"), iv("site_fst_sitepool"), iv("fst_clusters"), f"{q(d.rich, .5):.2f} [{q(d.rich, .1):.2f}, {q(d.rich, .9):.2f}]",
                f"{q(d.rho_0, .5):+.2f} [{q(d.rho_0, .1):+.2f}, {q(d.rho_0, .9):+.2f}]", f"{q(d.rho_2, .5):+.2f} [{q(d.rho_2, .1):+.2f}, {q(d.rho_2, .9):+.2f}]"]
               + cells + [str(covered)])
        L.append("| " + " | ".join(row) + " |")
    return L


def summarize_c3m(per_seed: pd.DataFrame, real: dict, cls: list[str], cands: list) -> list[str]:
    """Per-cell bands with eligibility and status (inside / outside / unavailable), accounting, profile F medians, context summaries."""
    K = len(cls)
    hdr = (["candidate", "F_W", "kappa", "total F over S", "site F_ST", "site F_ST (site pool)", "2-cluster F_ST", "rich", "undefined pool"]
           + [f"F {cls[c]} [band] n_finite status" for c in range(COMMON)] + [f"sites {cls[c]} [band] status" for c in range(COMMON, K)]
           + ["lost from target (count per class)"] + [f"profile F target {cls[c]}" for c in range(COMMON)]
           + [f"profile F deposit {cls[c]}" for c in range(COMMON)])
    L = ["| " + " | ".join(hdr) + " |", "|" + "---|" * len(hdr)]
    rec = (["**record**", "", "", "", f"{real['site_fst']:.4f}", f"{real['site_fst_sitepool']:.4f}", f"{real['fst_clusters']:.4f}", f"{real['rich']:.2f}", ""]
           + [fmt_f(real[f"F_{c}"]) for c in range(COMMON)] + [str(real[f"occ_{c}"]) for c in range(COMMON, K)] + ["", "", "", "", "", "", ""])
    L.append("| " + " | ".join(rec) + " |")
    for c, (label, _point, _op, _g) in enumerate(cands):
        d = per_seed[per_seed.candidate == c]; cells = []
        for k in range(COMMON):
            b = rd.band(d[f"F_{k}"]); st = rd.status(real[f"F_{k}"], b)
            cells.append(f"{fmt_f(b['median'])} [{fmt_f(b['lo'])}, {fmt_f(b['hi'])}] {b['n_finite']} {st}" + ("" if b["eligible"] else " (ineligible)"))
        for k in range(COMMON, K):
            b = rd.band(d[f"occ_{k}"]); st = rd.status(real[f"occ_{k}"], b)
            cells.append(f"{b['median']:.0f} [{b['lo']:.0f}, {b['hi']:.0f}] {st}")
        lost = ",".join(str(int(d[f"lost_{k}"].astype(bool).sum())) for k in range(K))
        def iv(col, d=d):
            return f"{q(d[col], .5):.4f} [{q(d[col], .1):.4f}, {q(d[col], .9):.4f}]"
        row = ([label, f"{d.F_W.iloc[0]:g}", f"{d.kappa.iloc[0]:.0f}" if d.kappa.iloc[0] > 0 else "none", f"{d.total_f_S.iloc[0]:.4f}",
                iv("site_fst"), iv("site_fst_sitepool"), iv("fst_clusters"), f"{q(d.rich, .5):.2f} [{q(d.rich, .1):.2f}, {q(d.rich, .9):.2f}]",
                str(int((d.reason.astype(str) == "undefined_pool").sum()))]
               + cells + [lost] + [fmt_f(q(d[f"pfT_{k}"], .5)) for k in range(COMMON)] + [fmt_f(q(d[f"pfB_{k}"], .5)) for k in range(COMMON)])
        L.append("| " + " | ".join(row) + " |")
    return L


def provenance(counts: pd.DataFrame, labels: np.ndarray) -> dict:
    rev = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()
    h = hashlib.sha256()
    for p in TS_FILES.values():
        h.update(p.read_bytes())
    scales_hash = h.hexdigest()[:16]
    h = hashlib.sha256(); h.update(np.ascontiguousarray(counts.to_numpy(int)).tobytes())
    h.update("|".join(map(str, counts.index)).encode()); h.update("|".join(map(str, counts.columns)).encode())
    h.update(np.asarray(labels, int).tobytes())
    return {"git_rev": rev, "numpy_version": np.__version__, "bit_generator": np.random.default_rng(0).bit_generator.__class__.__name__,
            "scales_hash": scales_hash, "data_hash": h.hexdigest()[:16], "initial": "uniform (every site at the pooled shares)"}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--generator", choices=GENERATORS, required=True)
    ap.add_argument("--settings", default="", help="comma-separated; C1/C2: Ld or multiples of L_R; R1: rates; R2: mean clump sizes; C3: F_W; M: F_W/m")
    ap.add_argument("--records", type=int, default=50)
    ap.add_argument("--points", default="primary,secondary")
    ap.add_argument("--marker-class", default=None, help="M only: the class planted, by its column name")
    ap.add_argument("--marker-cluster", default=None, choices=("small", "large"), help="M only: the spatial cluster planted")
    ap.add_argument("--scales", default="auto", choices=("auto", "frozen", "measured"),
                    help="auto: frozen for C3 and M, measured for the step 0 generators (the only combinations supported)")
    a = ap.parse_args(argv)
    settings = [s for s in a.settings.split(",") if s] or ["neutral"]
    drift = a.generator in DRIFT_GENERATORS
    points = a.points.split(",")
    if drift and len(points) != 1:
        raise SystemExit("C3 and M runs take exactly one --points value (one artifact per operating point)")
    if a.scales != "auto" and (a.scales == "frozen") != drift:
        raise SystemExit("frozen scales are the C3 and M convention; measured scales are the step 0 convention (spec 4.3)")
    if (a.generator == "M") != (a.marker_class is not None and a.marker_cluster is not None):
        raise SystemExit("--marker-class and --marker-cluster are required for M and allowed only for M")
    m47 = importlib.import_module("47_revision_analysis")
    from a0_predictive_summaries import spatial_clusters
    counts, coords = load_curated(); xy = coords.to_numpy(float); D = geo_km(xy); cls = list(counts.columns); y_real = counts.to_numpy(int)
    labels = spatial_clusters(xy)
    G.update(N=len(counts), K=counts.shape[1], D=D, pooled=y_real.astype(float).sum(0) / y_real.sum(), totals=y_real.sum(1),
             labels=labels, diversity=m47.diversity, records=a.records, marker=None)
    if a.generator == "M":
        if a.marker_class not in cls:
            raise SystemExit(f"unknown class {a.marker_class!r}; the columns are {cls}")
        sizes = {int(g): int((labels == g).sum()) for g in np.unique(labels)}
        if len(sizes) != 2 or len(set(sizes.values())) != 2:
            raise SystemExit(f"the partition must have two clusters of different size, got {sizes}")
        small = min(sizes, key=sizes.get); large = max(sizes, key=sizes.get)
        G.update(marker={"c": cls.index(a.marker_class), "cluster": small if a.marker_cluster == "small" else large},
                 marker_class=a.marker_class, marker_cluster=a.marker_cluster)
    if drift:
        G["prov"] = provenance(counts, labels)
    screen.G.update(G); screen.G["weights"] = {}
    G["observed"] = rs.summaries(y_real, G["labels"], m47.diversity)
    real = per_class_rows(y_real, G["labels"], m47.diversity)
    pts = json.loads(POINTS.read_text())
    cands = []
    for j, st in enumerate(settings):
        for point in points:
            op = pts[point]
            label = f"{a.generator} {st}" if a.generator != "neutral" else "D22 " + point
            cands.append((label, point, {k: op[k] for k in OP_KEYS}, {"kind": a.generator, "setting": st, "setting_index": j}))
    args = [(c, *cand) for c, cand in enumerate(cands)]
    t0 = time.time()
    with mp.get_context("fork").Pool(min(28, mp.cpu_count(), len(args))) as pool:
        rows = [r for rr in pool.map(run_candidate, args) for r in rr]
    secs = time.time() - t0
    per_seed = pd.DataFrame(rows); OUTD.mkdir(parents=True, exist_ok=True)
    stem = output_stem(a.generator, points[0], a.marker_class, a.marker_cluster)
    per_seed.to_csv(OUTD / f"{stem}.csv", index=False)
    cmd = (f"--generator {a.generator} --settings {a.settings or 'neutral'} --records {a.records} --points {a.points}"
           + (f" --marker-class {a.marker_class} --marker-cluster {a.marker_cluster}" if a.generator == "M" else ""))
    L = [f"# Per-class table: generator {a.generator} against the record" + (f" ({points[0]} point)" if drift else ""), "",
         (f"Produced by `analyses/per_type/per_class_table.py {cmd}`: {len(cands)} candidates, {a.records} records each, {secs:.0f} s. "
          "Cells: median [10th, 90th percentile] across records. F is the site-weighted method-of-moments F per class with binomial "
          "sampling variance removed; sites is the number of the 28 sites at which the class is present; site F_ST is 1 - hs/ht with the "
          "sherd-weighted pool and with the site-mean pool; rho is the Spearman correlation of a class's share with collection size "
          "(PP = Parkin_Punctated). Exploratory design validation (rule 1), not evidence about the record."), ""]
    if drift:
        sc = frozen_scales(points[0])
        L += [(f"Regional drift of the innovation target (D35; regional-drift spec draft 3). Time scales frozen at the step 0 values: "
               f"L_d {sc['L_d']}, L_R {sc['L_R']}, W = 8 L_R = {sc['W']}, S {sc['S']}, burn-in {sc['burnin']}; initial state every site at the "
               f"pooled shares. Bands use NumPy's linear quantile rule with inclusive endpoints and are eligible when at least {rd.MIN_FINITE} "
               "of the records give a finite value; status is inside, outside or unavailable (ineligible band). 'lost from target' counts "
               "records in which the class's target share is exactly zero at the last stored step (numerical target loss). Profile F is the "
               "sampling-free var/(mean(1-mean)) of the target snapshots at the sampled ends and of the span profiles before sherds are drawn. "
               + (f"Marker: {a.marker_class} planted at the deposit stage, odds times m, in the {a.marker_cluster} cluster "
                  f"(label {G['marker']['cluster']}, {int((labels == G['marker']['cluster']).sum())} sites); planted records share their "
                  "neutral record's streams, so the separation readings are made by drift_sensitivity_figure.py against the C3 cells. " if a.generator == "M" else "")
               + f"Provenance per row: git {G['prov']['git_rev'][:12]}, NumPy {G['prov']['numpy_version']} ({G['prov']['bit_generator']}), "
               f"scales hash {G['prov']['scales_hash']}, data hash {G['prov']['data_hash']}."), ""]
        L += summarize_c3m(per_seed, real, cls, cands)
    else:
        L += [("Reading rule (Beau, 2026-10-19): a common class is covered when the record's F lies inside the interval, a rare class when "
               "its site count does; the last column counts covered classes and no selection is made."), ""]
        L += summarize(per_seed, real, cls, cands)
        L += ["", "A trailing \"c\" marks a covered class."]
    L += ["", f"Per-record values: `output/probe_latent_groups/{stem}.csv`.", ""]
    out = FINDD / f"per_type_generators_{stem.replace('per_type_table_', '')}.md"; out.write_text("\n".join(L)); print("\n".join(L)); print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

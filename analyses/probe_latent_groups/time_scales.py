"""Simulator time scales L_d and L_R, burn-in, and the stationarity check.

Latent-group probe, pilot step 5 (spec 3.3, 6.5; Finding 8). Exploratory
design validation (rule 1: committed even though exploratory).

S0 at each operating point named in operating_points.json (D22: the primary
exponential 6 km point and the secondary Gaussian 12 km point; the
inherited 24 km point is measured too, for the record), 28 sites,
straight-line distance, innovation from the pooled profile, 5 seeds, first
stored length 5,000 steps after a burn-in of 5,000; the estimator and
doubling rule are in mls_emergence.probe.timescales. After L_R is measured, if 5 L_R exceeds the
burn-in used, the measurement is repeated once with burn-in 5 L_R. The
stationarity rule is then applied to S0 and to S2 (Phillips's areas at leak
0.03) at the stored length the probe will use, and the regional-mean
trajectories are plotted alongside.

Usage: .venv/bin/python analyses/probe_latent_groups/time_scales.py
Outputs per point: output/probe_latent_groups/time_scales_<point>.json (the
  primary also as time_scales.json, which the screen scripts read),
  figures/probe_latent_groups/time_scales_<point>_acf.png and _regional.png,
  output/findings/simulator_time_scales.md (all points)
"""
from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "analyses"))

from mls_emergence.dataio.basin import data_fingerprint, geo_km, load_curated
from mls_emergence.probe import generators as gen
from mls_emergence.probe import timescales as ts
from mls_emergence.transmission.spatial import drift_record

OUTD = ROOT / "output" / "probe_latent_groups"
FIG = ROOT / "figures" / "probe_latent_groups"
FIND = ROOT / "output" / "findings" / "simulator_time_scales.md"
POINTS_FILE = Path(__file__).resolve().parent / "operating_points.json"
SEEDS = [1, 2, 3, 4, 5]
S0 = 5000
BURNIN0 = 5000
STRONG_LEAK = 0.03
plt.rcParams.update({"font.family": "sans-serif", "font.size": 8})


def make_run(weights, pooled, burnin, op):
    def run(seed, S):
        return drift_record(weights, k=len(pooled), n_ind=op["n_ind"], innovation=op["innovation"], mixing=op["mixing"],
                            steps=S, n_records=S, burnin=burnin, seed=seed, initial="uniform", target=pooled)
    return run


def main() -> int:
    points = {k: v for k, v in json.loads(POINTS_FILE.read_text()).items() if not k.startswith("_")}
    counts, coords = load_curated()
    fp = data_fingerprint(counts, coords)
    D = geo_km(coords.to_numpy(float))
    pooled = counts.to_numpy(float).sum(0); pooled /= pooled.sum()
    m36 = importlib.import_module("36_canonical_phase_map")
    labels, _ = m36.assign_primary_phases([str(i) for i in counts.index], coords.to_numpy(float))
    labels = np.asarray(labels)
    OUTD.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
    sections = []; status = 0
    for name, op in points.items():
        rc, sec = measure_point(name, op, D, pooled, labels, fp)
        sections.extend(sec); status = max(status, rc)
    L = ["# Simulator time scales, burn-in and stationarity (pilot step 5; D22 operating points)", "",
         ("Produced by `analyses/probe_latent_groups/time_scales.py`; estimator in `mls_emergence.probe.timescales` "
          "(spec 3.3, with the S/10 window rule). Exploratory design validation, not evidence about the record. "
          f"Points from `operating_points.json`; seeds {SEEDS}; first stored length {S0} after burn-in {BURNIN0}.")] + sections
    FIND.parent.mkdir(parents=True, exist_ok=True); FIND.write_text("\n".join(L) + "\n")
    print("\n".join(sections)); print(f"\nwrote {FIND}")
    return status


def measure_point(name, op, D, pooled, labels, fp):
    fam = gen.family_name(op["kernel"], op["length"])
    w0 = gen.family_weights(op["kernel"], op["length"], D)
    w2 = gen.family_weights(op["kernel"], op["length"], D, labels=labels, leak=STRONG_LEAK)
    burnin = BURNIN0
    res = ts.measure_time_scales(make_run(w0, pooled, burnin, op), SEEDS, S0=S0)
    repeated = False
    if res["status"] == "ok" and 5 * res["L_R"] > burnin:
        burnin = 5 * res["L_R"]; repeated = True
        res = ts.measure_time_scales(make_run(w0, pooled, burnin, op), SEEDS, S0=S0)
    if res["status"] != "ok":
        msg = f"STOP ({name}, {fam}): time-scale measurement status {res['status']} at S = {res['S']} after {res['doublings']} doublings"
        print(msg)
        (OUTD / f"time_scales_{name}.json").write_text(json.dumps({k: v for k, v in res.items() if not k.startswith("curves")}, indent=1))
        return 2, ["", f"## {name}: {fam}, {op}", "", msg]
    L_d, L_R = res["L_d"], res["L_R"]
    S_probe = int(max(2 * L_d + L_R, 10 * L_R))
    burnin_probe = int(max(5 * L_R, burnin))

    # Stationarity at the probe's stored length, S0 and S2.
    stat = {}
    traj = {}
    for sc, w in (("S0", w0), ("S2", w2)):
        runs = [drift_record(w, k=len(pooled), n_ind=op["n_ind"], innovation=op["innovation"], mixing=op["mixing"],
                             steps=S_probe, n_records=S_probe, burnin=burnin_probe, seed=100 + s, initial="uniform",
                             target=pooled) for s in SEEDS]
        stat[sc] = ts.stationarity(runs)
        traj[sc] = np.stack([r.mean(axis=1) for r in runs])    # seeds x S x K regional means

    out = {"point": name, "family": fam, "operating_point": op, "seeds": SEEDS, "first_stored_length": S0, "burnin_first": BURNIN0,
           "burnin_used": burnin, "measurement_repeated_with_5LR_burnin": repeated,
           "stored_length_measured": res["S"], "doublings": res["doublings"],
           "L_d_per_seed": res["local"], "L_R_per_seed": res["regional"], "L_d": L_d, "L_R": L_R,
           "n_included_local": res["n_included_local"], "n_excluded_local": res["n_excluded_local"],
           "n_excluded_regional": res["n_excluded_regional"],
           "probe_stored_length": S_probe, "probe_burnin": burnin_probe, "stationarity": stat,
           "data_fingerprint": fp, "labels_for_S2": labels.tolist()}
    (OUTD / f"time_scales_{name}.json").write_text(json.dumps(out, indent=1))
    if name == "primary":
        (OUTD / "time_scales.json").write_text(json.dumps(out, indent=1))
    np.save(OUTD / f"time_scales_{name}_curves_local.npy", np.array(res["curves_local"]))
    np.save(OUTD / f"time_scales_{name}_curves_regional.npy", np.array(res["curves_regional"]))

    # Figure: the five local and five regional curves with the 1/e line.
    fig, axes = plt.subplots(1, 2, figsize=(7, 3), sharey=True)
    for ax, key, title in ((axes[0], "curves_local", "local (site, class) mean"), (axes[1], "curves_regional", "regional mean")):
        for c in res[key]:
            ax.plot(np.arange(len(c)), c, lw=0.8, alpha=0.8)
        ax.axhline(1 / np.e, color="black", lw=0.6, ls="--")
        ax.set_xlabel("lag (simulator steps)"); ax.set_xlim(0, min(len(res[key][0]) - 1, 4 * max(L_d, L_R)))
        ax.text(0.98, 0.95, title, transform=ax.transAxes, ha="right", va="top")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("autocorrelation"); axes[0].axvline(L_d, color="grey", lw=0.6); axes[1].axvline(L_R, color="grey", lw=0.6)
    fig.tight_layout(); fig.savefig(FIG / f"time_scales_{name}_acf.png", dpi=300); plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(7, 4.5), sharex=True)
    for ax, sc in zip(axes, ("S0", "S2")):
        for s_ in range(len(SEEDS)):
            for k in range(traj[sc].shape[2]):
                ax.plot(traj[sc][s_, :, k], lw=0.5, alpha=0.6, color=plt.cm.viridis(k / 9))
        ax.axvline(S_probe // 2, color="black", lw=0.6, ls="--")
        ax.text(0.98, 0.95, f"{sc}: regional mean class frequency, {len(SEEDS)} seeds", transform=ax.transAxes, ha="right", va="top")
        ax.spines[["top", "right"]].set_visible(False)
    axes[1].set_xlabel("stored step"); fig.tight_layout(); fig.savefig(FIG / f"time_scales_{name}_regional.png", dpi=300); plt.close(fig)

    L = ["", f"## {name}: {fam}, n_ind {op['n_ind']}, innovation {op['innovation']:g}, mixing {op['mixing']:g}", "",
         (f"Burn-in used {burnin} (repeated with 5 L_R: {repeated}); stored length at which every curve was accepted: {res['S']} "
          f"after {res['doublings']} doublings."), "",
         "| scale | per seed | median | series included / excluded |", "|---|---|---|---|",
         f"| L_d (local) | {res['local']} | **{L_d}** | {res['n_included_local']} / {res['n_excluded_local']} |",
         f"| L_R (regional) | {res['regional']} | **{L_R}** | classes excluded {res['n_excluded_regional']} |", "",
         f"Probe settings that follow: stored length max(2 L_d + L_R, 10 L_R) = {S_probe}; burn-in max(5 L_R, {burnin}) = {burnin_probe}.", "",
         "| scenario | first-half dispersion per seed | second-half | difference of means | range across seeds | stationary |",
         "|---|---|---|---|---|---|"]
    for sc in ("S0", "S2"):
        st = stat[sc]
        L.append(f"| {sc} | {[round(v, 5) for v in st['first_half']]} | {[round(v, 5) for v in st['second_half']]} | "
                 f"{st['difference_of_means']:.5f} | {st['range_across_seeds']:.5f} | {'yes' if st['stationary'] else 'NO: stop for discussion'} |")
    L += ["", f"Figures: `time_scales_{name}_acf.png`, `time_scales_{name}_regional.png`."]
    return (0 if all(stat[n]["stationary"] for n in stat) else 3), L



if __name__ == "__main__":
    raise SystemExit(main())

"""Calibrated finite-class drift comparison; regenerate with --force.

Canonical spatial pipeline for the September 2026 revision. Two neutral
innovation models are run, each calibrated to the observed within-assemblage
diversity before any between-group statistic is compared:

- ``pooled``: innovation events draw a class in proportion to the regional
  pooled class profile, representing uneven aggregation of variants into the
  measured archaeological classes (main baseline).
- ``uniform``: the standard K-allele model, innovation uniform over classes
  (alternative null).

Calibration evaluates a grid of node population size, innovation rate, and
mixing rate against three within-assemblage summaries: mean Gini-Simpson
diversity, mean class richness at the observed sherd totals, and pooled
Gini-Simpson diversity. A cell is diversity-matched when all three fall within
``tolerance`` (relative) of the observed values. The baseline cell is the
diversity-matched cell whose median partition F_ST is highest, so that the
drift comparison is the one most favorable to drift; if no cell matches, the
minimum-loss cell is used and flagged. Between-group statistics are recorded
for every cell but never enter the matching criterion. Outputs include complete replicate tables and an
input/code/config fingerprint; a cache is reused only when it matches.
"""
from pathlib import Path
import argparse
import ast
from concurrent.futures import ProcessPoolExecutor
import hashlib
import importlib
from importlib.metadata import version
import json
import os
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analyses"))
sys.path.insert(0, str(ROOT / "src"))
import matplotlib
matplotlib.use("Agg")
import make_figures as mf
import make_map as mm
from mls_emergence.transmission.spatial import copying_weights, drift_record, sample_record
from mls_emergence.signatures.variance import cultural_fst
from mls_emergence.signatures.assortativity import similarity_matrix, mantel

sd = importlib.import_module("23_phases_vs_spatial_drift")
basin = importlib.import_module("17_basin_results")
full = importlib.import_module("35_basin_pullout")
phases = importlib.import_module("36_canonical_phase_map")
regions = importlib.import_module("31_within_region_structure")
OUT = ROOT / "output" / "revision_2026_09"
MODELS = ["pooled", "uniform"]
CONFIG = dict(k="observed class count", steps=1800, n_records=90,
              burnin=1200, window=8, models=MODELS, tolerance=0.10,
              calib_n_ind=[120, 2000, 10000],
              calib_innovations=[.0002, .0005, .001, .002, .004, .008, .012, .024, .048, .1, .2],
              calib_mixings=[.001, .002, .005, .01, .02, .05, .1, .2, .4], calib_reps=6, length=24.,
              baseline_reps=500, grid_reps=30, sensitivity_reps=50,
              lengths=[12., 24.], leaks=[1., .5, .1, .03])
REGION_LABEL = {"basin": "St. Francis basin", "valley": "Parkin partition", "cmv": "SE Missouri"}
MODEL_LABEL = {"pooled": "pooled-profile innovation", "uniform": "uniform innovation"}


def fst_by(m, labels):
    return cultural_fst(np.array([m[labels == g].sum(0) for g in np.unique(labels)]))


def diversity(m):
    """Within-assemblage Gini-Simpson, class richness, and pooled Gini-Simpson."""
    m = np.asarray(m, float)
    p = m / m.sum(1, keepdims=True)
    pooled = m.sum(0) / m.sum()
    return dict(hs=float((1 - (p ** 2).sum(1)).mean()),
                rich=float((m > 0).sum(1).mean()),
                ht=float(1 - (pooled ** 2).sum()))


def load_sets():
    b, bc = mf._load_curated()
    v, vc = full.load_full()
    c, xy = regions.load_cmv_miss()
    sets = {}
    for name, m, coords, names in [
        ("basin", b.to_numpy(int), bc.to_numpy(float), list(b.index)),
        ("valley", v.to_numpy(int), vc.to_numpy(float), list(v.index)),
        ("cmv", c.astype(int), xy, [str(i) for i in range(len(c))]),
    ]:
        coords = np.asarray(coords)
        d = sd.geo_km(coords) if name == "cmv" else mm.river_distance_matrix(coords)[0]
        ca = mf.correspondence_axis(m)[0]
        if name == "basin":
            ca = basin.oriented_ca(b)[0].reindex(names).to_numpy()
        order = np.argsort(ca)
        ranks = np.empty(len(m)); ranks[order] = np.linspace(0, 1, len(m))
        centered = coords - coords.mean(0)
        labs = {k: mf._kmeans_labels(centered, k, seed=7) for k in range(2, 5)}
        k = max(labs, key=lambda k: mf.silhouette_mean(centered, labs[k]))
        parkin = phases.assign_phases(names, coords) == "Parkin" if name == "valley" else None
        sets[name] = dict(m=m, coords=coords, names=names, d=d, ranks=ranks,
                          labels=labs[k], parkin=parkin, pooled=m.sum(0) / m.sum())
    return sets


def fingerprint(sets):
    h = hashlib.sha256(json.dumps(CONFIG, sort_keys=True).encode())
    h.update(json.dumps({name: version(name) for name in ["numpy", "scipy", "pandas", "networkx"]}, sort_keys=True).encode())
    # Hash the numerical pipeline, not unrelated numbered experiments or plotting.
    # Loaded arrays below bind source data, geographic construction, and assignments.
    for path in [ROOT / "src/mls_emergence/transmission/spatial.py",
                 ROOT / "src/mls_emergence/signatures/variance.py",
                 ROOT / "src/mls_emergence/signatures/assortativity.py",
                 ROOT / "analyses/23_phases_vs_spatial_drift.py"]:
        h.update(path.read_bytes())
    tree = ast.parse(Path(__file__).read_text())
    numerical = {"fst_by", "diversity", "simulate", "stats", "calibrate", "worker", "run"}  # rates dict carries n_ind
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name in numerical:
            h.update(ast.dump(node, include_attributes=False).encode())
    for name, data in sets.items():
        h.update(name.encode())
        for key in ["m", "coords", "d", "ranks", "labels", "parkin", "pooled"]:
            h.update(np.asarray(data[key]).tobytes() if data[key] is not None else b"none")
        h.update(json.dumps(data["names"]).encode())
    return h.hexdigest()


def simulate(data, seed, model, *, innovation, mixing, n_ind=120, length=24., leak=1.,
             burnin=1200, initial="uniform"):
    labels = data["parkin"] if data["parkin"] is not None else data["labels"]
    w = copying_weights(data["d"], length, labels=labels, leak=leak)
    target = data["pooled"] if model == "pooled" else None
    return drift_record(w, k=data["m"].shape[1], n_ind=n_ind, seed=seed, innovation=innovation,
                        mixing=mixing, burnin=burnin, initial=initial, target=target)


def stats(m, data):
    """Community and partition F_ST, modularity, distance decay, boundary excess, diversity."""
    labels, q, n = sd.ceramic_communities(m)
    dd_r, _ = mantel(similarity_matrix(m), data["d"], n_perm=1)
    # Between-community F_ST and boundary excess are undefined when the
    # modularity partition returns a single community (rule 5: cultural_fst
    # refuses fewer than two groups rather than returning 0.0, because 0.0 is
    # exactly the "no differentiation" value and would be indistinguishable
    # from a real null). They are recorded as NaN and every summary below
    # reports how many realizations defined them.
    undefined = n < 2
    result = dict(ncom=n, community_fst=float("nan") if undefined else fst_by(m, labels),
                  spatial_fst=fst_by(m, data["labels"]), Q=q, dd_r=dd_r,
                  be_community=float("nan") if undefined else sd.boundary_excess_labeled(m, data["d"], labels),
                  be_spatial=sd.boundary_excess_labeled(m, data["d"], data["labels"]),
                  **diversity(m))
    if data["parkin"] is not None:
        result["parkin_fst"] = fst_by(m, data["parkin"])
    return result, labels


def tail(obs, null):
    a = np.asarray(null, float); a = a[np.isfinite(a)]
    return (1 + int((a >= obs).sum())) / (len(a) + 1)


def calibrate(name, data, model):
    """Grid over population size, innovation, and mixing; pick the diversity-matched
    cell most favorable to drift (highest partition F_ST), else the minimum-loss cell."""
    obs = diversity(data["m"]); totals = data["m"].sum(1); rows = []
    part = data["parkin"] if data["parkin"] is not None else data["labels"]
    for n_ind in CONFIG["calib_n_ind"]:
        for innovation in CONFIG["calib_innovations"]:
            for mixing in CONFIG["calib_mixings"]:
                sims, fst = [], []
                for seed in range(CONFIG["calib_reps"]):
                    f = simulate(data, 60000 + seed, model, innovation=innovation, mixing=mixing, n_ind=n_ind)
                    m = sample_record(f, data["ranks"], totals, np.random.default_rng(61000 + seed))
                    sims.append(diversity(m)); fst.append(fst_by(m, part))
                mean = {key: float(np.mean([s[key] for s in sims])) for key in obs}
                rel = {key: abs(mean[key] - obs[key]) / obs[key] for key in obs}
                rows.append(dict(region=name, model=model, n_ind=n_ind, innovation=innovation, mixing=mixing,
                                 loss=sum(v ** 2 for v in rel.values()), matched=max(rel.values()) <= CONFIG["tolerance"],
                                 fst_median=float(np.median(fst)),
                                 **{f"sim_{k}": v for k, v in mean.items()}, **{f"obs_{k}": v for k, v in obs.items()}))
    frame = pd.DataFrame(rows)
    matched = frame[frame.matched]
    best = matched.loc[matched.fst_median.idxmax()] if len(matched) else frame.loc[frame.loss.idxmin()]
    frame["selected"] = (frame.n_ind == best.n_ind) & (frame.innovation == best.innovation) & (frame.mixing == best.mixing)
    rates = dict(innovation=float(best.innovation), mixing=float(best.mixing), n_ind=int(best.n_ind))
    return frame, rates, dict(matched=bool(best.matched), n_matched=int(len(matched)), n_cells=int(len(frame)))


def worker(args):
    name, data, model, config = args
    CONFIG.update(config)  # spawned workers re-import this module with defaults
    calibration, rates, match = calibrate(name, data, model)
    print(f"calibrated {name}/{model}: {rates} {match}", flush=True)
    totals = data["m"].sum(1)
    rows, influence, co_members = [], [], []
    for seed in range(CONFIG["baseline_reps"]):
        f = simulate(data, seed, model, **rates)
        rng = np.random.default_rng(100000 + seed)
        for kind, ranks in [("contemporaneous", np.ones(len(totals))),
                            ("time_transgressive", data["ranks"]),
                            ("reversed", 1 - data["ranks"])]:
            m = sample_record(f, ranks, totals, rng)
            st, labels = stats(m, data)
            rows.append(dict(region=name, model=model, seed=seed, sampling=kind, **st))
            if kind == "time_transgressive" and "Parkin" in data["names"]:
                pk = data["names"].index("Parkin")
                co_members.append(labels == labels[pk])
            if name == "valley" and kind == "time_transgressive":
                for j in np.where(data["parkin"])[0]:
                    keep = np.arange(len(m)) != j
                    influence.append(dict(model=model, seed=seed, site=data["names"][j],
                                          fst=fst_by(m[keep], data["parkin"][keep])))
        if (seed + 1) % 100 == 0:
            print(f"baseline {name}/{model}: {seed+1}/{CONFIG['baseline_reps']}", flush=True)
    co = None
    if co_members:
        co = pd.DataFrame(dict(region=name, model=model, site=data["names"], probability=np.mean(co_members, axis=0)))

    grid = []
    for length in CONFIG["lengths"]:
        for leak in CONFIG["leaks"]:
            for seed in range(CONFIG["grid_reps"]):
                f = simulate(data, 20000 + seed, model, length=length, leak=leak, **rates)
                m = sample_record(f, data["ranks"], totals, np.random.default_rng(30000 + seed))
                st, _ = stats(m, data)
                grid.append(dict(region=name, model=model, length=length, leak=leak, seed=seed, **st))
        print(f"grid {name}/{model}: length={length}", flush=True)

    sensitivity = []
    if name in ("basin", "valley"):
        cases = [("longer_burnin", dict(burnin=2400)),
                 ("monomorphic_burnin", dict(initial="monomorphic")),
                 ("monomorphic_no_burnin", dict(initial="monomorphic", burnin=0)),
                 ("innovation_halved", dict(innovation=rates["innovation"] / 2)),
                 ("innovation_doubled", dict(innovation=rates["innovation"] * 2)),
                 ("mixing_halved", dict(mixing=rates["mixing"] / 2)),
                 ("mixing_doubled", dict(mixing=rates["mixing"] * 2))]
        for case, kw in cases:
            for seed in range(CONFIG["sensitivity_reps"]):
                f = simulate(data, 40000 + seed, model, **{**rates, **kw})
                m = sample_record(f, data["ranks"], totals, np.random.default_rng(50000 + seed))
                st, _ = stats(m, data)
                sensitivity.append(dict(region=name, model=model, case=case, seed=seed, **st))
        print(f"sensitivity {name}/{model} done", flush=True)
    return dict(calibration=calibration, rates=dict(region=name, model=model, **rates, **match),
                baseline=pd.DataFrame(rows), influence=pd.DataFrame(influence),
                co_membership=co, grid=pd.DataFrame(grid), sensitivity=pd.DataFrame(sensitivity))


def run(sets):
    tasks = [(name, data, model, dict(CONFIG)) for name, data in sets.items() for model in MODELS]
    workers = min(len(tasks), os.cpu_count() or 1)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(worker, tasks))
    pd.concat([r["calibration"] for r in results]).to_csv(OUT / "calibration.csv", index=False)
    pd.DataFrame([r["rates"] for r in results]).to_csv(OUT / "calibrated_rates.csv", index=False)
    pd.concat([r["baseline"] for r in results]).to_csv(OUT / "baseline.csv", index=False)
    pd.concat([r["influence"] for r in results if len(r["influence"])]).to_csv(OUT / "leave_one_out.csv", index=False)
    pd.concat([r["co_membership"] for r in results if r["co_membership"] is not None]).to_csv(OUT / "co_membership.csv", index=False)
    pd.concat([r["grid"] for r in results]).to_csv(OUT / "boundary_grid.csv", index=False)
    pd.concat([r["sensitivity"] for r in results if len(r["sensitivity"])]).to_csv(OUT / "initialization_sensitivity.csv", index=False)


def _pct(a):
    """2.5, 50, 97.5 percentiles over the DEFINED (finite) values, with their count."""
    a = np.asarray(a, float); a = a[np.isfinite(a)]
    if len(a) == 0:
        return dict(lo=float("nan"), median=float("nan"), hi=float("nan"), n_defined=0)
    lo, med, hi = np.percentile(a, [2.5, 50, 97.5])
    return dict(lo=float(lo), median=float(med), hi=float(hi), n_defined=int(len(a)))


def summarize(sets):
    base = pd.read_csv(OUT / "baseline.csv")
    grid = pd.read_csv(OUT / "boundary_grid.csv")
    loo = pd.read_csv(OUT / "leave_one_out.csv")
    sens = pd.read_csv(OUT / "initialization_sensitivity.csv")
    calib = pd.read_csv(OUT / "calibration.csv")
    rates = pd.read_csv(OUT / "calibrated_rates.csv")
    observed = {name: stats(d["m"], d)[0] for name, d in sets.items()}
    result = dict(config=CONFIG, observed=observed, calibration=[], comparisons=[],
                  boundary_grid=[], leave_one_out=[], sensitivity=[])
    for _, r in rates.iterrows():
        sel = calib[(calib.region == r.region) & (calib.model == r.model) & calib.selected].iloc[0]
        frame = base[(base.region == r.region) & (base.model == r.model) & (base.sampling == "time_transgressive")]
        result["calibration"].append(dict(region=r.region, model=r.model, innovation=r.innovation,
            mixing=r.mixing, n_ind=int(r.n_ind), matched=bool(r.matched), n_matched=int(r.n_matched),
            n_cells=int(r.n_cells), loss=float(sel.loss),
            observed={k: float(sel[f"obs_{k}"]) for k in ["hs", "rich", "ht"]},
            achieved={k: _pct(frame[k]) for k in ["hs", "rich", "ht"]}))
    for (model, name, kind), frame in base.groupby(["model", "region", "sampling"]):
        for metric in ["community_fst", "spatial_fst", "ncom"] + (["parkin_fst"] if name == "valley" else []):
            a = frame[metric].to_numpy()
            result["comparisons"].append(dict(model=model, region=name, sampling=kind, metric=metric,
                observed=observed[name][metric], **_pct(a), p_upper=tail(observed[name][metric], a),
                p_lower=tail(-observed[name][metric], -a)))
    for (model, name, length, leak), frame in grid.groupby(["model", "region", "length", "leak"]):
        metric = "parkin_fst" if name == "valley" else "spatial_fst"
        a = frame[metric]; obs = observed[name][metric]; q = _pct(a)
        result["boundary_grid"].append(dict(model=model, region=name, length=length, leak=leak,
            metric=metric, observed=obs, **q, inside=bool(q["lo"] <= obs <= q["hi"]), p_upper=tail(obs, a)))
    data = sets["valley"]
    for model in MODELS:
        for j in np.where(data["parkin"])[0]:
            keep = np.arange(len(data["m"])) != j
            obs = fst_by(data["m"][keep], data["parkin"][keep])
            a = loo.loc[(loo.model == model) & (loo.site == data["names"][j]), "fst"]
            result["leave_one_out"].append(dict(model=model, site=data["names"][j], observed=obs,
                                                **_pct(a), p_upper=tail(obs, a)))
    for (model, name, case), frame in sens.groupby(["model", "region", "case"]):
        metric = "parkin_fst" if name == "valley" else "spatial_fst"
        result["sensitivity"].append(dict(model=model, region=name, case=case, metric=metric,
            **_pct(frame[metric]), p_upper=tail(observed[name][metric], frame[metric]),
            hs=float(frame.hs.mean()), rich=float(frame.rich.mean())))
    (OUT / "summary.json").write_text(json.dumps(result, indent=2))

    lines = ["# Calibrated spatial drift comparison", "",
             "Two neutral innovation models, each calibrated to observed within-assemblage "
             "Gini-Simpson diversity, richness, and pooled diversity before comparison. "
             "Fixed observed repertoire; full eight-slice windows; 1,200-generation burn-in.", "",
             "## Calibration", "",
             "| Region | Model | N | Innovation | Mixing | Matched cells | Obs H_S / rich / H_T | Achieved (baseline median) |",
             "|---|---|---:|---:|---:|---:|---|---|"]
    for r in result["calibration"]:
        o, a = r["observed"], r["achieved"]
        lines.append(f"| {r['region']} | {r['model']} | {r['n_ind']} | {r['innovation']:.4f} | {r['mixing']:.3f} | "
                     f"{r['n_matched']}/{r['n_cells']}{'' if r['matched'] else ' (unmatched; min loss)'} | "
                     f"{o['hs']:.3f} / {o['rich']:.2f} / {o['ht']:.3f} | "
                     f"{a['hs']['median']:.3f} / {a['rich']['median']:.2f} / {a['ht']['median']:.3f} |")
    lines += ["", "## Baseline comparisons (upper-tail Monte Carlo p, add-one)", "",
              "| Model | Region | Sampling | Statistic | Observed | Median [95%] | p upper | p lower |",
              "|---|---|---|---|---:|---|---:|---:|"]
    for r in result["comparisons"]:
        lines.append(f"| {r['model']} | {r['region']} | {r['sampling']} | {r['metric']} | {r['observed']:.4f} | "
                     f"{r['median']:.4f} [{r['lo']:.4f}, {r['hi']:.4f}] | {r['p_upper']:.4f} | {r['p_lower']:.4f} |")
    lines += ["", "## Boundary grid (cells bracketing observed, of 2 lengths)", "",
              "| Model | Region | Leak | Bracketing / total |", "|---|---|---:|---:|"]
    for model in MODELS:
        for name in sets:
            for leak in CONFIG["leaks"]:
                rows = [r for r in result["boundary_grid"] if r["model"] == model and r["region"] == name and r["leak"] == leak]
                lines.append(f"| {model} | {name} | {leak:g} | {sum(r['inside'] for r in rows)} / {len(rows)} |")
    lines += ["", "## Sensitivity (50 realizations per case)", "",
              "| Model | Region | Case | Median [95%] | p upper | mean H_S | mean richness |", "|---|---|---|---|---:|---:|---:|"]
    for r in result["sensitivity"]:
        lines.append(f"| {r['model']} | {r['region']} | {r['case']} | {r['median']:.4f} [{r['lo']:.4f}, {r['hi']:.4f}] | "
                     f"{r['p_upper']:.3f} | {r['hs']:.3f} | {r['rich']:.2f} |")
    lines += ["", "## Parkin leave-one-out (observation removed from observed and simulated counts)", "",
              "| Model | Site | Observed | Median [95%] | p upper |", "|---|---|---:|---|---:|"]
    for r in result["leave_one_out"]:
        lines.append(f"| {r['model']} | {r['site']} | {r['observed']:.4f} | {r['median']:.4f} [{r['lo']:.4f}, {r['hi']:.4f}] | {r['p_upper']:.3f} |")
    (OUT / "report.md").write_text("\n".join(lines) + "\n")
    print("wrote calibrated spatial results; manuscript figures come from analysis 50", flush=True)


def main(force=False):
    OUT.mkdir(parents=True, exist_ok=True)
    sets = load_sets(); key = fingerprint(sets)
    manifest = OUT / "manifest.json"
    products = ["calibration.csv", "calibrated_rates.csv", "baseline.csv", "leave_one_out.csv",
                "co_membership.csv", "boundary_grid.csv", "initialization_sensitivity.csv"]
    valid = manifest.exists() and json.loads(manifest.read_text()).get("fingerprint") == key
    if force or not valid or not all((OUT / p).exists() for p in products):
        run(sets)
        manifest.write_text(json.dumps(dict(fingerprint=key, config=CONFIG), indent=2))
    else:
        print("validated cache fingerprint; rendering stored results", flush=True)
    summarize(sets)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--force", action="store_true")
    main(parser.parse_args().force)

"""91_turnover_sensitivity.py - does the turnover result survive the smoothing scale, and can it see a boundary?

Analysis 75 maps the local rate of compositional change (turnover) with one
Gaussian bandwidth, 15 km, over 28 assemblages, and reports that the
reconstructed phase boundaries run through slower-changing country than the
boundaries of 200 same-size alternative divisions built the same way. The
review of 2026-09-24 (item 5) asks two things of that result before it can
carry the reading "the lines fall where composition changes slowly":

A. BANDWIDTH. Is the percentile a property of the record or of 15 km? The
   comparison is rerun at bandwidths from 10 to 30 km, with the same 200
   alternatives (same seed, same construction), on the same grid.

B. RECOVERY UNDER THE SAME SAMPLING GEOMETRY. With 28 points and a smoothing
   scale near the site spacing, would a copying boundary at the phase lines
   put them on a ridge at all? Simulated records from analysis 84's grid, with
   and without copying cut to 3 percent across the phase lines (no local
   innovation, and local innovation 0.2 for the strongest cell), are drawn at
   the real sherd counts and put through the same comparison. If the boundary
   cells do not raise the phase lines' percentile, the turnover comparison
   cannot see a boundary on this geometry and should be read as description.

The boundaries are this study's RECONSTRUCTION: Voronoi territories of the
assemblages merged by their Phillips (1970) phases, exactly as
36_canonical_phase_map builds them (analysis 75's territory_geometry). They are not published lines.

Output: output/findings/turnover_sensitivity.md
Usage: .venv/bin/python analyses/91_turnover_sensitivity.py [--grid 160] [--reps 100]
"""
from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analyses"))
sys.path.insert(0, str(ROOT / "src"))
OUT_MD = ROOT / "output" / "findings" / "turnover_sensitivity.md"
BANDWIDTHS = [10.0, 12.5, 15.0, 20.0, 25.0, 30.0]
N_ALT = 200
SIM_CELLS = [(1.0, 0.0), (0.1, 0.0), (0.03, 0.0), (0.03, 0.2)]   # (copying factor, local innovation)
SIM_BANDWIDTH = 15.0


def percentile_of(props, pts_km, gx, gy, bandwidth, b_pts, alt_pts, g75):
    """Median turnover along the phase boundaries and its percentile among the alternatives'."""
    field, _ = g75.turnover_field(pts_km, props, gx / 1000.0, gy / 1000.0, bandwidth, g75.MIN_WEIGHT)
    on_b = g75.read_field_at(field, gx, gy, b_pts)
    obs = float(np.median(on_b)) if len(on_b) else float("nan")
    alt = np.array([np.median(v) for v in (g75.read_field_at(field, gx, gy, P) for P in alt_pts) if len(v)])
    return obs, float(np.median(alt)), 100.0 * float(np.mean(alt < obs)), len(alt)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", type=int, default=160)
    ap.add_argument("--reps", type=int, default=100)
    args = ap.parse_args()
    g75 = importlib.import_module("75_groupness_surface")
    t74 = importlib.import_module("74_phase_partition_test")
    mf = importlib.import_module("make_figures")
    mm = importlib.import_module("make_map")
    ph = importlib.import_module("36_canonical_phase_map")
    rev = importlib.import_module("47_revision_analysis")
    a84 = importlib.import_module("84_phases_as_groups")
    from pyproj import Transformer
    from mls_emergence.transmission.spatial import copying_weights, drift_record, sample_record

    counts, coords = mf._load_curated()
    names = [str(i) for i in counts.index]
    m = counts.to_numpy(float)
    xy = coords.to_numpy(float)
    tr = Transformer.from_crs("EPSG:4326", mm.UTM15N, always_xy=True)
    E, N = (np.asarray(v, float) for v in tr.transform(xy[:, 1], xy[:, 0]))
    pts_km = np.column_stack([E, N]) / 1000.0
    labels_ph, _ = ph.assign_primary_phases(names, xy)
    pad = 18_000.0
    ext = (E.min() - pad, E.max() + pad, N.min() - pad, N.max() + pad)
    gx, gy = np.meshgrid(np.linspace(ext[0], ext[1], args.grid), np.linspace(ext[2], ext[3], args.grid))
    cells, owner, env = g75.territory_geometry(E, N, ext)
    b_pts = g75.sample_line(g75.territory_boundaries(list(labels_ph), cells, owner, env))

    # The same 200 alternatives as analysis 75 (same generator, same seed).
    km = t74.km_xy(xy)
    plist = sorted(set(labels_ph))
    sizes = np.bincount(np.array([plist.index(l) for l in labels_ph]))
    slot = np.repeat(np.arange(len(sizes)), sizes)
    rng = np.random.default_rng(75)
    alt_pts = []
    for _ in range(N_ALT):
        lab = t74.assign_exact(km, km[rng.choice(len(km), size=len(sizes), replace=False)], slot)
        gb = g75.territory_boundaries(list(lab), cells, owner, env)
        if gb is not None:
            alt_pts.append(g75.sample_line(gb))

    L = ["# Turnover along the reconstructed phase boundaries: bandwidth and recovery", "",
         f"Produced by `analyses/91_turnover_sensitivity.py`. Basin set, {len(names)} assemblages; "
         f"{args.grid} x {args.grid} grid (analysis 75 uses 220); {len(alt_pts)} same-size alternative "
         "divisions around random centers, built as analysis 75 builds them (seed 75). The boundaries are "
         "this study's reconstruction of territories from the assemblages' Phillips (1970) phases, not his drawn lines. "
         "Percentile = share of alternatives whose median turnover along their boundaries is below the phase "
         "boundaries'; a boundary between communities would sit high.", "",
         "## A. Bandwidth", "",
         "| bandwidth (km) | phase boundaries, median turnover per km | alternatives, median | percentile of the phase boundaries |",
         "|---:|---:|---:|---:|"]
    props = m / m.sum(1, keepdims=True)
    for bw in BANDWIDTHS:
        o, a, pct, n = percentile_of(props, pts_km, gx, gy, bw, b_pts, alt_pts, g75)
        L.append(f"| {bw:g} | {o:.5f} | {a:.5f} | {pct:.1f} |")
        print(L[-1], flush=True)
    L.append("")

    data = rev.load_sets()["basin"]
    if [str(n) for n in data["names"]] != names:
        raise RuntimeError("assemblage order differs between loaders")
    phase = np.array([plist.index(l) for l in labels_ph])
    totals = data["m"].sum(1)
    rates = pd.read_csv(ROOT / "output" / "revision_2026_09" / "calibrated_rates.csv")
    row = rates[(rates.region == "basin") & (rates.model == "pooled")].iloc[0]
    L += ["## B. Recovery of a copying boundary at the phase lines", "",
          f"{args.reps} simulated records per cell from analysis 84's grid and seeds (N {int(row['n_ind'])}, "
          f"innovation {row['innovation']}, mixing {row['mixing']}, 24 km along the rivers), sampled at the real "
          f"sherd counts, bandwidth {SIM_BANDWIDTH:g} km. Entries are percentiles of the phase boundaries among "
          "the same alternatives.", "",
          "| copying factor | local innovation | percentile, 5% / 25% / 50% / 75% / 95% over records | records at or below the observed percentile |",
          "|---|---|---|---:|"]
    o_obs, _, pct_obs, _ = percentile_of(props, pts_km, gx, gy, SIM_BANDWIDTH, b_pts, alt_pts, g75)
    for leak, strength in SIM_CELLS:
        li, si = a84.LEAKS.index(leak), a84.STRENGTHS.index(strength)
        w = copying_weights(data["d"], 24.0, labels=phase, leak=leak)
        pcts = []
        for r in range(args.reps):
            seed = 84000 + 10000 * li + 1000 * si + r
            rng = np.random.default_rng(seed + 5_000_000)
            tgt = a84.group_targets_by(data["pooled"], phase, strength, rng)
            rec = drift_record(w, k=m.shape[1], n_ind=int(row["n_ind"]), seed=seed,
                               innovation=float(row["innovation"]), mixing=float(row["mixing"]),
                               burnin=1200, initial="uniform", target=tgt)
            ms = sample_record(rec, data["ranks"], totals, rng).astype(float)
            pcts.append(percentile_of(ms / ms.sum(1, keepdims=True), pts_km, gx, gy, SIM_BANDWIDTH,
                                      b_pts, alt_pts, g75)[2])
        pcts = np.array(pcts)
        q = " / ".join(f"{x:.0f}" for x in np.percentile(pcts, [5, 25, 50, 75, 95]))
        L.append(f"| {leak:g} | {strength:g} | {q} | {np.mean(pcts <= pct_obs + 1e-9):.2f} |")
        print(L[-1], flush=True)
    L += ["", f"Observed at this grid and bandwidth: median {o_obs:.5f} per km, percentile {pct_obs:.1f}.", "",
          "## Reading", "",
          "If the percentile stays low across bandwidths, the result is not an artifact of 15 km. If the "
          "boundary cells in B do not lift the phase lines' percentile above the no-boundary cell, the "
          "comparison cannot see a copying boundary on this sampling geometry, and the turnover result is "
          "description rather than evidence against a boundary.", ""]
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print(f"wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

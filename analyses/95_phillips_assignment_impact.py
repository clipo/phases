#!/usr/bin/env python3
"""What the primary phase tests return under Phillips's (1970) own phase areas.

WHY. The manuscript says the primary test concerns the phases as Phillips
formulated them. Until 2026-10-01 the pipeline assigned phases from Mainfort
(1996: Figure 1) as transcribed in 36_canonical_phase_map, and that scheme
differs from Phillips's at five of the 28 analyzed assemblages: Beck, Belle
Meade, Commerce and Hollywood are Kent for Phillips and were Walls, and Castile
Landing is Parkin for Phillips and was Kent. Mainfort (2003:177) states the
first two in his own words ("two of Phillips's Kent phase sites, Belle Meade
and Beck") and lists Hollywood and Commerce among Phillips's Kent sites
(2003:176). Author ruling, 2026-10-01: the primary analysis is the phases as
Phillips assigned them, on the survey counts and Lipo's (2001) compilation,
independent of Mainfort.

WHAT THIS DOES. It runs the phase-line tests that carry the primary argument on
the same 28 assemblages under BOTH assignments, through the functions the
pipeline itself uses, so the two columns differ in the labels and in nothing
else:

  1. the partition comparison of analysis 86 (F_ST under the scheme against
     same-size divisions around random centers and compact ones, as posterior
     probabilities; boundary excess at the lines), for the three phases and for
     Parkin against the rest;
  2. the recovery test of analysis 76 (what reproduces the scheme: the site
     map alone, map plus composition, composition alone), scored as there on
     the assemblages the scheme assigns directly.

The Phillips column must reproduce output/findings/partition_posterior.md,
since it uses that analysis's seeds; the script raises if it does not, so a
difference between the columns cannot be an artifact of how this script calls
the tests.

STATUS. Written before the switch, when it measured what switching would
change. The pipeline was switched on 2026-10-01
(`36_canonical_phase_map.assign_primary_phases`), so the Phillips column is now
the pipeline's and the 1996 column is the superseded scheme, kept as the record
of what the change did. It runs after analysis 86 and checks its Phillips
column against that analysis's output.

Usage: PYTHONPATH=src .venv/bin/python analyses/95_phillips_assignment_impact.py
"""
from __future__ import annotations

import argparse
import importlib
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))
sys.path.insert(0, str(ROOT / "src"))

OUT_MD = ROOT / "output" / "findings" / "phillips_assignment_impact.md"
PIPELINE_MD = ROOT / "output" / "findings" / "partition_posterior.md"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=2000)
    ap.add_argument("--alt", type=int, default=300)
    ap.add_argument("--seeds", type=int, default=40)
    ap.add_argument("--random", type=int, default=400)
    args = ap.parse_args()

    rev = importlib.import_module("47_revision_analysis")
    mf = importlib.import_module("make_figures")
    ph = importlib.import_module("36_canonical_phase_map")
    t74 = importlib.import_module("74_phase_partition_test")
    t76 = importlib.import_module("76_phase_recovery")
    t86 = importlib.import_module("86_partition_posterior")

    counts, coords = mf._load_curated()
    names = [str(i) for i in counts.index]
    m = counts.to_numpy(float)
    xy = coords.to_numpy(float)
    pts = t74.km_xy(xy)
    data = rev.load_sets()["basin"]
    if [str(n) for n in data["names"]] != names:
        raise RuntimeError("assemblage order differs between loaders")

    schemes = {
        "Mainfort 1996 (superseded)": ph.assign_phases_by_territory(names, xy),
        "Phillips 1970 (pipeline)": ph.assign_primary_phases(names, xy),
    }
    for key, (labs, _) in schemes.items():
        if "unassigned" in set(labs):
            raise RuntimeError(f"{key} leaves analyzed assemblages unassigned: "
                               + ", ".join(np.array(names)[np.asarray(labs) == "unassigned"]))

    G = t76.unit_scale(pts)
    p = m / m.sum(1, keepdims=True)
    C = t76.unit_scale(t76.composition_features(p, "chisq"))
    weights = [0.0, 0.1, 0.25, 0.5, 1.0, np.inf]

    res = {}
    for key, (labs, derived) in schemes.items():
        labs = np.asarray(labs)
        direct = ~np.asarray(derived, dtype=bool)
        plist = sorted(set(labs))
        pidx = np.array([plist.index(l) for l in labs])
        k = len(plist)
        r = dict(sizes={pn: int((labs == pn).sum()) for pn in plist}, n_direct=int(direct.sum()))
        r["basin"] = t86.compare(m, pidx, pts, rev, t74, args.draws, args.alt, 86000, dist=data["d"])
        parkin = np.array([1 if l == "Parkin" else 0 for l in labs])
        r["parkin"] = t86.compare(m, parkin, pts, rev, t74, args.draws, args.alt, 86200, dist=data["d"])
        rec = {}
        for w in weights:
            F = G if w == 0 else (C if not np.isfinite(w) else np.column_stack([G, np.sqrt(w) * C]))
            a = [t76.ari(pidx[direct], t76.cluster(F, k, "kmeans", sd, mf)[direct])
                 for sd in range(args.seeds)]
            rec[w] = (float(np.median(a)), float(np.min(a)), float(np.max(a)))
        r["recovery"] = rec
        rng = np.random.default_rng(76)
        slot = np.repeat(np.arange(k), np.bincount(pidx))
        rnd = [t76.ari(pidx[direct], t74.assign_exact(
            pts, pts[rng.choice(len(pts), size=k, replace=False)], slot)[direct])
            for _ in range(args.random)]
        r["random"] = float(np.median(rnd))
        res[key] = r

    # The pipeline column has to be the pipeline's own numbers.
    a = res["Phillips 1970 (pipeline)"]["basin"]
    if args.draws == 2000 and args.alt == 300:
        txt = PIPELINE_MD.read_text()
        got = re.search(r"\| F_ST under the phase scheme \| ([0-9.]+) \|", txt)
        if got is None or f"{a['phase'][0]:.4f}" != got.group(1):
            raise RuntimeError("the pipeline column does not reproduce partition_posterior.md "
                               f"({a['phase'][0]:.4f} here, {got.group(1) if got else 'no value'} there)")

    keys = list(schemes)
    A, B = (res[k] for k in keys)
    la, lb = (np.asarray(schemes[k][0]) for k in keys)
    da, db = (np.asarray(schemes[k][1], dtype=bool) for k in keys)
    km = ph.phillips_outline_distance(xy, lb)
    moved = [i for i in range(len(names)) if la[i] != lb[i]]

    def row(label, f):
        return f"| {label} | {f(A)} | {f(B)} |"

    q = lambda t: f"{t[0]:.4f} ({t[1]:.4f} to {t[2]:.4f})"
    qb = lambda t: f"{t[0]:+.1f} ({t[1]:+.1f} to {t[2]:+.1f})"
    wname = {0.0: "site map alone", np.inf: "composition alone"}
    L = ["# The primary phase tests under Phillips's (1970) phase areas", "",
         "Produced by `analyses/95_phillips_assignment_impact.py`. Same 28 assemblages, same counts, same",
         "coordinates, same seeds; the two columns differ only in the phase labels. The superseded column",
         "is Mainfort (1996: Figure 1) with unmapped assemblages given the territory they fall in, which the",
         "pipeline used until 2026-10-01; the Phillips column, the pipeline's since then, is the outline each",
         "assemblage lies in (`data/raw/phillips1970_phase_outlines.csv`, Lipo 2001: Figure 2.6), with an",
         "assemblage inside no outline given the nearest.", "",
         "## The assignments", "",
         "| scheme | " + " | ".join(sorted(set(la) | set(lb))) + " | assigned directly | by territory or nearness |",
         "|---|" + "---|" * (len(set(la) | set(lb)) + 2)]
    for key, r, d in ((keys[0], A, da), (keys[1], B, db)):
        L.append(f"| {key} | " + " | ".join(str(r["sizes"].get(pn, 0)) for pn in sorted(set(la) | set(lb)))
                 + f" | {r['n_direct']} | {int(d.sum())} |")
    L += ["",
          f"{len(moved)} assemblages change phase: "
          + "; ".join(f"{names[i]} ({la[i]} to {lb[i]})" for i in moved) + ".",
          "Assigned by nearness under Phillips: "
          + "; ".join(f"{names[i]} ({lb[i]}, {km[i]:.1f} km outside the line)" for i in np.where(db)[0]) + ".",
          "",
          f"## Do the phase lines separate the pottery better than other divisions? ({args.draws} posterior draws, "
          f"{args.alt} alternatives of each kind)", "",
          f"| three phases | {keys[0]} | {keys[1]} |", "|---|---|---|",
          row("plug-in F_ST", lambda r: f"{r['basin']['plug_in']:.4f}"),
          row("F_ST under the scheme, posterior median (95% CrI)", lambda r: q(r["basin"]["phase"])),
          row("median F_ST, same sizes around random centers", lambda r: q(r["basin"]["alt"])),
          row("median F_ST, same sizes made compact", lambda r: q(r["basin"]["opt"])),
          row("P(scheme separates better than random centers)", lambda r: f"{r['basin']['p_alt']:.2f}"),
          row("P(scheme separates better than compact)", lambda r: f"{r['basin']['p_opt']:.2f}"),
          row("boundary excess at the lines, posterior median (95% CrI)", lambda r: qb(r["basin"]["be_phase"])),
          row("P(larger boundary excess than random centers)", lambda r: f"{r['basin']['be_p_alt']:.2f}"),
          row("P(larger boundary excess than compact)", lambda r: f"{r['basin']['be_p_opt']:.2f}"),
          "",
          f"| Parkin against the rest | {keys[0]} | {keys[1]} |", "|---|---|---|",
          row("plug-in F_ST", lambda r: f"{r['parkin']['plug_in']:.4f}"),
          row("F_ST, posterior median (95% CrI)", lambda r: q(r["parkin"]["phase"])),
          row("median F_ST, same sizes around random centers", lambda r: q(r["parkin"]["alt"])),
          row("median F_ST, same sizes made compact", lambda r: q(r["parkin"]["opt"])),
          row("P(separates better than random centers)", lambda r: f"{r['parkin']['p_alt']:.2f}"),
          row("P(separates better than compact)", lambda r: f"{r['parkin']['p_opt']:.2f}"),
          row("boundary excess, posterior median (95% CrI)", lambda r: qb(r["parkin"]["be_phase"])),
          row("P(larger boundary excess than random centers)", lambda r: f"{r['parkin']['be_p_alt']:.2f}"),
          row("P(larger boundary excess than compact)", lambda r: f"{r['parkin']['be_p_opt']:.2f}"),
          "",
          f"## What recovers the scheme? (adjusted Rand index over the directly assigned assemblages, k-means, "
          f"{args.seeds} starts, chi-square composition; median, with the range over starts)", "",
          f"| weight on composition | {keys[0]} | {keys[1]} |", "|---|---|---|"]
    for w in weights:
        L.append(row(wname.get(w, f"{w:g}"), lambda r: "{:.3f} ({:.3f} to {:.3f})".format(*r["recovery"][w])))
    L += [row(f"same-size divisions around random centers (median of {args.random})",
              lambda r: f"{r['random']:.3f}"), ""]
    OUT_MD.write_text("\n".join(L))
    print("\n".join(L))
    print(f"\nwrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

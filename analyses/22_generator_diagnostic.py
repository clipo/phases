"""22_generator_diagnostic.py — is "only F_ST recovers emergence" a record property or a generator artifact?

The recovery experiment (script 21) injects emergence as between-spatial-cluster divergence,
which is exactly what cultural F_ST measures, so it could be near-tautological that an F_ST
detector recovers it. This script tests that objection directly by delivering emergence through
three different channels and asking which signature recovers each, all at the real record's
configuration and rarefied to a common count.

Generators (each parameterized by strength s in [0,1], ramped along the sequence):
  A  spatial divergence + conformity: each SPATIAL cluster develops a distinct, conformist
     repertoire. Channel = between-spatial-cluster divergence (the F_ST channel). [= script 21]
  B  single-pool conformity: ONE shared repertoire that grows more conformist over the
     sequence, with NO between-group divergence. Channel = within-group conformity (the
     neutral-departure channel). Should light neutral, not F_ST.
  C  social (non-spatial) divergence + conformity: assemblages assort into social groups that
     are scrambled relative to the spatial clusters; each social group diverges. Channel =
     between-SOCIAL-group divergence (the seriation channel). Should fragment seriation while
     leaving the spatial-cluster F_ST flat.

If each signature recovers its own channel, the four-signature criterion is vindicated in
principle: you need four because emergence can arrive through any channel and a single
signature is blind to the others. If B and C recover nothing, then at this record's resolution
the signatures genuinely fail regardless of channel, and "F_ST-led" is a property of the record.

Also reports: rarefaction-depth sensitivity (does neutral/seriation recover at NRARE 100/150?)
and a k-cluster sweep (does the spatial signature recover at k=4,5?).

Writes output/generator_diagnostic.md.

Usage: .venv/bin/python analyses/22_generator_diagnostic.py
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import make_figures as mf  # noqa: E402
m21 = importlib.import_module("21_signal_recovery")
res = importlib.import_module("17_basin_results")

OUT = ROOT / "output" / "generator_diagnostic.md"
N_BINS = m21.N_BINS
KAPPA = m21.KAPPA
SIG = m21.SIG
N_SEEDS = 60
S_SHOW = [0.0, 0.4, 0.8]


def _ramp(s, t):
    return s * t / (N_BINS - 1)


def profiles_A(k, K, s, rng, clusters, bins_arr):
    grid = m21.emergence_profiles(k, K, s, rng)            # (k, N_BINS, K) spatial divergence+conf
    return np.array([grid[clusters[i], bins_arr[i]] for i in range(len(clusters))])


def profiles_B(K, s, rng, bins_arr):
    p0 = m21.zipf_base(K)
    pool = np.empty((N_BINS, K))
    for t in range(N_BINS):
        pr = p0 ** (1.0 + KAPPA * _ramp(s, t))
        pool[t] = pr / pr.sum()
    return np.array([pool[bins_arr[i]] for i in range(len(bins_arr))])


def profiles_C(K, s, rng, bins_arr, n):
    S = 3
    social = rng.integers(0, S, n)                          # social groups scrambled vs space
    p0 = m21.zipf_base(K)
    shifts = rng.permutation(K)[:S]
    grid = np.empty((S, N_BINS, K))
    for t in range(N_BINS):
        a = _ramp(s, t)
        for g in range(S):
            base = (1 - a) * p0 + a * np.roll(p0, int(shifts[g]))
            pr = base ** (1.0 + KAPPA * a)
            grid[g, t] = pr / pr.sum()
    return np.array([grid[social[i], bins_arr[i]] for i in range(n)])


def score(profiles, clusters, bins_arr, coords_c, N_arr, nrare, rng):
    M = np.zeros((len(N_arr), profiles.shape[1]))
    for i in range(len(N_arr)):
        M[i] = rng.multinomial(int(N_arr[i]), profiles[i])
    return m21.sig_rhos(m21.rarefy(M, nrare, rng), clusters, bins_arr, coords_c)


def mean_rhos(gen, s, clusters, bins_arr, coords_c, N_arr, k, K, n, nrare=m21.NRARE):
    acc = {sg: [] for sg in SIG}
    for seed in range(N_SEEDS):
        rng = np.random.default_rng(7000 + seed)
        if gen == "A":
            prof = profiles_A(k, K, s, rng, clusters, bins_arr)
        elif gen == "B":
            prof = profiles_B(K, s, rng, bins_arr)
        else:
            prof = profiles_C(K, s, rng, bins_arr, n)
        r = score(prof, clusters, bins_arr, coords_c, N_arr, nrare, rng)
        for sg in SIG:
            acc[sg].append(r[sg])
    return {sg: float(np.nanmean(acc[sg])) for sg in SIG}


def build_config(k_override=None):
    counts, coords = mf._load_curated()
    ca, _ = res.oriented_ca(counts)
    K = counts.shape[1]
    Ni = dict(zip(counts.index, counts.to_numpy(float).sum(1)))
    cdf = coords.dropna()
    have = list(cdf.index)
    cc = cdf[["Latitude", "Longitude"]].to_numpy(float)
    cc_c = cc - cc.mean(0)
    if k_override is None:
        sil = {kk: mf.silhouette_mean(cc_c, mf._kmeans_labels(cc_c, kk, seed=7)) for kk in range(2, 7)}
        k = max(sil, key=sil.get)
    else:
        k = k_override
    cl = mf._kmeans_labels(cc_c, k, seed=7)
    cluster_of = dict(zip(have, cl))
    bins = pd.qcut(ca.reindex(have), N_BINS, labels=False, duplicates="drop")
    keep = [i for i in have if not pd.isna(bins[i])]
    clusters = np.array([cluster_of[i] for i in keep])
    bins_arr = np.array([int(bins[i]) for i in keep])
    coords_c = cc_c[[have.index(i) for i in keep]]
    N_arr = np.array([Ni[i] for i in keep])
    return clusters, bins_arr, coords_c, N_arr, k, K, len(keep)


def main():
    clusters, bins_arr, coords_c, N_arr, k, K, n = build_config()

    L = ["# Generator diagnostic: is 'only F_ST' a record property or a generator artifact?", "",
         f"Real config: n = {n}, K = {K} types, k = {k} spatial clusters, {N_BINS} bins, "
         f"rarefied to {m21.NRARE}. {N_SEEDS} seeds. Entries are mean Spearman rho of the "
         f"signature vs ordinal position (emergence should drive the channel's signature toward +1).",
         ""]

    GEN_DESC = {"A": "spatial divergence+conformity (F_ST channel)",
                "B": "single-pool conformity (neutral channel)",
                "C": "social non-spatial divergence (seriation channel)"}
    for gen in ("A", "B", "C"):
        L += [f"## Generator {gen}: {GEN_DESC[gen]}", "",
              "| s | " + " | ".join(SIG) + " |", "|---|" + "---|" * len(SIG)]
        rows = {}
        for s in S_SHOW:
            r = mean_rhos(gen, s, clusters, bins_arr, coords_c, N_arr, k, K, n)
            rows[s] = r
            L.append(f"| {s:.1f} | " + " | ".join(f"{r[sg]:+.2f}" for sg in SIG) + " |")
        # which signature responds most strongly to this generator's emergence (s=0.8 minus s=0)
        resp = {sg: rows[0.8][sg] - rows[0.0][sg] for sg in SIG}
        winner = max(resp, key=resp.get)
        L += ["", f"- Largest response (rho at s=0.8 minus s=0): **{winner}** "
              f"({resp[winner]:+.2f}). Per-signature response: " +
              ", ".join(f"{sg} {resp[sg]:+.2f}" for sg in SIG) + ".", ""]

    # rarefaction-depth sensitivity (generator A, s=0.6): do neutral/seriation recover deeper?
    L += ["## Rarefaction-depth sensitivity (generator A, s = 0.6)", "",
          "| NRARE | " + " | ".join(SIG) + " |", "|---|" + "---|" * len(SIG)]
    for nr in (50, 100, 150):
        r = mean_rhos("A", 0.6, clusters, bins_arr, coords_c, N_arr, k, K, n, nrare=nr)
        L.append(f"| {nr} | " + " | ".join(f"{r[sg]:+.2f}" for sg in SIG) + " |")
    L += ["", "Larger NRARE keeps more rare types; if neutral/seriation climb with depth, their "
          "weakness at NRARE=50 is partly a rarefaction-floor artifact, not a record property.", ""]

    # k-cluster sweep (generator A, s=0.6): does the spatial signature recover at higher k?
    L += ["## Cluster-count sweep (generator A, s = 0.6)", "",
          "| k | " + " | ".join(SIG) + " |", "|---|" + "---|" * len(SIG)]
    for kk in (3, 4, 5):
        cl2, ba2, cc2, N2, kused, K2, n2 = build_config(k_override=kk)
        r = mean_rhos("A", 0.6, cl2, ba2, cc2, N2, kused, K2, n2)
        L.append(f"| {kused} | " + " | ".join(f"{r[sg]:+.2f}" for sg in SIG) + " |")
    L += ["", "At k=3 the spatial boundary and F_ST summarize one partition; if the spatial "
          "signature climbs at k=4,5 its weakness is a low-cluster-count artifact.", ""]

    # verdict
    L += ["## Verdict", "",
          "Read the three generator tables together. If each signature shows its largest response "
          "to the generator whose channel it measures (A->F_ST, B->neutral, C->seriation), the "
          "four signatures are each empirically sufficient for a distinct emergence channel, and "
          "the four-signature criterion is vindicated in principle even though, for the spatial "
          "bounded-group emergence the Parkin hypothesis predicts, F_ST is the cleanest at this "
          "resolution. If B and C fail to light any signature, the signatures genuinely fail at "
          "this resolution regardless of channel, and the F_ST-led reading is a record property."]
    OUT.write_text("\n".join(L), encoding="utf-8")
    print(f"wrote {OUT}")
    print("\n".join(L))


if __name__ == "__main__":
    main()

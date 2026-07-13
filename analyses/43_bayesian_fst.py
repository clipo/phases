"""43_bayesian_fst.py - Bayesian credible interval on the basin cultural F_ST.

Reuses analysis 07's spatial clusters (prepare_inputs) to build the observed
St. Francis basin between-cluster counts, fits the hierarchical
Dirichlet-multinomial model, and reports a credible interval on the Gini-Simpson
F_ST alongside the frequentist plug-in. The credible interval quantifies
estimation uncertainty; the separate drift-null test (analyses 21/33/35/37) asks
whether drift alone could produce the value.

Usage: .venv/bin/python analyses/43_bayesian_fst.py [--fast]
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from figstyle import OI_BLUE, OI_VERMIL, save  # noqa: E402

OUT_MD = ROOT / "output" / "bayesian_fst.md"
OUT_NC = ROOT / "output" / "bayesian_fst.nc"

# target_accept=0.99 + a larger budget tames the low-F_ST Dirichlet concentration
# funnel (the basin F_ST ~0.03 drives large concentration, a Neal-funnel analogue):
# verified R-hat < 1.01, adequate ESS, and 0 divergences across seeds. Sampler-only;
# the posterior is unchanged.
FULL = dict(draws=3000, tune=3000, chains=4, target_accept=0.99)
FAST = dict(draws=300, tune=500, chains=2, target_accept=0.99)


def basin_group_counts(inp):
    """Observed basin between-cluster counts (n_clusters, K) and sizes."""
    ids = list(inp.have_coords_ids)
    labels = np.array([inp.cluster_of[a] for a in ids])
    counts_have = inp.counts.loc[ids].to_numpy(float)
    uniq = np.unique(labels)
    gc = np.array([counts_have[labels == c].sum(axis=0) for c in uniq])
    return gc, gc.sum(axis=1)


def main(fast=False):
    from mls_emergence.inference import sample_fst, fst_summary

    a07 = importlib.import_module("07_refined_empirical")

    cfg = FAST if fast else FULL
    inp = a07.prepare_inputs()
    gc, sizes = basin_group_counts(inp)

    idata = sample_fst(gc, random_seed=0, **cfg)
    s = fst_summary(idata, sizes=sizes, group_counts=gc)

    import arviz as az
    def _flat(diag):
        ds = diag.dataset if hasattr(diag, "dataset") else diag
        return np.concatenate([np.atleast_1d(v.values).ravel() for v in ds.data_vars.values()])
    rhat = float(np.max(_flat(az.rhat(idata))))
    ess = float(np.min(_flat(az.ess(idata))))
    ndiv = int(idata.sample_stats["diverging"].sum())

    L = ["# Bayesian credible interval on the basin cultural F_ST", "",
         f"Observed St. Francis basin, {gc.shape[0]} spatial clusters, "
         f"{gc.shape[1]} decorated types. Hierarchical Dirichlet-multinomial, "
         f"{cfg['draws']} draws x {cfg['chains']} chains "
         f"({'FAST' if fast else 'full'}).", "",
         "## Cultural F_ST (Gini-Simpson estimator)", "",
         f"- Posterior mean F_ST = {s['fst_mean']:.4f}, 95% credible interval "
         f"[{s['fst_hdi95'][0]:.4f}, {s['fst_hdi95'][1]:.4f}].",
         f"- Frequentist plug-in (variance.cultural_fst) = {s['plugin_fst']:.4f}.",
         f"- P(F_ST > 0) = {s['p_fst_gt0']:.3f}.", "",
         "The credible interval is the estimation uncertainty on the reported "
         "F_ST. Whether the value is distinguishable from stochastic drift is a "
         "separate test (analyses 21/33/35/37).", "",
         "## MCMC diagnostics", "",
         f"- max R-hat = {rhat:.4f} (want < 1.01); min ESS = {ess:.0f}; "
         f"divergences = {ndiv}.", ""]

    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    idata.to_netcdf(OUT_NC)

    fig, ax = plt.subplots(figsize=(3.5, 3))
    ax.hist(s["fst_samples"], bins=40, density=True, color=OI_BLUE, alpha=0.7)
    ax.axvline(s["plugin_fst"], color=OI_VERMIL, lw=1.4, label="frequentist plug-in")
    ax.set_xlabel("cultural $F_{ST}$")
    ax.set_ylabel("posterior density")
    ax.legend(frameon=False, fontsize=7)
    save(fig, "bayesian_fst")
    plt.close(fig)

    print(f"wrote {OUT_MD}")
    print("\n".join(L))


if __name__ == "__main__":
    main(fast="--fast" in sys.argv)

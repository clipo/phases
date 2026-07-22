"""43_bayesian_fst.py - Bayesian cultural F_ST for the basin (Balding-Nichols).

Reuses analysis 07's spatial clusters (prepare_inputs) to build the observed
St. Francis basin between-cluster counts, fits the Balding-Nichols model
(mirroring the hyperlocality project's src/bayes.py), and reports three
model-based quantities from one fit:

  1. the BN cultural F_ST parameter (F ~ Uniform(0,1), estimated directly);
  2. the Gini-Simpson F_ST the manuscript reports, reconstructed conjugately for
     the observed clusters (a credible interval on the exact reported estimator);
  3. the structure-vs-panmixia Bayes factor.

The headline F_ST is reported under a second prior (Beta(1,3)) as a sensitivity
check. The credible interval quantifies estimation uncertainty; the Bayes factor
weighs structure against panmixia. Neither replaces the separate stochastic-drift
null (analyses 21/33/35/37), which asks whether spatial drift alone could produce
the value.

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
from figstyle import OI_BLUE, OI_VERMIL, OI_GREEN, save  # noqa: E402

OUT_MD = ROOT / "output" / "bayesian_fst.md"
OUT_NC = ROOT / "output" / "bayesian_fst.nc"

# The marginalized Balding-Nichols model is well conditioned, so the default
# target_accept=0.9 suffices; no funnel-taming is needed (unlike the earlier
# non-centered Dirichlet-multinomial).
FULL = dict(draws=2000, tune=2000, chains=4)
FAST = dict(draws=300, tune=500, chains=2)
BF_FULL = dict(draws=2000, chains=4)
BF_FAST = dict(draws=500, chains=2)

ALT_PRIOR = ("beta", 1.0, 3.0)   # a sterner prior (more mass on low F_ST)


def basin_group_counts(inp):
    """Observed basin between-cluster counts (n_clusters, K) and sizes."""
    ids = list(inp.have_coords_ids)
    labels = np.array([inp.cluster_of[a] for a in ids])
    counts_have = inp.counts.loc[ids].to_numpy(float)
    uniq = np.unique(labels)
    gc = np.array([counts_have[labels == c].sum(axis=0) for c in uniq])
    return gc, gc.sum(axis=1)


def main(fast=False):
    import arviz as az
    from mls_emergence.inference import (
        sample_fst, fst_summary, gini_simpson_summary, bayes_factor_structure,
    )

    a07 = importlib.import_module("07_refined_empirical")

    cfg = FAST if fast else FULL
    bf_cfg = BF_FAST if fast else BF_FULL
    inp = a07.prepare_inputs()
    gc, _sizes = basin_group_counts(inp)

    # (1) BN cultural F_ST parameter, uniform prior
    idata = sample_fst(gc, random_seed=0, **cfg)
    s = fst_summary(idata)

    # (2) Gini-Simpson readout (the manuscript estimator), conjugate reconstruction
    g = gini_simpson_summary(idata, gc, seed=0)

    # prior sensitivity on the BN parameter
    idata_alt = sample_fst(gc, random_seed=0, f_prior=ALT_PRIOR, **cfg)
    s_alt = fst_summary(idata_alt)

    # (3) structure vs panmixia Bayes factor
    bf = bayes_factor_structure(gc, seed=0, **bf_cfg)

    def _flat(diag):
        ds = diag.dataset if hasattr(diag, "dataset") else diag
        return np.concatenate([np.atleast_1d(v.values).ravel() for v in ds.data_vars.values()])
    rhat = float(np.max(_flat(az.rhat(idata))))
    ess = float(np.min(_flat(az.ess(idata))))
    ndiv = int(idata.sample_stats["diverging"].sum())

    L = ["# Bayesian cultural F_ST for the basin (Balding-Nichols model)", "",
         f"Observed St. Francis basin, {gc.shape[0]} spatial clusters, "
         f"{gc.shape[1]} decorated types. Balding-Nichols Dirichlet-multinomial "
         f"(F ~ Uniform(0,1), per-cluster frequencies marginalized), "
         f"{cfg['draws']} draws x {cfg['chains']} chains "
         f"({'FAST' if fast else 'full'}). This mirrors the hyperlocality "
         f"project's Bayesian F_ST model.", "",
         "## 1. Cultural F_ST parameter (Balding-Nichols theta)", "",
         f"- Posterior median F_ST = {s['fst_median']:.4f} "
         f"(mean {s['fst_mean']:.4f}), 95% credible interval "
         f"[{s['fst_hdi95'][0]:.4f}, {s['fst_hdi95'][1]:.4f}].",
         "- Flat prior on F_ST itself, so the estimate is not pulled toward the "
         "drift level by the prior.",
         f"- Prior sensitivity (Beta(1,3) prior on F): median "
         f"{s_alt['fst_median']:.4f}, 95% CI "
         f"[{s_alt['fst_hdi95'][0]:.4f}, {s_alt['fst_hdi95'][1]:.4f}].", "",
         "## 2. Gini-Simpson F_ST (the manuscript estimator), conjugate readout", "",
         f"- Posterior median = {g['gst_median']:.4f} (mean {g['gst_mean']:.4f}), "
         f"95% credible interval [{g['gst_hdi95'][0]:.4f}, {g['gst_hdi95'][1]:.4f}].",
         f"- Frequentist plug-in (variance.cultural_fst) = {g['plugin_fst']:.4f}, "
         f"which falls "
         f"{'inside' if g['gst_hdi95'][0] <= g['plugin_fst'] <= g['gst_hdi95'][1] else 'outside'} "
         f"the 95% interval.", "",
         "This readout is a credible interval on the exact estimator the "
         "manuscript reports, reconstructed from the same fit as "
         "p_g | x_g ~ Dirichlet((1-F)/F * pi + x_g).", "",
         "## 3. Structure vs panmixia (Bayes factor)", "",
         f"- log marginal likelihood: structure (M1) = {bf['logml_structure']:.2f}, "
         f"panmixia (M0) = {bf['logml_panmixia']:.2f}.",
         f"- 2 ln BF10 = {bf['two_ln_bf']:.2f} "
         f"(+/- {bf['two_ln_bf_chain_sd']:.2f} across chains); log10 BF10 = "
         f"{bf['log10_bf']:.2f}.",
         f"- Kass & Raftery (1995): {bf['evidence']}.", "",
         "CAVEAT: this Bayes factor tests structure against EXACT panmixia (every "
         "cluster sharing one frequency vector). With thousands of sherds, exact "
         "panmixia is rejected trivially for any real assemblage set, so a large "
         "value here is expected and is NOT evidence of bounded groups. It is a "
         "much weaker null than drift: spatially structured neutral drift itself "
         "produces F_ST > 0 and would also reject panmixia. The paper's operative "
         "test is structure vs DRIFT, addressed by the stochastic-drift null "
         "(analyses 21/33/35/37) and the drift-versus-groups comparison, not by "
         "this panmixia Bayes factor. It is reported here for parity with the "
         "hyperlocality analysis and as a model-adequacy check, not as the "
         "structure test.", "",
         "## MCMC diagnostics (uniform-prior fit)", "",
         f"- max R-hat = {rhat:.4f} (want < 1.01); min ESS = {ess:.0f}; "
         f"divergences = {ndiv}.", ""]

    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    idata.to_netcdf(OUT_NC)

    fig, ax = plt.subplots(figsize=(3.5, 3))
    ax.hist(s["fst_samples"], bins=40, density=True, color=OI_BLUE, alpha=0.65,
            label="BN F_ST parameter")
    ax.hist(g["gst_samples"], bins=40, density=True, color=OI_GREEN, alpha=0.55,
            label="Gini-Simpson readout")
    ax.axvline(g["plugin_fst"], color=OI_VERMIL, lw=1.4, label="frequentist plug-in")
    ax.set_xlabel("cultural $F_{ST}$")
    ax.set_ylabel("posterior density")
    ax.legend(frameon=False, fontsize=7)
    save(fig, "bayesian_fst")
    plt.close(fig)

    print(f"wrote {OUT_MD}")
    print("\n".join(L))


if __name__ == "__main__":
    main(fast="--fast" in sys.argv)

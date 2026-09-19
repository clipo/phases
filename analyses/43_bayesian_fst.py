"""43_bayesian_fst.py - Bayesian cultural F_ST for the basin (Balding-Nichols).

Reuses analysis 07's spatial clusters (prepare_inputs) to build the observed
St. Francis basin between-cluster counts, fits the Balding-Nichols model
(the estimator is shared with a sibling project, ``../mataa``; method
provenance only, no Rapa Nui data is involved), and reports three
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

# PRIMARY is Beta(1,10), not Uniform(0,1). Changed 2026-09-02 on the strength of
# the prior predictive in analyses/57_fst_prior_predictive.py: at this design the
# flat prior expects a median Gini-Simpson F_ST of 0.27 and puts 47 percent of
# its mass above 0.30, against an observed 0.018. Nobody holds that belief about
# assemblages a few tens of km apart in one drainage, so on the scale that
# matters the flat prior is strongly and wrongly informative rather than
# uninformative (F2, rule 20a).
#
# The family must be Beta(1, b): any Beta(a, b) with a > 1 has zero density at
# F = 0 and so asserts that panmixia is impossible, which is disqualifying for a
# comparison against panmixia.
#
# Measured effect of the change on this fit: the BN parameter moves 0.0738 ->
# 0.0626, and the Gini-Simpson readout the manuscript reports moves 0.0179 ->
# 0.0178. The reported quantity is robust to the prior; the parameter is not.
PRIMARY_PRIOR = ("beta", 1.0, 10.0)
ALT_PRIOR = ("uniform",)          # reported alongside, per rule 20(b)


def basin_group_counts(inp, scope: str = "basin"):
    """Observed between-cluster counts (n_clusters, K) and sizes.

    ``scope="basin"`` (default) restricts to the canonical drainage-basin
    membership in ``data/processed/basin_members_curated.txt`` (29 assemblages)
    and re-selects the number of spatial clusters on the BASIN's coordinates by
    the same silhouette rule analysis 07 uses. That gives k = 3, which is the
    "three spatial clusters" the manuscript describes.

    ``scope="region"`` reproduces the pre-2026-08-31 behavior: the whole curated
    set with coordinates (55 assemblages) grouped by ``inp.cluster_of``, whose k
    is selected on the wider set and comes out at 5. It is retained because it
    is a legitimate quantity at a different spatial scale, not because it was
    the intent here.

    HISTORY. Until 2026-08-31 this function ignored the basin membership
    entirely and always returned the 55-assemblage grouping, while being named
    ``basin_group_counts`` and feeding a report headed "Observed St. Francis
    basin". Restricting properly moves the Gini-Simpson readout from 0.0327 to
    0.0179. Measured in ``analyses/47_basin_scope_check.py``; see
    ``output/findings/basin_scope_check.md`` and
    ``docs/CODE_REVIEW_2026-08-31.md`` F18. No manuscript number was affected,
    because the manuscript figures reach the data through
    ``make_figures._load_curated()``, which does apply the membership.
    """
    if scope not in ("basin", "region"):
        raise ValueError(f"scope must be 'basin' or 'region', got {scope!r}")

    ids_all = [str(a) for a in inp.have_coords_ids]
    counts_all = inp.counts.loc[list(inp.have_coords_ids)].to_numpy(float)

    if scope == "region":
        ids = ids_all
        labels = np.array([inp.cluster_of[a] for a in inp.have_coords_ids])
        counts_use = counts_all
    else:
        import make_figures as mf
        members = mf._basin_members("curated")
        keep = [i for i, a in enumerate(ids_all) if a in members]
        if not keep:
            raise RuntimeError(
                "no assemblage in prepare_inputs() is in the canonical basin "
                "membership; data/processed/basin_members_curated.txt is stale "
                "or the id conventions have diverged")
        ids = [ids_all[i] for i in keep]
        counts_use = counts_all[keep]
        coords_b = np.asarray(inp.cc_c)[keep]
        labels = _basin_cluster_labels(coords_b)

    uniq = np.unique(labels)
    gc = np.array([counts_use[labels == c].sum(axis=0) for c in uniq])
    return gc, gc.sum(axis=1)


def _basin_cluster_labels(coords, seed: int = 7, kmin: int = 2, kmax: int = 6):
    """Analysis 07's cluster rule, applied to whatever coordinates it is given.

    k is the value in [kmin, kmax] maximizing mean silhouette. Kept identical to
    ``07_refined_empirical`` so the basin grouping is produced by the same rule
    as the regional one, differing only in the point set it sees.
    """
    from mls_emergence.signatures.assortativity import _kmeans_labels
    silhouette_mean = importlib.import_module("07_refined_empirical").silhouette_mean
    sil = {k: silhouette_mean(coords, _kmeans_labels(coords, k, seed=seed))
           for k in range(kmin, kmax + 1)}
    return _kmeans_labels(coords, max(sil, key=sil.get), seed=seed)


def fitted_basin_ids(inp) -> list:
    """The assemblage ids the default (basin) fit actually uses. For tests."""
    import make_figures as mf
    members = mf._basin_members("curated")
    return sorted(str(a) for a in inp.have_coords_ids if str(a) in members)


def main(fast=False):
    import arviz as az
    from mls_emergence.inference import (
        sample_fst, fst_summary, gini_simpson_summary, bayes_factor_structure,
    )

    a07 = importlib.import_module("07_refined_empirical")

    cfg = FAST if fast else FULL
    bf_cfg = BF_FAST if fast else BF_FULL
    inp = a07.prepare_inputs()
    gc, _sizes = basin_group_counts(inp)          # scope="basin" by default
    a43_ids = fitted_basin_ids(inp)

    # (1) BN cultural F_ST parameter, uniform prior
    idata = sample_fst(gc, random_seed=0, f_prior=PRIMARY_PRIOR, **cfg)
    s = fst_summary(idata)

    # (2) Gini-Simpson readout (the manuscript estimator), conjugate reconstruction
    g = gini_simpson_summary(idata, gc, seed=0)

    # prior sensitivity on the BN parameter
    idata_alt = sample_fst(gc, random_seed=0, f_prior=ALT_PRIOR, **cfg)
    s_alt = fst_summary(idata_alt)

    # (3) The structure-vs-panmixia Bayes factor is DEMOTED (F1) and no longer
    # computed here. See the docstring of inference.bayesian_fst.
    # bayes_factor_structure remains runnable for provenance.

    def _flat(diag):
        ds = diag.dataset if hasattr(diag, "dataset") else diag
        return np.concatenate([np.atleast_1d(v.values).ravel() for v in ds.data_vars.values()])

    def _diagnostics(id_, label):
        """Full rule-16 set. Applied to EVERY fit, including sensitivity refits.

        Until 2026-08-31 only the uniform-prior fit was diagnosed and the
        Beta(1,3) refit -- the stated defense of the prior -- reported an
        interval with nothing behind it. See docs/CODE_REVIEW_2026-08-31.md F11.
        """
        rhat = float(np.max(_flat(az.rhat(id_))))
        ess_bulk = float(np.min(_flat(az.ess(id_))))
        ess_tail = float(np.min(_flat(az.ess(id_, method="tail"))))
        div = int(id_.sample_stats["diverging"].sum())
        ntot = int(id_.sample_stats["diverging"].size)
        td = id_.sample_stats.get("tree_depth")
        n_td = int((td >= 10).sum()) if td is not None else -1
        eb = id_.sample_stats.get("energy")
        if eb is not None:
            e = np.asarray(eb.values)
            de = np.diff(e, axis=1)
            ebfmi = float(np.min((de ** 2).mean(axis=1) / e.var(axis=1)))
        else:
            ebfmi = float("nan")
        return (f"- {label}: max R-hat = {rhat:.4f} (want < 1.01); "
                f"min bulk ESS = {ess_bulk:.0f}; min tail ESS = {ess_tail:.0f}; "
                f"divergences = {div}/{ntot} ({100 * div / ntot:.2f}%); "
                f"treedepth >= 10 = {n_td}; min E-BFMI = {ebfmi:.3f} "
                f"(want > 0.3).")

    diag_lines = [_diagnostics(idata, "Beta(1,10) prior (primary)"),
                  _diagnostics(idata_alt, "Uniform(0,1) prior (sensitivity refit)")]

    L = ["# Bayesian cultural F_ST for the basin (Balding-Nichols model)", "",
         f"Observed St. Francis basin (canonical drainage membership, "
         f"{len(a43_ids)} assemblages), {gc.shape[0]} spatial clusters, "
         f"{gc.shape[1]} decorated types. Balding-Nichols Dirichlet-multinomial "
         f"(F ~ Beta(1,10), per-cluster frequencies marginalized), "
         f"{cfg['draws']} draws x {cfg['chains']} chains "
         f"({'FAST' if fast else 'full'}). The estimator is shared with a "
         f"sibling project for consistency of method; no data is shared.", "",
         "## 1. Cultural F_ST parameter (Balding-Nichols theta)", "",
         f"- Posterior median F_ST = {s['fst_median']:.4f} "
         f"(mean {s['fst_mean']:.4f}), 95% credible interval "
         f"[{s['fst_hdi95'][0]:.4f}, {s['fst_hdi95'][1]:.4f}].",
         "- Beta(1,10) prior on F_ST, adopted on prior-predictive grounds "
         "(analyses/57). It leans toward SMALL F_ST, which is the direction "
         "of this paper's conclusion, so the Uniform(0,1) refit below is the "
         "check that the prior is not doing the work.",
         f"- Prior sensitivity (Uniform(0,1) prior on F): median "
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
         "## 3. Structure vs panmixia: not reported", "",
         "The Bayes factor against exact panmixia that this section used to "
         "carry is **demoted and no longer computed** (F1). Two reasons.", "",
         "First, its evidence band is prior-unstable exactly where it would "
         "matter. `../mataa` computed the same comparison exactly and measured "
         "the Kass-Raftery band moving across a boundary, 4.48 log units from "
         "Uniform(0,1) to Beta(1,10), for a near-null case. That is the "
         "Lindley-Bartlett effect on a nested boundary comparison and no "
         "estimator fixes it, so the quantity is robust where it is not needed "
         "and fragile where it would be load-bearing.", "",
         "Second, exact panmixia is not the alternative anyone entertains. With "
         "thousands of sherds it is rejected trivially by any real assemblage "
         "set, and spatially structured neutral drift, which IS the alternative "
         "this paper tests, produces F_ST > 0 and would reject it too. The "
         "operative comparison is structure against drift, addressed by the "
         "drift-versus-groups analysis, not by a panmixia Bayes factor.", "",
         "What stands in its place is section 1: the posterior for F itself, "
         "under a prior justified by prior predictive (analyses/57) and reported "
         "under two priors leaning opposite ways.", "",
         "## MCMC diagnostics", "",
         *diag_lines, ""]

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

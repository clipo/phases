import numpy as np

from mls_emergence.inference import (
    fst_from_frequencies, sample_fst, fst_summary,
    gini_simpson_summary, gini_simpson_posterior, bayes_factor_structure,
)


def test_fst_zero_for_identical_groups():
    p = np.array([[0.5, 0.5], [0.5, 0.5]])
    assert abs(fst_from_frequencies(p, sizes=[1, 1])) < 1e-12


def test_fst_one_for_disjoint_groups():
    p = np.array([[1.0, 0.0], [0.0, 1.0]])
    # H_within = 0 for both; H_T = Gini-Simpson of [0.5,0.5] = 0.5 -> F_ST = 1
    assert abs(fst_from_frequencies(p, sizes=[1, 1]) - 1.0) < 1e-12


def test_fst_intermediate_and_size_weighting():
    p = np.array([[0.8, 0.2], [0.2, 0.8]])
    val = fst_from_frequencies(p, sizes=[1, 1])
    assert 0.0 < val < 1.0
    # a size-imbalanced version differs (pool shifts toward the larger group)
    val2 = fst_from_frequencies(p, sizes=[9, 1])
    assert val2 != val


def test_bn_recovers_high_fst():
    # near-disjoint groups -> high F_ST parameter
    gc = np.array([[90, 5, 5], [5, 5, 90]])
    idata = sample_fst(gc, draws=400, tune=600, chains=2, random_seed=1)
    s = fst_summary(idata)
    assert s["fst_median"] > 0.3
    assert s["fst_hdi95"][0] < s["fst_median"] < s["fst_hdi95"][1]
    assert s["rhat"] < 1.1


def test_bn_recovers_low_fst():
    # nearly identical groups -> low F_ST parameter
    gc = np.array([[400, 300, 300], [380, 310, 310]])
    idata = sample_fst(gc, draws=400, tune=600, chains=2, random_seed=2)
    s = fst_summary(idata)
    assert s["fst_median"] < 0.15


def test_fst_summary_keys():
    gc = np.array([[60, 20, 20], [20, 20, 60]])
    idata = sample_fst(gc, draws=300, tune=400, chains=2, random_seed=3)
    s = fst_summary(idata)
    for k in ["fst_median", "fst_mean", "fst_hdi95", "rhat", "ess", "fst_samples"]:
        assert k in s


def test_gini_simpson_readout_matches_plugin_and_brackets_it():
    from mls_emergence.signatures.variance import cultural_fst
    gc = np.array([[80, 10, 10], [10, 10, 80]])
    idata = sample_fst(gc, draws=500, tune=600, chains=2, random_seed=4)
    g = gini_simpson_summary(idata, gc)
    # the derived readout reports the SAME frequentist estimator as the plug-in
    assert abs(g["plugin_fst"] - cultural_fst(gc.astype(float))) < 1e-12
    # the plug-in should fall inside the 95% credible interval of the readout
    assert g["gst_hdi95"][0] <= g["plugin_fst"] <= g["gst_hdi95"][1]
    # readout is bounded in [0, 1]
    samp = gini_simpson_posterior(idata, gc)
    assert samp.min() >= 0.0 and samp.max() <= 1.0


def test_gini_simpson_readout_is_reproducible():
    gc = np.array([[70, 20, 10], [10, 20, 70]])
    idata = sample_fst(gc, draws=300, tune=400, chains=2, random_seed=5)
    a = gini_simpson_posterior(idata, gc, seed=0)
    b = gini_simpson_posterior(idata, gc, seed=0)
    assert np.allclose(a, b)


def test_bayes_factor_favors_structure_when_groups_differ():
    gc = np.array([[95, 3, 2], [2, 3, 95]])
    bf = bayes_factor_structure(gc, draws=800, chains=2, seed=1)
    # strongly differentiated groups -> positive evidence for structure
    assert bf["two_ln_bf"] > 2.0
    assert "structure" in bf["evidence"]


def test_bayes_factor_favors_panmixia_when_groups_identical():
    # two clusters drawn from one shared distribution -> evidence for panmixia
    gc = np.array([[40, 30, 30], [41, 29, 30]])
    bf = bayes_factor_structure(gc, draws=800, chains=2, seed=2)
    assert bf["two_ln_bf"] < 2.0

import numpy as np

from mls_emergence.inference import (
    fst_from_frequencies, sample_fst, fst_summary,
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


def test_recovers_high_fst():
    rng = np.random.default_rng(0)
    # near-disjoint groups -> high F_ST
    gc = np.array([[90, 5, 5], [5, 5, 90]])
    idata = sample_fst(gc, draws=400, tune=600, chains=2, random_seed=1)
    s = fst_summary(idata, sizes=gc.sum(1), group_counts=gc)
    assert s["fst_mean"] > 0.4
    assert s["fst_hdi95"][0] < s["fst_mean"] < s["fst_hdi95"][1]


def test_recovers_low_fst():
    # nearly identical groups -> low F_ST
    gc = np.array([[40, 30, 30], [38, 31, 31]])
    idata = sample_fst(gc, draws=400, tune=600, chains=2, random_seed=2)
    s = fst_summary(idata, sizes=gc.sum(1), group_counts=gc)
    assert s["fst_mean"] < 0.15


def test_summary_keys_and_plugin():
    from mls_emergence.signatures.variance import cultural_fst
    gc = np.array([[60, 20, 20], [20, 20, 60]])
    idata = sample_fst(gc, draws=300, tune=400, chains=2, random_seed=3)
    s = fst_summary(idata, sizes=gc.sum(1), group_counts=gc)
    for k in ["fst_mean", "fst_hdi95", "p_fst_gt0", "plugin_fst", "fst_samples"]:
        assert k in s
    assert abs(s["plugin_fst"] - cultural_fst(gc.astype(float))) < 1e-12

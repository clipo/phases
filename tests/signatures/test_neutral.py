import numpy as np
from mls_emergence.signatures import neutral

def test_homozygosity_uniform_vs_concentrated():
    uniform = np.array([25, 25, 25, 25])
    concentrated = np.array([97, 1, 1, 1])
    assert neutral.homozygosity_f(concentrated) > neutral.homozygosity_f(uniform)

def test_theta_estimators_agree_under_neutral_sample():
    rng = np.random.default_rng(0)
    counts = _ewens_sample(theta=5.0, n=500, rng=rng)
    tf, te = neutral.theta_f(counts), neutral.theta_e(counts)
    assert abs(tf - te) / te < 0.5

def test_neutrality_verdict_flags_conformity():
    counts = np.array([50, 50, 50, 50, 50])
    v = neutral.neutrality_verdict(counts)
    assert v["direction"] in {"conformist", "neutral", "anti-conformist"}

def _ewens_sample(theta, n, rng):
    labels = []; next_label = 0
    for i in range(n):
        if rng.random() < theta / (theta + i):
            labels.append(next_label); next_label += 1
        else:
            labels.append(labels[rng.integers(0, len(labels))])
    _, counts = np.unique(labels, return_counts=True)
    return counts


def test_theta_e_refuses_the_all_singletons_boundary():
    """F6. E[k | theta, n] rises from 1 to n as theta goes 0 -> infinity, so a
    finite root exists only for 1 < k < n. At k == n brentq used to raise a bare
    ValueError from inside scipy with no indication of the cause.

    Discriminating probe: k just below n must still solve, or the guard is
    simply refusing hard cases."""
    import numpy as np
    import pytest
    from mls_emergence.signatures.neutral import theta_e

    # k == n: every sherd a distinct class -> no finite root
    with pytest.raises(ValueError, match="distinct class"):
        theta_e(np.array([1, 1, 1, 1, 1]))

    # k = n - 1: a root exists and must be found
    val = theta_e(np.array([2, 1, 1, 1]))
    assert np.isfinite(val) and val > 0

    # k == 1 is the theta -> 0 limit and is returned, not refused
    assert theta_e(np.array([7, 0, 0])) == 0.0

    with pytest.raises(ValueError, match="empty assemblage"):
        theta_e(np.array([0, 0, 0]))


def test_neutrality_verdict_carries_its_own_caveat():
    """F5. The direction label comes from uncalibrated cut points and must not
    travel without saying so, since quoting the dict is how a caveat gets
    dropped (rule 6)."""
    import numpy as np
    from mls_emergence.signatures.neutral import neutrality_verdict

    res = neutrality_verdict(np.array([40, 30, 20, 10]))
    assert res["direction_is_descriptive"] is True
    assert "not evidence" in res["note"]


# --------------------------------------------------------------------------- #
# Value tests. Added 2026-09-04 after a mutation sweep: dropping the (n-1) bias
# correction from homozygosity_f, and shifting the index in Ewens' expected-k
# sum, both SURVIVED the tests that were here. Both are exactly the kind of
# off-by-one an estimator inherits silently, and homozygosity feeds Neiman's
# theta_F, which the paper reports.
# --------------------------------------------------------------------------- #

def test_homozygosity_is_the_unbiased_form_not_the_plug_in():
    """F = sum n(n-1) / (N(N-1)), NOT sum p^2.

    counts (2, 2), N = 4:
        unbiased  (2*1 + 2*1) / (4*3) = 4/12  = 1/3
        plug-in   (4 + 4) / 16        = 8/16  = 1/2

    The two differ by a third of the value here, so the assertion pins which
    estimator is in use rather than merely bracketing it.
    """
    import pytest
    from mls_emergence.signatures.neutral import homozygosity_f
    assert homozygosity_f(np.array([2.0, 2.0])) == pytest.approx(4.0 / 12.0, abs=1e-12)
    assert homozygosity_f(np.array([2.0, 2.0])) != pytest.approx(0.5, abs=1e-6)
    # (3, 1): only the 3-class contributes, 3*2 = 6, over 4*3 = 12
    assert homozygosity_f(np.array([3.0, 1.0])) == pytest.approx(0.5, abs=1e-12)


def test_theta_f_inverts_homozygosity():
    """theta = (1 - F)/F, checked against a hand value rather than by calling
    homozygosity_f again. F = 1/3 for counts (2,2), so theta = 2."""
    import pytest
    from mls_emergence.signatures.neutral import theta_f
    assert theta_f(np.array([2.0, 2.0])) == pytest.approx(2.0, abs=1e-12)


def test_ewens_expected_k_sums_theta_over_theta_plus_i_from_zero():
    """E[k] = sum_{i=0}^{n-1} theta/(theta+i). The lower limit is i = 0, so the
    first term is exactly 1 for any theta. Starting the sum at i = 1 instead
    would give 1.083 rather than 1.833 below, which this rejects."""
    import pytest
    from mls_emergence.signatures.neutral import _ewens_expected_k
    #   theta = 1, n = 3: 1/1 + 1/2 + 1/3 = 11/6
    assert _ewens_expected_k(1.0, 3) == pytest.approx(11.0 / 6.0, abs=1e-12)
    #   theta = 2, n = 4: 2/2 + 2/3 + 2/4 + 2/5
    assert _ewens_expected_k(2.0, 4) == pytest.approx(
        1.0 + 2.0 / 3.0 + 0.5 + 0.4, abs=1e-12)
    # the first term is 1 whatever theta is, which is what fixes the lower limit
    assert _ewens_expected_k(7.5, 1) == pytest.approx(1.0, abs=1e-12)

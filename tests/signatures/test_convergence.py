import pandas as pd
from mls_emergence.signatures import convergence as C

def test_convergence_high_when_all_signatures_extreme():
    panel = pd.DataFrame({
        "neutral_departure": [0.0, 1.0, 2.0],
        "seriation_groups":  [1.0, 2.0, 3.0],
        "fst":               [0.0, 0.3, 0.6],
        "assortativity":     [0.1, 0.4, 0.8],
    }, index=["early", "mid", "late"])
    score = C.convergence_score(panel)
    assert score["late"] > score["early"]

def test_time_derivative_positive_for_rising_panel():
    s = pd.Series([0.1, 0.2, 0.5, 0.9], index=[1, 2, 3, 4])
    assert C.time_derivative(s) > 0


# --------------------------------------------------------------------------- #
# Value tests. Added 2026-09-04 after a mutation sweep: scaling time_derivative
# by 1.5 and halving convergence_score both SURVIVED. That matters because
# validation/harness.discriminates compares the standardized slope against
# deriv_threshold, so a scale error moves the convergence verdict itself.
# --------------------------------------------------------------------------- #

def test_time_derivative_is_the_ols_slope_in_units_per_index_step():
    """A straight line of slope m must return exactly m. The scale is the whole
    point: discriminates() thresholds this value, so a constant factor is not a
    cosmetic difference."""
    import pandas as pd, pytest
    from mls_emergence.signatures.convergence import time_derivative
    assert time_derivative(pd.Series([0.0, 1.0, 2.0, 3.0, 4.0])) == pytest.approx(1.0, abs=1e-12)
    assert time_derivative(pd.Series([0.0, 2.0, 4.0, 6.0, 8.0])) == pytest.approx(2.0, abs=1e-12)
    assert time_derivative(pd.Series([5.0, 4.0, 3.0, 2.0, 1.0])) == pytest.approx(-1.0, abs=1e-12)
    # flat is exactly zero, and an intercept shift changes nothing
    assert time_derivative(pd.Series([3.0, 3.0, 3.0])) == pytest.approx(0.0, abs=1e-12)
    assert time_derivative(pd.Series([100.0, 101.0, 102.0])) == pytest.approx(1.0, abs=1e-12)


def test_time_derivative_uses_the_index_not_the_position():
    """Unevenly spaced index values must be honoured, since the seriation axis
    is not always 0..n-1. y = 2x on x = (0, 1, 4) has slope 2 either way only if
    the index is read; treating it as 0,1,2 would give a different answer."""
    import pandas as pd, pytest
    from mls_emergence.signatures.convergence import time_derivative
    s = pd.Series([0.0, 2.0, 8.0], index=[0.0, 1.0, 4.0])
    assert time_derivative(s) == pytest.approx(2.0, abs=1e-12)


def test_convergence_score_is_the_mean_of_z_scores():
    """Population z-scores (ddof=0) averaged across columns. Hand-computed.

    A column (1, 2, 3) has mean 2 and population SD sqrt(2/3), so its z-scores
    are (-c, 0, +c) with c = 1/sqrt(2/3). Two identical columns average to the
    same thing, which pins the scale that halving the result would break.
    """
    import numpy as np, pandas as pd, pytest
    from mls_emergence.signatures.convergence import convergence_score
    panel = pd.DataFrame({"a": [1.0, 2.0, 3.0], "b": [1.0, 2.0, 3.0]})
    c = 1.0 / np.sqrt(2.0 / 3.0)
    got = convergence_score(panel)
    assert got.iloc[0] == pytest.approx(-c, abs=1e-12)
    assert got.iloc[1] == pytest.approx(0.0, abs=1e-12)
    assert got.iloc[2] == pytest.approx(+c, abs=1e-12)


def test_convergence_score_averages_rather_than_letting_one_column_dominate():
    """The docstring's claim, made testable: one extreme column must not carry
    the score. With one rising and one falling column the mean is zero."""
    import pandas as pd, pytest
    from mls_emergence.signatures.convergence import convergence_score
    panel = pd.DataFrame({"up": [1.0, 2.0, 3.0], "down": [3.0, 2.0, 1.0]})
    for v in convergence_score(panel):
        assert v == pytest.approx(0.0, abs=1e-12)

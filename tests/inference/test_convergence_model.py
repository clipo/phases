import numpy as np

from mls_emergence.inference import (
    sample_convergence, convergence_summary,
)


def _synthetic_panel(slope_panel, slope_ser, T=6, n_panel=3, noise=0.15, seed=0):
    """Panel z-values with a common linear slope per signature plus noise."""
    rng = np.random.default_rng(seed)
    t = (np.arange(T) - (T - 1) / 2)
    t = t / t.std()
    y = np.array([slope_panel * t + rng.normal(0, noise, T) for _ in range(n_panel)])
    # standardize each row to z (mean 0, sd 1) and rescale se accordingly
    se = np.full((n_panel, T), noise)
    b_ser_obs = slope_ser
    se_ser = 0.15
    return y, se, t, b_ser_obs, se_ser


def test_recovers_joint_positive_convergence():
    y, se, t, bso, ses = _synthetic_panel(slope_panel=1.0, slope_ser=1.0, seed=1)
    idata = sample_convergence(y, se, t, bso, ses,
                               draws=500, tune=500, chains=2, random_seed=1)
    s = convergence_summary(idata)
    assert s["p_convergence"] > 0.85          # all four slopes clearly positive
    assert all(b > 0 for b in s["b_mean"])
    assert s["cscore_slope_p_pos"] > 0.9


def test_low_convergence_when_slopes_negative():
    y, se, t, bso, ses = _synthetic_panel(slope_panel=-1.0, slope_ser=-1.0, seed=2)
    idata = sample_convergence(y, se, t, bso, ses,
                               draws=500, tune=500, chains=2, random_seed=2)
    s = convergence_summary(idata)
    assert s["p_convergence"] < 0.10          # not converging


def test_mixed_slopes_give_intermediate_pconv():
    # three panel signatures up, seriation down -> NOT all four positive
    y, se, t, bso, ses = _synthetic_panel(slope_panel=1.0, slope_ser=-1.0, seed=3)
    idata = sample_convergence(y, se, t, bso, ses,
                               draws=500, tune=500, chains=2, random_seed=3)
    s = convergence_summary(idata)
    assert s["p_convergence"] < 0.5           # the down signature blocks convergence
    assert s["b_mean"][-1] < 0                 # seriation slope recovered negative


def test_non_centered_model_samples_without_divergences():
    """At the library default this model must sample cleanly.

    Previously asserted ndiv <= 5, which tolerated the symptom the test is named
    for and passed until a library version bump moved the count from 5 to 7.
    The default is now 0.99 and the assertion is zero. The other side: a nonzero
    count here means the tau funnel has returned and the default is too loose
    again, not that the tolerance should be raised.
    """
    import arviz as az
    y, se, t, bso, ses = _synthetic_panel(slope_panel=1.0, slope_ser=1.0, seed=5)
    idata = sample_convergence(y, se, t, bso, ses,
                               draws=500, tune=1000, chains=2, random_seed=5)
    ndiv = int(idata.sample_stats["diverging"].sum())
    assert ndiv == 0, f"{ndiv} divergences at the library default; funnel not resolved"
    rhat = float(az.rhat(idata, var_names=["tau"])["tau"].values)
    assert rhat < 1.01, f"tau R-hat {rhat:.4f} exceeds 1.01"


def test_summary_keys_and_shapes():
    y, se, t, bso, ses = _synthetic_panel(0.5, 0.5, seed=4)
    idata = sample_convergence(y, se, t, bso, ses,
                               draws=300, tune=300, chains=2, random_seed=4)
    s = convergence_summary(idata)
    for k in ["p_convergence", "b_mean", "b_p_pos", "b_hdi95", "mu_mean",
              "mu_p_pos", "cscore_slope_mean", "cscore_slope_p_pos"]:
        assert k in s
    assert len(s["b_mean"]) == 4 and len(s["b_hdi95"]) == 4

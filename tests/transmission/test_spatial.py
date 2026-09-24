import numpy as np
import pytest
from mls_emergence.transmission.spatial import copying_weights, drift_record, sample_record


def test_leak_one_is_exactly_the_drift_network():
    d = np.array([[0, 2, 3], [2, 0, 1], [3, 1, 0]])
    np.testing.assert_array_equal(copying_weights(d), copying_weights(d, labels=[0, 0, 1], leak=1))
    w = copying_weights(d, labels=[0, 0, 1], leak=0)
    assert w[0, 2] == 0 and w[2, 2] == 1


def test_neutral_absorbing_state_without_innovation():
    w = copying_weights(np.ones((3, 3)))
    f = drift_record(w, initial="monomorphic", innovation=0, steps=20,
                     n_records=10, burnin=10, k=4)
    assert (f[:, :, 0] == 1).all()
    assert (f[:, :, 1:] == 0).all()


def test_full_window_even_for_earliest_sample_and_totals_preserved():
    f = np.zeros((5, 2, 2)); f[:, :, 0] = 1
    f[0, :, :] = [0, 1]
    totals = np.array([10000, 23])
    out = sample_record(f, [0, 1], totals, np.random.default_rng(2), window=2)
    # Early sample includes BOTH first slices rather than just the first.
    assert .47 < out[0, 0] / totals[0] < .53
    assert out[1, 1] == 0
    np.testing.assert_array_equal(out.sum(1), totals)


def test_uniform_target_matches_default_and_skewed_target_shifts_pooled_profile():
    w = copying_weights(np.ones((4, 4)))
    a = drift_record(w, k=5, seed=3, steps=30, n_records=5, burnin=5)
    b = drift_record(w, k=5, seed=3, steps=30, n_records=5, burnin=5, target=np.ones(5))
    np.testing.assert_array_equal(a, b)
    target = np.array([0.7, 0.2, 0.05, 0.03, 0.02])
    f = drift_record(w, k=5, seed=3, steps=400, n_records=40, burnin=200,
                     innovation=0.05, mixing=0.05, target=target)
    pooled = f.mean((0, 1))
    assert pooled[0] > 0.5 and pooled[0] > pooled[1] > pooled[2:].max()
    with pytest.raises(ValueError):
        drift_record(w, k=5, target=np.ones(4))


# --------------------------------------------------------------------------- #
# Heterogeneity across nodes. Added 2026-09-16 with the three departures tested
# in analyses/64 and 65. Each of these was a scalar before, and each scalar was
# a substantive assumption rather than a convenience: equal populations, equal
# accumulation spans, and one region-wide innovation profile that no boundary
# interrupts.
# --------------------------------------------------------------------------- #

def _ring(n=6):
    import numpy as np
    w = np.ones((n, n))
    np.fill_diagonal(w, 0)
    return w / w.sum(1, keepdims=True)


def test_per_node_population_changes_drift_rate_not_just_shape():
    """Drift goes as 1/N, so a small node must lose diversity faster than a large
    one in the SAME run. The probe: compare the smallest and largest node under a
    per-node vector. If the vector were ignored and one scalar used, the two
    would be exchangeable and this ordering would appear only by chance."""
    import numpy as np
    from mls_emergence.transmission.spatial import drift_record
    w = _ring(6)
    n_ind = np.array([20, 20, 20, 5000, 5000, 5000])
    rec = drift_record(w, k=8, n_ind=n_ind, steps=400, n_records=40,
                       burnin=200, seed=3)
    final = rec[-1]
    gs = 1.0 - (final ** 2).sum(1)
    small, large = gs[:3].mean(), gs[3:].mean()
    assert small < large, (
        f"small populations should hold less diversity, got {small:.3f} vs {large:.3f}")


def test_per_node_innovation_target_separates_groups():
    """Two groups innovating from different profiles must end up differing. With
    one shared profile they must not, beyond drift. The second assertion is the
    probe that makes the first mean something."""
    import numpy as np
    from mls_emergence.transmission.spatial import drift_record
    w = _ring(6)
    shared = np.array([.4, .3, .2, .1])
    split = np.vstack([np.tile(shared, (3, 1)),
                       np.tile(shared[::-1], (3, 1))])
    kw = dict(k=4, n_ind=4000, steps=600, n_records=20, burnin=300,
              innovation=0.05, mixing=0.01, seed=5)
    apart = drift_record(w, target=split, **kw)[-1]
    together = drift_record(w, target=shared, **kw)[-1]
    d_apart = abs(apart[:3].mean(0) - apart[3:].mean(0)).sum()
    d_together = abs(together[:3].mean(0) - together[3:].mean(0)).sum()
    assert d_apart > 3 * d_together, (
        f"group-specific innovation should separate the groups: "
        f"{d_apart:.3f} against {d_together:.3f} when shared")


def test_per_site_window_averages_different_spans():
    """A long window averages more of the sequence than a short one, so on a
    record that changes through time the two must differ. Built on a record with
    a deliberate trend so that averaging has something to average over."""
    import numpy as np
    from mls_emergence.transmission.spatial import sample_record
    t, n, k = 20, 4, 3
    rec = np.zeros((t, n, k))
    for i in range(t):                      # class 0 falls, class 2 rises
        frac = i / (t - 1)
        rec[i, :, 0] = 1 - frac
        rec[i, :, 2] = frac
        rec[i, :, 1] = 0.0
    rec += 1e-9
    rec /= rec.sum(-1, keepdims=True)
    ranks = np.full(n, 0.5)
    totals = np.full(n, 4000)
    rng = np.random.default_rng(0)
    short = sample_record(rec, ranks, totals, rng, window=np.array([1, 1, 1, 1]))
    rng = np.random.default_rng(0)
    long_ = sample_record(rec, ranks, totals, rng, window=np.array([15, 15, 15, 15]))
    # The long window reaches further back, so it carries more of the early class.
    assert long_[:, 0].mean() > short[:, 0].mean(), (
        "a longer accumulation span should retain more of the earlier repertoire")


@pytest.mark.parametrize("bad, match", [
    ("pop", "one value per node"),
    ("target", "per-node target"),
    ("window", "one value per node"),
])
def test_heterogeneous_inputs_guard_their_length(bad, match):
    """Rule 5. A vector of the wrong length is the likely mistake here, and it
    would otherwise broadcast or truncate into a plausible-looking result."""
    import numpy as np
    from mls_emergence.transmission.spatial import drift_record, sample_record
    w = _ring(6)
    if bad == "pop":
        with pytest.raises(ValueError, match=match):
            drift_record(w, k=4, n_ind=np.array([10, 20]), steps=40,
                         n_records=4, burnin=10, seed=0)
    elif bad == "target":
        with pytest.raises(ValueError, match=match):
            drift_record(w, k=4, n_ind=100, steps=40, n_records=4, burnin=10,
                         seed=0, target=np.ones((2, 4)) / 4)
    else:
        rec = drift_record(w, k=4, n_ind=100, steps=40, n_records=4, burnin=10, seed=0)
        with pytest.raises(ValueError, match=match):
            sample_record(rec, np.full(6, .5), np.full(6, 100),
                          np.random.default_rng(0), window=np.array([2, 3]))

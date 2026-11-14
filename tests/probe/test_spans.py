"""Paired spans and named streams (spec 3.3, 3.4; R2-4)."""
import numpy as np
import pytest

from mls_emergence.probe import spans


def test_base_key_and_named_streams_are_independent_of_consumption_order():
    k = spans.base_key("S2", 3)
    assert k == 2026101600 + 10000 * 2 + 3
    a1 = spans.stream(k, "lengths").integers(0, 10**9, 5)
    # Consuming another stream first does not change this one.
    spans.stream(k, "history").integers(0, 10**9, 1000)
    a2 = spans.stream(k, "lengths").integers(0, 10**9, 5)
    assert np.array_equal(a1, a2)
    # Different names, same key: different draws. Same name, different setting: different draws.
    assert not np.array_equal(a1, spans.stream(k, "history").integers(0, 10**9, 5))
    e0 = spans.stream(k, "ends", setting=0).integers(0, 10**9, 5)
    e1 = spans.stream(k, "ends", setting=1).integers(0, 10**9, 5)
    assert not np.array_equal(e0, e1)
    with pytest.raises(ValueError):
        spans.stream(k, "lengths", setting=1)
    with pytest.raises(ValueError):
        spans.stream(k, "ends")
    with pytest.raises(ValueError):
        spans.stream(k, "nonsense")


def test_t0_and_t1_share_lengths_and_differ_in_ends():
    k = spans.base_key("S0", 0); S, L_d, n = 4000, 300, 28
    w0 = spans.draw_lengths(spans.stream(k, "lengths"), n, L_d)
    w1 = spans.draw_lengths(spans.stream(k, "lengths"), n, L_d)
    assert np.array_equal(w0, w1) and w0.min() >= 1 and w0.max() <= L_d
    e0 = spans.ends_t0(n, S); e1 = spans.ends_t1(spans.stream(k, "ends", setting=1), n, S, L_d)
    assert (e0 == S - 1).all()
    assert e1.min() >= S - L_d and e1.max() <= S - 1 and not np.array_equal(e0, e1)


def test_area_offsets_are_one_per_area():
    rng = np.random.default_rng(1)
    areas = np.array(["Parkin"] * 12 + ["Kent"] * 11 + ["Walls"] * 5)
    off = spans.area_offsets(rng, areas, L_R=1000)
    for a in set(areas):
        assert len(set(off[areas == a].tolist())) == 1
    assert off.min() >= 0 and off.max() <= 1000


def test_sample_spans_realizes_the_drawn_ends_exactly():
    # A record whose value at step t is t itself, so the sampled mean reveals
    # the window: mean over [e - w + 1, e] equals e - (w - 1) / 2.
    S, n = 500, 7
    t = np.arange(S, dtype=float)
    record = np.stack([np.stack([t / S, 1 - t / S], axis=1)] * n, axis=1)   # S x n x 2
    rng = np.random.default_rng(0)
    lengths = rng.integers(1, 101, size=n)
    ends = rng.integers(100, S, size=n)
    ranks = spans.ranks_for_ends(ends, lengths, S)
    realized = lengths - 1 + np.rint(ranks * (S - lengths)).astype(int)
    assert np.array_equal(realized, ends)
    # Large totals make the sampled share close to the window mean.
    totals = np.full(n, 10**6)
    y = spans.sample_spans(record, ends, lengths, totals, np.random.default_rng(1))
    share = y[:, 0] / totals
    expected = (ends - (lengths - 1) / 2) / S
    assert np.allclose(share, expected, atol=2e-3)
    # Spans that do not fit are refused.
    with pytest.raises(ValueError):
        spans.ranks_for_ends(np.array([5]), np.array([10]), S)
    # T0 guard: a span as long as the whole record ending at S - 1 has rank 0 and still realizes.
    r = spans.ranks_for_ends(np.array([S - 1]), np.array([S]), S)
    assert r[0] == 0.0


def test_ends_window_range_and_domain():
    # C1 (spec 6.1): each end uniform over the last W stored steps; a window longer than the record is refused.
    rng = np.random.default_rng(0)
    e = spans.ends_window(rng, 200000, S=2920, W=2336)
    assert e.min() >= 2920 - 2336 and e.max() < 2920
    assert abs((e < 2920 - 1168).mean() - 0.5) < 0.01      # uniform over the window, not over its last half
    with pytest.raises(ValueError, match="window exceeds"):
        spans.ends_window(rng, 5, S=100, W=101)


def test_span_profile_matches_the_window_sample_record_uses():
    # R2 draws clumps from the span-averaged profile; it must be the same window sample_record averages.
    rec = np.zeros((10, 2, 3)); rec[:, 0, 0] = np.arange(10); rec[:, 0, 1] = 1.0; rec[:, 1, 2] = 1.0
    p = spans.span_profile(rec, ends=np.array([9, 4]), lengths=np.array([2, 5]))
    assert np.allclose(p[0], np.array([8.5, 1.0, 0.0]) / 9.5) and np.allclose(p[1], [0, 0, 1])
    with pytest.raises(ValueError, match="inside the stored record"):
        spans.span_profile(rec, ends=np.array([9, 3]), lengths=np.array([2, 5]))


def test_ends_window_rejects_a_non_integer_window():
    import pytest
    with pytest.raises(ValueError):
        spans.ends_window(np.random.default_rng(0), 5, 100, 2.5)
